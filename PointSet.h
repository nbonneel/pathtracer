#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "Geometry.h"
#include "nanoflann.hpp"
#include "TriangleMesh.h"

// only for NanoFlann
template <typename T>
struct PointCloud
{
	struct Point {
		T  x, y, z;
	};

	std::vector<Point>  pts;

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return pts.size(); }

	// Returns the dim'th component of the idx'th point in the class:
	// Since this is inlined and the "dim" argument is typically an immediate value, the
	//  "if/else's" are actually solved at compile time.
	inline T kdtree_get_pt(const size_t idx, int dim) const
	{
		if (dim == 0) return pts[idx].x;
		else if (dim == 1) return pts[idx].y;
		else return pts[idx].z;
	}

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }

};

class PointSet : public Object {
public:
	PointSet() {};
	PointSet(const char* filename, double radius, bool mirror = false, bool transp = false, bool normal_swapped = false) {
		init(filename, radius, mirror, transp, normal_swapped);
	};

	void init(const char* filename, double radius, bool mirror = false, bool transp = false, bool normal_swapped = false) {
		miroir = mirror;
		transparent = transp;
		this->radius = radius;
		flip_normals = normal_swapped;


		FILE* f = fopen(filename, "r");
		char* line = new char[1000];
		while (!feof(f)) {
			if (!fgets(line, 1000, f)) break;
			std::string s(line);

			Vector p;
			int r, g, b;
			int tmp1, tmp2;
			sscanf(line, "%u %u %lf %lf %lf %u %u %u", &tmp1, &tmp2, &p[0], &p[1], &p[2], &r, &g, &b);
			vertices.push_back(p);
			colors.push_back(Vector(r / 255., g / 255., b / 255.));
			Vector n = Vector(uniform(engine), uniform(engine), uniform(engine));
			n.normalize();
			normals.push_back(n);
		}
		fclose(f);

		estimate_normals();

		BBox b = build_centers_bbox(0, vertices.size());
		Vector center(0., 0., 0.);
		double s = sqrt((b.bounds[1] - b.bounds[0]).getNorm2());
		for (int i = 0; i < vertices.size(); i++) {
			vertices[i] = (vertices[i] - b.bounds[0]) / s * 50;
			center += vertices[i];
		}
		center = center / vertices.size();
		build_bvh(0, vertices.size());


		this->rotation_center = center;
		name = "PointSet";
	}

	void estimate_normals() {

		int N = vertices.size();
		PointCloud<double> cloud;
		cloud.pts.resize(N);
		for (size_t i = 0; i < N; i++) {
			cloud.pts[i].x = vertices[i][0];
			cloud.pts[i].y = vertices[i][1];
			cloud.pts[i].z = vertices[i][2];
		}
		typedef nanoflann::KDTreeSingleIndexAdaptor<
			nanoflann::L2_Simple_Adaptor<double, PointCloud<double> >,
			PointCloud<double>, 3> my_kd_tree_t;

		my_kd_tree_t   index(3 /*dim*/, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
		index.buildIndex();

		std::vector<double> allNNDist(N);
		for (int i = 0; i < N; i++) {

			const size_t num_results = 5;
			size_t ret_index[num_results];
			double out_dist_sqr[num_results];
			nanoflann::KNNResultSet<double> resultSet(num_results);
			resultSet.init(ret_index, out_dist_sqr);
			index.findNeighbors(resultSet, &vertices[i][0], nanoflann::SearchParams());

			Vector center(0., 0., 0.);
			for (int j = 0; j < num_results; j++) {
				center += vertices[ret_index[j]];
			}
			allNNDist[i] = sqrt(out_dist_sqr[1]);
			center *= 1. / num_results;
			double cov[9] = { 0,0,0,0,0,0,0,0,0 };
			for (int l = 0; l < num_results; l++) {
				Vector v = vertices[ret_index[l]] - center;
				for (int j = 0; j < 3; j++) {
					for (int k = 0; k < 3; k++) {
						cov[j * 3 + k] += v[j] * v[k];
					}
				}
			}
			cimg_library::CImg<double> mat(cov, 3, 3), val(3), vec(3, 3);
			mat.symmetric_eigen(val, vec);
			normals[i] = Vector(vec(2, 0), vec(2, 1), vec(2, 2));
		}

		std::sort(allNNDist.begin(), allNNDist.end());

		radius = allNNDist[allNNDist.size() / 2] * 0.3;

	}

	void save_to_file(FILE* f) {

		fprintf(f, "NEW POINTSET\n");

		Object::save_to_file(f);

		fprintf(f, "radius: %lf\n", radius);

	}

	static PointSet* create_from_file(FILE* f) {
		PointSet* result = new PointSet();
		result->Object::load_from_file(f);

		fscanf(f, "radius: %lf\n", &result->radius);

		result->init(result->name.c_str(), result->radius, result->miroir, result->transparent, result->flip_normals);
		return result;
	}


	bool intersection(const Ray& d, Vector& P, double &t, MaterialValues &mat, double cur_best_t, int &triangle_id) const;
	bool intersection_shadow(const Ray& d, double &t, double cur_best_t, double dist_light) const;

	BBox build_bbox(int i0, int i1);
	BBox build_centers_bbox(int i0, int i1);
	void build_bvh(int i0, int i1);
	void build_bvh_recur(BVH* b, int node, int i0, int i1);

	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> colors;
	double radius;
	BVH bvh;
};