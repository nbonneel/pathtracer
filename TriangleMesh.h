#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "Geometry.h"


template<typename T>
class BVHNodesT {
public:
	//int i0, i1;
	bool isleaf;
	int fg, fd;
	BBoxT<T> bbox;
};

template<typename T>
class BVHT {
public:
	BBoxT<T> bbox;
	std::vector<BVHNodesT<T> > nodes;
};

typedef BVHT<double> BVHd;
typedef BVHT<float> BVHf;
typedef BVHf BVH;

typedef BVHNodesT<double> BVHNodesd;
typedef BVHNodesT<float> BVHNodesf;
typedef BVHNodesf BVHNodes;


class Edge {
public:
	Edge(int i = 0, int j = 0) {
		if (i<j) {
			a = i; b = j;
		} else {
			a = j; b = i;
		}
	}
	int a, b;
	bool operator<(const Edge& r) const {
		if (a < r.a) { return true; }
		if (a > r.a) { return false; }
		return (b < r.b);
	}
	bool operator==(const Edge& r) const {
		if (a != r.a) return false;
		if (b != r.b) return false;
		return true;
	}
};

class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
		showEdges[0] = true;
		showEdges[1] = true;
		showEdges[2] = true;
	};
	int vtxi, vtxj, vtxk;
	int uvi, uvj, uvk;
	int ni, nj, nk;
	int group;
	bool showEdges[3];
};

class Triangle {
public:
	Triangle() {};
	Triangle(const Vector& A, const Vector &B, const Vector& C) : A(A) {
		u = B - A;
		v = C - A;
		N = cross(u, v);
		m11 = u.getNorm2();
		m22 = v.getNorm2();
		m12 = dot(u, v);
		invdetm = 1. / (m11*m22 - m12*m12);
	};
	float area() {
		return sqrt(N.getNorm2())*0.5;
	}
	bool intersection(const Ray& d, Vector& P, float &t, float &alpha, float &beta, float &gamma) const {
		t = dot(A - d.origin, N) / dot(d.direction, N);
		if (t < 0 || t != t) return false;  //isnan

		P = d.origin + t*d.direction;
		Vector w(P - A);
		float b11 = dot(w, u);
		float b21 = dot(w, v);
		float detb = b11*m22 - b21*m12;
		beta = detb * invdetm;    // coord barycentrique w.r.t � B
		if (beta < 0) return false;

		float detg = b21*m11 - b11*m12;
		gamma = detg * invdetm;   // coord barycentrique w.r.t � C
		if (gamma < 0) return false;

		alpha = 1 - beta - gamma;
		if (alpha < 0) return false;

		//N.normalize(); // only normalized at the last moment, for the actual nearest triangle being intersected (not all intersected triangles)
		return true;

	}

	Vector A, u, v, N;
	float m11, m12, m22, invdetm;
	float uvs[3][2];
	Vector normals[3];

};

Vector TransformH(const Vector &in, float H); // transform hue of a color

class TriMesh : public Object {
public:
  ~TriMesh() {}
	TriMesh() { type = OT_TRIMESH; };
	TriMesh(Scene* scene, const char* obj, float scaling, const Vector& offset, bool mirror = false, const char* colors_csv_filename = NULL, bool preserve_input = false, bool center = true, Vector rot_center = Vector(std::nan(""), std::nan(""), std::nan("")));

	void init(Scene* scene, const char* obj, float scaling, const Vector& offset, bool mirror = false, const char* colors_csv_filename = NULL, bool load_textures = true, bool preserve_input = false, bool center = true, Vector rot_center = Vector(std::nan(""), std::nan(""), std::nan("")));

	void readOBJ(const char* filename, bool load_textures);
	void readVRML(const char* filename);
	void readOFF(const char* filename);
	void exportMTL(const char* mtlfile);
	void saveOBJ(const char* obj);

	int getNbConnected(int &alsoReturnsNbEdges, int &nonManifoldFaces, int &nbBoundaryEdges);
	void findQuads(int &nbTriangles, int &nbOthers, int &nbRealEdges); //nbOthers = nb of quads or higher order polygons
	//void load_normal_map(const char* filename);  
	void save_to_file(FILE* f) {

		fprintf(f, "NEW MESH\n");

		Object::save_to_file(f);

		fprintf(f, "is_centered: %u\n", is_centered ? 1 : 0);
		fprintf(f, "has_csv: %u\n", (csv_file.size() != 0) ? 1 : 0);
		fprintf(f, "csv_file: %s\n", csv_file.c_str());
	}

	static TriMesh* create_from_file(FILE* f, Scene* scene, const char* replacedNames = NULL) {
		TriMesh* result = new TriMesh();
		result->Object::load_from_file(f, replacedNames);
		int hascsv;
		char line[512];
		fscanf(f, "%[^\n]\n", line);
		bool is_centered = true;
		if (line[0] == 'i' && line[1] == 's') { //for backward compatibility
			int c;
			sscanf(line, "is_centered: %u\n", &c);
			is_centered = (c == 1);
			fscanf(f, "has_csv: %u\n", &hascsv);
		} else {
			sscanf(line, "has_csv: %u\n", &hascsv);
		}
		
		if (hascsv)
			fscanf(f, "csv_file: %[^\n]\n", line);
		else
			fscanf(f, "csv_file: \n");

		result->init(scene, result->name.c_str(), 1., Vector(0, 0, 0), result->miroir, hascsv ? line : NULL, false, false, is_centered, result->rotation_center);
		return result;
	}

	virtual void colorAnisotropy() {

		facecolors.resize(indices.size());
		for (int i = 0; i < indices.size(); i++) {
			float angle1 = std::abs(dot(vertices[indices[i].vtxj] - vertices[indices[i].vtxi], vertices[indices[i].vtxk] - vertices[indices[i].vtxi]))/sqrt((vertices[indices[i].vtxj] - vertices[indices[i].vtxi]).getNorm2()*(vertices[indices[i].vtxk] - vertices[indices[i].vtxi]).getNorm2());
			float angle2 = std::abs(dot(vertices[indices[i].vtxi] - vertices[indices[i].vtxj], vertices[indices[i].vtxk] - vertices[indices[i].vtxj])) / sqrt((vertices[indices[i].vtxi] - vertices[indices[i].vtxj]).getNorm2()*(vertices[indices[i].vtxk] - vertices[indices[i].vtxj]).getNorm2());
			float angle3 = std::abs(dot(vertices[indices[i].vtxi] - vertices[indices[i].vtxk], vertices[indices[i].vtxj] - vertices[indices[i].vtxk])) / sqrt((vertices[indices[i].vtxi] - vertices[indices[i].vtxk]).getNorm2()*(vertices[indices[i].vtxj] - vertices[indices[i].vtxk]).getNorm2());
			float aniso = std::acos(std::max(angle1, std::max(angle2, angle3)))*180/M_PI;
			facecolors[i] = TransformH(Vector(1,0,0), std::min(240.f, std::max(0.f, aniso/60*240)));
			//facecolors[i] = Vector(std::pow(facecolors[i][0], 2.2), std::pow(facecolors[i][1], 2.2), std::pow(facecolors[i][2], 2.2)); // just here, assumes a 2.2 gamma
		}

		// save color map to file
		std::vector<unsigned char> img(240*30 * 3, 0);
		for (int i=0; i<240; i++)
			for (int j = 0; j < 30; j++) {
				Vector col = TransformH(Vector(1, 0, 0), i);
				img[i * 3 * 30 + j * 3 + 0] = std::pow(col[0], 1./2.2) * 255;  //beware albedo values were gamma corrected. Same on the legend
				img[i * 3 * 30 + j * 3 + 1] = std::pow(col[1], 1. / 2.2) * 255;
				img[i * 3 * 30 + j * 3 + 2] = std::pow(col[2], 1. / 2.2) * 255;
			}
		save_image("legend.png", &img[0], 30, 240);
	}

	virtual void randomColors() {

		int r1 = rand() % 10000+1;
		int r2 = rand() % 10000 + 1;
		int r3 = rand() % 10000 + 1;
		for (int i = 0; i < facecolors.size(); i++) {
			Vector c = facecolors[i];
			long long ir = (long long)(c[0] * 1024); ir = (ir*r1 + ir*ir*(r1 + 1) + ir*15 +r1+ 3)%1024;
			long long ig = (long long)(c[1] * 1024); ig = (ig*r2 + ig*ig*(r2 + 9) + ig * 7 +r2+ 3) % 1024;
			long long ib = (long long)(c[2] * 1024); ib = (ib*r3 + ib*ib*(r3 + 3) + ib * 18 +r3+ 3) % 1024;
			facecolors[i] = Vector(ir / 1024., ig / 1024., ib / 1024.);
		}
	}

	std::map<std::string, int> groupNames;

	std::vector<TriangleIndices> indices;
	std::vector<Triangle> triangleSoup;
	std::vector<Vector> tangentSoup, bitangentSoup;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs; // Vector en 3D mais on n'utilise que 2 composantes
	std::vector<Vector> facecolors;
	std::vector<std::map<int, Vector> > edgecolor; // [id_vertex1][idvertex2]
	std::vector<Vector> vertexcolors;
	std::vector<int> permuted_triangle_index;

	std::vector<Vector> tangents;
	std::vector<Vector> bitangents;
	std::string csv_file;

#ifdef USE_EMBREE
	bool intersection(const Ray& d, Vector& P, float &t, MaterialValues &mat, float cur_best_t, int &triangle_id) const { return false; };
	bool intersection_shadow(const Ray& d, float &t, float cur_best_t, float dist_light) const { return false;  };
	bool reservoir_sampling_intersection(const Ray& r, Vector& P, float &t, MaterialValues &mat, int &triangle_id, int &current_nb_intersections, float min_t, float max_t) const { return false; };
#else
	bool intersection(const Ray& d, Vector& P, float &t, MaterialValues &mat, float cur_best_t, int &triangle_id) const;
	bool intersection_shadow(const Ray& d, float &t, float cur_best_t, float dist_light) const;
	bool reservoir_sampling_intersection(const Ray& r, Vector& P, float &t, MaterialValues &mat, int &triangle_id, int &current_nb_intersections, float min_t, float max_t) const;
#endif

	MaterialValues getMaterial(int triId, float alpha, float beta, float gamma, MaterialValues &mat) const;

	BBox build_bbox(int i0, int i1);
	BBox bbox;

#ifdef USE_EMBREE
	RTCScene embree_scene_for_instance;
	RTCGeometry instance_geom;
#else
	BVH bvh;	
	BBox build_centers_bbox(int i0, int i1);

	void build_bvh(BVH* node, int i0, int i1);
	void build_bvh_recur(BVH* b, int node, int i0, int i1, int depth);
#endif

	bool is_centered;


	int max_bvh_triangles;
	int bvh_depth;
	float bvh_avg_depth;
	int bvh_nb_nodes;
private:

	void load_edge_colors(const char* csvfilename); // beware, has to be called before the build bvh
	void setup_tangents();   // has to be called after the bvh construction
};




class Yarns : public Object {
public:
	Yarns() { type = OT_YARNS; };
	Yarns(const char* filename) {
		type = OT_YARNS;
		FILE* f = fopen(filename, "r+");
		int nbyarns;
		fscanf(f, "%u\n", &nbyarns);
		for (int i = 0; i < nbyarns; i++) {
			int nbsegments;
			fscanf(f, "%u\n", &nbsegments);
			float x, y, z;
			fscanf(f, "%f %f %f\n", &x, &y, &z);
			for (int j = 0; j < nbsegments - 1; j++) {
				float x2, y2, z2;
				fscanf(f, "%f %f %f\n", &x2, &y2, &z2);
				Cylinder* c = new Cylinder(Vector(x, y, z)*50.f, Vector(x2, y2, z2)*50.f, 0.1f);
				cyls.push_back(c);
				x = x2;
				y = y2;
				z = z2;
			}
		}
		fclose(f);
		build_bvh(&bvh, 0, cyls.size());
	}

	bool intersection(const Ray& d, Vector& P, float &t, MaterialValues &mat, float cur_best_t, int &triangle_id) const;
	bool intersection_shadow(const Ray& d, float &t, float cur_best_t, float dist_light) const {
		MaterialValues mat;
		int tid;
		Vector P;
		return intersection(d, P, t, mat, cur_best_t, tid);
	}
	bool reservoir_sampling_intersection(const Ray& r, Vector& P, float &t, MaterialValues &mat, int &triangle_id, int &current_nb_intersections, float min_t, float max_t) const;

	BBox build_bbox(int i0, int i1);
	BBox build_centers_bbox(int i0, int i1);

	void build_bvh(BVH* node, int i0, int i1);
	void build_bvh_recur(BVH* b, int node, int i0, int i1, int depth);

	std::vector<Cylinder*> cyls;
	BVH bvh;
};