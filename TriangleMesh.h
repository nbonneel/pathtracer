#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "Geometry.h"

class Geometry : public Object {
public:
	Geometry() {};
	Geometry(const char* obj, double scaling, const Vector& offset, bool mirror = false, bool transp = false, const char* colors_csv_filename = NULL);

	void init(const char* obj, double scaling, const Vector& offset, bool mirror = false, bool transp = false, const char* colors_csv_filename = NULL, bool load_textures = true);

	void readOBJ(const char* obj, bool load_textures);
	void readVRML(const char* obj);
	void exportMTL(const char* mtlfile);

	//void load_normal_map(const char* filename);  
	void save_to_file(FILE* f) {

		fprintf(f, "NEW MESH\n");

		Object::save_to_file(f);

		fprintf(f, "has_csv: %u\n", (csv_file.size() != 0) ? 1 : 0);
		fprintf(f, "csv_file: %s\n", csv_file.c_str());
	}

	static Geometry* create_from_file(FILE* f) {
		Geometry* result = new Geometry();
		result->Object::load_from_file(f);
		int hascsv;
		char line[512];
		fscanf(f, "has_csv: %u\n", &hascsv);
		if (hascsv)
			fscanf(f, "csv_file: %[^\n]\n", line);
		else
			fscanf(f, "csv_file: \n");

		result->init(result->name.c_str(), 1., Vector(0, 0, 0), result->miroir, result->transparent, hascsv ? line : NULL, false);
		return result;
	}


	std::map<std::string, int> groupNames;

	std::vector<TriangleIndices> indices;
	std::vector<Triangle> triangleSoup;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs; // Vector en 3D mais on n'utilise que 2 composantes
	std::vector<Vector> facecolors;
	std::vector<std::map<int, Vector> > edgecolor; // [id_vertex1][idvertex2]
	std::vector<Vector> vertexcolors;

	std::vector<Vector> tangents;
	std::vector<Vector> bitangents;
	std::string csv_file;

	bool intersection(const Ray& d, Vector& P, double &t, MaterialValues &mat, double cur_best_t, int &triangle_id) const;
	bool intersection_shadow(const Ray& d, double &t, double cur_best_t, double dist_light) const;

	BBox build_bbox(int i0, int i1);
	BBox build_centers_bbox(int i0, int i1);

	void build_bvh(BVH* node, int i0, int i1);
	void build_bvh_recur(BVH* b, int node, int i0, int i1, int depth);

	BVH bvh;
	int max_bvh_triangles;
	int bvh_depth;
	double bvh_avg_depth;
	int bvh_nb_nodes;
private:

	void load_edge_colors(const char* csvfilename); // beware, has to be called before the build bvh
	void setup_tangents();   // has to be called after the bvh construction
};