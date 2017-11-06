#pragma once

#include <iostream>
#include <vector>
#include "Vector.h"
#include <algorithm>
#include <map>
#include "utils.h"


#include "CImg.h"
#include "utils.h"

#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif


void load_bmp(const char* filename, std::vector<unsigned char> &tex, int &W, int& H);

class PointSet;

class MaterialValues {
public:
	MaterialValues() {
		shadingN = Vector(0, 1., 0);
		Kd = Vector(0.5, 0.5, 0.5);
		Ne = Vector(100, 100, 100);
		Ks = Vector(0., 0., 0.);
		Ke = Vector(0., 0., 0.);
	}
	Vector shadingN, Kd, Ks, Ne, Ke;
};

class Texture {
public:
	Texture() { multiplier = Vector(1., 1., 1.); W = 0; H = 0; type = 0; };
	Texture(const char* filename, int texType, const Vector &multiplier): type(texType) {
		W = 0;
		H = 0;		
		switch (texType) {
		case 0: loadColors(filename); break;   // diffuse defaults to white
		case 1: loadColors(filename); break;  // specular defaults to black
		case 2: loadNormals(filename); break;   // normals defaults to 0,0,1
		case 3: loadColors(filename); break;   // alpha defaults to white
		case 4: loadColors(filename); break;   // roughness defaults to 1
		}
		this->multiplier = multiplier;
	}

	static double wrap(double u) {
		u = fmod(u, 1);
		if (u < 0) u += 1;
		return u;
	}
	static Texture defaultNormal() {
		return Texture("Null", true, Vector(0.,0.,1.));
	}
	static Texture defaultSpecular() {
		return Texture("Null", false, Vector(0., 0., 0.));
	}
	static Texture defaultDiffuse() {
		return Texture("Null", false, Vector(1., 1., 1.));
	}

	void clear_texture() {
		W = 0;
		H = 0;
		filename = "Null";
		values.clear();
	}

	Vector getVec(double u, double v) const {
		if (W > 0) {
			// no assert on u and v ; assume they are btw 0  and 1
			int x = u * (W - 1);
			int y = v * (H - 1);
			int idx = (y*W + x)*3;
			double cr = values[idx]*multiplier[0];
			double cg = values[idx + 1]*multiplier[1];
			double cb = values[idx + 2]*multiplier[2];
			return Vector(cr, cg, cb);
		} else {
			return multiplier;
		}
	}
	Vector getNormal(double u, double v) const {
		if (W > 0) {
			// no assert on u and v ; assume they are btw 0  and 1
			int x = u * (W - 1);
			int y = v * (H - 1);
			int idx = (y*W + x) * 3;
			double cr = values[idx];
			double cg = values[idx + 1];
			double cb = values[idx + 2];
			return Vector(cr, cg, cb);
		} else {
			return Vector(0.,0.,1.);
		}
	}
	double getValRed(double u, double v) const {
		if (W > 0) {
			// no assert on u and v ; assume they are btw 0  and 1
			int x = u * (W - 1);
			int y = v * (H - 1);
			int idx = (y*W + x) * 3;
			double cr = values[idx] * multiplier[0];
			return cr;
		} else {
			return multiplier[0];
		}
	}
	void loadColors(const char* file) {
		if (file) {
			filename = std::string(file);
			if (load_image(file, values, W, H)) {
				for (int i = 0; i < values.size(); i++) {
					values[i] /= 255.;
				}
			} 
		} 
	}

	void loadNormals(const char* file) {
		if (file) {
			filename = std::string(file);
			if (load_image(file, values, W, H)) {
				for (int i = 0; i < values.size() / 3; i++) {
					Vector v(values[i * 3] - 128, values[i * 3 + 1] - 128, values[i * 3 + 2] - 128);
					v.normalize();
					values[i * 3] = v[0];
					values[i * 3 + 1] = v[1];
					values[i * 3 + 2] = v[2];
				}
			}
		} 
	}

	size_t W, H;
	std::vector<double> values;
	Vector multiplier;
	std::string filename;
	unsigned char type; // 0: diffuse, 1: specular, 2: normal, 3:alpha, 4: roughness
};

class BBox {
public:
	BBox() {};
	BBox(const Vector& bmin, const Vector& bmax) {
		bounds[0] = bmin;
		bounds[1] = bmax;	
	};

	double area() {
		Vector sides = bounds[1] - bounds[0];
		return 2 * (sides[0] * sides[1] + sides[0] * sides[2] + sides[1] * sides[2]);
	}

	bool intersection(const Ray& d) const {

		double t_1_x = (bounds[0][0] - d.origin[0]) / d.direction[0];
		double t_2_x = (bounds[1][0] - d.origin[0]) / d.direction[0];
		double t_min_x = std::min(t_1_x, t_2_x);
		double t_max_x = std::max(t_1_x, t_2_x);

		double t_1_y = (bounds[0][1] - d.origin[1]) / d.direction[1];
		double t_2_y = (bounds[1][1] - d.origin[1]) / d.direction[1];
		double t_min_y = std::min(t_1_y, t_2_y);
		double t_max_y = std::max(t_1_y, t_2_y);

		double t_1_z = (bounds[0][2] - d.origin[2]) / d.direction[2];
		double t_2_z = (bounds[1][2] - d.origin[2]) / d.direction[2];
		double t_min_z = std::min(t_1_z, t_2_z);
		double t_max_z = std::max(t_1_z, t_2_z);

		double tlast = std::min(std::min(t_max_x, t_max_y), t_max_z);
		if (tlast < 0) return false;

		if (tlast - std::max(std::max(t_min_x, t_min_y), t_min_z) >= -0.000001) return true;
		return false;
	} 
	 
	bool intersection(const Ray& d, double &t) const { 
		 
		double invdx = 1. / d.direction[0];
		double invdy = 1. / d.direction[1];
		double invdz = 1. / d.direction[2];
		double t_1_x = (bounds[0][0] - d.origin[0]) *invdx;
		double t_2_x = (bounds[1][0] - d.origin[0]) *invdx;
		double t_min_x = std::min(t_1_x, t_2_x);
		double t_max_x = std::max(t_1_x, t_2_x);

		double t_1_y = (bounds[0][1] - d.origin[1]) *invdy;
		double t_2_y = (bounds[1][1] - d.origin[1]) *invdy;
		double t_min_y = std::min(t_1_y, t_2_y);
		double t_max_y = std::max(t_1_y, t_2_y);

		double t_1_z = (bounds[0][2] - d.origin[2]) *invdz;
		double t_2_z = (bounds[1][2] - d.origin[2]) *invdz;
		double t_min_z = std::min(t_1_z, t_2_z);
		double t_max_z = std::max(t_1_z, t_2_z);

		double tlast = std::min(std::min(t_max_x, t_max_y), t_max_z);
		if (tlast < 0) return false;

		double tfirst = std::max(std::max(t_min_x, t_min_y), t_min_z);

		if (tfirst > 0) t = tfirst; else t = 0;

		if (tlast - tfirst >= -0.000001) return true;
		return false;
	}

	bool intersection_invd(const Ray& invd, const char signs[3], double &t) const {

		double t_max;
		t_max = (bounds[signs[0]][0] - invd.origin[0]) *invd.direction[0];
		if (t_max < 0) return false;
		t = (bounds[1-signs[0]][0] - invd.origin[0]) *invd.direction[0];

		double t_min_y, t_max_y;
		t_max_y = (bounds[signs[1]][1] - invd.origin[1]) *invd.direction[1];
		if (t_max_y < 0) return false;
		t_min_y = (bounds[1-signs[1]][1] - invd.origin[1]) *invd.direction[1];

		if (t_min_y > t_max  ||  t_max_y < t) return false;
		if (t_min_y > t) t = t_min_y;
		if (t_max_y < t_max) t_max = t_max_y;
		
		double t_min_z, t_max_z;
		t_max_z = (bounds[signs[2]][2] - invd.origin[2]) *invd.direction[2];
		if (t_max_z < 0) return false;
		t_min_z = (bounds[1-signs[2]][2] - invd.origin[2]) *invd.direction[2];

		if (t > t_max_z || t_min_z > t_max) return false;
		if (t_min_z > t) t = t_min_z;
		//if (t_max_z < t_max) tmax = t_max_z;
		if (t < 0) t = 0;


		return true;
	}

	bool intersection_invd_positive_x(const Ray& invd, const char signs[3], double &t) const {

		double t_max;
		t_max = (bounds[1][0] - invd.origin[0]);
		if (t_max < 0) return false;
		t_max *= invd.direction[0];
		t = (bounds[0][0] - invd.origin[0]) *invd.direction[0];

		double t_min_y, t_max_y;
		t_max_y = (bounds[signs[1]][1] - invd.origin[1]) *invd.direction[1];
		if (t_max_y < 0) return false;
		t_min_y = (bounds[1 - signs[1]][1] - invd.origin[1]) *invd.direction[1];

		if (t_min_y > t_max || t_max_y < t) return false;
		if (t_min_y > t) t = t_min_y;
		if (t_max_y < t_max) t_max = t_max_y;

		double t_min_z, t_max_z;
		t_max_z = (bounds[signs[2]][2] - invd.origin[2]) *invd.direction[2];
		if (t_max_z < 0) return false;
		t_min_z = (bounds[1 - signs[2]][2] - invd.origin[2]) *invd.direction[2];

		if (t > t_max_z || t_min_z > t_max) return false;
		if (t_min_z > t) t = t_min_z;
		//if (t_max_z < t_max) tmax = t_max_z;
		if (t < 0) t = 0;


		return true;
	}

	bool intersection_invd_negative_x(const Ray& invd, const char signs[3], double &t) const {

		double t_max;
		t_max = (bounds[0][0] - invd.origin[0]);
		if (t_max > 0) return false;
		t_max *= invd.direction[0];
		t = (bounds[1][0] - invd.origin[0]) *invd.direction[0];

		double t_min_y, t_max_y;
		t_max_y = (bounds[signs[1]][1] - invd.origin[1]) *invd.direction[1];
		if (t_max_y < 0) return false;
		t_min_y = (bounds[1 - signs[1]][1] - invd.origin[1]) *invd.direction[1];

		if (t_min_y > t_max || t_max_y < t) return false;
		if (t_min_y > t) t = t_min_y;
		if (t_max_y < t_max) t_max = t_max_y;

		double t_min_z, t_max_z;
		t_max_z = (bounds[signs[2]][2] - invd.origin[2]) *invd.direction[2];
		if (t_max_z < 0) return false;
		t_min_z = (bounds[1 - signs[2]][2] - invd.origin[2]) *invd.direction[2];

		if (t > t_max_z || t_min_z > t_max) return false;
		if (t_min_z > t) t = t_min_z;
		//if (t_max_z < t_max) tmax = t_max_z;
		if (t < 0) t = 0;


		return true;
	}


	bool intersection_invd(const Ray& invd, const char signs[3], double &t, double &t_far) const {

		t_far = (bounds[signs[0]][0] - invd.origin[0]) *invd.direction[0];
		if (t_far < 0) return false;
		t = (bounds[1 - signs[0]][0] - invd.origin[0]) *invd.direction[0];

		double t_min_y, t_max_y;
		t_max_y = (bounds[signs[1]][1] - invd.origin[1]) *invd.direction[1];
		if (t_max_y < 0) return false;
		t_min_y = (bounds[1 - signs[1]][1] - invd.origin[1]) *invd.direction[1];

		if (t_min_y > t_far || t_max_y < t) return false;
		if (t_min_y > t) t = t_min_y;
		if (t_max_y < t_far) t_far = t_max_y;

		double t_min_z, t_max_z;
		t_max_z = (bounds[signs[2]][2] - invd.origin[2]) *invd.direction[2];
		if (t_max_z < 0) return false;
		t_min_z = (bounds[1 - signs[2]][2] - invd.origin[2]) *invd.direction[2];

		if (t > t_max_z || t_min_z > t_far) return false;
		if (t_min_z > t) t = t_min_z;
		if (t_max_z < t_far) t_far = t_max_z;
		if (t < 0) t = 0;


		return true;
	}

	Vector bounds[2];

};

class Object {
public:
	Object() { scale = 1; flip_normals = false; refr_index = 1.3; };
	virtual bool intersection(const Ray& d, Vector& P, double &t, MaterialValues &mat, double cur_best_t, int &triangle_id) const = 0;
	virtual bool intersection_shadow(const Ray& d, double &t, double cur_best_t, double dist_light) const = 0;

	virtual Vector get_translation(double time) const {
		return max_translation;// *time;
	}
	virtual Vector get_rotation(double time) const {
		return max_rotation; // *time;
	}
	virtual Vector get_scale(double time) const {
		return Vector(scale, scale, scale); // *time;
	}
	virtual void build_matrix(double time) {

		Vector v2;
		for (int i = 0; i < 3; i++) {
			Vector v(0, 0, 0);
			v[i] = 1;
			v2 = rotate_dir(v, get_rotation(time));
			trans_matrix[0 * 4 + i] = v2[0] * scale;
			trans_matrix[1 * 4 + i] = v2[1] * scale;
			trans_matrix[2 * 4 + i] = v2[2] * scale;

			rot_matrix[0 * 3 + i] = v2[0];
			rot_matrix[1 * 3 + i] = v2[1];
			rot_matrix[2 * 3 + i] = v2[2];
						
			v2 = inverse_rotate_dir(v2, get_rotation(time));
			inv_trans_matrix[0 * 4 + i] = v2[0] / scale;
			inv_trans_matrix[1 * 4 + i] = v2[1] / scale;
			inv_trans_matrix[2 * 4 + i] = v2[2] / scale;
		}

		v2 = rotate_dir(-rotation_center, get_rotation(time));
		trans_matrix[0 * 4 + 3] = v2[0]*scale + rotation_center[0] + get_translation(time)[0];
		trans_matrix[1 * 4 + 3] = v2[1] * scale + rotation_center[1] + get_translation(time)[1];
		trans_matrix[2 * 4 + 3] = v2[2] * scale + rotation_center[2] + get_translation(time)[2];

		v2 = inverse_rotate_dir(-rotation_center-get_translation(time), get_rotation(time));
		inv_trans_matrix[0 * 4 + 3] = v2[0] /scale + rotation_center[0] ;
		inv_trans_matrix[1 * 4 + 3] = v2[1] / scale + rotation_center[1];
		inv_trans_matrix[2 * 4 + 3] = v2[2] / scale + rotation_center[2];		

	}

	Vector apply_transformation(const Vector& v) {
		Vector v2;
		v2[0] = trans_matrix[0] * v[0] + trans_matrix[1] * v[1] + trans_matrix[2] * v[2] + trans_matrix[3];
		v2[1] = trans_matrix[4] * v[0] + trans_matrix[5] * v[1] + trans_matrix[6] * v[2] + trans_matrix[7];
		v2[2] = trans_matrix[8] * v[0] + trans_matrix[9] * v[1] + trans_matrix[10] * v[2] + trans_matrix[11];
		return v2;
	}
	Vector apply_rotation_scaling(const Vector& v) {
		Vector v2;
		v2[0] = trans_matrix[0] * v[0] + trans_matrix[1] * v[1] + trans_matrix[2] * v[2];
		v2[1] = trans_matrix[4] * v[0] + trans_matrix[5] * v[1] + trans_matrix[6] * v[2];
		v2[2] = trans_matrix[8] * v[0] + trans_matrix[9] * v[1] + trans_matrix[10] * v[2];
		return v2;
	}
	Vector apply_rotation(const Vector& v) {
		Vector v2;
		v2[0] = rot_matrix[0] * v[0] + rot_matrix[1] * v[1] + rot_matrix[2] * v[2];
		v2[1] = rot_matrix[3] * v[0] + rot_matrix[4] * v[1] + rot_matrix[5] * v[2];
		v2[2] = rot_matrix[6] * v[0] + rot_matrix[7] * v[1] + rot_matrix[8] * v[2];
		return v2;
	}
	Vector apply_inverse_transformation(const Vector& v) {
		Vector v2;
		v2[0] = inv_trans_matrix[0] * v[0] + inv_trans_matrix[1] * v[1] + inv_trans_matrix[2] * v[2] + inv_trans_matrix[3];
		v2[1] = inv_trans_matrix[4] * v[0] + inv_trans_matrix[5] * v[1] + inv_trans_matrix[6] * v[2] + inv_trans_matrix[7];
		v2[2] = inv_trans_matrix[8] * v[0] + inv_trans_matrix[9] * v[1] + inv_trans_matrix[10] * v[2] + inv_trans_matrix[11];
		return v2;
	}
	Vector apply_inverse_rotation_scaling(const Vector& v) {
		Vector v2;
		v2[0] = inv_trans_matrix[0] * v[0] + inv_trans_matrix[1] * v[1] + inv_trans_matrix[2] * v[2];
		v2[1] = inv_trans_matrix[4] * v[0] + inv_trans_matrix[5] * v[1] + inv_trans_matrix[6] * v[2];
		v2[2] = inv_trans_matrix[8] * v[0] + inv_trans_matrix[9] * v[1] + inv_trans_matrix[10] * v[2];
		return v2;
	}

	MaterialValues queryMaterial(int idx, double u, double v) const {
		
		MaterialValues mat;
		u = Texture::wrap(u);
		v = Texture::wrap(v);
		
		if (idx >= textures.size()) {
			mat.Kd = Vector(1., 1., 1.);
		} else {
			mat.Kd = textures[idx].getVec(u, v);
		}
		if (idx >= specularmap.size()) {
			mat.Ks = Vector(0., 0., 0.);
		} else {
			mat.Ks = specularmap[idx].getVec(u, v);
		}
		if (idx >= roughnessmap.size()) {
			mat.Ne = Vector(1., 1., 1.);
		} else {
			mat.Ne = roughnessmap[idx].getVec(u, v);
		}
		mat.Ke = Vector(0., 0., 0.);
		return mat;
	}

	virtual void save_to_file(FILE* f) {
		fprintf(f, "name: %s\n", name.c_str());
		fprintf(f, "refr_index: %lf\n", refr_index);
		fprintf(f, "miroir: %u\n", miroir?1:0);
		fprintf(f, "transparent: %u\n", transparent ? 1 : 0);
		fprintf(f, "translation: (%lf, %lf, %lf)\n", max_translation[0], max_translation[1], max_translation[2]);
		fprintf(f, "rotation: (%lf, %lf, %lf)\n", max_rotation[0], max_rotation[1], max_rotation[2]);
		fprintf(f, "center: (%lf, %lf, %lf)\n", rotation_center[0], rotation_center[1], rotation_center[2]);
		fprintf(f, "scale: %lf\n", scale);
		fprintf(f, "display_edges: %u\n", display_edges ? 1 : 0);
		fprintf(f, "interp_normals: %u\n", interp_normals ? 1 : 0);
		fprintf(f, "flip_normals: %u\n", flip_normals ? 1 : 0);
		fprintf(f, "nb_textures: %u\n", textures.size());
		for (int i = 0; i < textures.size(); i++) {
			fprintf(f, "texture: %s\n", textures[i].filename.c_str());
			fprintf(f, "multiplier: (%lf, %lf, %lf)\n", textures[i].multiplier[0], textures[i].multiplier[1], textures[i].multiplier[2]);
		}
		fprintf(f, "nb_normalmaps: %u\n", normal_map.size());
		for (int i = 0; i < normal_map.size(); i++) {
			fprintf(f, "texture: %s\n", normal_map[i].filename.c_str());
			fprintf(f, "multiplier: (%lf, %lf, %lf)\n", normal_map[i].multiplier[0], normal_map[i].multiplier[1], normal_map[i].multiplier[2]);
		}
		fprintf(f, "nb_specularmaps: %u\n", specularmap.size());
		for (int i = 0; i < specularmap.size(); i++) {
			fprintf(f, "texture: %s\n", specularmap[i].filename.c_str());
			fprintf(f, "multiplier: (%lf, %lf, %lf)\n", specularmap[i].multiplier[0], specularmap[i].multiplier[1], specularmap[i].multiplier[2]);
		}
		fprintf(f, "nb_alphamaps: %u\n", alphamap.size());
		for (int i = 0; i < alphamap.size(); i++) {
			fprintf(f, "texture: %s\n", alphamap[i].filename.c_str());
			fprintf(f, "multiplier: (%lf, %lf, %lf)\n", alphamap[i].multiplier[0], alphamap[i].multiplier[1], alphamap[i].multiplier[2]);
		}
		fprintf(f, "nb_expmaps: %u\n", roughnessmap.size());
		for (int i = 0; i < roughnessmap.size(); i++) {
			fprintf(f, "texture: %s\n", roughnessmap[i].filename.c_str());
			fprintf(f, "multiplier: (%lf, %lf, %lf)\n", roughnessmap[i].multiplier[0], roughnessmap[i].multiplier[1], roughnessmap[i].multiplier[2]);
		}
	}

	virtual void load_from_file(FILE* f) {
		char line[512];
		int mybool;
		fscanf(f, "name: %[^\n]\n", line);
		name = std::string(line);
		fscanf(f, "refr_index: %lf\n", &refr_index);
		fscanf(f, "miroir: %u\n", &mybool); miroir = mybool;
		fscanf(f, "transparent: %u\n", &mybool); transparent = mybool;
		fscanf(f, "translation: (%lf, %lf, %lf)\n", &max_translation[0], &max_translation[1], &max_translation[2]);
		fscanf(f, "rotation: (%lf, %lf, %lf)\n", &max_rotation[0], &max_rotation[1], &max_rotation[2]);
		fscanf(f, "center: (%lf, %lf, %lf)\n", &rotation_center[0], &rotation_center[1], &rotation_center[2]);
		fscanf(f, "scale: %lf\n", &scale);
		fscanf(f, "display_edges: %u\n", &mybool); display_edges = mybool;
		fscanf(f, "interp_normals: %u\n", &mybool); interp_normals = mybool;
			fscanf(f, "flip_normals: %u\n", &mybool); flip_normals = mybool;

		int n; 
		fscanf(f, "nb_textures: %u\n", &n);
		for (int i = 0; i < n; i++) {
			fscanf(f, "texture: %[^\n]\n", line);
			if (line[0] == 'C' && line[1] == 'o') {  //Color
				Vector col;
				sscanf(line, "Color: (%lf, %lf, %lf)\n", &col[0], &col[1], &col[2]);				
				add_col_texture(col/255.);				
			} else {
				add_texture(line);
			}			
			fscanf(f, "multiplier: (%lf, %lf, %lf)\n", &textures[textures.size() - 1].multiplier[0], &textures[textures.size() - 1].multiplier[1], &textures[textures.size() - 1].multiplier[2]);
		}
		fscanf(f, "nb_normalmaps: %u\n", &n);
		for (int i = 0; i < n; i++) {
			fscanf(f, "texture: %[^\n]\n", line);
			if (line[0] == 'N' && line[1] == 'u') {  //Null
				add_null_normalmap();
			} else {
				add_normalmap(line);
			}
			fscanf(f, "multiplier: (%lf, %lf, %lf)\n", &normal_map[normal_map.size() - 1].multiplier[0], &normal_map[normal_map.size() - 1].multiplier[1], &normal_map[normal_map.size() - 1].multiplier[2]);
		}
		fscanf(f, "nb_specularmaps: %u\n", &n);
		for (int i = 0; i < n; i++) {
			fscanf(f, "texture: %[^\n]\n", line);
			if (line[0] == 'C' && line[1] == 'o') {  //Color
				Vector col;
				sscanf(line, "Color: (%lf, %lf, %lf)\n", &col[0], &col[1], &col[2]);
				add_col_specular(col/255.);
			} else {
				add_specularmap(line);
			}
			fscanf(f, "multiplier: (%lf, %lf, %lf)\n", &specularmap[specularmap.size() - 1].multiplier[0], &specularmap[specularmap.size() - 1].multiplier[1], &specularmap[specularmap.size() - 1].multiplier[2]);
		}
		fscanf(f, "nb_alphamaps: %u\n", &n);
		for (int i = 0; i < n; i++) {
			fscanf(f, "texture: %[^\n]\n", line);
			double col;
			if (sscanf(line, "%lf\n", &col)==1) {  //Color			
				add_col_alpha(col);
			} else {
				add_alphamap(line);
			}
			fscanf(f, "multiplier: (%lf, %lf, %lf)\n", &alphamap[alphamap.size() - 1].multiplier[0], &alphamap[alphamap.size() - 1].multiplier[1], &alphamap[alphamap.size() - 1].multiplier[2]);
		}

		fscanf(f, "nb_expmaps: %u\n", &n);
		for (int i = 0; i < n; i++) {
			fscanf(f, "texture: %[^\n]\n", line);
			if (line[0] == 'C' && line[1] == 'o') {  //Color
				Vector col;
				sscanf(line, "Color: (%lf, %lf, %lf)\n", &col[0], &col[1], &col[2]);
				add_col_roughness(col);
			} else {
				add_roughnessmap(line);
			}
			fscanf(f, "multiplier: (%lf, %lf, %lf)\n", &roughnessmap[roughnessmap.size() - 1].multiplier[0], &roughnessmap[roughnessmap.size() - 1].multiplier[1], &roughnessmap[roughnessmap.size() - 1].multiplier[2]);
		}
	}

	static Object* create_from_file(FILE* f);

	virtual void add_texture(const char* filename);
	virtual void set_texture(const char* filename, int idx);
	virtual void add_specularmap(const char* filename);
	virtual void set_specularmap(const char* filename, int idx);
	virtual void add_normalmap(const char* filename);
	virtual void set_normalmap(const char* filename, int idx);
	virtual void add_roughnessmap(const char* filename);
	virtual void set_roughnessmap(const char* filename, int idx);
	virtual void add_col_roughness(const Vector &col);
	virtual void set_col_roughness(const Vector &col, int idx);
	virtual void add_alphamap(const char* filename);
	virtual void set_alphamap(const char* filename, int idx);
	virtual void set_col_alpha(double col, int idx);
	virtual void remove_texture(int id);
	virtual void remove_normal(int id);
	virtual void remove_specular(int id);
	virtual void remove_alpha(int id);
	virtual void remove_roughness(int id);
	virtual void swap_alpha(int id1, int id2);
	virtual void swap_textures(int id1, int id2);
	virtual void swap_normal(int id1, int id2);
	virtual void swap_specular(int id1, int id2);
	virtual void swap_roughness(int id1, int id2);
	virtual void add_col_texture(const Vector& col);
	virtual void set_col_texture(const Vector& col, int idx);
	virtual void add_col_specular(const Vector& col);
	virtual void set_col_specular(const Vector& col, int idx);
	virtual void add_null_normalmap();
	virtual void set_null_normalmap(int idx);
	virtual void add_col_alpha(double col);

	std::string name;
	double refr_index;
	bool miroir;
	bool transparent;
	Vector max_translation;
	Vector max_rotation; // Euler angles
	Vector rotation_center;
	double scale;
	bool display_edges, interp_normals, flip_normals;

	/*std::vector<std::vector<unsigned char> > textures;
	std::vector<std::vector<unsigned char> > specularmap;
	std::vector<std::vector<unsigned char> > alphamap;
	std::vector<std::vector<unsigned char> > roughnessmap;
	std::vector<std::vector<double> > normal_map;

	std::vector<std::string> textures_filenames, normalmap_filenames, alphamap_filenames, specularmap_filenames, roughnessmap_filenames;
	std::vector<int> w, h, normalW, normalH, alphamapW, alphamapH, specularmapW, specularmapH, roughnessmapW, roughnessmapH;*/


	std::vector<Texture> textures, specularmap, alphamap, roughnessmap, normal_map;
	double trans_matrix[12], inv_trans_matrix[12], rot_matrix[9];
};




class Sphere : public Object {
public:
	Sphere() {};
	Sphere(const Vector &origin, double rayon, bool mirror = false, bool transp = false, bool normal_swapped = false, const char* envmap_file = NULL) {
		init(origin, rayon, mirror, transp, normal_swapped, envmap_file);
	};

	void init(const Vector &origin, double rayon, bool mirror = false, bool transp = false, bool normal_swapped = false, const char* envmap_file = NULL) {
		O = origin;
		R = rayon;
		miroir = mirror;
		transparent = transp;
		flip_normals = normal_swapped;
		has_envmap = (envmap_file != NULL);
		envmapfilename = "";
		if (has_envmap) {
			load_envmap(envmap_file);
		}
		this->rotation_center = origin;
		name = "Sphere";
		//bbox.bounds[0] = O - Vector(R, R, R);
		//bbox.bounds[1] = O + Vector(R, R, R);
	}

	void save_to_file(FILE* f) {

		fprintf(f, "NEW SPHERE\n");

		Object::save_to_file(f);

		fprintf(f, "is_envmap: %u\n", has_envmap ? 1 : 0);
		fprintf(f, "envmapfilename: %s\n", envmapfilename.c_str());
		fprintf(f, "O: (%lf, %lf, %lf)\n", O[0], O[1], O[2]);
		fprintf(f, "R: %f\n", R);
	}

	static Sphere* create_from_file(FILE* f) {
		Sphere* result = new Sphere();
		result->Object::load_from_file(f);

		char line[512];
		int has_envmap;
		fscanf(f, "is_envmap: %u\n", &has_envmap);

		std::string envmapfilename;
		if (has_envmap) {
			fscanf(f, "envmapfilename: %[^\n]\n", line);
			envmapfilename = std::string(line);
		} else {
			fscanf(f, "envmapfilename: \n");
		}

		Vector O;
		double R;
		fscanf(f, "O: (%lf, %lf, %lf)\n", &O[0], &O[1], &O[2]);
		fscanf(f, "R: %lf\n", &R);

		result->init(O, R, result->miroir, result->transparent, result->flip_normals, has_envmap ? envmapfilename.c_str() : NULL);
		return result;
	}

	void load_envmap(const char* filename) {
		load_bmp(filename, envtex, envW, envH);
		envmapfilename = std::string(filename);
		has_envmap = true;
	}

	bool intersection(const Ray& d, Vector& P, double &t, MaterialValues &mat, double cur_best_t, int &triangle_id) const {
		//if (!bbox.intersection(d)) return false;

		// resout a*t^2 + b*t + c = 0
		double a = d.direction.getNorm2();
		double b = 2 * dot(d.direction, d.origin - O);
		double c = (d.origin - O).getNorm2() - R*R;

		double delta = b*b - 4 * a*c;
		if (delta < 0) return false;
		
		double sqDelta = sqrt(delta);
		double inv2a = 1. / (2. * a);
		double t2 = (-b + sqDelta) * inv2a;

		if (t2 < 0) return false;

		double t1 = (-b - sqDelta) * inv2a;

		if (t1 > 0)
			t = t1;
		else
			t = t2;

		P = d.origin + t*d.direction;
		Vector N = (P - O).getNormalized();

		double theta = 1 - acos(N[1]) / M_PI;
		double phi = (atan2(-N[2], N[0]) + M_PI) / (2.*M_PI);
		mat = queryMaterial(0, theta, phi);

		mat.shadingN = N;
		if (has_envmap) {
			int idx = (int)(theta*(envH-1.))*envW + (int)(phi*(envW-1.));		
			if (idx<0 || idx>=envW*envH) mat.Ke = Vector(0., 0., 0.); else
			mat.Ke = Vector(envtex[idx*3+0], envtex[idx * 3 + 1], envtex[idx * 3 + 2])/255. *100000.;
		} else {
			mat.Ke = Vector(0.,0.,0.);
		}

		if (flip_normals) mat.shadingN = -mat.shadingN;
		triangle_id = -1;
		return true;
	}

	bool intersection_shadow(const Ray& d, double &t, double cur_best_t, double dist_light) const {
		double a = d.direction.getNorm2();
		double b = 2 * dot(d.direction, d.origin - O);
		double c = (d.origin - O).getNorm2() - R*R;

		double delta = b*b - 4 * a*c;
		if (delta < 0) return false;

		double sqDelta = sqrt(delta);
		double inv2a = 1. / (2. * a);
		double t2 = (-b + sqDelta) * inv2a;

		if (t2 < 0) return false;

		double t1 = (-b - sqDelta) * inv2a;

		if (t1 > 0)
			t = t1;
		else
			t = t2;
		return true;
	}

	Vector O;
	double R;
	bool has_envmap;
	std::vector<unsigned char> envtex;
	std::string envmapfilename;
	int envW, envH;
	//BBox bbox;
};


class Disk {
public:
	Disk(const Vector& center, const Vector& normal, double radius): c(center),n(normal), r(radius) {
	}
	bool intersection(const Ray& d, Vector& P, Vector &N, double &t) {
		N = n;
		t = dot(c - d.origin, N) / dot(d.direction, N);
		if (t < 0 || t != t) return false;  //isnan

		P = d.origin + t*d.direction;
		double r2 = (P - c).getNorm2();
		return (r2 <= r*r);
	}
	const Vector& c;
	const Vector& n;
	double r;
};




class Plane : public Object {
public:
	Plane() {};
	Plane(const Vector& A, const Vector &N, bool mirror = false, bool transp = false) {
		init(A, N, mirror, transp);
	};

	void init(const Vector& A, const Vector &N, bool mirror = false, bool transp = false) {
		this->A = A;
		this->vecN = N;
		miroir = mirror;
		transparent = transp;
		name = "Plane";
	}

	bool intersection(const Ray& d, Vector& P, double &t, MaterialValues &mat, double cur_best_t, int &triangle_id) const {
		mat.shadingN = vecN;
		double ddot = dot(d.direction, vecN);
		if (abs(ddot) < 1E-15) return false;
		t = dot(A - d.origin, vecN) / ddot;
		if (t <= 0.) return false;

		P = d.origin + t*d.direction;
		triangle_id = -1;

		double u = P[0]*0.1;
		double v = P[2]*0.1; 
		mat = queryMaterial(0, u, v);

		return true;
	}

	bool intersection_shadow(const Ray& d, double &t, double cur_best_t, double dist_light) const {
		double ddot = dot(d.direction, vecN);
		if (abs(ddot) < 1E-15) return false;
		t = dot(A - d.origin, vecN) / ddot;
		if (t <= 0.) return false;
		return true;
	}

	void save_to_file(FILE* f) {

		fprintf(f, "NEW PLANE\n");

		Object::save_to_file(f);

		fprintf(f, "Point: (%lf, %lf, %lf)\n", A[0], A[1], A[2]);
		fprintf(f, "N: (%lf, %lf, %lf)\n", vecN[0], vecN[1], vecN[2]);
	}

	static Plane* create_from_file(FILE* f) {
		Plane* result = new Plane();
		result->Object::load_from_file(f);

		char line[512];

		Vector Point, N;
		fscanf(f, "Point: (%lf, %lf, %lf)\n", &Point[0], &Point[1], &Point[2]);
		fscanf(f, "N: (%lf, %lf, %lf)\n", &N[0], &N[1], &N[2]);
		

		result->init(Point, N, result->transparent);
		return result;
	}

	Vector A, vecN;
};




class Scene {
public:
	Scene() {};

	void addObject(Object* o) { objects.push_back(o); }

	void clear() {
		for (int i = 0; i < objects.size(); i++) {
			delete objects[i];
		}
		objects.clear();
		lumiere = NULL;

	}
	void deleteObject(int id) {
		if (id < 0 || id >= objects.size()) return;
		delete objects[id];
		objects.erase(objects.begin() + id);
	}


	bool intersection(const Ray& d, Vector& P, int &sphere_id, double &min_t, MaterialValues &mat, int &triangle_id) const {

		bool has_inter = false;
		min_t = 1E99;

		for (int i = 0; i < objects.size(); i++) {
			Vector localP;
			MaterialValues localmat;
			double t;

			Vector transformed_dir = objects[i]->apply_inverse_rotation_scaling(d.direction);
			Vector new_origin = objects[i]->apply_inverse_transformation(d.origin); 
			Ray transformed_ray(new_origin, transformed_dir, d.time);

			bool local_has_inter = objects[i]->intersection(transformed_ray, localP, t, localmat, min_t, triangle_id);

			if (local_has_inter) {		
				if (t < min_t) {
					has_inter = true;
					min_t = t;
					
					P = objects[i]->apply_transformation(localP);
					sphere_id = i;
					mat = localmat;
					mat.shadingN = objects[i]->apply_rotation(localmat.shadingN); 
				}
			}
		}

		return has_inter;
	}


	bool intersection_shadow(const Ray& d, double &min_t, double dist_light) const {

		bool has_inter = false;
		min_t = 1E99;

		for (int i = 0; i < objects.size(); i++) {
			Vector transformed_dir = objects[i]->apply_inverse_rotation_scaling(d.direction);
			Vector new_origin = objects[i]->apply_inverse_transformation(d.origin); 
			Ray transformed_ray(new_origin, transformed_dir, d.time);

			double t;

			bool local_has_inter = objects[i]->intersection_shadow(transformed_ray, t, min_t, dist_light);

			if (local_has_inter) {
				if (t < dist_light*0.999) {
					return true;
				}
			}
		}

		return false;
	}

	std::vector<Object*> objects;
	Sphere *lumiere;
	double intensite_lumiere;
	double envmap_intensity;
	double fog_density;
	int fog_type; // 0 : uniform, 1: exponential
};