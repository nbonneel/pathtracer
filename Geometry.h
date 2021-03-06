#pragma once

#include <iostream>
#include <vector>
#include "Vector.h"
#include <algorithm>
#include <map>
#include "utils.h"
#include "BRDF.h"

#define cimg_display 0
#include "CImg.h"
#include "utils.h"

#include <omp.h>

#ifdef USE_EMBREE
#include "embree3/rtcore.h"

#endif

#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

enum ObjectType {OT_TRIMESH, OT_SPHERE, OT_PLANE, OT_POINTSET, OT_CYLINDER, OT_YARNS, OT_FLUID};

void load_bmp(const char* filename, std::vector<unsigned char> &tex, int &W, int& H);

class Scene;

class PointSet;




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
		t_max = (bounds[(int)signs[0]][0] - invd.origin[0]) *invd.direction[0];
		if (t_max < 0) return false;
		t = (bounds[1-signs[0]][0] - invd.origin[0]) *invd.direction[0];

		double t_min_y, t_max_y;
		t_max_y = (bounds[(int)signs[1]][1] - invd.origin[1]) *invd.direction[1];
		if (t_max_y < 0) return false;
		t_min_y = (bounds[1-signs[1]][1] - invd.origin[1]) *invd.direction[1];

		if (t_min_y > t_max  ||  t_max_y < t) return false;
		if (t_min_y > t) t = t_min_y;
		if (t_max_y < t_max) t_max = t_max_y;
		
		double t_min_z, t_max_z;
		t_max_z = (bounds[(int)signs[2]][2] - invd.origin[2]) *invd.direction[2];
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
		t_max_y = (bounds[(int)signs[1]][1] - invd.origin[1]) *invd.direction[1];
		if (t_max_y < 0) return false;
		t_min_y = (bounds[1 - signs[1]][1] - invd.origin[1]) *invd.direction[1];

		if (t_min_y > t_max || t_max_y < t) return false;
		if (t_min_y > t) t = t_min_y;
		if (t_max_y < t_max) t_max = t_max_y;

		double t_min_z, t_max_z;
		t_max_z = (bounds[(int)signs[2]][2] - invd.origin[2]) *invd.direction[2];
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
		t_max_y = (bounds[(int)signs[1]][1] - invd.origin[1]) *invd.direction[1];
		if (t_max_y < 0) return false;
		t_min_y = (bounds[1 - signs[1]][1] - invd.origin[1]) *invd.direction[1];

		if (t_min_y > t_max || t_max_y < t) return false;
		if (t_min_y > t) t = t_min_y;
		if (t_max_y < t_max) t_max = t_max_y;

		double t_min_z, t_max_z;
		t_max_z = (bounds[(int)signs[2]][2] - invd.origin[2]) *invd.direction[2];
		if (t_max_z < 0) return false;
		t_min_z = (bounds[1 - signs[2]][2] - invd.origin[2]) *invd.direction[2];

		if (t > t_max_z || t_min_z > t_max) return false;
		if (t_min_z > t) t = t_min_z;
		//if (t_max_z < t_max) tmax = t_max_z;
		if (t < 0) t = 0;


		return true;
	}


	bool intersection_invd(const Ray& invd, const char signs[3], double &t, double &t_far) const {

		t_far = (bounds[(int)signs[0]][0] - invd.origin[0]) *invd.direction[0];
		if (t_far < 0) return false;
		t = (bounds[1 - signs[0]][0] - invd.origin[0]) *invd.direction[0];

		double t_min_y, t_max_y;
		t_max_y = (bounds[(int)signs[1]][1] - invd.origin[1]) *invd.direction[1];
		if (t_max_y < 0) return false;
		t_min_y = (bounds[1 - signs[1]][1] - invd.origin[1]) *invd.direction[1];

		if (t_min_y > t_far || t_max_y < t) return false;
		if (t_min_y > t) t = t_min_y;
		if (t_max_y < t_far) t_far = t_max_y;

		double t_min_z, t_max_z;
		t_max_z = (bounds[(int)signs[2]][2] - invd.origin[2]) *invd.direction[2];
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
	virtual ~Object() {};
	Object() { 
		scale = 1; 
		flip_normals = false; 
		miroir = false; 
		ghost = false; 
		scene = NULL;
		brdf = new PhongBRDF();
	};
	virtual bool intersection(const Ray& d, Vector& P, double &t, MaterialValues &mat, double cur_best_t, int &triangle_id) const = 0;
	virtual bool intersection_shadow(const Ray& d, double &t, double cur_best_t, double dist_light) const = 0;

	// uniformly random choice of intersection (see reservoir sampling) between min_t and max_t (may return no intersection if current_nb_intersections>0 even though there are intersections ; they will just not be randomly chosen).
	// Should NOT change any value if there is no intersection. Should ONLY change current_nb_intersections if there are intersections in the interval but they were not chosen due to randomness.
	virtual bool reservoir_sampling_intersection(const Ray& d, Vector& P, double &t, MaterialValues &mat, int &triangle_id, int &current_nb_intersections, double min_t, double max_t) const = 0;

	virtual Vector get_translation(double frame, bool is_recording) const {

		if (is_recording) return max_translation;

		std::map<double, Vector>::const_iterator it1 = translation_keyframes.upper_bound(frame); 
		std::map<double, Vector>::const_iterator it2 = it1;
		if (it1 == translation_keyframes.end()) {
			if (translation_keyframes.size() != 0) {
				return translation_keyframes.rbegin()->second;
			} else
				return max_translation;
		}
		if (it1 == translation_keyframes.begin()) { 
			return it1->second;
		}
		it1--;		
		double t = (frame - it1->first) / (it2->first - it1->first);
		return (1 - t)*it1->second + t*it2->second;		
	}
	virtual Matrix<3,3> get_rotation(double frame, bool is_recording) const {
		if (is_recording) return mat_rotation;

		std::map<double, Matrix<3, 3> >::const_iterator it1 = rotation_keyframes.upper_bound(frame);
		std::map<double, Matrix<3, 3> >::const_iterator it2 = it1;
		if (it1 == rotation_keyframes.end()) {
			if (rotation_keyframes.size() != 0) {
				return rotation_keyframes.rbegin()->second;
			} else
				return mat_rotation;
		}
		if (it1 == rotation_keyframes.begin()) {
			return it1->second;
		}
		it1--;
		double t = (frame - it1->first) / (it2->first - it1->first);
		return Slerp(it1->second, it2->second, t);
	}
	virtual double get_scale(double frame, bool is_recording) const {
		if (is_recording) return scale;

		std::map<double, double >::const_iterator it1 = scale_keyframes.upper_bound(frame);
		std::map<double, double >::const_iterator it2 = it1;
		if (it1 == scale_keyframes.end()) {
			if (scale_keyframes.size() != 0) {
				return scale_keyframes.rbegin()->second;
			} else
				return scale;
		}
		if (it1 == scale_keyframes.begin()) {
			return it1->second;
		}
		it1--;
		double t = (frame - it1->first) / (it2->first - it1->first);
		return (1 - t)*it1->second + t*it2->second;
	}
	virtual void add_keyframe(int frame) {
		rotation_keyframes[frame] = mat_rotation;
		translation_keyframes[frame] = max_translation;
		scale_keyframes[frame] = scale;
	}
	std::map<double, Matrix<3, 3> > rotation_keyframes; // should do one Transform class.  Frames are double in case I do motion blur later.
	std::map<double, Vector > translation_keyframes;
	std::map<double, double > scale_keyframes;

	virtual void build_matrix(double frame, bool is_recording) {

		Matrix<3, 3> m = get_rotation(frame, is_recording);
		Matrix<3, 3> mt = m.getTransposed();
		Vector v2;
		double s = get_scale(frame, is_recording);
		Vector tr = get_translation(frame, is_recording);
		for (int i = 0; i < 3; i++) {
			Vector v(0, 0, 0);
			v[i] = 1;
			v2 = Vector(m[0*3+i], m[1 * 3 + i], m[2 * 3 + i]);// rotate_dir(v, get_rotation(time));
			trans_matrix[0 * 4 + i] = v2[0] * s;
			trans_matrix[1 * 4 + i] = v2[1] * s;
			trans_matrix[2 * 4 + i] = v2[2] * s;

			rot_matrix[0 * 3 + i] = v2[0];
			rot_matrix[1 * 3 + i] = v2[1];
			rot_matrix[2 * 3 + i] = v2[2];
						
			//v2 = inverse_rotate_dir(v2, get_rotation(time));
			v2 = Vector(mt[0 * 3 + i], mt[1 * 3 + i], mt[2 * 3 + i]);
			inv_trans_matrix[0 * 4 + i] = v2[0] / s;
			inv_trans_matrix[1 * 4 + i] = v2[1] / s;
			inv_trans_matrix[2 * 4 + i] = v2[2] / s;
		}

		//v2 = rotate_dir(-rotation_center, get_rotation(time));
		v2 = m*(-rotation_center);
		trans_matrix[0 * 4 + 3] = v2[0] * s + rotation_center[0] + tr[0];
		trans_matrix[1 * 4 + 3] = v2[1] * s + rotation_center[1] + tr[1];
		trans_matrix[2 * 4 + 3] = v2[2] * s + rotation_center[2] + tr[2];

		//v2 = inverse_rotate_dir(-rotation_center-get_translation(time), get_rotation(time));
		v2 = mt*(-rotation_center - tr);
		inv_trans_matrix[0 * 4 + 3] = v2[0] / s + rotation_center[0] ;
		inv_trans_matrix[1 * 4 + 3] = v2[1] / s + rotation_center[1];
		inv_trans_matrix[2 * 4 + 3] = v2[2] / s + rotation_center[2];		

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
		if (idx >= subsurface.size()) {
			mat.Ksub = Vector(0., 0., 0.);
		} else {
			mat.Ksub = subsurface[idx].getVec(u, v);
		}
		if (idx >= roughnessmap.size()) {
			mat.Ne = Vector(1., 1., 1.);
		} else {
			mat.Ne = roughnessmap[idx].getVec(u, v);
		}
		if (idx >= transparent_map.size()) {
			mat.transp = false;
		} else {
			mat.transp = transparent_map[idx].getBool(u, v);
		}
		if (idx >= refr_index_map.size()) {
			mat.refr_index = 1.3;
		} else {
			mat.refr_index = refr_index_map[idx].getValRed(u, v);
		}
		mat.Ke = Vector(0., 0., 0.);
		return mat;
	}

	virtual void colorAnisotropy() {
		// see tri meshes
	}

	virtual void randomColors() {
		// see tri meshes
	}

	virtual void save_to_file(FILE* f) {
		fprintf(f, "name: %s\n", name.c_str());		
		fprintf(f, "miroir: %u\n", miroir?1:0);		
		fprintf(f, "ghost: %u\n", ghost ? 1 : 0);
		fprintf(f, "translation: (%lf, %lf, %lf)\n", max_translation[0], max_translation[1], max_translation[2]);
		fprintf(f, "rotation: (%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf)\n", mat_rotation[0], mat_rotation[1], mat_rotation[2], mat_rotation[3], mat_rotation[4], mat_rotation[5], mat_rotation[6], mat_rotation[7], mat_rotation[8]);
		fprintf(f, "center: (%lf, %lf, %lf)\n", rotation_center[0], rotation_center[1], rotation_center[2]);
		fprintf(f, "scale: %lf\n", scale); // TODO SAVE ANIMATIONS
		fprintf(f, "display_edges: %u\n", display_edges ? 1 : 0);
		fprintf(f, "interp_normals: %u\n", interp_normals ? 1 : 0);
		fprintf(f, "flip_normals: %u\n", flip_normals ? 1 : 0);
		fprintf(f, "nb_transforms: %u\n", static_cast<unsigned int>(translation_keyframes.size()));
		for (std::map<double, double>::iterator it = scale_keyframes.begin(); it != scale_keyframes.end(); ++it) {
			fprintf(f, "%lf %lf\n", it->first, it->second);
		}
		for (std::map<double, Vector>::iterator it = translation_keyframes.begin(); it != translation_keyframes.end(); ++it) {
			fprintf(f, "%lf %lf, %lf, %lf\n", it->first, it->second[0], it->second[1], it->second[2]);
		}
		for (std::map<double, Matrix<3, 3> >::iterator it = rotation_keyframes.begin(); it != rotation_keyframes.end(); ++it) {
			fprintf(f, "%lf %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", it->first, it->second[0], it->second[1], it->second[2], it->second[3], it->second[4], it->second[5], it->second[6], it->second[7], it->second[8]);
		}
		fprintf(f, "nb_textures: %u\n", static_cast<unsigned int>(textures.size()));
		for (int i = 0; i < textures.size(); i++) {
			fprintf(f, "texture: %s\n", textures[i].filename.c_str());
			fprintf(f, "multiplier: (%lf, %lf, %lf)\n", textures[i].multiplier[0], textures[i].multiplier[1], textures[i].multiplier[2]);
		}
		fprintf(f, "nb_normalmaps: %u\n", static_cast<unsigned int>(normal_map.size()));
		for (int i = 0; i < normal_map.size(); i++) {
			fprintf(f, "texture: %s\n", normal_map[i].filename.c_str());
			fprintf(f, "multiplier: (%lf, %lf, %lf)\n", normal_map[i].multiplier[0], normal_map[i].multiplier[1], normal_map[i].multiplier[2]);
		}
		fprintf(f, "nb_subsurfaces: %u\n", static_cast<unsigned int>(subsurface.size()));
		for (int i = 0; i < subsurface.size(); i++) {
			fprintf(f, "texture: %s\n", subsurface[i].filename.c_str());
			fprintf(f, "multiplier: (%lf, %lf, %lf)\n", subsurface[i].multiplier[0], subsurface[i].multiplier[1], subsurface[i].multiplier[2]);
		}
		fprintf(f, "nb_specularmaps: %u\n", static_cast<unsigned int>(specularmap.size()));
		for (int i = 0; i < specularmap.size(); i++) {
			fprintf(f, "texture: %s\n", specularmap[i].filename.c_str());
			fprintf(f, "multiplier: (%lf, %lf, %lf)\n", specularmap[i].multiplier[0], specularmap[i].multiplier[1], specularmap[i].multiplier[2]);
		}
		fprintf(f, "nb_alphamaps: %u\n", static_cast<unsigned int>(alphamap.size()));
		for (int i = 0; i < alphamap.size(); i++) {
			fprintf(f, "texture: %s\n", alphamap[i].filename.c_str());
			fprintf(f, "multiplier: (%lf, %lf, %lf)\n", alphamap[i].multiplier[0], alphamap[i].multiplier[1], alphamap[i].multiplier[2]);
		}
		fprintf(f, "nb_expmaps: %u\n", static_cast<unsigned int>(roughnessmap.size()));
		for (int i = 0; i < roughnessmap.size(); i++) {
			fprintf(f, "texture: %s\n", roughnessmap[i].filename.c_str());
			fprintf(f, "multiplier: (%lf, %lf, %lf)\n", roughnessmap[i].multiplier[0], roughnessmap[i].multiplier[1], roughnessmap[i].multiplier[2]);
		}
		fprintf(f, "nb_transpmaps: %u\n", static_cast<unsigned int>(transparent_map.size()));
		for (int i = 0; i < transparent_map.size(); i++) {
			fprintf(f, "texture: %s\n", transparent_map[i].filename.c_str());
			fprintf(f, "multiplier: %lf)\n", transparent_map[i].multiplier[0]);
		}
		fprintf(f, "nb_refrindexmaps: %u\n", static_cast<unsigned int>(refr_index_map.size()));
		for (int i = 0; i < refr_index_map.size(); i++) {
			fprintf(f, "texture: %s\n", refr_index_map[i].filename.c_str());
			fprintf(f, "multiplier: %lf)\n", refr_index_map[i].multiplier[0]);
		}
	}

	virtual void load_from_file(FILE* f, const char* replacedNames = NULL) {
		char line[512];
		int mybool;
		fscanf(f, "name: %[^\n]\n", line);
		name = std::string(line);

		if (replacedNames) {
			name.replace(name.find("#"), std::string("#").length(), std::string(replacedNames));
		}

		fscanf(f, "miroir: %u\n", &mybool); miroir = mybool;
		fscanf(f, "%[^\n]\n", line);
		if (line[0] == 'g') {
			sscanf(line, "ghost: %u\n", &mybool); ghost = mybool;
			fscanf(f, "translation: (%lf, %lf, %lf)\n", &max_translation[0], &max_translation[1], &max_translation[2]);
		} else {
			sscanf(line, "translation: (%lf, %lf, %lf)\n", &max_translation[0], &max_translation[1], &max_translation[2]);
		}
		fscanf(f, "rotation: (%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf)\n", &mat_rotation[0], &mat_rotation[1], &mat_rotation[2], &mat_rotation[3], &mat_rotation[4], &mat_rotation[5], &mat_rotation[6], &mat_rotation[7], &mat_rotation[8]);
		fscanf(f, "center: (%lf, %lf, %lf)\n", &rotation_center[0], &rotation_center[1], &rotation_center[2]);
		fscanf(f, "scale: %lf\n", &scale);
		fscanf(f, "display_edges: %u\n", &mybool); display_edges = mybool;
		fscanf(f, "interp_normals: %u\n", &mybool); interp_normals = mybool;
		fscanf(f, "flip_normals: %u\n", &mybool); flip_normals = mybool;

		fscanf(f, " %[^\n]\n", line);
		scale_keyframes.clear();
		translation_keyframes.clear();
		rotation_keyframes.clear();
		int n;
		if (line[4] == 'r') {
			sscanf(line, "nb_transforms: %u\n", &n);
			for (int i = 0; i < n; i++) {
				double s;
				double frame;
				fscanf(f, "%lf %lf\n", &frame, &s);
				scale_keyframes[frame] = s;
			}
			for (int i = 0; i < n; i++) {
				double frame;
				Vector t;
				fscanf(f, "%lf %lf, %lf, %lf\n", &frame, &t[0], &t[1], &t[2]);
				translation_keyframes[frame] = t;
			}
			for (int i = 0; i < n; i++) {
				double frame;
				Matrix<3,3> m;
				fscanf(f, "%lf %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", &frame, &m[0], &m[1], &m[2], &m[3], &m[4], &m[5], &m[6], &m[7], &m[8]);
				rotation_keyframes[frame] = m;
			}
			fscanf(f, "nb_textures: %u\n", &n);
		} else {
			sscanf(line, "nb_textures: %u\n", &n);
		}
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

		fscanf(f, "nb_subsurfaces: %u\n", &n);
		for (int i = 0; i < n; i++) {
			fscanf(f, "texture: %[^\n]\n", line);
			if (line[0] == 'C' && line[1] == 'o') {  //Color
				Vector col;
				sscanf(line, "Color: (%lf, %lf, %lf)\n", &col[0], &col[1], &col[2]);
				add_col_subsurface(col / 255.);
			} else {
				add_subsurface(line);
			}
			fscanf(f, "multiplier: (%lf, %lf, %lf)\n", &subsurface[subsurface.size() - 1].multiplier[0], &subsurface[subsurface.size() - 1].multiplier[1], &subsurface[subsurface.size() - 1].multiplier[2]);
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
		
		fscanf(f, "nb_transpmaps: %u\n", &n);
		for (int i = 0; i < n; i++) {
			fscanf(f, "texture: %[^\n]\n", line);
			add_transp_map(line);			
			fscanf(f, "multiplier: %lf)\n", &transparent_map[i].multiplier[0]);
		}
		fscanf(f, "nb_refrindexmaps: %u\n", &n);
		for (int i = 0; i < n; i++) {
			fscanf(f, "texture: %[^\n]\n", line);
			add_refr_map(line);
			fscanf(f, "multiplier: %lf)\n", &refr_index_map[i].multiplier[0]);
		}
	}

	static Object* create_from_file(FILE* f, Scene* s, const char* replacedNames = NULL);

	virtual void add_texture(const char* filename);
	virtual void set_texture(const char* filename, int idx);
	virtual void add_specularmap(const char* filename);
	virtual void set_specularmap(const char* filename, int idx);
	virtual void add_normalmap(const char* filename);
	virtual void set_normalmap(const char* filename, int idx);
	virtual void add_roughnessmap(const char* filename);	
	virtual void set_roughnessmap(const char* filename, int idx);
	virtual void add_transp_map(const char* filename);
	virtual void set_transp_map(const char* filename, int idx);
	virtual void add_subsurface(const char* filename);
	virtual void set_subsurface(const char* filename, int idx);
	virtual void set_col_transp(double col, int idx);
	virtual void add_col_transp(double col);
	virtual void add_refr_map(const char* filename);
	virtual void set_refr_map(const char* filename, int idx);
	virtual void set_col_refr(double col, int idx);
	virtual void add_col_refr(double col);
	virtual void add_col_roughness(const Vector &col);
	virtual void set_col_roughness(const Vector &col, int idx);
	virtual void add_alphamap(const char* filename);
	virtual void set_alphamap(const char* filename, int idx);
	virtual void set_col_alpha(double col, int idx);
	virtual void remove_refr(int id);
	virtual void remove_transp(int id);
	virtual void remove_texture(int id);
	virtual void remove_normal(int id);
	virtual void remove_specular(int id);
	virtual void remove_alpha(int id);
	virtual void remove_roughness(int id);
	virtual void remove_subsurface(int id);
	virtual void swap_refr(int id1, int id2);
	virtual void swap_transp(int id1, int id2);
	virtual void swap_alpha(int id1, int id2);
	virtual void swap_textures(int id1, int id2);
	virtual void swap_normal(int id1, int id2);
	virtual void swap_subsurface(int id1, int id2);
	virtual void swap_specular(int id1, int id2);
	virtual void swap_roughness(int id1, int id2);
	virtual void add_col_texture(const Vector& col);
	virtual void set_col_texture(const Vector& col, int idx);
	virtual void add_col_specular(const Vector& col);
	virtual void set_col_specular(const Vector& col, int idx);
	virtual void add_col_subsurface(const Vector& col);
	virtual void set_col_subsurface(const Vector& col, int idx);
	virtual void add_null_normalmap();
	virtual void set_null_normalmap(int idx);
	virtual void add_col_alpha(double col);

	std::string name;
	bool miroir;	
	Vector max_translation;
	Matrix<3, 3> mat_rotation; 
	Vector rotation_center;
	double scale;
	bool display_edges, interp_normals, flip_normals, ghost;
	ObjectType type;

	std::vector<Texture> textures, specularmap, alphamap, roughnessmap, normal_map, subsurface, transparent_map, refr_index_map;
	double trans_matrix[12], inv_trans_matrix[12], rot_matrix[9];
	BRDF* brdf;
	Scene* scene;
};


class Cylinder : public Object {
public:
	Cylinder() { type = OT_CYLINDER; };
	Cylinder(const Vector& A, const Vector&B, double R):A(A),B(B), R(R) {
		d = (B - A).getNormalized();
		len = sqrt((B - A).getNorm2());
		type = OT_CYLINDER;
	}

	bool intersection(const Ray& r, Vector& P, double &t, MaterialValues &mat, double cur_best_t, int &triangle_id) const {

		// norm2(X*t+Y) = R^2
		Vector X = r.direction - dot(r.direction, d)*d;
		Vector Y = r.origin - A - dot(r.origin - A, d)*d;
		double a = X.getNorm2();
		double b = 2 * dot(X, Y);
		double c = Y.getNorm2() - R*R;
		double delta = b*b - 4 * a*c;
		if (delta < 0) return false;
		double sdelta = sqrt(delta);
		double t2 = (-b + sdelta) / (2 * a);
		if (t2 < 0) return false;
		double t1 = (-b - sdelta) / (2 * a);
		if (t1 > 0) t = t1; else t = t2;
		P = r.origin + t*r.direction;
		double dP = dot(P - A, d);
		if (dP < 0 || dP > len) return false;
		
		Vector proj = A + dP*d;
		mat = queryMaterial(0, dP/len, 0.5);
		mat.shadingN = (P - proj); // .getNormalized();
		if (flip_normals) mat.shadingN = -mat.shadingN;
		triangle_id = -1;
		return true;
	}

	virtual bool reservoir_sampling_intersection(const Ray& r, Vector& P, double &t, MaterialValues &mat, int &triangle_id, int &current_nb_intersections, double min_t, double max_t) const {
		// norm2(X*t+Y) = R^2
		Vector X = r.direction - dot(r.direction, d)*d;
		Vector Y = r.origin - A - dot(r.origin - A, d)*d;
		double a = X.getNorm2();
		double b = 2 * dot(X, Y);
		double c = Y.getNorm2() - R * R;
		double delta = b * b - 4 * a*c;
		if (delta < 0) return false;
		double sdelta = sqrt(delta);
		double t2 = (-b + sdelta) / (2 * a);
		if (t2 < min_t) return false;
		double t1 = (-b - sdelta) / (2 * a);
		int threadid = omp_get_thread_num();
		const float invmax = 1.f / engine[threadid].max();
		bool has_inter = false;
		double dP;

		if (t1 >= min_t) { // 2 intersections with /infinite/ cylinder

			Vector trialP = r.origin + t1 * r.direction;
			dP = dot(P - A, d);
			if (dP >= 0 && dP <= len && t1 < max_t) {
				current_nb_intersections++;
				float r1 = engine[threadid]()*invmax;
				if (r1 < 1. / current_nb_intersections) {
					t = t1;
					P = trialP;
					has_inter = true;
				}
			}

			trialP = r.origin + t2 * r.direction;
			double dP2 = dot(P - A, d);
			if (dP2 >= 0 && dP2 <= len && t2 < max_t) {
				current_nb_intersections++;
				float r1 = engine[threadid]()*invmax;
				if (r1 < 1. / current_nb_intersections) {
					t = t2;
					P = trialP;
					has_inter = true;
					dP = dP2;
				}
			}
						  
		} else { // 1 intersection
			Vector trialP = r.origin + t2 * r.direction;
			dP = dot(P - A, d);
			if (dP >= 0 && dP <= len && t2 < max_t) {
				current_nb_intersections++;
				float r1 = engine[threadid]()*invmax;
				if (r1 < 1. / current_nb_intersections) {
					t = t2;
					has_inter = true;
					P = trialP;
				}
				
			}
		}
		if (!has_inter) return false;


		Vector proj = A + dP * d;
		mat = queryMaterial(0, dP / len, 0.5);
		mat.shadingN = (P - proj); // .getNormalized();
		if (flip_normals) mat.shadingN = -mat.shadingN;
		triangle_id = -1;
		return true;
	}

	bool intersection_shadow(const Ray& r, double &t, double cur_best_t, double dist_light) const {
		Vector P;
		MaterialValues mat;
		int tri;
		return intersection(r, P, t, mat, cur_best_t, tri);
	}

	Vector A, B, d;
	double R, len;
};


class Sphere : public Object {
public:
	Sphere() { type = OT_SPHERE; };
	Sphere(const Vector &origin, double rayon, bool mirror = false, bool normal_swapped = false, const char* envmap_file = NULL) {
		init(origin, rayon, mirror, normal_swapped, envmap_file);
	};

	void init(const Vector &origin, double rayon, bool mirror = false, bool normal_swapped = false, const char* envmap_file = NULL) {
		type = OT_SPHERE;
		O = origin;
		R = rayon;
		miroir = mirror;
		flip_normals = normal_swapped;
		has_envmap = (envmap_file != NULL);
		envmapfilename = "";
		if (has_envmap) {
			load_envmap(envmap_file);
			flip_normals = true;
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

		result->init(O, R, result->miroir, result->flip_normals, has_envmap ? envmapfilename.c_str() : NULL);
		return result;
	}

	void load_envmap(const char* filename) {
		load_bmp(filename, envtex, envW, envH);
		envmapfilename = std::string(filename);
		has_envmap = true;
	}

	bool intersection(const Ray& d, Vector& P, double &t, MaterialValues &mat, double cur_best_t, int &triangle_id) const {
		/*if (has_envmap) {
			P = d.origin + R * d.direction;
			Vector N = d.direction;
			double theta = 1 - acos(N[1]) / M_PI;
			double phi = (atan2(-N[2], N[0]) + M_PI) / (2.*M_PI);
			mat = queryMaterial(0, theta, phi);
			mat.shadingN = N;
			int idx = (int)(theta*(envH - 1.))*envW + (int)(phi*(envW - 1.));
			if (idx < 0 || idx >= envW * envH) mat.Ke = Vector(0., 0., 0.); else
				mat.Ke = Vector(envtex[idx * 3 + 0], envtex[idx * 3 + 1], envtex[idx * 3 + 2]) / 255. *100000.;

			if (flip_normals) mat.shadingN = -mat.shadingN;
			triangle_id = -1;
			return true;
		}*/
		
		//if (!bbox.intersection(d)) return false;

		// resout a*t^2 + b*t + c = 0
		double a =  d.direction.getNorm2(); // can't suppose unit norm as objects are scaled by scaling ray length
		double b =  dot(d.direction, d.origin - O);
		double c = (d.origin - O).getNorm2() - R*R;

		double delta = b*b - a*c;
		if (delta < 0) return false;
		
		double sqDelta = sqrt(delta);
		double inva = 1. / a;		
		double t2 = (-b + sqDelta)*inva;

		if (t2 < 0) return false;
		
		double t1 = (-b - sqDelta)*inva;

		if (t1 > 0)
			t = t1;
		else
			t = t2;

		P = d.origin + t*d.direction;
		Vector N = (P - O);

		if (has_envmap) {
			//P = d.origin + R * d.direction;
			//Vector N = d.direction;
			N.fast_normalize();
			double theta = 1 - acos(N[1]) / M_PI;
			double phi = (atan2(-N[2], N[0]) + M_PI) / (2.*M_PI);
			mat = queryMaterial(0, theta, phi);
			mat.shadingN = -N;
			int idx = 3*((int)(theta*(envH - 1.))*envW + (int)(phi*(envW - 1.)));
			if (idx < 0 || idx >= 3*envW * envH) mat.Ke = Vector(0., 0., 0.); else
				mat.Ke = Vector(envtex[idx + 0], envtex[idx + 1], envtex[idx + 2]) / 255. *100000.;
			//if (flip_normals) mat.shadingN = -mat.shadingN;
			triangle_id = -1;
			return true;
		}

		if (textures.size() != 0 || specularmap.size() != 0 || roughnessmap.size() != 0 || transparent_map.size() != 0 || refr_index_map.size() != 0) {
			N.fast_normalize();
			double theta = 1 - acos(N[1]) / M_PI;
			double phi = (atan2(-N[2], N[0]) + M_PI) / (2.*M_PI);
			mat = queryMaterial(0, theta, phi);
		}

		mat.shadingN = N;
		mat.Ke = Vector(0.,0.,0.);

		if (flip_normals) mat.shadingN = -mat.shadingN;
		triangle_id = -1;
		return true;
	}

	virtual bool reservoir_sampling_intersection(const Ray& d, Vector& P, double &t, MaterialValues &mat, int &triangle_id, int &current_nb_intersections, double min_t, double max_t) const {

		if (has_envmap) { // in our context (used for subsurface scattering), no intersection with an envmap since max_t is finite.
			return false;
		}

		double a = d.direction.getNorm2(); // can't suppose unit norm as objects are scaled by scaling ray length
		double b = dot(d.direction, d.origin - O);
		double c = (d.origin - O).getNorm2() - R * R;

		double delta = b * b - a * c;
		if (delta < 0) return false;

		double sqDelta = sqrt(delta);
		double inva = 1. / a;
		double t2 = (-b + sqDelta)*inva;

		if (t2 < min_t) return false;

		double t1 = (-b - sqDelta)*inva;

		int threadid = omp_get_thread_num();
		const float invmax = 1.f / engine[threadid].max();

		bool has_inter = false;
		if (t1 >= min_t) { // 2 intersections

			if (t1 < max_t) {
				current_nb_intersections++;
				float r1 = engine[threadid]()*invmax;
				if (r1 < 1. / current_nb_intersections) {
					t = t1;
					has_inter = true;
				}
			}
			if (t2 < max_t) {
				current_nb_intersections++;
				float r1 = engine[threadid]()*invmax;
				if (r1 < 1. / current_nb_intersections) {
					t = t2;
					has_inter = true;
				}
			}

		} else { // 1 intersection
			if (t2 < max_t) {
				current_nb_intersections++;
				float r1 = engine[threadid]()*invmax;
				if (r1 < 1. / current_nb_intersections) {
					t = t2;
					has_inter = true;
				}
			}
		}

		if (!has_inter) return false;

		P = d.origin + t * d.direction;
		Vector N = (P - O);


		if (textures.size() != 0 || specularmap.size() != 0 || roughnessmap.size() != 0 || transparent_map.size() != 0 || refr_index_map.size() != 0) {
			N.fast_normalize();
			double theta = 1 - acos(N[1]) / M_PI;
			double phi = (atan2(-N[2], N[0]) + M_PI) / (2.*M_PI);
			mat = queryMaterial(0, theta, phi);
		}

		mat.shadingN = N;
		mat.Ke = Vector(0., 0., 0.);

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
	Plane() { type = OT_PLANE; };
	Plane(const Vector& A, const Vector &N, bool mirror = false) {
		init(A, N, mirror);
	};

	void init(const Vector& A, const Vector &N, bool mirror = false) {
		this->A = A;
		this->vecN = N;
		miroir = mirror;
		name = "Plane";
		type = OT_PLANE;
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

	virtual bool reservoir_sampling_intersection(const Ray& d, Vector& P, double &t, MaterialValues &mat, int &triangle_id, int &current_nb_intersections, double min_t, double max_t) const {
		mat.shadingN = vecN;
		double ddot = dot(d.direction, vecN);
		if (abs(ddot) < 1E-15) return false;
		double curt = dot(A - d.origin, vecN) / ddot;
		if (curt < min_t || curt>= max_t) return false;		

		int threadid = omp_get_thread_num();
		const float invmax = 1.f / engine[threadid].max();
		current_nb_intersections++;
		float r1 = engine[threadid]()*invmax;
		if (r1 >= 1. / current_nb_intersections) {
			return false;
		}

		t = curt;
		P = d.origin + t * d.direction;
		triangle_id = -1;

		double u = P[0] * 0.1;
		double v = P[2] * 0.1;
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

		Vector Point, N;
		fscanf(f, "Point: (%lf, %lf, %lf)\n", &Point[0], &Point[1], &Point[2]);
		fscanf(f, "N: (%lf, %lf, %lf)\n", &N[0], &N[1], &N[2]);
		

		result->init(Point, N);
		return result;
	}

	Vector A, vecN;
};


#ifdef USE_EMBREE

struct MyEmbreeIntersection
{
	RTCIntersectContext context;
	int random_inter_count;
	RTCHit random_hit;
	float random_t;
};


/*static thread_local int embree_random_inter_count[64];
static thread_local RTCHit embree_random_hit[64];
static thread_local float embree_random_t[64];*/

void random_intersectionFilter(const RTCFilterFunctionNArguments* args);
#endif

class Scene {
public:
	Scene() {
		backgroundW = 0; backgroundH = 0; current_time = 0; current_frame = 0; duration = 1;
		double_frustum_start_t = 0;

#ifdef USE_EMBREE
#if _DEBUG
		const char * config = "verbose=3";
#else
		const char * config = "";
#endif		
		// embree init
		embree_device = rtcNewDevice(config);
		embree_scene = rtcNewScene(embree_device);
		//rtcSetSceneFlags(embree_scene, RTC_SCENE_FLAG_DYNAMIC | RTC_SCENE_FLAG_ROBUST);
		rtcSetSceneFlags(embree_scene, RTC_SCENE_FLAG_CONTEXT_FILTER_FUNCTION);
		rtcSetSceneBuildQuality(embree_scene, RTC_BUILD_QUALITY_HIGH);
		for (int i = 0; i < 64; i++) {
			embree_coherent[i] = { RTC_INTERSECT_CONTEXT_FLAG_COHERENT,    0, {} };
			embree_incoherent[i] = { RTC_INTERSECT_CONTEXT_FLAG_INCOHERENT,  0, {} };
			embree_random[i].context = { RTC_INTERSECT_CONTEXT_FLAG_INCOHERENT,  0, {} };			
			rtcInitIntersectContext(&embree_coherent[i]);
			rtcInitIntersectContext(&embree_incoherent[i]);
			rtcInitIntersectContext(&embree_random[i].context);
			embree_coherent[i].flags = RTC_INTERSECT_CONTEXT_FLAG_COHERENT;
			embree_random[i].context.filter = random_intersectionFilter;
			
		}
		embree_bvh_up_to_date = false;
#endif
	};

	~Scene() {
#ifdef USE_EMBREE
		rtcReleaseScene(embree_scene);
		rtcReleaseDevice(embree_device);
#endif
	}

	void prepare_render(bool is_recording);

	void addObject(Object* o);

	void clear() {
#ifdef USE_EMBREE
		for (int i = 0; i < embree_objects.size(); i++) {
			rtcDetachGeometry(embree_scene, embree_objects[i]);
		}
		embree_to_real_objects.clear();
		embree_bvh_up_to_date = false;
		/*for (int i = 0; i < objects.size(); i++) {
			if (objects[i]->type == OT_TRIMESH) {
				TriMesh* g = dynamic_cast<TriMesh*>(objects[i]);
			}
		}*/
#endif

		for (int i = 0; i < objects.size(); i++) {
			delete objects[i];
		}
		objects.clear();
		lumiere = NULL;

	}
	void deleteObject(int id) {
		if (id < 0 || id >= objects.size()) return;

#ifdef USE_EMBREE
		int embreeID = -1;
		for (int i = 0; i < embree_to_real_objects.size(); i++) {
			int realID = embree_to_real_objects[i];
			if (realID == id) {
				embreeID = i;
				rtcDetachGeometry(embree_scene, embree_objects[i]);				
			}
		}
		embree_to_real_objects.erase(embree_to_real_objects.begin() + embreeID);
		for (int i = 0; i < embree_to_real_objects.size(); i++) {
			int realID = embree_to_real_objects[i];
			if (realID >= id) {
				embree_to_real_objects[i]--;
			}
		}
#endif

		delete objects[id];
		objects.erase(objects.begin() + id);
	}


	bool intersection(const Ray& d, Vector& P, int &sphere_id, double &min_t, MaterialValues &mat, int &triangle_id, bool avoid_ghosts = false, bool isCoherent = false) const;


	bool intersection_shadow(const Ray& d, double &min_t, double dist_light, bool avoid_ghosts = false) const;

	bool get_random_intersection(const Ray& d, Vector& P, int &sphere_id, double &min_t, MaterialValues &mat, int &triangle_id, double tmin, double tmax, bool avoid_ghosts = false) const;

	void clear_background() {
		background.clear();
		backgroundW = 0;
		backgroundH = 0;
		backgroundfilename = std::string("");
	}

	void load_background(const char* filename, double gamma) {
		std::vector<unsigned char> bg;
		backgroundfilename = std::string(filename);
		load_bmp(filename, bg, backgroundW, backgroundH);
		background.resize(bg.size());
		for (int i = 0; i < backgroundW*backgroundH * 3; i++)
			background[i] = std::pow(bg[i]/255., gamma)* 196964.699; // will be gamma corrected during rendering ; beware: values are kept between [0,255^2.2] for consistency with the renderer
	}

	std::vector<Object*> objects;
	std::vector<double> background;
	int backgroundW, backgroundH;
	std::string backgroundfilename;
	Sphere *lumiere;
	double intensite_lumiere;
	double envmap_intensity;
	double fog_density, fog_absorption, fog_density_decay, fog_absorption_decay;
	int nbframes, current_frame;
	double duration, current_time;
	int fog_type; // 0 : uniform, 1: exponential
	int fog_phase_type; // 0 isotropic ; 1 Schlick 2 rayleigh
	double double_frustum_start_t;
	double phase_aniso;

#ifdef USE_EMBREE
	RTCDevice embree_device;
	RTCScene embree_scene;
	mutable RTCIntersectContext embree_coherent[64];   //  primary rays
	mutable RTCIntersectContext embree_incoherent[64]; // secondary/visibility rays
	mutable MyEmbreeIntersection embree_random[64]; // reservoir sampling of intersections
	
	bool embree_bvh_up_to_date;
	std::vector<int> embree_objects;
	std::vector<int> embree_to_real_objects;
#endif
};

