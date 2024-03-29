#define _CRT_SECURE_NO_WARNINGS
#include "Geometry.h"
#include <list>
#include <set>
#include <string> 
#include <algorithm>
#include "utils.h"
#include "TriangleMesh.h"
#include "PointSet.h"

Object* Object::create_from_file(FILE* f, Scene* scene, const char* replacedNames) {

	char line[255];
	fscanf(f, "%[^\n]\n", line);
	if (line[4] == 'M') { // mesh
		return TriMesh::create_from_file(f, scene, replacedNames);
	}
	if (line[4] == 'S') { // sphere
		return Sphere::create_from_file(f);
	}
	if (line[4] == 'P' && line[5] == 'L') { // Plane
		return Plane::create_from_file(f);
	}
	if (line[4] == 'P' && line[5] == 'O') { // PointSet
		return PointSet::create_from_file(f, replacedNames);
	}
  return nullptr;
}


void load_bmp(const char* filename, std::vector<unsigned char> &tex, int &W, int& H) {
	size_t w, h;
	load_image(filename, tex, w, h);
	W = w;
	H = h;
	
	/*FILE* f;
	f = fopen(filename, "rb");
	unsigned char info[54];
	fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

	W = *(int*)&info[18]; // extract image height and width from header
	H = *(int*)&info[22];

	int size = 3 * W * H;
	tex.resize(size); // allocate 3 bytes per pixel
	fread(&tex[0], sizeof(unsigned char), size, f); // read the rest of the data at once
	fclose(f);

	for (int i = 0; i < size; i += 3) {
		std::swap(tex[i], tex[i + 2]);
	}*/
}


void Object::add_normalmap(const char* filename) {
	normal_map.push_back(Texture(filename, 2, Vector(0., 0., 1.)));
}

void Object::set_normalmap(const char* filename, int idx) {
	if (idx >= normal_map.size()) return;
	//normal_map[idx] = Texture(filename, 2, Vector(0., 0., 1.));
	normal_map[idx].loadNormals(filename);
}

void Object::add_alphamap(const char* filename) {
	alphamap.push_back(Texture(filename, 3, Vector(1., 1., 1.)));
}

void Object::set_texture(const char* filename, int idx) {
	if (idx >= textures.size()) return;
	//textures[idx] = Texture(filename, 0, Vector(1., 1., 1.));
	textures[idx].loadColors(filename);
}

void Object::add_texture(const char* filename) {
	textures.push_back(Texture(filename, 0, Vector(1., 1., 1.)));
}
void Object::add_specularmap(const char* filename) {
	specularmap.push_back(Texture(filename, 1, Vector(0., 0., 0.)));
}
void Object::set_specularmap(const char* filename, int idx) {
	if (idx >= specularmap.size()) return;
	//specularmap[idx] = Texture(filename, 1, Vector(0., 0., 0.));
	specularmap[idx].loadColors(filename);
}
void Object::add_subsurface(const char* filename) {
	subsurface.push_back(Texture(filename, 1, Vector(0., 0., 0.)));
}
void Object::set_subsurface(const char* filename, int idx) {
	if (idx >= subsurface.size()) return;
	//specularmap[idx] = Texture(filename, 1, Vector(0., 0., 0.));
	subsurface[idx].loadColors(filename);
}
void Object::add_transp_map(const char* filename) {
	transparent_map.push_back(Texture(filename, 5, Vector(1., 1., 1.)));
}
void Object::set_transp_map(const char* filename, int idx) {
	if (idx >= transparent_map.size()) return;
	transparent_map[idx].loadColors(filename);
}
void Object::set_col_transp(float col, int idx) {
	if (idx >= transparent_map.size()) return;
	transparent_map[idx] = Texture("Null", 5, Vector(col, col, col));
}
void Object::remove_transp(int id) {
	transparent_map.erase(transparent_map.begin() + id);
}
void Object::swap_transp(int id1, int id2) {
	std::swap(transparent_map[id1], transparent_map[id2]);
}
void Object::add_col_transp(float col) {
	transparent_map.push_back(Texture("Null", 5, Vector(col, col, col)));
}
void Object::add_refr_map(const char* filename) {
	refr_index_map.push_back(Texture(filename, 6, Vector(1., 1., 1.)));
}
void Object::set_refr_map(const char* filename, int idx) {
	if (idx >= refr_index_map.size()) return;
	refr_index_map[idx].loadColors(filename);
}
void Object::set_col_refr(float col, int idx) {
	if (idx >= refr_index_map.size()) return;
	refr_index_map[idx] = Texture("Null", 6, Vector(col, col, col));
}
void Object::remove_refr(int id) {
	refr_index_map.erase(refr_index_map.begin() + id);
}
void Object::swap_refr(int id1, int id2) {
	std::swap(refr_index_map[id1], refr_index_map[id2]);
}
void Object::add_col_refr(float col) {
	refr_index_map.push_back(Texture("Null", 6, Vector(col, col, col)));
}
void Object::add_roughnessmap(const char* filename) {
	roughnessmap.push_back(Texture(filename, 4, Vector(1., 1., 1.)));
}
void Object::set_roughnessmap(const char* filename, int idx) {
	if (idx >= roughnessmap.size()) return;
	roughnessmap[idx] = Texture(filename, 4, Vector(1., 1., 1.));
}

void Object::set_alphamap(const char* filename, int idx) {
	if (idx >= alphamap.size()) return;
	alphamap[idx] = Texture(filename, 3, Vector(1., 1., 1.));
}
void Object::set_col_alpha(float col, int idx) {
	if (idx >= alphamap.size()) return;
	alphamap[idx] = Texture("Null", 3, Vector(col, col, col));
}

void Object::add_col_specular(const Vector& col) {
//	std::string name = std::string("Color: (") + std::to_string((int)(col[0] * 255)) + std::string(", ") + std::to_string((int)(col[1] * 255)) + std::string(", ") + std::to_string((int)(col[2] * 255)) + std::string(")");
	specularmap.push_back(Texture("Null", 1, col));
}
void Object::add_col_subsurface(const Vector& col) {
	subsurface.push_back(Texture("Null", 1, col));
}

void Object::add_col_roughness(const Vector& col) {	
	//std::string name = std::string("Color: (") + std::to_string((int)(col[0])) + std::string(", ") + std::to_string((int)(col[1])) + std::string(", ") + std::to_string((int)(col[2])) + std::string(")");
	roughnessmap.push_back(Texture("Null", 4, col));
}

void Object::set_col_roughness(const Vector& col, int idx) {
	if (idx >= roughnessmap.size()) return;

	//std::string name = std::string(std::string("Color: (") + std::to_string((int)(col[0])) + std::string(", ") + std::to_string((int)(col[1])) + std::string(", ") + std::to_string((int)(col[2])) + std::string(")"));
	//roughnessmap[idx] = Texture(name.c_str(), 4, col);
	roughnessmap[idx].multiplier = col;
}

void Object::set_col_specular(const Vector& col, int idx) {
	if (idx >= specularmap.size()) return;

	//std::string name = std::string("Color: (") + std::to_string((int)(col[0] * 255)) + std::string(", ") + std::to_string((int)(col[1] * 255)) + std::string(", ") + std::to_string((int)(col[2] * 255)) + std::string(")");
	specularmap[idx].multiplier = col;// = Texture(name.c_str(), 1, col);
}
void Object::set_col_subsurface(const Vector& col, int idx) {
	if (idx >= subsurface.size()) return;

	//std::string name = std::string("Color: (") + std::to_string((int)(col[0] * 255)) + std::string(", ") + std::to_string((int)(col[1] * 255)) + std::string(", ") + std::to_string((int)(col[2] * 255)) + std::string(")");
	subsurface[idx].multiplier = col;// = Texture(name.c_str(), 1, col);
}


void Object::remove_texture(int id) {
	textures.erase(textures.begin() + id);
}
void Object::remove_alpha(int id) {
	alphamap.erase(alphamap.begin() + id);
}
void Object::remove_specular(int id) {
	specularmap.erase(specularmap.begin() + id);
}
void Object::remove_normal(int id) {
	normal_map.erase(normal_map.begin() + id);
}
void Object::remove_roughness(int id) {
	roughnessmap.erase(roughnessmap.begin() + id);
}
void Object::remove_subsurface(int id) {
	subsurface.erase(subsurface.begin() + id);
}
void Object::swap_roughness(int id1, int id2) {
	std::swap(roughnessmap[id1], roughnessmap[id2]);
}
void Object::swap_textures(int id1, int id2) {
	std::swap(textures[id1], textures[id2]);
}
void Object::swap_normal(int id1, int id2) {
	std::swap(normal_map[id1], normal_map[id2]);
}
void Object::swap_specular(int id1, int id2) {
	std::swap(specularmap[id1], specularmap[id2]);
}
void Object::swap_subsurface(int id1, int id2) {
	std::swap(subsurface[id1], subsurface[id2]);
}
void Object::swap_alpha(int id1, int id2) {
	std::swap(alphamap[id1], alphamap[id2]);
}
void Object::add_col_texture(const Vector& col) {
		
	//std::string name = (std::string("Color: (") + std::to_string((int)(col[0] * 255)) + std::string(", ") + std::to_string((int)(col[1] * 255)) + std::string(", ") + std::to_string((int)(col[2] * 255)) + std::string(")"));
	textures.push_back(Texture("Null", 0, col));
}

void Object::set_col_texture(const Vector& col, int idx) {
	if (idx >= textures.size()) return;

	//std::string name = std::string("Color: (") + std::to_string((int)(col[0] * 255)) + std::string(", ") + std::to_string((int)(col[1] * 255)) + std::string(", ") + std::to_string((int)(col[2] * 255)) + std::string(")");
	textures[idx].multiplier = col;// = Texture("Null", 0, col);
}

void Object::add_null_normalmap() {
	normal_map.push_back(Texture::defaultNormal());
}
void Object::set_null_normalmap(int idx) {
	if (idx >= textures.size()) return;
	normal_map[idx].clear_texture();
}
void Object::add_col_alpha(float col) {
	alphamap.push_back(Texture("Null", 3, Vector(col, col, col)));
}



void Scene::addObject(Object* o) {
	objects.push_back(o);

	TriMesh* g = dynamic_cast<TriMesh*>(o);
	castToMesh.push_back(g);

#ifdef USE_EMBREE
	if (g) {
		g->instance_geom = rtcNewGeometry(embree_device, RTC_GEOMETRY_TYPE_INSTANCE);
		rtcSetGeometryInstancedScene(g->instance_geom, g->embree_scene_for_instance);
		rtcSetGeometryTimeStepCount(g->instance_geom, 1);
		int geomID = rtcAttachGeometry(embree_scene, g->instance_geom);
		rtcReleaseGeometry(g->instance_geom);
		embree_objects.push_back(geomID);
		embree_to_real_objects.resize(geomID + 1);
		embree_to_real_objects[geomID] = objects.size() - 1;
		float trans[16];
		for (int i = 0; i < 12; i++) {
			trans[i] = g->trans_matrix[i];
		}
		trans[12] = 0; trans[13] = 0; trans[14] = 0; trans[15] = 1; 
		rtcSetGeometryTransform(g->instance_geom, 0, RTC_FORMAT_FLOAT3X4_ROW_MAJOR, trans);
		rtcCommitGeometry(g->instance_geom);
		rtcSetSceneBuildQuality(embree_scene, RTC_BUILD_QUALITY_HIGH);
		rtcCommitScene(embree_scene);
		embree_bvh_up_to_date = false;
	}
	prepare_render(false);
#endif
}

void Scene::prepare_render(bool is_recording) {

	for (int i = 0; i < objects.size(); i++) {
		objects[i]->build_matrix(current_frame, is_recording); // time = 0 here
	}

#ifdef USE_EMBREE
	if (!embree_bvh_up_to_date) {
		for (int i = 0; i < objects.size(); i++) {

			if (objects[i]->type == OT_TRIMESH) {
				TriMesh* g = castToMesh[i]; // dynamic_cast<TriMesh*>(objects[i]);
				float trans[16];
				for (int j = 0; j < 12; j++) {
					trans[j] = g->trans_matrix[j];
				}
				trans[12] = 0; trans[13] = 0; trans[14] = 0; trans[15] = 1;
				rtcSetGeometryTransform(g->instance_geom, 0, RTC_FORMAT_FLOAT3X4_ROW_MAJOR, trans);
				rtcCommitGeometry(g->instance_geom);
			}
		}
		rtcCommitScene(embree_scene);
		embree_bvh_up_to_date = true;
	}


#endif
}


#ifdef USE_EMBREE

void random_intersectionFilter(const RTCFilterFunctionNArguments* args) {
	// avoid crashing when debug visualizations are used 
	if (args->context == nullptr) return;

	assert(args->N == 1);
	int* valid = args->valid;
	
	//ignore inactive rays 
	if (valid[0] != -1) return;

	int rayid = omp_get_thread_num(); // ((RTCRay*)args->ray)->id;
	int threadid = rayid;
	const float invmax = 1.f / engine[rayid].max();

	MyEmbreeIntersection* inter = (MyEmbreeIntersection*)(args->context);
	inter->random_inter_count++;
	
	float r1 = engine[rayid]()*invmax;
	if (r1 <= 1. / inter->random_inter_count) {
		inter->random_hit = *(((RTCHit*)args->hit)); //RTCHitN_instID(args->hit, args->N, 0, 0);
		inter->random_t = ((RTCRay*)args->ray)->tfar;
	}
	valid[0] = 0;	
}
#endif 

// super slow version ; returns a random intersection in tmin tmax. Uses reservoir sampling. If sphere_id!=-1, considers only intersections with this object
bool Scene::get_random_intersection(const Ray& d, Vector& P, int &sphere_id, float &min_t, MaterialValues &mat, int &triangle_id, float tmin, float tmax, bool avoid_ghosts) const {
	/*bool hasinter = false;
	double startt = tmin;
	int ninter = 0;
	for (int i = 0; i < 10; i++) {
		Ray start(d.origin + startt * d.direction, d.direction, d.time);
		Vector localP;
		MaterialValues localmat;
		double localt;
		int localid = sphere_id, localtri;
		bool localhasinter = intersection(start, localP, localid, localt, localmat, localtri, avoid_ghosts, false);
		if (!localhasinter) return hasinter;		
		if (localt > tmax - startt) return hasinter;
		startt += localt + 1E-5;
		if ((sphere_id != -1) && (localid != sphere_id)) continue;
		ninter++;
		hasinter = true;
		if (ninter == 1) {
			P = localP;
			sphere_id = localid;
			min_t = startt;
			mat = localmat;
			triangle_id = localtri;

		} else {
			int threadid = omp_get_thread_num();
			float invmax = 1.f / engine[threadid].max();
			float r1 = engine[threadid]()*invmax;
			if (r1 < 1.f / ninter) {
				P = localP;
				sphere_id = localid;
				min_t = startt;
				mat = localmat;
				triangle_id = localtri;
			}
		}
		

	}
	return hasinter;
	*/

	bool has_inter = false;
	min_t = 1E99;


	float t;
	int nb_intersections = 0;

	for (int i = 0; i < objects.size(); i++) {
		if (sphere_id != -1 && sphere_id != i) continue;
		if (avoid_ghosts && objects[i]->ghost) continue;

#ifdef USE_EMBREE
		if (objects[i]->type == OT_TRIMESH) {
			continue;
		}
#endif

		Vector transformed_dir = objects[i]->apply_inverse_rotation_scaling(d.direction);
		Vector new_origin = objects[i]->apply_inverse_transformation(d.origin);
		Ray transformed_ray(new_origin, transformed_dir, d.time);

		has_inter = has_inter || objects[i]->reservoir_sampling_intersection(transformed_ray, P, min_t, mat, triangle_id, nb_intersections, tmin, tmax);
		
	}


#ifdef USE_EMBREE
	RTCRayHit embree_ray;

	embree_ray =
	{
		{
			(float)d.origin[0],  (float)d.origin[1], (float)d.origin[2],         // origin - Visual Studio bug requires explicitly casting to float
			(float)tmin,                                              // tnear    
			(float)d.direction[0], (float)d.direction[1], (float)d.direction[2], // direction
			0 ,                                             // time
			(float)tmax,              // tfar
			avoid_ghosts ? ((~0u) - 1) : (~0u),                    // mask
			0,                      // ray id 
			0                       // ray flags
		},

		{
			0, 0, 0,                    // intersection normal
			0, 0,                       // intersection u, v
			RTC_INVALID_GEOMETRY_ID,    // primitive ID
			RTC_INVALID_GEOMETRY_ID,    // geometry ID
			{ RTC_INVALID_GEOMETRY_ID } // instance ID
		}
	};

	int threadid = omp_get_thread_num();
	
	embree_random[threadid].random_hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
	embree_random[threadid].random_inter_count = nb_intersections;	

	rtcIntersect1(embree_scene, &embree_random[threadid].context, &embree_ray);
	embree_ray.hit = embree_random[threadid].random_hit;
	nb_intersections = embree_random[threadid].random_inter_count;
	int embreeObjectID = embree_ray.hit.instID[0];
	has_inter = has_inter || (embreeObjectID != RTC_INVALID_GEOMETRY_ID);

	if (has_inter) {

		if (embreeObjectID != RTC_INVALID_GEOMETRY_ID) {  // embree
			min_t = embree_random[threadid].random_t;
			triangle_id = embree_ray.hit.primID;
			P = d.origin + min_t * d.direction;
			int systemObjectID = embree_to_real_objects[embreeObjectID];
			sphere_id = systemObjectID;
			TriMesh* g = castToMesh[systemObjectID];//  dynamic_cast<TriMesh*>(objects[systemObjectID]);
			g->getMaterial(triangle_id, 1 - embree_ray.hit.u - embree_ray.hit.v, embree_ray.hit.u, embree_ray.hit.v, mat);
		} else {
			P = objects[sphere_id]->apply_transformation(P);
			mat.shadingN.fast_normalize();
		}
		mat.shadingN = objects[sphere_id]->apply_rotation(mat.shadingN);
	}
#else

	if (has_inter) {
		P = objects[sphere_id]->apply_transformation(P);
		mat.shadingN = objects[sphere_id]->apply_rotation(mat.shadingN);
	}
	mat.shadingN.fast_normalize();
#endif

	
	

	return has_inter;
}

void Scene::first_intersection_batch(int batch_size) { // dont forget to set W and H

	int threadid = omp_get_thread_num();
	firstIntersection_P[threadid].resize(batch_size);
	firstIntersection_sphere_id[threadid].resize(batch_size);
	firstIntersection_min_t[threadid].resize(batch_size);
	firstIntersection_mat[threadid].resize(batch_size);
	firstIntersection_triangle_id[threadid].resize(batch_size);
	firstIntersection_has_inter[threadid].resize(batch_size);
	std::fill(firstIntersection_has_inter[threadid].begin(), firstIntersection_has_inter[threadid].end(), false);
	std::fill(firstIntersection_min_t[threadid].begin(), firstIntersection_min_t[threadid].end(), 1E9);
		
#ifdef USE_EMBREE
	std::vector<RTCRayHit16> allhits(batch_size / 16);

	for (int id = 0; id < batch_size /16; id++) {

		for (int k = 0; k < 16; k++) {
			allhits[id].ray.org_x[k] = firstIntersection_Ray[threadid][id * 16 + k].origin[0];
			allhits[id].ray.org_y[k] = firstIntersection_Ray[threadid][id * 16 + k].origin[1];
			allhits[id].ray.org_z[k] = firstIntersection_Ray[threadid][id * 16 + k].origin[2];

			allhits[id].ray.tnear[k] = 0;

			allhits[id].ray.dir_x[k] = firstIntersection_Ray[threadid][id * 16 + k].direction[0];
			allhits[id].ray.dir_y[k] = firstIntersection_Ray[threadid][id * 16 + k].direction[1];
			allhits[id].ray.dir_z[k] = firstIntersection_Ray[threadid][id * 16 + k].direction[2];
			
			allhits[id].ray.time[k] = 0;
			allhits[id].ray.tfar[k] = 1E8;
			allhits[id].ray.mask[k] = ~0u;
			allhits[id].ray.id[k] = id*16+k;
			allhits[id].ray.flags[k] = 0;
			allhits[id].hit.primID[k] = RTC_INVALID_GEOMETRY_ID;
			allhits[id].hit.geomID[k] = RTC_INVALID_GEOMETRY_ID;
			allhits[id].hit.instID[0][k] = RTC_INVALID_GEOMETRY_ID;
		}
	}

	rtcIntersectNM(embree_scene, &embree_coherent[threadid], (RTCRayHitN*)&allhits[0], 16, ceil(batch_size / 16.), sizeof(RTCRayHit16));

	for (int id = 0; id < batch_size / 16; id++) {
		for (int k = 0; k < 16; k++) {		
			const RTCRayHit16 &embree_ray = allhits[id];
				int embreeObjectID = embree_ray.hit.instID[0][k];
				if (embreeObjectID != RTC_INVALID_GEOMETRY_ID) {  // embree
					Ray d(Vector(embree_ray.ray.org_x[k], embree_ray.ray.org_y[k], embree_ray.ray.org_z[k]), Vector(embree_ray.ray.dir_x[k], embree_ray.ray.dir_y[k], embree_ray.ray.dir_z[k]), 0);
					int rayId = embree_ray.ray.id[k];
					firstIntersection_has_inter[threadid][rayId] = true;
					firstIntersection_min_t[threadid][rayId] = embree_ray.ray.tfar[k];
					firstIntersection_triangle_id[threadid][rayId] = embree_ray.hit.primID[k];
					firstIntersection_P[threadid][rayId] = d.origin + embree_ray.ray.tfar[k] * d.direction;
					int systemObjectID = embree_to_real_objects[embreeObjectID];
					firstIntersection_sphere_id[threadid][rayId] = systemObjectID;
					TriMesh* g = castToMesh[systemObjectID]; // dynamic_cast<TriMesh*>(objects[systemObjectID]);
					g->getMaterial(embree_ray.hit.primID[k], 1 - embree_ray.hit.u[k] - embree_ray.hit.v[k], embree_ray.hit.u[k], embree_ray.hit.v[k], firstIntersection_mat[threadid][rayId]);
				} 
		}
	}


	for (int id = 0; id < batch_size; id++) {
		
		bool found_not_mesh = false;
		for (int i = 0; i < objects.size(); i++) {
			if (objects[i]->type == OT_TRIMESH) continue;
			if (i == 1 && firstIntersection_has_inter[threadid][id]) continue; // we won't do better instesecting the envmap ;)

			const Ray &d = firstIntersection_Ray[threadid][id];
			Vector transformed_dir = objects[i]->apply_inverse_rotation_scaling(d.direction);
			Vector new_origin = objects[i]->apply_inverse_transformation(d.origin);
			Ray transformed_ray(new_origin, transformed_dir, d.time);

			Vector localP;
			int triangle_id;
			float min_t;
			MaterialValues localmat;
			float t;

			bool local_has_inter = objects[i]->intersection(transformed_ray, localP, t, localmat, firstIntersection_min_t[threadid][id], triangle_id);

			if (local_has_inter) {
				if (t < firstIntersection_min_t[threadid][id]) {
					firstIntersection_has_inter[threadid][id] = true;
					firstIntersection_min_t[threadid][id] = t;

					//P = objects[sphere_id]->apply_transformation(localP);
					firstIntersection_P[threadid][id] = localP;
					firstIntersection_sphere_id[threadid][id] = i;
					firstIntersection_mat[threadid][id] = localmat;
					found_not_mesh = true;
					//mat.shadingN = objects[i]->apply_rotation(localmat.shadingN);
				}
			}
		}

		if (found_not_mesh) {
			firstIntersection_P[threadid][id] = objects[firstIntersection_sphere_id[threadid][id]]->apply_transformation(firstIntersection_P[threadid][id]);
		}

		firstIntersection_mat[threadid][id].shadingN = objects[firstIntersection_sphere_id[threadid][id]]->apply_rotation(firstIntersection_mat[threadid][id].shadingN);
		firstIntersection_mat[threadid][id].shadingN.fast_normalize();
	}

#else

//#pragma omp parallel for schedule(dynamic, 1)
	for (int id = 0; id < batch_size; id++) {
		firstIntersection_has_inter[threadid][id] = intersection(firstIntersection_Ray[threadid][id], firstIntersection_P[threadid][id], firstIntersection_sphere_id[threadid][id], firstIntersection_min_t[threadid][id], firstIntersection_mat[threadid][id], firstIntersection_triangle_id[threadid][id], false, true);
	}

#endif
}


bool Scene::intersection(const Ray& d, Vector& P, int &sphere_id, float &min_t, MaterialValues &mat, int &triangle_id, bool avoid_ghosts, bool isCoherent) const {

	bool has_inter = false;
	min_t = 1E99;


	Vector localP;
	MaterialValues localmat;
	float t;

	for (int i = 0; i < objects.size(); i++) {
		if (avoid_ghosts && objects[i]->ghost) continue;

#ifdef USE_EMBREE
		if (objects[i]->type == OT_TRIMESH) continue;
#endif

		Vector transformed_dir = objects[i]->apply_inverse_rotation_scaling(d.direction);
		Vector new_origin = objects[i]->apply_inverse_transformation(d.origin);
		Ray transformed_ray(new_origin, transformed_dir, d.time);

		bool local_has_inter = objects[i]->intersection(transformed_ray, localP, t, localmat, min_t, triangle_id);

		if (local_has_inter) {
			if (t < min_t) {
				has_inter = true;
				min_t = t;

				//P = objects[sphere_id]->apply_transformation(localP);
				P = localP;
				sphere_id = i;
				mat = localmat;
				//mat.shadingN = objects[i]->apply_rotation(localmat.shadingN);
			}
		}
	}


#ifdef USE_EMBREE
	RTCRayHit embree_ray;

	embree_ray =
	{
		{
			(float)d.origin[0],  (float)d.origin[1], (float)d.origin[2],         // origin - Visual Studio bug requires explicitly casting to float
			0,                                              // tnear    
			(float)d.direction[0], (float)d.direction[1], (float)d.direction[2], // direction
			0 ,                                             // time
			(float)min_t,              // tfar
			avoid_ghosts?((~0u)-1): (~0u),                    // mask
			0,                      // ray id 
			0                       // ray flags
		},

		{
			0, 0, 0,                    // intersection normal
			0, 0,                       // intersection u, v
			RTC_INVALID_GEOMETRY_ID,    // primitive ID
			RTC_INVALID_GEOMETRY_ID,    // geometry ID
			{ RTC_INVALID_GEOMETRY_ID } // instance ID
		}
	};

	if (isCoherent) 
		rtcIntersect1(embree_scene, &embree_coherent[omp_get_thread_num()], &embree_ray);  
	else
		rtcIntersect1(embree_scene, &embree_incoherent[omp_get_thread_num()], &embree_ray);

	int embreeObjectID = embree_ray.hit.instID[0];
	has_inter = has_inter || (embreeObjectID != RTC_INVALID_GEOMETRY_ID);

	if (has_inter) {

		if (embreeObjectID != RTC_INVALID_GEOMETRY_ID) {  // embree
			min_t = embree_ray.ray.tfar;
			triangle_id = embree_ray.hit.primID;
			P = d.origin + min_t * d.direction;
			int systemObjectID = embree_to_real_objects[embreeObjectID];
			sphere_id = systemObjectID;
			TriMesh* g = castToMesh[systemObjectID]; // dynamic_cast<TriMesh*>(objects[systemObjectID]);
			g->getMaterial(triangle_id, 1 - embree_ray.hit.u - embree_ray.hit.v, embree_ray.hit.u, embree_ray.hit.v, mat);			
		} else {
			P = objects[sphere_id]->apply_transformation(P);			
		}
		mat.shadingN = objects[sphere_id]->apply_rotation(mat.shadingN);
	}
#else

	if (has_inter) {
		P = objects[sphere_id]->apply_transformation(P);
		mat.shadingN = objects[sphere_id]->apply_rotation(mat.shadingN);

    }
#endif

	mat.shadingN.fast_normalize();

	return has_inter;

}


bool Scene::intersection_shadow(const Ray& d, float &min_t, float dist_light, bool avoid_ghosts, bool isCoherent) const {

	min_t = 1E99;

#ifdef USE_EMBREE
	RTCRay embree_ray;

	embree_ray =
	{	
			(float)d.origin[0],  (float)d.origin[1], (float)d.origin[2],         // origin - Visual Studio bug requires explicitly casting to float
			0,                                              // tnear    
			(float)d.direction[0], (float)d.direction[1], (float)d.direction[2], // direction
			0 ,                                             // time
			(float)(dist_light*0.9991),              // tfar
			avoid_ghosts ? ((~0u) - 1) : (~0u),                    // mask
			0,                      // ray id 
			0                       // ray flags
	};

	if (isCoherent)
		rtcOccluded1(embree_scene, &embree_coherent[omp_get_thread_num()], &embree_ray);  
	else
		rtcOccluded1(embree_scene, &embree_incoherent[omp_get_thread_num()], &embree_ray);  

	if (embree_ray.tfar<0) {
		return true;
	}	

#endif

	for (int i = 0; i < objects.size(); i++) {
		if (avoid_ghosts && objects[i]->ghost) continue;

#ifdef USE_EMBREE
		if (objects[i]->type == OT_TRIMESH) continue;
#endif

		Vector transformed_dir = objects[i]->apply_inverse_rotation_scaling(d.direction);
		Vector new_origin = objects[i]->apply_inverse_transformation(d.origin);
		Ray transformed_ray(new_origin, transformed_dir, d.time);

		float t;

		bool local_has_inter = objects[i]->intersection_shadow(transformed_ray, t, min_t, dist_light);

		if (local_has_inter) {
			if (t < dist_light*0.999) {
				return true;
			}
		}
	}

	return false;
}