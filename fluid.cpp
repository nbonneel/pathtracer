#include "fluid.h"
#include <omp.h>

BBoxd Fluid::build_bbox(int frame, int i0, int i1) {

	BBoxd result;
	result.bounds[1] = particles[frame][i0] + Vectord(radius, radius, radius);
	result.bounds[0] = particles[frame][i0] - Vectord(radius, radius, radius);
	for (int i = i0; i < i1; i++) { // indice de triangle
		result.bounds[0] = min(result.bounds[0], particles[frame][i] - Vectord(radius, radius, radius));
		result.bounds[1] = max(result.bounds[1], particles[frame][i] + Vectord(radius, radius, radius));
	}
	return result;
}

BBoxd Fluid::build_centers_bbox(int frame, int i0, int i1) {

	BBoxd result;
	result.bounds[1] = particles[frame][i0];
	result.bounds[0] = particles[frame][i0];
	for (int i = i0; i < i1; i++) { // indice de triangle
		result.bounds[0] = min(result.bounds[0], particles[frame][i]);
		result.bounds[1] = max(result.bounds[1], particles[frame][i]);
	}
	return result;
}

void Fluid::build_bvh(int frame, int i0, int i1) {
	return;
	bvh.nodes.clear();
	bvh.bbox = build_bbox(frame, i0, i1);
	bvh.nodes.reserve(Nparticles * 2);
	build_bvh_recur(frame, &bvh, 0, i0, i1);	
}

void Fluid::build_grid(int frame) {
	accelgrid.clear();
	accelgridindices.clear();
	accelgrid.resize(Nx*Ny*Nz, 0);
	accelgridindices.resize(Nx*Ny*Nz);
	int ni = std::ceil(radius / dx[0]);
	int nj = std::ceil(radius / dx[1]);
	int nk = std::ceil(radius / dx[2]);

	for (int id = 0; id < particles[frame].size(); id++) {
		int minx = std::max(0, (int)((particles[frame][id][0] - radius - extent.bounds[0][0]) / dx[0]));
		int miny = std::max(0, (int)((particles[frame][id][1] - radius - extent.bounds[0][1]) / dx[1]));
		int minz = std::max(0, (int)((particles[frame][id][2] - radius - extent.bounds[0][2]) / dx[2]));
		int maxx = std::min(Nx - 1, (int)((particles[frame][id][0] + radius - extent.bounds[0][0]) / dx[0]));
		int maxy = std::min(Ny - 1, (int)((particles[frame][id][1] + radius - extent.bounds[0][1]) / dx[1]));
		int maxz = std::min(Nz - 1, (int)((particles[frame][id][2] + radius - extent.bounds[0][2]) / dx[2]));

		for (int i = minz; i <= maxz; i++) {
			for (int j = miny; j <= maxy; j++) {
				for (int k = minx; k <= maxx; k++) {
					accelgrid[i*Nx*Ny + j * Nx + k] = 1;
					accelgridindices[i*Nx*Ny + j * Nx + k].push_back(id);
				}
			}
		}
	}
}


bool Fluid::intersection_transparent2(const Ray& d, Vector& P, float &t, MaterialValues &mat, float cur_best_t, int &triangle_id) const {
	t = cur_best_t;
	bool has_inter = false;
	float t_box_left, t_box_right;
	int best_index = -1;
	bool goleft, goright;
	Vector localP, localN;
	float localt;
	int frame = d.time;

	Ray invd(d.origin, Vector(1. / d.direction[0], 1. / d.direction[1], 1. / d.direction[2]), d.time);
	char signs[3];
	signs[0] = (invd.direction[0] >= 0) ? 1 : 0;
	signs[1] = (invd.direction[1] >= 0) ? 1 : 0;
	signs[2] = (invd.direction[2] >= 0) ? 1 : 0;

	if (!extent.intersection_invd(invd, signs, t_box_left)) return false;
	
	int threadid = omp_get_thread_num();
	int interid = 0;
	Vectord startingP = toVecd(d.origin) + (t_box_left+1E-9) * toVecd(d.direction);
	int boxid[3] = { static_cast<int>(std::round((startingP[0] - extent.bounds[0][0]) / dx[0])), static_cast<int>(std::round((startingP[1] - extent.bounds[0][1]) / dx[1] )), static_cast<int>(std::round((startingP[2] - extent.bounds[0][2]) / dx[2] ))};
	if (boxid[0] < 0 || boxid[1] < 0 || boxid[2] < 0 || boxid[0] >= Nx || boxid[1] >= Ny || boxid[2] >= Nz) return false;

	char stepX = signs[0] * 2 - 1;
	char stepY = signs[1] * 2 - 1;
	char stepZ = signs[2] * 2 - 1;
	Vectord step(stepX, stepY, stepZ);

	Vectord tMax = (((Vectord(boxid[0], boxid[1], boxid[2])+Vectord(signs[0], signs[1], signs[2]))*dx+extent.bounds[0])-startingP) / toVecd(d.direction); // next ray-plane intersection
	Vectord tDelta = (step * dx) / toVecd(d.direction);
	
	while (true) {		
		int voxid = boxid[2] * Nx*Ny + boxid[1] * Nx + boxid[0];
		if (accelgrid[voxid]) {
			int nids = accelgridindices[voxid].size();
			const int* ids = &accelgridindices[voxid][0];
			bool curhasinter = false;
			for (int i = 0; i < nids; i++) {
				float t1, t2;
				if (SimplerSphere(particles[frame][ids[i]], radius).both_intersections(d, t1, t2)) {
					has_inter = true;
					curhasinter = true;
					allts[threadid][interid] = std::make_pair(t1, std::make_pair(t2, ids[i]));
					interid++;
				}
			}
			if (!curhasinter && has_inter) break;
		}
		if (tMax[0] < tMax[1]){
			if (tMax[0] < tMax[2]) {
				boxid[0] += step[0];
				if (boxid[0] >= Nx || boxid[0]<0) break;
				tMax[0] += tDelta[0];
			} else {
				boxid[2] += step[2];
				if (boxid[2] >= Nz || boxid[2] < 0) break;
				tMax[2] += tDelta[2];
			}
		} else {
			if (tMax[1] < tMax[2]) {
				boxid[1] += step[1];
				if (boxid[1] >= Ny || boxid[1] < 0) break;
				tMax[1] += tDelta[1];
			} else {
				boxid[2] += step[2];
				if (boxid[2] >= Nz || boxid[2] < 0) break;
				tMax[2] += tDelta[2];
			}
		}
	}

	if (has_inter) {
		std::sort(allts[threadid], allts[threadid] + interid + 1);
		if (allts[threadid][0].first > 0) {
			best_index = allts[threadid][0].second.second;
			t = allts[threadid][0].first;
		} else {
			int k = 0;
			while ((k < interid) && (allts[threadid][k].second.first < 0)) {
				k++;
			}
			t = allts[threadid][k].second.first;
			best_index = k;
			k++;
			while ((k < interid) && (allts[threadid][k].first <= t)) {
				t = std::max(t, allts[threadid][k].second.first);
				best_index = k;
				k++;
			}
		}
		triangle_id = best_index;

		P = d.origin + t * d.direction;
		mat.shadingN = P - toVecf(particles[frame][best_index]); mat.shadingN.normalize();

		//if (dot(mat.shadingN, d.direction) > 0 && !mat.transp) mat.shadingN = -mat.shadingN;
		if (flip_normals) mat.shadingN = -mat.shadingN;

		//mat.Kd = Vector(0.5, 0.5, 0.5);
		//mat.Kd = Vector(localN[0]*0.5+0.5, localN[1] * 0.5 + 0.5, localN[2] * 0.5 + 0.5);
		//mat.Kd = Vector(std::abs(localN[0]), std::abs(localN[1]), std::abs(localN[2]));
	}


	return has_inter;
}



bool Fluid::intersection_opaque2(const Ray& d, Vector& P, float &t, MaterialValues &mat, float cur_best_t, int &triangle_id) const {
	t = cur_best_t;
	bool has_inter = false;
	float t_box_left, t_box_right;
	int best_index = -1;
	bool goleft, goright;
	Vector localP, localN;
	double localt;
	int frame = d.time;

	Ray invd(d.origin, Vector(1.f / d.direction[0], 1.f / d.direction[1], 1.f / d.direction[2]), d.time);
	char signs[3];
	signs[0] = (invd.direction[0] >= 0) ? 1 : 0;
	signs[1] = (invd.direction[1] >= 0) ? 1 : 0;
	signs[2] = (invd.direction[2] >= 0) ? 1 : 0;

	if (!extent.intersection_invd(invd, signs, t_box_left)) return false;

	int threadid = omp_get_thread_num();
	int interid = 0;
	Vectord startingP = toVecd(d.origin) + (t_box_left + 1E-9) * toVecd(d.direction);
	int boxid[3] = { static_cast<int>(std::round((startingP[0] - extent.bounds[0][0]) / dx[0])), static_cast<int>(std::round((startingP[1] - extent.bounds[0][1]) / dx[1])) , static_cast<int>(std::round((startingP[2] - extent.bounds[0][2]) / dx[2])) };
	if (boxid[0] < 0 || boxid[1] < 0 || boxid[2] < 0 || boxid[0] >= Nx || boxid[1] >= Ny || boxid[2] >= Nz) return false;

	char stepX = signs[0] * 2 - 1;
	char stepY = signs[1] * 2 - 1;
	char stepZ = signs[2] * 2 - 1;
	Vectord step(stepX, stepY, stepZ);

	Vectord tMax = (((Vectord(boxid[0], boxid[1], boxid[2]) + Vectord(signs[0], signs[1], signs[2]))*dx + extent.bounds[0]) - startingP) / toVecd(d.direction); // next ray-plane intersection
	Vectord tDelta = (step * dx) / toVecd(d.direction);

	t = 1E30;	
	while (true) {
		int voxid = boxid[2] * Nx*Ny + boxid[1] * Nx + boxid[0];
		if (accelgrid[voxid]) {
			int nids = accelgridindices[voxid].size();
			const int* ids = &accelgridindices[voxid][0];
			for (int i = 0; i < nids; i++) {
				float t1;
				Vector curP, curN;
				if (SimplerSphere(particles[frame][ids[i]], radius).intersection(d, curP, curN, t1)) {					
					if (t1 < t) {
						has_inter = true;
						t = t1;
						P = curP;
						mat.shadingN = curN;
						triangle_id = ids[i];
					}					
				}
			}
			if (has_inter) break;
		}
		if (tMax[0] < tMax[1]) {
			if (tMax[0] < tMax[2]) {
				boxid[0] += step[0];
				if (boxid[0] >= Nx || boxid[0] < 0) break;
				tMax[0] += tDelta[0];
			} else {
				boxid[2] += step[2];
				if (boxid[2] >= Nz || boxid[2] < 0) break;
				tMax[2] += tDelta[2];
			}
		} else {
			if (tMax[1] < tMax[2]) {
				boxid[1] += step[1];
				if (boxid[1] >= Ny || boxid[1] < 0) break;
				tMax[1] += tDelta[1];
			} else {
				boxid[2] += step[2];
				if (boxid[2] >= Nz || boxid[2] < 0) break;
				tMax[2] += tDelta[2];
			}
		}
	}

	if (has_inter) {	
		//if (dot(mat.shadingN, d.direction) > 0 && !mat.transp) mat.shadingN = -mat.shadingN;
		if (flip_normals) mat.shadingN = -mat.shadingN;
		//mat.Kd = Vector(0.5, 0.5, 0.5);
		//mat.Kd = Vector(localN[0]*0.5+0.5, localN[1] * 0.5 + 0.5, localN[2] * 0.5 + 0.5);
		//mat.Kd = Vector(std::abs(localN[0]), std::abs(localN[1]), std::abs(localN[2]));
		mat.Kd = toVecf(visualparticlescolor[triangle_id]);
	}


	return has_inter;
}

bool Fluid::intersection_shadow2(const Ray& d, float &t, float cur_best_t, float dist_light) const {
	t = cur_best_t;
	float t_box_left, t_box_right;
	int best_index = -1;
	bool goleft, goright;
	Vector localP, localN;
	float localt;
	int frame = d.time;

	Ray invd(d.origin, Vector(1. / d.direction[0], 1. / d.direction[1], 1. / d.direction[2]), d.time);
	char signs[3];
	signs[0] = (invd.direction[0] >= 0) ? 1 : 0;
	signs[1] = (invd.direction[1] >= 0) ? 1 : 0;
	signs[2] = (invd.direction[2] >= 0) ? 1 : 0;

	if (!extent.intersection_invd(invd, signs, t_box_left)) return false;

	int threadid = omp_get_thread_num();
	int interid = 0;
	Vectord startingP = toVecd(d.origin) + (t_box_left + 1E-9) * toVecd(d.direction);
	int boxid[3] = { static_cast<int>(std::round((startingP[0] - extent.bounds[0][0]) / dx[0])), static_cast<int>(std::round((startingP[1] - extent.bounds[0][1]) / dx[1])) , static_cast<int>(std::round((startingP[2] - extent.bounds[0][2]) / dx[2])) };
	if (boxid[0] < 0 || boxid[1] < 0 || boxid[2] < 0 || boxid[0] >= Nx || boxid[1] >= Ny || boxid[2] >= Nz) return false;
	char stepX = signs[0] * 2 - 1;
	char stepY = signs[1] * 2 - 1;
	char stepZ = signs[2] * 2 - 1;
	Vectord step(stepX, stepY, stepZ);

	Vectord tMax = (((Vectord(boxid[0], boxid[1], boxid[2]) + Vectord(signs[0], signs[1], signs[2]))*dx + extent.bounds[0]) - startingP) / toVecd(d.direction); // next ray-plane intersection
	Vectord tDelta = (step * dx) / toVecd(d.direction);

	t = 1E30;
	while (true) {
		int voxid = boxid[2] * Nx*Ny + boxid[1] * Nx + boxid[0];
		if (accelgrid[voxid]) {
			int nids = accelgridindices[voxid].size();
			const int* ids = &accelgridindices[voxid][0];
			for (int i = 0; i < nids; i++) {
				float t1;
				Vector curP, curN;
				if (SimplerSphere(particles[frame][ids[i]], radius).intersection_shadow(d, t1)) {
					if (t1 < dist_light) {
						t = t1;
						return true;
					}
				}
			}
		}
		if (tMax[0] < tMax[1]) {
			if (tMax[0] < tMax[2]) {
				boxid[0] += step[0];
				if (boxid[0] >= Nx || boxid[0] < 0) break;
				tMax[0] += tDelta[0];
			} else {
				boxid[2] += step[2];
				if (boxid[2] >= Nz || boxid[2] < 0) break;
				tMax[2] += tDelta[2];
			}
		} else {
			if (tMax[1] < tMax[2]) {
				boxid[1] += step[1];
				if (boxid[1] >= Ny || boxid[1] < 0) break;
				tMax[1] += tDelta[1];
			} else {
				boxid[2] += step[2];
				if (boxid[2] >= Nz || boxid[2] < 0) break;
				tMax[2] += tDelta[2];
			}
		}
	}


	return false;
}

void Fluid::build_bvh_recur(int frame, BVHd* b, int node, int i0, int i1) {

	BVHNodesd n;
	//n.i0 = i0;
	//n.i1 = i1;
	n.bbox = build_bbox(frame, i0, i1);
	n.fg = i0;
	n.fd = i1;
	n.isleaf = true;
	b->nodes.push_back(n);

	BBoxd centerBB = build_centers_bbox(frame, i0, i1);

	Vectord diag = centerBB.bounds[1] - centerBB.bounds[0];
	int split_dim;
	if ((diag[0] >= diag[1]) && (diag[0] >= diag[2])) {
		split_dim = 0;
	} else {
		if ((diag[1] >= diag[0]) && (diag[1] >= diag[2])) {
			split_dim = 1;
		} else {
			split_dim = 2;
		}
	}


	double best_split_factor = 0.5;
	double best_area_bb = 1E50;
#ifdef _DEBUG
	int max_tests = 1;
#else
	int max_tests = 3;
#endif

	for (int test_split = 0; test_split < max_tests; test_split++) {

		double cur_split_factor = (test_split + 1) / (double)(max_tests + 1);
		double split_val = centerBB.bounds[0][split_dim] + diag[split_dim] * cur_split_factor;
		BBoxd bb_left(Vectord(1E10, 1E10, 1E10), Vectord(-1E10, -1E10, -1E10)), bb_right(Vectord(1E10, 1E10, 1E10), Vectord(-1E10, -1E10, -1E10));
		int nl = 0, nr = 0;

		for (int i = i0; i < i1; i++) {
			double center_split_dim = particles[frame][i][split_dim];

			if (center_split_dim <= split_val) {
				bb_left.bounds[0] = min(bb_left.bounds[0], particles[frame][i] - Vectord(radius, radius, radius));
				bb_left.bounds[1] = max(bb_left.bounds[1], particles[frame][i] + Vectord(radius, radius, radius));
				nl++;
			} else {
				bb_right.bounds[0] = min(bb_right.bounds[0], particles[frame][i] - Vectord(radius, radius, radius));
				bb_right.bounds[1] = max(bb_right.bounds[1], particles[frame][i] + Vectord(radius, radius, radius));
				nr++;
			}
		}
		double sum_area_bb = bb_left.area()*nl + bb_right.area()*nr;
		if (sum_area_bb < best_area_bb) {
			best_split_factor = cur_split_factor;
			best_area_bb = sum_area_bb;
		}
	}

	double split_val = centerBB.bounds[0][split_dim] + diag[split_dim] * best_split_factor;
	int pivot = i0 - 1;
	for (int i = i0; i < i1; i++) {
		double center_split_dim = particles[frame][i][split_dim];

		if (center_split_dim <= split_val) {
			pivot++;
			std::swap(particles[frame][i], particles[frame][pivot]);
			/*std::swap(colors[i], colors[pivot]);
			std::swap(normals[i], normals[pivot]);
			std::swap(radius[i], radius[pivot]);*/
		}
	}


	if (pivot < i0 || pivot >= i1 - 1 || i1 <= i0 + 1 /**/) { // 1 spheres per leaf
		return;
	}

	b->nodes[node].isleaf = false;
	b->nodes[node].fg = b->nodes.size();
	build_bvh_recur(frame, b, b->nodes[node].fg, i0, pivot + 1);

	b->nodes[node].fd = b->nodes.size();
	build_bvh_recur(frame, b, b->nodes[node].fd, pivot + 1, i1);

}


bool Fluid::intersection_opaque(const Ray& d, Vector& P, float &t, MaterialValues &mat, float cur_best_t, int &triangle_id) const {
	t = cur_best_t;
	bool has_inter = false;
	float t_box_left, t_box_right;
	int best_index = -1;
	bool goleft, goright;
	Vector localP, localN;
	float localt;
	int frame = d.time;

	Ray invd(d.origin, Vector(1. / d.direction[0], 1. / d.direction[1], 1. / d.direction[2]), d.time);
	char signs[3];
	signs[0] = (invd.direction[0] >= 0) ? 1 : 0;
	signs[1] = (invd.direction[1] >= 0) ? 1 : 0;
	signs[2] = (invd.direction[2] >= 0) ? 1 : 0;

	if (!bvh.bbox.intersection_invd(invd, signs, t_box_left)) return false;
	if (t_box_left > cur_best_t) return false;

	int tid = omp_get_thread_num();
	int* l = nodelist[tid];
	float* tnear = tnearlist[tid];
	int idx_back = -1;

	l[++idx_back] = 0;
	tnear[idx_back] = t_box_left;

	while (idx_back >= 0) {

		if (tnear[idx_back] > t) {
			idx_back--;
			continue;
		}
		const int current = l[idx_back--];

		const int fg = bvh.nodes[current].fg;
		const int fd = bvh.nodes[current].fd;

		if (!bvh.nodes[current].isleaf) {
			goleft = (bvh.nodes[fg].bbox.intersection_invd(invd, signs, t_box_left) && t_box_left < t);
			goright = (bvh.nodes[fd].bbox.intersection_invd(invd, signs, t_box_right) && t_box_right < t);

			if (goleft&&goright) {
				if (t_box_left < t_box_right) {
					l[++idx_back] = fd;  tnear[idx_back] = t_box_right;
					l[++idx_back] = fg;  tnear[idx_back] = t_box_left;
				} else {
					l[++idx_back] = fg;  tnear[idx_back] = t_box_left;
					l[++idx_back] = fd;  tnear[idx_back] = t_box_right;
				}
			} else {
				if (goleft) { l[++idx_back] = fg; tnear[idx_back] = t_box_left; }
				if (goright) { l[++idx_back] = fd; tnear[idx_back] = t_box_right; }
			}
		} else {  // feuille

			for (int i = fg; i < fd; i++) {
				if (SimplerSphere(particles[frame][i], radius).intersection(d, localP, localN, localt)) {
					if (localt < t) {
						has_inter = true;
						best_index = i;
						t = localt;
					}
				}
			}
		}

	}

	if (has_inter) {
		int i = best_index;
		triangle_id = best_index;
		SimplerSphere(particles[frame][i], radius).intersection(d, localP, localN, localt);
		P = localP;

		queryMaterial(0, 0, 0, mat);
		mat.shadingN = localN;

		//if (dot(mat.shadingN, d.direction) > 0 && !mat.transp) mat.shadingN = -mat.shadingN;
		if (flip_normals) mat.shadingN = -mat.shadingN;

		mat.Kd = Vector(0.5, 0.5, 0.5);
		//mat.Kd = Vector(localN[0]*0.5+0.5, localN[1] * 0.5 + 0.5, localN[2] * 0.5 + 0.5);
		//mat.Kd = Vector(std::abs(localN[0]), std::abs(localN[1]), std::abs(localN[2]));
	}

	return has_inter;
}



bool Fluid::intersection_transparent(const Ray& d, Vector& P, float &t, MaterialValues &mat, float cur_best_t, int &triangle_id) const {
	t = cur_best_t;
	bool has_inter = false;
	float t_box_left, t_box_right;
	int best_index = -1;
	bool goleft, goright;
	Vector localP, localN;
	float localt;
	int frame = d.time;

	Ray invd(d.origin, Vector(1. / d.direction[0], 1. / d.direction[1], 1. / d.direction[2]), d.time);
	char signs[3];
	signs[0] = (invd.direction[0] >= 0) ? 1 : 0;
	signs[1] = (invd.direction[1] >= 0) ? 1 : 0;
	signs[2] = (invd.direction[2] >= 0) ? 1 : 0;

	if (!bvh.bbox.intersection_invd(invd, signs, t_box_left)) return false;
	if (t_box_left > cur_best_t) return false;
	if (t_box_left > 0) return intersection_opaque(d, P, t, mat, cur_best_t, triangle_id); // ray coming from outside the global bbox : same thing as opaque intersections.

	int threadid = omp_get_thread_num();
	int* l = nodelist[threadid];
	float* tnear = tnearlist[threadid];
	int idx_back = -1;

	l[++idx_back] = 0;
	tnear[idx_back] = t_box_left;


	int interid = 0;
	while (idx_back >= 0) {

		if (tnear[idx_back] > cur_best_t) {
			idx_back--;
			continue;
		}

		const int current = l[idx_back--];

		const int fg = bvh.nodes[current].fg;
		const int fd = bvh.nodes[current].fd;

		if (!bvh.nodes[current].isleaf) {
			if (signs[0] == 1) {
				goleft = bvh.nodes[fg].bbox.intersection_invd_positive_x(invd, signs, t_box_left);
				goright = bvh.nodes[fd].bbox.intersection_invd_positive_x(invd, signs, t_box_right);
			} else {
				goleft = bvh.nodes[fg].bbox.intersection_invd_negative_x(invd, signs, t_box_left);
				goright = bvh.nodes[fd].bbox.intersection_invd_negative_x(invd, signs, t_box_right);
			}

			if (goleft&&goright) {				
				l[++idx_back] = fd;  tnear[idx_back] = t_box_right;
				l[++idx_back] = fg;  tnear[idx_back] = t_box_left;
			} else {
				if (goleft) { l[++idx_back] = fg; tnear[idx_back] = t_box_left; }
				if (goright) { l[++idx_back] = fd; tnear[idx_back] = t_box_right; }
			}
		} else {  // feuille

			for (int i = fg; i < fd; i++) {
				float t1, t2;
				if (SimplerSphere(particles[frame][i], radius).both_intersections(d, t1, t2)) {					
					has_inter = true;					
					allts[threadid][interid] = std::make_pair(t1, std::make_pair(t2, i));
					interid++;
				}
			}
		}

	}

	if (has_inter) {
		std::sort(allts[threadid], allts[threadid]+ interid+1);
		if (allts[threadid][0].first > 0) {
			best_index = allts[threadid][0].second.second;
			t = allts[threadid][0].first;
		} else {
			int k = 0;
			while ((k < interid) && (allts[threadid][k].second.first < 0)) {
				k++;
			}
			t = allts[threadid][k].second.first;
			best_index = k;		
			k++;
			while ((k < interid) && (allts[threadid][k].first <= t)) {
				t = std::max(t, allts[threadid][k].second.first);
				best_index = k;
				k++;
			}
		}
		triangle_id = best_index;
		
		P = d.origin + t * d.direction;
		mat.shadingN = P - toVecf(particles[frame][best_index]); mat.shadingN.normalize();

		//if (dot(mat.shadingN, d.direction) > 0 && !mat.transp) mat.shadingN = -mat.shadingN;
		if (flip_normals) mat.shadingN = -mat.shadingN;

		//mat.Kd = Vector(0.5, 0.5, 0.5);
		//mat.Kd = Vector(localN[0]*0.5+0.5, localN[1] * 0.5 + 0.5, localN[2] * 0.5 + 0.5);
		//mat.Kd = Vector(std::abs(localN[0]), std::abs(localN[1]), std::abs(localN[2]));
	}

	
	return has_inter;
}


bool Fluid::intersection(const Ray& d, Vector& P, float &t, MaterialValues &mat, float cur_best_t, int &triangle_id) const {
	queryMaterial(0, 0, 0, mat);
	if (mat.transp) {
		return intersection_transparent2(d, P, t, mat, cur_best_t, triangle_id);
	} else {
		return intersection_opaque2(d, P, t, mat, cur_best_t, triangle_id);
	}
}

bool Fluid::reservoir_sampling_intersection(const Ray& d, Vector& P, float &t, MaterialValues &mat, int &triangle_id, int &current_nb_intersections, float min_t, float max_t) const {
	
	// NOT IMPLEMENTED (correctly).
	float localt;
	Vector localP;
	int local_tri_id;
	MaterialValues localmat;
	if (intersection(d, localP, localt, localmat, max_t, local_tri_id)) {
		if (t > min_t && t > max_t) {
			int threadid = omp_get_thread_num();
			const float invmax = 1.f / engine[threadid].max();
			float r1 = engine[threadid]()*invmax;
			current_nb_intersections++;
			if (r1 < 1. / current_nb_intersections) {	
				P = localP;
				t = localt;
				triangle_id = local_tri_id;
				mat = localmat;
				return true;
			}
		}
	}
	return false;

}


bool Fluid::intersection_shadow(const Ray& d, float &t, float cur_best_t, float dist_light) const { // not correct, but no big deal

	return intersection_shadow2(d, t, cur_best_t, dist_light);

	t = cur_best_t;
	bool has_inter = false;
	float t_box_left, t_box_right;
	int best_index = -1;
	bool goleft, goright;
	Vector localP, localN;
	float localt;
	int frame = d.time;

	Ray invd(d.origin, Vector(1. / d.direction[0], 1. / d.direction[1], 1. / d.direction[2]), d.time);
	char signs[3];
	signs[0] = (invd.direction[0] >= 0) ? 1 : 0;
	signs[1] = (invd.direction[1] >= 0) ? 1 : 0;
	signs[2] = (invd.direction[2] >= 0) ? 1 : 0;

	if (!bvh.bbox.intersection_invd(invd, signs, t_box_left)) return false;
	if (t_box_left > cur_best_t || t_box_left > dist_light) return false;

	int tid = omp_get_thread_num();
	int* l = nodelist[tid];
	float* tnear = tnearlist[tid];
	int idx_back = -1;

	l[++idx_back] = 0;
	tnear[idx_back] = t_box_left;

	while (idx_back >= 0) {

		if (tnear[idx_back] > t) {
			idx_back--;
			continue;
		}
		const int current = l[idx_back--];

		const int fg = bvh.nodes[current].fg;
		const int fd = bvh.nodes[current].fd;

		if (!bvh.nodes[current].isleaf) {
			goleft = (bvh.nodes[fg].bbox.intersection_invd(invd, signs, t_box_left) && (t_box_left < t) && (t_box_left < dist_light));
			goright = (bvh.nodes[fd].bbox.intersection_invd(invd, signs, t_box_right) && (t_box_right < t) && (t_box_right < dist_light));

			if (goleft&&goright) {
				if (t_box_left < t_box_right) {
					l[++idx_back] = fd;  tnear[idx_back] = t_box_right;
					l[++idx_back] = fg;  tnear[idx_back] = t_box_left;
				} else {
					l[++idx_back] = fg;  tnear[idx_back] = t_box_left;
					l[++idx_back] = fd;  tnear[idx_back] = t_box_right;
				}
			} else {
				if (goleft) { l[++idx_back] = fg; tnear[idx_back] = t_box_left; }
				if (goright) { l[++idx_back] = fd; tnear[idx_back] = t_box_right; }
			}
		} else {  // feuille

			for (int i = fg; i < fd; i++) {
				if (SimplerSphere(particles[frame][i], radius).intersection_shadow(d, localt)) {
					if (localt < t) {
						has_inter = true;
						best_index = i;
						t = localt;
					}
				}
			}
		}

	}

	return has_inter;
}
