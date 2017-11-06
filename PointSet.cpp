#include "PointSet.h"


BBox PointSet::build_bbox(int i0, int i1) {

	BBox result;
	result.bounds[1] = vertices[i0] + Vector(radius, radius, radius);
	result.bounds[0] = vertices[i0] - Vector(radius, radius, radius);
	for (int i = i0; i < i1; i++) { // indice de triangle
		result.bounds[0] = min(result.bounds[0], vertices[i] - Vector(radius, radius, radius));
		result.bounds[1] = max(result.bounds[1], vertices[i] + Vector(radius, radius, radius));
	}
	return result;
}

BBox PointSet::build_centers_bbox(int i0, int i1) {

	BBox result;
	result.bounds[1] = vertices[i0];
	result.bounds[0] = vertices[i0];
	for (int i = i0; i < i1; i++) { // indice de triangle
		result.bounds[0] = min(result.bounds[0], vertices[i]);
		result.bounds[1] = max(result.bounds[1], vertices[i]);
	}
	return result;
}

void PointSet::build_bvh(int i0, int i1) {
	bvh.bbox = build_bbox(i0, i1);
	bvh.nodes.reserve(vertices.size() * 2);
	build_bvh_recur(&bvh, 0, i0, i1);
}

void PointSet::build_bvh_recur(BVH* b, int node, int i0, int i1) {

	BVHNodes n;
	//n.i0 = i0;
	//n.i1 = i1;
	n.bbox = build_bbox(i0, i1);
	n.fg = i0;
	n.fd = i1;
	n.isleaf = true;
	b->nodes.push_back(n);

	BBox centerBB = build_centers_bbox(i0, i1);

	Vector diag = centerBB.bounds[1] - centerBB.bounds[0];
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
	int max_tests = 16;
#endif

	for (int test_split = 0; test_split < max_tests; test_split++) {

		double cur_split_factor = (test_split + 1) / (double)(max_tests + 1);
		double split_val = centerBB.bounds[0][split_dim] + diag[split_dim] * cur_split_factor;
		BBox bb_left(Vector(1E10, 1E10, 1E10), Vector(-1E10, -1E10, -1E10)), bb_right(Vector(1E10, 1E10, 1E10), Vector(-1E10, -1E10, -1E10));
		int nl = 0, nr = 0;

		for (int i = i0; i < i1; i++) {
			double center_split_dim = vertices[i][split_dim];

			if (center_split_dim <= split_val) {
				bb_left.bounds[0] = min(bb_left.bounds[0], vertices[i] - Vector(radius, radius, radius));
				bb_left.bounds[1] = max(bb_left.bounds[1], vertices[i] + Vector(radius, radius, radius));
				nl++;
			} else {
				bb_right.bounds[0] = min(bb_right.bounds[0], vertices[i] - Vector(radius, radius, radius));
				bb_right.bounds[1] = max(bb_right.bounds[1], vertices[i] + Vector(radius, radius, radius));
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
		double center_split_dim = vertices[i][split_dim];

		if (center_split_dim <= split_val) {
			pivot++;
			std::swap(vertices[i], vertices[pivot]);
			std::swap(colors[i], colors[pivot]);
			std::swap(normals[i], normals[pivot]);
		}
	}


	if (pivot < i0 || pivot >= i1 - 1 || i1 <= i0 + 4) { // 4 triangles per leaf
		return;
	}

	b->nodes[node].isleaf = false;
	b->nodes[node].fg = b->nodes.size();
	build_bvh_recur(b, b->nodes[node].fg, i0, pivot + 1);

	b->nodes[node].fd = b->nodes.size();
	build_bvh_recur(b, b->nodes[node].fd, pivot + 1, i1);

}


bool PointSet::intersection(const Ray& d, Vector& P, double &t, MaterialValues &mat, double cur_best_t, int &triangle_id) const
{
	t = cur_best_t;
	bool has_inter = false;
	double t_box_left, t_box_right;
	int best_index = -1;
	bool goleft, goright;
	Vector localP, localN;
	double localt;
	double alpha, beta, gamma;

	Ray invd(d.origin, Vector(1. / d.direction[0], 1. / d.direction[1], 1. / d.direction[2]), d.time);
	char signs[3];
	signs[0] = (invd.direction[0] >= 0) ? 1 : 0;
	signs[1] = (invd.direction[1] >= 0) ? 1 : 0;
	signs[2] = (invd.direction[2] >= 0) ? 1 : 0;

	if (!bvh.bbox.intersection_invd(invd, signs, t_box_left)) return false;
	if (t_box_left > cur_best_t) return false;

	int l[50];
	double tnear[50];
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
				if (Disk(vertices[i], normals[i], radius).intersection(d, localP, localN, localt)) {
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
		Disk(vertices[i], normals[i], radius).intersection(d, localP, localN, localt);
		localN.normalize();
		P = localP;
		mat.shadingN = localN;

		if (dot(mat.shadingN, d.direction) > 0 && !transparent) mat.shadingN = -mat.shadingN;
		if (flip_normals) mat.shadingN = -mat.shadingN;

		if (colors.size() > i)
			mat.Kd = colors[i];
		else
			mat.Kd = Vector(0.5, 0.5, 0.5);
		if (display_edges) {
			double r2 = (localP - vertices[i]).getNorm2();
			if (r2 > (radius*radius*0.95*0.95))
				mat.Kd = Vector(0., 0., 0.);
		}


	}



	/*if (!bb.intersection(d)) return false;

	t = 1E99;
	bool has_inter = false;
	for (int i = 0; i < faces.size() / 3; i++) {
	int i0 = faces[i * 3];
	int i1 = faces[i * 3+1];
	int i2 = faces[i * 3+2];
	Triangle tri(vertices[i0], vertices[i1], vertices[i2], albedo, miroir, transparent);
	Vector localP, localN;
	double localt;
	if (tri.intersection(d, localP, localN, localt)) {
	has_inter = true;
	if (localt < t) {
	t = localt;
	P = localP;
	N = localN;
	}
	}
	}*/

	return has_inter;
}


bool PointSet::intersection_shadow(const Ray& d, double &t, double cur_best_t, double dist_light) const
{
	t = cur_best_t;
	bool has_inter = false;
	double t_box_left, t_box_right;
	int best_index = -1;
	bool goleft, goright;
	Vector localP, localN;
	double localt;
	double alpha, beta, gamma;

	Ray invd(d.origin, Vector(1. / d.direction[0], 1. / d.direction[1], 1. / d.direction[2]), d.time);
	char signs[3];
	signs[0] = (invd.direction[0] >= 0) ? 1 : 0;
	signs[1] = (invd.direction[1] >= 0) ? 1 : 0;
	signs[2] = (invd.direction[2] >= 0) ? 1 : 0;

	if (!bvh.bbox.intersection_invd(invd, signs, t_box_left)) return false;
	if (t_box_left > cur_best_t || t_box_left > dist_light) return false;

	int l[50];
	double tnear[50];
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
				if (Disk(vertices[i], normals[i], radius).intersection(d, localP, localN, localt)) {
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
