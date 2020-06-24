#pragma once

#define _CRT_SECURE_NO_WARNINGS
#include "Geometry.h"
#include "Vector.h"

#include "TriangleMesh.h"

#include <random>

class SimplerSphere {
public:
	SimplerSphere() {};
	SimplerSphere(const Vector &origin, double rayon) {
		init(origin, rayon);
	};

	void init(const Vector &origin, double rayon) {
		O = origin;
		R2 = rayon* rayon;
	}
	bool intersection(const Ray& d, Vector& P, Vector &N, double &t) const {
		// resout a*t^2 + b*t + c = 0
		double a = 1; // d.direction.getNorm2();
		double b = 2 * dot(d.direction, d.origin - O);
		double c = (d.origin - O).getNorm2() - R2;

		double delta = b * b - 4 /** a*/ *c;
		if (delta < 0) return false;

		double sqDelta = sqrt(delta);
		double inv2a = 0.5; // 1. / (2. * a);
		double t2 = (-b + sqDelta) * inv2a;

		if (t2 < 0) return false;

		double t1 = (-b - sqDelta) * inv2a;

		if (t1 > 0)
			t = t1;
		else
			t = t2;

		P = d.origin + t * d.direction;
		N = (P - O).getNormalized();
		return true;
	}

	bool both_intersections(const Ray& d, double &t1, double&t2) const {
		// solves t^2 - 2*b*t + c = 0
		double b = -dot(d.direction, d.origin - O);
		double c = (d.origin - O).getNorm2() - R2;

		double delta = b * b - c;
		if (delta < 0) return false;

		double sqDelta = sqrt(delta);		
		t2 = b + sqDelta;
		if (t2 < 0) return false;

		t1 = b - sqDelta;
		//if (t1 < 0) t1 = t2;
		return true;
	}

	bool intersection_shadow(const Ray& d, double &t) const {
		double b = -dot(d.direction, d.origin - O);
		double c = (d.origin - O).getNorm2() - R2;

		double delta = b * b - c;
		if (delta < 0) return false;

		double sqDelta = sqrt(delta);	
		double t2 = b + sqDelta;

		if (t2 < 0) return false;
		double t1 = b - sqDelta;

		if (t1 > 0)
			t = t1;
		else
			t = t2;
		return true;
	}
	Vector O;
	double R2;
};

class Fluid : public Object {	
public:
	Fluid(Scene& s, BBox extent, int Nx, int Ny, int Nz, int Nparticles, double rho, double sphere_radius, int nbframes, double dt, int nsubsteps) : scene(s), extent(extent), Nx(Nx), Ny(Ny), Nz(Nz), Nparticles(Nparticles), rho(rho), radius(sphere_radius), nbframes(nbframes), dt(dt/nsubsteps), nsubsteps(nsubsteps){
		velX.resize((Nx + 1)*Ny*Nz, 0.);
		velY.resize(Nx*(Ny + 1)*Nz, 0.);
		velZ.resize(Nx*Ny*(Nz + 1), 0.);
		pressure.resize(Nx*Ny*Nz, 0.);
		celltypes.resize(Nx*Ny*Nz, 0);
		this->brdf = new LambertBRDF(Vector(0.5, 0.6, 0.9));
		opaque = true;

		newVelX.resize((Nx + 1)*Ny*Nz, 0.);
		newVelY.resize(Nx*(Ny + 1)*Nz, 0.);
		newVelZ.resize(Nx*Ny*(Nz + 1), 0.);
		
		Nv = Vector(Nx, Ny, Nz);
		dx = (extent.bounds[1] - extent.bounds[0]) / Vector(Nx, Ny, Nz);
		dx2 = dx * dx;
		invdx2 = Vector(1.,1.,1.) / dx2;
		invdx = Vector(1., 1., 1.) / dx;
		name = "Fluid1";


		rasterize();
		

	}

	void rasterize() {

#pragma omp parallel for schedule(static)
		for (int i = 0; i < Nz; i++)
			for (int j = 0; j < Ny; j++)
				for (int k = 0; k < Nx; k++) {
					Vector center = extent.bounds[0] + Vector(k+0.5, j+0.5, i+0.5)*dx;
					Vector dir(0.5, 0., 0.5); dir.normalize();
					Vector P;
					int id, tid;
					double mint;
					MaterialValues mat0, mat1;
					bool i1 = scene.intersection(Ray(center, dir, 0), P, id, mint, mat0, tid);
					bool i2 = scene.intersection(Ray(center, -dir, 0), P, id, mint, mat1, tid);

					if (i1&&i2&&dot(mat0.shadingN, dir) > 0 && dot(mat1.shadingN, dir) < 0) {
						celltypes[i*Ny*Nx + j * Nx + k] = 2;
					}
				}
	}


	void jfa(int k, int Nx, int Ny, int Nz, int* sites) {

		int neigh[27][3] = { {-k, -k, -k}, {-k, -k, 0}, {-k, -k, k}, {-k, 0, -k}, {-k, 0, 0}, {-k, 0, k}, {-k, k, -k}, {-k, k, 0}, {-k, k, k},
		                     {0, -k, -k}, {0, -k, 0}, {0, -k, k}, {0, 0, -k},  {0, 0, 0}, {0, 0, k}, {0, k, -k}, {0, k, 0}, {0, k, k},
		                     {k, -k, -k}, {k, -k, 0}, {k, -k, k}, {k, 0, -k}, {k, 0, 0}, {k, 0, k}, {k, k, -k}, {k, k, 0}, {k, k, k} };

		std::vector<int> newsites(Nx*Ny*Nz);		
#pragma omp parallel for
		for (int i = 0; i < Nz; i++) {
			for (int j = 0; j < Ny; j++) {
				for (int k = 0; k < Nx; k++) {
					int bestdist = std::numeric_limits<int>::max();
					int bestneigh = -1;
					for (int n = 0; n < 27; n++) {
						int ni = i + neigh[n][0]; if (ni < 0) continue; if (ni >= Nz) continue;
						int nj = j + neigh[n][1]; if (nj < 0) continue; if (nj >= Ny) continue;
						int nk = k + neigh[n][2]; if (nk < 0) continue; if (nk >= Nx) continue;
						int neigh = sites[ni*Nx*Ny + nj * Nx + nk];
						if (neigh < 0) continue;
						int neighi = neigh / (Ny*Nx);
						int neighj = (neigh / Nx)%Ny;
						int neighk = neigh % Nx;
						int distneigh = (neighi - i)*(neighi - i) + (neighj - j)*(neighj - j) + (neighk - k)*(neighk - k);
						if (distneigh < bestdist) {
							bestdist = distneigh;
							bestneigh = neigh;
						}
					}
					newsites[i*Nx*Ny + j * Nx + k] = bestneigh;
				}
			}
		}

		memcpy(sites, &newsites[0], Nx*Ny*Nz * sizeof(int));
	}

	void jfa(int Nx, int Ny, int Nz, int* sites) {
		int N = std::max(Nx, std::max(Ny, Nz));
		for (int i = N / 2; i >= 1; i /= 2) {
			jfa(i, Nx, Ny, Nz, sites);
		}
		jfa(1, Nx, Ny, Nz, sites);
	}

	void NNextrapolate(int Nx, int Ny, int Nz, int* sites, double* vel) {
		jfa(Nx, Ny, Nz, sites);
		for (int i = 0; i < Nx*Ny*Nz; i++) {
			vel[i] = vel[sites[i]];
		}
	}

	void extrapolateVel2() {

		std::vector<int> computedVelX((Nx + 1)*Ny*Nz);
		/*for (int i = 0; i < Nz; i++) {
			for (int j = 0; j < Ny; j++) {
				computedVelX[i*(Nx + 1)*Ny + j * (Nx+1) + 0] = i * (Nx + 1)*Ny + j * (Nx+1) + 0;
				computedVelX[i*(Nx + 1)*Ny + j * (Nx+1) + Nx] = i * (Nx + 1)*Ny + j * (Nx+1) + Nx;
			}
		}*/
		for (int i = 0; i < Nz; i++) {
			for (int j = 0; j < Ny; j++) {
				for (int k = 0; k < Nx+1; k++) {
					if (k==0 || k==Nx || celltypes[i*Nx*Ny + j * Nx + k] == 1 || celltypes[i*Nx*Ny + j * Nx + k - 1] == 1) {
						computedVelX[i*(Nx + 1)*Ny + j * (Nx+1) + k] = i * (Nx + 1)*Ny + j * (Nx+1) + k;
					} else {
						computedVelX[i*(Nx + 1)*Ny + j * (Nx+1) + k] = -1;
					}
				}
			}
		}
		
		NNextrapolate(Nx + 1, Ny, Nz, &computedVelX[0], &velX[0]);

		std::vector<int> computedVelY(Nx*(Ny+1)*Nz);
		for (int i = 0; i < Nz; i++) {
			for (int j = 0; j < Ny+1; j++) {
				for (int k = 0; k < Nx; k++) {
					if (j == 0 || j == Ny || celltypes[i*Nx*Ny + j * Nx + k] == 1 || celltypes[i*Nx*Ny + (j-1) * Nx + k] == 1) {
						computedVelY[i*Nx*(Ny+1) + j * Nx + k] = i * Nx*(Ny+1) + j * Nx + k;
					} else {
						computedVelY[i*Nx*(Ny+1) + j * Nx + k] = -1;
					}
				}
			}
		}

		NNextrapolate(Nx, Ny+1, Nz, &computedVelY[0], &velY[0]);

		std::vector<int> computedVelZ(Nx*Ny*(Nz+1));
		for (int i = 0; i < Nz+1; i++) {
			for (int j = 0; j < Ny; j++) {
				for (int k = 0; k < Nx; k++) {
					if (i == 0 || i == Nz || celltypes[i*Nx*Ny + j * Nx + k] == 1 || celltypes[(i-1)*Nx*Ny + j* Nx + k] == 1) {
						computedVelZ[i*Nx*Ny + j * Nx + k] = i * Nx*Ny + j * Nx + k;
					} else {
						computedVelZ[i*Nx*Ny + j * Nx + k] = -1;
					}
				}
			}
		}

		NNextrapolate(Nx, Ny, Nz+1, &computedVelZ[0], &velZ[0]);
	}

	void init_particles(bool initwithshape, int selected_object) {

		std::vector<Vector> initpart;
		std::random_device dev;
		std::mt19937 rng(dev());
		std::uniform_real_distribution<double> uniform(0, 1);

		int nbcells = 0;

		if (initwithshape) {
			cellcol.resize(Nx*Ny*Nz);
			Vector dir(0.5, 0., 0.5); dir.normalize();
			Vector transformed_dir = scene.objects[selected_object]->apply_inverse_rotation_scaling(dir);
			Vector randdirs[5];
			for (int i = 0; i < 5; i++) {
				randdirs[i] = Vector(uniform(rng) - 0.5, uniform(rng) - 0.5, uniform(rng) - 0.5);
				randdirs[i].normalize();
				randdirs[i] = scene.objects[selected_object]->apply_inverse_rotation_scaling(randdirs[i]);
			}


			for (int i = 0; i < Nz; i++)
				for (int j = 0; j < Ny; j++)
					for (int k = 0; k < Nx; k++) {
						Vector center = extent.bounds[0] + Vector(k + 0.5, j + 0.5, i + 0.5)*dx;
						Vector new_origin = scene.objects[selected_object]->apply_inverse_transformation(center);
						Vector P;
						int id, tid;
						double mint1, mint2, mint;
						Vector col(1., 1., 1.);
						MaterialValues mat0, mat1;
						bool i1 = scene.objects[selected_object]->intersection(Ray(new_origin, transformed_dir, 0), P, mint1, mat0, 1E9, tid);
						bool i2 = scene.objects[selected_object]->intersection(Ray(new_origin, -transformed_dir, 0), P, mint2, mat1, 1E9, tid);
						mint = mint1;
						col = mat0.Kd;
						if (mint2 < mint) {
							mint = mint2;
							col = mat1.Kd;
						}

						if (i1&&i2&&dot(scene.objects[selected_object]->apply_rotation(mat0.shadingN), dir) > 0 && dot(scene.objects[selected_object]->apply_rotation(mat1.shadingN), dir) < 0) {
							celltypes[i*Ny*Nx + j * Nx + k] = 1;
							nbcells++;


							for (int testray = 0; testray < 5; testray++) {
								bool i1 = scene.objects[selected_object]->intersection(Ray(new_origin, randdirs[i], 0), P, mint, mat0, 1E9, tid);
								bool i2 = scene.objects[selected_object]->intersection(Ray(new_origin, -randdirs[i], 0), P, mint, mat1, 1E9, tid);
								if (mint1 < mint) {
									mint = mint1;
									col = mat0.Kd;
								}
								if (mint2 < mint) {
									mint = mint2;
									col = mat1.Kd;
								}
							}
							cellcol[i*Ny*Nx + j * Nx + k] = col;
						}
					}

		} else {

			for (int i = (Nz*0.35); i < (Nz*0.65); i++)
				for (int j = (Ny*0.3); j < (Ny*0.8); j++)
					for (int k = (Nx*0.35); k < (Nx*0.65); k++) {
						if (celltypes[i*Ny*Nx + j * Nx + k] == 0) {
							celltypes[i*Ny*Nx + j * Nx + k] = 1;
							nbcells++;
						}
					}
		}
		initpart.reserve(Nparticles);
		visualparticlescolor.reserve(Nparticles);
		double nbpartpercell = Nparticles / (double)nbcells;
		int ipartpercell = std::ceil(nbpartpercell);
		for (int id = 0; id < celltypes.size(); id++) {
			if (celltypes[id] == 1) {
				for (int jd = 0; jd < ipartpercell; jd++) {
					if (uniform(rng) <= nbpartpercell / ipartpercell) {
						int k = id % Nx;
						int j = (id / Nx) % Ny;
						int i = (id / (Nx*Ny)) % Nz;
						initpart.push_back(Vector(k + uniform(rng), j + uniform(rng), i + uniform(rng)) * dx + extent.bounds[0]);
						if(cellcol.size()>0)
							visualparticlescolor.push_back(cellcol[id]);
						else
							visualparticlescolor.push_back(Vector(1.,1.,1.));
					}
				}
			}
		}


		Nparticles = initpart.size();
		particles.resize(nbframes + 1, std::vector<Vector>(Nparticles));
		particles[0] = initpart;

		ghostparticles.reserve(Nx*Ny*Nz * 1);
		for (int i = 0; i < Nz; i++) {
			for (int j = 0; j < Ny; j++) {
				for (int k = 0; k < Nx; k++) {
					if (celltypes[i*Nx*Ny + j * Nx + k] == 1) {
						for (int l = 0; l < 4; l++) { // not jittered
							Vector randP = extent.bounds[0] + Vector(k + uniform(rng), j + uniform(rng), i + uniform(rng))*dx;
							ghostparticles.push_back(randP);
						}
					}
				}
			}
		}
		ghostparticlesnew.resize(ghostparticles.size());


		/*initpart = ghostparticles;
		Nparticles = initpart.size();
		particles[0] = initpart;*/
	}

	template<typename T>
	T interp(const Vector& pos, const std::vector<T> &field, int Nx, int Ny, int Nz) {
		Vector posNormalized = (pos - extent.bounds[0]) / dx; // (extent.bounds[1] - extent.bounds[0])*Vector(Nx, Ny, Nz);
		posNormalized[0] = std::max(0., std::min(Nx - 1.000001, posNormalized[0]));
		posNormalized[1] = std::max(0., std::min(Ny - 1.000001, posNormalized[1]));
		posNormalized[2] = std::max(0., std::min(Nz - 1.000001, posNormalized[2]));
		int icoords[3] = {posNormalized[0], posNormalized[1] , posNormalized[2]};
		//if (icoords[0] >= Nx - 1 || icoords[1] >= Ny - 1 || icoords[2] >= Nz - 1) return T(0);
		//if (icoords[0] < 0 || icoords[1] <0 || icoords[2] <0) return T(0);

		Vector fcoords(posNormalized[0] - icoords[0], posNormalized[1] - icoords[1], posNormalized[2] - icoords[2]);
		/*if (icoords[0] < 0) { icoords[0] = 0; fcoords[0] = 0; }
		if (icoords[0] >= Nx-1) { icoords[0] = Nx-2; fcoords[0] = 1; }
		if (icoords[1] < 0) { icoords[1] = 0; fcoords[1] = 0; }
		if (icoords[1] >= Ny - 1) { icoords[1] = Ny - 2; fcoords[1] = 1; }
		if (icoords[2] < 0) { icoords[2] = 0; fcoords[2] = 0; }
		if (icoords[2] >= Nz - 1) { icoords[2] = Nz - 2; fcoords[2] = 1; }*/


		T result = (((1 - fcoords[0])*field[icoords[2] * Ny*Nx + icoords[1] * Nx + icoords[0]] + fcoords[0] * field[icoords[2] * Ny*Nx + icoords[1] * Nx + (icoords[0] + 1)])*(1 - fcoords[1])
			+ ((1 - fcoords[0])*field[icoords[2] * Ny*Nx + (icoords[1] + 1) * Nx + icoords[0]] + fcoords[0] * field[icoords[2] * Ny*Nx + (icoords[1] + 1) * Nx + (icoords[0] + 1)])*(fcoords[1])) * (1 - fcoords[2])
			+ (((1 - fcoords[0])*field[(icoords[2] + 1) * Ny*Nx + icoords[1] * Nx + icoords[0]] + fcoords[0] * field[(icoords[2] + 1) * Ny*Nx + icoords[1] * Nx + (icoords[0] + 1)])*(1 - fcoords[1])
				+ ((1 - fcoords[0])*field[(icoords[2] + 1) * Ny*Nx + (icoords[1] + 1) * Nx + icoords[0]] + fcoords[0] * field[(icoords[2] + 1) * Ny*Nx + (icoords[1] + 1) * Nx + (icoords[0] + 1)])*(fcoords[1])) * (fcoords[2]);
		return result;
	}
	
	void advect() {	
#pragma omp parallel for schedule(static)
		for (int i = 0; i < Nz; i++) {
			for (int j = 0; j < Ny; j++) {
				for (int k = 1; k < Nx; k++) {

					// get vector coordinates at staggered grid location by bilinear interpolation
					Vector curVelX;
					curVelX[0] = velX[i*Ny*(Nx + 1) + j * (Nx + 1) + k];
					curVelX[1] = (velY[i*(Ny + 1)*Nx + j *Nx + k] + velY[i*(Ny + 1)*Nx + (j+1) *Nx + k] + velY[i*(Ny + 1)*Nx + j * Nx + k-1] + velY[i*(Ny + 1)*Nx + (j+1) *Nx + k-1])*0.25;
					curVelX[2] = (velZ[i*Ny*Nx + j * Nx + k] + velZ[(i + 1)*Ny*Nx + j * Nx + k] + velZ[i*Ny*Nx + j * Nx + k-1] + velZ[(i + 1)*Ny*Nx + j * Nx + k-1])*0.25;

					Vector newVelPosition = Vector(k, j+0.5, i+0.5) * dx + extent.bounds[0] - curVelX*dt; // previous position in world space
					newVelPosition[0] = std::max(std::min(newVelPosition[0], extent.bounds[1][0]), extent.bounds[0][0]);  // stick to walls ; should change
					newVelPosition[1] = std::max(std::min(newVelPosition[1] - 0.5*dx[1], extent.bounds[1][1] - dx[1]), extent.bounds[0][1]);
					newVelPosition[2] = std::max(std::min(newVelPosition[2] - 0.5*dx[2], extent.bounds[1][2] - dx[2]), extent.bounds[0][2]);

					newVelX[i*Ny*(Nx + 1) + j * (Nx + 1) + k] = interp(newVelPosition, velX, Nx+1, Ny, Nz);
				}
			}
		}

#pragma omp parallel for schedule(static)
		for (int i = 0; i < Nz; i++) {
			for (int j = 1; j < Ny; j++) {
				for (int k = 0; k < Nx; k++) {

					// get vector coordinates at staggered grid location by bilinear interpolation
					Vector curVelY;
					curVelY[0] = (velX[i*Ny*(Nx+1) + j * (Nx + 1) + k] + velX[i*Ny*(Nx + 1) + (j - 1) *(Nx + 1) + k] + velX[i*Ny*(Nx + 1) + j * (Nx + 1) + k+1] + velX[i*Ny*(Nx + 1) + (j - 1) *(Nx + 1) + k+1])*0.25;
					curVelY[1] = velY[i*(Ny+1)*Nx + j * Nx + k];
					curVelY[2] = (velZ[i*Ny*Nx + j * Nx + k] + velZ[i*Ny*Nx + (j-1) * Nx + k] + velZ[(i + 1)*Ny*Nx + j * Nx + k] + velZ[(i+1)*Ny*Nx + (j-1) * Nx + k])*0.25;

					Vector newVelPosition = Vector(k+0.5, j, i + 0.5) * dx + extent.bounds[0] - curVelY * dt; // previous position in world space
					newVelPosition[0] = std::max(std::min(newVelPosition[0] - 0.5*dx[0], extent.bounds[1][0] - dx[0]), extent.bounds[0][0]);  // stick to walls ; should change
					newVelPosition[1] = std::max(std::min(newVelPosition[1], extent.bounds[1][1]), extent.bounds[0][1]);
					newVelPosition[2] = std::max(std::min(newVelPosition[2] - 0.5*dx[2], extent.bounds[1][2] - dx[2]), extent.bounds[0][2]);

					newVelY[i*(Ny+1)*Nx + j * Nx + k] = interp(newVelPosition, velY, Nx, Ny+1, Nz);
				}
			}
		}

#pragma omp parallel for schedule(static)
		for (int i = 1; i < Nz; i++) {
			for (int j = 0; j < Ny; j++) {
				for (int k = 0; k < Nx; k++) {

					// get vector coordinates at staggered grid location by bilinear interpolation
					Vector curVelZ;
					curVelZ[0] = (velX[i*Ny*(Nx + 1) + j * (Nx + 1) + k] + velX[(i-1)*Ny*(Nx + 1) + j *(Nx + 1) + k] + velX[(i - 1)*Ny*(Nx + 1) + j * (Nx + 1) + k+1] + velX[i*Ny*(Nx + 1) + j *(Nx + 1) + k+1])*0.25;					
					curVelZ[1] = (velY[i*(Ny+1)*Nx + j * Nx + k] + velY[i*(Ny+1)*Nx + (j + 1) * Nx + k] + velY[(i - 1)*(Ny+1)*Nx + j * Nx + k] + velY[(i - 1)*(Ny+1)*Nx + (j + 1) * Nx + k])*0.25;
					curVelZ[2] = velZ[i*Ny*Nx + j * Nx + k];

					Vector newVelPosition = Vector(k + 0.5, j+0.5, i) * dx + extent.bounds[0] - curVelZ * dt; // previous position in world space
					newVelPosition[0] = std::max(std::min(newVelPosition[0] - 0.5*dx[0], extent.bounds[1][0] - dx[0]), extent.bounds[0][0]);  // stick to walls ; should change
					newVelPosition[1] = std::max(std::min(newVelPosition[1] - 0.5*dx[1], extent.bounds[1][1] - dx[1]), extent.bounds[0][1]);
					newVelPosition[2] = std::max(std::min(newVelPosition[2], extent.bounds[1][2]), extent.bounds[0][2]);

					newVelZ[i*Ny*Nx + j * Nx + k] = interp(newVelPosition, velZ, Nx, Ny, Nz+1);
				}
			}
		}

		velX = newVelX;
		velY = newVelY;
		velZ = newVelZ;
	}

	void pressure_update() {

#pragma omp parallel for schedule(static)
		for (int i = 0; i < Nz; i++) {
			for (int j = 0; j < Ny; j++) {
				velX[i*Ny*(Nx + 1) + j * (Nx + 1) + 0] = 0;
				for (int k = 1; k < Nx; k++) {
					velX[i*Ny*(Nx + 1) + j * (Nx + 1) + k] -= dt / rho * (pressure[i*Ny*Nx + j * Nx + k] - pressure[i*Ny*Nx + j * Nx + k-1]) / dx[0];
				}
				velX[i*Ny*(Nx + 1) + j * (Nx + 1) + Nx] = 0;
			}
		}
#pragma omp parallel for schedule(static)
		for (int i = 0; i < Nz; i++) {
			for (int k = 0; k < Nx; k++) {
				velY[i*(Ny + 1)*Nx + 0 + k] = 0;
			}
			for (int j = 1; j < Ny; j++) {
				for (int k = 0; k < Nx; k++) {
					velY[i*(Ny+1)*Nx + j * Nx + k] -= dt / rho * (pressure[i*Ny*Nx + j * Nx + k] - pressure[i*Ny*Nx + (j-1) * Nx + k]) / dx[1];
				}
			}
			for (int k = 0; k < Nx; k++) {
				velY[i*(Ny + 1)*Nx + Ny*Nx + k] = 0;
			}
		}

		for (int j = 0; j < Ny; j++) {
			for (int k = 0; k < Nx; k++) {
				velZ[0 + j * Nx + k] = 0;
			}
		}
#pragma omp parallel for schedule(static)
		for (int i = 1; i < Nz; i++) {
			for (int j = 0; j < Ny; j++) {
				for (int k = 0; k < Nx; k++) {
					velZ[i*Ny*Nx + j * Nx + k] -= dt / rho * (pressure[i*Ny*Nx + j * Nx + k] - pressure[(i-1)*Ny*Nx + j * Nx + k]) / dx[2];
				}
			}
		}
		for (int j = 0; j < Ny; j++) {
			for (int k = 0; k < Nx; k++) {
				velZ[Nz*Ny*Nx + j * Nx + k] = 0;
			}
		}
	}

	void applyA(const std::vector<double> &b, std::vector<double> &r) {

		memset(&r[0], 0, r.size() * sizeof(double));
		double ncInit = 2*invdx2[0] + 2*invdx2[1] + 2*invdx2[2]; // default diagonal term

#pragma omp parallel for schedule(static)
		for (int id = 0; id < fluidindices.size(); id++) {
			int pixcell = fluidindices[id];
			int k = pixcell % Nx;
			int j = (pixcell / Nx)%Ny;
			int i = (pixcell/(Nx*Ny)) % Nz;

			double  bx, bmx, by, bmy, bz, bmz;
			double nc = ncInit;
			if (k == Nx - 1 || (celltypes[pixcell + 1] == 2)) {
				bx = 0;
				nc -= invdx2[0];
			} else {
				if (celltypes[pixcell + 1] == 0) {
					bx = 0;
				} else {
					bx = b[pixcell + 1];
				}
			}
			if (k == 0 || (celltypes[pixcell - 1] == 2)) {
				bmx = 0;
				nc -= invdx2[0];
			} else {
				if (celltypes[pixcell - 1] == 0) {
					bmx = 0;
				} else {
					bmx = b[pixcell - 1];
				}
			}

			/*if (j == Ny - 1 || (celltypes[pixcell + Nx] == 0)) {  // air interface
				by = 0;
			} else {
				if (celltypes[pixcell + Nx] == 2) {
					by = 0;
					nc -= invdx2[1];
				} else {
					by = b[pixcell + Nx];
				}
			}*/
			if (j == Ny - 1 || (celltypes[pixcell + Nx] == 2)) {  // beware : for ceiling, interface should be air, not solid
				by = 0;
				nc -= invdx2[1];
			} else {
				if (celltypes[pixcell + Nx] == 0) {
					by = 0;
				} else {
					by = b[pixcell + Nx];
				}
			}
			if (j == 0 || (celltypes[pixcell - Nx] == 2)) {
				bmy = 0;
				nc -= invdx2[1];
			} else {
				if (celltypes[pixcell - Nx] == 0) {
					bmy = 0;
				} else {
					bmy = b[pixcell - Nx];
				}
			}
			if (i == Nz - 1 || (celltypes[pixcell + Nx * Ny] == 2)) {
				bz = 0;
				nc -= invdx2[2];
			} else {
				if (celltypes[pixcell + Nx * Ny] == 0) {
					bz = 0;
				} else {
					bz = b[pixcell + Nx * Ny];
				}
			}
			if (i == 0 || (celltypes[pixcell - Nx * Ny] == 2)) {
				bmz = 0;
				nc -= invdx2[2];
			} else {
				if (celltypes[pixcell - Nx * Ny] == 0) {
					bmz = 0;
				} else {
					bmz = b[pixcell - Nx * Ny];
				}
			}
			r[pixcell] = nc * b[pixcell] - (bmx + bx)*invdx2[0] - (bmy + by)*invdx2[1] - (bmz + bz)*invdx2[2];
		}
	}

	void rhs(std::vector<double> &rhs) {
		memset(&rhs[0], 0, rhs.size() * sizeof(double));

#pragma omp parallel for
		for (int i = 0; i < Nz; i++) {
			for (int j = 0; j < Ny; j++) {
				for (int k = 0; k < Nx; k++) {

					int pcell = i * Ny*Nx + j * Nx + k;

					if (celltypes[pcell] != 1) {
						rhs[pcell] = 0;
						continue;
					}

					int pvelX = i * Ny*(Nx + 1) + j * (Nx + 1) + k;
					int pvelY = i * (Ny+1)*Nx + j * Nx + k;
					int pvelZ = pcell;

					double mdivu = (velX[pvelX] - velX[pvelX+1])*invdx[0] + (velY[pvelY] - velY[pvelY + Nx])*invdx[1] + (velZ[pvelZ] - velZ[pvelZ + Ny*Nx])*invdx[2];
					
					double bc = 0;
					if (k == Nx - 1 || (celltypes[pcell + 1] == 2)) {
						double a = rho / (dt*dx[0]);
						bc += a * velX[pvelX + 1];
					}
					if (k == 0 || (celltypes[pcell - 1] == 2)) {
						double a = rho / (dt*dx[0]);
						bc -= a * velX[pvelX];
					}
					//if ( (j != Ny - 1) && (celltypes[pcell+Nx] == 2)) { // ceiling should be air, not solid
					if (j == Ny - 1 || (celltypes[pcell + Nx] == 2)) {
						double a = rho / (dt*dx[1]);
						bc += a * velY[pvelY + Nx];
					}
					if (j == 0 || (celltypes[pcell - Nx] == 2)) {
						double a = rho / (dt*dx[1]);
						bc -= a * velY[pvelY];
					}
					if (i == Nz - 1 || (celltypes[pcell + Ny*Nx] == 2)) {
						double a = rho / (dt*dx[2]);
						bc += a * velZ[pvelZ + Ny*Nx];
					}
					if (i == 0 || (celltypes[pcell - Ny*Nx] == 2)) {
						double a = rho / (dt*dx[2]);
						bc -= a * velZ[pvelZ];
					}

					rhs[pcell] = rho/dt*mdivu + bc;
				}
			}
		}
	}

	void precond(const std::vector<double> &r, std::vector<double> &z) {

		/*z = r;
		return;*/

		// Jacobi to start with
		for (int i = 0; i < Nz; i++) {
			for (int j = 0; j < Ny; j++) {
				for (int k = 0; k < Nx; k++) {
					double nc = 2 / dx2[0] + 2 / dx2[1] + 2 / dx2[2];
					if (celltypes[i*Ny*Nx + j * Nx + k] == 0) {
						nc = 1;
					} else {
						if (k == Nx - 1 || (celltypes[i*Ny*Nx + j * Nx + k + 1] == 2)) {
							nc -= 1/dx2[0];
						}
						if (k == 0 || (celltypes[i*Ny*Nx + j * Nx + k - 1] == 2)) {
							nc -= 1 / dx2[0];
						}
						//if (/*j == Ny - 1 ||*/  (j != Ny - 1) &&  (celltypes[i*Ny*Nx + (j + 1) * Nx + k] == 2)) { // should be air interface at the ceiling
						if (j == Ny - 1 || (celltypes[i*Ny*Nx + (j + 1) * Nx + k] == 2)) { // should be air interface at the ceiling
							nc -= 1 / dx2[1];
						}
						if (j == 0 || (celltypes[i*Ny*Nx + (j - 1) * Nx + k] == 2)) {
							nc -= 1 / dx2[1];
						}
						if (i == Nz - 1 || (celltypes[(i + 1)*Ny*Nx + j * Nx + k] == 2)) {
							nc -= 1 / dx2[2];
						}
						if (i == 0 || (celltypes[(i - 1)*Ny*Nx + j * Nx + k] == 2)) {
							nc -= 1 / dx2[2];
						}
					}
					z[i*Ny*Nx + j * Nx + k] = r[i*Ny*Nx + j * Nx + k] / nc;

				}
			}
		}
	}

	void conjGrad() {
		int N = Nx * Ny*Nz;
		std::vector<double> x(N, 0.);  //x = pressure;
		std::vector<double> Ax(N, 0.);
		std::vector<double> Ap(N, 0.);
		std::vector<double> r(N, 0.);
		std::vector<double> z(N, 0.);
		std::vector<double> p(N, 0.);
		std::vector<double> b(N, 0.);

		fluidindices.clear();
		fluidindices.reserve(N);
		for (int i = 0; i < N; i++) {
			if (celltypes[i] == 1) fluidindices.push_back(i);
		}
		int Nfluid = fluidindices.size();

		rhs(b);

		applyA(x, Ax);
		double rr = 0;
		for (int i=0; i < Nfluid; i++) {
			r[fluidindices[i]] = b[fluidindices[i]] - Ax[fluidindices[i]];
			rr += r[fluidindices[i]] * r[fluidindices[i]];
		}
		if (rr < 1E-14) {
			pressure = x;
			return;
		}
		precond(r, z);
		p = z;

		for (int iter = 0; iter < 10000; iter++) {

			applyA(p, Ap);

			double rz = 0, pAp = 0, rr = 0;
			for (int i = 0; i < Nfluid;  i++) {
				rz += r[fluidindices[i]]*z[fluidindices[i]];
				rr += r[fluidindices[i]] * r[fluidindices[i]];
				pAp += p[fluidindices[i]] * Ap[fluidindices[i]];
			}
			/*FILE* f = fopen("iter.txt", "a+");
			fprintf(f, "%u %f\n", iter, rr);
			fclose(f);*/
			if (std::abs(pAp) < 1E-14) break;

			double ak = rz / pAp;					
			for (int i = 0; i < Nfluid; i++) {
				x[fluidindices[i]] += ak * p[fluidindices[i]];
				r[fluidindices[i]] -= ak * Ap[fluidindices[i]];
			}
			precond(r, z);

			double rz2 = 0;
			for (int i = 0; i < Nfluid; i++) {
				rz2 += r[fluidindices[i]] * z[fluidindices[i]];
			}
			double bk = rz2 / rz;
			if (std::abs(rz) < 1E-14) break;

			for (int i = 0; i < Nfluid; i++) {
				p[fluidindices[i]] = z[fluidindices[i]] + bk * p[fluidindices[i]];
			}
			if (rr < 1E-14) break;
		}

		pressure = x;
	}

	void addforces() {
#pragma omp parallel for schedule(static)
		for (int i = 0; i < Nz; i++) {
			for (int j = 1; j < Ny+1; j++) {
				for (int k = 0; k < Nx; k++) {
					velY[i*(Ny + 1)*Nx + j * Nx + k] -= 9.81*dt;
				}
			}
		}
	}

	void extrapolateVel() {
		// naive NN extrapolation ; should do jump flooding instead

		std::vector<bool> computed(Nz*Ny*Nz, false);
		for (int i = 0; i < Nx*Ny*Nz; i++)
			if (celltypes[i] == 1) computed[i] = true;

		for (int iter = 0; iter < 30; iter++) {
			for (int i = 0; i < Nz; i++) {
				for (int j = 0; j < Ny; j++) {
					for (int k = 1; k < Nx; k++) {
						if (computed[i*Nx*Ny + j * Nx + k - 1] || computed[i*Nx*Ny + j * Nx + k]) continue; //fluid-fluid, fluid-air, fluid-solid: already computed

						//air-air, air-solid, solid-solid boundary : no velocity computed at Vx_ijk => extrapolate

						if (k >= 2 && (computed[i*Nx*Ny + j * Nx + k - 2])) {
							velX[i*(Nx + 1)*Ny + j * (Nx + 1) + k] = velX[i*(Nx + 1)*Ny + j * (Nx + 1) + k - 1];
							computed[i*Nx*Ny + j * Nx + k - 1] = true;
							continue;
						}
						if (k <= Nx - 2 && (computed[i*Nx*Ny + j * Nx + k + 1])) {
							velX[i*(Nx + 1)*Ny + j * (Nx + 1) + k] = velX[i*(Nx + 1)*Ny + j * (Nx + 1) + k + 1];
							computed[i*Nx*Ny + j * Nx + k] = true;
							continue;
						}
					}
				}
			}
			for (int i = 0; i < Nz; i++) {
				for (int j = 1; j < Ny; j++) {
					for (int k = 0; k < Nx; k++) {
						if (computed[i*Nx*Ny + (j - 1) * Nx + k] || computed[i*Nx*Ny + j * Nx + k]) continue; //fluid-fluid, fluid-air, fluid-solid: already computed

						//air-air, air-solid, solid-solid boundary : no velocity computed at Vx_ijk => extrapolate
						if (j >= 2 && (computed[i*Nx*Ny + (j - 2) * Nx + k])) {
							velY[i*Nx*(Ny + 1) + j * Nx + k] = velY[i*Nx*(Ny + 1) + (j - 1) * Nx + k];
							computed[i*Nx*Ny + (j - 1) * Nx + k] = true;
							continue;
						}
						if (j <= Ny - 2 && (computed[i*Nx*Ny + (j + 1) * Nx + k])) {
							velY[i*Nx*(Ny + 1) + (j)* Nx + k] = velY[i*Nx*(Ny + 1) + (j + 1) * Nx + k];
							computed[i*Nx*Ny + j * Nx + k] = true;
							continue;
						}
					}
				}
			}
			for (int i = 1; i < Nz; i++) {
				for (int j = 0; j < Ny; j++) {
					for (int k = 0; k < Nx; k++) {
						if (computed[(i - 1)*Nx*Ny + j * Nx + k] || computed[i*Nx*Ny + j * Nx + k]) continue; //fluid-fluid, fluid-air, fluid-solid: already computed

						//air-air, air-solid, solid-solid boundary : no velocity computed at Vx_ijk => extrapolate

						if (i >= 2 && (computed[(i - 2)*Nx*Ny + j * Nx + k])) {
							velZ[i*Nx*Ny + j * Nx + k] = velZ[(i - 1)*Nx*Ny + j * Nx + k];
							computed[(i - 1)*Nx*Ny + j * Nx + k] = true;
							continue;
						}
						if (i <= Nz - 2 && (computed[(i + 1)*Nx*Ny + j * Nx + k])) {
							velZ[i*Nx*Ny + j * Nx + k] = velZ[(i + 1)*Nx*Ny + j * Nx + k];
							computed[i*Nx*Ny + j * Nx + k] = true;
							continue;
						}
					}
				}
			}
		}


	}

	void moveparticles(std::vector<Vector>& particles_new, const std::vector<Vector>& particles_old) {
#pragma omp parallel for schedule(static)
		for (int i = 0; i < particles_old.size(); i++) {
			double curdt = dt;
			for (int j = 0; j < 30; j++) {
				//particles_new[i] = particles_old[i] + curdt * Vector(interp(particles_old[i], velX, Nx + 1, Ny, Nz), interp(particles_old[i], velY, Nx, Ny + 1, Nz), interp(particles_old[i], velZ, Nx, Ny, Nz + 1));				
				//particles_new[i] = particles[frame-1][i] + 2*curdt * Vector(interp(particles_old[i], velX, Nx + 1, Ny, Nz), interp(particles_old[i], velY, Nx, Ny + 1, Nz), interp(particles_old[i], velZ, Nx, Ny, Nz + 1));
				/*Vector halfwayPos = particles_old[i] + curdt / 2 * Vector(interp(particles_old[i], velX, Nx + 1, Ny, Nz), interp(particles_old[i], velY, Nx, Ny + 1, Nz), interp(particles_old[i], velZ, Nx, Ny, Nz + 1));
				Vector halfwayVel = Vector(interp(halfwayPos, velX, Nx + 1, Ny, Nz), interp(halfwayPos, velY, Nx, Ny + 1, Nz), interp(halfwayPos, velZ, Nx, Ny, Nz + 1));
				particles_new[i] = particles_old[i] + curdt * halfwayVel; // RK2*/

				//particles_new[i] = particles_old[i] + curdt * Vector(interp(particles_old[i] - Vector(0, 0.5*dx[1], 0.5*dx[2]), velX, Nx + 1, Ny, Nz), interp(particles_old[i] - Vector(0.5*dx[0], 0, 0.5*dx[2]), velY, Nx, Ny + 1, Nz), interp(particles_old[i] - Vector(0.5*dx[0], 0.5*dx[1], 0), velZ, Nx, Ny, Nz + 1));

			    Vector k1 = Vector(interp(particles_old[i] - Vector(0, 0.5*dx[1], 0.5*dx[2]), velX, Nx + 1, Ny, Nz), interp(particles_old[i] - Vector(0.5*dx[0], 0, 0.5*dx[2]), velY, Nx, Ny + 1, Nz), interp(particles_old[i] - Vector(0.5*dx[0], 0.5*dx[1], 0), velZ, Nx, Ny, Nz + 1));
				Vector k2 = Vector(interp(particles_old[i] - Vector(0, 0.5*dx[1], 0.5*dx[2]) + curdt * 0.5 * k1, velX, Nx + 1, Ny, Nz), interp(particles_old[i] - Vector(0.5*dx[0], 0, 0.5*dx[2]) + curdt * 0.5*k1, velY, Nx, Ny + 1, Nz), interp(particles_old[i] - Vector(0.5*dx[0], 0.5*dx[1], 0) + curdt * 0.5*k1, velZ, Nx, Ny, Nz + 1));
				Vector k3 = Vector(interp(particles_old[i] - Vector(0, 0.5*dx[1], 0.5*dx[2]) + curdt * 0.5 * k2, velX, Nx + 1, Ny, Nz), interp(particles_old[i] - Vector(0.5*dx[0], 0, 0.5*dx[2]) + curdt * 0.5*k2, velY, Nx, Ny + 1, Nz), interp(particles_old[i] - Vector(0.5*dx[0], 0.5*dx[1], 0) + curdt * 0.5*k2, velZ, Nx, Ny, Nz + 1));
				Vector k4 = Vector(interp(particles_old[i] - Vector(0, 0.5*dx[1], 0.5*dx[2]) + curdt * k3, velX, Nx + 1, Ny, Nz), interp(particles_old[i] - Vector(0.5*dx[0], 0, 0.5*dx[2]) + curdt *k3, velY, Nx, Ny + 1, Nz), interp(particles_old[i] - Vector(0.5*dx[0], 0.5*dx[1], 0) + curdt *k3, velZ, Nx, Ny, Nz + 1));
				particles_new[i] = particles_old[i] + curdt/6.*(k1+2*k2+2*k3+k4);

				Vector pgrid = (particles_new[i] - extent.bounds[0]) / dx;
				int pi[3] = { pgrid[0], pgrid[1] , pgrid[2] };
				if ((pi[0] < 0 || pi[1] < 0 || pi[2] < 0) || (pi[0] >= Nx || pi[1] >= Ny || pi[2] >= Nz) || (celltypes[pi[2] * Ny*Nx + pi[1] * Nx + pi[0]] == 2)) {
					curdt *= 0.75;
				} else break;
			}
		}
	}

	void timestep() {
		advect();
		extrapolateVel2();
		addforces();
		conjGrad();
		pressure_update();
		//extrapolateVel2();
		//extrapolateVel();
		
		moveparticles(visualparticlesnew, visualparticles);
		visualparticles = visualparticlesnew;

		moveparticles(ghostparticlesnew, ghostparticles);
		ghostparticles = ghostparticlesnew;

#pragma omp parallel for schedule(static)
		for (int i = 0; i < Nx*Ny*Nz; i++)
			if (celltypes[i] == 1) celltypes[i] = 0;

		int fluidcount = 0;
		for (int i = 0; i < ghostparticles.size(); i++) {
			Vector pgrid = (ghostparticles[i] - extent.bounds[0]) / dx;
			int pi[3] = {pgrid[0], pgrid[1] , pgrid[2] };
			if (pi[0] < 0 || pi[1] < 0 || pi[2] < 0) continue;
			if (pi[0] >= Nx || pi[1] >= Ny || pi[2] >= Nz) continue;
			if (celltypes[pi[2] * Ny*Nx + pi[1] * Nx + pi[0]] != 0) continue;
			celltypes[pi[2]*Ny*Nx+pi[1]*Nx+pi[0]] = 1;
			fluidcount++;
		}
		FILE* f = fopen("fluidcount.txt", "a+");
		fprintf(f, "%u\n", fluidcount);
		fclose(f);
		for (int i = 0; i < visualparticles.size(); i++) {
			Vector pgrid = (visualparticles[i] - extent.bounds[0]) / dx;
			int pi[3] = { pgrid[0], pgrid[1] , pgrid[2] };
			if (pi[0] < 0 || pi[1] < 0 || pi[2] < 0) continue;
			if (pi[0] >= Nx || pi[1] >= Ny || pi[2] >= Nz) continue;
			if (celltypes[pi[2] * Ny*Nx + pi[1] * Nx + pi[0]] != 0) continue;
			celltypes[pi[2] * Ny*Nx + pi[1] * Nx + pi[0]] = 1;
		}

		/*std::vector<char> newcelltypes = celltypes;
#pragma omp parallel for schedule(static)
		for (int i = 0; i < Nz; i++) {
			for (int j = 0; j < Ny; j++) {
				for (int k = 0; k < Nx; k++) {

					// get vector coordinates at staggered grid location by bilinear interpolation
					Vector curVel;
					curVel[0] = (velX[i*Ny*(Nx + 1) + j * (Nx + 1) + k] + velX[i*Ny*(Nx + 1) + j * (Nx + 1) + k+1] )*0.5;
					curVel[1] = (velY[i*(Ny + 1)*Nx + j * Nx + k] + velY[i*(Ny + 1)*Nx + (j + 1) * Nx + k])*0.5;
					curVel[2] = (velZ[i*Ny*Nx + j * Nx + k] + velZ[(i+1)*Ny*Nx + j * Nx + k]) * 0.5;

					int id[3] = { (int)(k + 0.5 - curVel[0] * dt / dx[0]), (int)(j + 0.5 - curVel[1] * dt / dx[1]), (int)(i + 0.5 - curVel[2] * dt / dx[2]) };
					id[0] = std::max(0, std::min(Nx-1, id[0]));
					id[1] = std::max(0, std::min(Ny-1, id[1]));
					id[2] = std::max(0, std::min(Nz-1, id[2]));
					if (celltypes[i*Ny*Nx + j * Nx + k] != 2)
						newcelltypes[i*Ny*Nx + j * Nx + k] = celltypes[id[2]*Ny*Nx + id[1] * Nx + id[0]];
				}
			}
		}
		celltypes = newcelltypes;*/

	}

	void run() {
		for (int i = 0; i < nbframes; i++) {
			visualparticles = particles[i];
			visualparticlesnew.resize(particles[i].size());
			for (int j = 0; j < nsubsteps; j++) {				
				timestep();
			}
			particles[i + 1] = visualparticlesnew;
			FILE* f = fopen("state.txt", "a+");
			fprintf(f, "step %u done\n", i);
			fclose(f);
		}
		build_bvh(0, 0, Nparticles);
		build_grid(0);
	}



	bool intersection(const Ray& d, Vector& P, double &t, MaterialValues &mat, double cur_best_t, int &triangle_id) const;
	bool intersection_shadow(const Ray& d, double &t, double cur_best_t, double dist_light) const;
	bool intersection_transparent(const Ray& d, Vector& P, double &t, MaterialValues &mat, double cur_best_t, int &triangle_id) const;
	bool intersection_opaque(const Ray& d, Vector& P, double &t, MaterialValues &mat, double cur_best_t, int &triangle_id) const;

	BBox build_bbox(int frame, int i0, int i1);
	BBox build_centers_bbox(int frame, int i0, int i1);
	void build_bvh(int frame, int i0, int i1);
	void build_bvh_recur(int frame, BVH* b, int node, int i0, int i1);

	void build_grid(int frame);
	bool intersection_transparent2(const Ray& d, Vector& P, double &t, MaterialValues &mat, double cur_best_t, int &triangle_id) const;
	bool intersection_opaque2(const Ray& d, Vector& P, double &t, MaterialValues &mat, double cur_best_t, int &triangle_id) const;
	bool intersection_shadow2(const Ray& d, double &t, double cur_best_t, double dist_light) const;
		
	BBox extent;
	int Nx, Ny, Nz, Nparticles;
	Vector Nv;
	std::vector<double> velX, velY, velZ, pressure;
	std::vector<std::vector<Vector> > particles;  // decouples particles used for simulation and for rendering
	std::vector<Vector> ghostparticles, ghostparticlesnew, visualparticles, visualparticlesnew, visualparticlescolor;
	std::vector<int> fluidindices;
	std::vector<char> celltypes;
	std::vector<Vector> cellcol;
	Vector dx, dx2, invdx, invdx2;
	double rho, radius, dt;
	BVH bvh;
	int nbframes, nsubsteps;
	std::vector<double> newVelX, newVelY, newVelZ;
	bool opaque;
	Scene& scene;

	mutable std::pair<double, std::pair<double, int> >  allts[64][500000];
	mutable int nodelist[64][500000];
	mutable double tnearlist[64][500000];

	std::vector<char> accelgrid;
	std::vector<std::vector<int> > accelgridindices;
};