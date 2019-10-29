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
		R = rayon;
	}
	bool intersection(const Ray& d, Vector& P, Vector &N, double &t) const {
		// resout a*t^2 + b*t + c = 0
		double a = d.direction.getNorm2();
		double b = 2 * dot(d.direction, d.origin - O);
		double c = (d.origin - O).getNorm2() - R * R;

		double delta = b * b - 4 * a*c;
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

		P = d.origin + t * d.direction;
		N = (P - O).getNormalized();
		return true;
	}

	bool both_intersections(const Ray& d, double &t1, double&t2) const {
		// resout a*t^2 + b*t + c = 0
		double a = d.direction.getNorm2();
		double b = 2 * dot(d.direction, d.origin - O);
		double c = (d.origin - O).getNorm2() - R * R;

		double delta = b * b - 4 * a*c;
		if (delta < 0) return false;

		double sqDelta = sqrt(delta);
		double inv2a = 1. / (2. * a);
		t2 = (-b + sqDelta) * inv2a;
		t1 = (-b - sqDelta) * inv2a;
		return true;
	}

	bool intersection_shadow(const Ray& d, double &t) const {
		double a = d.direction.getNorm2();
		double b = 2 * dot(d.direction, d.origin - O);
		double c = (d.origin - O).getNorm2() - R * R;

		double delta = b * b - 4 * a*c;
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
};

class Fluid : public Object {	
public:
	Fluid(BBox extent, int Nx, int Ny, int Nz, int Nparticles, double rho, double sphere_radius, int nbframes, double dt) :extent(extent), Nx(Nx), Ny(Ny), Nz(Nz), Nparticles(Nparticles), rho(rho), radius(sphere_radius), nbframes(nbframes), dt(dt){
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
		particles.resize(nbframes + 1, std::vector<Vector>(Nparticles));

		std::random_device dev;
		std::mt19937 rng(dev());
		std::uniform_real_distribution<double> uniform(0, 1);
		for (int i = 0; i < Nparticles; i++) {
			particles[0][i] = Vector(uniform(rng)*1./3.+1./3., uniform(rng)/2+0.4, uniform(rng)*1./3.+1./3.)*(extent.bounds[1] - extent.bounds[0]) + extent.bounds[0];
		}
#pragma omp parallel for schedule(static)
		for (int i = Nz/3; i < (2*Nz)/3; i++)
			for (int j =(Ny*0.4); j < (Ny*0.9); j++)
				for (int k = Nx/3; k < (2*Nx)/3; k++) {
					celltypes[i*Ny*Nx + j * Nx + k] = 1;
				}
	
	}

	template<typename T>
	T interp(const Vector& pos, const std::vector<T> &field, int Nx, int Ny, int Nz) {
		Vector posNormalized = (pos - extent.bounds[0])/dx;
		int icoords[3] = {posNormalized[0], posNormalized[1] , posNormalized[2]};
		if (icoords[0] >= Nx - 1 || icoords[1] >= Ny - 1 || icoords[2] >= Nz - 1) return T(0);
		if (icoords[0] < 0 || icoords[1] <0 || icoords[2] <0) return T(0);

		Vector fcoords(posNormalized[0] - icoords[0], posNormalized[1] - icoords[1], posNormalized[2] - icoords[2]);

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
					newVelPosition[1] = std::max(std::min(newVelPosition[1] - 0.5*dx[1], extent.bounds[1][1]), extent.bounds[0][1]);
					newVelPosition[2] = std::max(std::min(newVelPosition[2] - 0.5*dx[2], extent.bounds[1][2]), extent.bounds[0][2]);

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
					newVelPosition[0] = std::max(std::min(newVelPosition[0] - 0.5*dx[0], extent.bounds[1][0]), extent.bounds[0][0]);  // stick to walls ; should change
					newVelPosition[1] = std::max(std::min(newVelPosition[1], extent.bounds[1][1]), extent.bounds[0][1]);
					newVelPosition[2] = std::max(std::min(newVelPosition[2] - 0.5*dx[2], extent.bounds[1][2]), extent.bounds[0][2]);

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
					newVelPosition[0] = std::max(std::min(newVelPosition[0] - 0.5*dx[0], extent.bounds[1][0]), extent.bounds[0][0]);  // stick to walls ; should change
					newVelPosition[1] = std::max(std::min(newVelPosition[1] - 0.5*dx[1], extent.bounds[1][1]), extent.bounds[0][1]);
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
				for (int k = 1; k < Nx; k++) {
					velX[i*Ny*(Nx + 1) + j * (Nx + 1) + k] -= dt / rho * (pressure[i*Ny*Nx + j * Nx + k] - pressure[i*Ny*Nx + j * Nx + k-1]) / dx[0];
				}
			}
		}
#pragma omp parallel for schedule(static)
		for (int i = 0; i < Nz; i++) {
			for (int j = 1; j < Ny; j++) {
				for (int k = 0; k < Nx; k++) {
					velY[i*(Ny+1)*Nx + j * Nx + k] -= dt / rho * (pressure[i*Ny*Nx + j * Nx + k] - pressure[i*Ny*Nx + (j-1) * Nx + k]) / dx[1];
				}
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
	}

	void applyA(const std::vector<double> &b, std::vector<double> &r) {

		memset(&r[0], 0, r.size() * sizeof(double));
		double ncInit = 2*invdx2[0] + 2*invdx2[1] + 2*invdx2[2]; // default diagonal term

#pragma omp parallel for
		for (int i = 0; i < Nz; i++) {
			for (int j = 0; j < Ny; j++) {
				for (int k = 0; k < Nx; k++) {
					
					int pixcell = i * Ny*Nx + j * Nx + k;					
					if (celltypes[pixcell] == 0) {
						r[pixcell] = b[pixcell];
						continue;
					}

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

					if (j == Ny-1 || (celltypes[pixcell + Nx] == 0)) {  // beware : for ceiling, interface should be air, not solid
						by = 0;						
					} else {
						if (celltypes[pixcell + Nx] == 2) {
							by = 0;		
							nc -= invdx2[1];
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
					if (i == Nz - 1 || (celltypes[pixcell + Nx*Ny] == 2)) {
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
					r[pixcell] = (nc * b[pixcell] - (bmx + bx)*invdx2[0] - (bmy + by)*invdx2[1]- (bmz+bz)*invdx2[2]);
				}
			}
		}
	}

	void rhs(std::vector<double> &rhs) {
		memset(&rhs[0], 0, rhs.size() * sizeof(double));

#pragma omp parallel for
		for (int i = 0; i < Nz; i++) {
			for (int j = 0; j < Ny; j++) {
				for (int k = 0; k < Nx; k++) {

					int pcell = i * Ny*Nx + j * Nx + k;

					if (celltypes[pcell] == 0) {
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
					if (/*j == Ny - 1 ||*/ (celltypes[pcell+Nx] == 2)) { // ceiling should be air, not solid
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

		z = r;
		return;

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
						if (/*j == Ny - 1 ||*/ (celltypes[i*Ny*Nx + (j + 1) * Nx + k] == 2)) { // should be air interface at the ceiling
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
		std::vector<double> x(N, 0.);
		std::vector<double> Ax(N, 0.);
		std::vector<double> Ap(N, 0.);
		std::vector<double> r(N, 0.);
		std::vector<double> z(N, 0.);
		std::vector<double> p(N, 0.);
		std::vector<double> b(N, 0.);
		rhs(b);

		applyA(x, Ax);
		double rr = 0;
		for (int i=0; i < N; i++) {
			r[i] = b[i] - Ax[i];
			rr += r[i] * r[i];
		}
		if (rr < 1E-10) {
			pressure = x;
			return;
		}
		precond(r, z);
		p = z;

		for (int iter = 0; iter < 1000; iter++) {

			applyA(p, Ap);

			double rz = 0, pAp = 0, rr = 0;
			for (int i = 0; i < N;  i++) {
				rz += r[i]*z[i];
				rr += r[i] * r[i];
				pAp += p[i] * Ap[i];
			}
			/*FILE* f = fopen("iter.txt", "a+");
			fprintf(f, "%u %f\n", iter, rr);
			fclose(f);*/
			double ak = rz / pAp;			
			for (int i = 0; i < N; i++) {
				x[i] += ak * p[i];
				r[i] -= ak * Ap[i];				
			}
			precond(r, z);

			double rz2 = 0;
			for (int i = 0; i < N; i++) {
				rz2 += r[i] * z[i];
			}
			double bk = rz2 / rz;
			for (int i = 0; i < N; i++) {
				p[i] = z[i] + bk * p[i];
			}
			if (rr < 1E-10) break;
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

	void timestep(int frame) {
		advect();
		addforces();
		conjGrad();
		pressure_update();

#pragma omp parallel for schedule(static)
		for (int i = 0; i < Nparticles; i++) {
			particles[frame+1][i] = particles[frame][i] + dt * Vector(interp(particles[frame][i], velX, Nx + 1, Ny, Nz), interp(particles[frame][i], velY, Nx, Ny + 1, Nz), interp(particles[frame][i], velZ, Nx, Ny, Nz + 1));
		}

#pragma omp parallel for schedule(static)
		for (int i = 0; i < Nx*Ny*Nz; i++)
			if (celltypes[i] == 1) celltypes[i] = 0;

		for (int i = 0; i < Nparticles; i++) {
			Vector pgrid = (particles[frame+1][i] - extent.bounds[0]) / (extent.bounds[1] - extent.bounds[0]) * Nv;
			int pi[3] = {pgrid[0], pgrid[1] , pgrid[2] };
			if (pi[0] < 0 || pi[1] < 0 || pi[2] < 0) continue;
			if (pi[0] >= Nx || pi[1] >= Ny || pi[2] >= Nz) continue;
			if (celltypes[pi[2] * Ny*Nx + pi[1] * Nx + pi[0]] != 0) continue;
			celltypes[pi[2]*Ny*Nx+pi[1]*Nx+pi[0]] = 1;
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
			timestep(i);
		}
		build_bvh(0, 0, Nparticles);
	}



	bool intersection(const Ray& d, Vector& P, double &t, MaterialValues &mat, double cur_best_t, int &triangle_id) const;
	bool intersection_shadow(const Ray& d, double &t, double cur_best_t, double dist_light) const;
	bool intersection_transparent(const Ray& d, Vector& P, double &t, MaterialValues &mat, double cur_best_t, int &triangle_id) const;
	bool intersection_opaque(const Ray& d, Vector& P, double &t, MaterialValues &mat, double cur_best_t, int &triangle_id) const;

	BBox build_bbox(int frame, int i0, int i1);
	BBox build_centers_bbox(int frame, int i0, int i1);
	void build_bvh(int frame, int i0, int i1);
	void build_bvh_recur(int frame, BVH* b, int node, int i0, int i1);

	BBox extent;
	int Nx, Ny, Nz, Nparticles;
	Vector Nv;
	std::vector<double> velX, velY, velZ, pressure;
	std::vector<std::vector<Vector> > particles;
	std::vector<char> celltypes;
	Vector dx, dx2, invdx, invdx2;
	double rho, radius, dt;
	BVH bvh;
	int nbframes;
	std::vector<double> newVelX, newVelY, newVelZ;
	bool opaque;
};