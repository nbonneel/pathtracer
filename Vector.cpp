#include "Vector.h"
#include <cmath>
#include <algorithm>
#include <omp.h>

Vector operator+(const Vector& a, const Vector &b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector &b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(double a, const Vector &b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector &b, double a) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator/(const Vector& a, double b) {
	return Vector(a[0]/b, a[1]/b, a[2]/b);
}
Vector operator/(const Vector& a, const Vector &b) {
	return Vector(a[0] / b[0], a[1] / b[1], a[2] / b[2]);
}
Vector operator-(const Vector& a) {
	return Vector(-a[0], -a[1], -a[2]);
}
Vector operator*(const Vector &a, const Vector &b) {
	return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}
Vector min(const Vector& a, const Vector& b) {
	return Vector(std::min(a[0], b[0]), std::min(a[1], b[1]), std::min(a[2], b[2]));
}
Vector max(const Vector& a, const Vector& b) {
	return Vector(std::max(a[0], b[0]), std::max(a[1], b[1]), std::max(a[2], b[2]));
}
Vector pow(const Vector& a, const Vector& b) {
	return Vector(std::pow(a[0], b[0]), std::pow(a[1], b[1]), std::pow(a[2], b[2]));
}
double dot(const Vector&a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vector cross(const Vector&a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

Vector getTangent(const Vector& N) {
	Vector tangent1;
	Vector absN(abs(N[0]), abs(N[1]), abs(N[2]));
	if (absN[0] <= absN[1] && absN[0] <= absN[2]) {
		tangent1 = Vector(0, -N[2], N[1]);
	} else
		if (absN[1] <= absN[0] && absN[1] <= absN[2]) {
			tangent1 = Vector(-N[2], 0, N[0]);
		} else
			tangent1 = Vector(-N[1], N[0], 0);
	tangent1.normalize();
	return tangent1;
}

Vector random_cos(const Vector &N, double r1, double r2) {
	double sr2 = sqrt(1. - r2);

	Vector direction_aleatoire_repere_local(cos(2. * M_PI*r1)*sr2, sin(2. * M_PI*r1)*sr2, sqrt(r2));
	Vector tangent1 = getTangent(N);
	Vector tangent2 = cross(tangent1, N);
	return direction_aleatoire_repere_local[2] * N + direction_aleatoire_repere_local[0] * tangent1 + direction_aleatoire_repere_local[1] * tangent2;
}

Vector random_cos(const Vector &N) {
	int threadid = omp_get_thread_num();
	/*float r1 = uniformf(engine[threadid]);
	float r2 = uniformf(engine[threadid]);*/
	const float invmax = 1.f / engine[threadid].max();
	float r1 = engine[threadid]()*invmax;
	float r2 = engine[threadid]()*invmax;
	return random_cos(N, r1, r2);
}



Vector random_uniform_sphere() {
	int threadid = omp_get_thread_num();
	float invmax = 1.f/engine[threadid].max();
	float r1 = engine[threadid]()*invmax; // uniform(engine[threadid]);
	float r2 = engine[threadid]()*invmax; // uniform(engine[threadid]);
	Vector result;
	result[0] = 2.*cos(2.*M_PI*r1)*sqrt(r2*(1 - r2));
	result[1] = 2.*sin(2.*M_PI*r1)*sqrt(r2*(1 - r2));
	result[2] = 1 - 2 * r2;
	return result;
}
Vector random_uniform_hemisphere(const Vector& N) {
	int threadid = omp_get_thread_num();
	float invmax = 1.f / engine[threadid].max();
	float r1 = engine[threadid]()*invmax; // uniform(engine[threadid]);
	float r2 = engine[threadid]()*invmax; // uniform(engine[threadid]);
	Vector resultlocal;
	resultlocal[0] = cos(2.*M_PI*r1)*sqrt(1 - r2*r2);
	resultlocal[1] = sin(2.*M_PI*r1)*sqrt(1 - r2 * r2);
	resultlocal[2] = r2;
	Vector tangent1 = getTangent(N);
	Vector tangent2 = cross(tangent1, N);
	return resultlocal[2] * N + resultlocal[0] * tangent1 + resultlocal[1] * tangent2;
}
Vector random_uniform_ball() {
	int threadid = omp_get_thread_num();
	float invmax = 1.f / engine[threadid].max();
	Vector result = std::pow(engine[threadid]()*invmax, 1./3.)*random_uniform_sphere();
	return result;
}
Vector random_uniform_hemiball(const Vector &N) {
	int threadid = omp_get_thread_num();
	float invmax = 1.f / engine[threadid].max();
	Vector result = std::pow(engine[threadid]()*invmax, 1. / 3.)*random_uniform_hemisphere(N);
	return result;
}
Vector boxMuller() { // returns radius in third component
	int threadid = omp_get_thread_num();
	float invmax = 1.f / engine[threadid].max();
	float r1 = engine[threadid]()*invmax;
	float r2 = engine[threadid]()*invmax;
	float s1 = sqrt(-2 * log(r1));
	float s2 = 2 * M_PI*r2;
	return Vector(s1*cos(s2), s1*sin(s2), s1);
}



Vector rotate_dir(const Vector&v, const Vector &angles) {

	Vector sinangles(sin(angles[0]), sin(angles[1]), sin(angles[2]));
	Vector cosangles(cos(angles[0]), cos(angles[1]), cos(angles[2]));

	Vector v1;
	v1[0] = v[0];
	v1[1] = cosangles[0] * v[1] - sinangles[0] * v[2];
	v1[2] = sinangles[0] * v[1] + cosangles[0] * v[2];

	Vector v2;	
	v2[0] = sinangles[1] * v1[2] + cosangles[1] * v1[0];
	v2[1] = v1[1];	
	v2[2] = cosangles[1] * v1[2] - sinangles[1] * v1[0];

	Vector v3;
	v3[0] = cosangles[2] * v2[0] - sinangles[2] * v2[1];
	v3[1] = sinangles[2] * v2[0] + cosangles[2] * v2[1];
	v3[2] = v2[2];

	return v3;
}

Vector inverse_rotate_dir(const Vector&v, const Vector &angles) {

	Vector sinangles(sin(-angles[0]), sin(-angles[1]), sin(-angles[2]));
	Vector cosangles(cos(-angles[0]), cos(-angles[1]), cos(-angles[2]));

	Vector v1;
	v1[0] = cosangles[2] * v[0] - sinangles[2] * v[1];
	v1[1] = sinangles[2] * v[0] + cosangles[2] * v[1];
	v1[2] = v[2];

	Vector v2;
	v2[0] = sinangles[1] * v1[2] + cosangles[1] * v1[0];
	v2[1] = v1[1];
	v2[2] = cosangles[1] * v1[2] - sinangles[1] * v1[0];

	Vector v3;
	v3[0] = v2[0];
	v3[1] = cosangles[0] * v2[1] - sinangles[0] * v2[2];
	v3[2] = sinangles[0] * v2[1] + cosangles[0] * v2[2];


	return v3;

}