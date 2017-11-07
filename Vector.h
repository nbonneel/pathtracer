#pragma once
#include <math.h>
#include <vector>
#include <sstream>
#include <random>

#ifndef M_PI
#define M_PI 3.1415926535897932
#endif
#ifdef min
#undef min
#undef max
#endif

/*static __declspec(thread) std::default_random_engine engine;
static __declspec(thread) std::uniform_real_distribution<double> uniform(0, 1);*/

static thread_local std::default_random_engine engine;
static thread_local std::uniform_real_distribution<double> uniform(0, 1);

class Vector;

Vector operator+(const Vector& a, const Vector &b);
Vector operator-(const Vector& a, const Vector &b);
Vector operator*(double a, const Vector &b);
Vector operator*(const Vector &b, double a);
Vector operator*(const Vector &a, const Vector &b);
Vector operator/(const Vector& a, double b);
Vector operator-(const Vector& a);
double dot(const Vector&a, const Vector& b);
Vector cross(const Vector&a, const Vector& b);

template<int M, int N>
class Matrix {
public:
	Matrix() {
		memset(values, 0, M*N * sizeof(double));
		for (int i = 0; i < N; i++) {
			values[i*N + i] = 1;
		}
	}
	Matrix<N, M> getTransposed() {
		Matrix<N, M> res;
		for (int i=0; i<M; i++)
			for (int j = 0; j < N; j++) {
				res[j*M + i] = values[i*N+j];
			}
		return res;
	}
	double operator[](int i) const { return values[i]; }
	double& operator[](int i)  { return values[i]; }
	double values[M*N];
};

static inline Matrix<3, 3> createRotationMatrixX(double angleX) {
	Matrix<3, 3> R;
	R[0] = cos(angleX);
	R[1] = -sin(angleX);
	R[3] = sin(angleX);
	R[4] = cos(angleX);
	return R;
}
static inline Matrix<3, 3> createRotationMatrixY(double angleY) {
	Matrix<3, 3> R;
	R[0] = cos(angleY);
	R[2] = sin(angleY);
	R[6] = -sin(angleY);
	R[8] = cos(angleY);
	return R;
}
static inline Matrix<3, 3> createRotationMatrixZ(double angleZ) {
	Matrix<3, 3> R;
	R[4] = cos(angleZ);
	R[5] = -sin(angleZ);
	R[7] = sin(angleZ);
	R[8] = cos(angleZ);
	return R;
}

class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {
		coord[0] = x;
		coord[1] = y;
		coord[2] = z;
	}
	const double& operator[](int i) const { return coord[i]; }
	double& operator[](int i) { return coord[i]; }

	double getNorm2() const {
		return coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2];
	}
	void normalize() {
		double norm = sqrt(getNorm2());
		coord[0] /= norm;
		coord[1] /= norm;
		coord[2] /= norm;
	}
	Vector getNormalized() const {
		Vector result(*this);
		result.normalize();
		return result;
	}
	Vector reflect(const Vector &N) const {
		Vector result = *this - 2.*dot(*this, N)*N;
		return result;
	}
	std::string toColorStr() {
		std::ostringstream os;
		os << "Color(" << coord[0] << ", " << coord[1] << ", " << coord[2] << ")";
		return os.str();
	}

	Vector& operator+=(const Vector& b) {
		coord[0] += b[0];
		coord[1] += b[1];
		coord[2] += b[2];
		return *this;
	}
	Vector& operator*=(double b) {
		coord[0] *= b;
		coord[1] *= b;
		coord[2] *= b;
		return *this;
	}
	Vector& operator/=(double b) {
		coord[0] /= b;
		coord[1] /= b;
		coord[2] /= b;
		return *this;
	}

private:
	double coord[3];
};

static inline  Vector operator*(const Matrix<3, 3> &mat, const Vector &b) {

	Vector res(0, 0, 0);
	for (int i = 0; i < 3; i++) {
		double v = 0;
		for (int j = 0; j < 3; j++) {
			v += mat[i*3 + j] * b[j];
		}
		res[i] = v;
	}
	return res;
}

static inline Vector operator*(const Matrix<3, 4> &mat, const Vector &b) {
	Vector res(0, 0, 0);
	for (int i = 0; i < 3; i++) {
		double v = 0;
		for (int j = 0; j < 3; j++) {
			v += mat[i*4 + j] * b[j];
		}
		res[i] = v;
	}
	for (int i = 0; i < 3; i++) {
		res[i] += mat[i*4 + 3];
	}
	return res;
}

template<int M, int N, int O>
Matrix<M, O> operator*(const Matrix<M, N> &matA, const Matrix<N, O> &matB) {

	Matrix<M, O> result;
	for (int i = 0; i < M; i++)
		for (int j = 0; j < O; j++) {
			double v = 0;
			for (int k = 0; k < N; k++) {
				v += matA[i*N + k] * matB[k*O + j];
			}
			result[i*O + j] = v;
		}
	return result;
}


Vector min(const Vector& a, const Vector& b);
Vector max(const Vector& a, const Vector& b);
Vector pow(const Vector& a, const Vector& b);
Vector random_cos(const Vector &N);
Vector random_Phong(const Vector &R, double phong_exponent);
Vector random_uniform();

Vector rotate_dir(const Vector&v, const Vector &angles);
Vector inverse_rotate_dir(const Vector&v, const Vector &angles); // performs the rotation by -angle, but in the reversed order 

class Ray {
public:
	Ray() {};
	Ray(const Vector& o, const Vector& d, double time) : origin(o), direction(d), time(time) {};
	Vector origin, direction;
	double time;
};


class Camera {
public:
	Camera() {};
	Camera(const Vector& position, const Vector &direction, const Vector& up): position(position), direction(direction), up(up), initial_position(position), initial_direction(direction), initial_up(up) {};

	void translate(const Vector& translation, double time) {
		position = initial_position + time*translation;
	}

	void rotate(double angle_x, double angle_y, double time) {
		double cur_angle_x = time*angle_x;
		double cur_angle_y = time*angle_y;

		Vector dirtmp;

		dirtmp[0] = direction[0];
		dirtmp[1] = cos(cur_angle_y)*direction[1] - sin(cur_angle_y)*direction[2];
		dirtmp[2] = sin(cur_angle_y)*direction[1] + cos(cur_angle_y)*direction[2];

		direction[0] = cos(cur_angle_x)*dirtmp[0] - sin(cur_angle_x)*dirtmp[2];
		direction[1] = dirtmp[1];
		direction[2] = sin(cur_angle_x)*dirtmp[0] + cos(cur_angle_x)*dirtmp[2];


		Vector uptmp;
		uptmp[0] = up[0];
		uptmp[1] = cos(cur_angle_y)*up[1] - sin(cur_angle_y)*up[2];
		uptmp[2] = sin(cur_angle_y)*up[1] + cos(cur_angle_y)*up[2];
		
		up[0] = cos(cur_angle_x)*uptmp[0] - sin(cur_angle_x)*uptmp[2];
		up[1] = uptmp[1];
		up[2] = sin(cur_angle_x)*uptmp[0] + cos(cur_angle_x)*uptmp[2];


	}

	void rotateAroundRight(double angle) {

		Vector dirtmp, uptmp;
		dirtmp =  sin(angle)*up  + cos(angle)*direction;		
		uptmp = cos(angle)*up - sin(angle)*direction;

		direction = dirtmp;
		up = uptmp;
	}
	void rotateAroundUp(double angle) {

		Vector right = cross(up, direction);
		direction = -sin(angle)*right + cos(angle)*direction;

	}

	void rotateAroundAxes(double angle_x, double angle_y, double time) {
		double cur_angle_x = time*angle_x;
		double cur_angle_y = time*angle_y;

		rotateAroundUp(cur_angle_x);
		rotateAroundRight(cur_angle_y);
	}
	

	Ray generateDirection(int i, int j, double time, double dx_sensor, double dy_sensor, double dx_aperture, double dy_aperture, int W, int H) {
		Vector direction_vec(j - W / 2 + 0.5 + dx_sensor, i - H / 2 + 0.5 + dy_sensor, W / (2 * tan(fov / 2)));
		direction_vec.normalize();
		Vector camera_right = cross(direction, up);
		direction_vec = camera_right * direction_vec[0] + up*direction_vec[1] + direction*direction_vec[2];

		Vector destination = position + focus_distance * direction_vec;
		Vector new_origin = position + Vector(dx_aperture, dy_aperture, 0);
		return Ray(new_origin, (destination - new_origin).getNormalized(), time);

	}

	Vector position, direction, up;
	double fov;             // field of view 
	double focus_distance; // distance mise au point
	double aperture;       // ouverture, pour depth of field (faible valeur = net partout)

private:
	Vector initial_position, initial_direction, initial_up;
};


