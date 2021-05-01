#pragma once
#include <math.h>
#include <vector>
#include <sstream>
#include <cstring>
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

static thread_local std::default_random_engine engine[64];
static thread_local std::uniform_real_distribution<double> uniform(0, 1);

class Vector;

Vector operator+(const Vector& a, const Vector &b);
Vector operator-(const Vector& a, const Vector &b);
Vector operator*(double a, const Vector &b);
Vector operator*(const Vector &b, double a);
Vector operator*(const Vector &a, const Vector &b);
Vector operator/(const Vector& a, double b);
Vector operator/(const Vector& a, const Vector &b);
Vector operator-(const Vector& a);
double dot(const Vector&a, const Vector& b);
Vector cross(const Vector&a, const Vector& b);

template<typename T>
class Quaternion
{
public:
	Quaternion() { }
	Quaternion(double wVal, double xVal, double yVal, double zVal)
	{
		w = wVal; x = xVal; y = yVal; z = zVal;
	}

	T w, x, y, z;
};
template<typename T>
Quaternion<T> operator*(const Quaternion<T>& q1, const Quaternion<T>& q2) {
	T w1, x1, y1, z1, w2, x2, y2, z2, w3, x3, y3, z3;

	w1 = q1.w; x1 = q1.x; y1 = q1.y; z1 = q1.z;
	w2 = q2.w; x2 = q2.x; y2 = q2.y; z2 = q2.z;

	w3 = w1*w2 - x1*x2 - y1*y2 - z1*z2;
	x3 = w1*x2 + x1*w2 + y1*z2 - z1*y2;
	y3 = w1*y2 + y1*w2 + z1*x2 - x1*z2;
	z3 = w1*z2 + z1*w2 + x1*y2 - y1*x2;

	return Quaternion<T>(w3, x3, y3, z3);
}


template<int M, int N>
class Matrix {
public:
	Matrix() {
		memset(values, 0, M*N * sizeof(double));
		for (int i = 0; i < N; i++) {
			values[i*N + i] = 1;
		}
	}
	Matrix<N, M> getTransposed() const {
		Matrix<N, M> res;
		for (int i=0; i<M; i++)
			for (int j = 0; j < N; j++) {
				res[j*M + i] = values[i*N+j];
			}
		return res;
	}
	void fromQuaternion(Quaternion<double> q) {
		double w, x, y, z;
		w = q.w; x = q.x; y = q.y; z = q.z;

		values[0] = w*w + x*x - y*y - z*z;
		values[1] = 2.0*x*y + 2.0*w*z;
		values[2] = 2.0*x*z - 2.0*y*w;
		values[3] = 2.0*x*y - 2.0*w*z;
		values[4] = w*w - x*x + y*y - z*z;
		values[5] = 2.0*y*z + 2.0*w*x;
		values[6] = 2.0*x*z + 2.0*w*y;
		values[7] = 2.0*y*z - 2.0*w*x;
		values[8] = w*w - x*x - y*y + z*z;
	}

	Quaternion<double> toQuaternion() const {
		double m00 = (*this)(0, 0);
		double m01 = (*this)(1, 0);
		double m02 = (*this)(2, 0);
		double m10 = (*this)(0, 1);
		double m11 = (*this)(1, 1);
		double m12 = (*this)(2, 1);
		double m20 = (*this)(0, 2);
		double m21 = (*this)(1, 2);
		double m22 = (*this)(2, 2);

		double tr = m00 + m11 + m22;
		double qw, qx, qy, qz;
		if (tr > 0) {
			double S = sqrt(tr + 1.0) * 2; // S=4*qw 
			qw = 0.25 * S;
			qx = (m21 - m12) / S;
			qy = (m02 - m20) / S;
			qz = (m10 - m01) / S;
		} else if ((m00 > m11)&(m00 > m22)) {
			double S = sqrt(1.0 + m00 - m11 - m22) * 2; // S=4*qx 
			qw = (m21 - m12) / S;
			qx = 0.25 * S;
			qy = (m01 + m10) / S;
			qz = (m02 + m20) / S;
		} else if (m11 > m22) {
			double S = sqrt(1.0 + m11 - m00 - m22) * 2; // S=4*qy
			qw = (m02 - m20) / S;
			qx = (m01 + m10) / S;
			qy = 0.25 * S;
			qz = (m12 + m21) / S;
		} else {
			double S = sqrt(1.0 + m22 - m00 - m11) * 2; // S=4*qz
			qw = (m10 - m01) / S;
			qx = (m02 + m20) / S;
			qy = (m12 + m21) / S;
			qz = 0.25 * S;
		}
		return Quaternion<double>(qw, qx, qy, qz);
	}
	double det() const {
		return (*this)(0, 0) * ((*this)(1, 1) * (*this)(2, 2) - (*this)(2, 1) * (*this)(1, 2)) -
			(*this)(0, 1) * ((*this)(1, 0) * (*this)(2, 2) - (*this)(1, 2) * (*this)(2, 0)) +
			(*this)(0, 2) * ((*this)(1, 0) * (*this)(2, 1) - (*this)(1, 1) * (*this)(2, 0));
	}
	double &operator()(int i, int j) {
		return values[i*N + j];
	}
	double operator()(int i, int j) const {
		return values[i*N + j];
	}
	Matrix<3,3> inverse() const {  // ONLY FOR 3x3

		double determinant = det();

		double invdet = 1 / determinant;

		Matrix<3,3> minv; // inverse of matrix m
		minv(0, 0) = ((*this)(1, 1) * (*this)(2, 2) - (*this)(2, 1) * (*this)(1, 2)) * invdet;
		minv(0, 1) = ((*this)(0, 2) * (*this)(2, 1) - (*this)(0, 1) * (*this)(2, 2)) * invdet;
		minv(0, 2) = ((*this)(0, 1) * (*this)(1, 2) - (*this)(0, 2) * (*this)(1, 1)) * invdet;
		minv(1, 0) = ((*this)(1, 2) * (*this)(2, 0) - (*this)(1, 0) * (*this)(2, 2)) * invdet;
		minv(1, 1) = ((*this)(0, 0) * (*this)(2, 2) - (*this)(0, 2) * (*this)(2, 0)) * invdet;
		minv(1, 2) = ((*this)(1, 0) * (*this)(0, 2) - (*this)(0, 0) * (*this)(1, 2)) * invdet;
		minv(2, 0) = ((*this)(1, 0) * (*this)(2, 1) - (*this)(2, 0) * (*this)(1, 1)) * invdet;
		minv(2, 1) = ((*this)(2, 0) * (*this)(0, 1) - (*this)(0, 0) * (*this)(2, 1)) * invdet;
		minv(2, 2) = ((*this)(0, 0) * (*this)(1, 1) - (*this)(1, 0) * (*this)(0, 1)) * invdet;
		return minv;
	}	
	double operator[](int i) const { return values[i]; }
	double& operator[](int i)  { return values[i]; }
	double values[M*N];
};

template<int M, int N> Matrix<M, N> operator+(const Matrix<M, N>& a, const Matrix<M, N>& b) {
	Matrix<M, N> result;
	for (int i = 0; i < M*N; i++) {
		result[i] = a[i] + b[i];
	}
	return result;
}

template<int M, int N, int O> Matrix<M, O> operator*(const Matrix<M, N>& a, const Matrix<N, O>& b) {
	Matrix<M, O> result;	
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < O; j++) {
			double s = 0;
			for (int k = 0; k < N; k++) {
				s += a[i*N+k]*b[k*O+j];
			}
			result.values[i*O + j] = s;
		}
	}
	return result;
}
template<int M, int N> Matrix<M, N> operator*(double a, const Matrix<M, N>& b) {
	Matrix<M, N> result;
	for (int i = 0; i < M*N; i++) {
		result[i] = a*b[i];
	}
	return result;
}

template<typename T>
Quaternion<T> Slerp(Quaternion<T> q1, Quaternion<T> q2, double t) {
	T w1, x1, y1, z1, w2, x2, y2, z2, w3, x3, y3, z3;
	Quaternion<T> q2New;
	double theta, mult1, mult2;

	w1 = q1.w; x1 = q1.x; y1 = q1.y; z1 = q1.z;
	w2 = q2.w; x2 = q2.x; y2 = q2.y; z2 = q2.z;

	// Reverse the sign of q2 if q1.q2 < 0.
	if (w1*w2 + x1*x2 + y1*y2 + z1*z2 < 0)
	{
		w2 = -w2; x2 = -x2; y2 = -y2; z2 = -z2;
	}

	theta = acos(w1*w2 + x1*x2 + y1*y2 + z1*z2);

	if (theta > 0.000001)
	{
		mult1 = sin((1 - t)*theta) / sin(theta);
		mult2 = sin(t*theta) / sin(theta);
	}

	// To avoid division by 0 and by very small numbers the approximation of sin(angle)
	// by angle is used when theta is small (0.000001 is chosen arbitrarily).
	else
	{
		mult1 = 1 - t;
		mult2 = t;
	}

	w3 = mult1*w1 + mult2*w2;
	x3 = mult1*x1 + mult2*x2;
	y3 = mult1*y1 + mult2*y2;
	z3 = mult1*z1 + mult2*z2;
	return Quaternion<T>(w3, x3, y3, z3);
}

template<int M, int N> Matrix<M, N> Slerp(const Matrix<M, N>& a, const Matrix<M, N>& b, double t) {

	Quaternion<double> qa = a.toQuaternion();
	Quaternion<double> qb = b.toQuaternion();
	Quaternion<double> qm = Slerp(qa, qb, t);

	Matrix<M, N> result;
	result.fromQuaternion(qm);
	return result;
}
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
	explicit Vector(const double* x) {
		memcpy(coord, x, 3 * sizeof(double));
	}
	explicit Vector(const float* x) {
		coord[0] = (double)x[0];
		coord[1] = (double)x[1];
		coord[2] = (double)x[2];
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
	std::string toRedValueStr() {
		std::ostringstream os;
		os << coord[0];
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




Vector min(const Vector& a, const Vector& b);
Vector max(const Vector& a, const Vector& b);
Vector pow(const Vector& a, const Vector& b);
Vector random_cos(const Vector &N);
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
	Camera(const Vector& position, const Vector &direction, const Vector& up): position(position), direction(direction), up(up), initial_position(position), initial_direction(direction), initial_up(up) {
		lenticular_max_angle = 35 * M_PI / 180.*0.25;
		lenticular_nb_images = 10;
		lenticular_pixel_width = 1;
		is_lenticular = false;
		isArray = false;
		current_viewX = 0; 
		current_viewY = 0; 
		nbviewX = 1;
		nbviewY = 1;
	};

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
		double k = W / (2 * tan(fov / 2));
		Vector camera_right = cross(direction, up);		

		Vector direction_vec;
		Vector C1;
		if (is_lenticular) {
			double L = focus_distance*tan(lenticular_max_angle / 2) / (lenticular_nb_images / 2.0);

			int offset = -((j / lenticular_pixel_width) % lenticular_nb_images - lenticular_nb_images / 2);

			Vector P = position + focus_distance*Vector(0, 0, 1);
			C1 = position + offset*L * camera_right;
			Vector v1 = (P - C1).getNormalized(); // direction we want central in offseted camera 
			Vector PprojCam1 = (k / dot(v1, direction)) * v1 + C1;
			double pixProjCam1_j = PprojCam1[0] + W / 2 - 0.5;
			double pixProjCam1_i = PprojCam1[1] + H / 2 - 0.5;  // new offset : this pixel should be the central pixel in camera 


			//Vector direction_vec(j - W / 2 + 0.0 + dx_sensor, i - H / 2 + 0.5 + dy_sensor, k);			
			direction_vec = Vector((j - pixProjCam1_j) + dx_sensor, (i - pixProjCam1_i) + dy_sensor, k);
		} else {
			C1 = position;
			direction_vec = Vector(j - W / 2 + 0.5 + dx_sensor, i - H / 2 + 0.5 + dy_sensor, k);

		}
		direction_vec.normalize();
		direction_vec = camera_right * direction_vec[0] + up*direction_vec[1] + direction*direction_vec[2];
		Vector destination = C1 + focus_distance * direction_vec;
		Vector new_origin = C1 + dx_aperture*camera_right + dy_aperture*up;
		return Ray(new_origin, (destination - new_origin).getNormalized(), time);
	}

	Vector position, direction, up;
	double fov;             // field of view 
	double focus_distance; // distance mise au point
	double aperture;       // ouverture, pour depth of field (faible valeur = net partout)
	double lenticular_max_angle;
	int lenticular_nb_images, lenticular_pixel_width;
	bool is_lenticular;
	int current_viewX, current_viewY, nbviewX, nbviewY;
	bool isArray;

private:
	Vector initial_position, initial_direction, initial_up;
};


