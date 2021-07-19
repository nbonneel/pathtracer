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

#ifndef M_TWO_PI
#define M_TWO_PI 6.28318530718
#endif

/*static __declspec(thread) std::default_random_engine engine;
static __declspec(thread) std::uniform_real_distribution<double> uniform(0, 1);*/

/*static thread_local std::minstd_rand  engine[64];
static thread_local std::uniform_real_distribution<double> uniform(0, 1);
static thread_local std::uniform_real_distribution<float> uniformf(0, 1);*/

#include "pcg_random.hpp"

static thread_local pcg32  engine[64];
//static thread_local std::uniform_real_distribution<double> uniform(0, 1);
static  std::uniform_real_distribution<float> uniformf(0, 1);

template<typename T> class VectorT;
typedef VectorT<double> Vectord;
typedef VectorT<float> Vectorf;

typedef Vectorf Vector;

template<int M, int N, typename T> class Matrix;
typedef Matrix<3,3,double> Matrix33d;
typedef Matrix<3, 3, float> Matrix33f;

typedef Matrix33f Matrix33;

template<typename T> inline T sqr(T x) { return x * x; }


/*Vector operator+(const Vector& a, const Vector &b);
Vector operator-(const Vector& a, const Vector &b);
Vector operator*(double a, const Vector &b);
Vector operator*(const Vector &b, double a);
Vector operator*(const Vector &a, const Vector &b);
Vector operator/(const Vector& a, double b);
Vector operator/(const Vector& a, const Vector &b);
Vector operator-(const Vector& a);
double dot(const Vector&a, const Vector& b);
Vector cross(const Vector&a, const Vector& b);*/

template<typename T>
class Quaternion
{
public:
	Quaternion() { }
	Quaternion(T wVal, T xVal, T yVal, T zVal)
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


template<int M, int N, typename T>
class Matrix {
public:
	Matrix() {
		memset(values, 0, M*N * sizeof(T));
		for (int i = 0; i < N; i++) {
			values[i*N + i] = 1;
		}
	}
	Matrix<N, M, T> getTransposed() const {
		Matrix<N, M, T> res;
		for (int i=0; i<M; i++)
			for (int j = 0; j < N; j++) {
				res[j*M + i] = values[i*N+j];
			}
		return res;
	}
	void fromQuaternion(Quaternion<T> q) {
		T w, x, y, z;
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

	Quaternion<T> toQuaternion() const {
		T m00 = (*this)(0, 0);
		T m01 = (*this)(1, 0);
		T m02 = (*this)(2, 0);
		T m10 = (*this)(0, 1);
		T m11 = (*this)(1, 1);
		T m12 = (*this)(2, 1);
		T m20 = (*this)(0, 2);
		T m21 = (*this)(1, 2);
		T m22 = (*this)(2, 2);

		T tr = m00 + m11 + m22;
		T qw, qx, qy, qz;
		if (tr > 0) {
			T S = sqrt(tr + 1.0) * 2; // S=4*qw 
			qw = 0.25 * S;
			qx = (m21 - m12) / S;
			qy = (m02 - m20) / S;
			qz = (m10 - m01) / S;
		} else if ((m00 > m11)&(m00 > m22)) {
			T S = sqrt(1.0 + m00 - m11 - m22) * 2; // S=4*qx 
			qw = (m21 - m12) / S;
			qx = 0.25 * S;
			qy = (m01 + m10) / S;
			qz = (m02 + m20) / S;
		} else if (m11 > m22) {
			T S = sqrt(1.0 + m11 - m00 - m22) * 2; // S=4*qy
			qw = (m02 - m20) / S;
			qx = (m01 + m10) / S;
			qy = 0.25 * S;
			qz = (m12 + m21) / S;
		} else {
			T S = sqrt(1.0 + m22 - m00 - m11) * 2; // S=4*qz
			qw = (m10 - m01) / S;
			qx = (m02 + m20) / S;
			qy = (m12 + m21) / S;
			qz = 0.25 * S;
		}
		return Quaternion<T>(qw, qx, qy, qz);
	}
	T det() const {
		return (*this)(0, 0) * ((*this)(1, 1) * (*this)(2, 2) - (*this)(2, 1) * (*this)(1, 2)) -
			(*this)(0, 1) * ((*this)(1, 0) * (*this)(2, 2) - (*this)(1, 2) * (*this)(2, 0)) +
			(*this)(0, 2) * ((*this)(1, 0) * (*this)(2, 1) - (*this)(1, 1) * (*this)(2, 0));
	}
	T &operator()(int i, int j) {
		return values[i*N + j];
	}
	T operator()(int i, int j) const {
		return values[i*N + j];
	}
	Matrix<3,3, T> inverse() const {  // ONLY FOR 3x3

		T determinant = det();

		T invdet = 1 / determinant;

		Matrix<3,3, T> minv; // inverse of matrix m
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
	T operator[](int i) const { return values[i]; }
	T& operator[](int i)  { return values[i]; }
	T values[M*N];
};

template<int M, int N, typename T> Matrix<M, N, T> operator+(const Matrix<M, N, T>& a, const Matrix<M, N, T>& b) {
	Matrix<M, N, T> result;
	for (int i = 0; i < M*N; i++) {
		result[i] = a[i] + b[i];
	}
	return result;
}

template<int M, int N, int O, typename T> Matrix<M, O, T> operator*(const Matrix<M, N, T>& a, const Matrix<N, O, T>& b) {
	Matrix<M, O, T> result;	
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < O; j++) {
			T s = 0;
			for (int k = 0; k < N; k++) {
				s += a[i*N+k]*b[k*O+j];
			}
			result.values[i*O + j] = s;
		}
	}
	return result;
}
template<int M, int N, typename T> Matrix<M, N, T> operator*(T a, const Matrix<M, N, T>& b) {
	Matrix<M, N, T> result;
	for (int i = 0; i < M*N; i++) {
		result[i] = a*b[i];
	}
	return result;
}

template<typename T>
Quaternion<T> Slerp(Quaternion<T> q1, Quaternion<T> q2, T t) {
	T w1, x1, y1, z1, w2, x2, y2, z2, w3, x3, y3, z3;
	Quaternion<T> q2New;
	T theta, mult1, mult2;

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

template<int M, int N, typename T> Matrix<M, N, T> Slerp(const Matrix<M, N, T>& a, const Matrix<M, N, T>& b, T t) {

	Quaternion<T> qa = a.toQuaternion();
	Quaternion<T> qb = b.toQuaternion();
	Quaternion<T> qm = Slerp(qa, qb, t);

	Matrix<M, N, T> result;
	result.fromQuaternion(qm);
	return result;
}
static inline Matrix<3, 3, float> createRotationMatrixX(float angleX) {
	Matrix<3, 3, float> R;
	R[0] = cos(angleX);
	R[1] = -sin(angleX);
	R[3] = sin(angleX);
	R[4] = cos(angleX);
	return R;
}
static inline Matrix<3, 3, float> createRotationMatrixY(float angleY) {
	Matrix<3, 3, float> R;
	R[0] = cos(angleY);
	R[2] = sin(angleY);
	R[6] = -sin(angleY);
	R[8] = cos(angleY);
	return R;
}
static inline Matrix<3, 3, float> createRotationMatrixZ(float angleZ) {
	Matrix<3, 3, float> R;
	R[4] = cos(angleZ);
	R[5] = -sin(angleZ);
	R[7] = sin(angleZ);
	R[8] = cos(angleZ);
	return R;
}
static inline float invSqRoot(float n) {

	const float threehalfs = 1.5F;
	float y = n;

	long i = *(long *)&y;

	i = 0x5f3759df - (i >> 1);
	//i = 0x5f37e75a - (i >> 1);  // arxiv 1603.04483 ; 0x5f37e75a or 0x5f37add5
	y = *(float *)&i;

	y = y * (threehalfs - ((n * 0.5F) * y * y));
	y = y * (threehalfs - ((n * 0.5F) * y * y));

	return y;
}

static inline double fastPow(double a, double b) {
	union {
		double d;
		int x[2];
	} u = { a };
	u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
	u.x[0] = 0;
	return u.d;
}
//https://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
inline double fastPrecisePow(double a, double b) {
	// calculate approximation with fraction of the exponent
	int e = (int)b;
	union {
		double d;
		int x[2];
	} u = { a };
	u.x[1] = (int)((b - e) * (u.x[1] - 1072632447) + 1072632447);
	u.x[0] = 0;

	// exponentiation by squaring with the exponent's integer part
	// double r = u.d makes everything much slower, not sure why
	double r = 1.0;
	while (e) {
		if (e & 1) {
			r *= a;
		}
		a *= a;
		e >>= 1;
	}

	return r * u.d;
}

template<typename T>
class VectorT {
public:
	explicit VectorT<T>(T x = 0, T y = 0, T z = 0) {
		coord[0] = x;
		coord[1] = y;
		coord[2] = z;
	}
	explicit VectorT<T>(const double* x) {
		//memcpy(coord, x, 3 * sizeof(double));
		coord[0] = (T)x[0];
		coord[1] = (T)x[1];
		coord[2] = (T)x[2];
	}
	explicit VectorT<T>(const float* x) {
		coord[0] = (T)x[0];
		coord[1] = (T)x[1];
		coord[2] = (T)x[2];
	}
	const T& operator[](int i) const { return coord[i]; }
	T& operator[](int i) { return coord[i]; }

	T getNorm2() const {
		return coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2];
	}
	void normalize() {
		T norm = sqrt(getNorm2());
		coord[0] /= norm;
		coord[1] /= norm;
		coord[2] /= norm;
	}
	void fast_normalize() {
		T norm2 = getNorm2();
		T invnorm = invSqRoot(norm2);  // 1. / sqrt(norm2); // 
		coord[0] *= invnorm;
		coord[1] *= invnorm;
		coord[2] *= invnorm;
	}
	VectorT<T> getNormalized() const {
		VectorT<T> result(*this);
		result.normalize();
		return result;
	}
	VectorT<T> reflect(const VectorT<T> &N) const {
		VectorT<T> result = *this - T(2.)*dot(*this, N)*N;
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

	VectorT<T>& operator+=(const VectorT<T>& b) {
		coord[0] += b[0];
		coord[1] += b[1];
		coord[2] += b[2];
		return *this;
	}
	VectorT<T>& operator-=(const VectorT<T>& b) {
		coord[0] -= b[0];
		coord[1] -= b[1];
		coord[2] -= b[2];
		return *this;
	}
	VectorT<T>& operator*=(T b) {
		coord[0] *= b;
		coord[1] *= b;
		coord[2] *= b;
		return *this;
	}
	VectorT<T>& operator*=(const VectorT<T> &b) {
		coord[0] *= b[0];
		coord[1] *= b[1];
		coord[2] *= b[2];
		return *this;
	}
	VectorT<T>& operator/=(T b) {
		coord[0] /= b;
		coord[1] /= b;
		coord[2] /= b;
		return *this;
	}

private:
	T coord[3];
};

template<typename T>
static inline  VectorT<T> operator*(const Matrix<3, 3, T> &mat, const VectorT<T> &b) {

	VectorT<T> res(0, 0, 0);
	for (int i = 0; i < 3; i++) {
		T v = 0;
		for (int j = 0; j < 3; j++) {
			v += mat[i*3 + j] * b[j];
		}
		res[i] = v;
	}
	return res;
}

template<typename T>
static inline VectorT<T> operator*(const Matrix<3, 4, T> &mat, const VectorT<T> &b) {
	VectorT<T> res(0, 0, 0);
	for (int i = 0; i < 3; i++) {
		T v = 0;
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




/*Vector min(const Vector& a, const Vector& b);
Vector max(const Vector& a, const Vector& b);
Vector pow(const Vector& a, const Vector& b);
Vector fastPow(const Vector& a, const Vector& b);
Vector fastPrecisePow(const Vector& a, const Vector& b);
Vector exp(const Vector& a);
Vector sqrt(const Vector& a);
Vector sqr(const Vector& a);
Vector random_cos(const Vector &N);
Vector random_cos(const Vector &N, double r1, double r2);
Vector random_uniform_sphere();
Vector random_uniform_ball();
Vector random_uniform_hemiball(const Vector &N);
Vector getTangent(const Vector& N);
Vector boxMuller();

Vector rotate_dir(const Vector&v, const Vector &angles);
Vector inverse_rotate_dir(const Vector&v, const Vector &angles); // performs the rotation by -angle, but in the reversed order */



template<typename T>
VectorT<T> operator+(const VectorT<T>& a, const VectorT<T> &b) {
	return VectorT<T>(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
template<typename T>
VectorT<T> operator-(const VectorT<T>& a, const VectorT<T> &b) {
	return VectorT<T>(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
template<typename T>
VectorT<T> operator*(T a, const VectorT<T> &b) {
	return VectorT<T>(a*b[0], a*b[1], a*b[2]);
}
template<typename T>
VectorT<T> operator*(const VectorT<T> &b, T a) {
	return VectorT<T>(a*b[0], a*b[1], a*b[2]);
}
template<typename T>
VectorT<T> operator/(const VectorT<T>& a, T b) {
	return VectorT<T>(a[0] / b, a[1] / b, a[2] / b);
}
template<typename T>
VectorT<T> operator/(const VectorT<T>& a, const VectorT<T> &b) {
	return VectorT<T>(a[0] / b[0], a[1] / b[1], a[2] / b[2]);
}
template<typename T>
VectorT<T> operator-(const VectorT<T>& a) {
	return VectorT<T>(-a[0], -a[1], -a[2]);
}
template<typename T>
VectorT<T> operator*(const VectorT<T> &a, const VectorT<T> &b) {
	return VectorT<T>(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}
template<typename T, typename T2>
VectorT<T> min(const VectorT<T>& a, const VectorT<T2>& b) {
	return VectorT<T>(std::min(a[0], (T)b[0]), std::min(a[1],(T) b[1]), std::min(a[2], (T)b[2]));
}
template<typename T, typename T2>
VectorT<T> max(const VectorT<T>& a, const VectorT<T2>& b) {
	return VectorT<T>(std::max(a[0], (T)b[0]), std::max(a[1], (T)b[1]), std::max(a[2],(T) b[2]));
}
template<typename T>
VectorT<T> pow(const VectorT<T>& a, const VectorT<T>& b) {
	return VectorT<T>(std::pow(a[0], b[0]), std::pow(a[1], b[1]), std::pow(a[2], b[2]));
}
template<typename T>
VectorT<T> fastPow(const VectorT<T>& a, const VectorT<T>& b) {
	return VectorT<T>(fastPow(a[0], b[0]), fastPow(a[1], b[1]), fastPow(a[2], b[2]));
}
template<typename T>
VectorT<T> fastPrecisePow(const VectorT<T>& a, const VectorT<T>& b) {
	return VectorT<T>(fastPrecisePow(a[0], b[0]), fastPrecisePow(a[1], b[1]), fastPrecisePow(a[2], b[2]));
}
template<typename T>
T dot(const VectorT<T>&a, const VectorT<T>& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
template<typename T>
VectorT<T> exp(const VectorT<T>& a) {
	return VectorT<T>(exp(a[0]), exp(a[1]), exp(a[2]));
}
template<typename T>
VectorT<T> sqrt(const VectorT<T>& a) {
	return VectorT<T>(sqrt(a[0]), sqrt(a[1]), sqrt(a[2]));
}
template<typename T>
VectorT<T> sqr(const VectorT<T>& a) {
	return VectorT<T>(a[0] * a[0], a[1] * a[1], a[2] * a[2]);
}

template<typename T>
VectorT<T> cross(const VectorT<T>&a, const VectorT<T>& b) {
	return VectorT<T>(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

template<typename T>
VectorT<T> getTangent(const VectorT<T>& N) {
	VectorT<T> tangent1;
	VectorT<T> absN(abs(N[0]), abs(N[1]), abs(N[2]));
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

template<typename T>
VectorT<T> random_cos(const VectorT<T> &N, T r1, T r2) {
	T sr2 = sqrt(T(1) - r2);

	VectorT<T> direction_aleatoire_repere_local(cos(T(2. * M_PI)*r1)*sr2, sin(T(2. * M_PI)*r1)*sr2, sqrt(r2));
	VectorT<T> tangent1 = getTangent(N);
	VectorT<T> tangent2 = cross(tangent1, N);
	return direction_aleatoire_repere_local[2] * N + direction_aleatoire_repere_local[0] * tangent1 + direction_aleatoire_repere_local[1] * tangent2;
}

template<typename T>
VectorT<T> random_cos(const VectorT<T> &N) {
	int threadid = omp_get_thread_num();
	/*float r1 = uniformf(engine[threadid]);
	float r2 = uniformf(engine[threadid]);*/
	const T invmax = T(1) / engine[threadid].max();
	T r1 = engine[threadid]()*invmax;
	T r2 = engine[threadid]()*invmax;
	return random_cos(N, r1, r2);
}



template<typename T>
VectorT<T> random_uniform_sphere() {
	int threadid = omp_get_thread_num();
	T invmax = T(1) / engine[threadid].max();
	T r1 = engine[threadid]()*invmax; // uniform(engine[threadid]);
	T r2 = engine[threadid]()*invmax; // uniform(engine[threadid]);
	VectorT<T> result;
	result[0] = T(2.)*cos(T(2.*M_PI)*r1)*sqrt(r2*(1 - r2));
	result[1] = T(2.)*sin(T(2.*M_PI)*r1)*sqrt(r2*(1 - r2));
	result[2] = T(1) - T(2) * r2;
	return result;
}

template<typename T>
VectorT<T> random_uniform_hemisphere(const VectorT<T>& N) {
	int threadid = omp_get_thread_num();
	T invmax = T(1) / engine[threadid].max();
	T r1 = engine[threadid]()*invmax; // uniform(engine[threadid]);
	T r2 = engine[threadid]()*invmax; // uniform(engine[threadid]);
	VectorT<T> resultlocal;
	resultlocal[0] = cos(T(2.*M_PI)*r1)*sqrt(1 - r2 * r2);
	resultlocal[1] = sin(T(2.*M_PI)*r1)*sqrt(1 - r2 * r2);
	resultlocal[2] = r2;
	VectorT<T> tangent1 = getTangent(N);
	VectorT<T> tangent2 = cross(tangent1, N);
	return resultlocal[2] * N + resultlocal[0] * tangent1 + resultlocal[1] * tangent2;
}

template<typename T>
VectorT<T> random_uniform_ball() {
	int threadid = omp_get_thread_num();
	T invmax = T(1) / engine[threadid].max();
	VectorT<T> result = std::pow(engine[threadid]()*invmax, T(1. / 3.))*random_uniform_sphere();
	return result;
}
template<typename T>
VectorT<T> random_uniform_hemiball(const VectorT<T> &N) {
	int threadid = omp_get_thread_num();
	T invmax = T(1) / engine[threadid].max();
	VectorT<T> result = std::pow(engine[threadid]()*invmax, T(1. / 3.))*random_uniform_hemisphere(N);
	return result;
}
template<typename T>
VectorT<T> boxMuller() { // returns radius in third component
	int threadid = omp_get_thread_num();
	T invmax = T(1) / engine[threadid].max();
	T r1 = engine[threadid]()*invmax;
	T r2 = engine[threadid]()*invmax;
	T s1 = sqrt(T(-2) * log(r1));
	T s2 = T(2 * M_PI)*r2;
	return VectorT<T>(s1*cos(s2), s1*sin(s2), s1);
}



template<typename T>
VectorT<T> rotate_dir(const VectorT<T>&v, const VectorT<T> &angles) {

	VectorT<T> sinangles(sin(angles[0]), sin(angles[1]), sin(angles[2]));
	VectorT<T> cosangles(cos(angles[0]), cos(angles[1]), cos(angles[2]));

	VectorT<T> v1;
	v1[0] = v[0];
	v1[1] = cosangles[0] * v[1] - sinangles[0] * v[2];
	v1[2] = sinangles[0] * v[1] + cosangles[0] * v[2];

	VectorT<T> v2;
	v2[0] = sinangles[1] * v1[2] + cosangles[1] * v1[0];
	v2[1] = v1[1];
	v2[2] = cosangles[1] * v1[2] - sinangles[1] * v1[0];

	VectorT<T> v3;
	v3[0] = cosangles[2] * v2[0] - sinangles[2] * v2[1];
	v3[1] = sinangles[2] * v2[0] + cosangles[2] * v2[1];
	v3[2] = v2[2];

	return v3;
}

template<typename T>
VectorT<T> inverse_rotate_dir(const VectorT<T>&v, const VectorT<T> &angles) {

	VectorT<T> sinangles(sin(-angles[0]), sin(-angles[1]), sin(-angles[2]));
	VectorT<T> cosangles(cos(-angles[0]), cos(-angles[1]), cos(-angles[2]));

	VectorT<T> v1;
	v1[0] = cosangles[2] * v[0] - sinangles[2] * v[1];
	v1[1] = sinangles[2] * v[0] + cosangles[2] * v[1];
	v1[2] = v[2];

	VectorT<T> v2;
	v2[0] = sinangles[1] * v1[2] + cosangles[1] * v1[0];
	v2[1] = v1[1];
	v2[2] = cosangles[1] * v1[2] - sinangles[1] * v1[0];

	VectorT<T> v3;
	v3[0] = v2[0];
	v3[1] = cosangles[0] * v2[1] - sinangles[0] * v2[2];
	v3[2] = sinangles[0] * v2[1] + cosangles[0] * v2[2];


	return v3;

}

Vectorf toVecf(const Vectord& v);
Vectord toVecd(const Vectorf& v);

class Ray {
public:
	Ray() {};
	Ray(const Vector& o, const Vector& d, float time) : origin(o), direction(d), time(time) {};
	Vector origin, direction;
	float time;
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

	void translate(const Vector& translation, float time) {
		position = initial_position + time*translation;
	}

	void rotate(float angle_x, float angle_y, float time) {
		float cur_angle_x = time*angle_x;
		float cur_angle_y = time*angle_y;

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

	void rotateAroundRight(float angle) {

		Vector dirtmp, uptmp;
		dirtmp =  sin(angle)*up  + cos(angle)*direction;		
		uptmp = cos(angle)*up - sin(angle)*direction;

		direction = dirtmp;
		up = uptmp;
	}
	void rotateAroundUp(float angle) {

		Vector right = cross(up, direction);
		direction = -sin(angle)*right + cos(angle)*direction;

	}

	void rotateAroundAxes(float angle_x, float angle_y, float time) {
		float cur_angle_x = time*angle_x;
		float cur_angle_y = time*angle_y;

		rotateAroundUp(cur_angle_x);
		rotateAroundRight(cur_angle_y);
	}
	

	Ray generateDirection(float init_t, int i, int j, float time, float dx_sensor, float dy_sensor, float dx_aperture, float dy_aperture, int W, int H) {
		float k = W / (2 * tan(fov / 2));
		Vector camera_right = cross(direction, up);		

		Vector direction_vec;
		Vector C1;
		if (is_lenticular) {
			float L = focus_distance*tan(lenticular_max_angle / 2) / (lenticular_nb_images / 2.0);

			int offset = -((j / lenticular_pixel_width) % lenticular_nb_images - lenticular_nb_images / 2);

			Vector P = position + focus_distance*Vector(0, 0, 1);
			C1 = position + offset*L * camera_right;
			Vector v1 = (P - C1).getNormalized(); // direction we want central in offseted camera 
			Vector PprojCam1 = (k / dot(v1, direction)) * v1 + C1;
			float pixProjCam1_j = PprojCam1[0] + W / 2 - 0.5;
			float pixProjCam1_i = PprojCam1[1] + H / 2 - 0.5;  // new offset : this pixel should be the central pixel in camera 


			//Vector direction_vec(j - W / 2 + 0.0 + dx_sensor, i - H / 2 + 0.5 + dy_sensor, k);			
			direction_vec = Vector((j - pixProjCam1_j) + dx_sensor, (i - pixProjCam1_i) + dy_sensor, k);
		} else {
			C1 = position;
			direction_vec = Vector(j - W / 2 + 0.5 + dx_sensor, i - H / 2 + 0.5 + dy_sensor, k);

		}
		//direction_vec.normalize();
		direction_vec.normalize();
		direction_vec = camera_right * direction_vec[0] + up*direction_vec[1] + direction*direction_vec[2];
		Vector destination = C1 + focus_distance/std::abs(dot(direction_vec, direction)) * direction_vec;
		Vector new_origin = C1 + dx_aperture*camera_right + dy_aperture*up;
		Vector new_direction = (destination - new_origin).getNormalized();
		return Ray(new_origin + init_t*new_direction/dot(new_direction, direction), new_direction, time);
	}

	Vector position, direction, up;
	float fov;             // field of view 
	float focus_distance; // distance mise au point
	float aperture;       // ouverture, pour depth of field (faible valeur = net partout)
	float lenticular_max_angle;
	float maxSpacingX, maxSpacingY;
	int lenticular_nb_images, lenticular_pixel_width;
	bool is_lenticular;
	int current_viewX, current_viewY, nbviewX, nbviewY;
	bool isArray;

private:
	Vector initial_position, initial_direction, initial_up;
};


