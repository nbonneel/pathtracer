#include "Vector.h"
#include <cmath>
#include <algorithm>
#include <omp.h>

Vectorf toVecf(const Vectord& v) {
	return Vectorf(v[0], v[1], v[2]);
}
Vectord toVecd(const Vectorf& v) {
	return Vectord(v[0], v[1], v[2]);
}