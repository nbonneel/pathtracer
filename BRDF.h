#pragma once

#include "Vector.h"
#include <omp.h>
#include "MERLBRDFRead.h"

class MaterialValues {
public:
	MaterialValues() {
		shadingN[0] = 0;  shadingN[1] = 1; shadingN[2] = 0;
		Kd[0] = 0.5; Kd[1] = 0.5; Kd[2] = 0.5;
		Ne[0] = 100; Ne[1] = 100; Ne[2] = 100;
		Ks[0] = 0.; Ks[1] = 0.; Ks[2] = 0.;
		Ke[0] = 0.; Ke[1] = 0.; Ke[2] = 0.;
		Ksub[0] = 0.; Ksub[1] = 0.; Ksub[2] = 0.;				
	}
	Vector shadingN, Kd, Ks, Ne, Ke, Ksub;
	bool transp;
	double refr_index;
};

class BRDF {
public:
	BRDF() {};
	virtual Vector sample(const MaterialValues& mat, const Vector& wi, const Vector& N, double &pdf, double r1, double r2) const = 0;
	virtual Vector eval(const MaterialValues& mat, const Vector& wi, const Vector& wo, const Vector& N) const = 0;

	Vector sample(const MaterialValues& mat, const Vector& wo, const Vector& N, double &pdf) const {  // performs MIS Kd-Ks
		int threadid = omp_get_thread_num();
		double invmax = 1.f / engine[threadid].max();
		return sample(mat, wo, N, pdf, engine[threadid]()*invmax, engine[threadid]()*invmax);
	}
};



class PhongBRDF : public BRDF {
public:
	PhongBRDF() {};

	static Vector random_Phong(const Vector &R, double phong_exponent, double r1, double r2) {
		double facteur = sqrt(1 - std::pow(r2, 2. / (phong_exponent + 1)));
		//double facteur = sqrt(1 - fastPrecisePow(r2, 2. / (phong_exponent + 1)));
		Vector direction_aleatoire_repere_local(cos(2 * M_PI*r1)*facteur, sin(2 * M_PI*r1)*facteur, std::pow(r2, 1. / (phong_exponent + 1)));
		//Vector aleatoire(uniform(engine) - 0.5, uniform(engine) - 0.5, uniform(engine) - 0.5);
		//Vector tangent1 = cross(R, aleatoire); tangent1.normalize();
		Vector tangent1;
		Vector absR(abs(R[0]), abs(R[1]), abs(R[2]));
		if (absR[0] <= absR[1] && absR[0] <= absR[2]) {
			tangent1 = Vector(0, -R[2], R[1]);
		} else
			if (absR[1] <= absR[0] && absR[1] <= absR[2]) {
				tangent1 = Vector(-R[2], 0, R[0]);
			} else
				tangent1 = Vector(-R[1], R[0], 0);
			tangent1.normalize();

			Vector tangent2 = cross(tangent1, R);

			return direction_aleatoire_repere_local[2] * R + direction_aleatoire_repere_local[0] * tangent1 + direction_aleatoire_repere_local[1] * tangent2;
	}

	Vector sample(const MaterialValues& mat, const Vector& wo, const Vector& N, double &pdf, double r1, double r2) const {  // performs MIS Kd-Ks
		int threadid = omp_get_thread_num();
		double avgNe = (mat.Ne[0] + mat.Ne[1] + mat.Ne[2]) / 3.;
		Vector direction_aleatoire;
		double p = 1 - (mat.Ks[0] + mat.Ks[1] + mat.Ks[2]) / 3.;
		/*if (s.objects[sphere_id]->ghost) {
			p = std::max(0.2, p);
		}*/
		bool sample_diffuse;
		Vector R = (-wo).reflect(N); // reflect takes a ray going towards a surface and reflects it outwards of it
		if (engine[threadid]() / (double)engine[threadid].max() < p) {
			sample_diffuse = true;
			direction_aleatoire = random_cos(N, r1, r2);
		} else {
			sample_diffuse = false;
			direction_aleatoire = random_Phong(R, avgNe, r1, r2);
		}

		double proba_phong = (avgNe + 1) / (2.*M_PI) * pow(dot(R, direction_aleatoire), avgNe);
		double proba_globale = p * dot(N, direction_aleatoire) / (M_PI)+(1. - p)*proba_phong;
		pdf = proba_globale;

		return direction_aleatoire;
	}

	Vector eval(const MaterialValues& mat, const Vector& wi, const Vector& wo, const Vector& N) const {
		Vector reflechi = (-wo).reflect(N);
		double d = dot(reflechi, wi);
		if (d < 0) return mat.Kd / M_PI;
		//Vector lobe = pow(Vector(d, d, d), mat.Ne) * ((mat.Ne + Vector(2., 2., 2.)) / (2.*M_PI));
		Vector lobe = Vector(pow(d, mat.Ne[0])*(mat.Ne[0]+2.)/M_TWO_PI, pow(d, mat.Ne[1])*(mat.Ne[1] + 2.) / M_TWO_PI, pow(d, mat.Ne[2])*(mat.Ne[2] + 2.) / M_TWO_PI);
		//Vector lobe = fastPrecisePow(Vector(d, d, d), mat.Ne) * ((mat.Ne + Vector(2., 2., 2.)) / (2.*M_PI));
		return mat.Kd / M_PI + lobe * mat.Ks;
	}
};

class LambertBRDF : public BRDF {
public:
	LambertBRDF() {};

	Vector sample(const MaterialValues& mat, const Vector& wo, const Vector& N, double &pdf, double r1, double r2) const {  // performs MIS Kd-Ks		
		Vector direction_aleatoire = random_cos(N, r1, r2);
		pdf = dot(N, direction_aleatoire) / (M_PI);
		return direction_aleatoire;
	}
	Vector eval(const MaterialValues& mat, const Vector& wi, const Vector& wo, const Vector& N) const {
		return mat.Kd / M_PI;
	}
	Vector Kd;
};


class TitopoBRDF : public BRDF { //our internal thetai,thetao, phio format
public:
	TitopoBRDF(std::string filename, int Nthetai, int Nthetao, int Nphid) :Nthetai(Nthetai), Nthetao(Nthetao), Nphid(Nphid) {
		data.resize(Nthetai * Nthetao * Nphid * 3);
		FILE* f = fopen(filename.c_str(), "rb");
		fread(&data[0], sizeof(float), Nthetai * Nthetao * Nphid * 3, f); // half
		fclose(f);

	}

	Vector sample(const MaterialValues& mat, const Vector& wi, const Vector& N, double &pdf, double r1, double r2) const {
		Vector d = random_cos(N, r1, r2);
		pdf = dot(N, d) / (M_PI);
		return d;
	}
	Vector eval(const MaterialValues& mat, const Vector& wi, const Vector& wo, const Vector& N) const {

		/*Vector sumVal(0., 0., 0.);
		for (int i = 0; i < Nthetai*Nthetao*Nphio; i++) {
			sumVal += Vector(&data[i * 3]);
		}
		std::cout << sumVal[0] << " " << sumVal[1] << " " << sumVal[2] << std::endl;*/

		Vector tangent1;
		Vector absN(abs(N[0]), abs(N[1]), abs(N[2]));
		if (absN[0] <= absN[1] && absN[0] <= absN[2]) {
			tangent1 = Vector(0, -N[2], N[1]);
		} else {
			if (absN[1] <= absN[0] && absN[1] <= absN[2]) {
				tangent1 = Vector(-N[2], 0, N[0]);
			} else
				tangent1 = Vector(-N[1], N[0], 0);
		}
		tangent1.normalize();
		Vector tangent2 = cross(tangent1, N);

		Vector wilocal = Vector(dot(wi, tangent1), dot(wi, tangent2), dot(wi, N));
		Vector wolocal = Vector(dot(wo, tangent1), dot(wo, tangent2), dot(wo, N));

		double thetai = acos(wilocal[2]);
		if (thetai >= M_PI / 2) return Vector(0., 0., 0.);
		double thetao = acos(wolocal[2]);
		if (thetao >= M_PI / 2) return Vector(0., 0., 0.);
		double phid = atan2(wolocal[1], wolocal[0]) - atan2(wilocal[1], wilocal[0]); // in -pi..pi 
		if (phid < 0) phid += 2 * M_PI;  // => 0..2pi
		if (phid < 0) phid += 2 * M_PI;  // => 0..2pi
		if (phid >= 2 * M_PI) phid -= 2 * M_PI;  // => 0..2pi
		if (phid >= 2 * M_PI) phid -= 2 * M_PI;  // => 0..2pi
		int idthetai = (int)(thetai / (M_PI / 2.)*Nthetai);
		int idthetao = (int)(thetao / (M_PI / 2.)*Nthetao);
		int idphid = (int)(phid / (M_PI * 2.)*Nphid);
		int idthetai2 = ((idthetai < Nthetai - 1) ? (idthetai + 1) : idthetai);
		int idthetao2 = ((idthetao < Nthetao - 1) ? (idthetao + 1) : idthetao);
		int idphid2 = ((idphid < Nphid - 1) ? (idphid + 1) : idphid);
		double fthetai = thetai / (M_PI / 2.)*Nthetai - idthetai;
		double fthetao = thetao / (M_PI / 2.)*Nthetao - idthetao;
		double fphid = phid / (M_PI * 2.)*Nphid - idphid;
		Vector val000 = Vector(&data[(idthetai*Nthetao*Nphid + idthetao * Nphid + idphid) * 3]);
		Vector val001 = Vector(&data[(idthetai*Nthetao*Nphid + idthetao * Nphid + idphid2) * 3]);
		Vector val010 = Vector(&data[(idthetai*Nthetao*Nphid + idthetao2 * Nphid + idphid) * 3]);
		Vector val100 = Vector(&data[(idthetai2*Nthetao*Nphid + idthetao * Nphid + idphid) * 3]);
		Vector val101 = Vector(&data[(idthetai2*Nthetao*Nphid + idthetao * Nphid + idphid2) * 3]);
		Vector val110 = Vector(&data[(idthetai2*Nthetao*Nphid + idthetao2 * Nphid + idphid) * 3]);
		Vector val011 = Vector(&data[(idthetai*Nthetao*Nphid + idthetao2 * Nphid + idphid2) * 3]);
		Vector val111 = Vector(&data[(idthetai2*Nthetao*Nphid + idthetao2 * Nphid + idphid2) * 3]);
		Vector lerpVal = ((val000*(1. - fphid) + val001 * fphid)*(1. - fthetao) + (val010*(1. - fphid) + val011 * fphid)*fthetao) * (1. - fthetai)
			+ ((val100*(1. - fphid) + val101 * fphid)*(1. - fthetao) + (val110*(1. - fphid) + val111 * fphid)*fthetao) * fthetai;
		return lerpVal /*/ Vector(1.0, 1.15, 1.66)*/;
	}
	std::vector<float> data;
	int Nthetai;
	int Nthetao;
	int Nphid;
};

class IsoMERLBRDF : public BRDF { // todo
public:
	IsoMERLBRDF(std::string filename) {
		read_brdf(filename.c_str(), data);
	}
	void setParameters(const MaterialValues &m) {};
	Vector sample(const MaterialValues& mat, const Vector& wi, const Vector& N, double &pdf, double r1, double r2) const {
		Vector d = random_cos(N, r1, r2);
		pdf = dot(N, d) / (M_PI);
		return d;
	}
	Vector eval(const MaterialValues& mat, const Vector& wi, const Vector& wo, const Vector& N) const {

		/*Vector sumVal(0., 0., 0.);
		for (int i = 0; i < 90*90*360/2; i++) {
			sumVal += Vector(data[i] * RED_SCALE, data[i + BRDF_SAMPLING_RES_THETA_H * BRDF_SAMPLING_RES_THETA_D*BRDF_SAMPLING_RES_PHI_D / 2] * GREEN_SCALE, data[i + BRDF_SAMPLING_RES_THETA_H * BRDF_SAMPLING_RES_THETA_D*BRDF_SAMPLING_RES_PHI_D] * BLUE_SCALE);
		}
		std::cout << sumVal[0] << " " << sumVal[1] << " " << sumVal[2] << std::endl;*/

		Vector tangent1;
		Vector absN(abs(N[0]), abs(N[1]), abs(N[2]));
		if (absN[0] <= absN[1] && absN[0] <= absN[2]) {
			tangent1 = Vector(0, -N[2], N[1]);
		} else {
			if (absN[1] <= absN[0] && absN[1] <= absN[2]) {
				tangent1 = Vector(-N[2], 0, N[0]);
			} else
				tangent1 = Vector(-N[1], N[0], 0);
		}
		tangent1.normalize();
		Vector tangent2 = cross(tangent1, N);

		Vector wilocal = Vector(dot(wi, tangent1), dot(wi, tangent2), dot(wi, N));
		Vector wolocal = Vector(dot(wo, tangent1), dot(wo, tangent2), dot(wo, N));


		double thetai = acos(wilocal[2]);
		if (thetai >= M_PI / 2) return Vector(0., 0., 0.);
		double thetao = acos(wolocal[2]);
		if (thetao >= M_PI / 2) return Vector(0., 0., 0.);
		double phio = atan2(wolocal[1], wolocal[0]); // in -pi..pi 
		if (phio < 0) phio += 2 * M_PI;  // => 0..2pi
		double phii = atan2(wilocal[1], wilocal[0]); // in -pi..pi 
		if (phii < 0) phii += 2 * M_PI;  // => 0..2pi

		Vector result;
		lookup_brdf_val(data, thetai, phii, thetao, phio, result[0], result[1], result[2]);
		return result;

	}
	double* data;
};



class Texture {
public:
	Texture() { multiplier = Vector(1., 1., 1.); W = 0; H = 0; type = 0; };
	Texture(const char* filename, int texType, const Vector &multiplier) : type(texType) {
		W = 0;
		H = 0;
		switch (texType) {
		case 0: loadColors(filename); break;   // diffuse defaults to white
		case 1: loadColors(filename); break;  // specular defaults to black
		case 2: loadNormals(filename); break;   // normals defaults to 0,0,1
		case 3: loadColors(filename); break;   // alpha defaults to white
		case 4: loadColors(filename); break;   // roughness defaults to 1
		case 5: loadColors(filename); break;   // transparency defaults to 0
		case 6: loadColors(filename); break;   // refr index defaults to 1.3
		}
		this->multiplier = multiplier;
	}

	static double wrap(double u) {
		u = fmod(u, 1);
		if (u < 0) u += 1;
		return u;
	}
	static Texture defaultNormal() {
		return Texture("Null", true, Vector(0., 0., 1.));
	}
	static Texture defaultSpecular() {
		return Texture("Null", false, Vector(0., 0., 0.));
	}
	static Texture defaultDiffuse() {
		return Texture("Null", false, Vector(1., 1., 1.));
	}

	void clear_texture() {
		W = 0;
		H = 0;
		filename = "Null";
		values.clear();
	}

	Vector getVec(double u, double v) const {
		if (W > 0) {
			// no assert on u and v ; assume they are btw 0  and 1
			int x = u * (W - 1);
			int y = v * (H - 1);
			int idx = (y*W + x) * 3;
			double cr = values[idx] * multiplier[0];
			double cg = values[idx + 1] * multiplier[1];
			double cb = values[idx + 2] * multiplier[2];
			return Vector(cr, cg, cb);
		} else {
			return multiplier;
		}
		/*if (W > 0) {
			// no assert on u and v ; assume they are btw 0  and 1
			float x = u * (W - 1);
			float y = v * (H - 1);
			int ix = x;
			int iy = y;
			if (x > W - 2 || y > H - 2) {
				int idx = (iy*W + ix) * 3;
				return Vector(values[idx]* multiplier[0], values[idx + 1]* multiplier[1], values[idx + 2]* multiplier[2]);
			} else {
				float fracx = x - ix;
				float fracy = y - iy;
				int idx = (iy*W + ix) * 3;
				int idxX = idx + 3;
				int idxY = idx + 3 * W;
				int idxXY = idxY + 3;
				return multiplier*Vector((values[idx] * (1 - fracx) + values[idxX] * fracx)*(1 - fracy) + (values[idxY] * (1 - fracx) + values[idxXY] * fracx)*fracy,
					(values[idx + 1] * (1 - fracx) + values[idxX + 1] * fracx)*(1 - fracy) + (values[idxY + 1] * (1 - fracx) + values[idxXY + 1] * fracx)*fracy,
					(values[idx + 2] * (1 - fracx) + values[idxX + 2] * fracx)*(1 - fracy) + (values[idxY + 2] * (1 - fracx) + values[idxXY + 2] * fracx)*fracy);

			}
		} else {
			return multiplier;
		}*/
	}
	bool getBool(double u, double v) const {
		if (W > 0) {
			// no assert on u and v ; assume they are btw 0  and 1
			int x = u * (W - 1);
			int y = v * (H - 1);
			int idx = (y*W + x) * 3;
			double cr = values[idx] * multiplier[0];
			return cr < 0.5;
		} else {
			return (multiplier[0] < 0.5);
		}
	}
	Vector getNormal(double u, double v) const {
		if (W > 0) {
			// no assert on u and v ; assume they are btw 0  and 1
			int x = u * (W - 1);
			int y = v * (H - 1);
			int idx = (y*W + x) * 3;
			return Vector(values[idx], values[idx + 1], values[idx + 2]);
		} else {
			return Vector(0., 0., 1.);
		}
		/*if (W > 0) {
			// no assert on u and v ; assume they are btw 0  and 1
			float x = u * (W - 1);
			float y = v * (H - 1);
			int ix = x;
			int iy = y;
			if (x > W - 2 || y > H - 2) {
				int idx = (iy*W + ix) * 3;
				return Vector(values[idx], values[idx + 1], values[idx + 2]);
			} else {
				float fracx = x - ix;
				float fracy = y - iy;
				int idx = (iy*W + ix) * 3;
				int idxX = idx + 3;
				int idxY = idx + 3 * W;
				int idxXY = idxY + 3;
				return Vector((values[idx] * (1 - fracx) + values[idxX] * fracx)*(1 - fracy) + (values[idxY] * (1 - fracx) + values[idxXY] * fracx)*fracy,
					(values[idx+1] * (1 - fracx) + values[idxX+1] * fracx)*(1 - fracy) + (values[idxY+1] * (1 - fracx) + values[idxXY+1] * fracx)*fracy,
					(values[idx+2] * (1 - fracx) + values[idxX+2] * fracx)*(1 - fracy) + (values[idxY+2] * (1 - fracx) + values[idxXY+2] * fracx)*fracy);

			}
		} else {
			return Vector(0., 0., 1.);
		}*/
	}
	double getValRed(double u, double v) const {
		if (W > 0) {
			// no assert on u and v ; assume they are btw 0  and 1
			int x = u * (W - 1);
			int y = v * (H - 1);
			int idx = (y*W + x) * 3;
			double cr = values[idx] * multiplier[0];
			return cr;
		} else {
			return multiplier[0];
		}
	}
	void loadColors(const char* file) {
		if (file) {
			filename = std::string(file);
			if (load_image(file, values, W, H)) {
#pragma omp parallel for
				for (int i = 0; i < values.size(); i++) {
					values[i] /= 255.;
					values[i] = std::pow(values[i], 2.2); // values are then gamma corrected. 
				}
			}
		}
	}

	void loadNormals(const char* file) {
		if (file) {
			filename = std::string(file);
			if (load_image(file, values, W, H)) {
				for (int i = 0; i < values.size() / 3; i++) {
					Vector v(values[i * 3] - 128, values[i * 3 + 1] - 128, values[i * 3 + 2] - 128);
					v.normalize();
					values[i * 3] = v[0];
					values[i * 3 + 1] = v[1];
					values[i * 3 + 2] = v[2];
				}
			}
		}
	}

	size_t W, H;
	std::vector<double> values;
	Vector multiplier;
	std::string filename;
	unsigned char type; // 0: diffuse, 1: specular, 2: normal, 3:alpha, 4: roughness
};