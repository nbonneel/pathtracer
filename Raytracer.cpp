// Raytracer.cpp : Defines the entry point for the console application.
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <string>
#include "chrono.h"
#include "Raytracer.h"
#undef max
#undef min



inline double sqr(double x) { return x*x; }

void save_img(const char* filename, const unsigned char* pixels, int W, int H) // stored as RGB
{
	unsigned char bmpfileheader[14] = { 'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0 };
	unsigned char bmpinfoheader[40] = { 40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0 };
	unsigned char bmppad[3] = { 0,0,0 };

	int filesize = 54 + 3 * W*H;  //w is your image width, h is image height, both int
	bmpfileheader[2] = (unsigned char)(filesize);
	bmpfileheader[3] = (unsigned char)(filesize >> 8);
	bmpfileheader[4] = (unsigned char)(filesize >> 16);
	bmpfileheader[5] = (unsigned char)(filesize >> 24);

	bmpinfoheader[4] = (unsigned char)(W);
	bmpinfoheader[5] = (unsigned char)(W >> 8);
	bmpinfoheader[6] = (unsigned char)(W >> 16);
	bmpinfoheader[7] = (unsigned char)(W >> 24);
	bmpinfoheader[8] = (unsigned char)(H);
	bmpinfoheader[9] = (unsigned char)(H >> 8);
	bmpinfoheader[10] = (unsigned char)(H >> 16);
	bmpinfoheader[11] = (unsigned char)(H >> 24);

	FILE* f;
	f = fopen(filename, "wb");
	fwrite(bmpfileheader, 1, 14, f);
	fwrite(bmpinfoheader, 1, 40, f);
	std::vector<unsigned char> bgr_pixels(W*H * 3);
	for (int i = 0; i < W*H; i++) {
		bgr_pixels[i * 3] = pixels[i * 3 + 2];
		bgr_pixels[i * 3+1] = pixels[i * 3 + 1];
		bgr_pixels[i * 3+2] = pixels[i * 3];
	}
	for (int i = 0; i < H; i++)
	{
		fwrite(&bgr_pixels[0] + (W*(H - i - 1) * 3), 3, W, f);
		fwrite(bmppad, 1, (4 - (W * 3) % 4) % 4, f);
	}
	fclose(f);
}


Vector Phong_BRDF(const Vector& wi, const Vector& wo, const Vector& N, const Vector& phong_exponent) {

	Vector reflechi = wo.reflect(N);
	double d = dot(reflechi, wi);
	if (d < 0) return Vector(0., 0., 0.);
	Vector lobe = pow(Vector(d,d,d), phong_exponent) * ((phong_exponent + Vector(2.,2.,2.)) / (2.*M_PI));
	return lobe;
}


double int_exponential(double y0, double ysol, double beta, double s, double uy) {
	double 	result = 0.1 * exp((-y0 + ysol)*beta) * (1 - exp(-s*uy*beta)) / (uy*beta);
	return result;
}

Vector Raytracer::getColor(const Ray &r, const Scene &s, int nbrebonds, bool show_lights) {

	if (nbrebonds == 0) return Vector(0, 0, 0);

	Vector centerLight = s.lumiere->O + s.lumiere->get_translation(r.time);

	Vector P;
	MaterialValues mat;
	int sphere_id, tri_id;
	double t;
	bool has_inter = s.intersection(r, P, sphere_id, t, mat, tri_id);
	const Vector &N = mat.shadingN;

	/*Vector randV(uniform(engine), uniform(engine), uniform(engine));
	Vector T1 = cross(randV, N); T1.normalize();
	Vector T2 = cross(T1, N);*/

	Vector newOrigin = P+0.1*N;// +-log(uniform(engine))*0.15*T1 + -log(uniform(engine))*0.15*T2 + 0.1*N;   // BEWARE: not used for normal maps!

	Vector intensite_pixel(0, 0, 0);
	if (has_inter) {

		if (sphere_id == 0) {
			intensite_pixel = show_lights ? (/*s.lumiere->albedo **/Vector(1.,1.,1.)* (s.intensite_lumiere / sqr(s.lumiere->scale))) : Vector(0., 0., 0.);
		} else {

			intensite_pixel += mat.Ke*s.envmap_intensity;

			if (s.objects[sphere_id]->miroir) {
				Vector direction_miroir = r.direction.reflect(N);
				Ray rayon_miroir(P + 0.001*N, direction_miroir, r.time);

				//if (uniform(engine) < 0.9)
				intensite_pixel += getColor(rayon_miroir, s, nbrebonds - 1);// / 0.9;

			} else
				if (mat.transp) {
					double n1 = 1;
					double n2 = mat.refr_index;
					Vector normale_pour_transparence(N);
					Ray new_ray;
					bool entering = true;
					if (dot(r.direction, N) > 0) {  // on sort de la sphere
						n1 = mat.refr_index;
						n2 = 1;
						normale_pour_transparence = -N;
						entering = false;
					}
					double radical = 1 - sqr(n1 / n2)*(1 - sqr(dot(normale_pour_transparence, r.direction)));
					if (radical > 0) {
						Vector direction_refracte = (n1 / n2)*(r.direction - dot(r.direction, normale_pour_transparence)*normale_pour_transparence) - normale_pour_transparence * sqrt(radical);
						Ray rayon_refracte(P - 0.001*normale_pour_transparence, direction_refracte, r.time);

						double R0 = sqr((n1 - n2) / (n1 + n2));
						double R;
						if (entering) {
							R = R0 + (1 - R0)*std::pow(1 + dot(r.direction, N), 5);
						} else {
							R = R0 + (1 - R0)*std::pow(1 - dot(direction_refracte, N), 5);
						}

						if (uniform(engine) < R) {
							new_ray = Ray(P + 0.001*normale_pour_transparence, r.direction.reflect(N), r.time);
						} else {
							new_ray = Ray(P - 0.001*normale_pour_transparence, direction_refracte, r.time);
						}
					} else {
						new_ray = Ray(P + 0.001*normale_pour_transparence, r.direction.reflect(N), r.time);
						//return Vector(0, 0, 0);
					}
					intensite_pixel += getColor(new_ray, s, nbrebonds - 1);
				} else {

					//if (dot(N, r.direction) > 0) N = -N;

			// contribution de l'eclairage direct
				/*Ray ray_light(P + 0.01*N, (s.lumiere->O - P).getNormalized(), r.time);
				Vector P_light, N_light, color_light;
				int sphere_id_light;
				double t_light;
				bool has_inter_light = s.intersection(ray_light, P_light, N_light, sphere_id_light, t_light, color_light);
				double d_light2 = (s.lumiere->O - P).getNorm2();

				if (has_inter_light && t_light*t_light < d_light2) {
					intensite_pixel = Vector(0, 0, 0);
				}
				else {
					intensite_pixel = albedo / M_PI * s.intensite_lumiere * std::max(0., dot((s.lumiere->O - P).getNormalized(), N)) / d_light2;
				}*/
					Vector axeOP = (P - centerLight).getNormalized();
					Vector dir_aleatoire = random_cos(axeOP);
					Vector point_aleatoire = dir_aleatoire * s.lumiere->R*s.lumiere->scale + centerLight;
					Vector wi = (point_aleatoire - P).getNormalized();
					double d_light2 = (point_aleatoire - P).getNorm2();
					Vector Np = dir_aleatoire;

					Ray ray_light(P + 0.01*wi, wi, r.time);
					double t_light;
					bool has_inter_light = s.intersection_shadow(ray_light, t_light, sqrt(d_light2));

					if ((has_inter_light) || (dot(N, wi) < 0)) {
						intensite_pixel += Vector(0, 0, 0);
					} else {
						//intensite_pixel = (s.intensite_lumiere / (4.*M_PI*d_light2) * std::max(0., dot(N, wi)) * dot(Np, -wi) / dot(axeOP, dir_aleatoire))*(M_PI) * ((1. - s.objects[sphere_id]->ks)*albedo/ (M_PI) + Phong_BRDF(wi, r.direction, N, s.objects[sphere_id]->phong_exponent)*s.objects[sphere_id]->ks * albedo);

						//Vector BRDF = albedo / M_PI   * (1. - s.objects[sphere_id]->ks)   + s.objects[sphere_id]->ks*Phong_BRDF(wi, r.direction, N, s.objects[sphere_id]->phong_exponent)*albedo;
						Vector BRDF = mat.Kd / M_PI + Phong_BRDF(wi, r.direction, N, mat.Ne)*mat.Ks;
						double J = 1.*dot(Np, -wi) / d_light2;
						double proba = dot(axeOP, dir_aleatoire) / (M_PI * s.lumiere->R*s.lumiere->R*s.lumiere->scale*s.lumiere->scale);
						if (proba > 0)
							intensite_pixel += s.intensite_lumiere / sqr(s.lumiere->scale) * std::max(0., dot(N, wi)) * J * BRDF / proba;
					}

					// ajout de la contribution indirecte
					double avgNe = (mat.Ne[0] + mat.Ne[1] + mat.Ne[2]) / 3.;
					Vector direction_aleatoire;
					double p = 1 - (mat.Ks[0] + mat.Ks[1] + mat.Ks[2]) / 3.;
					bool sample_diffuse;
					Vector R = r.direction.reflect(N);
					if (uniform(engine) < p) {
						sample_diffuse = true;
						direction_aleatoire = random_cos(N);
					} else {
						sample_diffuse = false;
						direction_aleatoire = random_Phong(R, avgNe);
					}

					if (dot(direction_aleatoire, N) < 0) return intensite_pixel;
					if (dot(direction_aleatoire, R) < 0) return intensite_pixel;

					Ray	rayon_aleatoire(P + 0.01*direction_aleatoire, direction_aleatoire, r.time);

					double proba_phong = (avgNe + 1) / (2.*M_PI) * pow(dot(R, direction_aleatoire), avgNe);
					double proba_globale = p * dot(N, direction_aleatoire) / (M_PI)+(1. - p)*proba_phong;
					if (proba_globale <= 0) return intensite_pixel;
					if (sample_diffuse)
						intensite_pixel += getColor(rayon_aleatoire, s, nbrebonds - 1, false) * mat.Kd * (dot(N, direction_aleatoire) / (M_PI) / proba_globale);
					else
						intensite_pixel += getColor(rayon_aleatoire, s, nbrebonds - 1, false)  * (dot(N, direction_aleatoire) * Phong_BRDF(direction_aleatoire, r.direction, N, mat.Ne) / proba_globale) *mat.Ks;

					//intensite_pixel += getColor(rayon_aleatoire, s, nbrebonds - 1, false)  * dot(N, direction_aleatoire) * Phong_BRDF(direction_aleatoire, r.direction, N, s.objects[sphere_id]->phong_exponent)*s.objects[sphere_id]->ks * albedo / proba_globale;
				}
		}
	}
	if (s.fog_density==0)
		return intensite_pixel;

	Vector Lv(0., 0., 0.);
	
	double p_uniform = 0.5;
	bool is_uniform_fog = (s.fog_type==0);
	double beta = s.fog_density;// is_uniform_fog ? 0.04 : 0.1;
	bool uniform_sampling_ray = true;
	int phase = 0; // 0 : uniform, 1: Schlick, 2: Rayleigh

	double int_ext;
	if (is_uniform_fog) {
		int_ext = beta * t;
	} else {
		int_ext = int_exponential(r.origin[1], s.objects[2]->get_translation(1)[1], beta, t, r.direction[1]);
	}
	double T = exp(-int_ext);


	double proba_t, random_t;
	double clamped_t = std::min(1000., t);
	if (uniform_sampling_ray) {		
		random_t = uniform(engine)*clamped_t;
		proba_t = 1. / clamped_t;
	} else {
		double alpha = 5. / clamped_t;

		do {
			random_t = -log(uniform(engine)) / alpha;
		} while (random_t > clamped_t);

		double normalization = 1. / alpha * (1 - exp(-alpha*clamped_t));
		proba_t = exp(-alpha*random_t) / normalization;
	}

	double int_ext_partielle;
	if (is_uniform_fog) {
		int_ext_partielle = beta * random_t;
	} else {
		int_ext_partielle = int_exponential(r.origin[1], s.objects[2]->get_translation(1)[1], beta, random_t, r.direction[1]);
	}
		

	Vector random_P = r.origin + random_t*r.direction;

	Vector random_dir;
	double proba_dir;

	Vector point_aleatoire;
	Vector axeOP = (random_P - centerLight).getNormalized();
	bool is_uniform;
	if (uniform(engine) < p_uniform) {
		random_dir = random_uniform();				
		is_uniform = true;
	} else {		
		Vector dir_aleatoire = random_cos(axeOP);
		point_aleatoire = dir_aleatoire * s.lumiere->R*s.lumiere->scale + centerLight;
		random_dir = (point_aleatoire - random_P).getNormalized();
		is_uniform = false;
	}


	double phase_func;
	double k = 0.4;
	switch (phase) {
	case 0:
		phase_func = 0.3 / (4.*M_PI);
		break;
	case 1:		
		phase_func = (1 - k*k) / (4.*M_PI*(1 + k*dot(random_dir, -r.direction)));
		break;
	case 2:
		phase_func = 3 / (16 * M_PI)*(1 + sqr(dot(random_dir, r.direction)));
		break;
	}

	Ray L_ray(random_P, random_dir, r.time);
	MaterialValues interMat;
	Vector interP;
	int interid, intertri;
	double intert;
	bool interinter = s.intersection(L_ray, interP, interid, intert, interMat, intertri);
	const Vector& interN = interMat.shadingN;

	double V;
	if (is_uniform) {
		V = 1;
	} else {
		double d_light2 = (point_aleatoire - random_P).getNorm2();
		if (interinter && intert*intert < d_light2*0.99) {
			V = 0;
		} else {
			V = 1;
		}
	}

	if (V == 0) {
		Lv = Vector(0., 0., 0.);
	} else {
		double pdf_uniform = 1. / (4.*M_PI);
		double J = dot(interN, -random_dir) / (interP - random_P).getNorm2();
		double pdf_light = (interinter && interid==0) ? (dot((interP - centerLight).getNormalized(), axeOP) / (M_PI * sqr(s.lumiere->R*s.lumiere->scale)) / J) : 0.;
		proba_dir = p_uniform * pdf_uniform + (1 - p_uniform)*pdf_light;

		Vector L = getColor(L_ray, s, nbrebonds - 1);

		double ext;
		if (is_uniform_fog) {
			ext = beta;
		} else {
			ext = 0.1 * exp(-beta*(random_P[1] - s.objects[2]->get_translation(1)[1]));
		}
		Lv = L * phase_func * ext * exp(-int_ext_partielle) / (proba_t * proba_dir);
	}

	return intensite_pixel * T + Lv;
}

void Raytracer::save_scene(const char* filename) {
	
	FILE* f = fopen(filename, "w+");
	fprintf(f, "W,H: %u, %u\n", W, H);
	fprintf(f, "nrays: %u\n", nrays);
	fprintf(f, "Cam: (%lf, %lf, %lf), (%lf, %lf, %lf), (%lf, %lf, %lf)\n", cam.position[0], cam.position[1], cam.position[2], cam.direction[0], cam.direction[1], cam.direction[2], cam.up[0], cam.up[1], cam.up[2]);
	fprintf(f, "fov: %lf\n", cam.fov);
	fprintf(f, "focus: %lf\n", cam.focus_distance);
	fprintf(f, "aperture: %lf\n", cam.aperture);
	fprintf(f, "sigma_filter: %lf\n", sigma_filter);
	fprintf(f, "bounces: %u\n", nb_bounces);
	fprintf(f, "intensite_lum: %le\n", s.intensite_lumiere);
	fprintf(f, "intensite_envmap: %lf\n", s.envmap_intensity);


	fprintf(f, "nbobjects: %u\n", s.objects.size());
	for (int i = 0; i < s.objects.size(); i++) {
		s.objects[i]->save_to_file(f);
	}

	fprintf(f, "fog_density: %lf\n", s.fog_density);
	fprintf(f, "fog_type: %u\n", s.fog_type);	
	fclose(f);
}


void Raytracer::load_scene(const char* filename) {

	s.clear();

	FILE* f = fopen(filename, "r");
	fscanf(f, "W,H: %u, %u\n", &W, &H);
	fscanf(f, "nrays: %u\n", &nrays);
	fscanf(f, "Cam: (%lf, %lf, %lf), (%lf, %lf, %lf), (%lf, %lf, %lf)\n", &cam.position[0], &cam.position[1], &cam.position[2], &cam.direction[0], &cam.direction[1], &cam.direction[2], &cam.up[0], &cam.up[1], &cam.up[2]);
	fscanf(f, "fov: %lf\n", &cam.fov);
	fscanf(f, "focus: %lf\n", &cam.focus_distance);
	fscanf(f, "aperture: %lf\n", &cam.aperture);
	fscanf(f, "sigma_filter: %lf\n", &sigma_filter);
	fscanf(f, "bounces: %u\n", &nb_bounces);
	fscanf(f, "intensite_lum: %le\n", &s.intensite_lumiere);
	fscanf(f, "intensite_envmap: %lf\n", &s.envmap_intensity);

	int nbo;
	fscanf(f, "nbobjects: %u\n", &nbo);

	for (int i = 0; i < nbo; i++) {
		s.objects.push_back(Object::create_from_file(f));
	}

	fscanf(f, "fog_density: %lf\n", &s.fog_density);
	fscanf(f, "fog_type: %u\n", &s.fog_type);

	fclose(f);

	s.lumiere = dynamic_cast<Sphere*>(s.objects[0]);
}

void Raytracer::loadScene() {

	//W = 1000;
	//H = 1200;
	W = 600;
	H = 720;
	nrays = 1000;   //10 pour debugguer, 1000 pour les rendus finaux	
	cam = Camera(Vector(0., 0., 50.), Vector(0, 0, -1), Vector(0, 1, 0));
	cam.fov = 35 * M_PI / 180;   // field of view
	cam.focus_distance = 50; // distance mise au point	
	cam.aperture = 0.1;     // ouverture, pour depth of field (faible valeur = net partout)
	sigma_filter = 0.5;  // antialiasing
	nb_bounces = 3;

	Sphere* slum = new Sphere(Vector(10, 23, 15), 10);   //sphere light: position, rayon, albedo (qui sert à rien)
	Sphere* s2 = new Sphere(Vector(0, 0, 0), 1000000); // envmap  

	Plane* plane = new Plane(Vector(0, 0, 0), Vector(0., 1., 0.));
	plane->max_translation = Vector(0.,-27.3,0.);

	
	s.addObject(slum);   // sphere light
	s.addObject(s2); // envmap
	s.addObject(plane);
	

	s.lumiere = slum;
	s.intensite_lumiere = 1000000000 * 4.*M_PI / (4.*M_PI*s.lumiere->R*s.lumiere->R*M_PI);
	s.envmap_intensity = 1;

	cam.rotate(0, -22 * M_PI / 180, 1); // helmet
}

double sum_area_table(double* sat, int sat_width, int i0, int i1, int j0, int j1) {

	double term1 = 0;
	if (i0 > 0) {
		term1 = sat[(i0 - 1)*sat_width + j1];
	}
	double term2 = 0;
	if (j0 > 0) {
		term2 = sat[i1*sat_width + j0-1];
	}
	double term3 = 0;
	if (i0>0 && j0 > 0) {
		term3 = sat[(i0-1)*sat_width + j0 - 1];
	}
	return sat[i1*sat_width + j1] - term1 - term2 + term3;
}

//A Fast, Compact Approximation of the Exponential Function
double fast_exp(double y) {
 double d;
* ((int*)(&d) + 0) = 0;
* ((int*)(&d) + 1) = (int)(1512775 * y + 1072632447);
return d;
}
void Raytracer::render_image()
{

	int fstart = 0;
	int fend = 1;
	Wlr = ceil(W / 16.);
	Hlr = ceil(H / 16. );
	computed.resize(W*H, false);
	imagedouble.resize(W*H * 3, 0);  
	imagedouble_lowres.resize(Wlr*Hlr * 3, 0);
	sample_count.resize(W*H, 0);
	
	for (int i = 0; i < W*H; i++) {
		computed[i] = false;
		sample_count[i] = 0;
	}
	for (int i = 0; i < Wlr*Hlr; i++) {
		imagedouble_lowres[i * 3] = 0;
		imagedouble_lowres[i * 3 + 1] = 0;
		imagedouble_lowres[i * 3 + 2] = 0;
	}

	for (int i = 0; i < s.objects.size(); i++) {
		s.objects[i]->build_matrix(0); // time = 0 here
	}

	
	//static int nbcalls = 0;
	std::vector<std::pair<int, int> > permut(64);
	for (int i = 0; i < 8; i++)
		for (int j = 0; j < 8; j++)
			permut[i * 8 + j] = std::make_pair(i, j);
	std::random_shuffle(permut.begin(), permut.end());

	int filter_size = std::ceil(sigma_filter * 2);
	int filter_total_width = 2 * filter_size + 1;
	filter_integral.resize(filter_total_width*filter_total_width);	
	filter_value.resize(filter_total_width*filter_total_width);
	for (int i = -filter_size; i <= filter_size; i++) {
		for (int j = -filter_size; j <= filter_size; j++) {
			double integ = 0;
			for (int i2 = -filter_size; i2 <= i; i2++) {
				for (int j2 = -filter_size; j2 <= j; j2++) {
					double w = fast_exp(-(i2*i2+j2*j2) / (2.*sigma_filter*sigma_filter)) / (sigma_filter*sigma_filter*2.*M_PI);
					integ += w;
				}
			}
			filter_integral[(i + filter_size)*filter_total_width + (j + filter_size)] = integ;
			filter_value[(i + filter_size)*filter_total_width + (j + filter_size)] = exp(-(i*i + j*j) / (2.*sigma_filter*sigma_filter)) / (sigma_filter*sigma_filter*2.*M_PI);
		}
	}


	/*double test = sum_area_table(&filter_integral[0], filter_total_width, filter_size- filter_size, filter_size+15, filter_size+8, filter_size+11);
	double integ = 0;
	for (int i2 = -filter_size; i2 <= 15; i2++) {
		for (int j2 = 8; j2 <= 11; j2++) {
			double w = exp(-(i2*i2 + j2*j2) / (2.*sigma_filter*sigma_filter)) / (sigma_filter*sigma_filter*2.*M_PI);
			integ += w;
		}
	}*/



	for (int time_step = fstart; time_step < fend; time_step++) {
		for (int k = 0; k < nrays; k++) {
			current_nb_rays = k;
			chrono.Start();
			for(int pix_in_block=0; pix_in_block<64; pix_in_block++) {
				int i1 = permut[pix_in_block].first;
				int j1 = permut[pix_in_block].second;
			//for (int i1 = 0; i1 < 8; i1++) {
				//for (int j1 = 0; j1 < 8; j1++) {
					if (stopped) {						
						return;
					}
#pragma omp parallel for schedule(dynamic, 2)
					for (int i = i1; i < H; i += 8) {

						for (int j = j1; j < W; j += 8) {

							double dx = uniform(engine) - 0.5;
							double dy = uniform(engine) - 0.5;

							double dx_aperture = (uniform(engine) - 0.5) * cam.aperture;
							double dy_aperture = (uniform(engine) - 0.5) * cam.aperture;

							double time = time_step / 60. + uniform(engine) / 60.;



							Ray r = cam.generateDirection(i, j, time, dx, dy, dx_aperture, dy_aperture, W, H);

							Vector color = getColor(r, s, nb_bounces);

							int bmin_i = std::max(0, i - filter_size);
							int bmax_i = std::min(i + filter_size, H - 1);
							int bmin_j = std::max(0, j - filter_size);
							int bmax_j = std::min(j + filter_size, W - 1);
							double ratio = 1. / sum_area_table(&filter_integral[0], filter_total_width, bmin_i - i + filter_size, bmax_i - i + filter_size, bmin_j - j + filter_size, bmax_j - j + filter_size);
							double denom1 = ratio / (sigma_filter*sigma_filter*2.*M_PI);
							double denom2 = 1. / (2.*sigma_filter*sigma_filter);
							double sw = 0;
							for (int i2 = bmin_i; i2 <= bmax_i; i2++) {
								for (int j2 = bmin_j; j2 <= bmax_j; j2++) {
									//double w = filter_value[(i2 - i + filter_size)*filter_total_width + j2-j + filter_size] * ratio; 

									double w = fast_exp(-(sqr(i2-i-dy)+ sqr(j2 - j - dx)) *denom2) *denom1;
									imagedouble[((H - i2 - 1)*W + j2) * 3 + 0] += color[0] * w;
									imagedouble[((H - i2 - 1)*W + j2) * 3 + 1] += color[1] * w;
									imagedouble[((H - i2 - 1)*W + j2) * 3 + 2] += color[2] * w;																		
									//sw += w;
									sample_count[(H - i2 - 1)*W + j2] += w; // 1. / sqr(filter_total_width);
								}
							}
							//sample_count[(H - i - 1)*W + j] += sw;

							computed[(H - i - 1)*W + j] = true;

							imagedouble_lowres[((Hlr - (i)/16 - 1)*Wlr + (j)/16) * 3 + 0] += color[0]/256.;
							imagedouble_lowres[((Hlr - (i)/16 - 1)*Wlr + (j)/16) * 3 + 1] += color[1] / 256.;
							imagedouble_lowres[((Hlr - (i)/16 - 1)*Wlr + (j)/16) * 3 + 2] += color[2] / 256.;

							/*bmin_i = std::max(0, i/16 - filter_size);
							bmax_i = std::min(i/16 + filter_size, Hlr - 1);
							bmin_j = std::max(0, j/16 - filter_size);
							bmax_j = std::min(j/16 + filter_size, Wlr - 1);
							ratio = 1. / sum_area_table(&filter_integral[0], filter_total_width, bmin_i - i/16 + filter_size, bmax_i - i/16 + filter_size, bmin_j - j/16 + filter_size, bmax_j - j/16 + filter_size);
							double sw = 0;
							for (int i2 = bmin_i; i2 <= bmax_i; i2++) {
								for (int j2 = bmin_j; j2 <= bmax_j; j2++) {
									double w = filter_value[(i2 - i/16 + filter_size)*filter_total_width + j2 - j/16 + filter_size] * ratio/256.;
									imagedouble_lowres[((Hlr - i2 - 1)*Wlr + j2) * 3 + 0] += color[0] * w;
									imagedouble_lowres[((Hlr - i2 - 1)*Wlr + j2) * 3 + 1] += color[1] * w;
									imagedouble_lowres[((Hlr - i2 - 1)*Wlr + j2) * 3 + 2] += color[2] * w;
									sw += w;
								}
							}*/

						}
					}
				}
			//}

			curTimePerFrame = chrono.GetDiffMs();
			/*nbcalls++;
			if (nbcalls%12==0)
				std::cout << nbcalls / (double)H * 100. << std::endl;*/
		}


		image.resize(W*H * 3, 0);
#pragma omp parallel for 
		for (int i = 0; i < H; i++) {
			for (int j = 0; j < W; j++) {
				image[((H - i - 1)*W + j) * 3 + 0] = std::min(255., std::max(0., std::pow(imagedouble[((H - i - 1)*W + j) * 3 + 0] / (nrays + 1), 1 / 2.2)));   // rouge
				image[((H - i - 1)*W + j) * 3 + 1] = std::min(255., std::max(0., std::pow(imagedouble[((H - i - 1)*W + j) * 3 + 1] / (nrays + 1), 1 / 2.2))); // vert
				image[((H - i - 1)*W + j) * 3 + 2] = std::min(255., std::max(0., std::pow(imagedouble[((H - i - 1)*W + j) * 3 + 2] / (nrays + 1), 1 / 2.2))); // bleu
			}
		}

		std::ostringstream os;
		//os << "testFog_" << time_step << ".bmp";
		os << "export" << time_step << ".bmp";
		save_img(os.str().c_str(), &image[0], W, H);
	}

	//std::cout << chrono.GetDiffMs() / (double)1000. << std::endl;
    return;
}

