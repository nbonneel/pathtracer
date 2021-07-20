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
#include <list>

#undef max
#undef min




float int_exponential(float y0, float ysol, float beta, float s, float uy) {
	//double 	result =  exp((-y0 + ysol)*beta) * (1 - exp(-s*uy*beta)) / (uy*beta); // not ok
	//double 	result = (exp(-beta * (y0 - ysol)) - exp(-beta*(y0 - ysol +s*uy)))/(uy*beta);
	
	//double 	result = -1.0 / (uy*beta) * (exp(-beta * (y0 + s * uy - ysol)) - exp(-beta * (y0 - ysol))); //ok
	float 	result;
	if (abs(uy*beta) < 0.0001) { // I may have issues with smalle values
		result = exp(-beta * (y0 - ysol)) * (s /*- s * s*uy*beta*0.5*/);
		//double result2 = (exp(-beta * (y0 - ysol)) - exp(-beta * (y0 + s * uy - ysol))) / (uy*beta); // ok
		//std::cout << result << " " << result2 << std::endl;
	} else {
		result = (exp(-beta * (y0 - ysol)) - exp(-beta * (y0 + s * uy - ysol))) / (uy*beta); // ok
	}
	//double 	result2 = exp(-beta * (y0 - ysol)) * (1 - exp(-beta * (s * uy))) / (uy*beta); //not ok
	//double 	result = (-beta * (y0 - ysol)) + log(1 - exp(-beta * (s * uy))) - log(beta);
	//result = exp(result)/ (uy);
	//std::cout << result1 << "  " << result2 << "  " << result << std::endl;
	return result;
}





bool Raytracer::fogContribution(const Ray &r, const Vector& sampleLightPos, float t, Vector curWeight, int nbrebonds, bool showLight, bool hadSS, Contrib& newContrib, float &attenuationFactor) {
	if (curWeight.getNorm2() < 1E-12) return false;
	Vector rayDirection = r.direction; 

	Vector Lv(0., 0., 0.);

	float p_uniform = 0.5f;
	bool is_uniform_fog = (s.fog_type == 0);
	float alpha = s.fog_absorption;
	float sigmaT = s.fog_absorption_decay;
	bool uniform_sampling_ray = false;
	int phase = s.fog_phase_type; // 0 : uniform, 1: Schlick, 2: Rayleigh
	float groundLevel = s.objects[2]->get_translation(r.time, is_recording)[1];

	float int_ext;
	if (is_uniform_fog) {
		int_ext = alpha * t *0.05;
	} else {
		int_ext = alpha*int_exponential(r.origin[1], groundLevel, sigmaT, t, rayDirection[1]);
	}
	float T = exp(-int_ext);

	int threadid = omp_get_thread_num();
	float proba_t, random_t;
	float clamped_t = std::min(1000.f, t);
	
	float a = dot(sampleLightPos - r.origin, r.direction);
	if (a > 0) {
		// http://library.imageworks.com/pdfs/imageworks-library-importance-sampling-of-area-lights-in-participating-media.pdf
		Vector projP = r.origin + a * r.direction;
		float D = sqrt((sampleLightPos - projP).getNorm2());
		float thetaA = -atan2(a, D);
		float b = t - a;
		float thetaB = atan2(b, D); // can be negative. should be ok.
		Vector farP = r.origin + t * r.direction;		
		float x = engine[threadid]()*invmax;
		random_t = D * tan((1 - x)*thetaA + x * thetaB);
		proba_t = D / ((thetaB - thetaA)*(D*D + random_t * random_t));
		random_t += a;

	} else {
		if (uniform_sampling_ray) {
			random_t = engine[threadid]()*invmax*clamped_t;
			proba_t = 1. / clamped_t;
		} else {
			float alpha = 5.f / clamped_t;

			do {
				random_t = -log(engine[threadid]()*invmax) / alpha;
			} while (random_t > clamped_t);

			float normalization = 1.f / alpha * (1.f - exp(-alpha * clamped_t));
			proba_t = exp(-alpha * random_t) / normalization;
		}
	}
	

	float int_ext_partielle;
	if (is_uniform_fog) {
		int_ext_partielle = alpha * random_t *0.05;
	} else {
		int_ext_partielle = alpha * int_exponential(r.origin[1], groundLevel, sigmaT, random_t, rayDirection[1]);
	}


	Vector random_P = r.origin + random_t * rayDirection;
	if (random_P[1] < groundLevel) return false;

	Vector random_dir;
	float proba_dir;

	Vector point_aleatoire;
	Vector axeOP = (random_P - centerLight).getNormalized();
	bool is_uniform;
	if (engine[threadid]()*invmax < p_uniform) {
		random_dir = random_uniform_sphere<float>();
		is_uniform = true;
	} else {
		Vector dir_aleatoire = random_cos(axeOP);
		point_aleatoire = dir_aleatoire * radiusLight + centerLight;
		random_dir = (point_aleatoire - random_P).getNormalized();
		is_uniform = false;
	}


	float phase_func;
	float k = s.phase_aniso;
	switch (phase) {
	case 0:
		phase_func = 1. / (4.*M_PI);
		break;
	case 1:
		phase_func = (1 - k * k) / (4.*M_PI*(1 + k * dot(random_dir, -rayDirection)));
		break;
	case 2:
		phase_func = 3 / (16 * M_PI)*(1 + sqr(dot(random_dir, rayDirection)));
		break;
	}

	Ray L_ray(random_P, random_dir, r.time);
	MaterialValues interMat;
	Vector interP;
	int interid, intertri;
	float intert;
	bool interinter = s.intersection(L_ray, interP, interid, intert, interMat, intertri);
	const Vector& interN = interMat.shadingN;

	float V;
	if (is_uniform) {
		V = 1;
	} else {
		float d_light2 = (point_aleatoire - random_P).getNorm2();
		if (interinter && intert*intert < d_light2*0.99) {
			V = 0;
		} else {
			V = 1;
		}
	}

	attenuationFactor = T;

	if (V == 0) {
		Lv = Vector(0.f, 0.f, 0.f);
		return false;
	} else {
		float pdf_uniform = 1. / (4.*M_PI);
		float J = dot(interN, -random_dir) / (interP - random_P).getNorm2();
		float pdf_light = (interinter && interid == 0) ? (dot((interP - centerLight).getNormalized(), axeOP) / (M_PI * sqr(radiusLight)) / J) : 0.;
		proba_dir = p_uniform * pdf_uniform + (1 - p_uniform)*pdf_light;

		//Vector L = getColor(L_ray, sampleID, nbrebonds - 1, screenI, screenJ, normalValue, albedoValue);

		float ext;
		if (is_uniform_fog) {
			ext = s.fog_density*0.05;
		} else {
			//ext = 0.1 * exp(-beta * (random_P[1] - s.objects[2]->get_translation(r.time, is_recording)[1]));
			ext = s.fog_density * exp(-s.fog_density_decay * (random_P[1] - groundLevel));
		}
		//Lv = L * phase_func * ext * exp(-int_ext_partielle) / (proba_t * proba_dir);
		Vector newweight = curWeight * (phase_func * ext * exp(-int_ext_partielle) / (proba_t * proba_dir));
		if (isnan(newweight[0])) {
			std::cout << "nan" << std::endl;
		}
		newContrib = Contrib(newweight, L_ray, nbrebonds - 1, showLight, hadSS);
		return true;
	}
	//return intensite_pixel * T + Lv;
}



Vector Raytracer::getColor(const Ray &r, int sampleID, int nbrebondsss, int screenI, int screenJ, Vector &normalValue, Vector &albedoValue, bool no_envmap, bool has_precomputed_rays, int rayID) {


	int threadid = omp_get_thread_num();

	Vector P;
	MaterialValues mat;
	int sphere_id, tri_id;
	float t;
	Contrib newContrib;
	float attenuationFactor;
	bool has_fog = (s.fog_density > 1E-8);

	Vector color(0, 0, 0);
	//Vector pathWeight(1, 1, 1);

	//Ray currentRay = r;
	//std::list<Contrib> contribs;	
	Contrib* contribs = &contribsArray[threadid][0];
	int contribIndexStart = 0;
	int contribIndexEnd = 1;
	//contribs.push_front(Contrib(Vector(1,1,1), r, nb_bounces, true, false));
	contribs[contribIndexStart] = Contrib(Vector(1.f, 1.f, 1.f), r, nb_bounces, true, false);
	bool has_dome = sphereEnv;
	bool has_backgroundimage = (s.backgroundW > 0) && (s.background.size() == s.backgroundW*s.backgroundH * 3);


	/*for (int nbrebonds = nb_bounces; nbrebonds > 0; nbrebonds--) {*/ // finally, not just a linear "tree" due to multiple scattering.
	while (contribIndexStart!= contribIndexEnd)
	{
		const Contrib& curContrib = contribs[contribIndexStart]; //contribs.front();
		Ray currentRay = curContrib.r;
		int nbrebonds = curContrib.depth;	
		Vector pathWeight = curContrib.weight;
		

		bool has_had_subsurface_interaction = curContrib.has_had_subsurface_interaction;
		bool show_lights = curContrib.show_lights;
		//contribs.pop_front();
		contribIndexStart++;
		if (contribIndexStart >= sizeCircArray) contribIndexStart = 0;
		//if (contribIndexStart < 0) contribIndexStart = sizeCircArray-1;

		if (nbrebonds == 0) continue;
		if (pathWeight.getNorm2() < sqr(0.01f)) continue;

		bool has_inter;
		if (has_precomputed_rays && nbrebonds == nb_bounces) {
			has_inter = s.firstIntersection_has_inter[threadid][rayID];
			P = s.firstIntersection_P[threadid][rayID];
			sphere_id = s.firstIntersection_sphere_id[threadid][rayID];
			mat = s.firstIntersection_mat[threadid][rayID];
			tri_id = s.firstIntersection_triangle_id[threadid][rayID];
		} else {
			has_inter = s.intersection(currentRay, P, sphere_id, t, mat, tri_id, false, nbrebonds == nb_bounces);
		}
		Vector N = mat.shadingN;
		if (has_inter && nbrebonds == nb_bounces) {
			normalValue = N;
			albedoValue = mat.Kd;
		}
		

		if ((nbrebonds == nb_bounces) && has_backgroundimage && (!has_inter || (has_inter && sphere_id == 1 && has_dome))) {
			int i = std::min(s.backgroundH - 1, std::max(0, (int)(screenI / (float)H*s.backgroundH)));
			int j = std::min(s.backgroundW - 1, std::max(0, (int)(screenJ / (float)W*s.backgroundW)));
			float r = s.background[i*s.backgroundW * 3 + j * 3];
			float g = s.background[i*s.backgroundW * 3 + j * 3 + 1];
			float b = s.background[i*s.backgroundW * 3 + j * 3 + 2];
			color += pathWeight*Vector(r, g, b);
			continue;
		}

		Vector rayDirection = currentRay.direction;
		bool is_subsurface = (mat.Ksub.getNorm2() > 1E-8);


		if (has_inter) {
			if (sphere_id == 1) {
				if (no_envmap) {
					if (has_fog) {
						bool hasContrib = fogContribution(currentRay, centerLight, t, pathWeight, nbrebonds, show_lights, has_had_subsurface_interaction, newContrib, attenuationFactor);
						if (hasContrib) {
							//contribs.push_back(newContrib);
							contribs[contribIndexEnd] = newContrib;
							contribIndexEnd++; if (contribIndexEnd >= sizeCircArray) contribIndexEnd = 0;
						}
					}
					continue; // we do not break anymore: there can be secondary rays due to multiple scattering still being computed
				} else {
					if (sphereEnv) {
						if (has_fog) {
							bool hasContrib = fogContribution(currentRay, centerLight, t, pathWeight, nbrebonds, show_lights, has_had_subsurface_interaction, newContrib, attenuationFactor);
							if (hasContrib) {
								//contribs.push_back(newContrib);
								contribs[contribIndexEnd] = newContrib;
								contribIndexEnd++; if (contribIndexEnd >= sizeCircArray) contribIndexEnd = 0;
							}
							color += attenuationFactor * pathWeight * s.envmap_intensity*mat.Ke;
						} else {
							color += pathWeight * s.envmap_intensity*mat.Ke;					
						}
						continue;
					}
				}
			}
			if (sphere_id == 0) {
				Vector currentContrib = show_lights ? (/*s.lumiere->albedo **/Vector(lightPower, lightPower, lightPower)) : Vector(0.f, 0.f, 0.f);
				if (has_fog) {
					bool hasContrib = fogContribution(currentRay, centerLight, t, pathWeight, nbrebonds, show_lights, has_had_subsurface_interaction, newContrib, attenuationFactor);
					if (hasContrib) {
						//contribs.push_back(newContrib);
						contribs[contribIndexEnd] = newContrib;
						contribIndexEnd++; if (contribIndexEnd >= sizeCircArray) contribIndexEnd = 0;
					}
					color += attenuationFactor * pathWeight * currentContrib;
				} else {
					color += pathWeight * currentContrib;		
				}
				continue;
			} else {
				const float subsProba = (has_had_subsurface_interaction || !is_subsurface)?0.f:0.6f;// sqrt(mat.Ksub.getNorm2());
				
				const float inv1MSubsProba = 1.f / (1.f - subsProba);
				Vector subsW(inv1MSubsProba, inv1MSubsProba, inv1MSubsProba);

				bool sub_interaction = false;
				if (/*!has_had_subsurface_interaction &&*/ is_subsurface && (engine[threadid]()*invmax < subsProba)) {
					sub_interaction = true;
					const float invSubsProba = 1.f / subsProba;
					subsW = Vector(invSubsProba, invSubsProba, invSubsProba);

					//const double diskR = 6;
					const float sigmasub = 1.5f;
					const float diskR = sqrt(12.46f) * sigmasub; // sigmasub / sqrt(12.46);


					Vector gauss(0.f, 0.f, 1000.f);

					float integ = 1.f - exp(-diskR * diskR / (2.f * sigmasub*sigmasub));
					float randR = sigmasub * sqrt(-2.f * log(1.f - engine[threadid]()*invmax*integ));
					float randangle = engine[threadid]()*invmax * 2.f * (float)M_PI;
					gauss[0] = randR * sin(randangle);
					gauss[1] = randR * cos(randangle);
					gauss[2] = randR;
					float gaussval = (1. / (sigmasub*sigmasub * 2.f * (float)M_PI))*exp(-(gauss[2] * gauss[2]) / (2.f * sigmasub*sigmasub));
					float pdfgauss = gaussval / integ;

					Vector Tg = getTangent(N);
					Vector Tg2 = cross(N, Tg);
					Vector PtaboveP = P + gauss[0] * Tg + gauss[1] * Tg2 + N * diskR;
					float r1 = engine[threadid]()*invmax;
					Vector& axis = -N;
					float tmax;
					float h = sqrt(diskR*diskR - gauss[2] * gauss[2]);
					Vector subsOrigin = PtaboveP + (diskR - h)*(-N);
					float wAxis;
					if (r1 < 0.5f) {
						wAxis = 0.5f;
						//axis = -N;
						tmax = 2.f * h;
					} else {
						wAxis = 0.25f;
						tmax = 2.f * gauss[2];
						if (r1 < 0.75f) {
							axis = Tg;
						} else
							axis = Tg2;
						float r2 = engine[threadid]()*invmax;
						if (r2 < 0.5f) {
							subsOrigin -= h * N;
						}
					}

					MaterialValues subsmat;
					int subsid, substriid;
					float subst;
					Vector localP2;

					subsid = sphere_id;
					bool subsinter = s.get_random_intersection(Ray(subsOrigin, axis, r.time), localP2, subsid, subst, subsmat, substriid, 0, tmax, false);
					if (subsinter) {
						//subsW *= 1. / std::max(gaussval*std::abs(dot(mat.shadingN, subsmat.shadingN)),0.005) *exp(-(P - localP2).getNorm2() / (2.*sigmasub*sigmasub));
						float chris = exp(-(P - localP2).getNorm2() / (2.*sigmasub*sigmasub));
						//double dist = (P - localP2).getNorm2();
						//double chris = (exp(-sqrt(dist) / sigmasub) + exp(-sqrt(dist / (3*sigmasub)))) / (8*M_PI*dist*sigmasub);
						//double surfA = 1;// sqrt(mat.Ksub.getNorm2());
						//double sval = 1.9 - surfA + 3.5*sqr(surfA - 0.8);
						//double chris = 30*surfA * sval*  (exp(-sval * dist / sigmasub) + exp(-sval * dist / (3 * sigmasub))) / (8 * M_PI*sigmasub*dist);
						float sumpdfs = sqr(0.5*dot(subsmat.shadingN, N)) + sqr(0.25*dot(subsmat.shadingN, Tg)) + sqr(0.25*dot(subsmat.shadingN, Tg2));
						float pdfdisk = wAxis * std::abs(dot(axis, subsmat.shadingN)) / sumpdfs;
						//Vector A = sqrt(mat.Ksub * subsmat.Ksub);
						//Vector sval = Vector(1.9, 1.9, 1.9) - A + 3.5*sqr(A - Vector(0.8, 0.8, 0.8));
						//Vector chris = sval*A*(exp(-sval * dist / sigmasub) + exp(-sval * dist / (3 * sigmasub))) / (8 * M_PI*sigmasub*dist);
						subsW *= pdfdisk / std::max(pdfgauss, 0.05f) *chris;
						rayDirection = (localP2 - P).getNormalized();
						P = localP2 + 0.005f*subsmat.shadingN;
						if (r1 < 0.5f) {
							subsW *= 2.f;
						} else
							subsW *= 4.f;
						subsW = subsW * ((/*Vector(1.,1.,1.)-*/mat.Ksub/*-mat.Ks*/) / float(M_PI)); // Kd subsurface...
						mat = subsmat;
						N = mat.shadingN;
						tri_id = substriid;

					}


				}

				BRDF* brdf = s.objects[sphere_id]->brdf;
				//brdf->setParameters(mat);

				color += pathWeight* mat.Ke*s.envmap_intensity;

				if (s.objects[sphere_id]->miroir) {
					Vector direction_miroir = rayDirection.reflect(N);
					Ray rayon_miroir(P + 0.001f*N, direction_miroir, r.time);

					//if (uniform(engine) < 0.9)
					//currentContrib += getColor(rayon_miroir, sampleID, nbrebonds - 1, screenI, screenJ, normalValue, albedoValue, show_lights, no_envmap);// / 0.9;
					//currentRay = rayon_miroir;
					if (has_fog) {
						bool hasContrib = fogContribution(currentRay, centerLight, t, pathWeight, nbrebonds, show_lights, has_had_subsurface_interaction, newContrib, attenuationFactor);
						if (hasContrib) {
							//contribs.push_back(newContrib);
							contribs[contribIndexEnd] = newContrib;
							contribIndexEnd++; if (contribIndexEnd >= sizeCircArray) contribIndexEnd = 0;
						}
						contribs[contribIndexEnd] = Contrib(attenuationFactor*pathWeight, rayon_miroir, nbrebonds - 1, show_lights, has_had_subsurface_interaction);
						contribIndexEnd++; if (contribIndexEnd >= sizeCircArray) contribIndexEnd = 0;
						//contribs.push_back(Contrib(attenuationFactor*pathWeight, rayon_miroir, nbrebonds - 1, show_lights, has_had_subsurface_interaction));
					} else {
						//contribs.push_back(Contrib(pathWeight, rayon_miroir, nbrebonds - 1, show_lights, has_had_subsurface_interaction));
						contribs[contribIndexEnd] = Contrib(pathWeight, rayon_miroir, nbrebonds - 1, show_lights, has_had_subsurface_interaction);
						contribIndexEnd++; if (contribIndexEnd >= sizeCircArray) contribIndexEnd = 0;

					}
					continue;
				} else
					if (mat.transp) {
						float n1 = 1.f;
						float n2 = mat.refr_index;
						Vector normale_pour_transparence(N);
						Ray new_ray;
						bool entering = true;
						if (dot(rayDirection, N) > 0) {  // on sort de la sphere
							n1 = mat.refr_index;
							n2 = 1;
							normale_pour_transparence = -N;
							entering = false;
						}
						float radical = 1.f - sqr(n1 / n2)*(1.f - sqr(dot(normale_pour_transparence, rayDirection)));
						if (radical > 0) {
							Vector direction_refracte = (n1 / n2)*(rayDirection - dot(rayDirection, normale_pour_transparence)*normale_pour_transparence) - normale_pour_transparence * sqrt(radical);
							Ray rayon_refracte(P - 0.001f*normale_pour_transparence, direction_refracte, r.time);

							float R0 = sqr((n1 - n2) / (n1 + n2));
							float R;
							if (entering) {
								R = R0 + (1 - R0)*std::pow(1.f + dot(rayDirection, N), 5.f);
							} else {
								R = R0 + (1 - R0)*std::pow(1.f - dot(direction_refracte, N), 5.f);
							}

							if (engine[threadid]()*invmax < R) {
								new_ray = Ray(P + 0.001f*normale_pour_transparence, rayDirection.reflect(N), r.time);
							} else {
								new_ray = Ray(P - 0.001f*normale_pour_transparence, direction_refracte, r.time);
							}
						} else {
							new_ray = Ray(P + 0.001f*normale_pour_transparence, rayDirection.reflect(N), r.time);
							//return Vector(0, 0, 0);
						}
						//currentContrib += getColor(new_ray, sampleID, nbrebonds - 1, screenI, screenJ, normalValue, albedoValue, show_lights, no_envmap);
						//currentRay = new_ray;
						if (has_fog) {
							bool hasContrib = fogContribution(currentRay, centerLight, t, pathWeight, nbrebonds, show_lights, has_had_subsurface_interaction, newContrib, attenuationFactor);
							if (hasContrib) {
								//contribs.push_back(newContrib);
								contribs[contribIndexEnd] = newContrib;
								contribIndexEnd++; if (contribIndexEnd >= sizeCircArray) contribIndexEnd = 0;
							}
							contribs[contribIndexEnd] = Contrib(attenuationFactor*pathWeight, new_ray, nbrebonds - 1, show_lights, has_had_subsurface_interaction);
							contribIndexEnd++; if (contribIndexEnd >= sizeCircArray) contribIndexEnd = 0;
							//contribs.push_back(Contrib(attenuationFactor*pathWeight, new_ray, nbrebonds - 1, show_lights, has_had_subsurface_interaction));
						} else {
							contribs[contribIndexEnd] = Contrib(pathWeight, new_ray, nbrebonds - 1, show_lights, has_had_subsurface_interaction);
							contribIndexEnd++; if (contribIndexEnd >= sizeCircArray) contribIndexEnd = 0;
							//contribs.push_back(Contrib(pathWeight, new_ray, nbrebonds - 1, show_lights, has_had_subsurface_interaction));
						}
						continue;
					} else {

						//if (dot(N, rayDirection) > 0) N = -N;	

						Vector axeOP = (P - centerLight); axeOP.fast_normalize();//.getNormalized();
						Vector dir_aleatoire;
						if (nbrebonds == nb_bounces && no_envmap)
							dir_aleatoire = random_cos(axeOP, samples2d[sampleID][0], samples2d[sampleID][1]);
						else
							dir_aleatoire = random_cos(axeOP);
						Vector point_aleatoire = dir_aleatoire * radiusLight + centerLight;
						Vector wi = (point_aleatoire - P); wi.fast_normalize();// .getNormalized();
						float d_light2 = (point_aleatoire - P).getNorm2();
						Vector Np = dir_aleatoire;
						//Vector kd =  mat.Kd;
						Vector tr;
						if (s.objects[sphere_id]->ghost) {
							Vector offset;   // try to be robust here since we don't do nbrebonds-1
							if (dot(N, rayDirection) > 0)
								offset = N;
							else
								offset = -N; 
								//tr = getColor(Ray(P + rayDirection * 0.001 + offset * 0.001, rayDirection, r.time), sampleID, nbrebonds, screenI, screenJ, normalValue, albedoValue, show_lights, no_envmap);
							currentRay = Ray(P + rayDirection * 0.001f + offset * 0.001f, rayDirection, r.time);
							//contribs.push_back(Contrib(pathWeight, currentRay, nbrebonds, show_lights, has_had_subsurface_interaction));
							contribs[contribIndexEnd] = Contrib(pathWeight, currentRay, nbrebonds, show_lights, has_had_subsurface_interaction);
							contribIndexEnd++; if (contribIndexEnd >= sizeCircArray) contribIndexEnd = 0;
							continue;
						}
						Ray ray_light(P + 0.01f*wi, wi, r.time);
						float t_light;
						bool isShadowed;
						if (dot(mat.shadingN, wi) < 0)
							isShadowed = true;
						else
							isShadowed = s.intersection_shadow(ray_light, t_light, sqrt(d_light2) - 0.01f, true, nbrebonds == nb_bounces);
						//bool direct_visible = true;
						Vector currentContrib(0, 0, 0);

						if (!isShadowed) {
							//if (dot(N, wi) < 0) N = -N; // two sided shading...

							Vector BRDF;
							if (sub_interaction) {
								BRDF = (/*Vector(1., 1., 1.) -*/ mat.Ksub /*- mat.Ks*/) / (float)M_PI;
							} else {
								BRDF = brdf->eval(mat, wi, -rayDirection, N);
							}
							float J = dot(Np, -wi) / d_light2;
							float proba = dot(axeOP, dir_aleatoire) / (M_PI * radiusLight*radiusLight);
							if (s.objects[sphere_id]->ghost) {
								//currentContrib += tr;
							} else {
								if (proba > 0.f) {
									currentContrib += subsW * (lightPower* std::max(0.f, dot(N, wi)) * J / proba)  * BRDF;
								}
							}


						}
						if (has_fog) {
							bool hasContrib = fogContribution(currentRay, point_aleatoire, t, pathWeight, nbrebonds, show_lights, has_had_subsurface_interaction, newContrib, attenuationFactor);
							if (hasContrib) {
								//contribs.push_back(newContrib);
								contribs[contribIndexEnd] = newContrib;
								contribIndexEnd++; if (contribIndexEnd >= sizeCircArray) contribIndexEnd = 0;
							}
							color += attenuationFactor * pathWeight * currentContrib;
						} else {
							color += pathWeight * currentContrib;
						}

						// ajout de la contribution indirecte
						float proba_globale;
						Vector direction_aleatoire;
						if (nbrebonds == nb_bounces && no_envmap) {
							if (sub_interaction) {
								direction_aleatoire = random_cos(N);
								proba_globale = dot(N, direction_aleatoire);
							} else
								direction_aleatoire = brdf->sample(mat, -rayDirection, N, proba_globale);
						} else {
							float tmp;
							float r1 = modf(randomPerPixel[screenI*W + screenJ][0] + samples2d[sampleID][0], &tmp); // cranley patterson
							float r2 = modf(randomPerPixel[screenI*W + screenJ][1] + samples2d[sampleID][1], &tmp);
							if (sub_interaction) {
								direction_aleatoire = random_cos(mat.shadingN, r1, r2);
								proba_globale = dot(N, direction_aleatoire) / (float)M_PI;
							} else
								direction_aleatoire = brdf->sample(mat, -rayDirection, N, proba_globale, r1, r2);
						}


						if (dot(direction_aleatoire, N) < 0 || dot(direction_aleatoire, rayDirection.reflect(N)) < 0 || proba_globale <= 0) {
							//delete brdf;
							//return currentContrib;
							continue;
						}


						Ray	rayon_aleatoire(P + 0.01f*direction_aleatoire, direction_aleatoire, r.time);

						Vector BRDFindirect;
						if (sub_interaction) {
							BRDFindirect = mat.Ksub / (float)M_PI;
						} else {
							BRDFindirect = brdf->eval(mat, direction_aleatoire, -rayDirection, N);
						}

						//currentContrib += subsW * getColor(rayon_aleatoire, sampleID, nbrebonds - 1, screenI, screenJ, normalValue, albedoValue, false, no_envmap, sub_interaction ? true : has_had_subsurface_interaction)  * ((dot(N, direction_aleatoire) / proba_globale)) * BRDFindirect;

						Vector newpathWeight = pathWeight * subsW * BRDFindirect* ((dot(N, direction_aleatoire) / proba_globale));

						if (has_fog) {					
							//contribs.push_back(Contrib(attenuationFactor*newpathWeight, rayon_aleatoire, nbrebonds - 1, false, sub_interaction?true:has_had_subsurface_interaction));
							contribs[contribIndexEnd] = Contrib(attenuationFactor*newpathWeight, rayon_aleatoire, nbrebonds - 1, false, sub_interaction ? true : has_had_subsurface_interaction);
							contribIndexEnd++; if (contribIndexEnd >= sizeCircArray) contribIndexEnd = 0;
						} else {
							contribs[contribIndexEnd] = Contrib(newpathWeight, rayon_aleatoire, nbrebonds - 1, false, sub_interaction ? true : has_had_subsurface_interaction);
							contribIndexEnd++; if (contribIndexEnd >= sizeCircArray) contribIndexEnd = 0;
//							contribs.push_back(Contrib(newpathWeight, rayon_aleatoire, nbrebonds - 1, false, sub_interaction ? true : has_had_subsurface_interaction));
						}


						/*if (s.objects[sphere_id]->ghost) {
							//if (nbrebonds != nb_bounces)
							//	return intensite_pixel;
							if (sample_diffuse)
								intensite_pixel += tr / 196964.699 / M_PI *getColor(rayon_aleatoire, nbrebonds - 1, screenI, screenJ, false, direct_visible)  / proba_globale;
							else
								intensite_pixel += getColor(rayon_aleatoire, nbrebonds - 1, screenI, screenJ, false, true)  * Phong_BRDF(direction_aleatoire, rayDirection, N, mat.Ne) / proba_globale *mat.Ks;
						} else {
							if (sample_diffuse)
								intensite_pixel += getColor(rayon_aleatoire, nbrebonds - 1, screenI, screenJ, false, no_envmap) * ((dot(N, direction_aleatoire) / (M_PI) / proba_globale)) * kd;
							else
								intensite_pixel += getColor(rayon_aleatoire, nbrebonds - 1, screenI, screenJ, false, no_envmap)  * ((dot(N, direction_aleatoire) / proba_globale)) *mat.Ks * Phong_BRDF(direction_aleatoire, rayDirection, N, mat.Ne);
						}*/
						//intensite_pixel += getColor(rayon_aleatoire, nbrebonds - 1, false)  * dot(N, direction_aleatoire) * Phong_BRDF(direction_aleatoire, rayDirection, N, s.objects[sphere_id]->phong_exponent)*s.objects[sphere_id]->ks * albedo / proba_globale;
					}
				//delete brdf;
			}
			} 

			if (s.fog_density == 0)
				continue;

			if (!has_inter) break;


		}


	return color;
}

#if 0
Vector Raytracer::getColor(const Ray &r, int sampleID, int nbrebonds, int screenI, int screenJ, Vector &normalValue, Vector &albedoValue, bool show_lights, bool no_envmap, bool has_had_subsurface_interaction) {

	if (nbrebonds == 0) return Vector(0, 0, 0);

	Vector centerLight = s.lumiere->apply_transformation(s.lumiere->O);// s.lumiere->O + s.lumiere->get_translation(r.time, is_recording);
	double lum_scale = s.lumiere->get_scale(r.time, is_recording);
	double radiusLight = lum_scale *s.lumiere->R;	
	double lightPower = s.intensite_lumiere / sqr(lum_scale);

	Vector P;
	MaterialValues mat;
	int sphere_id, tri_id;
	double t;
	bool has_inter = s.intersection(r, P, sphere_id, t, mat, tri_id, false, nbrebonds==nb_bounces);
	Vector N = mat.shadingN;	
	if (has_inter && nbrebonds==nb_bounces) {
		normalValue = N;
		albedoValue = mat.Kd;
	}

	/*Vector randV(uniform(engine), uniform(engine), uniform(engine));
	Vector T1 = cross(randV, N); T1.normalize();
	Vector T2 = cross(T1, N);*/

//	Vector newOrigin = P+0.1*N;// +-log(uniform(engine))*0.15*T1 + -log(uniform(engine))*0.15*T2 + 0.1*N;   // BEWARE: not used for normal maps!

	// we just replace the envmap by the background image for direct rays
	if ( (nbrebonds == nb_bounces) && (s.backgroundW > 0) && (s.background.size() == s.backgroundW*s.backgroundH * 3) && (!has_inter || (has_inter && sphere_id==1 && dynamic_cast<Sphere*>(s.objects[1]) ))) {
		int i = std::min(s.backgroundH - 1, std::max(0, (int)(screenI / (double)H*s.backgroundH)));
		int j = std::min(s.backgroundW - 1, std::max(0, (int)(screenJ / (double)W*s.backgroundW)));
		double r = s.background[i*s.backgroundW * 3 + j * 3];
		double g = s.background[i*s.backgroundW * 3 + j * 3 + 1];
		double b = s.background[i*s.backgroundW * 3 + j * 3 + 2];
		return Vector(r, g, b);
	}

	Vector rayDirection = r.direction;
	bool is_subsurface = (mat.Ksub.getNorm2() > 1E-8);

	Vector intensite_pixel(0, 0, 0);
	if (has_inter) {
		int threadid = omp_get_thread_num();

		if (sphere_id == 1) {
			if (no_envmap)
				return Vector(0, 0, 0);
			else {
				Sphere* env = dynamic_cast<Sphere*>(s.objects[1]);
				if (env) {
					return s.envmap_intensity*mat.Ke;
				}
			}				
		}
		if (sphere_id == 0) {
			intensite_pixel = show_lights ? (/*s.lumiere->albedo **/Vector(1.,1.,1.)* lightPower) : Vector(0., 0., 0.);
		} else {
			
			double subsProba = 0.6;// sqrt(mat.Ksub.getNorm2());
			if (has_had_subsurface_interaction || !is_subsurface) subsProba = 0;
			Vector subsW = Vector(1. / (1. - subsProba), 1. / (1. - subsProba), 1. / (1. - subsProba));

			bool sub_interaction = false;
			if (/*!has_had_subsurface_interaction &&*/ is_subsurface && (engine[threadid]()*invmax < subsProba) ) {
				sub_interaction = true;
				subsW = Vector(1. / (subsProba), 1. / (subsProba), 1. / (subsProba));

				//const double diskR = 6;
				const double sigmasub = 1.5;
				const double diskR = sqrt(12.46) * sigmasub; // sigmasub / sqrt(12.46);
				/*float r1 = engine[threadid]()*invmax;
				float sr2 = sqrt(engine[threadid]()*invmax);
				float x = diskR*sin(2 * M_PI*r1)*sr2;
				float y = diskR*cos(2 * M_PI*r1)*sr2;*/

				Vector gauss(0, 0, 1000);
				/*while (gauss[2] > diskR) {
					gauss = boxMuller()*sigmasub;
				}*/
				double integ = 1 - exp(-diskR * diskR / (2 * sigmasub*sigmasub));
				double randR = sigmasub * sqrt(-2 * log(1 - engine[threadid]()*invmax*integ));

				double randangle = engine[threadid]()*invmax * 2 * M_PI;
				gauss[0] = randR * sin(randangle);
				gauss[1] = randR * cos(randangle);
				gauss[2] = randR;
				double gaussval = (1. / (sigmasub*sigmasub * 2 * M_PI))*exp(-(gauss[2] * gauss[2]) / (2 * sigmasub*sigmasub));
				double pdfgauss = gaussval / integ;

				Vector Tg = getTangent(N);
				Vector Tg2 = cross(N, Tg);
				Vector PtaboveP = P + gauss[0] * Tg + gauss[1] * Tg2 + N * diskR;
				float r1 = engine[threadid]()*invmax;
				Vector axis;
				double tmax;
				double h = sqrt(diskR*diskR - gauss[2] * gauss[2]);
				Vector subsOrigin = PtaboveP + (diskR - h)*(-N);
				double wAxis;
				if (r1 < 0.5) {
					wAxis = 0.5;
					axis = -N;
					tmax = 2 * h;
				} else {
					wAxis = 0.25;
					tmax = 2 * gauss[2];
					if (r1 < 0.75) {
						axis = Tg;
					} else
						axis = Tg2;
					float r2 = engine[threadid]()*invmax;
					if (r2 < 0.5) {
						subsOrigin -= h * N;
					}
				}
			
			MaterialValues subsmat;
				int subsid, substriid;
				double subst;
				Vector localP2;

				subsid = sphere_id;
				bool subsinter = s.get_random_intersection(Ray(subsOrigin, axis, r.time), localP2, subsid, subst, subsmat, substriid, 0, tmax, false);
				if (subsinter) {					
					//subsW *= 1. / std::max(gaussval*std::abs(dot(mat.shadingN, subsmat.shadingN)),0.005) *exp(-(P - localP2).getNorm2() / (2.*sigmasub*sigmasub));
					double chris = exp(-(P - localP2).getNorm2() / (2.*sigmasub*sigmasub));
					//double dist = (P - localP2).getNorm2();
					//double chris = (exp(-sqrt(dist) / sigmasub) + exp(-sqrt(dist / (3*sigmasub)))) / (8*M_PI*dist*sigmasub);
					//double surfA = 1;// sqrt(mat.Ksub.getNorm2());
					//double sval = 1.9 - surfA + 3.5*sqr(surfA - 0.8);
					//double chris = 30*surfA * sval*  (exp(-sval * dist / sigmasub) + exp(-sval * dist / (3 * sigmasub))) / (8 * M_PI*sigmasub*dist);
					double sumpdfs = sqr(0.5*dot(subsmat.shadingN, N)) + sqr(0.25*dot(subsmat.shadingN, Tg)) + sqr(0.25*dot(subsmat.shadingN, Tg2));
					double pdfdisk = wAxis * std::abs(dot(axis, subsmat.shadingN)) / sumpdfs;
					//Vector A = sqrt(mat.Ksub * subsmat.Ksub);
					//Vector sval = Vector(1.9, 1.9, 1.9) - A + 3.5*sqr(A - Vector(0.8, 0.8, 0.8));
					//Vector chris = sval*A*(exp(-sval * dist / sigmasub) + exp(-sval * dist / (3 * sigmasub))) / (8 * M_PI*sigmasub*dist);
					subsW *= pdfdisk / std::max(pdfgauss, 0.05) *chris;
					rayDirection = (localP2 - P).getNormalized();
					P = localP2 + 0.005*subsmat.shadingN;
					if (r1 < 0.5) {
						subsW *= 2;
					} else
						subsW *= 4;
					subsW = subsW*((/*Vector(1.,1.,1.)-*/mat.Ksub/*-mat.Ks*/) / M_PI); // Kd subsurface...
					mat = subsmat;
					N = mat.shadingN;
					tri_id = substriid;					

				}


			}

			BRDF* brdf = s.objects[sphere_id]->brdf; 
			//brdf->setParameters(mat);

			intensite_pixel += mat.Ke*s.envmap_intensity;

			if (s.objects[sphere_id]->miroir) {
				Vector direction_miroir = rayDirection.reflect(N);
				Ray rayon_miroir(P + 0.001*N, direction_miroir, r.time);

				//if (uniform(engine) < 0.9)
				intensite_pixel += getColor(rayon_miroir, sampleID, nbrebonds - 1, screenI, screenJ, normalValue, albedoValue, show_lights, no_envmap);// / 0.9;

			} else
				if (mat.transp) {
					double n1 = 1;
					double n2 = mat.refr_index;
					Vector normale_pour_transparence(N);
					Ray new_ray;
					bool entering = true;
					if (dot(rayDirection, N) > 0) {  // on sort de la sphere
						n1 = mat.refr_index;
						n2 = 1;
						normale_pour_transparence = -N;
						entering = false;
					}
					double radical = 1 - sqr(n1 / n2)*(1 - sqr(dot(normale_pour_transparence, rayDirection)));
					if (radical > 0) {
						Vector direction_refracte = (n1 / n2)*(rayDirection - dot(rayDirection, normale_pour_transparence)*normale_pour_transparence) - normale_pour_transparence * sqrt(radical);
						Ray rayon_refracte(P - 0.001*normale_pour_transparence, direction_refracte, r.time);

						double R0 = sqr((n1 - n2) / (n1 + n2));
						double R;
						if (entering) {
							R = R0 + (1 - R0)*std::pow(1 + dot(rayDirection, N), 5);
						} else {
							R = R0 + (1 - R0)*std::pow(1 - dot(direction_refracte, N), 5);
						}

						if (engine[threadid]()*invmax < R) {
							new_ray = Ray(P + 0.001*normale_pour_transparence, rayDirection.reflect(N), r.time);
						} else {
							new_ray = Ray(P - 0.001*normale_pour_transparence, direction_refracte, r.time);
						}
					} else {
						new_ray = Ray(P + 0.001*normale_pour_transparence, rayDirection.reflect(N), r.time);
						//return Vector(0, 0, 0);
					}
					intensite_pixel += getColor(new_ray, sampleID, nbrebonds - 1, screenI, screenJ, normalValue, albedoValue, show_lights, no_envmap);
				} else {

					//if (dot(N, rayDirection) > 0) N = -N;	
	
					Vector axeOP = (P - centerLight); axeOP.fast_normalize();//.getNormalized();
					Vector dir_aleatoire;
					if (nbrebonds == nb_bounces && no_envmap)
						dir_aleatoire = random_cos(axeOP, samples2d[sampleID][0], samples2d[sampleID][1]);
					else
						dir_aleatoire = random_cos(axeOP);
					Vector point_aleatoire = dir_aleatoire * radiusLight + centerLight;
					Vector wi = (point_aleatoire - P); wi.fast_normalize();// .getNormalized();
					double d_light2 = (point_aleatoire - P).getNorm2();
					Vector Np = dir_aleatoire;
					//Vector kd =  mat.Kd;
					Vector tr;
					if (s.objects[sphere_id]->ghost) {
						Vector offset;   // try to be robust here since we don't do nbrebonds-1
						if (dot(N, rayDirection) > 0)
							offset = N;
						else
							offset = -N;
						tr = getColor(Ray(P + rayDirection *0.001 + offset*0.001, rayDirection, r.time), sampleID, nbrebonds, screenI, screenJ, normalValue, albedoValue, show_lights, no_envmap);
					}
					Ray ray_light(P + 0.01*wi, wi, r.time);
					double t_light;
					bool has_inter_light = s.intersection_shadow(ray_light, t_light, sqrt(d_light2)-0.01, true);
					bool direct_visible = true;
					if ((has_inter_light) ) {
						//intensite_pixel += Vector(0, 0, 0);
						direct_visible = false;
					} else {
						
						//if (dot(N, wi) < 0) N = -N; // two sided shading...
						
						Vector BRDF;
						if (sub_interaction) {
							BRDF = (/*Vector(1., 1., 1.) -*/ mat.Ksub /*- mat.Ks*/) / M_PI;
						} else {
							BRDF = brdf->eval(mat, wi, -rayDirection, N);
						}
						double J = 1.*dot(Np, -wi) / d_light2;
						double proba = dot(axeOP, dir_aleatoire) / (M_PI * radiusLight*radiusLight);
						if (s.objects[sphere_id]->ghost) {
							intensite_pixel += tr;
						} else {
							if (proba > 0) {
								intensite_pixel += subsW*(lightPower* std::max(0., dot(N, wi)) * J / proba)  * BRDF;
							}
						}
						
							
					}

					// ajout de la contribution indirecte
					double proba_globale;
					Vector direction_aleatoire;
					if (nbrebonds == nb_bounces && no_envmap) {
						if (sub_interaction) {
							direction_aleatoire = random_cos(N);
							proba_globale = dot(N, direction_aleatoire);
						} else
						direction_aleatoire = brdf->sample(mat, -rayDirection, N, proba_globale);
					}  else {		
						float tmp;
						float r1 = modf(randomPerPixel[screenI*W+screenJ][0] + samples2d[sampleID][0], &tmp); // cranley patterson
						float r2 = modf(randomPerPixel[screenI*W + screenJ][1] + samples2d[sampleID][1], &tmp);
						if (sub_interaction) {
							direction_aleatoire = random_cos(mat.shadingN, r1, r2);
							proba_globale = dot(N, direction_aleatoire) / M_PI;
						} else
							direction_aleatoire = brdf->sample(mat, -rayDirection, N, proba_globale, r1, r2);
					}
						

					if (dot(direction_aleatoire, N) < 0 || dot(direction_aleatoire, rayDirection.reflect(N)) < 0 || proba_globale <= 0) {
						//delete brdf;
						return intensite_pixel;
					}


					Ray	rayon_aleatoire(P + 0.01*direction_aleatoire, direction_aleatoire, r.time);

					Vector BRDFindirect;
					if (sub_interaction) {
						BRDFindirect = mat.Ksub / M_PI;
					} else {
						BRDFindirect = brdf->eval(mat, direction_aleatoire, -rayDirection, N);
					}

					intensite_pixel += subsW*getColor(rayon_aleatoire, sampleID, nbrebonds - 1, screenI, screenJ, normalValue, albedoValue, false, no_envmap, sub_interaction?true:has_had_subsurface_interaction)  * ((dot(N, direction_aleatoire) / proba_globale)) * BRDFindirect;
					
					/*if (s.objects[sphere_id]->ghost) {
						//if (nbrebonds != nb_bounces)
						//	return intensite_pixel;
						if (sample_diffuse)
							intensite_pixel += tr / 196964.699 / M_PI *getColor(rayon_aleatoire, nbrebonds - 1, screenI, screenJ, false, direct_visible)  / proba_globale;
						else
							intensite_pixel += getColor(rayon_aleatoire, nbrebonds - 1, screenI, screenJ, false, true)  * Phong_BRDF(direction_aleatoire, rayDirection, N, mat.Ne) / proba_globale *mat.Ks;
					} else {
						if (sample_diffuse)
							intensite_pixel += getColor(rayon_aleatoire, nbrebonds - 1, screenI, screenJ, false, no_envmap) * ((dot(N, direction_aleatoire) / (M_PI) / proba_globale)) * kd;
						else
							intensite_pixel += getColor(rayon_aleatoire, nbrebonds - 1, screenI, screenJ, false, no_envmap)  * ((dot(N, direction_aleatoire) / proba_globale)) *mat.Ks * Phong_BRDF(direction_aleatoire, rayDirection, N, mat.Ne);
					}*/
					//intensite_pixel += getColor(rayon_aleatoire, nbrebonds - 1, false)  * dot(N, direction_aleatoire) * Phong_BRDF(direction_aleatoire, rayDirection, N, s.objects[sphere_id]->phong_exponent)*s.objects[sphere_id]->ks * albedo / proba_globale;
				}
				//delete brdf;
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
		int_ext = int_exponential(r.origin[1], s.objects[2]->get_translation(r.time, is_recording)[1], beta, t, rayDirection[1]);
	}
	double T = exp(-int_ext);

	int threadid = omp_get_thread_num();
	double proba_t, random_t;
	double clamped_t = std::min(1000., t);
	if (uniform_sampling_ray) {		
		random_t = engine[threadid]()*invmax*clamped_t;
		proba_t = 1. / clamped_t;
	} else {
		double alpha = 5. / clamped_t;

		do {
			random_t = -log(engine[threadid]()*invmax) / alpha;
		} while (random_t > clamped_t);

		double normalization = 1. / alpha * (1 - exp(-alpha*clamped_t));
		proba_t = exp(-alpha*random_t) / normalization;
	}

	double int_ext_partielle;
	if (is_uniform_fog) {
		int_ext_partielle = beta * random_t;
	} else {
		int_ext_partielle = int_exponential(r.origin[1], s.objects[2]->get_translation(r.time, is_recording)[1], beta, random_t, rayDirection[1]);
	}
		

	Vector random_P = r.origin + random_t* rayDirection;

	Vector random_dir;
	double proba_dir;

	Vector point_aleatoire;
	Vector axeOP = (random_P - centerLight).getNormalized();
	bool is_uniform;
	if (engine[threadid]()*invmax < p_uniform) {
		random_dir = random_uniform_sphere();				
		is_uniform = true;
	} else {		
		Vector dir_aleatoire = random_cos(axeOP);
		point_aleatoire = dir_aleatoire * radiusLight + centerLight;
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
		phase_func = (1 - k*k) / (4.*M_PI*(1 + k*dot(random_dir, -rayDirection)));
		break;
	case 2:
		phase_func = 3 / (16 * M_PI)*(1 + sqr(dot(random_dir, rayDirection)));
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
		double pdf_light = (interinter && interid==0) ? (dot((interP - centerLight).getNormalized(), axeOP) / (M_PI * sqr(radiusLight)) / J) : 0.;
		proba_dir = p_uniform * pdf_uniform + (1 - p_uniform)*pdf_light;

		Vector L = getColor(L_ray, sampleID, nbrebonds - 1, screenI, screenJ, normalValue, albedoValue);

		double ext;
		if (is_uniform_fog) {
			ext = beta;
		} else {
			ext = 0.1 * exp(-beta*(random_P[1] - s.objects[2]->get_translation(r.time, is_recording)[1]));
		}
		Lv = L * phase_func * ext * exp(-int_ext_partielle) / (proba_t * proba_dir);
	}

	return intensite_pixel * T + Lv;
}
#endif

void Raytracer::save_scene(const char* filename) {
	
	FILE* f = fopen(filename, "w+");
	fprintf(f, "W,H: %u, %u\n", W, H);
	fprintf(f, "nrays: %u\n", nrays);
	fprintf(f, "nbframes: %u\n", s.nbframes);
	fprintf(f, "Cam: (%f, %f, %f), (%f, %f, %f), (%f, %f, %f)\n", cam.position[0], cam.position[1], cam.position[2], cam.direction[0], cam.direction[1], cam.direction[2], cam.up[0], cam.up[1], cam.up[2]);
	fprintf(f, "fov: %f\n", cam.fov);
	fprintf(f, "focus: %f\n", cam.focus_distance);
	fprintf(f, "aperture: %f\n", cam.aperture);			
	fprintf(f, "sigma_filter: %f\n", sigma_filter);
	fprintf(f, "gamma: %f\n", gamma);

	fprintf(f, "is_lenticular: %u\n", cam.is_lenticular);
	fprintf(f, "lenticular_nb_images: %u\n", cam.lenticular_nb_images);
	fprintf(f, "lenticular_max_angle: %f\n", cam.lenticular_max_angle);
	fprintf(f, "lenticular_pixel_width: %u\n", cam.lenticular_pixel_width);
	fprintf(f, "isArray: %u\n", cam.isArray);
	fprintf(f, "nbviewX: %u\n", cam.nbviewX);
	fprintf(f, "nbviewY: %u\n", cam.nbviewY);
	fprintf(f, "maxSpacingX: %f\n", cam.maxSpacingX);
	fprintf(f, "maxSpacingY: %f\n", cam.maxSpacingY);

	fprintf(f, "bounces: %u\n", nb_bounces);
	fprintf(f, "has_denoiser: %u\n", has_denoiser);
	

	fprintf(f, "intensite_lum: %f\n", s.intensite_lumiere);
	fprintf(f, "intensite_envmap: %f\n", s.envmap_intensity);
	if (s.backgroundfilename.size()>0)
		fprintf(f, "background: %s\n", s.backgroundfilename.c_str());
	

	fprintf(f, "nbobjects: %u\n", static_cast<unsigned int>(s.objects.size()));
	for (int i = 0; i < s.objects.size(); i++) {
		s.objects[i]->save_to_file(f);
	}

	fprintf(f, "fog_density: %f\n", s.fog_density);
	fprintf(f, "fog_absorption: %f\n", s.fog_absorption);
	fprintf(f, "fog_density_decay: %f\n", s.fog_density_decay);
	fprintf(f, "fog_absorption_decay: %f\n", s.fog_absorption_decay);
	fprintf(f, "fog_type: %u\n", s.fog_type);	
	fprintf(f, "fog_phase_type: %u\n", s.fog_phase_type);

	fprintf(f, "double_frustum_start_t: %f\n", s.double_frustum_start_t);

	
	
	fclose(f);
}


void Raytracer::load_scene(const char* filename, const char* replacedNames) {

	s.clear();
	char line[255], bg[255];

	FILE* f = fopen(filename, "r");
	fscanf(f, "W,H: %u, %u\n", &W, &H);
	fscanf(f, "nrays: %u\n", &nrays);
	fscanf(f, "%[^\n]\n", line);
	if (line[0] == 'n') {
		sscanf(line, "nbframes: %u\n", &s.nbframes);
		fscanf(f, "Cam: (%f, %f, %f), (%f, %f, %f), (%f, %f, %f)\n", &cam.position[0], &cam.position[1], &cam.position[2], &cam.direction[0], &cam.direction[1], &cam.direction[2], &cam.up[0], &cam.up[1], &cam.up[2]);
	} else {
		sscanf(line, "Cam: (%f, %f, %f), (%f, %f, %f), (%f, %f, %f)\n", &cam.position[0], &cam.position[1], &cam.position[2], &cam.direction[0], &cam.direction[1], &cam.direction[2], &cam.up[0], &cam.up[1], &cam.up[2]);
	}
	fscanf(f, "fov: %f\n", &cam.fov);
	fscanf(f, "focus: %f\n", &cam.focus_distance);
	fscanf(f, "aperture: %f\n", &cam.aperture);
	fscanf(f, "sigma_filter: %f\n", &sigma_filter);
	fscanf(f, "gamma: %f\n", &gamma);

	int nbo;
	
	fscanf(f, "%[^\n]\n", line);
		
	int read = sscanf(line, "is_lenticular: %u\n", &nbo);
	if (read == 1) {
		cam.is_lenticular = nbo;
		fscanf(f, "lenticular_nb_images: %u\n", &cam.lenticular_nb_images);
		fscanf(f, "lenticular_max_angle: %f\n", &cam.lenticular_max_angle);
		fscanf(f, "lenticular_pixel_width: %u\n", &cam.lenticular_pixel_width);
		fscanf(f, "isArray: %u\n", &cam.isArray);
		fscanf(f, "nbviewX: %u\n", &cam.nbviewX);
		fscanf(f, "nbviewY: %u\n", &cam.nbviewY);
		fscanf(f, "maxSpacingX: %f\n", &cam.maxSpacingX);
		fscanf(f, "maxSpacingY: %f\n", &cam.maxSpacingY);
		fscanf(f, "bounces: %u\n", &nb_bounces);
	} else {
		sscanf(line, "bounces: %u\n", &nb_bounces);
	}

	

	fscanf(f, "%[^\n]\n", line);
	read = sscanf(line, "has_denoiser: %u\n", &nbo);
	if (read == 1) {
		has_denoiser = nbo;
		fscanf(f, "intensite_lum: %f\n", &s.intensite_lumiere);
	} else {
		sscanf(line, "intensite_lum: %f\n", &s.intensite_lumiere);
	}
	fscanf(f, "intensite_envmap: %f\n", &s.envmap_intensity);


	fscanf(f, "%[^\n]\n", line);
	if (line[0] == 'n') {
		sscanf(line, "nbobjects: %u\n", &nbo);
		s.clear_background();
	} else {
		sscanf(line, "background: %[^\n]\n", bg);
		s.load_background(bg, gamma);
		fscanf(f, "nbobjects: %u\n", &nbo);
	}

	for (int i = 0; i < nbo; i++) {
		s.addObject(Object::create_from_file(f, &s, replacedNames));
	}

	fscanf(f, "fog_density: %f\n", &s.fog_density);
	fscanf(f, "%[^\n]\n", line);
	read = sscanf(line, "fog_absorption: %f\n", &s.fog_absorption);
	if (read == 1) {
		fscanf(f, "fog_density_decay: %f\n", &s.fog_density_decay);
		fscanf(f, "fog_absorption_decay: %f\n", &s.fog_absorption_decay);
		fscanf(f, "fog_type: %u\n", &s.fog_type);
	} else {
		sscanf(line, "fog_type: %u\n", &s.fog_type);
	}
	
	fscanf(f, "%[^\n]\n", line);
	read = sscanf(line, "fog_phase_type: %u\n", &s.fog_phase_type);
	if (read == 1)
		fscanf(f, "double_frustum_start_t: %f\n", &s.double_frustum_start_t);

	fclose(f);

	s.lumiere = dynamic_cast<Sphere*>(s.objects[0]);
}

void Raytracer::loadScene() {

	//W = 1000;
	//H = 1200;
#ifdef _DEBUG
	W = 320;
	H = 400;
#else
	W = 1000;
	H = 800;
#endif
	nrays = 100;   //10 pour debugguer, 100 pour les rendus finaux	
	cam = Camera(Vector(0., 0., 50.), Vector(0, 0, -1), Vector(0, 1, 0));
	cam.fov = 35 * M_PI / 180;   // field of view
	cam.focus_distance = 50; // distance mise au point	
	cam.aperture = 0.1;     // ouverture, pour depth of field (faible valeur = net partout)
	sigma_filter = 0.5;  // antialiasing
	nb_bounces = 3;

	Sphere* slum = new Sphere(Vector(10, 23, 15), 10);   //sphere light: position, rayon, albedo (qui sert à rien)
	Sphere* s2 = new Sphere(Vector(0, 0, 0), 1000000); s2->flip_normals = true; // envmap  

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

float sum_area_table(float* sat, int sat_width, int i0, int i1, int j0, int j1) {

	float term1 = 0;
	if (i0 > 0) {
		term1 = sat[(i0 - 1)*sat_width + j1];
	}
	float term2 = 0;
	if (j0 > 0) {
		term2 = sat[i1*sat_width + j0-1];
	}
	float term3 = 0;
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


inline uint32_t ReverseBits(uint32_t n) {
	n = (n << 16) | (n >> 16);
	n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
	n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
	n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
	n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);
	return n;
}

Vector extensibleLattice2d(uint32_t id) {
	uint32_t rid = ReverseBits(id);
	float phi_id = rid * pow(2.0, -32);

	float tmp;
	float x = modf(phi_id * 1 + 0.456789123, &tmp);       // sequence: 1, 182667, 469891, 498753, 110745, 446247, 250185, 118627, 245333, 283199, 408519, 391023, 246327
	float y = modf(phi_id * 182667 + 0.123456789, &tmp);  // see https://web.maths.unsw.edu.au/~fkuo/lattice/index.html for more ( lattice32001_order2 )
	return Vector(x, y, 0);
}

void Raytracer::prepare_render(float time) {

	sphereEnv = dynamic_cast<Sphere*>(s.objects[1]);

	for (int i = 0; i < omp_get_max_threads(); i++) {
		engine[i] = pcg32(i);
	}
	if (randomPerPixel.size() != W * H) {
		Wlr = ceil(W / 16.f);
		Hlr = ceil(H / 16.f);
		computed.resize(W*H, false);
		image.resize(W*H * 3, 0);
		imagedouble.resize(W*H * 3, 0);
		imagedouble_lowres.resize(Wlr*Hlr * 3, 0);
		filteredImage.resize(W*H * 3, 0);
		albedoImage.resize(W*H * 3, 0);
		normalImage.resize(W*H * 3, 0);
		sample_count.resize(W*H, 0);

		randomPerPixel.resize(W*H);
		for (int i = 0; i < W*H; i++) {
			randomPerPixel[i][0] = engine[0]()*invmax;
			randomPerPixel[i][1] = engine[0]()*invmax;
		}
	}
	if (nrays != last_nrays) {
		samples2d.resize(nrays);
		for (int i = 0; i < nrays; i++) {
			samples2d[i] = extensibleLattice2d(i);
		}
		last_nrays = nrays;
	}

	if (sigma_filter != lastfilter) {

		filter_size = std::ceil(sigma_filter * 2);
		filter_total_width = 2 * filter_size + 1;
		filter_integral.resize(filter_total_width*filter_total_width);
		filter_value.resize(filter_total_width*filter_total_width);
		for (int i = -filter_size; i <= filter_size; i++) {
			for (int j = -filter_size; j <= filter_size; j++) {
				float integ = 0;
				for (int i2 = -filter_size; i2 <= i; i2++) {
					for (int j2 = -filter_size; j2 <= j; j2++) {
						float w = fast_exp(-(i2*i2 + j2 * j2) / (2.*sigma_filter*sigma_filter)) / (sigma_filter*sigma_filter*2.*M_PI);
						integ += w;
					}
				}
				filter_integral[(i + filter_size)*filter_total_width + (j + filter_size)] = integ;
				filter_value[(i + filter_size)*filter_total_width + (j + filter_size)] = exp(-(i*i + j * j) / (2.*sigma_filter*sigma_filter)) / (sigma_filter*sigma_filter*2.*M_PI);
			}
		}
		lastfilter = sigma_filter;
	}
	s.prepare_render(is_recording);

	centerLight = s.lumiere->apply_transformation(s.lumiere->O);// s.lumiere->O + s.lumiere->get_translation(r.time, is_recording);
	lum_scale = s.lumiere->get_scale(time, is_recording);
	radiusLight = lum_scale * s.lumiere->R;
	lightPower = s.intensite_lumiere / sqr(lum_scale);

	std::fill(computed.begin(), computed.end(), false);
	std::fill(sample_count.begin(), sample_count.end(), 0);
	std::fill(imagedouble_lowres.begin(), imagedouble_lowres.end(), 0);
	std::fill(imagedouble.begin(), imagedouble.end(), 0);

	std::fill(filteredImage.begin(), filteredImage.end(), 0);
	std::fill(albedoImage.begin(), albedoImage.end(), 0);
	std::fill(normalImage.begin(), normalImage.end(), 0);
	
}

void Raytracer::precomputeRayBatch(int i1, int j1, int i2, int j2, int nspp) {  // i1 included to i2 excluded
	int threadid = omp_get_thread_num();
	int batchW = j2 - j1;
	int batchH = i2 - i1;
	int batchSize = batchW * batchH*nspp;

	// compute the first intersection with the whole set of first rays
	s.firstIntersection_Ray[threadid].resize(batchSize);
	s.firstIntersection_dx[threadid].resize(batchSize);
	s.firstIntersection_dy[threadid].resize(batchSize);

	for (int i = 0; i < batchH; i++) {
		for (int j = 0; j < batchW; j++) {
			for (int k = 0; k < nspp; k++) {
				int rayId = (i*batchW + j)*nspp + k;
				float dx = engine[threadid]()*invmax - 0.5f;  // not perfectly uniform but much faster than std::uniform
				float dy = engine[threadid]()*invmax - 0.5f;
				s.firstIntersection_dx[threadid][rayId] = dx;
				s.firstIntersection_dy[threadid][rayId] = dy;

				float dx_aperture = (engine[threadid]()*invmax - 0.5f) * cam.aperture;
				float dy_aperture = (engine[threadid]()*invmax - 0.5f) * cam.aperture;

				Ray r = cam.generateDirection(s.double_frustum_start_t, i+i1, j+j1, s.current_frame, dx, dy, dx_aperture, dy_aperture, W, H);
				s.firstIntersection_Ray[threadid][rayId] = r;
			}
		}
	}
	s.first_intersection_batch(batchSize);
}

void Raytracer::render_image()
{

	int fstart = 0;
	int fend = 1;

	float denom2 = 1.f / (2.*sigma_filter*sigma_filter);

	/*double test = sum_area_table(&filter_integral[0], filter_total_width, filter_size- filter_size, filter_size+15, filter_size+8, filter_size+11);
	double integ = 0;
	for (int i2 = -filter_size; i2 <= 15; i2++) {
		for (int j2 = 8; j2 <= 11; j2++) {
			double w = exp(-(i2*i2 + j2*j2) / (2.*sigma_filter*sigma_filter)) / (sigma_filter*sigma_filter*2.*M_PI);
			integ += w;
		}
	}*/

	prepare_render(s.current_frame);
	
	realtime_ray_iter = 0;
		for ( ; realtime_ray_iter < nrays; realtime_ray_iter++) {
			current_nb_rays = realtime_ray_iter;
			chrono.Start();
			for(int pix_in_block=0; pix_in_block<64; pix_in_block++) {
				int i1 = pix_in_block % 8;// permut[pix_in_block].first;
				int j1 = pix_in_block / 8; // permut[pix_in_block].second;
			//for (int i1 = 0; i1 < 8; i1++) {
				//for (int j1 = 0; j1 < 8; j1++) {
					if (stopped) {						
						return;
					}
#pragma omp parallel for schedule(dynamic, 1)
					for (int i = i1; i < H; i += 8) {
						int threadid = omp_get_thread_num();

						for (int j = j1; j < W; j += 8) {


							float dx = engine[threadid]()*invmax - 0.5f;  // not perfectly uniform but much faster than std::uniform
							float dy = engine[threadid]()*invmax - 0.5f;

							float dx_aperture = (engine[threadid]()*invmax - 0.5f) * cam.aperture;
							float dy_aperture = (engine[threadid]()*invmax - 0.5f) * cam.aperture;

							float time = s.current_frame;// + engine[threadid]()*invmax;



							Ray r = cam.generateDirection(s.double_frustum_start_t, i, j, time, dx, dy, dx_aperture, dy_aperture, W, H);

							Vector normal, albedo;
							Vector color = getColor(r, realtime_ray_iter, nb_bounces, i, j, normal, albedo, false, false);

							int bmin_i = std::max(0, i - filter_size);
							int bmax_i = std::min(i + filter_size, H - 1);
							int bmin_j = std::max(0, j - filter_size);
							int bmax_j = std::min(j + filter_size, W - 1);
							float ratio = 1.f / sum_area_table(&filter_integral[0], filter_total_width, bmin_i - i + filter_size, bmax_i - i + filter_size, bmin_j - j + filter_size, bmax_j - j + filter_size);
							float denom1 = ratio / (sigma_filter*sigma_filter*2.*M_PI);
							//double dx = s.firstIntersection_dx[i*W + j];
							//double dy = s.firstIntersection_dy[i*W + j];

							for (int i2 = bmin_i; i2 <= bmax_i; i2++) {
								for (int j2 = bmin_j; j2 <= bmax_j; j2++) {
									//double w = filter_value[(i2 - i + filter_size)*filter_total_width + j2-j + filter_size] * ratio; 

									float w =  fast_exp(-(sqr(i2 - i - dy) + sqr(j2 - j - dx)) *denom2) *denom1;
									imagedouble[((H - i2 - 1)*W + j2) * 3 + 0] += color[0] * w;
									imagedouble[((H - i2 - 1)*W + j2) * 3 + 1] += color[1] * w;
									imagedouble[((H - i2 - 1)*W + j2) * 3 + 2] += color[2] * w;																		
									//sw += w;
									sample_count[(H - i2 - 1)*W + j2] += w; // 1. / sqr(filter_total_width);
								}
							}

							/*imagedouble[((H - i - 1)*W + j) * 3 + 0] += color[0] ;
							imagedouble[((H - i - 1)*W + j) * 3 + 1] += color[1] ;
							imagedouble[((H - i - 1)*W + j) * 3 + 2] += color[2] ;
							sample_count[(H - i - 1)*W + j]++;*/
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


#pragma omp parallel for 
		for (int i = 0; i < H; i++) {
			for (int j = 0; j < W; j++) {
				image[((H - i - 1)*W + j) * 3 + 0] = std::min(255., std::max(0., 255.*std::pow(imagedouble[((H - i - 1)*W + j) * 3 + 0] / 196964.7 / std::max(sample_count[(H - i - 1)*W + j],1.f), 1 / gamma)));   // rouge
				image[((H - i - 1)*W + j) * 3 + 1] = std::min(255., std::max(0., 255.*std::pow(imagedouble[((H - i - 1)*W + j) * 3 + 1] / 196964.7 / std::max(sample_count[(H - i - 1)*W + j],1.f), 1 / gamma))); // vert
				image[((H - i - 1)*W + j) * 3 + 2] = std::min(255., std::max(0., 255.*std::pow(imagedouble[((H - i - 1)*W + j) * 3 + 2] / 196964.7 / std::max(sample_count[(H - i - 1)*W + j],1.f), 1 / gamma))); // bleu
			}
		}

		if (autosave) {
			std::ostringstream os;
			//os << "testFog_" << time_step << ".bmp";
			if (cam.isArray) {
				os << "exportD" << s.current_frame << "_" << cam.current_viewX << "_" << cam.nbviewX << "_" << cam.current_viewY << "_" << cam.nbviewY << ".jpg";
			} else {
				os << "exportD" << s.current_frame << ".jpg";
			}
			save_image(os.str().c_str(), &image[0], W, H);
		}
	stopped = true;

	//std::cout << chrono.GetDiffMs() / (double)1000. << std::endl;
    return;
}

void Raytracer::render_image_nopreviz() {

	prepare_render(s.current_frame);
	float denom2 = 1.f / (2.f*sigma_filter*sigma_filter);
	const int batchWidth = 4;
	const int batchHeight = 4;
	const int nbBatchX = ceil(W / (float)batchWidth);
	const int nbBatchY = ceil(H / (float)batchHeight);

	chrono.Start();
	int maxthreads = omp_get_max_threads();
	std::vector<float> imagedoublethreads(W*H * 3 * maxthreads, 0);
	std::vector<float> samplecountthreads(W*H  * maxthreads, 0);
	std::vector<float> albedodoublethreads(W*H * 3 * maxthreads, 0);
	std::vector<float> normaldoublethreads(W*H * 3 * maxthreads, 0);

#pragma omp parallel 
	{
		int threadid = omp_get_thread_num();
		float* curimagedouble = &imagedoublethreads[threadid*W*H * 3];
		float* cursamplecount = &samplecountthreads[threadid*W*H];
		float* curalbedodouble = &albedodoublethreads[threadid*W*H * 3];
		float* curnormaldouble = &normaldoublethreads[threadid*W*H * 3];


#pragma omp for schedule(dynamic, 1) 
		for (int batchid = 0; batchid < nbBatchX*nbBatchY; batchid++) {
			int batchi = batchid / nbBatchX;
			int batchj = batchid % nbBatchX;

			precomputeRayBatch(batchi*batchHeight, batchj*batchWidth, std::min(H, batchi*batchHeight + batchHeight), std::min(W, batchj*batchWidth + batchWidth), nrays);
			int batchW = std::min(W, batchj*batchWidth + batchWidth) - batchj* batchWidth;
			int batchH = std::min(H, batchi*batchHeight + batchHeight) - batchi* batchHeight;

			for (int id = 0; id < batchW*batchH; id++) {

				int i = batchi * batchHeight + id / batchW;
				int j = batchj * batchWidth + id % batchW;

				int bmin_i = std::max(0, i - filter_size);
				int bmax_i = std::min(i + filter_size, H - 1);
				int bmin_j = std::max(0, j - filter_size);
				int bmax_j = std::min(j + filter_size, W - 1);
				float ratio = 1.f / sum_area_table(&filter_integral[0], filter_total_width, bmin_i - i + filter_size, bmax_i - i + filter_size, bmin_j - j + filter_size, bmax_j - j + filter_size);
				float denom1 = ratio / (sigma_filter*sigma_filter*2.*M_PI);

				for (int k = 0; k < nrays; k++) {

					/*float dx = engine[threadid]()*invmax - 0.5f;
					float dy = engine[threadid]()*invmax - 0.5f;

					float dx_aperture = (engine[threadid]()*invmax - 0.5f) * cam.aperture;
					float dy_aperture = (engine[threadid]()*invmax - 0.5f) * cam.aperture;

					float time = s.current_frame;//  +engine[threadid]()*invmax; //current_time // beware, now everything is in number of frames! ; deactivated motion blur as the behavior is strange when the light radius is changing over time


					Ray r = cam.generateDirection(s.double_frustum_start_t, i, j, time, dx, dy, dx_aperture, dy_aperture, W, H);*/

					const Ray &r = s.firstIntersection_Ray[threadid][id*nrays+k];
					float dx = s.firstIntersection_dx[threadid][id*nrays + k];
					float dy = s.firstIntersection_dy[threadid][id*nrays + k];

					Vector normal, albedo;
					Vector color = getColor(r, k, nb_bounces, i, j, normal, albedo, false, true, id*nrays + k);

					if (has_denoiser) {
						int idx = ((H - i - 1)*W + j) * 3;
						curimagedouble[idx + 0] += color[0];
						curimagedouble[idx + 1] += color[1];
						curimagedouble[idx + 2] += color[2];
						//sw += w;
						cursamplecount[(H - i - 1)*W + j] += 1; // 1. / sqr(filter_total_width);

						curnormaldouble[idx + 0] += normal[0];
						curnormaldouble[idx + 1] += normal[1];
						curnormaldouble[idx + 2] += normal[2];

						curalbedodouble[idx + 0] += albedo[0];
						curalbedodouble[idx + 1] += albedo[1];
						curalbedodouble[idx + 2] += albedo[2];
					} else { // Intel Open Image Denoiser does not support splatting (the result appears unfiltered)
						for (int i2 = bmin_i; i2 <= bmax_i; i2++) {
							for (int j2 = bmin_j; j2 <= bmax_j; j2++) {
								//double w = filter_value[(i2 - i + filter_size)*filter_total_width + j2-j + filter_size] * ratio; 

								int idx = ((H - i2 - 1)*W + j2) * 3;
								float w = fast_exp(-(sqr(i2 - i - dy) + sqr(j2 - j - dx)) *denom2) *denom1;
								curimagedouble[idx + 0] += color[0] * w;
								curimagedouble[idx + 1] += color[1] * w;
								curimagedouble[idx + 2] += color[2] * w;
								//sw += w;
								cursamplecount[(H - i2 - 1)*W + j2] += w; // 1. / sqr(filter_total_width);
							}
						}
					}
				}



			}
		}
	}

	for (int th = 0; th < maxthreads; th++) {
		for (int i = 0; i < W*H; i++) {
			imagedouble[i * 3] += imagedoublethreads[th*W*H * 3 + i * 3];
			imagedouble[i * 3+1] += imagedoublethreads[th*W*H * 3 + i * 3+1];
			imagedouble[i * 3+2] += imagedoublethreads[th*W*H * 3 + i * 3+2];
			sample_count[i] += samplecountthreads[th*W*H + i];

			if (has_denoiser) {
				albedoImage[i * 3] += albedodoublethreads[th*W*H * 3 + i * 3];
				albedoImage[i * 3 + 1] += albedodoublethreads[th*W*H * 3 + i * 3 + 1];
				albedoImage[i * 3 + 2] += albedodoublethreads[th*W*H * 3 + i * 3 + 2];
				normalImage[i * 3] += imagedoublethreads[th*W*H * 3 + i * 3];
				normalImage[i * 3 + 1] += imagedoublethreads[th*W*H * 3 + i * 3 + 1];
				normalImage[i * 3 + 2] += imagedoublethreads[th*W*H * 3 + i * 3 + 2];
			}		
		}
	}

	for (int i = 0; i < W*H; i++) {
		float normnormal = sqrt(normalImage[i * 3 + 0] * normalImage[i * 3 + 0] + normalImage[i * 3 + 1] * normalImage[i * 3 + 1] + normalImage[i * 3 + 2] * normalImage[i * 3 + 2]);
		for (int j = 0; j < 3; j++) {
			imagedouble[i*3+j] /= sample_count[i];
			albedoImage[i * 3 + j] /= sample_count[i];
			normalImage[i * 3 + j] /= normnormal;
		}				
	}

	curTimePerFrame = chrono.GetDiffMs() / nrays;

	std::ostringstream os;

	if (!has_denoiser) {
		#pragma omp parallel for
			for (int i = 0; i < H; i++) {
				for (int j = 0; j < W; j++) {
					image[((H - i - 1)*W + j) * 3 + 0] = std::min(255., std::max(0., 255.*std::pow(imagedouble[((H - i - 1)*W + j) * 3 + 0] / 196964.7, 1 / gamma)));   // rouge
					image[((H - i - 1)*W + j) * 3 + 1] = std::min(255., std::max(0., 255.*std::pow(imagedouble[((H - i - 1)*W + j) * 3 + 1] / 196964.7, 1 / gamma))); // vert
					image[((H - i - 1)*W + j) * 3 + 2] = std::min(255., std::max(0., 255.*std::pow(imagedouble[((H - i - 1)*W + j) * 3 + 2] / 196964.7, 1 / gamma))); // bleu
				}
			}

			//os << "testFog_" << time_step << ".bmp";
			if (autosave) {
				if (cam.isArray) {
					os << "exportE" << s.current_frame << "_" << cam.current_viewX << "_" << cam.nbviewX << "_" << cam.current_viewY << "_" << cam.nbviewY << ".jpg";
				} else {
					os << "exportE" << s.current_frame << ".jpg";
				}
				save_image(os.str().c_str(), &image[0], W, H);
			}
	} else {

#ifdef USE_OPENIMAGEDENOISER
		chrono.Start();
		filter.setImage("color", &imagedouble[0], oidn::Format::Float3, W, H); // beauty
		filter.setImage("albedo", &albedoImage[0], oidn::Format::Float3, W, H); // auxiliary
		filter.setImage("normal", &normalImage[0], oidn::Format::Float3, W, H); // auxiliary
		filter.setImage("output", &filteredImage[0], oidn::Format::Float3, W, H); // denoised beauty
		filter.set("inputScale", 1.f / 196964.7f);
		//filter.set("cleanAux", true);
		filter.set("hdr", true); // beauty image is HDR
		filter.commit();
		std::cout << "commit time (ms) : " << chrono.GetDiffMs() << std::endl;;
		chrono.Start();
		// Filter the image
		filter.execute();
		const char* errorMessage;
		if (device.getError(errorMessage) != oidn::Error::None)
			std::cout << "Error: " << errorMessage << std::endl;
		std::cout << "filtering time (ms) : " << chrono.GetDiffMs() << std::endl;;
#pragma omp parallel for 
		for (int i = 0; i < H; i++) {
			for (int j = 0; j < W; j++) {
				image[((H - i - 1)*W + j) * 3 + 0] = std::min(255., std::max(0., 255.*std::pow(filteredImage[((H - i - 1)*W + j) * 3 + 0] / 196964.7, 1 / gamma)));   // rouge
				image[((H - i - 1)*W + j) * 3 + 1] = std::min(255., std::max(0., 255.*std::pow(filteredImage[((H - i - 1)*W + j) * 3 + 1] / 196964.7, 1 / gamma))); // vert
				image[((H - i - 1)*W + j) * 3 + 2] = std::min(255., std::max(0., 255.*std::pow(filteredImage[((H - i - 1)*W + j) * 3 + 2] / 196964.7, 1 / gamma))); // bleu
			}
		}
		os.str("");
		//os << "testFog_" << time_step << ".bmp";
		if (autosave) {
			if (cam.isArray) {
				os << "exportEFiltered" << s.current_frame << "_" << cam.current_viewX << "_" << cam.nbviewX << "_" << cam.current_viewY << "_" << cam.nbviewY << ".jpg";
			} else {
				os << "exportEFiltered" << s.current_frame << ".jpg";
			}
			save_image(os.str().c_str(), &image[0], W, H);
		}


	}


/*#pragma omp parallel for 
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			image[((H - i - 1)*W + j) * 3 + 0] = std::min(255., std::max(0., 255.*std::pow(albedoImage[((H - i - 1)*W + j) * 3 + 0], 1 / gamma)));   // rouge
			image[((H - i - 1)*W + j) * 3 + 1] = std::min(255., std::max(0., 255.*std::pow(albedoImage[((H - i - 1)*W + j) * 3 + 1], 1 / gamma))); // vert
			image[((H - i - 1)*W + j) * 3 + 2] = std::min(255., std::max(0., 255.*std::pow(albedoImage[((H - i - 1)*W + j) * 3 + 2], 1 / gamma))); // bleu
		}
	}
	os.str("");
	//os << "testFog_" << time_step << ".bmp";
	if (cam.isArray) {
		os << "exportEAlbedo" << s.current_frame << "_" << cam.current_viewX << "_" << cam.nbviewX << "_" << cam.current_viewY << "_" << cam.nbviewY << ".jpg";
	} else {
		os << "exportEAlbedo" << s.current_frame << ".jpg";
	}
	save_image(os.str().c_str(), &image[0], W, H);



#pragma omp parallel for 
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			image[((H - i - 1)*W + j) * 3 + 0] = std::min(255., std::max(0., 255.*std::pow((1+normalImage[((H - i - 1)*W + j) * 3 + 0] )/2, 1 / gamma)));   // rouge
			image[((H - i - 1)*W + j) * 3 + 1] = std::min(255., std::max(0., 255.*std::pow((1+normalImage[((H - i - 1)*W + j) * 3 + 1] )/2, 1 / gamma))); // vert
			image[((H - i - 1)*W + j) * 3 + 2] = std::min(255., std::max(0., 255.*std::pow((1+normalImage[((H - i - 1)*W + j) * 3 + 2] )/2, 1 / gamma))); // bleu
		}
	}
	os.str("");
	//os << "testFog_" << time_step << ".bmp";
	if (cam.isArray) {
		os << "exportENormals" << s.current_frame << "_" << cam.current_viewX << "_" << cam.nbviewX << "_" << cam.current_viewY << "_" << cam.nbviewY << ".jpg";
	} else {
		os << "exportENormals" << s.current_frame << ".jpg";
	}
	save_image(os.str().c_str(), &image[0], W, H);*/
#endif
}