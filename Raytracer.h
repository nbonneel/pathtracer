#pragma once

#include <vector>
#include "Vector.h"
#include "Geometry.h"
#include "TriangleMesh.h"
#include "PointSet.h"



class Raytracer {
public:

	Raytracer() : invmax(1.f / engine[0].max())
	{
		stopped = false; current_nb_rays = 0; curTimePerFrame = 0; is_recording = false;  gamma = 2.2;
		for (int i = 0; i < omp_get_max_threads(); i++) {
			engine[i] = pcg32(i);
		}
	};
	void loadScene();
	void prepare_render();
	void render_image();
	void render_image_nopreviz();
	void clear_image() { 

		image.resize(W*H * 3, 0);
		imagedouble.resize(W*H * 3, 0.);
		computed.resize(W*H, false);
		std::fill(imagedouble.begin(), imagedouble.end(), 0);
		std::fill(image.begin(), image.end(), 0);
		std::fill(computed.begin(), computed.end(), false);

		Wlr = ceil(W / 16.);
		Hlr = ceil(H / 16.);
		imagedouble_lowres.resize(Wlr*Hlr * 3);

		permut.resize(64);
		for (int i = 0; i < 8; i++)
			for (int j = 0; j < 8; j++)
				permut[i * 8 + j] = std::make_pair(i, j);
		//std::random_shuffle(permut.begin(), permut.end());
	};
	void stopRender() { stopped = true; };
	void save_scene(const char* filename);
	void load_scene(const char* filename, const char* replacedNames = NULL);

	Vector getColor(const Ray &r, int sampleID, int nbrebonds, int screenI, int screenJ, bool show_lights = true, bool no_envmap = false, bool has_had_subsurface_interaction = false);

	int W, H, Wlr, Hlr;
	int nrays, last_nrays;  
	std::vector<std::pair<int, int> > permut;
	Camera cam;   // camera: position, direction de visée, up vector
	
	double sigma_filter, lastfilter;  // antialiasing
	int filter_size, filter_total_width;
	int nb_bounces;
	double gamma;
	bool is_recording;
	
	Scene s;

	bool stopped;
	int current_nb_rays;
	std::vector<unsigned char> image;
	std::vector<double> imagedouble;
	std::vector<bool> computed;
	double curTimePerFrame;
	PerfChrono chrono;
	const float invmax;
	std::vector<double> filter_value;
	std::vector<double> filter_integral;

	std::vector<double> filter_value_lowres;
	std::vector<double> filter_integral_lowres;


	std::vector<double> imagedouble_lowres;
	std::vector<double> sample_count;

	std::vector<Vector> samples2d;
	std::vector<Vector> randomPerPixel;
};

