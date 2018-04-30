#pragma once

#include <vector>
#include "Vector.h"
#include "Geometry.h"
#include "TriangleMesh.h"
#include "PointSet.h"



class Raytracer {
public:

	Raytracer() { stopped = false; current_nb_rays = 0; curTimePerFrame = 0; };
	void loadScene();
	void render_image();
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
	};
	void stopRender() { stopped = true; };
	void save_scene(const char* filename);
	void load_scene(const char* filename);

	Vector getColor(const Ray &r, const Scene &s, int nbrebonds, int screenI, int screenJ, bool show_lights = true);

	int W, H, Wlr, Hlr;
	int nrays;  
	Camera cam;   // camera: position, direction de visée, up vector
	
	double sigma_filter;  // antialiasing
	int nb_bounces;
	
	Scene s;

	bool stopped;
	int current_nb_rays;
	std::vector<unsigned char> image;
	std::vector<double> imagedouble;
	std::vector<bool> computed;
	double curTimePerFrame;
	PerfChrono chrono;
	std::vector<double> filter_value;
	std::vector<double> filter_integral;

	std::vector<double> filter_value_lowres;
	std::vector<double> filter_integral_lowres;


	std::vector<double> imagedouble_lowres;
	std::vector<double> sample_count;
};

void save_img(const char* filename, const unsigned char* pixels, int W, int H);