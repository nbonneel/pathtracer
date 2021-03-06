#pragma once

#include <vector>
#include "Vector.h"
#include "Geometry.h"
#include "TriangleMesh.h"
#include "PointSet.h"

#ifdef USE_OPENIMAGEDENOISER
#include <OpenImageDenoise/oidn.hpp>
#endif

class Contrib {
public:
	Contrib() {};
	Contrib(const Vector& w, const Ray& r, int d, bool showlights, bool hadSS) :weight(w), r(r), depth(d), show_lights(showlights), has_had_subsurface_interaction(hadSS) {};
	Vector weight;
	Ray r;
	int depth;
	bool show_lights, has_had_subsurface_interaction;
};

class Raytracer {
public:

	Raytracer() : invmax(1.f / engine[0].max())
	{
		stopped = false; current_nb_rays = 0; curTimePerFrame = 0; is_recording = false;  gamma = 2.2;
		for (int i = 0; i < omp_get_max_threads(); i++) {
			engine[i] = pcg32(i);
		}
		has_denoiser = false;
		autosave = true;
#ifdef USE_OPENIMAGEDENOISER
		device = oidn::newDevice();
		device.commit();
		filter = device.newFilter("RT");
#endif
	};
	void loadScene();
	void prepare_render(double time);
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

	Vector getColor(const Ray &r, int sampleID, int nbrebonds, int screenI, int screenJ, Vector &normalValue, Vector &albedoValue, bool no_envmap = false);
	bool fogContribution(const Ray &r, const Vector& sampleLightPos, double t, Vector curWeight, int nbrebonds, bool showLight, bool hadSS, Contrib& newContrib, double &attenuationFactor);

	int W, H, Wlr, Hlr;
	int nrays, last_nrays;  
	std::vector<std::pair<int, int> > permut;
	Camera cam;   // camera: position, direction de vis�e, up vector
	
	double sigma_filter, lastfilter;  // antialiasing
	int filter_size, filter_total_width;
	int nb_bounces;
	double gamma;
	bool is_recording;
	
	Scene s;

	bool stopped;
	bool has_denoiser, autosave;
	int current_nb_rays;
	std::vector<unsigned char> image;
	std::vector<float> imagedouble;
	std::vector<float> albedoImage, normalImage, filteredImage;	
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

	Vector centerLight;
	double lum_scale, radiusLight, lightPower;

#ifdef USE_OPENIMAGEDENOISER
	oidn::DeviceRef device;
	oidn::FilterRef filter;
#endif
};

