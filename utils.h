#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <string>
#include "hdrwriter.h"
#include <vector>
#include "CImg.h"

bool file_exists(const char* filename);
std::string& str_replace(const std::string &search, const std::string &replace, std::string &subject);
std::string extractFilePath (const std::string& str);
std::string extractFileName (const std::string& str);
std::string extractFilePathWithEndingSlash(const std::string& str);



template<typename T>
void save_image(const char* filename, int num, int W, int H, const T* val, T maxval, int nbchan) {

	if (nbchan == 3) {

		std::string s(filename);
		if (s.find(".hdr") != std::string::npos) {
			std::string ret = EncodeRadianceHDR(W, H, (const float*)val, nbchan);
			FILE* f = fopen(filename, "wb+");
			fwrite(&ret[0], sizeof(ret[0]), ret.size(), f);
			fclose(f);
			return;
		}


		std::vector<unsigned char> deinterleaved(W*H * 3);
		for (int i = 0; i<W*H; i++) {
			deinterleaved[i] = min(255., max(0., val[i * 3] * (255. / maxval)));
			deinterleaved[i + W*H] = min(255., max(0., val[i * 3 + 1] * (255. / maxval)));
			deinterleaved[i + 2 * W*H] = min(255., max(0., val[i * 3 + 2] * (255. / maxval)));
		}

		cimg_library::CImg<unsigned char> cimg(&deinterleaved[0], W, H, 1, 3);
		cimg.save(filename, num);
	} else {

		std::vector<unsigned char> deinterleaved(W*H * 3);
		for (int i = 0; i<W*H; i++) {
			deinterleaved[i] = min(255., max(0., val[i] * (255. / maxval)));
			deinterleaved[i + W*H] = deinterleaved[i];
			deinterleaved[i + 2 * W*H] = deinterleaved[i];
		}

		cimg_library::CImg<unsigned char> cimg(&deinterleaved[0], W, H, 1, 3);
		cimg.save(filename, num);
	}
}
template<typename T>
void save_image(const char* filename, int num, int W, int H, const std::vector<T> &val, T maxval) {

	if (val.size() == W*H * 3) {
		std::string s(filename);
		if (s.find(".hdr") != std::string::npos && sizeof(T) == sizeof(float)) {
			std::string ret = EncodeRadianceHDR(W, H, (const float*)&val[0], 3);
			FILE* f = fopen(filename, "wb+");
			fwrite(&ret[0], sizeof(ret[0]), ret.size(), f);
			fclose(f);
			return;
		}
		std::vector<unsigned char> deinterleaved(W*H * 3);
		for (int i = 0; i < W*H; i++) {
			deinterleaved[i] = min(255., max(0., val[i * 3] * (255. / maxval)));
			deinterleaved[i + W*H] = min(255., max(0., val[i * 3 + 1] * (255. / maxval)));
			deinterleaved[i + 2 * W*H] = min(255., max(0., val[i * 3 + 2] * (255. / maxval)));
		}

		cimg_library::CImg<unsigned char> cimg(&deinterleaved[0], W, H, 1, 3);
		cimg.save(filename, num);
	} else {
		if (val.size() == W*H * 2) {
			std::vector<unsigned char> deinterleaved(W*H * 3);
			for (int i = 0; i < W*H; i++) {
				deinterleaved[i] = min(255., max(0., val[i * 2] * (255. / maxval)));
				deinterleaved[i + W*H] = min(255., max(0., val[i * 2 + 1] * (255. / maxval)));
				deinterleaved[i + 2 * W*H] = min(255., max(0., 0 * (255. / maxval)));
			}

			cimg_library::CImg<unsigned char> cimg(&deinterleaved[0], W, H, 1, 3);
			cimg.save(filename, num);
		} else {
			std::vector<unsigned char> deinterleaved(W*H * 3);
			for (int i = 0; i < W*H; i++) {
				deinterleaved[i] = min(255., max(0., val[i] * (255. / maxval)));
				deinterleaved[i + W*H] = deinterleaved[i];
				deinterleaved[i + 2 * W*H] = deinterleaved[i];
			}

			cimg_library::CImg<unsigned char> cimg(&deinterleaved[0], W, H, 1, 3);
			cimg.save(filename, num);
		}
	}
}


int get_png_bitdepth(const char* filename);

template<typename T>
bool load_image(const char* filename, std::vector<T> &val, size_t &W, size_t &H, bool bits16 = false) {

	if (!file_exists(filename)) return false;
	if (bits16 || (std::string(filename).find('.png') != std::string::npos && get_png_bitdepth(filename) == 16)) {
		cimg_library::CImg<float> cimg(filename);
		W = cimg.width();
		H = cimg.height();
		val.resize(W*H * 3);

		for (int i = 0; i<W*H; i++) {
			val[i * 3] = cimg.data()[i]  / 150.;
			val[i * 3 + 1] = cimg.data()[i + W*H] / 150.;
			val[i * 3 + 2] = cimg.data()[i + W*H * 2] / 150.;
		}
		return true;
	}
	cimg_library::CImg<unsigned char> cimg(filename);
	W = cimg.width();
	H = cimg.height();
	if (cimg.size() == W*H) {
		val.resize(W*H * 3);
		for (int i = 0; i<W*H; i++) {
			val[i * 3] = cimg.data()[i] ;
			val[i * 3 + 1] = cimg.data()[i] ;
			val[i * 3 + 2] = cimg.data()[i] ;
		}
		for (int i = 0; i < H / 2; i++) {
			for (int j = 0; j < W; j++) {
				std::swap(val[((H - i - 1)*W + j) * 3 + 0], val[(i*W + j) * 3 + 0]);
				std::swap(val[((H - i - 1)*W + j) * 3 + 1], val[(i*W + j) * 3 + 1]);
				std::swap(val[((H - i - 1)*W + j) * 3 + 2], val[(i*W + j) * 3 + 2]);
			}
		}
		return true;
	}
	val.resize(W*H * 3);

	for (int i = 0; i<W*H; i++) {
		val[i * 3] = cimg.data()[i] ;
		val[i * 3 + 1] = cimg.data()[i + W*H] ;
		val[i * 3 + 2] = cimg.data()[i + W*H * 2];
	}
	for (int i = 0; i < H/2; i++) {
		for (int j = 0; j < W; j++) {
			std::swap(val[((H - i - 1)*W + j) * 3 + 0], val[(i*W + j) * 3 + 0]);
			std::swap(val[((H - i - 1)*W + j) * 3 + 1], val[(i*W + j) * 3 + 1]);
			std::swap(val[((H - i - 1)*W + j) * 3 + 2], val[(i*W + j) * 3 + 2]);
		}
	}
	return true;
}