#include "utils.h"
#include <fstream>
#include <sstream>
#define cimg_display 0
#include "CImg.h"
#include "hdrwriter.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


std::string extractFilePath (const std::string& str)
{
  size_t found;
  found=str.find_last_of("/\\");
  return str.substr(0,found);
}
std::string extractFilePathWithEndingSlash(const std::string& str)
{
	size_t found;
	found = str.find_last_of("/\\");
	return str.substr(0, found+1);
}

std::string extractFileName (const std::string& str)
{
  size_t found;
  found=str.find_last_of("/\\");
  return str.substr(found+1);
}


bool file_exists(const char* filename) {
	std::ifstream ifile(filename);
	return ifile.good();
}

int strToInt(const std::string &s) {

	int Result;
	if (!(std::istringstream(s) >> Result)) Result = 3;
	return Result;
}


std::string& str_replace(const std::string &search, const std::string &replace, std::string &subject)
{
	std::string buffer;

	int sealeng = search.length();
	int strleng = subject.length();

	if (sealeng == 0)
		return subject;//no change

	for (int i = 0, j = 0; i<strleng; j = 0)
	{
		while (i + j<strleng && j<sealeng && subject[i + j] == search[j])
			j++;
		if (j == sealeng)//found 'search'
		{
			buffer.append(replace);
			i += sealeng;
		} else
		{
			buffer.append(&subject[i++], 1);
		}
	}
	subject = buffer;
	return subject;
}

int get_png_bitdepth(const char* filename) {

	FILE* f = fopen(filename, "rb");
	char tmp[8];
	fread(tmp, sizeof(char), 8, f);
	unsigned int w, h;
	unsigned int d;
	char b, c;
	fread(&d, sizeof(d), 1, f);
	fread(tmp, sizeof(char), 4, f);
	fread(&w, sizeof(w), 1, f);
	fread(&h, sizeof(h), 1, f);
	fread(&b, sizeof(b), 1, f);
	fread(&c, sizeof(c), 1, f);
	fclose(f);
	int bb = b;
	return bb;
}
template bool load_image(const char* filename, std::vector<double> &val, size_t &W, size_t &H, bool bits16);
template bool load_image(const char* filename, std::vector<float> &val, size_t &W, size_t &H, bool bits16);
template bool load_image(const char* filename, std::vector<int> &val, size_t &W, size_t &H, bool bits16);
template bool load_image(const char* filename, std::vector<unsigned char> &val, size_t &W, size_t &H, bool bits16);

template<typename T>
bool load_image(const char* filename, std::vector<T> &val, size_t &W, size_t &H, bool bits16) {

	if (!file_exists(filename)) return false;

	unsigned char* data;
	int x, y, n;
	if (data = stbi_load(filename, &x, &y, &n, 3)) {
		W = x;
		H = y;
		val.resize(W*H * 3);
		for (int i = 0; i < W*H * 3; i++) {
			val[i] = data[i];
		}
		for (int i = 0; i < H / 2; i++) {
			for (int j = 0; j < W; j++) {
				std::swap(val[((H - i - 1)*W + j) * 3 + 0], val[(i*W + j) * 3 + 0]);
				std::swap(val[((H - i - 1)*W + j) * 3 + 1], val[(i*W + j) * 3 + 1]);
				std::swap(val[((H - i - 1)*W + j) * 3 + 2], val[(i*W + j) * 3 + 2]);
			}
		}
		stbi_image_free(data);
		return true;
	}

	if (bits16 || (std::string(filename).find(".png") != std::string::npos && get_png_bitdepth(filename) == 16)) {
		cimg_library::CImg<float> cimg(filename);
		W = cimg.width();
		H = cimg.height();
		val.resize(W*H * 3);

		for (int i = 0; i<W*H; i++) {
			val[i * 3] = cimg.data()[i] / 150.;
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
			val[i * 3] = cimg.data()[i];
			val[i * 3 + 1] = cimg.data()[i];
			val[i * 3 + 2] = cimg.data()[i];
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
		val[i * 3] = cimg.data()[i];
		val[i * 3 + 1] = cimg.data()[i + W*H];
		val[i * 3 + 2] = cimg.data()[i + W*H * 2];
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

template void save_image(const char* filename, const double* val, int W, int H, double maxval, int nbchan);
template void save_image(const char* filename, const float* val, int W, int H, float maxval, int nbchan);
template void save_image(const char* filename, const int* val, int W, int H, int maxval, int nbchan);
template void save_image(const char* filename, const unsigned char* val, int W, int H, unsigned char maxval, int nbchan);

template<typename T>
void save_image(const char* filename, const T* val, int W, int H, T maxval, int nbchan) {

	if (nbchan == 3) {

		std::string s(filename), ls(filename);
		std::transform(s.begin(), s.end(), ls.begin(), ::tolower);

		if (ls.find(".hdr") != std::string::npos) {
			std::string ret = EncodeRadianceHDR(W, H, (const float*)val, nbchan);
			FILE* f = fopen(filename, "wb+");
			fwrite(&ret[0], sizeof(ret[0]), ret.size(), f);
			fclose(f);
			return;
		}
		std::vector<unsigned char> scaled(W*H * 3); 
		for (int i = 0; i<3*W*H; i++) {
			scaled[i] = std::min(255., std::max(0., val[i] * (255. / maxval)));
		}

		// preferably use STB
		if (ls.find(".bmp") != std::string::npos) {
			stbi_write_bmp(filename, W, H, 3, &scaled[0]);
			return;
		}
		if (ls.find(".tga") != std::string::npos) {
			stbi_write_tga(filename, W, H, 3, &scaled[0]);
			return;
		}
		if (ls.find(".jpg") != std::string::npos) {
			stbi_write_jpg(filename, W, H, 3, &scaled[0], 100);
			return;
		}
		if (ls.find(".png") != std::string::npos) {
			stbi_write_png(filename, W, H, 3, &scaled[0], 0);
			return;
		}

		std::vector<unsigned char> deinterleaved(W*H * 3);
		for (int i = 0; i<W*H; i++) {
			deinterleaved[i] = std::min(255., std::max(0., val[i * 3] * (255. / maxval)));
			deinterleaved[i + W*H] = std::min(255., std::max(0., val[i * 3 + 1] * (255. / maxval)));
			deinterleaved[i + 2 * W*H] = std::min(255., std::max(0., val[i * 3 + 2] * (255. / maxval)));
		}

		cimg_library::CImg<unsigned char> cimg(&deinterleaved[0], W, H, 1, 3);
		cimg.save(filename, 0);
	} else {

		std::vector<unsigned char> deinterleaved(W*H * 3);
		for (int i = 0; i<W*H; i++) {
			deinterleaved[i] = std::min(255., std::max(0., val[i] * (255. / maxval)));
			deinterleaved[i + W*H] = deinterleaved[i];
			deinterleaved[i + 2 * W*H] = deinterleaved[i];
		}

		cimg_library::CImg<unsigned char> cimg(&deinterleaved[0], W, H, 1, 3);
		cimg.save(filename, 0);
	}
}
