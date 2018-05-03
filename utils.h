#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <string>

#include <vector>



bool file_exists(const char* filename);
std::string& str_replace(const std::string &search, const std::string &replace, std::string &subject);
std::string extractFilePath (const std::string& str);
std::string extractFileName (const std::string& str);
std::string extractFilePathWithEndingSlash(const std::string& str);

int get_png_bitdepth(const char* filename);

template<typename T> void save_image(const char* filename, const T* val, int W, int H, T maxval = 255, int nbchan = 3);
template<typename T> bool load_image(const char* filename, std::vector<T> &val, size_t &W, size_t &H, bool bits16 = false);