#pragma once
#include <string>

std::string EncodeFloatTIFF(unsigned int w, unsigned int h, float* RGBdata, unsigned int floatsPerPixel = 4);
std::string EncodeRadianceHDR(unsigned int w, unsigned int h, const float* RGBdata, unsigned int floatsPerPixel = 4);