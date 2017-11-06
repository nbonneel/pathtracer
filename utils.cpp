#include "utils.h"
#include <fstream>
#include <sstream>

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
	int cc = c;
	return bb;
}