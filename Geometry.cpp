#define _CRT_SECURE_NO_WARNINGS
#include "Geometry.h"
#include <list>
#include <set>
#include <string> 
#include <algorithm>
#include "utils.h"
#include "TriangleMesh.h"


Object* Object::create_from_file(FILE* f) {

	char line[255];
	fscanf(f, "%[^\n]\n", line);
	if (line[4] == 'M') { // mesh
		return Geometry::create_from_file(f);
	}
	if (line[4] == 'S') { // sphere
		return Sphere::create_from_file(f);
	}
	if (line[4] == 'P') { // Plane
		return Plane::create_from_file(f);
	}
}


void load_bmp(const char* filename, std::vector<unsigned char> &tex, int &W, int& H) {
	size_t w, h;
	load_image(filename, tex, w, h);
	W = w;
	H = h;
	
	/*FILE* f;
	f = fopen(filename, "rb");
	unsigned char info[54];
	fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

	W = *(int*)&info[18]; // extract image height and width from header
	H = *(int*)&info[22];

	int size = 3 * W * H;
	tex.resize(size); // allocate 3 bytes per pixel
	fread(&tex[0], sizeof(unsigned char), size, f); // read the rest of the data at once
	fclose(f);

	for (int i = 0; i < size; i += 3) {
		std::swap(tex[i], tex[i + 2]);
	}*/
}


void Object::add_normalmap(const char* filename) {
	normal_map.push_back(Texture(filename, 2, Vector(0., 0., 1.)));
}

void Object::set_normalmap(const char* filename, int idx) {
	if (idx >= normal_map.size()) return;
	//normal_map[idx] = Texture(filename, 2, Vector(0., 0., 1.));
	normal_map[idx].loadNormals(filename);
}

void Object::add_alphamap(const char* filename) {
	alphamap.push_back(Texture(filename, 3, Vector(1., 1., 1.)));
}

void Object::set_texture(const char* filename, int idx) {
	if (idx >= textures.size()) return;
	//textures[idx] = Texture(filename, 0, Vector(1., 1., 1.));
	textures[idx].loadColors(filename);
}

void Object::add_texture(const char* filename) {
	textures.push_back(Texture(filename, 0, Vector(1., 1., 1.)));
}
void Object::add_specularmap(const char* filename) {
	specularmap.push_back(Texture(filename, 1, Vector(0., 0., 0.)));
}
void Object::set_specularmap(const char* filename, int idx) {
	if (idx >= specularmap.size()) return;
	//specularmap[idx] = Texture(filename, 1, Vector(0., 0., 0.));
	specularmap[idx].loadColors(filename);
}


void Object::add_roughnessmap(const char* filename) {
	roughnessmap.push_back(Texture(filename, 4, Vector(1., 1., 1.)));
}
void Object::set_roughnessmap(const char* filename, int idx) {
	if (idx >= roughnessmap.size()) return;
	roughnessmap[idx] = Texture(filename, 4, Vector(1., 1., 1.));
}

void Object::set_alphamap(const char* filename, int idx) {
	if (idx >= alphamap.size()) return;
	alphamap[idx] = Texture(filename, 3, Vector(1., 1., 1.));
}
void Object::set_col_alpha(double col, int idx) {
	if (idx >= alphamap.size()) return;
	alphamap[idx] = Texture("Null", 3, Vector(col, col, col));
}
void Object::add_col_specular(const Vector& col) {
//	std::string name = std::string("Color: (") + std::to_string((int)(col[0] * 255)) + std::string(", ") + std::to_string((int)(col[1] * 255)) + std::string(", ") + std::to_string((int)(col[2] * 255)) + std::string(")");
	specularmap.push_back(Texture("Null", 1, col));
}

void Object::add_col_roughness(const Vector& col) {	
	//std::string name = std::string("Color: (") + std::to_string((int)(col[0])) + std::string(", ") + std::to_string((int)(col[1])) + std::string(", ") + std::to_string((int)(col[2])) + std::string(")");
	roughnessmap.push_back(Texture("Null", 4, col));
}

void Object::set_col_roughness(const Vector& col, int idx) {
	if (idx >= roughnessmap.size()) return;

	//std::string name = std::string(std::string("Color: (") + std::to_string((int)(col[0])) + std::string(", ") + std::to_string((int)(col[1])) + std::string(", ") + std::to_string((int)(col[2])) + std::string(")"));
	//roughnessmap[idx] = Texture(name.c_str(), 4, col);
	roughnessmap[idx].multiplier = col;
}

void Object::set_col_specular(const Vector& col, int idx) {
	if (idx >= specularmap.size()) return;

	//std::string name = std::string("Color: (") + std::to_string((int)(col[0] * 255)) + std::string(", ") + std::to_string((int)(col[1] * 255)) + std::string(", ") + std::to_string((int)(col[2] * 255)) + std::string(")");
	specularmap[idx].multiplier = col;// = Texture(name.c_str(), 1, col);
}



void Object::remove_texture(int id) {
	textures.erase(textures.begin() + id);
}
void Object::remove_alpha(int id) {
	alphamap.erase(alphamap.begin() + id);
}
void Object::remove_specular(int id) {
	specularmap.erase(specularmap.begin() + id);
}
void Object::remove_normal(int id) {
	normal_map.erase(normal_map.begin() + id);
}
void Object::remove_roughness(int id) {
	roughnessmap.erase(roughnessmap.begin() + id);
}
void Object::swap_roughness(int id1, int id2) {
	std::swap(roughnessmap[id1], roughnessmap[id2]);
}
void Object::swap_textures(int id1, int id2) {
	std::swap(textures[id1], textures[id2]);
}
void Object::swap_normal(int id1, int id2) {
	std::swap(normal_map[id1], normal_map[id2]);
}
void Object::swap_specular(int id1, int id2) {
	std::swap(specularmap[id1], specularmap[id2]);
}
void Object::swap_alpha(int id1, int id2) {
	std::swap(alphamap[id1], alphamap[id2]);
}
void Object::add_col_texture(const Vector& col) {
		
	//std::string name = (std::string("Color: (") + std::to_string((int)(col[0] * 255)) + std::string(", ") + std::to_string((int)(col[1] * 255)) + std::string(", ") + std::to_string((int)(col[2] * 255)) + std::string(")"));
	textures.push_back(Texture("Null", 0, col));
}

void Object::set_col_texture(const Vector& col, int idx) {
	if (idx >= textures.size()) return;

	//std::string name = std::string("Color: (") + std::to_string((int)(col[0] * 255)) + std::string(", ") + std::to_string((int)(col[1] * 255)) + std::string(", ") + std::to_string((int)(col[2] * 255)) + std::string(")");
	textures[idx].multiplier = col;// = Texture("Null", 0, col);
}

void Object::add_null_normalmap() {
	normal_map.push_back(Texture::defaultNormal());
}
void Object::set_null_normalmap(int idx) {
	if (idx >= textures.size()) return;
	normal_map[idx].clear_texture();
}
void Object::add_col_alpha(double col) {
	alphamap.push_back(Texture("Null", 3, Vector(col, col, col)));
}


