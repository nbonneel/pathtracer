#include "TriangleMesh.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream> 
#include <fstream>
#include <string>
#include <set>
#include <list>

void Geometry::readVRML(const char* obj) {
	bool shuffle_colors = true;

	FILE* f;
	f = fopen(obj, "r");
	int curGroup = -1;
	bool readpoints = false, readcolors = false, readfaces = false;

	int max_line_size = 50000000;
	char* line = new char[max_line_size];
	while (!feof(f)) {
		if (!fgets(line, max_line_size, f)) break;
		std::string s(line);

		if (s.find("point") != std::string::npos)
			readpoints = true;
		if (s.find("color") != std::string::npos)
			readcolors = true;
		if (s.find("coordIndex") != std::string::npos)
			readfaces = true;

		float f1, f2, f3;
		int offset;
		char* consumedline = line;
		if (readpoints) {
			int nread;
			while (true) {
				if (consumedline[0] == '\n') break;
				nread = sscanf(consumedline, "%f %f %f%n", &f1, &f2, &f3, &offset);
				if (nread == 3) {
					Vector vec(f1, f2, f3);
					vertices.push_back(vec);
					consumedline = consumedline + offset;
				} else {
					consumedline = consumedline + 1;
				}
			}
		}

		if (readcolors) {
			int nread;
			while (true) {
				if (consumedline[0] == '\n' || consumedline[0] == '\0') break;
				nread = sscanf(consumedline, "%f %f %f%n", &f3, &f2, &f1, &offset);
				if (nread == 3) {
					Vector vec(f3, f2, f1);
					if (shuffle_colors) {
						Vector shuff = vec*vec*255. * 987 + Vector(1234, 4567, 7891);
						//Vector shuff = vec*Vector(vec[2], vec[0], vec[1])*255. * 9871 + Vector(1234, 14567, 17891);
						//Vector shuff = vec*Vector(vec[2], vec[0], vec[1])*256. * 9871 + Vector(1234, 14567, 17891);						
						//Vector shuff = vec*Vector(vec[2], vec[0], vec[1])*258. * 9871 + Vector(1234, 14567, 17891);
						//Vector shuff = vec*255.;
						vec = Vector(((int)(shuff[0]) % 255) / 255., ((int)(shuff[1]) % 255) / 255., ((int)(shuff[2]) % 255) / 255.);
					}

					facecolors.push_back(vec);
					consumedline = consumedline + offset;
				} else {
					consumedline = consumedline + 1;
				}
			}
		}

		if (readfaces) {
			int nread;
			while (true) {
				if (consumedline[0] == '\n') break;
				int i1, i2, i3;
				nread = sscanf(consumedline, "%d,%d,%d,-1%n", &i1, &i2, &i3, &offset);
				if (nread == 3) {
					/*faces.push_back(i1);
					faces.push_back(i2);
					faces.push_back(i3);*/
					TriangleIndices t;
					t.vtxi = i1;
					t.vtxj = i2;
					t.vtxk = i3;
					indices.push_back(t);
					consumedline = consumedline + offset;
				} else {
					consumedline = consumedline + 1;
				}
			}
		}

		if (s.find("]") != std::string::npos) {
			readpoints = false;
			readcolors = false;
			readfaces = false;
		}

	}
	fclose(f);

	delete[] line;
}

void Geometry::load_edge_colors(const char* csvfilename) {
	if (!csvfilename) return;

	edgecolor.resize(vertices.size());
	//vertexcolors.resize(vertices.size());

	std::map<Edge, std::vector<int> > edges_to_faces;
	for (int i = 0; i < indices.size(); i++) {
		int idv1 = std::min(indices[i].vtxi, indices[i].vtxj);
		int idv2 = std::max(indices[i].vtxi, indices[i].vtxj);
		edges_to_faces[Edge(idv1, idv2)].push_back(i);

		idv1 = std::min(indices[i].vtxk, indices[i].vtxj);
		idv2 = std::max(indices[i].vtxk, indices[i].vtxj);
		edges_to_faces[Edge(idv1, idv2)].push_back(i);

		idv1 = std::min(indices[i].vtxi, indices[i].vtxk);
		idv2 = std::max(indices[i].vtxi, indices[i].vtxk);
		edges_to_faces[Edge(idv1, idv2)].push_back(i);
	}

	std::map<Edge, Edge > faces_to_edges;
	for (auto it = edges_to_faces.begin(); it != edges_to_faces.end(); ++it) {
		if (it->second.size() == 2) {
			int idf1 = std::min(it->second[0], it->second[1]);
			int idf2 = std::max(it->second[0], it->second[1]);
			int idv1 = it->first.a;
			int idv2 = it->first.b;
			faces_to_edges[Edge(idf1, idf2)] = Edge(idv1, idv2);
		}
	}


	FILE * f = fopen(csvfilename, "r+");
	int nbEdges = 0;
	int nbNodes = 0;
	std::vector<std::pair<int, int> > edges;
	std::vector<float> values;
	while (!feof(f)) {
		int cut;
		float val0, val1;
		int idNode0;
		float n0x, n0y, n0z;
		int idNode1;
		float n1x, n1y, n1z;
		char line[255];
		if (!fgets(line, 255, f)) break;

		int nsc = sscanf(line, "%u %f %f %u %f %f %f %u %f %f %f\n", &cut, &val0, &val1, &idNode0, &n0x, &n0y, &n0z, &idNode1, &n1x, &n1y, &n1z);
		if (nsc != 11) continue;

		nbEdges++;
		nbNodes = std::max(nbNodes, idNode0);
		nbNodes = std::max(nbNodes, idNode1);

		int idFace1 = std::min(idNode0, idNode1);
		int idFace2 = std::max(idNode0, idNode1);

		float n0 = sqrt(n0x*n0x + n0y*n0y + n0z*n0z);
		float n1 = sqrt(n1x*n1x + n1y*n1y + n1z*n1z);
		n0x /= n0; n0y /= n0; n0z /= n0;
		n1x /= n1; n1y /= n1; n1z /= n1;
		double scalaire = n0x*n1x + n0y*n1y + n0z*n1z;

		//scalaire = dot(&mesh.normals[idNode0], &mesh.normals[idNode1]);


		//float v = (std::min(1.f, std::max(0.f, val0)) + std::min(1.f, std::max(0.f, val1)))*0.5;		
		val0 = std::min(1.f, std::max(0.f, val0));
		val1 = std::min(1.f, std::max(0.f, val1));
		float v = (val0 + val1)*0.5;// sqrt(val0*val1);

		int idVertex1 = faces_to_edges[Edge(idFace1, idFace2)].a;
		int idVertex2 = faces_to_edges[Edge(idFace1, idFace2)].b;
		edgecolor[idVertex1][idVertex2] = Vector(1., 1, 1)*v + (1. - v)*Vector(1., 0., 0.);

		//vertexcolors[idVertex1] = Vector(1., 1, 1)*v + (1. - v)*Vector(1., 0., 0.);
		//vertexcolors[idVertex2] = Vector(1., 1, 1)*v + (1. - v)*Vector(1., 0., 0.);

		if (scalaire > cos(5 * 3.1416 / 180)) v = 1;
		//v = std::max(0., scalaire); v = sqrt(v);
		values.push_back(std::max(0.001, 1. - v));
	}
	fclose(f);
}

Vector TransformH(
	const Vector &in,  // color to transform
	float H
)
{
	float U = cos(H*M_PI / 180);
	float W = sin(H*M_PI / 180);

	Vector ret;
	ret[0] = (.299 + .701*U + .168*W)*in[0]
		+ (.587 - .587*U + .330*W)*in[1]
		+ (.114 - .114*U - .497*W)*in[2];
	ret[1] = (.299 - .299*U - .328*W)*in[0]
		+ (.587 + .413*U + .035*W)*in[1]
		+ (.114 - .114*U + .292*W)*in[2];
	ret[2] = (.299 - .3*U + 1.25*W)*in[0]
		+ (.587 - .588*U - 1.05*W)*in[1]
		+ (.114 + .886*U - .203*W)*in[2];
	return ret;
}


void Geometry::readOBJ(const char* obj, bool load_textures) {

	char matfile[255];
	char grp[255];

	std::set<int> marked_vertices;

	FILE* f;
	f = fopen(obj, "r");
	int curGroup = -1;
	while (!feof(f)) {
		char line[255];
		if (!fgets(line, 255, f)) break;

		std::string linetrim(line);
		linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
		strcpy(line, linetrim.c_str());

		if (line[0] == 'u' && line[1] == 's') {
			sscanf(line, "usemtl %[^\n]\n", grp);

			if (groupNames.find(std::string(grp)) != groupNames.end()) {
				curGroup = groupNames[std::string(grp)];
			} else {
				curGroup = groupNames.size();
				groupNames[std::string(grp)] = curGroup;
			}
		}
		if (line[0] == 'm' && line[1] == 't' && line[2] == 'l') {
			sscanf(line, "mtllib %[^\n]\n", matfile);
		}
		if (line[0] == 'v' && line[1] == ' ') {
			Vector vec;
			//sscanf(line, "v %lf %lf %lf\n", &vec[2], &vec[0], &vec[1]); // car
			//sscanf_s(line, "v %lf %lf %lf\n", &vec[0], &vec[2], &vec[1]); // girl
			//sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);  // dragon, car2

			Vector col;
			if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
				col[0] = std::min(1., std::max(0., col[0]));
				col[1] = std::min(1., std::max(0., col[1]));
				col[2] = std::min(1., std::max(0., col[2]));

				//if (col[1] > 0.6 && col[0]<0.3) { marked_vertices.insert(vertices.size()); }

				vertices.push_back(vec);
				vertexcolors.push_back(col);

				/*col[1] = 1. - col[1];
				col[2] = 1. - col[2];*/
				//col[1] = 1. ;
				//col[2] = 1.;
				//vertexcolors.push_back(TransformH(col,130));
			} else {
				sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);  // helmet
																			 //vec[2] = -vec[2]; //car2
				vertices.push_back(vec);
			}
		}
		if (line[0] == 'v' && line[1] == 'n') {
			Vector vec;
			//sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]); //dragon, car2
			sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]); // helmet
																		 //sscanf(line, "vn %lf %lf %lf\n", &vec[2], &vec[0], &vec[1]); //car
																		 //sscanf_s(line, "vn %lf %lf %lf\n", &vec[0], &vec[2], &vec[1]); //girl
																		 //vec[2] = -vec[2];  //car2
			normals.push_back(vec);
		}
		if (line[0] == 'v' && line[1] == 't') {
			Vector vec;
			sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
			uvs.push_back(vec);
		}
		if (line[0] == 'f') {
			TriangleIndices t;
			int i0, i1, i2, i3;
			int j0, j1, j2, j3;
			int k0, k1, k2, k3;
			int nn;
			//faceGroup.push_back(curGroup);
			t.group = curGroup;
			t.showEdges[0] = true;
			t.showEdges[1] = true;

			char* consumedline = line + 1;
			int offset;

			nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
			if (nn == 9) {
				if (marked_vertices.find(i0 - 1) != marked_vertices.end()) continue;
				if (marked_vertices.find(i1 - 1) != marked_vertices.end()) continue;
				if (marked_vertices.find(i2 - 1) != marked_vertices.end()) continue;
				if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
				if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
				if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
				if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
				if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
				if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
				if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
				if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
				if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
				t.showEdges[2] = (consumedline[offset] == '\n') || ((consumedline[offset] == ' ' && consumedline[offset + 1] == '\n')); // ugly hack cause ppl put spaces before \n sometimes.
				indices.push_back(t);
			} else {
				nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
				if (nn == 6) {
					if (marked_vertices.find(i0 - 1) != marked_vertices.end()) continue;
					if (marked_vertices.find(i1 - 1) != marked_vertices.end()) continue;
					if (marked_vertices.find(i2 - 1) != marked_vertices.end()) continue;
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					t.showEdges[2] = (consumedline[offset] == '\n') || ((consumedline[offset] == ' ' && consumedline[offset + 1] == '\n'));
					indices.push_back(t);
				} else {
					nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
					if (nn == 3) {
						if (marked_vertices.find(i0 - 1) != marked_vertices.end()) continue;
						if (marked_vertices.find(i1 - 1) != marked_vertices.end()) continue;
						if (marked_vertices.find(i2 - 1) != marked_vertices.end()) continue;
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						t.showEdges[2] = (consumedline[offset] == '\n') || ((consumedline[offset] == ' ' && consumedline[offset + 1] == '\n'));
						indices.push_back(t);
					} else {
						nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
						if (marked_vertices.find(i0 - 1) != marked_vertices.end()) continue;
						if (marked_vertices.find(i1 - 1) != marked_vertices.end()) continue;
						if (marked_vertices.find(i2 - 1) != marked_vertices.end()) continue;
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
						if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
						if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
						t.showEdges[2] = (consumedline[offset] == '\n') || ((consumedline[offset] == ' ' && consumedline[offset + 1] == '\n'));
						indices.push_back(t);
					}
				}
			}
			//}				


			consumedline = consumedline + offset;

			while (true) {
				if (consumedline[0] == '\n') break;
				if (consumedline[0] == '\0') break;
				nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
				TriangleIndices t2;
				t2.group = curGroup;
				t2.showEdges[0] = false;
				t2.showEdges[1] = true;
				if (nn == 3) {
					if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
					if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
					if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
					if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
					if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
					if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
					if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
					if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
					if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
					t2.showEdges[2] = (consumedline[offset] == '\n') || ((consumedline[offset] == ' ' && consumedline[offset+1] == '\n'));
					indices.push_back(t2);
					consumedline = consumedline + offset;
					i2 = i3;
					j2 = j3;
					k2 = k3;
				} else {
					nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
					if (nn == 2) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						t2.showEdges[2] = (consumedline[offset] == '\n') || ((consumedline[offset] == ' ' && consumedline[offset + 1] == '\n'));
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						indices.push_back(t2);
					} else {
						nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
							if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
							if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
							t2.showEdges[2] = (consumedline[offset] == '\n') || ((consumedline[offset] == ' ' && consumedline[offset + 1] == '\n'));
							consumedline = consumedline + offset;
							i2 = i3;
							k2 = k3;
							indices.push_back(t2);
						} else {
							nn = sscanf(consumedline, "%u%n", &i3, &offset);
							if (nn == 1) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								t2.showEdges[2] = (consumedline[offset] == '\n') || ((consumedline[offset] == ' ' && consumedline[offset + 1] == '\n'));
								consumedline = consumedline + offset;
								i2 = i3;
								indices.push_back(t2);
							} else {
								consumedline = consumedline + 1;
							}
						}
					}
				}
			}





		}
		//}

	}
	fclose(f);

	if (groupNames.size() == 0) {
		for (int i = 0; i < indices.size(); i++) {
			indices[i].group = 0;
		}
		groupNames["Default"] = 0;
	}

	//#ifndef _DEBUG
	if (load_textures) {


		for (int i = 0; i < groupNames.size(); i++) {
			add_col_texture(Vector(0.5, 0.5, 0.5));
			add_col_specular(Vector(0., 0., 0.));
			add_col_roughness(Vector(0., 0., 0.));
			add_null_normalmap();
			add_col_alpha(1.);
		}

		std::string filenamemat = extractFilePathWithEndingSlash(std::string(obj)) + std::string(matfile);
		f = fopen(filenamemat.c_str(), "r");
		int illum = -1;
		if (f) {
			while (!feof(f)) {
				char line[255];
				if (!fgets(line, 255, f)) break;
				if (line[0] == 'n' && line[1] == 'e' && line[2] == 'w') {
					sscanf(line, "newmtl %[^\n]\n", grp);
					illum = -1;
				}

				if (line[0] == 'm' && line[4] == 'K' && line[5] == 'd') {
					char texturefile[255];
					sscanf(line, "map_Kd %[^\n]\n", texturefile);
					//add_texture((std::string("BeautifulGirl/Texture/") + std::string(texturefile)).c_str());
					std::string filenamtext = extractFilePathWithEndingSlash(std::string(obj)) + std::string(texturefile);
					//set_texture(filenamtext.c_str(), groupNames[std::string(grp)]);
					textures[groupNames[std::string(grp)]].loadColors(filenamtext.c_str());
				}

				if (line[0] == 'm' && line[4] == 'K' && line[5] == 's') {
					char texturefile[255];
					sscanf(line, "map_Ks %[^\n]\n", texturefile);
					//add_texture((std::string("BeautifulGirl/Texture/") + std::string(texturefile)).c_str());
					std::string filenamtext = extractFilePathWithEndingSlash(std::string(obj)) + std::string(texturefile);
					//set_specularmap(filenamtext.c_str(), groupNames[std::string(grp)]);
					specularmap[groupNames[std::string(grp)]].loadColors(filenamtext.c_str());
				}
				if (line[0] == 'm' && line[4] == 'B' && line[5] == 'u') {
					char texturefile[255];
					sscanf(line, "map_Bump %[^\n]\n", texturefile);
					//add_texture((std::string("BeautifulGirl/Texture/") + std::string(texturefile)).c_str());
					std::string filenamtext = extractFilePathWithEndingSlash(std::string(obj)) + std::string(texturefile);
					//set_normalmap(filenamtext.c_str(), groupNames[std::string(grp)]);
					normal_map[groupNames[std::string(grp)]].loadNormals(filenamtext.c_str());
				}
				if (line[0] == 'm' && line[1] == 'a' && line[4] == 'd') {
					char texturefile[255];
					sscanf(line, "map_d %[^\n]\n", texturefile);
					std::string filenamtext = extractFilePathWithEndingSlash(std::string(obj)) + std::string(texturefile);
					//set_alphamap(filenamtext.c_str(), groupNames[std::string(grp)]);
					alphamap[groupNames[std::string(grp)]].loadColors(filenamtext.c_str());
				}

				if (line[0] == 'K' && line[1] == 'd') {
					Vector Kd;
					sscanf(line, "Kd %lf %lf %lf\n", &Kd[0], &Kd[1], &Kd[2]);
					//set_col_texture(Kd, groupNames[std::string(grp)]);
					textures[groupNames[std::string(grp)]].multiplier = Kd;
				}
				if (line[0] == 'K' && line[1] == 's') {
					Vector Ks;
					sscanf(line, "Ks %lf %lf %lf\n", &Ks[0], &Ks[1], &Ks[2]);
					//set_col_specular(Ks, groupNames[std::string(grp)]);
					specularmap[groupNames[std::string(grp)]].multiplier = Ks;
					if (illum == 0 || illum == 1) specularmap[groupNames[std::string(grp)]].multiplier = Vector(0., 0., 0.);
				}
				if (line[0] == 'N' && line[1] == 's') {
					Vector Ns;
					int ret = sscanf(line, "Ns %lf %lf %lf\n", &Ns[0], &Ns[1], &Ns[2]);
					if (ret == 1) Ns = Vector(Ns[0], Ns[0], Ns[0]);
					//set_col_roughness(Ns, groupNames[std::string(grp)]);
					roughnessmap[groupNames[std::string(grp)]].multiplier = Ns;
				}
				if (line[0] == 'i' && line[1] == 'l' && line[3] == 'u') {
					int ret = sscanf(line, "illum %u\n", &illum);
					if (illum == 0 || illum == 1) specularmap[groupNames[std::string(grp)]].multiplier = Vector(0., 0., 0.);

				}

			}
			fclose(f);
		}
	}
	//#endif
	//uvs.clear();
	//uvIds.clear();
}

void Geometry::exportMTL(const char* mtlfile) {

	FILE* f = fopen(mtlfile, "w+");

	for (std::map<std::string, int>::iterator it = groupNames.begin(); it != groupNames.end(); ++it) {

		std::string gname = it->first;
		int gid = it->second;

		fprintf(f, "newmtl %s\n", gname.c_str());
		fprintf(f, "Kd %f %f %f\n", textures[gid].multiplier[0], textures[gid].multiplier[1], textures[gid].multiplier[2]);
		if (textures[gid].W>0)
			fprintf(f, "map_Kd %s\n", textures[gid].filename.c_str());
		fprintf(f, "Ks %f %f %f\n", specularmap[gid].multiplier[0], specularmap[gid].multiplier[1], specularmap[gid].multiplier[2]);
		if (specularmap[gid].W>0)
			fprintf(f, "map_Ks %s\n", specularmap[gid].filename.c_str());
		fprintf(f, "Ns %f\n", roughnessmap[gid].multiplier[0]);
		if (roughnessmap[gid].W>0)
			fprintf(f, "map_Ns %s\n", roughnessmap[gid].filename.c_str());
		if (alphamap[gid].W>0)
			fprintf(f, "map_d %s\n", alphamap[gid].filename.c_str());
		if (normal_map[gid].W>0)
			fprintf(f, "map_bump %s\n", normal_map[gid].filename.c_str());
	}

	fclose(f);

}


void Geometry::setup_tangents() {

	std::vector<Vector> tan1(vertices.size(), Vector(0, 0, 0));
	std::vector<Vector> tan2(vertices.size(), Vector(0, 0, 0));
	for (int i = 0; i < indices.size(); i++) {
		int a = indices[i].vtxi;
		int b = indices[i].vtxj;
		int c = indices[i].vtxk;
		if (indices[i].uvi == -1) continue;
		if (indices[i].uvj == -1) continue;
		if (indices[i].uvk == -1) continue;

		Vector vA = vertices[b] - vertices[a];
		Vector vB = vertices[c] - vertices[a];
		Vector sA = uvs[indices[i].uvj] - uvs[indices[i].uvi];
		Vector sB = uvs[indices[i].uvk] - uvs[indices[i].uvi];
		double det = (sA[0] * sB[1] - sB[0] * sA[1]);
		Vector sdir, tdir;
		if (det != 0) {
			sdir = (sB[1] * vA - sA[1] * vB) / det; // tan 1
			tdir = (sA[0] * vB - sB[0] * vA) / det; // tan 2
		} else {
			sdir = vA*0.00001;
			tdir = vB*0.00001;
		}

		tan1[a] += sdir; tan1[b] += sdir; tan1[c] += sdir;
		tan2[a] += tdir; tan2[b] += tdir; tan2[c] += tdir;
	}
	/*std::vector<Vector> computedVtxNormals(vertices.size(), Vector(0.,0.,0.));
	for (int i = 0; i < indices.size(); i++) {
	int a = indices[i].vtxi;
	int b = indices[i].vtxj;
	int c = indices[i].vtxk;
	indices[i].ni = a;
	indices[i].nj = b;
	indices[i].nk = c;
	Vector vA = vertices[b] - vertices[a];
	Vector vB = vertices[c] - vertices[a];
	Vector faceN = cross(vA, vB).getNormalized();
	computedVtxNormals[a] += faceN; computedVtxNormals[b] += faceN; computedVtxNormals[c] += faceN;
	}
	for (int i = 0; i < computedVtxNormals.size(); i++) {
	double n = computedVtxNormals[i].getNorm2();
	if (n > 1E-12) {
	computedVtxNormals[i] /= sqrt(n);
	} else {
	computedVtxNormals[i] = Vector(0, 1, 0);
	}
	}*/
	// no vtx normals = per face normal (sharp edges)
	for (int i = 0; i < indices.size(); i++) {
		int a = indices[i].vtxi;
		int b = indices[i].vtxj;
		int c = indices[i].vtxk;
		Vector vA = vertices[b] - vertices[a];
		Vector vB = vertices[c] - vertices[a];
		Vector faceN = cross(vA, vB).getNormalized();
		if (indices[i].ni == -1) {
			normals.push_back(faceN);
			indices[i].ni = normals.size() - 1;
			if (indices[i].nj == -1) indices[i].nj = normals.size() - 1;
			if (indices[i].nk == -1) indices[i].nk = normals.size() - 1;
		}
		if (indices[i].nj == -1) {
			normals.push_back(faceN);
			indices[i].nj = normals.size() - 1;
			if (indices[i].nk == -1) indices[i].nk = normals.size() - 1;
		}
		if (indices[i].nk == -1) {
			normals.push_back(faceN);
			indices[i].nk = normals.size() - 1;
		}
	}

	std::vector<int> verticesToNormals(vertices.size());
	for (int i = 0; i < indices.size(); i++) {
		verticesToNormals[indices[i].vtxi] = indices[i].ni;
		verticesToNormals[indices[i].vtxj] = indices[i].nj;
		verticesToNormals[indices[i].vtxk] = indices[i].nk;
	}
	/*	for (int i = 0; i < vertices.size(); i++) {
	if (verticesToNormals[i] == -1) {
	normals.push_back(computedVtxNormals[i]);
	verticesToNormals[i] = normals.size() - 1;
	}
	}*/


	tangents.resize(vertices.size());
	bitangents.resize(vertices.size());
	for (int i = 0; i < vertices.size(); i++) {
		Vector N = normals[verticesToNormals[i]].getNormalized();
		tangents[i] = (tan1[i] - N*dot(tan1[i], N)).getNormalized();
		double wtangent = (dot(cross(N, tan1[i]), tan2[i]) < 0) ? -1.0 : 1.0;
		bitangents[i] = cross(N, tangents[i])*wtangent;
	}

}


Geometry::Geometry(const char* obj, double scaling, const Vector& offset, bool mirror, bool transp, const char* colors_csv_filename) {
	init(obj, scaling, offset, mirror, transp, colors_csv_filename, true);
}

void Geometry::init(const char* obj, double scaling, const Vector& offset, bool mirror, bool transp, const char* colors_csv_filename, bool load_textures) {
	display_edges = false;
	interp_normals = true;
	miroir = mirror;
	transparent = transp;
	name = obj;

	std::string filename(obj);
	std::string lowerFilename = filename;
	std::transform(filename.begin(), filename.end(), lowerFilename.begin(), ::tolower);

	if (lowerFilename.find(".wrl") != std::string::npos) {
		readVRML(obj);
	} else {
		if (lowerFilename.find(".obj") != std::string::npos) {
			readOBJ(obj, load_textures);
		}
	}

	for (int i = 0; i < vertices.size(); i++) {
		std::swap(vertices[i][0], vertices[i][2]);
		vertices[i][0] = -vertices[i][0];
	}
	for (int i = 0; i < normals.size(); i++) {
		std::swap(normals[i][0], normals[i][2]);
		normals[i][0] = -normals[i][0];
	}

	BBox bb(Vector(1E9, 1E9, 1E9), Vector(-1E9, -1E9, -1E9));
	for (int i = 0; i < vertices.size(); i++) {
		bb.bounds[0][0] = std::min(bb.bounds[0][0], vertices[i][0]); bb.bounds[1][0] = std::max(bb.bounds[1][0], vertices[i][0]);
		bb.bounds[0][1] = std::min(bb.bounds[0][1], vertices[i][1]); bb.bounds[1][1] = std::max(bb.bounds[1][1], vertices[i][1]);
		bb.bounds[0][2] = std::min(bb.bounds[0][2], vertices[i][2]); bb.bounds[1][2] = std::max(bb.bounds[1][2], vertices[i][2]);
	}

	double s = std::max(bb.bounds[1][0] - bb.bounds[0][0], std::max(bb.bounds[1][1] - bb.bounds[0][1], bb.bounds[1][2] - bb.bounds[0][2]));
	Vector c = (bb.bounds[0] + bb.bounds[1])*0.5;
	//rotation_center = Vector(0., 0., 0.);
	for (int i = 0; i < vertices.size(); i++) {
		vertices[i][0] = (vertices[i][0] - c[0]) / s * scaling + offset[0];
		vertices[i][1] = (vertices[i][1] - c[1]) / s * scaling + offset[1];
		vertices[i][2] = (vertices[i][2] - c[2]) / s * scaling + offset[2];
		//rotation_center += vertices[i];
	}
	//rotation_center *= 1. / vertices.size();
	if (colors_csv_filename)
		csv_file = std::string(colors_csv_filename);
	else
		csv_file = std::string("");
	load_edge_colors(colors_csv_filename);

	max_bvh_triangles = 0;
	build_bvh(&bvh, 0, indices.size());
	triangleSoup.resize(indices.size());
	for (int i = 0; i < indices.size(); i++) {
		triangleSoup[i] = Triangle(vertices[indices[i].vtxi], vertices[indices[i].vtxj], vertices[indices[i].vtxk]);
	}

	rotation_center = (bvh.bbox.bounds[0] + bvh.bbox.bounds[1])*0.5;

	setup_tangents();
}

BBox Geometry::build_bbox(int i0, int i1) {

	BBox result;
	result.bounds[1] = vertices[indices[i0].vtxi];
	result.bounds[0] = vertices[indices[i0].vtxi];
	for (int i = i0; i < i1; i++) { // indice de triangle
		for (int k = 0; k < 3; k++) { // indice de dimension
			result.bounds[0][k] = std::min(result.bounds[0][k], vertices[indices[i].vtxi][k]);
			result.bounds[1][k] = std::max(result.bounds[1][k], vertices[indices[i].vtxi][k]);
			result.bounds[0][k] = std::min(result.bounds[0][k], vertices[indices[i].vtxj][k]);
			result.bounds[1][k] = std::max(result.bounds[1][k], vertices[indices[i].vtxj][k]);
			result.bounds[0][k] = std::min(result.bounds[0][k], vertices[indices[i].vtxk][k]);
			result.bounds[1][k] = std::max(result.bounds[1][k], vertices[indices[i].vtxk][k]);
		}
	}
	return result;
}

BBox Geometry::build_centers_bbox(int i0, int i1) {

	BBox result;
	result.bounds[1] = (vertices[indices[i0].vtxi] + vertices[indices[i0].vtxj] + vertices[indices[i0].vtxk]) / 3.;
	result.bounds[0] = (vertices[indices[i0].vtxi] + vertices[indices[i0].vtxj] + vertices[indices[i0].vtxk]) / 3.;
	for (int i = i0; i < i1; i++) { // indice de triangle
		for (int k = 0; k < 3; k++) { // indice de dimension
			Vector center = (vertices[indices[i].vtxi] + vertices[indices[i].vtxj] + vertices[indices[i].vtxk]) / 3.;
			result.bounds[0][k] = std::min(result.bounds[0][k], center[k]);
			result.bounds[1][k] = std::max(result.bounds[1][k], center[k]);
		}
	}
	return result;
}

void Geometry::build_bvh(BVH* b, int i0, int i1) {
	bvh.bbox = build_bbox(i0, i1);
	bvh.nodes.reserve(indices.size() * 2);
	bvh_avg_depth = 0;
	bvh_nb_nodes = 0;
	build_bvh_recur(b, 0, i0, i1, 0);
	bvh_avg_depth /= (double)bvh_nb_nodes;
}

void Geometry::build_bvh_recur(BVH* b, int node, int i0, int i1, int depth) {

	BVHNodes n;
	//n.i0 = i0;
	//n.i1 = i1;
	n.bbox = build_bbox(i0, i1);
	n.fg = i0;
	n.fd = i1;
	n.isleaf = true;
	b->nodes.push_back(n);
	bvh_depth = std::max(bvh_depth, depth);
	bvh_avg_depth += depth;
	bvh_nb_nodes++;

	BBox centerBB = build_centers_bbox(i0, i1);

	Vector diag = centerBB.bounds[1] - centerBB.bounds[0];
	int split_dim;
	if ((diag[0] >= diag[1]) && (diag[0] >= diag[2])) {
		split_dim = 0;
	} else {
		if ((diag[1] >= diag[0]) && (diag[1] >= diag[2])) {
			split_dim = 1;
		} else {
			split_dim = 2;
		}
	}


	double best_split_factor = 0.5;
	double best_area_bb = 1E50;
#ifdef _DEBUG
	int max_tests = 1;
#else
	int max_tests = 16;
#endif

	for (int test_split = 0; test_split < max_tests; test_split++) {

		double cur_split_factor = (test_split + 1) / (double)(max_tests + 1);
		double split_val = centerBB.bounds[0][split_dim] + diag[split_dim] * cur_split_factor;
		BBox bb_left(Vector(1E10, 1E10, 1E10), Vector(-1E10, -1E10, -1E10)), bb_right(Vector(1E10, 1E10, 1E10), Vector(-1E10, -1E10, -1E10));
		int nl = 0, nr = 0;

		for (int i = i0; i < i1; i++) {
			double center_split_dim = (vertices[indices[i].vtxi][split_dim] + vertices[indices[i].vtxj][split_dim] + vertices[indices[i].vtxk][split_dim]) / 3.;

			if (center_split_dim <= split_val) {
				bb_left.bounds[0] = min(bb_left.bounds[0], vertices[indices[i].vtxi]);
				bb_left.bounds[0] = min(bb_left.bounds[0], vertices[indices[i].vtxj]);
				bb_left.bounds[0] = min(bb_left.bounds[0], vertices[indices[i].vtxk]);
				bb_left.bounds[1] = max(bb_left.bounds[1], vertices[indices[i].vtxi]);
				bb_left.bounds[1] = max(bb_left.bounds[1], vertices[indices[i].vtxj]);
				bb_left.bounds[1] = max(bb_left.bounds[1], vertices[indices[i].vtxk]);
				nl++;
			} else {
				bb_right.bounds[0] = min(bb_right.bounds[0], vertices[indices[i].vtxi]);
				bb_right.bounds[0] = min(bb_right.bounds[0], vertices[indices[i].vtxj]);
				bb_right.bounds[0] = min(bb_right.bounds[0], vertices[indices[i].vtxk]);
				bb_right.bounds[1] = max(bb_right.bounds[1], vertices[indices[i].vtxi]);
				bb_right.bounds[1] = max(bb_right.bounds[1], vertices[indices[i].vtxj]);
				bb_right.bounds[1] = max(bb_right.bounds[1], vertices[indices[i].vtxk]);
				nr++;
			}
		}
		double sum_area_bb = bb_left.area()*nl + bb_right.area()*nr;
		if (sum_area_bb < best_area_bb) {
			best_split_factor = cur_split_factor;
			best_area_bb = sum_area_bb;
		}
	}

	double split_val = centerBB.bounds[0][split_dim] + diag[split_dim] * best_split_factor;
	int pivot = i0 - 1;
	for (int i = i0; i < i1; i++) {
		double center_split_dim = (vertices[indices[i].vtxi][split_dim] + vertices[indices[i].vtxj][split_dim] + vertices[indices[i].vtxk][split_dim]) / 3.;

		if (center_split_dim <= split_val) {
			pivot++;
			std::swap(indices[i], indices[pivot]);

			if ((facecolors.size() != 0) && (i < facecolors.size()) && (pivot < facecolors.size())) {
				std::swap(facecolors[i], facecolors[pivot]);
			}
		}
	}


	if (pivot < i0 || pivot >= i1 - 1 || i1 <= i0 + 4) { // 1 triangles per leaf
		max_bvh_triangles = std::max(max_bvh_triangles, i1 - i0);
		return;
	}

	b->nodes[node].isleaf = false;
	b->nodes[node].fg = b->nodes.size();
	build_bvh_recur(b, b->nodes[node].fg, i0, pivot + 1, depth + 1);

	b->nodes[node].fd = b->nodes.size();
	build_bvh_recur(b, b->nodes[node].fd, pivot + 1, i1, depth + 1);

}


bool Geometry::intersection(const Ray& d, Vector& P, double &t, MaterialValues &mat, double cur_best_t, int &triangle_id) const
{
	t = cur_best_t;
	bool has_inter = false;
	double t_box_left, t_box_right;
	int best_index = -1;
	bool goleft, goright;
	Vector localP, localN;
	double localt;
	double alpha, beta, gamma;

	Ray invd(d.origin, Vector(1. / d.direction[0], 1. / d.direction[1], 1. / d.direction[2]), d.time);
	char signs[3];
	signs[0] = (invd.direction[0] >= 0) ? 1 : 0;
	signs[1] = (invd.direction[1] >= 0) ? 1 : 0;
	signs[2] = (invd.direction[2] >= 0) ? 1 : 0;

	if (!bvh.bbox.intersection_invd(invd, signs, t_box_left)) return false;
	if (t_box_left > cur_best_t) return false;

	int l[50];
	double tnear[50];
	int idx_back = -1;

	l[++idx_back] = 0;
	tnear[idx_back] = t_box_left;

	while (idx_back >= 0) {

		if (tnear[idx_back] > t) {
			idx_back--;
			continue;
		}
		const int current = l[idx_back--];

		const int fg = bvh.nodes[current].fg;
		const int fd = bvh.nodes[current].fd;

		if (!bvh.nodes[current].isleaf) {
			if (signs[0] == 1) {
				goleft = (bvh.nodes[fg].bbox.intersection_invd_positive_x(invd, signs, t_box_left) && t_box_left < t);
				goright = (bvh.nodes[fd].bbox.intersection_invd_positive_x(invd, signs, t_box_right) && t_box_right < t);
			} else {
				goleft = (bvh.nodes[fg].bbox.intersection_invd_negative_x(invd, signs, t_box_left) && t_box_left < t);
				goright = (bvh.nodes[fd].bbox.intersection_invd_negative_x(invd, signs, t_box_right) && t_box_right < t);
			}

			if (goleft&&goright) {
				if (t_box_left < t_box_right) {
					l[++idx_back] = fd;  tnear[idx_back] = t_box_right;
					l[++idx_back] = fg;  tnear[idx_back] = t_box_left;
				} else {
					l[++idx_back] = fg;  tnear[idx_back] = t_box_left;
					l[++idx_back] = fd;  tnear[idx_back] = t_box_right;
				}
			} else {
				if (goleft) { l[++idx_back] = fg; tnear[idx_back] = t_box_left; }
				if (goright) { l[++idx_back] = fd; tnear[idx_back] = t_box_right; }
			}
		} else {  // feuille
				  //bool go_in = (bvh.nodes[current].bbox.intersection_invd(invd, signs, t_box_left) && t_box_left < t);
				  //if (go_in)
			for (int i = fg; i < fd; i++) {
				if (triangleSoup[i].intersection(d, localP, localt, alpha, beta, gamma)) {
					if (localt < t) {
						int textureId = indices[i].group;
						if (uvs.size()>0 && alphamap.size() > textureId && indices[i].uvi >= 0 && indices[i].uvj >= 0 && indices[i].uvk >= 0) {
							double u = uvs[indices[i].uvi][0] * alpha + uvs[indices[i].uvj][0] * beta + uvs[indices[i].uvk][0] * gamma;
							double v = uvs[indices[i].uvi][1] * alpha + uvs[indices[i].uvj][1] * beta + uvs[indices[i].uvk][1] * gamma;
							u = Texture::wrap(u);
							v = Texture::wrap(v);
							if (alphamap[textureId].getValRed(u, v) < 0.5) continue;
						}
						has_inter = true;
						best_index = i;
						t = localt;
					}
				}
			}
		}

	}

	if (has_inter) {
		int i = best_index;
		triangle_id = best_index;
		triangleSoup[i].intersection(d, localP, localt, alpha, beta, gamma);
		localN.normalize();
		P = localP;

		double u = 0, v = 0;
		int textureId = indices[i].group;
		bool has_uv = false;

		if ((uvs.size() != 0) && (indices[i].group >= 0) && (indices[i].uvi >= 0) && (indices[i].uvj >= 0) && (indices[i].uvk >= 0)) {
			if (!(indices[i].uvi >= uvs.size() || indices[i].uvj >= uvs.size() || indices[i].uvk >= uvs.size())) {
				u = (uvs[indices[i].uvi][0] * alpha + uvs[indices[i].uvj][0] * beta + uvs[indices[i].uvk][0] * gamma);
				v = (uvs[indices[i].uvi][1] * alpha + uvs[indices[i].uvj][1] * beta + uvs[indices[i].uvk][1] * gamma);
				u = Texture::wrap(u);
				v = Texture::wrap(v);
				has_uv = true;
			}
		}
		mat = queryMaterial(textureId, u, v);

		if (!interp_normals || (indices[i].ni == -1)) {
			mat.shadingN = triangleSoup[i].N.getNormalized();
		} else {
			mat.shadingN = normals[indices[i].ni] * alpha + normals[indices[i].nj] * beta + normals[indices[i].nk] * gamma;
			mat.shadingN.normalize();
		}
		if (dot(mat.shadingN, d.direction) > 0 && !transparent) mat.shadingN = -mat.shadingN;
		if (flip_normals) mat.shadingN = -mat.shadingN;


		if ((normal_map.size() != 0) && has_uv && (textureId < normal_map.size())) {

			Vector tangent = (tangents[indices[i].vtxi] * alpha + tangents[indices[i].vtxj] * beta + tangents[indices[i].vtxk] * gamma).getNormalized();
			Vector bitangent = (bitangents[indices[i].vtxi] * alpha + bitangents[indices[i].vtxj] * beta + bitangents[indices[i].vtxk] * gamma).getNormalized();

			Vector NsLocal = normal_map[textureId].getNormal(u, v);
			Vector Ns = NsLocal[0] * tangent + NsLocal[1] * bitangent + NsLocal[2] * mat.shadingN;
			if (Ns[0] == 0. && Ns[1] == 0 && Ns[2] == 0)
				Ns = mat.shadingN;
			Ns.normalize();
			if (!(isnan(Ns[0]) || isnan(Ns[1]) || isnan(Ns[2])))
				mat.shadingN = Ns;
			if (dot(mat.shadingN, d.direction) > 0) mat.shadingN = -mat.shadingN;
		}


		if (vertexcolors.size() != 0) {
			Vector col = (vertexcolors[indices[i].vtxi] * alpha + vertexcolors[indices[i].vtxj] * beta + vertexcolors[indices[i].vtxk] * gamma);
			mat.Kd = col;
			if (display_edges) {
				if ((alpha < 0.05) && (indices[i].showEdges[1]))
					mat.Kd = Vector(0., 0., 0.);
				if ((beta < 0.05) && (indices[i].showEdges[2]))
					mat.Kd = Vector(0., 0., 0.);
				if ((gamma < 0.05) && (indices[i].showEdges[0]))
					mat.Kd = Vector(0., 0., 0.);
			}
		}
		if (facecolors.size() != 0) {
			mat.Kd = facecolors[i];
		}

		if (display_edges) {
			if (edgecolor.size() != 0) {
				if (alpha < 0.05 || beta < 0.05 || gamma < 0.05) {
					int id1, id2;
					if (alpha < 0.05) {
						id1 = std::min(indices[i].vtxj, indices[i].vtxk);
						id2 = std::max(indices[i].vtxj, indices[i].vtxk);
					}
					if (beta < 0.05) {
						id1 = std::min(indices[i].vtxi, indices[i].vtxk);
						id2 = std::max(indices[i].vtxi, indices[i].vtxk);
					}
					if (gamma < 0.05) {
						id1 = std::min(indices[i].vtxj, indices[i].vtxi);
						id2 = std::max(indices[i].vtxj, indices[i].vtxi);
					}

					const std::map<int, Vector> *curmap = &edgecolor[id1];
					const std::map<int, Vector>::iterator it = const_cast<std::map<int, Vector> *>(curmap)->find(id2);
					if (it != const_cast<std::map<int, Vector> *>(curmap)->end())
						mat.Kd = it->second;
					else
						mat.Kd = Vector(0., 0., 0.);
				}
			} else {
				if ((alpha < 0.05) && (indices[i].showEdges[1]))
					mat.Kd = Vector(0., 0., 0.);
				if ((beta < 0.05) && (indices[i].showEdges[2]))
					mat.Kd = Vector(0., 0., 0.);
				if ((gamma < 0.05) && (indices[i].showEdges[0]))
					mat.Kd = Vector(0., 0., 0.);
			}
		}
	}

	return has_inter;
}



bool Geometry::intersection_shadow(const Ray& d, double &t, double cur_best_t, double dist_light) const
{
	t = cur_best_t;
	bool has_inter = false;
	double t_box_left, t_box_right;
	int best_index = -1;
	bool goleft, goright;
	Vector localP, localN;
	double localt;
	double alpha, beta, gamma;

	Ray invd(d.origin, Vector(1. / d.direction[0], 1. / d.direction[1], 1. / d.direction[2]), d.time);
	char signs[3];
	signs[0] = (invd.direction[0] >= 0) ? 1 : 0;
	signs[1] = (invd.direction[1] >= 0) ? 1 : 0;
	signs[2] = (invd.direction[2] >= 0) ? 1 : 0;

	if (!bvh.bbox.intersection_invd(invd, signs, t_box_left)) return false;
	if (t_box_left > cur_best_t || t_box_left>dist_light) return false;

	int l[50];
	double tnear[50];
	int idx_back = -1;

	l[++idx_back] = 0;
	tnear[idx_back] = t_box_left;

	while (idx_back >= 0) {

		if (tnear[idx_back] > t) {
			idx_back--;
			continue;
		}
		const int current = l[idx_back--];

		const int fg = bvh.nodes[current].fg;
		const int fd = bvh.nodes[current].fd;

		if (!bvh.nodes[current].isleaf) {
			goleft = (bvh.nodes[fg].bbox.intersection_invd(invd, signs, t_box_left) && (t_box_left < t) && (t_box_left < dist_light));
			goright = (bvh.nodes[fd].bbox.intersection_invd(invd, signs, t_box_right) && (t_box_right < t) && (t_box_right < dist_light));

			if (goleft&&goright) {
				if (t_box_left < t_box_right) {
					l[++idx_back] = fd;  tnear[idx_back] = t_box_right;
					l[++idx_back] = fg;  tnear[idx_back] = t_box_left;
				} else {
					l[++idx_back] = fg;  tnear[idx_back] = t_box_left;
					l[++idx_back] = fd;  tnear[idx_back] = t_box_right;
				}
			} else {
				if (goleft) { l[++idx_back] = fg; tnear[idx_back] = t_box_left; }
				if (goright) { l[++idx_back] = fd; tnear[idx_back] = t_box_right; }
			}
		} else {  // feuille

			for (int i = fg; i < fd; i++) {
				if (triangleSoup[i].intersection(d, localP, localt, alpha, beta, gamma)) {
					if (localt < t) {
						int textureId = indices[i].group;
						if (uvs.size()>0 && alphamap.size() > textureId && indices[i].uvi >= 0 && indices[i].uvj >= 0 && indices[i].uvk >= 0) {
							double u = uvs[indices[i].uvi][0] * alpha + uvs[indices[i].uvj][0] * beta + uvs[indices[i].uvk][0] * gamma;
							double v = uvs[indices[i].uvi][1] * alpha + uvs[indices[i].uvj][1] * beta + uvs[indices[i].uvk][1] * gamma;
							u = Texture::wrap(u);
							v = Texture::wrap(v);
							if (alphamap[textureId].getValRed(u, v) < 0.5) continue;
						}
						has_inter = true;
						best_index = i;
						t = localt;
						if (t < dist_light*0.999) return true;
					}
				}
			}
		}

	}


	return has_inter;
}


void Geometry::findQuads(int &nbTriangles, int &nbOthers, int &nbRealEdges) {
	std::map<Edge, std::pair<bool, std::vector<int> > > edges_to_faces;

	int nb_hidden = 0;
	nbTriangles = 0;
	for (int i = 0; i < indices.size(); i++) {

		edges_to_faces[Edge(indices[i].vtxi, indices[i].vtxj)].second.push_back(i);
		edges_to_faces[Edge(indices[i].vtxj, indices[i].vtxk)].second.push_back(i);
		edges_to_faces[Edge(indices[i].vtxi, indices[i].vtxk)].second.push_back(i);

		edges_to_faces[Edge(indices[i].vtxi, indices[i].vtxj)].first = indices[i].showEdges[0];
		edges_to_faces[Edge(indices[i].vtxj, indices[i].vtxk)].first = indices[i].showEdges[1];
		edges_to_faces[Edge(indices[i].vtxi, indices[i].vtxk)].first = indices[i].showEdges[2];

		if (indices[i].showEdges[0] && indices[i].showEdges[1] && indices[i].showEdges[2])
			nbTriangles++;
	}
	for (auto it = edges_to_faces.begin(); it != edges_to_faces.end(); ++it) {
		if (!it->second.first) nb_hidden++;

	}
	nbRealEdges = edges_to_faces.size() - nb_hidden;
	int nbFacets = indices.size() - nb_hidden;
	nbOthers = nbFacets - nbTriangles;
}

int Geometry::getNbConnected(int &alsoReturnsNbEdges, int &nonManifoldFaces, int &nbBoundaryEdges) {

	std::map<Edge, std::vector<int> > edges_to_faces;
	for (int i = 0; i < indices.size(); i++) {		
		edges_to_faces[Edge(indices[i].vtxi, indices[i].vtxj)].push_back(i);		
		edges_to_faces[Edge(indices[i].vtxj, indices[i].vtxk)].push_back(i);		
		edges_to_faces[Edge(indices[i].vtxi, indices[i].vtxk)].push_back(i);
	}
	alsoReturnsNbEdges = edges_to_faces.size();

	nonManifoldFaces = 0;
	nbBoundaryEdges = 0;
	std::map<int, std::vector<int> > faces_to_neighbors;
	for (std::map<Edge, std::vector<int> >::iterator it = edges_to_faces.begin(); it != edges_to_faces.end(); ++it) {
		if (it->second.size() > 2) nonManifoldFaces++;
		if (it->second.size() == 1) nbBoundaryEdges++;
		for (int i = 0; i < it->second.size(); i++) {
			for (int j = 0; j < it->second.size(); j++) {
				if (it->second[i] == it->second[j]) continue;
				faces_to_neighbors[it->second[i]].push_back(it->second[j]);
			}
		}
	}

	std::list<int> visited_faces;
	std::vector<bool> visited(indices.size(), false);
	int nbcomp = 0, last_non_visited = 0, nb_visited = 0;

	while (nb_visited != indices.size()) {
		for (int i = last_non_visited; i < indices.size(); i++) {
			if (!visited[i]) {
				visited_faces.push_back(i);
				last_non_visited = i;
				break;
			}
		}

		nbcomp++;
		while (visited_faces.size() != 0) {
			int curface = visited_faces.front();
			visited_faces.pop_front();
			if (visited[curface]) continue;
			visited[curface] = true;
			nb_visited++;

			for (int i = 0; i < faces_to_neighbors[curface].size(); i++) {
				if (!visited[faces_to_neighbors[curface][i]]) {
					visited_faces.push_back(faces_to_neighbors[curface][i]);
				}
			}
		}
	}

	return nbcomp;
}

