#pragma once

#define BRDF_SAMPLING_RES_THETA_H       90
#define BRDF_SAMPLING_RES_THETA_D       90
#define BRDF_SAMPLING_RES_PHI_D         360
#define RED_SCALE (1.0/1500.0)
#define GREEN_SCALE (1.15/1500.0)
#define BLUE_SCALE (1.66/1500.0)

void lookup_brdf_val(double* brdf, double theta_in, double fi_in,
	double theta_out, double fi_out,
	double& red_val, double& green_val, double& blue_val);

bool read_brdf(const char *filename, double* &brdf);
