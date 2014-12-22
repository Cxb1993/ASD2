#ifndef _DATA_H_
#define _DATA_H_

#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "hdf5.h"
#include "../lib/common/jsonxx/jsonxx.h"

using namespace jsonxx;
using namespace std;

class Data {
public:
	Data();
	~Data();
	
	// Template is more elegant solution, 
	// however, I didn't use it to make easy to understand code for c++ newbie
	static int* Allocate1Di(int nx);
	static float* Allocate1Df(int nx);
	static double* Allocate1Dd(int nx);
	static long double* Allocate1Dld(int nx);

	static int Deallocate1Di(int *to_deallocate);
	static int Deallocate1Df(float *to_deallocate);
	static int Deallocate1Dd(double *to_deallocate);
	static int Deallocate1Dld(long double *to_deallocate);

	static int** Allocate2Di(int nx, int ny);
	static float** Allocate2Df(int nx, int ny);
	static double** Allocate2Dd(int nx, int ny);
	static long double** Allocate2Dld(int nx, int ny);

	static int Deallocate2Di(int **to_deallocate);
	static int Deallocate2Df(float **to_deallocate);
	static int Deallocate2Dd(double **to_deallocate);
	static int Deallocate2Dld(long double **to_deallocate);

	static int*** Allocate3Di(int nx, int ny, int nz);
	static float*** Allocate3Df(int nx, int ny, int nz);
	static double*** Allocate3Dd(int nx, int ny, int nz);
	static long double*** Allocate3Dld(int nx, int ny, int nz);

	static int Deallocate3Di(int ***to_deallocate);
	static int Deallocate3Df(float ***to_deallocate);
	static int Deallocate3Dd(double ***to_deallocate);
	static int Deallocate3Dld(long double ***to_deallocate);

	static hid_t OpenHDF5forRead(string file_name);
	static hid_t OpenHDF5forWrite(string file_fname);
	static int CloseHDF5(hid_t h5_file);

	static int OpenPLTforWrite(string file_name, std::ofstream& outfile);
	static int OpenPLTforAppend(string file_name, std::ofstream& outfile);
	static int OpenPLTforRead(string file_name, std::ifstream& infile);

	static int WritePLTBasic1D(std::ofstream& of_plt, int nx);
	static int WritePLTBasic2D(std::ofstream& of_plt, int nx, int ny);
	static int WritePLTBasic3D(std::ofstream& of_plt, int nx, int ny, int nz);
	
	static int ClosePLTWrite(std::ofstream& file);
	static int ClosePLTRead(std::ifstream& file);

	template <typename T>
	T getNumericValFromJson(jsonxx::Object io, std::string name);
	bool getBooleanValFromJson(jsonxx::Object io, std::string name);
	std::string getStringValFromJson(jsonxx::Object io, std::string name);
};

#endif
