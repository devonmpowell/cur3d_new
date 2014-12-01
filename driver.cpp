/***************************************************

	rasterization.cpp

	by Devon Powell

***************************************************/

#include <cstdio>
#include <vector>
#include <string>
#include <cstdint>
#include <cfloat>
#include <cmath>
#include "HDF_IO.hh"

#include "curas.ch"


// custom macros and typedefs
using namespace std;
typedef uint32_t luint;
#define LUINT_MAX UINT32_MAX


// forward declarations
//void read_hdf5(string filename, string fieldname);
void write_hdf5(string filename, string fieldname, vector<float> &data, luint nx, luint ny, luint nz);


int main(int argc, char **argv) {
	setbuf(stdout, NULL);

	// data arrays
	luint nx, ny, nz, ntot;
	vector<float> field;
	nx = ny = nz = 128;
	ntot = nx*ny*nz;
	field.assign(ntot, 0.0);


	curas_fill(field.data(), nz, ny, nz);


	write_hdf5("/lustre/ki/pfs/dmpowel1/curas_out/test.hdf5", "RHO", field, nx, ny, nz);



	return 0;    
}

void write_hdf5(string filename, string fieldname, vector<float> &data, luint nx, luint ny, luint nz) {
	HDFCreateFile(filename);
	luint dims[3] = {nx, ny, nz};
	HDFWriteDataset3D(filename, fieldname, dims, data);
	return;
}


/*
// read in a multidimensional hdf5 file
void read_hdf5(string filename, string fieldname) {
	vector<int> dims;
	HDFGetDatasetExtent(filename, fieldname, dims);
	nx = dims[0]; ny = dims[1]; nz = dims[2];
	ntot = nx*ny*nz;
	HDFReadDataset(filename, fieldname, field);
	return;
}
*/

