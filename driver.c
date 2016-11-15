/***************************************************

	rasterization.cpp

	by Devon Powell

***************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "r3d.h"
//#include "utils.h"

#define popc(x) __builtin_popcount(x)

typedef struct {
	r3d_real det[256];	
	r3d_int sgn[256];
} r3d_tet_data;

// position of the lowest set bit
r3d_int lbpos(r3d_int x) {
	// TODO: is there an even faster implementation??
	const unsigned lbtable[32] = {0,  1, 28,  2, 29, 14, 24,  3, 30, 22, 20, 15, 25, 17,  4,  8,
		31, 27, 13, 23, 21, 19, 16,  7, 26, 12, 18,  6, 11,  5, 10,  9};
	return lbtable[((x&-x)*0x077CB531)>>27];
}

r3d_int snoob(r3d_int x) {
	unsigned rightOne;
	unsigned nextHigherOneBit;
 	unsigned rightOnesPattern;
	rightOne = x&-x;    
	nextHigherOneBit = x + rightOne;
	rightOnesPattern = x ^ nextHigherOneBit;
	rightOnesPattern = (rightOnesPattern)/rightOne;
	rightOnesPattern >>= 2;
	return nextHigherOneBit | rightOnesPattern;
}

void print_binary(r3d_int bits, r3d_int n) {
	r3d_int i;
	//for(i = n-1; i >= 0; --i) 
	for(i = 0; i < n; ++i) 
		printf("%d", (bits>>i)&1);
	printf("\n");
}

r3d_int process_intersections(r3d_int bits, r3d_int parity, r3d_rvec3* verts, r3d_int nverts, r3d_tet_data* tdata) {

	static r3d_int num_int = 0;
	r3d_int i, v, tetmask, tet, tetpar, tmppar, tint, mask, nset, dim, intersection, tmpbits;
	r3d_int andcmp, orcmp;
	r3d_real vol, bary[2][4];
	r3d_rvec3 vpos;

	// TODO: this may not be great in CUDA...
	dim = popc(bits)-1;
	if(dim == 3) {
		return parity^tdata->sgn[bits];
	}

	// decoupling the recursion parity for the two tets
	// compute the barycentric coordinates of the intersection for both tets simultaneously
	intersection = 1;
	memset(bary, 0, sizeof(bary));
	for(tet = 0; tet < 2; ++tet) {
		tetmask = 0x0F<<(4*tet);
		tetpar =  tdata->sgn[tetmask]; // TODO: do we need tetpar??
		andcmp = 1; orcmp = 0; vol = 0.0;
		for(nset = 0, tmpbits = bits&tetmask; tmpbits; nset++, tmpbits ^= mask) {
			mask = tmpbits&-tmpbits;
			tmppar = (tetpar^!(nset&1)); // wtf // TODO: do we need parity really? 
			tint = process_intersections(bits^mask, parity^tmppar, verts, nverts, tdata);
			andcmp &= tint; orcmp |= tint;

			// TODO: only calculate barycentric coords for 3D!!!!

			// TODO: a bit more elegant...
			// TODO: fix the lbpos function!!
			v = lbpos(mask>>(4*tet)); // get the position of the lowest bit in the mask
			bary[tet][v] = (1-2*(parity^tmppar))*tdata->det[bits^mask]; 
			vol += bary[tet][v]; 
		}
		//if(vol[tet] == 0.0) printf("NOOOOO! zero volume!!!\n");
		for(v = 0; v < 4; ++v) bary[tet][v] /= vol; 
		intersection &= andcmp|~orcmp; // if the orientation is consistent, weaker constraint than positive
	}


	if(dim == 4 && intersection) {

		++num_int;

		printf("\n");

		// we have found some intersection between the tets
		printf("Found 4-intersection, dim = %d, bits = ", dim); 
		print_binary(bits, 8);

		for(tet = 0; tet < 2; ++tet) {
			printf(" tet %d, vol = , verts = ", tet);
			tetmask = 0x0F<<(4*tet);
			for(nset = 0, tmpbits = bits&tetmask; tmpbits; nset++, tmpbits ^= mask) {
				mask = tmpbits&-tmpbits;
				v = lbpos(mask>>(4*tet));
				printf("%d ", v);
			}
			printf("\n");


			for(i = 0; i < 3; ++i) {
				vpos.xyz[i] = 0.0;
				for(v = 0; v < 4; ++v) vpos.xyz[i] += bary[tet][v]*verts[v+4*tet].xyz[i];
			}
			printf(" intersection location = %f %f %f\n", vpos.x, vpos.y, vpos.z);
			printf("  # %d\n", num_int);

		}

	}

	return intersection;
}


void raster_test_det() {

	r3d_int v, i;

	r3d_rvec3 verts[8];
	memset(&verts, 0, sizeof(verts));
	for(v = 1; v < 4; ++v) verts[v].xyz[v-1] = 1.0;
	for(v = 5; v < 8; ++v) verts[v].xyz[v-5] = -1.0;
	verts[0].x = -0.3;
	verts[0].y = -0.3;
	verts[0].z = -0.3;
	verts[4].x = 0.2;
	verts[4].y = 0.2;
	verts[4].z = 0.2;
	//for(v = 0; v < 4; ++v) // experimental!
	//for(i = 0; i < 3; ++i) verts[v].xyz[i] = -0.1; 
	for(v = 0; v < 8; ++v)
		printf("verts[%d] = %f %f %f\n", v, verts[v].x, verts[v].y, verts[v].z);

	r3d_tet_data tdata;
	r3d_int nvb, allbits;
	r3d_rvec3 dv[4];

	// fill in all possible determinants, indexing by all numbers with 4 bits set 
	// get the four vertex indices, i.e. which four bits were set
	// build the determinant array
	memset(&tdata, 0, sizeof(tdata));
	for(allbits = 0x0F; allbits <= 0xF0; allbits = snoob(allbits)) {
		for (i = 0, nvb = 0; nvb < 4; ++i)
			if((allbits>>i)&1) dv[nvb++] = verts[i]; 
		tdata.det[allbits] = r3d_orient(dv);
		tdata.sgn[allbits] = (tdata.det[allbits] > 0);
	}

	process_intersections(0xFF, 1, verts, 8, &tdata);

	for(v = 0; v < 8; ++v)
		printf("verts[%d] = %f %f %f\n", v, verts[v].x, verts[v].y, verts[v].z);

}


int main(int argc, char **argv) {
	setbuf(stdout, NULL);
#ifdef CUDA
	printf("CUDA enabled.\n");

	r3d_int nside = 32;
	r3d_dvec3 sz = {{nside, nside, 1}}; 

	r3d_real* data = malloc(nside*nside*sizeof(r3d_real));
	cur3d_fill(data, sz, 2);

	FILE* f = fopen("out.dat", "wb");
	if(!f) return 0;

	fwrite(data, sizeof(r3d_real), nside*nside, f);



	close(f);
	free(data);
#else
	printf("CUDA disabled.\n");
	raster_test_det();
#endif



	return 0;    
}



