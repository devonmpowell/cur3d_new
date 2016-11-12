/***************************************************

	rasterization.cpp

	by Devon Powell

***************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "cur3d.h"
//#include "utils.h"

#define popc(x) __builtin_popcount(x)


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
	unsigned next = 0;
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

	r3d_int mask, nset, dim, intersection, tmpbits;
	r3d_int tint, v, bits0, bits1, tet, i;

	// TODO: this may not be great in CUDA...
	dim = popc(bits)-1;
	if(dim == 3) {
		return parity^tdata->sgn[bits];
	}

	// decoupling the recursion parity for the two tets
	intersection = 1;	
	for(tet = 0; tet < 2; ++tet) {
		for(nset = 0, tmpbits = bits&(0x0F<<(4*tet)); tmpbits; nset++, tmpbits ^= mask) {
			mask = tmpbits&-tmpbits;
			tint = process_intersections(bits^mask, parity^(nset&1), verts, nverts, tdata);
			intersection &= tint;
		}
	}

	if(dim == 4 & intersection) {

		// we have found some intersection between the tets


		printf("Found 4-intersection, bits = ", parity); 
		print_binary(bits, 8);

		for(tet = 0; tet < 2; ++tet) {
			printf(" tet %d, verts = ", tet); 
			for(nset = 0, tmpbits = bits&(0x0F<<(4*tet)); tmpbits; nset++, tmpbits ^= mask) {
				mask = tmpbits&-tmpbits;
				printf("%d ", lbpos(mask>>(4*tet)));
			}
			printf("\n");
		}
	}

	return intersection;
}


void raster_test_det() {

	r3d_int v, i, c, k, j, z, f, g, fsgn;
	r3d_int curbits, shift, mask, nint;
	r3d_real len, invlen, dotfac;

	r3d_rvec3 verts[8];
	memset(&verts, 0, sizeof(verts));
	for(v = 1; v < 4; ++v)
		verts[v].xyz[v-1] = 1.0;
	for(v = 5; v < 8; ++v)
		verts[v].xyz[v-5] = -1.0;
	verts[0].x = -0.2;
	verts[0].y = -0.2;
	verts[0].z = -0.2;
	verts[4].x = 0.2;
	verts[4].y = 0.2;
	verts[4].z = 0.2;


	// experimental!
	//for(v = 0; v < 4; ++v)
	//for(i = 0; i < 3; ++i) verts[v].xyz[i] = -0.1; 


	for(v = 0; v < 8; ++v)
		printf("verts[%d] = %f %f %f\n", v, verts[v].x, verts[v].y, verts[v].z);

	r3d_tet_data tdata;
	r3d_int vertinds[4];
	r3d_int nvb, allbits, bits, dim, bits0, bits1, nv0, nv1;
	r3d_int inout;
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
//#else
	printf("CUDA disabled.\n");
	raster_test_det();
#endif



	return 0;    
}



