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

#if 0

	// TODO: is there a more bit-hacky may to do this??
	typedef struct {
		r3d_int parity;
		r3d_int bits;
		//r3d_int vinds[4];
		//r3d_rvec3 pos;
		//r3d_rvec3 ltd[3];
	} r3d_ltd;
	r3d_ltd stack[128]; // no idea how big
	r3d_ltd curltd, tmpltd; 
	r3d_ltd tetltd[24];
	r3d_int ntl; 
	r3d_int nstack = 0;

	// initialize the base case, with both tets
	tmpltd.parity = 1;
	tmpltd.bits = 0xFF;
	stack[nstack++] = tmpltd;
	nint = 0, ntl = 0;
	while(nstack) {

		// pop the stack
		// stop if we are at a point (0-dimensional), quit
		curltd = stack[--nstack];
		bits = curltd.bits;
		dim = popc(bits)-1;
		//printf("found node with parity %d, bits = ", curltd.parity); print_binary(bits, 8);

		// TODO: this may not be great in CUDA...
		if(dim == 3) {
			printf("found intersection with parity %d, bits = ", curltd.parity); print_binary(bits, 8);
			nint++;
			continue;
		}

		// otherwise cycle through set bits to lower dimensionality 
		// push them to the stack 
		for(nvb = 0; bits; nvb++, bits ^= mask) {
			mask = bits&-bits; // get the rightmost set bit 
			tmpltd = curltd;
			tmpltd.bits ^= mask; 
			tmpltd.parity ^= nvb&1; 
			stack[nstack++] = tmpltd;
		}



		// TODO: this may not be great in CUDA...
		if(!dim) {
			v = lbpos(bits);
			curltd.pos = verts[v]; // index the set bit
			curltd.vinds[dim] = v;
			tetltd[ntl++] = curltd;
			continue;
		}

		// otherwise cycle through set bits to lower dimensionality 
		// push them to the stack 
		nvb = 0;
		while(bits) {
			mask = bits&-bits; // extract least significant bit on a 2s complement machine
			bits ^= mask;  // toggle the bit off
			v = lbpos(mask); // get the index of the set bit
			tmpltd = curltd;
			tmpltd.ltd[dim-1] = verts[v];
			tmpltd.vinds[dim] = v;
			tmpltd.bits ^= mask; 
			tmpltd.parity ^= (nvb++)&1; 
			stack[nstack++] = tmpltd;
		}
	}
	printf("num ltds = %d\n", nint);
#endif

#if 0
	for(j = 0; j < ntl; ++j) {
		curltd = tetltd[j];
		//pritnf("stack[%d] = %f %f %f\n", i, stack[i].pos.x, stack[i].pos.y, stack[i].pos.z);

			for(i = 0; i < 3; ++i) {
				for(j = 0; j < 3; ++j)
					curltd.ltd[i].xyz[j] -= curltd.pos.xyz[j];

				// gram-schmidt
				// TODO: is it dangerous to sweep upward in dimension??
				// May introduce errors to produce face normals via different calculations here..
				// TODO: can we count on normals to be robust if we calculate them top-down??
				for(j = 0; j < i; ++j) {
					dotfac = 0.0;
					for(k = 0; k < 3; ++k)
						dotfac += curltd.ltd[j].xyz[k]*curltd.ltd[i].xyz[k];
					for(k = 0; k < 3; ++k)
						curltd.ltd[i].xyz[k] -= dotfac*curltd.ltd[j].xyz[k];
				}
				len = 0.0;
				for(j = 0; j < 3; ++j)
					len += curltd.ltd[i].xyz[j]*curltd.ltd[i].xyz[j];
				invlen = 1.0/sqrt(len);
				for(j = 0; j < 3; ++j)
					curltd.ltd[i].xyz[j] *= invlen;
			}

			printf("Current ltd: parity = %d, vert order = %d %d %d %d, bits = ", curltd.parity, curltd.vinds[0], curltd.vinds[1], curltd.vinds[2], curltd.vinds[3]); print_binary(curltd.bits, 4);
			printf(" pos = %f %f %f\n", curltd.pos.x, curltd.pos.y, curltd.pos.z);
			for(i = 0; i < 3; ++i)
				printf("       %f %f %f\n", curltd.ltd[i].x, curltd.ltd[i].y, curltd.ltd[i].z);
	}
#endif

	return;
	
	// next, count through all 5-vertex permutations to get intersections
	nint = 0;
	r3d_int predicate;
	for(allbits = 0x1F; allbits <= 0xF8; allbits = snoob(allbits)) {
		print_binary(allbits, 8);
		//for(i = 0; i < 8; ++i) {
			//if(!((allbits>>i)&1)) continue;
			//bits = allbits^(1<<i);
			//predicate = tdata.sgn[bits];
			//printf(" predicate = %d, det = %.5e, bits = ", predicate, tdata.det[bits]); print_binary(bits, 8);
		//}
		
		// allbits set for each tet
		bits0 = allbits&0x0F;
		bits1 = allbits&0xF0;
		nv0 = popc(bits0);
		nv1 = popc(bits1);
		printf(" %d verts from tet 0, %d from tet 1\n", nv0, nv1);

		//if(nv0 == 1)
		
		++nint;
	}
	printf("num. inersection checks = %d\n", nint);

	return;



#if 0

	// step through all intersection possibilities
	r3d_int bitsmin, bitsmax, imask, jmask, idim, jdim, bits; 

	//bits = 0x0;
	r3d_int nset,  b, sgn;
	r3d_int dim, ndim = 3;

	for(dim = 0; dim < 3; ++dim) {

		printf("dim = %d of %d\n", dim, ndim);
	
		idim = dim;
		imask = (1<<(idim+1))-1;
		jdim = ndim-dim;
		jmask = (1<<jdim)-1;
		sgn = 0;
		r3d_int isgn, jsgn;
		isgn = 0;
		for(i = imask; i <= imask<<(ndim-idim); i = snoob(i)) {

			isgn ^= 1;

			//printf("nin")
			inout = 0;
			jsgn = 0;
			for(j = jmask; j <= jmask<<(ndim-jdim+1); j = snoob(j)) {
				//printf("popc(j) = %d\n", popc(j));
				jsgn ^= 1;

				
				bits = i|(j<<4);

				inout += isgn^jsgn^tdata.sgn[bits];

				printf("ij = "); 
				print_binary(bits, 8);
				//printf("  det = %f, sgn = %d\n", tdata.det[bits], tdata.sgn[bits]);
	
			}
				//printf("inout = %d\n", inout);
			if(inout == jdim+1 || inout == 0) {
				printf("Found intersection, dim = %d vs. %d\n", idim, jdim);
				printf("i = "); 
				print_binary(i&0x0F, 8);
				printf("inout = %d\n", inout);
			}
		} 


			//for(f = 0; f < 4; ++f) {
				//// compute the sign flip
				//fsgn = f&0x1;
				//inout |= (fsgn^tdata.sgn[i])<<f;
				//printf("vert[%d] = %f %f %f is %1d face[%d], flip = %d^%d, bits = ", v, verts[v].x, verts[v].y, verts[v].z, inout, f, fsgn, tdata.sgn[i]);
				//print_binary(i, 8);
			//}
			//if(inout == ((1<<4)-1) || inout == 0x0) {
				//printf("this vertex is inside the other tet!\n");
			//}
	
		
	}



#elif 0

	//const static r3d_int r3d_tet_faces_to_verts[4][3] = {{0, 1, 2}, {3, 0, 2}, {1, 3, 2}, {3, 1, 0}};
	const static r3d_int r3d_tet_faces_to_verts[4][3] = {{1, 3, 2}, {3, 0, 2}, {3, 1, 0}, {0, 1, 2}};
#define r3d_tet_face_sgn(x) (x&1)
	//const static r3d_int r3d_tet_face_sgn[4] = {1, 0, 1, 0};
	//const static r3d_int r3d_tet_faces_to_verts[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};
	//const static r3d_int r3d_tet_face_sgn[4] = {0, 1, 0, 1};

	// TODO: is th sign flip for faces a pattern or just coincidence? 

	// is v on the same side of all the faces of f? 
	r3d_int vbits, fbits, vsgn;
	for(v = 0; v < 4; ++v) {
		vbits = (1<<v);
		vsgn = v&1;
		inout = 0;
		for(f = 0; f < 4; ++f) {
			// compute the index for this determinant
			// compute the sign flip
			fbits = 0x0F^(1<<f);
			fsgn = r3d_tet_face_sgn(f);
			i = (fbits<<4)|vbits;
			inout += vsgn^fsgn^tdata.sgn[i];
			printf("vert[%d] = %f %f %f is %1d face[%d], flip = %d^%d, bits = ", v, verts[v].x, verts[v].y, verts[v].z, inout, f, fsgn, tdata.sgn[i]);
			print_binary(i, 8);
		}
		if(inout == 4 || inout == 0) {
			printf("this vertex is inside the other tet!\n");
		}
	}

	// do any edges of the tet intersect faces of the second?
	r3d_int eme, fme, fyou, e, flip;
	r3d_int sgn, ebits, bits, mesgn, yousgn;
	r3d_int pp0, pp1, p0, p1;
	for(fme = 0; fme < 4; ++fme)
	for(eme = 0; eme < 3; ++eme) {
		p0 = r3d_tet_faces_to_verts[fme][eme];
		p1 = r3d_tet_faces_to_verts[fme][(eme+1)%3];
		mesgn = (p0 < p1);
		for(f = 0; f < 4; ++f) {

			// first check to see if the edge crosses this face
			// skip if both verts are on the same side of f
			bits = 0x0;
			yousgn = r3d_tet_face_sgn(f);
			for(e = 0; e < 3; ++e) bits |= 1<<(4+r3d_tet_faces_to_verts[f][e]);
			bits &= 0xF0; bits |= 1<<p0; 
			if(tdata.sgn[bits]^yousgn) continue;
			bits &= 0xF0; bits |= 1<<p1;
			if(!(tdata.sgn[bits]^yousgn)) continue;
			
			// test whether the edge passes through the face polygon
			bits = ((1<<p0)|(1<<p1))<<4;
			inout = 0;
			for(e = 0; e < 3; ++e) {
				pp0 = r3d_tet_faces_to_verts[f][e];
				pp1 = r3d_tet_faces_to_verts[f][(e+1)%3];
				bits &= 0xF0; bits |= 1<<pp0; bits |= 1<<pp1;
				yousgn = (pp0 < pp1);
				inout += yousgn^tdata.sgn[bits];
			}
			if(inout != 3) continue;
				
			printf("face %d, edge %d-%d passes thru face %d!\n", fme, p0, p1, f);
		}
	}
#else
	printf(" ----------------- \n");

	//r3d_geometry_iterator gitme, gityou;
	//r3d_geometry_iterator_init(&gitme, &verts[0]);
	//r3d_geometry_iterator_init(&gityou, &verts[4]);

	//while(r3d_geometry_iterator_next(&gitme, 0)) {
	
		//printf("next!\n");
	
	
	//}


#endif
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



