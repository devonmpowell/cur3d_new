/***************************************************

	rasterization.cpp

	by Devon Powell

***************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "r3d.h"
//#include "utils.h"


#define R3D_TET 0
#define R3D_CUBE 1
#define FUZZ 1.0e-30

#define popc(x) __builtin_popcount(x)

#define zero(x) memset(&(x), 0, sizeof(x));

#define nmom(order) ((((order+1)*(order+2)*(order+3))/6))

// position of the lowest set bit
r3d_int lbpos(unsigned int x) {
	// TODO: seems not to work for bit positions higher than 4...
	// TODO: is there an even faster implementation??
	const unsigned lbtable[32] = {0,  1, 28,  2, 29, 14, 24,  3, 30, 22, 20, 15, 25, 17,  4,  8,
		31, 27, 13, 23, 21, 19, 16,  7, 26, 12, 18,  6, 11,  5, 10,  9};
	return lbtable[((x&-x)*0x077CB531U)>>27];
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
	for(i = 0; i < n; ++i) 
		printf("%d", (bits>>i)&1);
	printf("\n");
}

// TODO: compute exactly how big these arrays must be!
typedef struct {
	
	// determinant and sign predicates
	// using enough bits to capture all vertices
	r3d_rvec3 allverts[16];
	r3d_int poffs[2], ptypes[2], pmasks[2];
	r3d_real det[4096];	
	r3d_int sgn[4096];

	// candidate vertices
	r3d_real candidates[1024][16]; 
	r3d_int ncand;

	// all simplices and their indices
	r3d_int simpinds[512][4];
 	r3d_int nsimp; 

} r3d_intersection;

typedef struct {
	r3d_rvec3 verts[8];
	r3d_int ptype;
} r3d_tinypoly;


void r3d_reduce_intersection(r3d_real* moments, int order, r3d_intersection* tdata);
void r3d_process_intersection(r3d_tinypoly* polys, r3d_int npoly);
r3d_int r3d_recurse_intersection(r3d_int dim, r3d_int bits, r3d_int parity, r3d_int* vertinds, r3d_real* bpos, r3d_intersection* tdata);
//void r3d_recurse_setbits(r3d_int bits, r3d_tet_data* tdata);

int main(int argc, char **argv) {
	
	setbuf(stdout, NULL);
	r3d_int v, i;

	r3d_rvec3 *curverts; 
	r3d_tinypoly *curpoly; 
	r3d_tinypoly polys[2];
	memset(polys, 0, sizeof(polys));

#if 0
	// initialize one cube and one tetrahedron 
	// hacky way to initialize cube vertices
	// in row-major order
	curpoly = &polys[0];
	curverts = polys[0].verts;
	curpoly->ptype = R3D_CUBE;
	for(v = 0; v < 8; ++v) {
		for(i = 0; i < 3; ++i) {
			curverts[v].xyz[2-i] = 1.0*((v>>i)&1); 
		}
	}

	curpoly = &polys[1];
	curverts = polys[1].verts;
	curpoly->ptype = R3D_TET;	
	for(v = 0; v < 4; ++v)
	for(i = 0; i < 3; ++i)
		curverts[v].xyz[i] = 0.5-0.13;
	for(v = 1; v < 4; ++v) {
		curverts[v].xyz[v-1] += 4;
	}

#else
	// initialize a pair of tetrahedra
	r3d_int tet;
	for(tet = 0; tet < 2; ++tet) {
		curpoly = &polys[tet];
		curverts = polys[tet].verts;
		curpoly->ptype = R3D_TET;
		for(v = 0; v < 4; ++v)
		for(i = 0; i < 3; ++i)
			curverts[v].xyz[i] = (2*tet-1)*0.5-0.13;
		for(v = 1; v < 4; ++v) {
			curverts[v].xyz[v-1] += (1-2*tet)*4;
		}
	}
#endif



#ifdef CUDA
	printf("CUDA enabled.\n");
	//r3d_int nside = 32;
	//r3d_dvec3 sz = {{nside, nside, 1}}; 

	//r3d_real* data = malloc(nside*nside*sizeof(r3d_real));
	//cur3d_fill(data, sz, 2);

	//FILE* f = fopen("out.dat", "wb");
	//if(!f) return 0;

	//fwrite(data, sizeof(r3d_real), nside*nside, f);

	//close(f);
	//free(data);

#else
	printf("CUDA disabled.\n");
#endif

	r3d_process_intersection(polys, 2);
	
	return 0;    
}

void r3d_process_intersection(r3d_tinypoly* polys, r3d_int npoly) {

	r3d_int bits, allbits, minpred, maxpred;
	r3d_int istack[8];
	r3d_real ibpos[16];
	r3d_real moments[10];
	r3d_rvec3 dv[4];
	r3d_intersection tdata;
	memset(&tdata, 0, sizeof(tdata));


	// copy polyhedron data to a single buffer
	// simultaneously build an indexing array
	// TODO: Can we circumvent an indexing step here?? 
	// TODO: Can we save some comparisons here??
	// TODO: Can we not copy?
	r3d_int nv, cumnv, p, v, i, nvb; 
	bits = 0;
	for(p = 0, cumnv = 0; p < npoly; ++p, cumnv += nv) {

		// copy the vertices into a common buffer for fast bit indexing
		nv = 4+4*polys[p].ptype; // TODO: HACK
		memcpy(&tdata.allverts[cumnv], polys[p].verts, nv*sizeof(r3d_rvec3));
		tdata.ptypes[p] = polys[p].ptype;
		tdata.poffs[p] = cumnv;
		tdata.pmasks[p] = ((1<<(cumnv+nv))-1)^((1<<cumnv)-1);
		printf("Polyhedron %d, type = %d, verts = \n", p, polys[p].ptype);
		for(v = 0; v < nv; ++v) 
			printf(" %f %f %f\n", tdata.allverts[cumnv+v].x, tdata.allverts[cumnv+v].y, tdata.allverts[cumnv+v].z);
			//printf("  (%f %f %f)\n", tdata.allverts[cumnv+v].x, tdata.allverts[cumnv+v].y, tdata.allverts[cumnv+v].z);
	}
	
	// fill in all possible determinants, indexing by all numbers with 4 bits set 
	// get the four vertex indices, i.e. which four bits were set
	// build the determinant array
	// TODO: figure out a good way to do this without
	// ever visiting pairs of vertices that we never use!
	allbits = (1<<cumnv)-1; 
	minpred = (allbits^(allbits<<4))&allbits;
	maxpred = (allbits^(allbits>>4))&allbits;
	for(bits = minpred; bits <= maxpred; bits = snoob(bits)) {

		// TODO: bit-hackiness
		//for(iset = 0, tbits = bits&(0x0F<<(4*tet)); tbits; ipar ^= 1, tbits ^= imask4) {
		for (i = 0, nvb = 0; nvb < 4; ++i)
			if((bits>>i)&1) dv[nvb++] = tdata.allverts[i]; 
		tdata.det[bits] = r3d_orient(dv);

		// ensure that the determinant is never identically zero
		// TODO: can we do better here for degenerate tets?? I.e., an extended Cayley-menger determinant??
		// TODO: use edge lengths as a distance proxy??
		tdata.sgn[bits] = (tdata.det[bits] >= 0.0);
		tdata.det[bits] += FUZZ*tdata.sgn[bits]; 
	}

	// recurseively process the intersection, then compute the moments 
	// of the resulting polytope
	r3d_recurse_intersection(3, allbits, 1, istack, ibpos, &tdata);
	r3d_reduce_intersection(moments, 0, &tdata);
}


// get the next bit mask for the given dimension
void r3d_next_mask(r3d_int* mask, r3d_int* tbits, r3d_int* tpar, r3d_int ptype, r3d_int dim) {
	if(ptype == R3D_TET) {
		*mask = (*tbits)&-(*tbits);
	}
	else if(ptype == R3D_CUBE) {
		//printf("Cube!\n");
		*mask = (*tbits)&-(*tbits);
	}
	(*tpar) ^= 1;
	(*tbits) ^= (*mask);
}

r3d_int r3d_recurse_intersection(r3d_int dim, r3d_int bits, r3d_int parity, r3d_int* vertinds, r3d_real* bpos, r3d_intersection* tdata) {

	r3d_int v, p, poff, intersection;
	r3d_int cmp, andcmp, orcmp, candind;;
	r3d_int iset, tpar, ibits, tbits, tmask; // for iterating over set bits 
	r3d_real vol;
	r3d_real ibpos[16];

	// if we are at D+2 vertices, check for an intersection
	//dim = popc(bits)-5;
	memset(bpos, 0, 8*sizeof(r3d_real));
	candind = tdata->ncand++; // reserve this index for later 
	vertinds[dim] = candind;
	if(!dim) {
		// compute the mass coordinates of the intersection for both polytopes simultaneously
		// see if we have found some intersection
		intersection = 1;
		for(p = 0; p < 2; ++p) {
			poff = tdata->poffs[p];
			andcmp = 1; orcmp = 0; vol = 0.0;
			for(tpar = parity, tbits = (bits&tdata->pmasks[p])>>poff; tbits;) {
				r3d_next_mask(&tmask, &tbits, &tpar, tdata->ptypes[p], dim); 
				ibits = (tmask<<poff)^bits;
				cmp = tpar^tdata->sgn[ibits]; 
				andcmp &= cmp; orcmp |= cmp;

				// TODO: how to get the bary coordinates right for each tet??
				v = lbpos(tmask); // TODO: a bit more elegant, fix the lbpos function!!
				bpos[poff+v] = (2*tpar-1)*tdata->det[ibits]; 
				vol += bpos[poff+v]; 
			}

			// TODO: NOT 4, whatever nverts is!
			for(v = 0; v < 4; ++v) bpos[poff+v] /= vol;
			intersection &= andcmp|!orcmp; // if the orientation is consistent we have an intersection
		}
		if(intersection) 
			memcpy(tdata->simpinds[tdata->nsimp++], vertinds, 4*sizeof(r3d_int));
	}
	else {
		// otherwise, recurse on the remaining set bits
		// get the newest vertex position by averaging the children. 
		for(p = 0, iset = 0; p < 2; ++p) {
			poff = tdata->poffs[p];
			for(tpar = parity, tbits = (bits&tdata->pmasks[p])>>poff; tbits;) {
				r3d_next_mask(&tmask, &tbits, &tpar, tdata->ptypes[p], dim); 
				cmp = r3d_recurse_intersection(dim-1, (tmask<<poff)^bits, tpar, vertinds, ibpos, tdata);
				iset += cmp;
				if(cmp) for(v = 0; v < 8; ++v) bpos[v] += ibpos[v]; // TODO: slightly cleaner for CUDA?
			}
		}
		for(v = 0; v < 8; ++v) bpos[v] /= iset; 
		intersection = (iset > 0);
	}
	// save the position for this candidate vertex 
	if(intersection) for(v = 0; v < 8; ++v) 
		tdata->candidates[candind][v] = bpos[v];
	return intersection;
}


void r3d_reduce_intersection(r3d_real* moments, int order, r3d_intersection* tdata) {

	r3d_int i, j, v, s, tet;  
	r3d_real mass;
	r3d_rvec3 v0, v1, v2;

	r3d_rvec3 vpos[4];
	r3d_real vtot[2];
	zero(vtot);

	// zero the moments and sum over all simpleces in the buffer
	//memset(moments, 0, nmom(order)*sizeof(double));
	for(s = 0; s < tdata->nsimp; ++s) {

		r3d_int sgn =  1; //tdata->simplices[s].sgn;
		r3d_int* inds = tdata->simpinds[s];

		printf("Decomposed simplex %d, sign = %d, vertinds =", s, sgn);
		for(v = 0; v < 4; ++v) printf(" %d", inds[v]);
		printf("\n");
		// find the vertex location!

		for(tet = 0; tet < 2; ++tet) {

			// get the real-space positon and volume (just for testing)
			zero(vpos);
			for(v = 0; v < 4; ++v)
				for(j = 0; j < 4; ++j)
					for(i = 0; i < 3; ++i)
						vpos[v].xyz[i] += tdata->allverts[4*tet+j].xyz[i]*tdata->candidates[inds[v]][4*tet+j];
			for(i = 0; i < 3; ++i) {
				v0.xyz[i] = vpos[1].xyz[i] - vpos[0].xyz[i]; 
				v1.xyz[i] = vpos[2].xyz[i] - vpos[0].xyz[i]; 
				v2.xyz[i] = vpos[3].xyz[i] - vpos[0].xyz[i]; 
			}

			// compute the volume
			// TODO: Understand how to recover the proper sign without the fabs...
			mass = fabs(-v2.x*v1.y*v0.z + v1.x*v2.y*v0.z + v2.x*v0.y*v1.z
					   - v0.x*v2.y*v1.z - v1.x*v0.y*v2.z + v0.x*v1.y*v2.z); 

			mass *= 0.1666666666666666;
			vtot[tet] += mass; 

			for(v = 0; v < 4; ++v)
				printf("  tet %d, pos = %f %f %f\n", tet, vpos[v].x, vpos[v].y, vpos[v].z);
			printf("   vol[tet %d] = %.5e\n", tet, mass);


			// TODO: 4D reduction process for barycentric coordinates??

		}
	}
	printf("Vtot[0] = %.5e, vtot[1] = %.5e\n", vtot[0], vtot[1]);


}

