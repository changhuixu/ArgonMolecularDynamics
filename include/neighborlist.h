#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H

#define Rm 1.50
double Rm2 = (Rm*Rm);

int ** neighborlist (int NATOMS, float * BOX, int ** nlist, struct atom * ALL_ATOMS)
{
	int i, j, k;
	double Lx=0.5*BOX[0], Ly=0.5*BOX[1], Lz=0.5*BOX[2];
	double dx, dy, dz;
	double r2;			// square of distance between two atoms

   #ifdef _OPENMP
   #pragma omp parallel shared(nlist, ALL_ATOMS) private(i,j,k,r2,dx,dy,dz)
   #pragma omp for
   #endif 
	for(i=0; i<NATOMS; i++)
	{
		nlist[i][0] = 0;
		for(j=i+1; j<NATOMS; j++)
		{
			dx = ALL_ATOMS[i].position[0] - ALL_ATOMS[j].position[0];
			if(dx > Lx) dx -= BOX[0]; else if (dx < -Lx) dx += BOX[0];
			dy = ALL_ATOMS[i].position[1] - ALL_ATOMS[j].position[1];
			if(dy > Ly) dy -= BOX[1]; else if (dy < -Ly) dy += BOX[1];
			dz = ALL_ATOMS[i].position[2] - ALL_ATOMS[j].position[2];
			if(dz > Lz) dz -= BOX[2]; else if (dz < -Lz) dz += BOX[2];

			r2 = dx*dx + dy*dy + dz*dz;	// r2 = r*r faster than r = sqrt()
			if(r2 < Rm2)
			{
				k = ++nlist[i][0];	// save the number of neighbor atoms
				nlist[i][k] = j;	// save atom index to neighbor list
			}
		}
	}
	return nlist;
}


#endif
