#include <stdlib.h>
#include "atom.h"
#include <time.h>
#include <math.h>
#include <omp.h>
#include "neighborlist.h"

#ifndef ENERGY_H
#define ENERGY_H

// Lennard-Jones model for liquid argon
// epsilon/k_B (K) = 125.7 (some places give epsilon = 1.654 * e-21 J, or 0.0104 eV)
// sigma (nm) = 0.3345
// White(1999), copied from http://www.sklogwiki.org/SklogWiki/index.php/Argon

#define RCUT 1.20
#define sigma 0.3345
#define epsilon 0.0104
double rcut2 = (RCUT * RCUT);
double sigma2 = (sigma * sigma);
double factE = (4.0 * epsilon);
double factF = (48.0 * epsilon);

double Ekinetic (int NATOMS, struct atom *ALL_ATOMS)
{
	int i;

	double Ekin = 0.0;
   #ifdef _OPENMP
   #pragma omp parallel  shared(Ekin, ALL_ATOMS) private(i)
   #pragma omp for reduction(+:Ekin) 
   #endif
	for(i=0; i<NATOMS; i++)
		Ekin += ALL_ATOMS[i].velocity[0] * ALL_ATOMS[i].velocity[0] + \
			ALL_ATOMS[i].velocity[1] * ALL_ATOMS[i].velocity[1] + \
			ALL_ATOMS[i].velocity[2] * ALL_ATOMS[i].velocity[2];	

	Ekin = Ekin * MASS * 0.005184549;	// units (eV)
	return Ekin;
}

double Epotential (int NATOMS, float *BOX, int ** nlist, struct atom *ALL_ATOMS)
{
	int i,j,k;

	double dx, dy, dz;
	double Lx=0.5*BOX[0], Ly=0.5*BOX[1], Lz=0.5*BOX[2];
	double fxi, fyi, fzi;
	double r2;			// square of distance between two atoms
	
	double sr2, sr6, fpr;
	
	double Epot = 0.0;			// clear old value
   #ifdef _OPENMP
   #pragma omp parallel for
   #endif
	for(i=0; i<NATOMS; i++)
	{
		ALL_ATOMS[i].force[0] = 0.0; 
		ALL_ATOMS[i].force[1] = 0.0; 
		ALL_ATOMS[i].force[2] = 0.0; 
	}

   #ifdef _OPENMP
   #pragma omp parallel shared(Epot, nlist, ALL_ATOMS) private(i,j,k,r2,sr2,sr6,fpr,dx,dy,dz,fxi,fyi,fzi)
   #pragma omp for reduction(+:Epot) 
   #endif
	for(i=0; i<NATOMS-1; i++)
	{
		for(k=1; k<=nlist[i][0]; k++)
		{
			j = nlist[i][k];
			dx = ALL_ATOMS[i].position[0] - ALL_ATOMS[j].position[0];
			if(dx > Lx) dx -= BOX[0]; else if (dx < -Lx) dx += BOX[0];
			dy = ALL_ATOMS[i].position[1] - ALL_ATOMS[j].position[1];
			if(dy > Ly) dy -= BOX[1]; else if (dy < -Ly) dy += BOX[1];
			dz = ALL_ATOMS[i].position[2] - ALL_ATOMS[j].position[2];
			if(dz > Lz) dz -= BOX[2]; else if (dz < -Lz) dz += BOX[2];

			r2 = dx*dx + dy*dy + dz*dz;		// r2 = r*r faster than r = sqrt()
			if(r2 < rcut2)
			{
				sr2 = sigma2 / r2;
				sr6 = sr2 * sr2 * sr2;
				fpr = factF * sr6 * (sr6 - 0.5) / r2;

				fxi = fpr * dx;
				fyi = fpr * dy;
				fzi = fpr * dz;
				
				ALL_ATOMS[i].force[0] += fxi; ALL_ATOMS[j].force[0] -=fxi;
				ALL_ATOMS[i].force[1] += fyi; ALL_ATOMS[j].force[1] -=fyi;
				ALL_ATOMS[i].force[2] += fzi; ALL_ATOMS[j].force[2] -=fzi;
				
				Epot += factE * sr6 * (sr6 - 1.0);
			}
		}
	}

	return Epot;
}

#endif

