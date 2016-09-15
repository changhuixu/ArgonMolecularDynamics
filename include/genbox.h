#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "atom.h"

#ifndef GENBOX_H
#define GENBOX_H

#define TEMPERATURE 85

// velocity has units of nm/ps

void genvelocity (int NATOMS, struct atom *ALL_ATOMS)
{
	int i;
	srand(time(NULL));

	double sum_v2=0.0;		// temperature scaling
	double vx, vy, vz;
	vx=0.0; vy=0.0; vz=0.0;		// conserve monentum

   #ifdef _OPENMP
   #pragma omp parallel for
   #endif 
	for(i=0; i<NATOMS; i++)
	{
		ALL_ATOMS[i].velocity[0]=(double) rand()/ (double) RAND_MAX - 0.5;	// range [-0.5, 0.5]
		sum_v2 += ALL_ATOMS[i].velocity[0] * ALL_ATOMS[i].velocity[0];	
		vx += ALL_ATOMS[i].velocity[0];
		ALL_ATOMS[i].velocity[1]=(double) rand()/ (double) RAND_MAX - 0.5; 
		sum_v2 += ALL_ATOMS[i].velocity[1] * ALL_ATOMS[i].velocity[1];	
		vy += ALL_ATOMS[i].velocity[1];
		ALL_ATOMS[i].velocity[2]=(double) rand()/ (double) RAND_MAX - 0.5;
		sum_v2 += ALL_ATOMS[i].velocity[2] * ALL_ATOMS[i].velocity[2];	
		vz += ALL_ATOMS[i].velocity[2];
	}

	vx /= NATOMS; vy /= NATOMS; vz /= NATOMS;
	double KBM=2.49345e-2*TEMPERATURE*(NATOMS-1)/MASS;
	double scale = sqrt(KBM/sum_v2);	

   #ifdef _OPENMP
   #pragma omp parallel for
   #endif 
	for(i=0; i<NATOMS; i++)
	{
		ALL_ATOMS[i].velocity[0] *= scale;
		ALL_ATOMS[i].velocity[1] *= scale;
		ALL_ATOMS[i].velocity[2] *= scale;
		ALL_ATOMS[i].velocity[0] -= vx;
		ALL_ATOMS[i].velocity[1] -= vy;
		ALL_ATOMS[i].velocity[2] -= vz;
	}
	
}

struct atom * genbox (int *NATOMS, float *BOX, struct atom *ALL_ATOMS)
{
	int i;
	srand(time(NULL));

	*NATOMS=1000;
	BOX[0] = BOX[1] = BOX[2] = 0.362 * pow((float)(*NATOMS), 1.0/3.0);

	ALL_ATOMS = (struct atom *) realloc (ALL_ATOMS, *NATOMS * sizeof(struct atom));

   #ifdef _OPENMP
   #pragma omp parallel for
   #endif 
	for(i=0; i<*NATOMS; i++)
	{
		ALL_ATOMS[i].position[0] = BOX[0] * (double) rand()/ (double) (RAND_MAX-1);
		ALL_ATOMS[i].position[1] = BOX[1] * (double) rand()/ (double) (RAND_MAX-1);
		ALL_ATOMS[i].position[2] = BOX[2] * (double) rand()/ (double) (RAND_MAX-1);
		ALL_ATOMS[i].resid = i+1;
		strncpy(ALL_ATOMS[i].resname, "Ar", 2);
		strncpy(ALL_ATOMS[i].atom_name, "Ar", 2);
		ALL_ATOMS[i].atom_index = i+1;
		ALL_ATOMS[i].velocity[0] = 0.0;
		ALL_ATOMS[i].velocity[1] = 0.0;
		ALL_ATOMS[i].velocity[2] = 0.0;
		
	}

	return ALL_ATOMS;
}

#endif

