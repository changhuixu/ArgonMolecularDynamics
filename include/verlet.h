#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <math.h>
#include "./energy.h"

#ifndef VERLET_H
#define VERLET_H

#define TIMESTEP 0.001
// define verlet integration timestep as 1fs
double dt2 = (0.5 * 96.451 * TIMESTEP / MASS);
// converted to accelerator

double verlet (int NATOMS, float * BOX, int ** nlist, struct atom *ALL_ATOMS)
{
	int i;
   #ifdef _OPENMP
   #pragma omp parallel for
   #endif 
	for(i=0; i<NATOMS; i++)
	{
		ALL_ATOMS[i].velocity[0] += ALL_ATOMS[i].force[0] * dt2; 
		ALL_ATOMS[i].velocity[1] += ALL_ATOMS[i].force[1] * dt2; 
		ALL_ATOMS[i].velocity[2] += ALL_ATOMS[i].force[2] * dt2;
		ALL_ATOMS[i].position[0] += ALL_ATOMS[i].velocity[0] * TIMESTEP;
		ALL_ATOMS[i].position[1] += ALL_ATOMS[i].velocity[1] * TIMESTEP;
		ALL_ATOMS[i].position[2] += ALL_ATOMS[i].velocity[2] * TIMESTEP;

		if(ALL_ATOMS[i].position[0] > BOX[0])  ALL_ATOMS[i].position[0] -= BOX[0];	//PBC
		else if(ALL_ATOMS[i].position[0] < 0) ALL_ATOMS[i].position[0] += BOX[0];
		if(ALL_ATOMS[i].position[0] > BOX[0] || ALL_ATOMS[i].position[0] < 0)
		        {printf("\natom %d out of box\n", i); exit(1);}             // double check
		if(ALL_ATOMS[i].position[1] > BOX[1])  ALL_ATOMS[i].position[1] -= BOX[1];
		else if(ALL_ATOMS[i].position[1] < 0) ALL_ATOMS[i].position[1] += BOX[1];
		if(ALL_ATOMS[i].position[2] > BOX[2])  ALL_ATOMS[i].position[2] -= BOX[2];
		else if(ALL_ATOMS[i].position[2] < 0) ALL_ATOMS[i].position[2] += BOX[2];

	}

	double E_pot=0.0;
	E_pot = Epotential (NATOMS, BOX, nlist, ALL_ATOMS);

   #ifdef _OPENMP
   #pragma omp parallel for
   #endif 
	for(i=0; i<NATOMS; i++)
	{
		ALL_ATOMS[i].velocity[0] += ALL_ATOMS[i].force[0] * dt2; 
		ALL_ATOMS[i].velocity[1] += ALL_ATOMS[i].force[1] * dt2; 
		ALL_ATOMS[i].velocity[2] += ALL_ATOMS[i].force[2] * dt2;

	}

	return E_pot;
}

#endif

