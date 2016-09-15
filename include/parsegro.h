#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "atom.h"

#ifndef PARSEGRO_H
#define PARSEGRO_H

struct atom * parsegro (char **argv, int * NATOMS, float * BOX, struct atom *ALL_ATOMS)
{
	FILE *fp = fopen (argv[1], "r");

	if(fp==0) {printf("\n\tinput file error. Please check %s.\n", argv[1]); exit(1);}

	char buffer1[100];        
	fgets(buffer1, 100, fp);
	
	fscanf(fp, "%d", NATOMS);

	ALL_ATOMS = (struct atom *) realloc (ALL_ATOMS, *NATOMS * sizeof(struct atom));

	int i;
	for(i=0; i < *NATOMS; i++)		// this loop has to be sequential
	{
		fscanf(fp, "%d%c%c%c%c%c%c%c%c%c%c %d %lf %lf %lf %lf %lf %lf",
		&ALL_ATOMS[i].resid,
		&ALL_ATOMS[i].resname[0],&ALL_ATOMS[i].resname[1],&ALL_ATOMS[i].resname[2],&ALL_ATOMS[i].resname[3],
		&ALL_ATOMS[i].resname[4], &ALL_ATOMS[i].atom_name[0],&ALL_ATOMS[i].atom_name[1],
		&ALL_ATOMS[i].atom_name[2],&ALL_ATOMS[i].atom_name[3],&ALL_ATOMS[i].atom_name[4],
		&ALL_ATOMS[i].atom_index, 
		&ALL_ATOMS[i].position[0], &ALL_ATOMS[i].position[1], &ALL_ATOMS[i].position[2], 
		&ALL_ATOMS[i].velocity[0], &ALL_ATOMS[i].velocity[1], &ALL_ATOMS[i].velocity[2]);
		ALL_ATOMS[i].force[0]=0.0; ALL_ATOMS[i].force[1]=0.0; ALL_ATOMS[i].force[2]=0.0;
	}
	
	fscanf(fp, "%f %f %f", &BOX[0], &BOX[1], &BOX[2]);

	fclose(fp);
	fp=NULL;

	return ALL_ATOMS;

}

void printgro (char **argv, int NATOMS, float * BOX, struct atom *ALL_ATOMS, int steps)
{
	FILE *gp = fopen (argv[2], "a+");
	if(gp==0) printf("\n\toutput file error. Please check %s.\n", argv[2]);

	fprintf(gp, "t=%08d (fs)\t  AUTHOR  Changhui Xu @ University of Iowa 2014\n", steps);
	fprintf(gp, "%d\n", NATOMS);

	int i;
	for(i=0; i < NATOMS; i++)		// this loop has to be sequential
	{
		fprintf(gp, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", 
		ALL_ATOMS[i].resid,
		ALL_ATOMS[i].resname, 
		ALL_ATOMS[i].atom_name,
		ALL_ATOMS[i].atom_index, 
		ALL_ATOMS[i].position[0], ALL_ATOMS[i].position[1], ALL_ATOMS[i].position[2], 
		ALL_ATOMS[i].velocity[0], ALL_ATOMS[i].velocity[1], ALL_ATOMS[i].velocity[2]);
	}
	fprintf(gp, "%10.5f%10.5f%10.5f\n", BOX[0], BOX[1], BOX[2]);

	fclose(gp);
	gp=NULL;
}

#endif

