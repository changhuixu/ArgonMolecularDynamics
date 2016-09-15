#include <stdlib.h>
#include <math.h>
#include <omp.h>

#ifndef RDF_H
#define RDF_H

#define bin 0.002
#define cube(i) ((i)*(i)*(i))

void rdf (char **argv, int NATOMS, float *BOX, struct atom * ALL_ATOMS, int nsteps)
{
	FILE *fp = fopen (argv[2], "r");
        if(fp==0) printf("\n\toutput file error. Please check %s.\n", argv[2]);

	int i,j,k;
	int N = 0.5*BOX[0]/bin + 1;
	float *gofr = (float *) malloc(N*sizeof(float));
	for(i=0; i<N; i++) gofr[i]=0.0;
	int nbin;

	double Lx=0.5*BOX[0], Ly=0.5*BOX[1], Lz=0.5*BOX[2];
	double dx, dy, dz;
	double r;

	for(k=0; k<= nsteps && !feof(fp); k++)	// this loop has to be sequential
	{
		char buffer1[100];        
		if(fgets(buffer1, 100, fp) == NULL) break;
		fgets(buffer1, 100, fp);
		
		for(i=0; i < NATOMS; i++)		// this loop has to be sequential
		{
			fscanf(fp, "%d%c%c%c%c%c%c%c%c%c%c %d %lf %lf %lf %lf %lf %lf",
			&ALL_ATOMS[i].resid,
			&ALL_ATOMS[i].resname[0],&ALL_ATOMS[i].resname[1],&ALL_ATOMS[i].resname[2],&ALL_ATOMS[i].resname[3],
			&ALL_ATOMS[i].resname[4], &ALL_ATOMS[i].atom_name[0],&ALL_ATOMS[i].atom_name[1],
			&ALL_ATOMS[i].atom_name[2],&ALL_ATOMS[i].atom_name[3],&ALL_ATOMS[i].atom_name[4],
			&ALL_ATOMS[i].atom_index, 
			&ALL_ATOMS[i].position[0], &ALL_ATOMS[i].position[1], &ALL_ATOMS[i].position[2], 
			&ALL_ATOMS[i].velocity[0], &ALL_ATOMS[i].velocity[1], &ALL_ATOMS[i].velocity[2]);
		}
		fgets(buffer1, 100, fp);

	   #ifdef _OPENMP
	   #pragma omp parallel shared(gofr, ALL_ATOMS) private(i,j,nbin,r,dx,dy,dz)
	   #pragma omp for
	   #endif 
		for(i=0; i<NATOMS; i++)
		{
			for(j=i+1; j<NATOMS; j++)
			{
				dx = ALL_ATOMS[i].position[0] - ALL_ATOMS[j].position[0];
				if(dx > Lx) dx -= BOX[0]; else if (dx < -Lx) dx += BOX[0];
				dy = ALL_ATOMS[i].position[1] - ALL_ATOMS[j].position[1];
				if(dy > Ly) dy -= BOX[1]; else if (dy < -Ly) dy += BOX[1];
				dz = ALL_ATOMS[i].position[2] - ALL_ATOMS[j].position[2];
				if(dz > Lz) dz -= BOX[2]; else if (dz < -Lz) dz += BOX[2];
	
				r = sqrt(dx*dx + dy*dy + dz*dz);
				nbin = r/bin;
				if(nbin<N) gofr[nbin] += 2.0;
			}
		}
	}

	fclose(fp);
	fp=NULL;

	double Volume = BOX[0]*BOX[1]*BOX[2];
	double rho = NATOMS/Volume;
	double rho2= rho*nsteps*NATOMS*4.1867;
	
   #ifdef _OPENMP
   #pragma omp parallel shared(gofr) private(i)
   #pragma omp for
   #endif 
	for(i=0; i<N; i++) gofr[i] /= (cube(i+1)-cube(i))*cube(bin)*rho2;

	FILE *gp = fopen("rdf.xvg", "w");
	for(i=0; i<N; i++) fprintf(gp, "%lf\t%lf\n", i*bin, gofr[i]);

	fclose(gp);
	gp=NULL;
	free(gofr);
}

#endif

