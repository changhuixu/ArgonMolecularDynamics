#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "./include/atom.h"
#include "./include/parsegro.h"
#include "./include/energy.h"
#include "./include/verlet.h"
#include "./include/stddev.h"
#include "./include/genbox.h"
#include "./include/neighborlist.h"
#include "./include/rdf.h"

// MASS is defined in "atom.h" as 39.948 g/mol
// TEMPERATURE is defined in "genbox.h" as 85 K
// RCUT is defined in "energy.h" as 1.2 nm
// R_verlet is defined in "neighborlist.h" as Rm=1.5 nm
// TIMESTEP is defined in "verlet.h" as 1 fs

int main( int argc, char *argv[])
{

        if(argc!=6)
        {
                printf("\n\tUsage:\n\t./md in.gro out.gro nsteps n_print_steps n_threads\n");
                return 0;
        }

	int i=0;
	int NATOMS=1;
	int nsteps=atoi(argv[3]);
	int n_print_steps=atoi(argv[4]);
	int NUM_THREADS=atoi(argv[5]);
	omp_set_num_threads(NUM_THREADS);
	float BOX[3]={0.0, 0.0, 0.0};

	struct atom *ALL_ATOMS;
	ALL_ATOMS = (struct atom *) malloc (NATOMS * sizeof(struct atom));	// initial allocation of memory

	ALL_ATOMS = parsegro (argv, &NATOMS, BOX, ALL_ATOMS);			// read input starting geometry

	int **nlist = (int **) malloc (NATOMS * sizeof(int*));
	for(i=0; i<NATOMS; i++) nlist[i] = (int *) malloc ((NATOMS/2) * sizeof(int));

	FILE *foutp = fopen("out.dat", "w");
	fprintf(foutp, "####################################################################################\n");
	fprintf(foutp, "#\t\t\tMolecular Dynamics of Liquid Argon\n");
	fprintf(foutp, "#\tOpenMP parallelized C code\tAuthor: Changhui Xu @UIOWA 2014 FALL\n");
	fprintf(foutp, "#\tNumber of processors available = %d\n", omp_get_num_procs());
	fprintf(foutp, "#\tMax number of threads = %d\n", omp_get_max_threads());
	fprintf(foutp, "#\tNumber of threads used = %d\n", NUM_THREADS);
	fprintf(foutp, "#\tTotol number of atoms: %d \t\tTotal dynamics steps: %d\n", NATOMS, nsteps);
	fprintf(foutp, "#\tLJ cut-off: %lf (nm) \t\tVelet nlist cut-off: %lf (nm)\n", RCUT, Rm);
	fprintf(foutp, "####################################################################################\n");

	double start_time = omp_get_wtime();					// start timer and begin simulation

	double Epot, Ekin;
	double * Etot;
	nlist = neighborlist (NATOMS, BOX, nlist, ALL_ATOMS);
	Etot = (double *) malloc ((nsteps+1) * sizeof(double));
	Ekin = Ekinetic (NATOMS, ALL_ATOMS);				// track the energies
	Epot = Epotential (NATOMS, BOX, nlist, ALL_ATOMS);
	Etot[0] = Epot + Ekin;
	fprintf(foutp, "#\n#\tStep    \ttotal energy (eV)\tPotential Energy (eV)\tKinetic Energy (eV)\n");
	fprintf(foutp, "000000000\t%18.8lf\t%20.8lf\t%20.8lf\n", Etot[0], Epot, Ekin);		// output energies
	
	for(i=1; i<=nsteps; i++)						// this loop has to be sequential
	{
		if(i%5 == 0) nlist = neighborlist (NATOMS, BOX, nlist, ALL_ATOMS);  // update neighbor list every 5fs
		Epot = verlet (NATOMS, BOX, nlist, ALL_ATOMS);			// verlet algorithm---> move atoms
		Ekin = Ekinetic (NATOMS, ALL_ATOMS);
		Etot[i] = Epot + Ekin;
		fprintf(foutp, "%09d\t%18.8lf\t%20.8lf\t%20.8lf\n", i, Etot[i], Epot, Ekin);

		if(i%n_print_steps == 0)
		printgro (argv, NATOMS, BOX, ALL_ATOMS, i);			// print atom trajectory
	}

	fprintf(foutp, "####################################################################################\n");
	fprintf(foutp, "#\t\t\tMolecular Dynamics of Liquid Argon\n");
	double Etot_mean = mean (Etot, nsteps+1);
	double Etot_stddev = standard_deviation (Etot, Etot_mean, nsteps+1);
	double Etot_drift = drift (Etot, nsteps+1);
	fprintf(foutp, "#===================================================================================\n");
	fprintf(foutp, "#\taverage total energy (eV)\t\tRMSD (eV)\t\tdrift (eV)\n");
	fprintf(foutp, "#\t\t%.8lf\t\t\t%.8lf\t\t%.8lf\n", Etot_mean, Etot_stddev, Etot_drift);
	fprintf(foutp, "#===================================================================================\n");
	double criteria = 1.036427e-5*nsteps*TIMESTEP*NATOMS;
	if(fabs(Etot_drift) < criteria) fprintf(foutp, "#\tClassical MD simulation in NVE ensamble for %lf ps ... Success!\n", nsteps*TIMESTEP);
	fflush(foutp);
	rdf(argv, NATOMS, BOX, ALL_ATOMS, nsteps/n_print_steps + 1);
	double run_time = omp_get_wtime() - start_time;				// end timer
	fprintf(foutp, "#\tTotal run_time: %lf seconds.\n", run_time);
	fprintf(foutp, "####################################################################################\n");

	for(i=0; i<NATOMS; i++) free(nlist[i]);
	free(nlist);
	free(ALL_ATOMS);
	free(Etot);
	fclose(foutp);
	foutp=NULL;

	return 0;
}
