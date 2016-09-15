#ifndef ATOM_H
#define ATOM_H

/*The number density of liquid Argon at 85K is 20/nm^3*/
#define MASS 39.948

struct atom
{
	int resid;		// residue number (5 positions, integer)
	char resname[5];	// residue name (5 characters)
	char atom_name[5];	// atom name (5 characters)
	int atom_index;		// atom number (5 positions, integer)
	double position[3];	// position (in nm, x y z in 3 columns, 
				//	each 8 positions with 3 decimal places)
	double velocity[3];	// velocity (in nm/ps (or km/s), x y z in 3 columns, 
				//	each 8 positions with 4 decimal places)

				// above is the data structure for gro file
				// below is the data structure used in molecular dynamics

	double force[3];	// forces --------> calculated in "energy.h"
};


#endif

