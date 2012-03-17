#ifndef __MPI_MB_H
#define __MPI_MB_H

#include <mpi.h>
#include "common.h"
#include "microblock.h"
#include "plist.h"
#include "ppile.h"

#define GHOST_LENGTH cutoff

#define EMIGRANT_TAG    1
#define GHOST_TAG       2

#define GHOST_TRUE	1
#define GHOST_FALSE	0

const double micro_length = 2*cutoff;

extern MPI_Datatype PARTICLE;

typedef struct {
	int rank;
	int proc_x, proc_y;
	int neighbors[8];
	int num_neighbors;
	
	double left_x, right_x, bottom_y, top_y;
	int num_micro_x, num_micro_y;
	double mfactor_x, mfactor_y; // Multiply by these to get approx microblock loc from coords
	
	microblock* nw_ghostblock;
	microblock* ne_ghostblock;
	microblock* sw_ghostblock;
	microblock* se_ghostblock;
	microblock* n_ghostblocks;
	microblock* s_ghostblocks;
	microblock* e_ghostblocks;
	microblock* w_ghostblocks;
} mpi_cell;

//
// Functions to initialize each processor's set of particles at the beginning of the simulation, and to
// save all the particles at the end of the simulation.
//
mpi_cell* create_mpi_cell(double sim_size, int num_proc_x, int num_proc_y, int rank);
void setup_microblocks(microblock* microblocks, mpi_cell* mycell);
void distribute_particles(microblock* microblocks, mpi_cell* mycell, plist* local, particle_t* particles, int n);
void prepare_save(int rank, int n_proc, plist* local, particle_t* particles, int n);
bool compare_particles(particle_t left, particle_t right); // check if left_id < right_id

void process_particles(microblock* microblocks, int num_mb);

//
// Functions to handle ghost buffers
//
void setup_ghost_structure(int max_particles);
void clean_ghost_structure();

void prepare_ghost_packets(mpi_cell* mycell, microblock* microblocks);
void send_ghost_packets(mpi_cell* mycell);
void receive_ghost_packets(mpi_cell* mycell, ppile* ghost, int max_particles);

//
// Functions to handle emigrant buffers
//
void init_emigrant_buf(int n);
void free_emigrant_buf();

void prepare_emigrants(mpi_cell* mycell, microblock* microblocks, plist* local);
void send_emigrants(mpi_cell* mycell);
void receive_immigrants(mpi_cell* mycell, microblock* microblocks, plist* local, int max_particles);

#endif
