#ifndef __MPI_PARTICLES_H__
#define __MPI_PARTICLES_H__

#include <mpi.h>

#define EMIGRANT_TAG    1
#define GHOST_TAG       2

#define GHOST_TRUE	1
#define GHOST_FALSE	0

extern MPI_Datatype PARTICLE;

//
// Functions to initialize emigrant buffers
//
void init_emigrant_buf(int n);
void free_emigrant_buf();

//
// Functions to initialize ghost buffers
//
void setup_ghost_structure(int max_particles);
void clean_ghost_structure();

//
// Functions to initialize each processor's set of particles at the beginning of the simulation, and to
// save all the particles at the end of the simulation.
//
int select_particles(partition_t* part, int n, particle_t* particles, int* local_ids, double left_x, double right_x, double bottom_y, double top_y);
void prepare_save(partition_t* part, int rank, int n_proc, int* local_ids, int nlocal, particle_t* particles, int n);
bool compare_particles(particle_t left, particle_t right); // check if left_id < right_id

//
// Functions to send/receive emigrants
//
void prepare_emigrants(partition_t* part, int* local_ids, int* local_n, double left_x, double right_x, double bottom_y, double top_y, int* neighbors);
void send_emigrants(int* neighbors);
void receive_immigrants(int* neighbors, int num_neighbors, partition_t* part, int* local_ids, int* local_n, int buf_size);

//
// Functions to send/receive ghost particles
//
void prepare_ghost_packets(partition_t* part, int* local_ids, int nlocal, double left_x, double  right_x, double bottom_y, double top_y, int neighbors[]);
void send_ghost_packets(int neighbors[]);
void receive_ghost_packets(partition_t* part, int* ghost_ids, int* nghosts, int* neighbors, int num_neighbors, int buf_size);

#endif
