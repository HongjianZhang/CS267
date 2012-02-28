#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common.h"
#include "mpi_particles.h"

//
// Buffers for storing emigrant/immigrant particles
//
particle_t **emigrant_buf;
particle_t *immigrant_buf;
int *emigrant_cnt;


//
// Stores the state of asynchronous MPI sends
//
MPI_Request mpi_em_requests[8];

//
// Allocates memory for the emigrant/immigrant buffers
// 'buf_size' is in terms of number of particles, not bytes
//
void init_emigrant_buf(int buf_size)
{
    immigrant_buf = (particle_t*)malloc(buf_size * sizeof(particle_t));

    emigrant_buf = (particle_t**)malloc(8 * sizeof(particle_t*));
    emigrant_cnt = (int*)malloc(buf_size * sizeof(int));

    for(int i = 0; i < 8; i++)
    {
        emigrant_buf[i] = (particle_t*)malloc(buf_size * sizeof(particle_t));
    }
}


//
// Frees memory that was allocated to the emigrant/immigrant buffers
//
void free_emigrant_buf()
{
    for(int i = 0; i < 8; i++)
    {
        free(emigrant_buf[i]);
    }

    free(immigrant_buf);
    free(emigrant_buf);
    free(emigrant_cnt);
}


//
// Determines which of the particles on this processor are emigrants.
// The code searches through all the particles in 'particles', comparing their position to provided x,y processor boundaries.
// If a particle is outside the processor boundaries, it is REMOVED from the 'particles' array and added to the global
// emigrant buf. 'num_particles' is updated at the end with the number of particles currently remaining in the system.
// 'neighbors' is indexed by the direction values (e.g. N=0, S=1, etc) and gives the corresponding neighbor's rank.
// If a neighbor is not present, its direction is indicated as -1.
//
void prepare_emigrants(particle_t* particles, char* p_valid, int* num_particles, double left_x, double right_x, double bottom_y, double top_y, int* neighbors)
{
    int num_particles_checked = 0;
    int num_particles_removed = 0;
    int exit_dir;

    // Initialize all emigrant counts to zero
    for(int i = 0; i < 8; i++)
        emigrant_cnt[i] = 0;

    // Loop through all the particles, checking if they have left the bounds of this processor
    for(int i = 0; num_particles_checked < num_particles; i++)
    {
        if(p_valid[i] == INVALID)
            continue;

        exit_dir = -1;

        // If moved north-west
        if( (particles[i].y > top_y) && (particles[i].x < left_x) && (neighbors[P_NW] != -1))
        {
            exit_dir = P_NW;
        }

        // Else if moved north-east
        else if( (particles[i].y > top_y) && (particles[i].x > right_x) && (neighbors[P_NE] != -1))
        {
            exit_dir = P_NE;
        }

        // Else if moved north
        else if( (particles[i].y > top_y) && (neighbors[P_N] != -1))
        {
            exit_dir = P_N;
        }

        // Else if moved south-west
        else if( (particles[i].y < bottom_y) && (particles[i].x < left_x) && (neighbors[P_SW] != -1))
        {
            exit_dir = P_SW;
        }

        // Else if moved south-east
        else if( (particles[i].y < bottom_y) && (particles[i].x > right_x) && (neighbors[P_SE] != -1))
        {
            exit_dir = P_SE;
        }

        // Else if moved south
        else if( (particles[i].y < bottom_y) && (neighbors[P_S] != -1))
        {
            exit_dir = P_S;
        }

        // Else if moved west
        else if( (particles[i].x < left_x) && (neighbors[P_W] != -1))
        {
            exit_dir = P_W;
        }

        // Else if moved east
        else if( (particles[i].x > right_x) && (neighbors[P_E] != -1))
        {
            exit_dir = P_E;
        }

        // If the particle is an emigrant, remove it from the array and place it in the correct buffer
        if(exit_dir != -1)
        {
            emigrant_buf[exit_dir][emigrant_cnt[exit_dir]] = particles[i];
            emigrant_cnt[exit_dir] += 1;
            remove_particle(i, p_valid);
            num_particles_removed++;
        }

        num_particles_checked++;
    }

    // Update the count of active particles on this processor
    (*num_particles) -= num_particles_removed;
}


//
// Actually sends the immigrants from this processor to neighboring processors.
// The emigrants are determined in the 'prepare_emigrants' function. This function sends the 
// contents of the global emigrant_buf arrays to the neighbors.
// 'neighbors' is indexed by the direction values (e.g. N=0, S=1, etc) and gives the corresponding neighbor's rank.
// If a neighbor is not present, its direction is indicated as -1.
//
void send_emigrants(int* neighbors)
{
    int num_requests = 0;

	for(int i = 0; i < 8; ++i)
	{
		if(neighbors[i] != NONE)
		{
            MPI_Isend ((void*)(emigrant_buf[i]), emigrant_cnt[i], PARTICLE, neighbors[i], EMIGRANT_TAG, MPI_COMM_WORLD, &(mpi_em_requests[num_requests]));
            num_requests++;
		}
	}
}


//
// Receives immigrant particles from neighboring processors, and adds
// them to the list of local particles.
// 'neighbors' is indexed by the direction values (e.g. N=0, S=1, etc) and gives the corresponding neighbor's rank.
// If a neighbor is not present, its direction is indicated as -1.
// 'num_particles' is the number of valid particles on this processor, and is updated once all particles have been received.
// 'buf_size' is in terms of number of particles, not bytes
//
void receive_immigrants(int* neighbors, int num_neighbors, particle_t* particles, char* p_valid, int* num_particles, int array_sz, int buf_size)
{
    MPI_Status status;
    int num_particles_rcvd = 0;

    for(int i = 0; i < 8; i++)
    {
        // If no neighbor in this direction, skip over it
        if(neighbors[i] == -1)
            continue;

        // Perform blocking read from neighbor
        MPI_Recv ((void*)(immigrant_buf), buf_size, PARTICLE, neighbors[i], EMIGRANT_TAG, MPI_COMM_WORLD, &status); 

        MPI_Get_count(&status, PARTICLE, &num_particles_rcvd);

        // If the neighbor sent particles, add them to the local particle list
        for(int j = 0; j < num_particles_rcvd; j++)
        {
            if(add_particle(immigrant_buf[j], array_sz, particles, p_valid) == -1)
            {
                printf("Error: insufficient space to add particle to local array\n");
                exit(-1);
            }
        }

        // Update the number of particles on the local processor
        (*num_particles) += num_particles_rcvd;
    }

    // Make sure that all previous emigrant messages have been sent, as we need to reuse the buffers
    MPI_Waitall(num_neighbors, mpi_em_requests, MPI_STATUSES_IGNORE);
}
