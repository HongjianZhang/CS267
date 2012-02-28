#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common.h"
#include "mpi_particles.h"
#include <vector>

using std::vector;

//
// Buffers for storing emigrant/immigrant particles
//
particle_t **emigrant_buf;
particle_t *immigrant_buf;


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
}


//
// Sends the immigrants from this processor to neighboring processors
// 'emigrants_to_send' contains the indexes of the particles to send, along with the direction (N, S, W, etc)
// of the processor to which that particle must be sent.
// 'neighbors' is indexed by the direction values (e.g. N=0, S=1, etc) and gives the corresponding neighbor's rank.
// If a neighbor is not present, its direction is indicated as -1.
// The particles are marked as "invalid" in the 'p_valid' array once they have been sent.
//
void send_emigrants(particle_t* emigrants_to_send, int* neighbors, particle_t* particles, char* p_valid, int array_sz)
{
    int dest_rank, part_index;

    // Array of number of particles that must be sent to each neighbor
    char emigrant_cnt[8];

    // Send each buffer to the appropriate neighbors. If there are no emigrants,
    // then send a message to that neighbor indicating so.
    // All sends are asyncrhonous
    int num_requests = 0;

    for(int i = 0; i < 8; i++)
    {
        // Only send a message to a neighbor if that neighbor exists
        if(neighbors[i] != -1)
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
// 'buf_size' is in terms of number of particles, not bytes
//
void receive_immigrants(int* neighbors, int num_neighbors, particle_t* particles, char* p_valid, int array_sz, int buf_size)
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
    }

    // Make sure that all previous emigrant messages have been sent, as we need to reuse the buffers
    MPI_Waitall(num_neighbors, mpi_em_requests, MPI_STATUSES_IGNORE);
}
