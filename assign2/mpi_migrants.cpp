#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common.h"
#include "mpi_mb.h"
#include "plist.h"

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
// emigrant buf. 'local_n' is updated at the end with the number of particles currently remaining in the system.
// 'neighbors' is indexed by the direction values (e.g. N=0, S=1, etc) and gives the corresponding neighbor's rank.
// If a neighbor is not present, its direction is indicated as -1.
//
void prepare_emigrants(mpi_cell* mycell, microblock* microblocks, plist* local)
{
	int message_dir;

    // Initialize all emigrant counts to zero
    for(int i = 0; i < 8; i++) emigrant_cnt[i] = 0;

    // Loop through all the particles in each microblock, checking if they have left the bounds of this processor or their microblock
	for(int mb = 0; mb < (mycell->num_micro_x*mycell->num_micro_y); ++mb)
	{
		for(int i = 0; i < microblocks[mb].num_particles; ++i)
		{
			particle_t* migrant = microblocks[mb].particles[i];
			if(migrant->x < microblocks[mb].left_x   ||\
			   migrant->x > microblocks[mb].right_x  ||\
			   migrant->y < microblocks[mb].bottom_y ||\
			   migrant->y > microblocks[mb].top_y)
			{
				// Remove particle from cell
				mb_rm_particle(microblocks + mb, i);
				--i;
				
				int new_loc;
				if(migrant->x < mycell->left_x)
				{
					if(migrant->y < mycell->bottom_y)	new_loc = p_sw
					else if(migrant->y > mycell->top_y) new_loc = p_nw
					else new_loc = p_w
				}
				else if(migrant->x > mycell->right_x)
				{
					if(migrant->y < mycell->bottom_y)	new_loc = p_se
					else if(migrant->y > mycell->top_y) new_loc = p_ne
					else new_loc = p_e
				}
				else
				{
					if(migrant->y < mycell->bottom_y)	new_loc = p_s
					else if(migrant->y > mycell->top_y) new_loc = p_n
					else new_loc = -1
				}
				
				if(new_loc == -1) // Relocate particle in cell
				{
					int mb_x, mb_y;
					mb_x = migrant->x * mycell->mfactor_x;
					mb_y = migrant->y * mycell->mfactor_y;

					mb_x = min(max(0,mb_x), mycell->num_micro_x);
					mb_y = min(max(0,mb_y), mycell->num_micro_y);
					
					mb_add_particle(microblocks + mb_y*mycell->num_micro_x + mb_x, migrant);
				}
				else // Send to another cell
				{
					emigrant_buf[new_loc][emigrant_cnt[new_loc]] = *migrant;
					emigrant_cnt[new_loc] += 1;
					
					rm_particle(local, migrant); //Actually remove this particle from our local particle list
				}
			}
		}
	}
	
}


//
// Actually sends the immigrants from this processor to neighboring processors.
// The emigrants are determined in the 'prepare_emigrants' function. This function sends the 
// contents of the global emigrant_buf arrays to the neighbors.
// 'neighbors' is indexed by the direction values (e.g. N=0, S=1, etc) and gives the corresponding neighbor's rank.
// If a neighbor is not present, its direction is indicated as -1.
//
void send_emigrants(mpi_cell* mycell)
{
    int num_requests = 0;

	for(int i = 0; i < 8; ++i)
	{
		if(mycell->neighbors[i] != NONE)
		{
			MPI_Isend ((void*)(emigrant_buf[i]), emigrant_cnt[i], PARTICLE, mycell->neighbors[i], EMIGRANT_TAG, MPI_COMM_WORLD, &(mpi_em_requests[num_requests]));
			num_requests++;
		}
	}
}


//
// Receives immigrant particles from neighboring processors, and adds
// them to the list of local particles.
// 'neighbors' is indexed by the direction values (e.g. N=0, S=1, etc) and gives the corresponding neighbor's rank.
// If a neighbor is not present, its direction is indicated as -1.
// 'local_n' is the number of valid particles on this processor, and is updated once all particles have been received.
// 'buf_size' is in terms of number of particles, not bytes
//
void receive_immigrants(mpi_cell* mycell, microblock* microblocks, plist* local, int max_particles)
{
	MPI_Status status;
	int num_particles_rcvd = 0;

	for(int i = 0; i < 8; i++)
	{
		// If no neighbor in this direction, skip over it
		if(mycell->neighbors[i] == -1) continue;

		// Perform blocking read from neighbor
		MPI_Recv ((void*)(immigrant_buf), max_particles, PARTICLE, mycell->neighbors[i], EMIGRANT_TAG, MPI_COMM_WORLD, &status); 

		MPI_Get_count(&status, PARTICLE, &num_particles_rcvd);

		// If the neighbor sent particles, add them to the local particle list
		// Also update the number of particles stored locally on this processor
		for(int j = 0; j < num_particles_rcvd; j++)
		{
			particle_t* target = add_particle(local, immigrant_buf[j]);
			
			int mb_x, mb_y;
			mb_x = target->x * mycell->mfactor_x;
			mb_y = target->y * mycell->mfactor_y;

			mb_x = min(max(0,mb_x), mycell->num_micro_x);
			mb_y = min(max(0,mb_y), mycell->num_micro_y);
			
			mb_add_particle(microblocks + mb_y*mycell->num_micro_x + mb_x, target);
		}
	}

	// Make sure that all previous emigrant messages have been sent, as we need to reuse the buffers
	MPI_Waitall(mycell->num_neighbors, mpi_em_requests, MPI_STATUSES_IGNORE);
}
