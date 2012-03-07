#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "mpi_particles.h"

#define GHOST_LENGTH(cutoff)

int ghost_packet_length[8];
particle_t* ghost_packet_particles[8];
// In order, positions are SW, S, SE, W, E, NW, N, NE

particle_t* received_ghosts;

MPI_Request mpi_ghost_requests[8];

void setup_ghost_structure(int max_particles)
{
	for(int i = 0; i < 8; ++i)
	{
		ghost_packet_particles[i] = (particle_t *) malloc(max_particles * sizeof(particle_t));
	}
	
	received_ghosts = (particle_t *) malloc(max_particles * sizeof(particle_t));
}

void clean_ghost_structure()
{
	for(int i = 0; i < 8; ++i)
	{
		free(ghost_packet_particles[i]);
	}
	
	free(received_ghosts);
}

void prepare_ghost_packets(partition_t* part, int* local_ids, int nlocal, double left_x, double  right_x, double bottom_y, double top_y, int neighbors[])
{
	// Reset so will just overwrite old packet data
	ghost_packet_length[p_sw] = 0; ghost_packet_length[p_s ] = 0; ghost_packet_length[p_se] = 0; ghost_packet_length[p_w ] = 0;
	ghost_packet_length[p_e ] = 0; ghost_packet_length[p_nw] = 0; ghost_packet_length[p_n ] = 0; ghost_packet_length[p_ne] = 0;
	
	for(int i = 0; i < nlocal; ++i)
	{
		particle_t particle = *get_particle(part, local_ids[i]);
		
		if(particle.x <= (left_x + GHOST_LENGTH)) // check if in W, SW, or NW ghost zone by x
		{
			if((neighbors[p_w] != -1)) // Add to packet if W neighbor exists
			{
				ghost_packet_particles[p_w][ghost_packet_length[p_w]] = particle;
				++ghost_packet_length[p_w];
			}
			if(particle.y <= (bottom_y + GHOST_LENGTH)) // Add to packet if SW neighbor exists and y bounded
			{
				if((neighbors[p_sw] != -1)) 
				{
					ghost_packet_particles[p_sw][ghost_packet_length[p_sw]] = particle;
					++ghost_packet_length[p_sw];
				}
			}
			else if (particle.y >= (top_y - GHOST_LENGTH)) // Add to packet if NW neighbor exists and y bounded
			{
				if((neighbors[p_nw] != -1)) 
				{
					ghost_packet_particles[p_nw][ghost_packet_length[p_nw]] = particle;
					++ghost_packet_length[p_nw];
				}
			}
		}
		else if(particle.x >= (right_x - GHOST_LENGTH)) // check if in E, SE, or NE ghost zone by x
		{
			if((neighbors[p_e] != -1)) // Add to packet if E neighbor exists
			{
				ghost_packet_particles[p_e][ghost_packet_length[p_e]] = particle;
				++ghost_packet_length[p_e];
			}
			if(particle.y <= (bottom_y + GHOST_LENGTH)) // Add to packet if SE neighbor exists and y bounded
			{
				if((neighbors[p_se] != -1)) 
				{
					ghost_packet_particles[p_se][ghost_packet_length[p_se]] = particle;
					++ghost_packet_length[p_se];
				}
			}
			else if (particle.y >= (top_y - GHOST_LENGTH)) // Add to packet if NE neighbor exists and y bounded
			{
				if((neighbors[p_ne] != -1)) 
				{
					ghost_packet_particles[p_ne][ghost_packet_length[p_ne]] = particle;
					++ghost_packet_length[p_ne];
				}
			}
		}
		
		if(particle.y <= (bottom_y + GHOST_LENGTH)) // check if in S ghost zone by y (SW, SE already handled)
		{
			if((neighbors[p_s] != -1)) // Add to packet if S neighbor exists
			{
				ghost_packet_particles[p_s][ghost_packet_length[p_s]] = particle;
				++ghost_packet_length[p_s];
			}
		}
		else if(particle.y >= (top_y - GHOST_LENGTH)) // check if in N ghost zone by y (NW, NE already handled)
		{
			if((neighbors[p_n] != -1)) // Add to packet if S neighbor exists
			{
				ghost_packet_particles[p_n][ghost_packet_length[p_n]] = particle;
				++ghost_packet_length[p_n];
			}
		}
	}
}

void send_ghost_packets(int neighbors[])
{
	int rc = 0;
	for(int i = 0; i < 8; ++i)
	{
		if(neighbors[i] != NONE)
		{
			MPI_Isend(ghost_packet_particles[i], ghost_packet_length[i], PARTICLE, neighbors[i], GHOST_TAG, MPI_COMM_WORLD, &(mpi_ghost_requests[rc++]));
		}
	}
}

void receive_ghost_packets(partition_t* part, int* ghost_ids, int* nghosts, int* neighbors, int num_neighbors, int buf_size)
{
    MPI_Status status;
	
	// Remove old ghosts
	for(int i = 0; i < *nghosts; ++i)
	{
		remove_particle(part, ghost_ids[i]);
	}
	
	// Receive new ghosts
	int total_received = 0;
    for(int i = 0; i < 8; i++)
    {
		int num_particles_rcvd = 0;
		
        // If no neighbor in this direction, skip over it
        if(neighbors[i] == NONE) continue;

        // Perform blocking read from neighbor
        MPI_Recv (received_ghosts+total_received, (buf_size-total_received), PARTICLE, neighbors[i], GHOST_TAG, MPI_COMM_WORLD, &status); 

        MPI_Get_count(&status, PARTICLE, &num_particles_rcvd);
		total_received += num_particles_rcvd;
    }
	
	// Add new ghosts to system from packet
	for(int i = 0; i < total_received; ++i)
	{
		ghost_ids[i] = add_particle(part, received_ghosts[i]);
		set_ghost(part, ghost_ids[i], GHOST_TRUE);
	}

    *nghosts = total_received;

    // Make sure that all previous ghost messages have been sent, as we need to reuse the buffers
    MPI_Waitall(num_neighbors, mpi_ghost_requests, MPI_STATUSES_IGNORE);
}
