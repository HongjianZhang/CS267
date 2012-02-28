#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "mpi_particles.h"

#define GHOST_LENGTH (cutoff*2)

int ghost_packet_length[8];
particle_t* ghost_packet_particles[8];
// In order, positions are SW, S, SE, W, E, NW, N, NE

MPI_Request mpi_ghost_requests[8];

void setup_ghost_structure(int max_particles)
{
	for(int i = 0; i < 8; ++i)
	{
		ghost_packet_particles[i] = (particle_t *) malloc(max_particles * sizeof(particle_t));
	}
}

void clean_ghost_structure()
{
	for(int i = 0; i < 8; ++i)
	{
		free(ghost_packet_particles[i]);
	}
}

void prepare_ghost_packets(particle_t particles[], char p_valid[], int num_particles, double left_x, double  right_x, double bottom_y, double top_y, int neighbors[])
{
	// Reset so will just overwrite old packet data
	ghost_packet_length[p_sw] = 0; ghost_packet_length[p_s ] = 0; ghost_packet_length[p_se] = 0; ghost_packet_length[p_w ] = 0;
	ghost_packet_length[p_e ] = 0; ghost_packet_length[p_nw] = 0; ghost_packet_length[p_n ] = 0; ghost_packet_length[p_ne] = 0;
	
	int seen_particles = 0;
	for(int i = 0; seen_particles < num_particles; ++i)
	{
		if(p_valid[i] == INVALID) continue;
		seen_particles++;
		
		if(particles[i].x <= (left_x + GHOST_LENGTH)) // check if in W, SW, or NW ghost zone by x
		{
			if((neighbors[p_w] != -1)) // Add to packet if W neighbor exists
			{
				ghost_packet_particles[p_w][ghost_packet_length[p_w]] = particles[i];
				++ghost_packet_length[p_w];
			}
			if(particles[i].y <= (bottom_y + GHOST_LENGTH)) // Add to packet if SW neighbor exists and y bounded
			{
				if((neighbors[p_sw] != -1)) 
				{
					ghost_packet_particles[p_sw][ghost_packet_length[p_sw]] = particles[i];
					++ghost_packet_length[p_sw];
				}
			}
			else if (particles[i].y >= (top_y - GHOST_LENGTH)) // Add to packet if NW neighbor exists and y bounded
			{
				if((neighbors[p_nw] != -1)) 
				{
					ghost_packet_particles[p_nw][ghost_packet_length[p_nw]] = particles[i];
					++ghost_packet_length[p_nw];
				}
			}
		}
		else if(particles[i].x >= (right_x - GHOST_LENGTH)) // check if in E, SE, or NE ghost zone by x
		{
			if((neighbors[p_e] != -1)) // Add to packet if E neighbor exists
			{
				ghost_packet_particles[p_e][ghost_packet_length[p_e]] = particles[i];
				++ghost_packet_length[p_e];
			}
			if(particles[i].y <= (bottom_y + GHOST_LENGTH)) // Add to packet if SE neighbor exists and y bounded
			{
				if((neighbors[p_se] != -1)) 
				{
					ghost_packet_particles[p_se][ghost_packet_length[p_se]] = particles[i];
					++ghost_packet_length[p_se];
				}
			}
			else if (particles[i].y >= (top_y - GHOST_LENGTH)) // Add to packet if NE neighbor exists and y bounded
			{
				if((neighbors[p_ne] != -1)) 
				{
					ghost_packet_particles[p_ne][ghost_packet_length[p_ne]] = particles[i];
					++ghost_packet_length[p_ne];
				}
			}
		}
		
		if(particles[i].y <= (bottom_y + GHOST_LENGTH)) // check if in S ghost zone by y (SW, SE already handled)
		{
			if((neighbors[p_s] != -1)) // Add to packet if S neighbor exists
			{
				ghost_packet_particles[p_s][ghost_packet_length[p_s]] = particles[i];
				++ghost_packet_length[p_s];
			}
		}
		else if(particles[i].y >= (top_y - GHOST_LENGTH)) // check if in N ghost zone by y (NW, NE already handled)
		{
			if((neighbors[p_n] != -1)) // Add to packet if S neighbor exists
			{
				ghost_packet_particles[p_n][ghost_packet_length[p_n]] = particles[i];
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
			MPI_Isend(ghost_packet_particles[i], ghost_packet_length[i], PARTICLE, neighbors[i], GHOST_TAG, MPI_COMM_WORLD, &(mpi_ghost_requests[++rc]));
		}
	}
}

void receive_ghost_packets(int* num_ghost_particles, particle_t* ghost_particles, int* neighbors, int num_neighbors, int buf_size)
{
    MPI_Status status;
	
	*num_ghost_particles = 0;

    for(int i = 0; i < 8; i++)
    {
		int num_particles_rcvd = 0;
		
        // If no neighbor in this direction, skip over it
        if(neighbors[i] == NONE) continue;

        // Perform blocking read from neighbor
        MPI_Recv (ghost_particles+(*num_ghost_particles), (buf_size-(*num_ghost_particles)), PARTICLE, neighbors[i], GHOST_TAG, MPI_COMM_WORLD, &status); 

        MPI_Get_count(&status, PARTICLE, &num_particles_rcvd);
		*num_ghost_particles += num_particles_rcvd;
    }

    // Make sure that all previous emigrant messages have been sent, as we need to reuse the buffers
    MPI_Waitall(num_neighbors, mpi_ghost_requests, MPI_STATUSES_IGNORE);
}
