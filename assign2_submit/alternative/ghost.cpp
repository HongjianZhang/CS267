#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "mpi_particles.h"

//
// Constant Definitions
//

#define GHOST_LENGTH (cutoff)
extern int rank;
//
// Local buffers
//

int* ghost_ids;

// For sending ghost packets
int new_ghost_packet_len[8];
int del_ghost_packet_len[8];
int mod_ghost_packet_len[8];
particle_t* new_ghost_packet[8];
particle_t* del_ghost_packet[8];
particle_t* mod_ghost_packet[8];
MPI_Request mpi_ghost_requests[24];

// For receiving ghost packets
particle_t* new_received_ghosts;
particle_t* mod_received_ghosts;
particle_t* del_received_ghosts;

//
// Ghost functions
//

void setup_ghost_structure(int max_particles)
{
	ghost_ids = (int*) malloc(max_particles * sizeof(int));

	for(int i = 0; i < 8; ++i)
	{
		new_ghost_packet[i] = (particle_t *) malloc(max_particles * sizeof(particle_t));
		del_ghost_packet[i] = (particle_t *) malloc(max_particles * sizeof(particle_t));
		mod_ghost_packet[i] = (particle_t *) malloc(max_particles * sizeof(particle_t));

        new_ghost_packet_len[i] = 0;
        del_ghost_packet_len[i] = 0;
        mod_ghost_packet_len[i] = 0;
	}
	
	new_received_ghosts = (particle_t *) malloc(max_particles * sizeof(particle_t));
	mod_received_ghosts = (particle_t *) malloc(max_particles * sizeof(particle_t));
	del_received_ghosts = (particle_t *) malloc(max_particles * sizeof(particle_t));
}

void clean_ghost_structure()
{
    free(ghost_ids);

	for(int i = 0; i < 8; ++i)
	{
		free(new_ghost_packet[i]);
		free(del_ghost_packet[i]);
		free(mod_ghost_packet[i]);
	}
	
	free(new_received_ghosts);
	free(mod_received_ghosts);
	free(del_received_ghosts);
}


void prepare_initial_ghost_packets(partition_t* part, int* local_ids, int nlocal, double left_x, double right_x, double bottom_y, double top_y, int neighbors[])
{
    for(int i = 0; i < nlocal; ++i)
    {
        particle_t particle = *get_particle(part, local_ids[i]);

        if(particle.x <= (left_x + GHOST_LENGTH)) // check if in W, SW, or NW ghost zone by x
        {
            if((neighbors[p_w] != -1)) // Add to packet if W neighbor exists
            {
                new_ghost_packet[p_w][ new_ghost_packet_len[p_w] ] = particle;
                new_ghost_packet_len[p_w] += 1;
            }

            if(particle.y <= (bottom_y + GHOST_LENGTH)) // Add to packet if SW neighbor exists and y bounded
            {
                if((neighbors[p_sw] != -1))
                {
                    new_ghost_packet[p_sw][ new_ghost_packet_len[p_sw] ] = particle;
                    new_ghost_packet_len[p_sw] += 1;
                }
            }

            else if (particle.y >= (top_y - GHOST_LENGTH)) // Add to packet if NW neighbor exists and y bounded
            {
                if((neighbors[p_nw] != -1))
                {
                    new_ghost_packet[p_nw][ new_ghost_packet_len[p_nw] ] = particle;
                    new_ghost_packet_len[p_nw] += 1;
                }
            }
        }

        else if(particle.x >= (right_x - GHOST_LENGTH)) // check if in E, SE, or NE ghost zone by x
        {
            if((neighbors[p_e] != -1)) // Add to packet if E neighbor exists
            {
                new_ghost_packet[p_e][ new_ghost_packet_len[p_e] ] = particle;
                new_ghost_packet_len[p_e] += 1;
            }

            if(particle.y <= (bottom_y + GHOST_LENGTH)) // Add to packet if SE neighbor exists and y bounded
            {
                if((neighbors[p_se] != -1))
                {
                    new_ghost_packet[p_se][ new_ghost_packet_len[p_se] ] = particle;
                    new_ghost_packet_len[p_se] += 1;
                }
            }

            else if (particle.y >= (top_y - GHOST_LENGTH)) // Add to packet if NE neighbor exists and y bounded
            {
                if((neighbors[p_ne] != -1))
                {
                    new_ghost_packet[p_ne][ new_ghost_packet_len[p_ne] ] = particle;
                    new_ghost_packet_len[p_ne] += 1;
                }
            }
        }

        if(particle.y <= (bottom_y + GHOST_LENGTH)) // check if in S ghost zone by y (SW, SE already handled)
        {
            if((neighbors[p_s] != -1)) // Add to packet if S neighbor exists
            {
                new_ghost_packet[p_s][ new_ghost_packet_len[p_s] ] = particle;
                new_ghost_packet_len[p_s] += 1;
            }
        }

        else if(particle.y >= (top_y - GHOST_LENGTH)) // check if in N ghost zone by y (NW, NE already handled)
        {
            if((neighbors[p_n] != -1)) // Add to packet if S neighbor exists
            {
                new_ghost_packet[p_n][ new_ghost_packet_len[p_n] ] = particle;
                new_ghost_packet_len[p_n] += 1;
            }
        }
    }
}


//
// Determines if a particle was a ghost last cycle, and whether it
// is a ghost this cycle. If either is true, the particle is added
// to the relevant ghost list
//
void update_ghost(particle_t &particle, double oldx, double oldy, double left_x, double right_x, double bottom_y, double top_y)
{
    int old_ghost_dir[8];
    int new_ghost_dir[8];
    int old_ghost_dir_cnt = 0;
    int new_ghost_dir_cnt = 0;

    // Determine whether a new ghost

	if(particle.x <= (left_x + GHOST_LENGTH)) // check if in W, SW, or NW ghost zone by x
	{
        new_ghost_dir[new_ghost_dir_cnt] = P_W;
        new_ghost_dir_cnt++;

		if(particle.y <= (bottom_y + GHOST_LENGTH)) // Add to packet if SW neighbor exists and y bounded
		{
            new_ghost_dir[new_ghost_dir_cnt] = P_SW;
            new_ghost_dir_cnt++;
		}
		else if (particle.y >= (top_y - GHOST_LENGTH)) // Add to packet if NW neighbor exists and y bounded
		{
            new_ghost_dir[new_ghost_dir_cnt] = P_NW;
            new_ghost_dir_cnt++;
		}
	}
	else if(particle.x >= (right_x - GHOST_LENGTH)) // check if in E, SE, or NE ghost zone by x
	{
        new_ghost_dir[new_ghost_dir_cnt] = P_E;
        new_ghost_dir_cnt++;

		if(particle.y <= (bottom_y + GHOST_LENGTH)) // Add to packet if SE neighbor exists and y bounded
		{
            new_ghost_dir[new_ghost_dir_cnt] = P_SE;
            new_ghost_dir_cnt++;
		}
		else if (particle.y >= (top_y - GHOST_LENGTH)) // Add to packet if NE neighbor exists and y bounded
		{
            new_ghost_dir[new_ghost_dir_cnt] = P_NE;
            new_ghost_dir_cnt++;
		}
	}
	
	if(particle.y <= (bottom_y + GHOST_LENGTH)) // check if in S ghost zone by y (SW, SE already handled)
	{
        new_ghost_dir[new_ghost_dir_cnt] = P_S;
        new_ghost_dir_cnt++;
	}
	else if(particle.y >= (top_y - GHOST_LENGTH)) // check if in N ghost zone by y (NW, NE already handled)
	{
        new_ghost_dir[new_ghost_dir_cnt] = P_N;
        new_ghost_dir_cnt++;
	}

    // Determine whether was a ghost last cycle

	if(oldx <= (left_x + GHOST_LENGTH)) // check if in W, SW, or NW ghost zone by x
	{
        old_ghost_dir[old_ghost_dir_cnt] = P_W;
        old_ghost_dir_cnt++;

		if(oldy <= (bottom_y + GHOST_LENGTH)) // Add to packet if SW neighbor exists and y bounded
		{
            old_ghost_dir[old_ghost_dir_cnt] = P_SW;
            old_ghost_dir_cnt++;
		}
		else if (oldy >= (top_y - GHOST_LENGTH)) // Add to packet if NW neighbor exists and y bounded
		{
            old_ghost_dir[old_ghost_dir_cnt] = P_NW;
            old_ghost_dir_cnt++;
		}
	}
	else if(oldx >= (right_x - GHOST_LENGTH)) // check if in E, SE, or NE ghost zone by x
	{
        old_ghost_dir[old_ghost_dir_cnt] = P_E;
        old_ghost_dir_cnt++;

		if(oldy <= (bottom_y + GHOST_LENGTH)) // Add to packet if SE neighbor exists and y bounded
		{
            old_ghost_dir[old_ghost_dir_cnt] = P_SE;
            old_ghost_dir_cnt++;
		}
		else if (oldy >= (top_y - GHOST_LENGTH)) // Add to packet if NE neighbor exists and y bounded
		{
            old_ghost_dir[old_ghost_dir_cnt] = P_NE;
            old_ghost_dir_cnt++;
		}
	}
	
	if(oldy <= (bottom_y + GHOST_LENGTH)) // check if in S ghost zone by y (SW, SE already handled)
	{
        old_ghost_dir[old_ghost_dir_cnt] = P_S;
        old_ghost_dir_cnt++;
	}
	else if(oldy >= (top_y - GHOST_LENGTH)) // check if in N ghost zone by y (NW, NE already handled)
	{
        old_ghost_dir[old_ghost_dir_cnt] = P_N;
        old_ghost_dir_cnt++;
	}

    // If particle is about to leave this processor, we need to add it to the relevant remove lists.
    // We do this by making new_ghost_dir_cnt = 0.
    if( (particle.x < left_x) || (particle.x > right_x) || (particle.y < bottom_y) || (particle.y > top_y) )
    {
        new_ghost_dir_cnt = 0;
    }

    // Iterate through the new ghost regions to which this particle belongs
    for(int i = 0; i < new_ghost_dir_cnt; i++)
    {
        int match = 0;

        // Iterate through the old ghost regions to which this particle belonged last cycle
        for(int j = 0; j < old_ghost_dir_cnt; j++)
        {
            if(new_ghost_dir[i] == old_ghost_dir[j])
            {
                match = 1;
                break;
            }
        }

        // If the ghost belonged to this neighbor region last cycle
        if(match)
        {
            int ghost_dir = new_ghost_dir[i];
            mod_ghost_packet[ghost_dir][ mod_ghost_packet_len[ghost_dir] ] = particle;
            mod_ghost_packet_len[ghost_dir] += 1;
        }
        // Else a new ghost
        else
        {
            int ghost_dir = new_ghost_dir[i];
            new_ghost_packet[ghost_dir][ new_ghost_packet_len[ghost_dir] ] = particle;
            new_ghost_packet_len[ghost_dir] += 1;
        }
    }

    // Iterate through the old ghost regions to find deleted ghosts
    for(int i = 0; i < old_ghost_dir_cnt; i++)
    {
        int match = 0;

        // Iterate through the new ghost regions
        for(int j = 0; j < new_ghost_dir_cnt; j++)
        {
            if(old_ghost_dir[i] == new_ghost_dir[j])
            {
                match = 1;
                break;
            }
        }

        // If the ghost is no longer a ghost, we need to add it to the deleted list
        if(!match)
        {
            int ghost_dir = old_ghost_dir[i];
            del_ghost_packet[ghost_dir][ del_ghost_packet_len[ghost_dir] ] = particle;
            del_ghost_packet_len[ghost_dir] += 1;
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
			MPI_Isend(new_ghost_packet[i], new_ghost_packet_len[i], PARTICLE, neighbors[i], NEW_GHOST_TAG, MPI_COMM_WORLD, &(mpi_ghost_requests[rc++]));
			MPI_Isend(mod_ghost_packet[i], mod_ghost_packet_len[i], PARTICLE, neighbors[i], MOD_GHOST_TAG, MPI_COMM_WORLD, &(mpi_ghost_requests[rc++]));
			MPI_Isend(del_ghost_packet[i], del_ghost_packet_len[i], PARTICLE, neighbors[i], DEL_GHOST_TAG, MPI_COMM_WORLD, &(mpi_ghost_requests[rc++]));
		}
	}

    // Clear ghost counters
    for(int i = 0; i < 8; i++)
    {
        new_ghost_packet_len[i] = 0;
        del_ghost_packet_len[i] = 0;
        mod_ghost_packet_len[i] = 0;
    }
}


void receive_ghost_packets(partition_t* part, int* ghost_ids, int* nghosts, int* neighbors, int num_neighbors, int buf_size)
{
    MPI_Status status;
	
	// Receive new ghosts
	int total_new_received = 0;
	int total_del_received = 0;
	int total_mod_received = 0;

    for(int i = 0; i < 8; i++)
    {
		int num_particles_rcvd = 0;
		
        // If no neighbor in this direction, skip over it
        if(neighbors[i] == NONE) continue;

        // Perform blocking read from neighbor for new ghost packets
        MPI_Recv (new_received_ghosts+total_new_received, (buf_size-total_new_received), PARTICLE, neighbors[i], NEW_GHOST_TAG, MPI_COMM_WORLD, &status); 
        MPI_Get_count(&status, PARTICLE, &num_particles_rcvd);
		total_new_received += num_particles_rcvd;

        // Perform blocking read from neighbor for modified ghost packets
        MPI_Recv (mod_received_ghosts+total_mod_received, (buf_size-total_mod_received), PARTICLE, neighbors[i], MOD_GHOST_TAG, MPI_COMM_WORLD, &status); 
        MPI_Get_count(&status, PARTICLE, &num_particles_rcvd);
		total_mod_received += num_particles_rcvd;

        // Perform blocking read from neighbor for deleted ghost packets
        MPI_Recv (del_received_ghosts+total_del_received, (buf_size-total_del_received), PARTICLE, neighbors[i], DEL_GHOST_TAG, MPI_COMM_WORLD, &status); 
        MPI_Get_count(&status, PARTICLE, &num_particles_rcvd);
		total_del_received += num_particles_rcvd;
    }
	
	// Remove ghosts that are no longer ghosts
	for(int i = 0; i < total_del_received; ++i)
	{
        int global_id = del_received_ghosts[i].globalID;
		remove_particle(part, ghost_ids[global_id]);
	}

	// Add new ghosts to system from packet
	for(int i = 0; i < total_new_received; ++i)
	{
		int local_id = add_particle(part, new_received_ghosts[i]);
		set_ghost(part, local_id, GHOST_TRUE);
        ghost_ids[new_received_ghosts[i].globalID] = local_id;
	}

	// Modify existing ghosts
	for(int i = 0; i < total_mod_received; ++i)
	{
        particle_t ghost = mod_received_ghosts[i];
        set_state(part, ghost_ids[ghost.globalID], ghost.x, ghost.y, ghost.vx, ghost.vy); // global ID, not local ID
	}

    // *nghosts = total_received;

    // Make sure that all previous ghost messages have been sent, as we need to reuse the buffers
    MPI_Waitall(num_neighbors * 3, mpi_ghost_requests, MPI_STATUSES_IGNORE);
}
