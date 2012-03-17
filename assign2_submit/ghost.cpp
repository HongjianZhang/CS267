#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "mpi_mb.h"
#include "ppile.h"

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

void add_ghosts(microblock* basket, int packet_id)
{
	for(int i = 0; i < basket->num_particles; ++i)
	{
		ghost_packet_particles[packet_id][ghost_packet_length[packet_id]] = *(basket->particles[i]);
		++ghost_packet_length[packet_id];
	}
}

void prepare_ghost_packets(mpi_cell* mycell, microblock* microblocks)
{
	// Reset so will just overwrite old packet data
	ghost_packet_length[p_sw] = 0; ghost_packet_length[p_s ] = 0; ghost_packet_length[p_se] = 0; ghost_packet_length[p_w ] = 0;
	ghost_packet_length[p_e ] = 0; ghost_packet_length[p_nw] = 0; ghost_packet_length[p_n ] = 0; ghost_packet_length[p_ne] = 0;
/*
	for(int i = 0; i < 8; ++i)
	{
		for(int mb = 0; mb < mycell->num_micro_x*mycell->num_micro_y; ++mb)
		{
			add_ghosts(microblocks + mb, i);
		}
	}
*/
	for(int x = 0; x < mycell->num_micro_x; ++x)
	{
		for(int i = 0; i < 8; ++i)
		{
			add_ghosts(microblocks + x, i);
			add_ghosts(microblocks + (mycell->num_micro_y-1)*mycell->num_micro_x + x, i);
		}
	}
	for(int y = 1; y < mycell->num_micro_y-1; ++y)
	{
		for(int i = 0; i < 8; ++i)
		{
			add_ghosts(microblocks + y*mycell->num_micro_x, i);
			add_ghosts(microblocks + y*mycell->num_micro_x + mycell->num_micro_x - 1, i);
		}
	}
/*	
	// Prepare corners
	if(mycell->neighbors[p_sw]) add_ghosts(&microblocks[0                      *mycell->num_micro_x + 0]                      , p_sw);
	if(mycell->neighbors[p_se]) add_ghosts(&microblocks[0                      *mycell->num_micro_x + (mycell->num_micro_x-1)], p_se);
	if(mycell->neighbors[p_nw]) add_ghosts(&microblocks[(mycell->num_micro_y-1)*mycell->num_micro_x + 0]                      , p_nw);
	if(mycell->neighbors[p_ne]) add_ghosts(&microblocks[(mycell->num_micro_y-1)*mycell->num_micro_x + (mycell->num_micro_x-1)], p_ne);
	
	// Prepare sides
	for(int x = 0; x < mycell->num_micro_x; ++x)
	{
		if(mycell->neighbors[p_n]) add_ghosts(&microblocks[(mycell->num_micro_y-1)*mycell->num_micro_x + x], p_n);
		if(mycell->neighbors[p_s]) add_ghosts(&microblocks[(0                    )*mycell->num_micro_x + x], p_s);
	}
	for(int y = 0; y < mycell->num_micro_y; ++y)
	{
		if(mycell->neighbors[p_e]) add_ghosts(&microblocks[(y)*mycell->num_micro_x + (mycell->num_micro_x-1)], p_e);
		if(mycell->neighbors[p_w]) add_ghosts(&microblocks[(y)*mycell->num_micro_x + 0]                      , p_w);
	}*/
}

void send_ghost_packets(mpi_cell* mycell)
{
	int rc = 0;
	for(int i = 0; i < 8; ++i)
	{
		if(mycell->neighbors[i] != NONE)
		{
			MPI_Isend(ghost_packet_particles[i], ghost_packet_length[i], PARTICLE, mycell->neighbors[i], GHOST_TAG, MPI_COMM_WORLD, &(mpi_ghost_requests[rc++]));
		}
	}
}

void receive_ghost_packets(mpi_cell* mycell, ppile* ghost, int max_particles)
{
    MPI_Status status;
	
	// Remove old ghosts
	clear_particles(ghost);
	
	mb_clear(mycell->nw_ghostblock); mb_clear(mycell->ne_ghostblock);
	mb_clear(mycell->sw_ghostblock); mb_clear(mycell->se_ghostblock);
	
	for(int x = 0; x < mycell->num_micro_x; ++x)
	{
		mb_clear(mycell->n_ghostblocks+x);
		mb_clear(mycell->s_ghostblocks+x);
	}
	for(int y = 0; y < mycell->num_micro_y; ++y)
	{
		mb_clear(mycell->e_ghostblocks+y);
		mb_clear(mycell->w_ghostblocks+y);
	}
	
	// Receive new ghosts
	int total_received = 0;
    for(int i = 0; i < 8; i++)
    {
		int num_particles_rcvd = 0;
		
        // If no neighbor in this direction, skip over it
        if(mycell->neighbors[i] == NONE) continue;

        // Perform blocking read from neighbor
        MPI_Recv (received_ghosts+total_received, (max_particles-total_received), PARTICLE, mycell->neighbors[i], GHOST_TAG, MPI_COMM_WORLD, &status);

        MPI_Get_count(&status, PARTICLE, &num_particles_rcvd);
		total_received += num_particles_rcvd;
    }
	
	// Add new ghosts to system from packet
	for(int i = 0; i < total_received; ++i)
	{
		particle_t* target = add_particle(ghost, received_ghosts[i]);
		
		int mb_x, mb_y;
		mb_x = (target->x - mycell->left_x)   * mycell->mfactor_x;
		mb_y = (target->y - mycell->bottom_y) * mycell->mfactor_y;
		
		if(target->x < mycell->left_x)
		{
			if(target->y < mycell->bottom_y)   mb_add_particle(mycell->sw_ghostblock, target);
			else if(target->y > mycell->top_y) mb_add_particle(mycell->nw_ghostblock, target);
			else mb_add_particle(&mycell->w_ghostblocks[mb_y], target);
		}
		else if(target->x >= mycell->right_x)
		{
			if(target->y < mycell->bottom_y)   mb_add_particle(mycell->se_ghostblock, target);
			else if(target->y > mycell->top_y) mb_add_particle(mycell->ne_ghostblock, target);
			else mb_add_particle(&mycell->e_ghostblocks[mb_y], target);
		}
		else
		{
			if(target->y < mycell->bottom_y)   mb_add_particle(&mycell->s_ghostblocks[mb_x], target);
			else if(target->y > mycell->top_y) mb_add_particle(&mycell->n_ghostblocks[mb_x], target);
		}
	}

    // Make sure that all previous ghost messages have been sent, as we need to reuse the buffers
//    MPI_Waitall(mycell->num_neighbors, mpi_ghost_requests, MPI_STATUSES_IGNORE);
	MPI_Barrier(MPI_COMM_WORLD);
}
