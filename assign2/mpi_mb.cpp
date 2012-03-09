#include "mpi_mb.h"
#include <stdlib.h>

mpi_cell* create_mpi_cell(double sim_size, int num_proc_x, int num_proc_y, int rank)
{
	mpi_cell* target = (mpi_cell*) malloc(sizeof(mpi_cell));

	// Determine where this cell is
	target->rank = rank;
	target->proc_x = rank % num_proc_x;
	target->proc_y = rank/num_proc_x;
	
	// Determine my cell boundaries
	target->left_x   = (target->proc_x==0)            ? (0)        : ((sim_size/num_proc_x)* target->proc_x);
	target->right_x  = (target->proc_x==num_proc_x-1) ? (sim_size) : ((sim_size/num_proc_x)*(target->proc_x+1));
	target->bottom_y = (target->proc_y==0)            ? (0)        : ((sim_size/num_proc_y)* target->proc_y);
	target->top_y    = (target->proc_y==num_proc_y-1) ? (sim_size) : ((sim_size/num_proc_y)*(target->proc_y+1));

	// Determine the ranks of my neighbors for message passing, NONE means no neighbor
	target->neighbors[p_sw] = ((target->proc_x != 0)            && (target->proc_y != 0)           ) ? ((target->proc_y-1)*num_proc_x + (target->proc_x-1)) : (NONE);
	target->neighbors[p_s ] = (                                    (target->proc_y != 0)           ) ? ((target->proc_y-1)*num_proc_x + (target->proc_x  )) : (NONE);
	target->neighbors[p_se] = ((target->proc_x != num_proc_x-1) && (target->proc_y != 0)           ) ? ((target->proc_y-1)*num_proc_x + (target->proc_x+1)) : (NONE);
	
	target->neighbors[p_w ] = ((target->proc_x != 0)                                               ) ? ((target->proc_y  )*num_proc_x + (target->proc_x-1)) : (NONE);
	target->neighbors[p_e ] = ((target->proc_x != num_proc_x-1)                                    ) ? ((target->proc_y  )*num_proc_x + (target->proc_x+1)) : (NONE);
	
	target->neighbors[p_nw] = ((target->proc_x != 0)            && (target->proc_y != num_proc_y-1)) ? ((target->proc_y+1)*num_proc_x + (target->proc_x-1)) : (NONE);
	target->neighbors[p_n ] = (                                    (target->proc_y != num_proc_y-1)) ? ((target->proc_y+1)*num_proc_x + (target->proc_x  )) : (NONE);
	target->neighbors[p_ne] = ((target->proc_x != num_proc_x-1) && (target->proc_y != num_proc_y-1)) ? ((target->proc_y+1)*num_proc_x + (target->proc_x+1)) : (NONE);
    
	for(int i = 0 ; i < 8; ++i)
	{
		if(target->neighbors[i] != NONE) target->num_neighbors += 1;
	}
	
	// Setup microblocks init data
	target->num_micro_x = (target->right_x - target->left_x)/micro_length;
	target->num_micro_y = (target->top_y   - target->bottom_y)/micro_length;
	
	target->mfactor_x = target->num_micro_x/(target->right_x - target->left_x);
	target->mfactor_y = target->num_micro_y/(target->top_y   - target->bottom_y);
	
	// Allocate and initialize the ghost microblocks
	target->nw_ghostblock = (microblock*) malloc(sizeof(microblock));
	target->ne_ghostblock = (microblock*) malloc(sizeof(microblock));
	target->sw_ghostblock = (microblock*) malloc(sizeof(microblock));
	target->se_ghostblock = (microblock*) malloc(sizeof(microblock));
	
	target->n_ghostblocks = (microblock*) malloc(target->num_micro_x * sizeof(microblock));
	target->s_ghostblocks = (microblock*) malloc(target->num_micro_x * sizeof(microblock));
	target->e_ghostblocks = (microblock*) malloc(target->num_micro_y * sizeof(microblock));
	target->w_ghostblocks = (microblock*) malloc(target->num_micro_y * sizeof(microblock));
	
	mb_init(target->nw_ghostblock); mb_init(target->ne_ghostblock);
	mb_init(target->sw_ghostblock); mb_init(target->se_ghostblock);
	
	for(int x = 0; x < target->num_micro_x; ++x)
	{
		mb_init(&(target->n_ghostblocks[x]));
		mb_init(&(target->s_ghostblocks[x]));
	}
	for(int y = 0; y < target->num_micro_y; ++y)
	{
		mb_init(&(target->e_ghostblocks[y]));
		mb_init(&(target->w_ghostblocks[y]));
	}
	
	return target;
}

void setup_microblocks(microblock* microblocks, mpi_cell* mycell)
{
	for(int x = 0; x < mycell->num_micro_x; ++x)
	{
		for(int y = 0; y < mycell->num_micro_y; ++y)
		{
			microblock* current = microblocks+y*mycell->num_micro_x + x;
			mb_init(current);
			
			// Setup cell boundaries
			current->left_x   = (x==0)                     ? (mycell->left_x)   : (mycell->left_x   + ((mycell->right_x-mycell->left_x)/mycell->num_micro_x)*(x  ));
			current->right_x  = (x==mycell->num_micro_x-1) ? (mycell->right_x)  : (mycell->left_x   + ((mycell->right_x-mycell->left_x)/mycell->num_micro_x)*(x+1));
			current->bottom_y = (y==0)                     ? (mycell->bottom_y) : (mycell->bottom_y + ((mycell->top_y-mycell->bottom_y)/mycell->num_micro_y)*(y  ));
			current->top_y    = (y==mycell->num_micro_y-1) ? (mycell->top_y)    : (mycell->bottom_y + ((mycell->top_y-mycell->bottom_y)/mycell->num_micro_y)*(y+1));
			
			// Link microblock to neighboring microblocks
			if(y == 0)
			{
				current->neighbors[p_sw] = (x != 0)                     ? (&(mycell->s_ghostblocks[x-1])) : ((mycell->sw_ghostblock));
				current->neighbors[p_s ] =                                (&(mycell->s_ghostblocks[x  ]));
				current->neighbors[p_se] = (x != mycell->num_micro_x-1) ? (&(mycell->s_ghostblocks[x+1])) : ((mycell->se_ghostblock));
			}
			else
			{
				current->neighbors[p_sw] = (x != 0)                     ? (&microblocks[(y-1)*mycell->num_micro_x + (x-1)]) : (&(mycell->w_ghostblocks[y-1]));
				current->neighbors[p_s ] =                                (&microblocks[(y-1)*mycell->num_micro_x + (x  )]);
				current->neighbors[p_se] = (x != mycell->num_micro_x-1) ? (&microblocks[(y-1)*mycell->num_micro_x + (x+1)]) : (&(mycell->e_ghostblocks[y-1]));
			}
			
			current->neighbors[p_w] = ((x != 0)                      ) ? (&microblocks[(y  )*mycell->num_micro_x + (x-1)]) : (&(mycell->w_ghostblocks[y]));
			current->neighbors[p_e] = ((x != mycell->num_micro_x-1)  ) ? (&microblocks[(y  )*mycell->num_micro_x + (x+1)]) : (&(mycell->e_ghostblocks[y]));
			
			if(y == mycell->num_micro_y-1)
			{
				current->neighbors[p_nw] = (x != 0)                     ? (&(mycell->n_ghostblocks[x-1])) : ((mycell->nw_ghostblock));
				current->neighbors[p_n ] =                                (&(mycell->n_ghostblocks[x  ]));
				current->neighbors[p_ne] = (x != mycell->num_micro_x-1) ? (&(mycell->n_ghostblocks[x+1])) : ((mycell->ne_ghostblock));
			}
			else
			{
				current->neighbors[p_nw] = (x != 0)                     ? (&microblocks[(y+1)*mycell->num_micro_x + (x-1)]) : (&(mycell->w_ghostblocks[y+1]));
				current->neighbors[p_n ] =                                (&microblocks[(y+1)*mycell->num_micro_x + (x  )]);
				current->neighbors[p_ne] = (x != mycell->num_micro_x-1) ? (&microblocks[(y+1)*mycell->num_micro_x + (x+1)]) : (&(mycell->e_ghostblocks[y+1]));
			}
		}
	}
}
void distribute_particles(microblock* microblocks, mpi_cell* mycell, plist* local, particle_t* particles, int n)
{
	for(int i = 0; i < n; ++i)
	{
		particle_t target = particles[i];
		if(target.x >= mycell->left_x && target.x < mycell->right_x && target.y > mycell->bottom_y && target.y < mycell->top_y)
		{
			int mb_x, mb_y;
			mb_x = (target.x - mycell->left_x)   * mycell->mfactor_x;
			mb_y = (target.y - mycell->bottom_y) * mycell->mfactor_y;
			particle_t* added_part = add_particle(local, target);
			mb_add_particle(microblocks + mb_y*(mycell->num_micro_x) + mb_x, added_part);
		}
	}
}

void process_particles(microblock* microblocks, int num_mb)
{
	for(int mb = 0; mb < num_mb; ++mb)
	{
		for(int i = 0; i < microblocks[mb].num_particles; ++i)
		{
			// Reset forces for current particle
			microblocks[mb].particles[i]->ax = 0; microblocks[mb].particles[i]->ay = 0;
			
			// Collide with particles in my block
			for(int j = 0; j < microblocks[mb].num_particles; ++j)
			{
				if(i != j) // do not collide with self
				{
					apply_force(*(microblocks[mb].particles[i]), *(microblocks[mb].particles[j]));
				}
			}
			
			// Collide with particles in neighboring blocks
			for(int n = 0; n < 8; ++n)
			{
				if(microblocks[mb].neighbors[n] != NO_MB) // Make sure I have a neighbor!
				{
					for(int j = 0; j < microblocks[mb].neighbors[n]->num_particles; ++j)
					{
						apply_force(*(microblocks[mb].particles[i]), *(microblocks[mb].neighbors[n]->particles[j]));
					}
				}
			}
		}
	}

	//
	//  move particles
	//
	// Go through each particle in each microblock
	for(int mb = 0; mb < num_mb; ++mb)
	{
		for(int i = 0; i < microblocks[mb].num_particles; ++i)
		{
			move(*(microblocks[mb].particles[i]));
		}
	}
}
