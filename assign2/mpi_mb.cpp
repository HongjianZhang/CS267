#include "mpi_mb.h"

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
	
	mycell->nw_ghostblock = (microblock*) malloc(sizeof(microblock));
	mycell->ne_ghostblock = (microblock*) malloc(sizeof(microblock));
	mycell->sw_ghostblock = (microblock*) malloc(sizeof(microblock));
	mycell->se_ghostblock = (microblock*) malloc(sizeof(microblock));
	
	mycell->n_ghostblocks = (microblock*) malloc(mycell->num_micro_x * sizeof(microblock));
	mycell->s_ghostblocks = (microblock*) malloc(mycell->num_micro_x * sizeof(microblock));
	mycell->e_ghostblocks = (microblock*) malloc(mycell->num_micro_y * sizeof(microblock));
	mycell->w_ghostblocks = (microblock*) malloc(mycell->num_micro_y * sizeof(microblock));
	
	return target;
}

void setup_microblocks(microblock* microblocks, mpi_cell* mycell);
{
	for(int x = 0; x < mycell->num_micro_x; ++x)
	{
		for(int y = 0; y < mycell->num_micro_y; ++y)
		{
			microblock* current = microblocks+y*num_micro_x + x;
			current->particles = (particle_t**) malloc(default_mbuf_depth * sizeof(particle_t*));
			current->max_particles = default_mbuf_depth;
			current->num_particles = 0;
			
			current->neighbors[p_sw] = ((x != 0)                     && (y != 0)                    ) ? (&microblocks[(y-1)*mycell->num_micro_x + (x-1)]) : (&(mycell->sw_ghostblock   ));
			current->neighbors[p_s ] = (                                (y != 0)                    ) ? (&microblocks[(y-1)*mycell->num_micro_x + (x  )]) : (&(mycell->s_ghostblocks[x]));
			current->neighbors[p_se] = ((x != mycell->num_micro_x-1) && (y != 0)                    ) ? (&microblocks[(y-1)*mycell->num_micro_x + (x+1)]) : (&(mycell->se_ghostblock   ));
			
			current->neighbors[p_w ] = ((x != 0)                                                    ) ? (&microblocks[(y  )*mycell->num_micro_x + (x-1)]) : (&(mycell->w_ghostblocks[y]));
			current->neighbors[p_e ] = ((x != mycell->num_micro_x-1)                                ) ? (&microblocks[(y  )*mycell->num_micro_x + (x+1)]) : (&(mycell->e_ghostblocks[y]));
			
			current->neighbors[p_nw] = ((x != 0)                     && (y != mycell->num_micro_y-1)) ? (&microblocks[(y+1)*mycell->num_micro_x + (x-1)]) : (&(mycell->nw_ghostblock   ));
			current->neighbors[p_n ] = (                                (y != mycell->num_micro_y-1)) ? (&microblocks[(y+1)*mycell->num_micro_x + (x  )]) : (&(mycell->n_ghostblocks[x]));
			current->neighbors[p_ne] = ((x != mycell->num_micro_x-1) && (y != mycell->num_micro_y-1)) ? (&microblocks[(y+1)*mycell->num_micro_x + (x+1)]) : (&(mycell->ne_ghostblock   ));
			
			current->left_x   = (x==0)                     ? (left_x)   : (((mycell->right_x-mycell->left_x)/mycell->num_micro_x)*(x  ));
			current->right_x  = (x==mycell->num_micro_x-1) ? (right_x)  : (((mycell->right_x-mycell->left_x)/mycell->num_micro_x)*(x+1));
			current->bottom_y = (y==0)                     ? (bottom_y) : (((mycell->top_y-mycell->bottom_y)/mycell->num_micro_y)*(y  ));
			current->top_y    = (y==mycell->num_micro_y-1) ? (top_y)    : (((mycell->top_y-mycell->bottom_y)/mycell->num_micro_y)*(y+1));
		}
	}
}
void distribute_particles(microblock* microblocks, mpi_cell* mycell, plist* local, particle_t* particles, int n)
{
	for(int i = 0; i < n; ++i)
	{
		int mb_x, mb_y;
		mb_x = particles[i].x * mycell->mfactor_x;
		mb_y = particles[i].y * mycell->mfactor_y;
		
		if(mb_x < 0 || mb_x > mycell->num_micro_x || mb_y < 0 || mb_y > mycell->num_micro_y) continue;
		
		particle_t* added_part = add_particle(local, particles[i]);
		mb_add_particle(microblocks + mb_y*num_micro_x + mb_x, added_part);
	}
}