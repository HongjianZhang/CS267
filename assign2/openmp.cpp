#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "microblock.h"

const double micro_length  = cutoff*3;
const int default_mbuf_depth = 9;

void mb_add_particle(microblock* microblock, particle_t* particle_addr);
void mb_expand_particle(microblock* microblock, int new_max);
void mb_rm_particle(microblock* microblock, int pos);
void setup_microblocks(microblock* microblocks, int num_micro_x, int num_micro_y, double sim_x, double sim_y);
void distribute_particles(microblock* microblocks, int num_micro_x, int num_micro_y, double mfactor_x, double mfactor_y, particle_t* particles, int n);


//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    double sim_size = set_size( n );
    init_particles( n, particles );
	
	// Setup microblocks init data
	int num_micro_x, num_micro_y;
	num_micro_x = sim_size/micro_length;
	num_micro_y = sim_size/micro_length;
	
	double mfactor_x, mfactor_y; // Multiply by these to get approx microblock loc from coords
	mfactor_x = num_micro_x/sim_size;
	mfactor_y = num_micro_y/sim_size;
	
	// Initialize microblocks
	microblock* microblocks = (microblock*) malloc(num_micro_x*num_micro_y * sizeof(microblock));
	setup_microblocks(microblocks, num_micro_x, num_micro_y, sim_size, sim_size);

	// Distribute particles to microblocks
	distribute_particles(microblocks, num_micro_x, num_micro_y, mfactor_x, mfactor_y, particles, n);
	
    //
    //  simulate a number of time steps
    //
	
	double simulation_time = read_timer( );

	for( int step = 0; step < NSTEPS; step++ )
	{
		//
		//  compute all forces
		//
		// Go through each particle in each microblock
		for(int mb = 0; mb < num_micro_x*num_micro_y; ++mb)
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
					int n_id = microblocks[mb].neighbors[n];
					if(n_id != NONE) // Make sure I have a neighbor!
					{
						for(int j = 0; j < microblocks[n_id].num_particles; ++j)
						{
							apply_force(*(microblocks[mb].particles[i]), *(microblocks[n_id].particles[j]));
						}
					}
				}
			}
		}

		//
		//  move particles
		//
		// Go through each particle in each microblock
		for(int mb = 0; mb < num_micro_x*num_micro_y; ++mb)
		{
			for(int i = 0; i < microblocks[mb].num_particles; ++i)
			{
				move(*(microblocks[mb].particles[i]));
			}
		}
		
		//  migrate particles between microblocks as necessary
		for(int mb = 0; mb < num_micro_x*num_micro_y; ++mb)
		{
			for(int i = 0; i < microblocks[mb].num_particles; ++i)
			{
				particle_t* migrant = microblocks[mb].particles + i);
				if(migrant->x < microblocks[mb].left_x   ||\
				   migrant->x > microblocks[mb].right_x  ||\
				   migrant->y < microblocks[mb].bottom_y ||\
				   migrant->y < microblocks[mb].top_y)
				{
					int mb_x, mb_y;
					mb_x = migrant->x * mfactor_x;
					mb_y = migrant->y * mfactor_y;

					mb_x = min(max(0,mb_x), num_micro_x);
					mb_y = min(max(0,mb_y), num_micro_y);
				   
					mb_rm_particle(microblocks + mb, i);
					mb_add_particle(microblocks + mb_y*num_micro_x + mb_x, migrant);
				}
			}
		}

		//
		//  save if necessary
		//
		if( fsave && (step%SAVEFREQ) == 0 )
			save( fsave, n, particles );
	}
	simulation_time = read_timer( ) - simulation_time;
    
    printf("n = %d,\tsimulation time = %g seconds\n", n, simulation_time );*/
    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}

void mb_expand_particle(microblock* microblock, int new_max)
{
	// Create the new particle queue and copy old addresses into it
	particle_t** new_storage = (particle_t**) malloc(new_max * sizeof(particle_t*));
	for(int i = 0; i < microblock->num_particles; ++i)
	{
		new_storage[i] = microblock->particles[i];
	}
	
	// Free old storage
	free(microblock->particles);
	
	// Use the new particle queue
	microblock->particles = new_storage;
	microblock->max_particles = new_max;
}

void mb_add_particle(microblock* microblock, particle_t* particle_addr)
{
	// Expand queue if needed
	if(microblock->num_particles == microblock->max_particles)
	{
		mb_expand_particle(microblock, microblock->max_particles + 2);
	}
	
	// Add the new particle
	microblock->particles[microblock->num_particles] = particle_addr;
	microblock->num_particles += 1;
}
void mb_rm_particle(microblock* microblock, int pos)
{
	// Remove by overwriting target with last array value, then decrementing
	microblock->particles[pos] = microblock->particles[microblock->num_particles];
	microblock->num_particles -= 1;
}

void setup_microblocks(microblock* microblocks, int num_micro_x, int num_micro_y, double sim_x, double sim_y)
{
	for(int x = 0; x < num_micro_x; ++x)
	{
		for(int y = 0; y < num_micro_y; ++y)
		{
			microblock* current = microblocks+ y*num_micro_x + x;
			current->particles = (particle_t**) malloc(default_mbuf_depth * sizeof(particle_t*));
			current->max_particles = default_mbuf_depth;
			current->num_particles = 0;
			
			current->neighbors[p_sw] = ((x != 0)             && (y != 0)            ) ? ((y-1)*num_micro_x + (x-1)) : (NONE);
			current->neighbors[p_s ] = (                        (y != 0)            ) ? ((y-1)*num_micro_x + (x  )) : (NONE);
			current->neighbors[p_se] = ((x != num_micro_x-1) && (y != 0)            ) ? ((y-1)*num_micro_x + (x+1)) : (NONE);
			
			current->neighbors[p_w ] = ((x != 0)                                    ) ? ((y  )*num_micro_x + (x-1)) : (NONE);
			current->neighbors[p_e ] = ((x != num_micro_x-1)                        ) ? ((y  )*num_micro_x + (x+1)) : (NONE);
			
			current->neighbors[p_nw] = ((x != 0)             && (y != num_micro_y-1)) ? ((y+1)*num_micro_x + (x-1)) : (NONE);
			current->neighbors[p_n ] = (                        (y != num_micro_y-1)) ? ((y+1)*num_micro_x + (x  )) : (NONE);
			current->neighbors[p_ne] = ((x != num_micro_x-1) && (y != num_micro_y-1)) ? ((y+1)*num_micro_x + (x+1)) : (NONE);
			
			current->left_x   = (x==0)             ? (0)     : ((sim_x/num_micro_x)*(x  ));
			current->right_x  = (x==num_micro_x-1) ? (sim_x) : ((sim_x/num_micro_x)*(x+1));
			current->bottom_y = (y==0)             ? (0)     : ((sim_y/num_micro_y)*(y  ));
			current->top_y    = (y==num_micro_y-1) ? (sim_y) : ((sim_y/num_micro_y)*(y+1));
		}
	}
}

void distribute_particles(microblock* microblocks, int num_micro_x, int num_micro_y, double mfactor_x, double mfactor_y, particle_t* particles, int n)
{
	for(int i = 0; i < n; ++i)
	{
		int mb_x, mb_y;
		mb_x = particles[i].x * mfactor_x;
		mb_y = particles[i].y * mfactor_y;
		
		mb_x = min(max(0,mb_x), num_micro_x);
		mb_y = min(max(0,mb_y), num_micro_y);
		
		mb_add_particle(microblocks + mb_y*num_micro_x + mb_x, particles + i);
	}
}
