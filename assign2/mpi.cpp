#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <algorithm>
#include "common.h"
#include "mpi_particles.h"
#include "microblock.h"

MPI_Datatype PARTICLE;

void mb_add_particle(microblock* microblock, particle_t* particle_addr);
void mb_expand_particle(microblock* microblock, int new_max);
void mb_rm_particle(microblock* microblock, int pos);
void setup_microblocks(microblock* microblocks, int num_micro_x, int num_micro_y, double left_x, double right_x, double bottom_y, double top_y);
void distribute_particles(microblock* microblocks, int num_micro_x, int num_micro_y, double left_x, double bottom_y, double mfactor_x, double mfactor_y, plist* local, particle_t* particles, int n);

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    //
    //  process command line parameters
    //
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );
    
    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	
    double sim_size = set_size(n);
    
    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    
    MPI_Type_contiguous( 7, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );
    
    //
    //  set up the data partitioning across processors
    //
	// Determine how many cells are along each axis
	int num_proc_x, num_proc_y;
	for(int factor = (int)floor(sqrt((double)n_proc)); factor >= 1; --factor)
	{
		if(n_proc % factor == 0)
		{
			num_proc_x = factor;
			num_proc_y = n_proc/factor;
			break;
		}
	}
	
	// Determine where this cell is
	int proc_x, proc_y;
	proc_x = rank % num_proc_x;
	proc_y = rank/num_proc_x;
	
	// Determine my cell boundaries
	double left_x, right_x, bottom_y, top_y;
	left_x   = (proc_x==0)            ? (0)        : ((sim_size/num_proc_x)*proc_x);
	right_x  = (proc_x==num_proc_x-1) ? (sim_size) : ((sim_size/num_proc_x)*(proc_x+1));
	bottom_y = (proc_y==0)            ? (0)        : ((sim_size/num_proc_y)*proc_y);
	top_y    = (proc_y==num_proc_y-1) ? (sim_size) : ((sim_size/num_proc_y)*(proc_y+1));

	// Determine the ranks of my neighbors for message passing, NONE means no neighbor
	int neighbors[8];
	int num_neighbors = 0;
	// In order, positions are SW, S, SE, W, E, NW, N, NE
	neighbors[p_sw] = ((proc_x != 0)            && (proc_y != 0)           ) ? ((proc_y-1)*num_proc_x + (proc_x-1)) : (NONE);
	neighbors[p_s ] = (                            (proc_y != 0)           ) ? ((proc_y-1)*num_proc_x + (proc_x  )) : (NONE);
	neighbors[p_se] = ((proc_x != num_proc_x-1) && (proc_y != 0)           ) ? ((proc_y-1)*num_proc_x + (proc_x+1)) : (NONE);
	
	neighbors[p_w ] = ((proc_x != 0)                                       ) ? ((proc_y  )*num_proc_x + (proc_x-1)) : (NONE);
	neighbors[p_e ] = ((proc_x != num_proc_x-1)                            ) ? ((proc_y  )*num_proc_x + (proc_x+1)) : (NONE);
	
	neighbors[p_nw] = ((proc_x != 0)            && (proc_y != num_proc_y-1)) ? ((proc_y+1)*num_proc_x + (proc_x-1)) : (NONE);
	neighbors[p_n ] = (                            (proc_y != num_proc_y-1)) ? ((proc_y+1)*num_proc_x + (proc_x  )) : (NONE);
	neighbors[p_ne] = ((proc_x != num_proc_x-1) && (proc_y != num_proc_y-1)) ? ((proc_y+1)*num_proc_x + (proc_x+1)) : (NONE);
    
	for(int i = 0 ; i < 8; ++i)
	{
		if(neighbors[i] != NONE) num_neighbors++;
	}
	
	// Setup microblocks init data
	cell_xw = right_x - left_x;
	cell_yw = top_y   - bottom_y;
	
	int num_micro_x, num_micro_y;
	num_micro_x = (right_x - left_x)/micro_length;
	num_micro_y = (top_y   - bottom_y)/micro_length;
	
	double mfactor_x, mfactor_y; // Multiply by these to get approx microblock loc from coords
	mfactor_x = num_micro_x/(right_x - left_x);
	mfactor_y = num_micro_y/(top_y   - bottom_y);
	
	microblock* microblocks = (microblock*) malloc(num_micro_x*num_micro_y * sizeof(microblock));
	setup_microblocks(microblocks, num_micro_x, num_micro_y, left_x, right_x, bottom_y, top_y);
    
    //
    //  allocate storage for local particles, ghost particles, ids
    //
	plist* local = alloc_plist(n);
	
    particle_t *ghost = (particle_t*) malloc( n * sizeof(particle_t) );
	int nghost = 0;
	
	setup_ghost_structure(n);
	init_emigrant_buf(n);

    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //

	particle_t* particles = (particle_t*) malloc(n * sizeof(particle_t));
	
	if( rank == 0 )
		init_particles( n, particles);
	
	MPI_Bcast((void *) particles, n, PARTICLE, 0, MPI_COMM_WORLD);
	distribute_particles(microblocks, num_micro_x, num_micro_y, left_x, bottom_y, mfactor_x, mfactor_y, local, particles, n);

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer();
    for( int step = 0; step < NSTEPS; step++ )
    {
		//
		//  Handle ghosting
		//
		prepare_ghost_packets(part, local_ids, nlocal, left_x, right_x, bottom_y, top_y, neighbors);
		send_ghost_packets(neighbors);
		receive_ghost_packets(part, ghost_ids, &nghost, neighbors, num_neighbors, n);
	
        //  Compute all forces
        //
		update_particles(part);
		
		//
		//  Handle migration
		//
		prepare_emigrants(part, local_ids, &nlocal, left_x, right_x, bottom_y, top_y, neighbors);
		send_emigrants(neighbors);
		receive_immigrants(neighbors, num_neighbors, part, local_ids, &nlocal, n);
		
		//
        //  save current step if necessary
        //
		if(savename && (step%SAVEFREQ) == 0)
		{
			prepare_save(part, rank, n_proc, local_ids, nlocal, particles, n);
			
			if(fsave)
				save( fsave, n, particles );
		}
    }
    simulation_time = read_timer( ) - simulation_time;
    
    if( rank == 0 )
        printf( "n = %d, n_procs = %d, simulation time = %g s\n", n, n_proc, simulation_time );
    
    //
    //  release resources
    //
    free( particles );
	
	free_emigrant_buf();
	clean_ghost_structure();
	
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}

void prepare_save(partition_t* part, int rank, int n_proc, int* local_ids, int nlocal, particle_t* particles, int n)
{
	// First, get the number of particles in each node into node 0. Also prepare array placement offsets.
	int* node_particles_num    = (int *) malloc(n_proc*sizeof(int));
	int* node_particles_offset = (int *) malloc(n_proc*sizeof(int));
	
	MPI_Gather(&nlocal, 1, MPI_INT, node_particles_num, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if(rank == 0)
	{
		node_particles_offset[0] = 0;
		for(int i = 1; i < n_proc; ++i)
		{
			node_particles_offset[i] = node_particles_offset[i-1] + node_particles_num[i-1];
		}
	}
	
	// Now, each node prepares a collapsed list of all valid particles
	particle_t* collapsed_local = (particle_t *) malloc(nlocal * sizeof(particle_t));
	
	for(int i = 0; i < nlocal; ++i)
	{
		particle_t packing = *(get_particle(part, local_ids[i]));
		collapsed_local[i] = packing;
	}
	
	// Next, send the particles to node 0
	MPI_Gatherv(collapsed_local, nlocal, PARTICLE, particles, node_particles_num, node_particles_offset, PARTICLE, 0, MPI_COMM_WORLD);
	
	// Finally, sort the particles at node 0
	if(rank == 0)
	{
		sort(particles, particles + n, compare_particles);
	}
	
	
	// Clean up
	free(collapsed_local);
	free(node_particles_num);
	free(node_particles_offset);
}
bool compare_particles(particle_t left, particle_t right) // check if id < id
{
	return left.globalID < right.globalID;
}

int select_particles(partition_t* part, int n, particle_t* particles, int* local_ids, double left_x, double right_x, double bottom_y, double top_y)
// Note, does not take particles precisely on the top or right border.
{
	int current_particle = 0;
	
	for(int i = 0; i < n; ++i)
	{
		if((particles[i].x >= left_x) && (particles[i].x < right_x) && (particles[i].y >= bottom_y) && (particles[i].y < top_y))
		// Particle in my box, take it
		{
			local_ids[current_particle] = add_particle(part, particles[i]);
			
			current_particle++;
		}
	}
	
	// Make sure we know how many local particles we have.
	return current_particle;
}
void setup_microblocks(microblock* microblocks, int num_micro_x, int num_micro_y, double left_x, double right_x, double bottom_y, double top_y);
{
	for(int x = 0; x < num_micro_x; ++x)
	{
		for(int y = 0; y < num_micro_y; ++y)
		{
			microblock* current = microblocks+ y*num_micro_x + x;
			current->particles = (particle_t**) malloc(default_mbuf_depth * sizeof(particle_t*));
			current->max_particles = default_mbuf_depth;
			current->num_particles = 0;
			
			current->neighbors[p_sw] = ((x != 0)             && (y != 0)            ) ? (&microblocks[(y-1)*num_micro_x + (x-1)]) : (NO_MB);
			current->neighbors[p_s ] = (                        (y != 0)            ) ? (&microblocks[(y-1)*num_micro_x + (x  )]) : (NO_MB);
			current->neighbors[p_se] = ((x != num_micro_x-1) && (y != 0)            ) ? (&microblocks[(y-1)*num_micro_x + (x+1)]) : (NO_MB);
			
			current->neighbors[p_w ] = ((x != 0)                                    ) ? (&microblocks[(y  )*num_micro_x + (x-1)]) : (NO_MB);
			current->neighbors[p_e ] = ((x != num_micro_x-1)                        ) ? (&microblocks[(y  )*num_micro_x + (x+1)]) : (NO_MB);
			
			current->neighbors[p_nw] = ((x != 0)             && (y != num_micro_y-1)) ? (&microblocks[(y+1)*num_micro_x + (x-1)]) : (NO_MB);
			current->neighbors[p_n ] = (                        (y != num_micro_y-1)) ? (&microblocks[(y+1)*num_micro_x + (x  )]) : (NO_MB);
			current->neighbors[p_ne] = ((x != num_micro_x-1) && (y != num_micro_y-1)) ? (&microblocks[(y+1)*num_micro_x + (x+1)]) : (NO_MB);
			
			current->left_x   = (x==0)             ? (left_x)   : (((right_x-left_x)/num_micro_x)*(x  ));
			current->right_x  = (x==num_micro_x-1) ? (right_x)  : (((right_x-left_x)/num_micro_x)*(x+1));
			current->bottom_y = (y==0)             ? (bottom_y) : (((top_y-bottom_y)/num_micro_y)*(y  ));
			current->top_y    = (y==num_micro_y-1) ? (top_y)    : (((top_y-bottom_y)/num_micro_y)*(y+1));

			omp_init_lock(&(current->lock));
		}
	}
}

void distribute_particles(microblock* microblocks, int num_micro_x, int num_micro_y, double left_x, double bottom_y, double mfactor_x, double mfactor_y, plist* local, particle_t* particles, int n);
{
	for(int i = 0; i < n; ++i)
	{
		int mb_x, mb_y;
		mb_x = particles[i].x * mfactor_x;
		mb_y = particles[i].y * mfactor_y;
		
		if(mb_x < 0 || mb_x > num_micro_x || mb_y < 0 || mb_y > num_micro_y) continue;
		
		particle_t* added_part = add_particle(local, particles[i]);
		
		mb_add_particle(microblocks + mb_y*num_micro_x + mb_x, added_part);
	}
}
