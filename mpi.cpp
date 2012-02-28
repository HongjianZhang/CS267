#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "mpi_particles.h"

MPI_Datatype PARTICLE;

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
    
    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    
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
	left_x   = (proc_x==0)            ? (0)    : ((size/num_proc_x)*proc_x);
	right_x  = (proc_x==num_proc_x-1) ? (size) : ((size/num_proc_x)*(proc_x+1));
	bottom_y = (proc_y==0)            ? (0)    : ((size/num_proc_y)*proc_y);
	top_y    = (proc_y==num_proc_y-1) ? (size) : ((size/num_proc_y)*(proc_y+1));
	
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
	
    //
    //  allocate storage for local particles, ghost particles
    //
    particle_t *local = (particle_t*) malloc( n * sizeof(particle_t) );
	char *p_valid = (char*) malloc(n*sizeof(char));
	int nlocal;
	
	setup_ghost_structure(n);
	init_emigrant_buf(n);

    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //

	particle_t* particles = (particle_t*) malloc(n * sizeof(particle_t));
	
    set_size( n );
    if( rank == 0 )
        init_particles( n, particles);
	
	MPI_Bcast((void *) particles, n, PARTICLES, 0, MPI_COMM_WORLD);
	select_particles(n, particles, local, p_valid, &nlocal); //TODO WRITE ME
	
	ghost_particles = (particle_t *) malloc(max_particles * sizeof(particle_t));
	int nghosts = 0;
	
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
		//
		//  Handle ghosting
		//
		prepare_ghost_packets(local, p_valid, nlocal, left_x, right_x, bottom_y, top_y, neighbors);
		send_ghost_packets(neighbors);
		receive_ghost_packets(&nghosts, ghost_particles, neighbors, num_neighbors, n);
		
        //
        //  Compute all forces
        //
		compute_forces(local, p_valid, nlocal, ghost_particles, nghosts);
        
        //
        //  Move particles
        //
		int seen_particles = 0;
		for(int i = 0; seen_particles < num_particles; ++i)
		{
			if(p_valid[i] == INVALID) continue;
			seen_particles++;
			
            move( local[i] );
		}
		
		//
		//  Handle migration
		//

		prepare_emigrants(local, p_valid, &num_particles, left_x, right_x, bottom_y, top_y, neighbors);
		send_emigrants(neighbors);
		receive_immigrants(neighbors, num_neighbors, local, p_valid, n, n);
		
		//
        //  save current step if necessary
        //
		if(savename && (step%SAVEFREQ) == 0)
		{
			prepare_save(rank, local, p_valid, nlocal, particles, n);
			
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
    free( local );
    free( particles );
	free( ghost_particles );
	
	free_emigrant_buf();
	clean_ghost_structure();
	
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}

void compute_forces(particle_t local[], char p_valid[], int num_particles, particle_t ghosts[], int num_ghosts)
{
	int seen_particles = 0;
	for(int i = 0; seen_particles < num_particles; ++i)
	{
		if(p_valid[i] == INVALID) continue;
		seen_particles++;
		
		local[i].ax = local[i].ay = 0;
		int nearby_seen_particles = 0;
		for (int j = 0; nearby_seen_particles < num_particles; ++j)
		{
			if(p_valid[j] == INVALID) continue;
			nearby_seen_particles++;
			
			apply_force( local[i], local[j] );
		}
		
		for(int j; j < num_ghosts; ++j)
		{
			apply_force( local[i], ghosts[j]);
		}
	}
}