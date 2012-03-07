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
	
	mpi_cell* mycell = create_mpi_cell(sim_size, num_proc_x, num_proc_y, rank);
	
	microblock* microblocks = (microblock*) malloc(num_micro_x*num_micro_y * sizeof(microblock));
	setup_microblocks(microblocks, mycell);

	//
	//  allocate storage for local particles, ghost particles, ids
	//
	plist* local = alloc_plist(n);
	plist* ghost = alloc_ppile(n);
	
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
			prepare_save(rank, n_proc, &local, particles, n);
			
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
	free( mycell );
	
	free_plist(local);
	free_ppile(ghost);
	
	free_emigrant_buf();
	clean_ghost_structure();
	
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}

void prepare_save(int rank, int n_proc, plist* local, particle_t* particles, int n);
{
	// First, get the number of particles in each node into node 0. Also prepare array placement offsets.
	int* node_particles_num    = (int *) malloc(n_proc*sizeof(int));
	int* node_particles_offset = (int *) malloc(n_proc*sizeof(int));
	
	MPI_Gather(&local->num_used_ids, 1, MPI_INT, node_particles_num, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if(rank == 0)
	{
		node_particles_offset[0] = 0;
		for(int i = 1; i < n_proc; ++i)
		{
			node_particles_offset[i] = node_particles_offset[i-1] + node_particles_num[i-1];
		}
	}
	
	// Now, each node prepares a collapsed list of all valid particles
	particle_t* collapsed_local = (particle_t *) malloc(local->num_used_ids * sizeof(particle_t));
	
	int cp = 0;
	for(int i = 0; i < local->end_particle; ++i)
	{
		if(local->is_id_active[i] == 0) continue
		
		particle_t packing = *(local->particles[i]);
		collapsed_local[cp] = packing;
		
		cp++;
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