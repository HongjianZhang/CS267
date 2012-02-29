#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"

//================================================================================
//==================== Main Driver ===============================================
//================================================================================
void run_simulation(particle_t* ps, int n, FILE* fsave){
  //Create partition and set active
  partition_t* p1 = alloc_partition(n);
  for(int i=0; i<n; i++)
    int id = add_particle(p1, ps[i]);

  //For each step
  for(int step = 0; step < NSTEPS; step++ ){
    update_particles(p1);

    if(fsave && (step%SAVEFREQ) == 0){      
      for(int i=0; i<p1->num_particles; i++)
        if(p1->is_id_active[i])
          ps[p1->particles[i].globalID] = p1->particles[i];
    
      save(fsave, n, ps );
    }
  }
}

int main( int argc, char **argv )
{    
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
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    run_simulation(particles, n, fsave);
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
