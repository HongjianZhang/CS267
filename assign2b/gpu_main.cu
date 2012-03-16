#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <cuda.h>
#include "common.h"

// === microblock layout ===
int mb_rows;
int mb_cols;

//================================================================================
//======================== Main Driver ===========================================
//================================================================================

int main( int argc, char **argv ) {
  // Synchronize with the CUDA runtime
  cudaThreadSynchronize();

  // Display help menu if desired
  if( find_option( argc, argv, "-h" ) >= 0 ) {
    printf( "Options:\n" );
    printf( "-h to see this help\n" );
    printf( "-n <int> to set the number of particles\n" );
    printf( "-o <filename> to specify the output file name\n" );
    return 0;
  }

  // Output file parameter
  char *savename = read_string(argc, argv, "-o", NULL);
  FILE *fsave = savename ? fopen(savename, "w") : NULL;

  // Number of Particles Parameter
  int n = read_int(argc, argv, "-n", 1000);

  // Initialize CPU Particle List
  particle_t *particles;
  cudaMallocHost(&particles, n * sizeof(particle_t));
  double size = set_size(n);
  init_particles(n, particles);

  // Initialize GPU Particle List
  particle_t* gpu_particles;
  cudaMalloc(&gpu_particles, n * sizeof(particle_t));
  cudaThreadSynchronize();

  // Determine number of blocks
  double mb_size = 2 * cutoff;
  mb_rows = (int)(size / mb_size);
  mb_cols = (int)(size / mb_size);

  // Initialize GPU Microblock List
  microblock* gpu_microblocks;
  cudaMalloc(&gpu_microblocks, mb_rows * mb_cols * sizeof(microblock));

  // Initialize CPU Microblock List
  microblock* cpu_microblocks = (microblock*)malloc(mb_rows * mb_cols * sizeof(microblock));

  // Synchronize mallocs
  cudaThreadSynchronize();

  // Copy the particles to the GPU
  double copy_time = read_timer();
  cudaMemcpy(gpu_particles, particles, n * sizeof(particle_t), cudaMemcpyHostToDevice);
  cudaThreadSynchronize();
  copy_time = read_timer( ) - copy_time;

  // Thread Block Structure
  int blks_p = (n + NUM_THREADS - 1) / NUM_THREADS;
  int blks_mb = (mb_rows * mb_cols + NUM_THREADS - 1) / NUM_THREADS;

  // Distribute particles into microblocks
  distribute_gpu <<< blks_mb, NUM_THREADS >>> (gpu_microblocks, mb_rows, mb_cols, gpu_particles, n, size);
  cudaThreadSynchronize();

  // Start Simulation
  double simulation_time = read_timer();
  for(int step = 0; step < NSTEPS; step++)
  {    
    // Compute Forces (one kernel per particle)
    compute_forces_gpu <<< blks_p, NUM_THREADS >>> (gpu_particles, n, gpu_microblocks, mb_rows, mb_cols);
              
    // Move particles (one kernel per particle)
    move_gpu <<< blks_p, NUM_THREADS >>> (gpu_particles, n, size);

    // Migrate particles (one kernel per microblock)
    //migrate_particles_gpu <<< blks_mb, NUM_THREADS >>> (gpu_microblocks, mb_rows, mb_cols, gpu_particles, n, size);

    // If save desired      
    if( fsave && (step%SAVEFREQ) == 0 )
    {
      cudaThreadSynchronize();

      // Copy the particles back to the CPU
      cudaMemcpy(particles, gpu_particles, n * sizeof(particle_t), cudaMemcpyDeviceToHost);
      save(fsave, n, particles);
    }
  }

  // Calculate Time Passed
  cudaThreadSynchronize();
  simulation_time = read_timer() - simulation_time;

  // Output Statistics
  printf( "CPU-GPU copy time = %g seconds\n", copy_time);
  printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );

  // Free memory
  cudaFreeHost(particles);
  cudaFree(gpu_particles);
  cudaFree(gpu_microblocks);

  // Close File
  if(fsave) fclose(fsave);
    
  return 0;
}
