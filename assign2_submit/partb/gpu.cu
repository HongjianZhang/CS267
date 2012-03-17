#include <cuda.h>
#include "common.h"

//================================================================================
//=================== Lock Free Addition And Removal =============================
//================================================================================

__device__ void add_pidx(microblock* mb, int pidx){
  for(int i=0; i<max_particles_per_mb; i++){
    //Set particle to NEWLY_ADDED if it was INVALID
    int is_active = atomicCAS(&(mb->valid[i]), INVALID, NEWLY_ADDED);

    //If it was INVALID
    if(is_active == INVALID){
      mb->p_idx[i] = pidx;
      mb->valid[i] = VALID; 
      return;
    }
  }
}

__device__ void remove_idx(microblock* mb, int i){
  mb->valid[i] = INVALID;
}

//================================================================================
//================================================================================
//================================================================================

__device__ int within_bounds(particle_gpu_t* p, double left_x, double right_x, double bottom_y, double top_y){
  return
    p->pos.x >= left_x && p->pos.x < right_x &&
    p->pos.y >= bottom_y && p->pos.y < top_y;
}

//
// Responsible for scanning through the particle list and adding the appropriate particle indices to each microblock.
// Don't forget to set microblock.n also.
//
__global__ void distribute_gpu (microblock* mb_list, int mb_rows, int mb_cols, particle_gpu_t* particles, int n, double phys_size)
{
  // Get X and Y co-ordinate of microblock
  int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

  int mb_x = thread_id % mb_cols;
  int mb_y = thread_id / mb_cols;

  // Make sure that we have a valid microblock to process
  if((mb_x >= mb_cols) || (mb_y >= mb_rows))
    return;

  // Determine boundaries of microblock
  double mb_width = phys_size / mb_cols;
  double mb_height = phys_size / mb_rows;
  double left_x = mb_x * mb_width;
  double right_x = (mb_x + 1) * mb_width;
  double bottom_y = mb_y * mb_height;
  double top_y = (mb_y + 1) * mb_height;

  //Clear valid bits
  for(int i=0; i<max_particles_per_mb; i++)
    mb_list[thread_id].valid[i] = INVALID;

  //For each particle, i
  for(int i=0; i<n; i++){
    //If the particle is within my bounds
    if(within_bounds(&particles[i], left_x, right_x, bottom_y, top_y))
    {
      //Add the particle index, i, to my list
      add_pidx(&mb_list[thread_id], i);

      //Set particle's mb index
      particles[i].mb_idx = thread_id;
    }
  }
}

__device__ void apply_force_gpu(particle_gpu_t &particle, particle_gpu_t &neighbor)
{
  double dx = neighbor.pos.x - particle.pos.x;
  double dy = neighbor.pos.y - particle.pos.y;
  double r2 = dx * dx + dy * dy;
  if( r2 > cutoff*cutoff )
      return;
  //r2 = fmax( r2, min_r*min_r );
  r2 = (r2 > min_r*min_r) ? r2 : min_r*min_r;
  double r = sqrt( r2 );

  //
  //  very simple short-range repulsive force
  //
  double coef = ( 1 - cutoff / r ) / r2 / mass;
  particle.a.x += coef * dx;
  particle.a.y += coef * dy;
}


__global__ void compute_forces_gpu (particle_gpu_t* particles, int n, microblock* mb_list, int mb_rows, int mb_cols)
{
  // Get thread (particle) ID
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= n) return;

  particle_gpu_t p = particles[tid];

  // Get X and Y co-ordinate of particle's microblock
  int mb_x = p.mb_idx % mb_cols;
  int mb_y = p.mb_idx / mb_cols;

  // Make sure that we have a valid microblock to process
  if((mb_x >= mb_cols) || (mb_y >= mb_rows))
    return;
		
  int neighbours[8];
  neighbours[p_sw] = ((mb_x != 0)             && (mb_y != 0)        ) ? ((mb_y-1)*mb_cols + (mb_x-1)) : (NO_MB);
  neighbours[p_s ] = (                           (mb_y != 0)        ) ? ((mb_y-1)*mb_cols + (mb_x  )) : (NO_MB);
  neighbours[p_se] = ((mb_x != mb_cols-1)     && (mb_y != 0)        ) ? ((mb_y-1)*mb_cols + (mb_x+1)) : (NO_MB);
  
  neighbours[p_w ] = ((mb_x != 0)                                   ) ? ((mb_y  )*mb_cols + (mb_x-1)) : (NO_MB);
  neighbours[p_e ] = ((mb_x != mb_cols-1)                           ) ? ((mb_y  )*mb_cols + (mb_x+1)) : (NO_MB);
  
  neighbours[p_nw] = ((mb_x != 0)             && (mb_y != mb_rows-1)) ? ((mb_y+1)*mb_cols + (mb_x-1)) : (NO_MB);
  neighbours[p_n ] = (                           (mb_y != mb_rows-1)) ? ((mb_y+1)*mb_cols + (mb_x  )) : (NO_MB);
  neighbours[p_ne] = ((mb_x != mb_cols-1)     && (mb_y != mb_rows-1)) ? ((mb_y+1)*mb_cols + (mb_x+1)) : (NO_MB);

  p.a.x = p.a.y = 0;

  // Collide with active particles in my block
  for(int j=0; j<max_particles_per_mb; j++)
  {
    if(mb_list[p.mb_idx].valid[j] == VALID)
    {
      int nidx = mb_list[p.mb_idx].p_idx[j];

      // If not myself
      if(nidx != tid) {
        apply_force_gpu(p, particles[nidx]);
      }
    }
  }

  // For each active neighbouring block
  for(int k=0; k<8; k++)
  {
    if(neighbours[k] != NO_MB) {
      // For each active particle in the neighbouring block
      for(int j=0; j<max_particles_per_mb; j++) {
        //Get neighbour
        microblock* n = &mb_list[neighbours[k]];
        if(n->valid[j] == VALID) {
          int nidx = n->p_idx[j];
          apply_force_gpu(p, particles[nidx]);
        }
      }
    }
  }

  // Copy particle p back to global array
  particles[tid] = p;
}


__global__ void move_gpu (particle_gpu_t * particles, int n, double size)
{
  // Get thread (particle) ID
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= n) return;

  particle_gpu_t p = particles[tid];
  
  //
  //  slightly simplified Velocity Verlet integration
  //  conserves energy better than explicit Euler method
  //
  p.v.x += p.a.x * dt;
  p.v.y += p.a.y * dt;
  p.pos.x  += p.v.x * dt;
  p.pos.y  += p.v.y * dt;

  //
  //  bounce from walls
  //
  while( p.pos.x < 0 || p.pos.x > size )
  {
      p.pos.x  = p.pos.x < 0 ? -(p.pos.x) : 2*size-p.pos.x;
      p.v.x = -(p.v.x);
  }
  while( p.pos.y < 0 || p.pos.y > size )
  {
      p.pos.y  = p.pos.y < 0 ? -(p.pos.y) : 2*size-p.pos.y;
      p.v.y = -(p.v.y);
  }

  // copy result back to global memory
  particles[tid] = p;
}


__global__ void migrate_particles_gpu (microblock* mb_list, int mb_rows, int mb_cols, particle_gpu_t* particles, int n, double phys_size)
{
  // Get X and Y co-ordinate of microblock
  int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

  int mb_x = thread_id % mb_cols;
  int mb_y = thread_id / mb_cols;

  // Make sure that we have a valid microblock to process
  if((mb_x >= mb_cols) || (mb_y >= mb_rows))
    return;

  // Determine boundaries of microblock
  double mb_width = phys_size / mb_cols;
  double mb_height = phys_size / mb_rows;
  double left_x = mb_x * mb_width;
  double right_x = (mb_x + 1) * mb_width;
  double bottom_y = mb_y * mb_height;
  double top_y = (mb_y + 1) * mb_height;

  // For each index i in my microblock	
  for(int i = 0; i < max_particles_per_mb; i++){
    microblock* mb = &mb_list[thread_id];
    //If this index is valid
    if(mb->valid[i] == VALID){
      //Get current particle
      int pidx = mb->p_idx[i];
      particle_gpu_t* p = &particles[pidx];

      //If particle is no longer within my bounds
      if(!within_bounds(p, left_x, right_x, bottom_y, top_y)){
        
        //Remove from current microblock
        remove_idx(mb, i);
        
        //and add to the appropriate neighbour
        //neighbour to add it to.
        int neighbour_x, neighbour_y;
        neighbour_x = p->pos.x / mb_width;
        neighbour_y = p->pos.y / mb_height;

        //Make sure we don't leave the simulation space
        neighbour_x = min(max(0,neighbour_x), mb_cols-1);
        neighbour_y = min(max(0,neighbour_y), mb_rows-1);

        //Neighbour to add to:
        microblock* neighbour = &mb_list[neighbour_y * mb_cols + neighbour_x];
        add_pidx(neighbour, pidx);
        
        //Update particle's own record of its mb
        p->mb_idx = neighbour_y * mb_cols + neighbour_x;
      }
    }
  }
}
