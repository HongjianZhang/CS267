#include <cuda.h>
#include "common.h"

//
// Responsible for scanning through the particle list and adding the appropriate particle indices to each microblock.
// Don't forget to set microblock.n also.
//
__global__ void distribute_gpu (microblock* mb_list, int mb_rows, int mb_cols, particle_t* particles, int n, double phys_size)
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

    // Iterate through array of particles
    for(int i = 0; i < n; i++)
    {
        // If particle lies within this microblock, add it
        if( (particles[i].x >= left_x) && (particles[i].x < right_x) &&
            (particles[i].y >= bottom_y) && (particles[i].y < top_y) )
        {
            mb_list[thread_id].p_idx[mb_list[thread_id].n] = i;
            mb_list[thread_id].n += 1;
        }
    }
}

__device__ void apply_force_gpu(particle_t &particle, particle_t &neighbor)
{
  double dx = neighbor.x - particle.x;
  double dy = neighbor.y - particle.y;
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
  particle.ax += coef * dx;
  particle.ay += coef * dy;

}

__global__ void compute_forces_gpu (microblock* mb_list, int mb_rows, int mb_cols, particle_t* particles)
{
    // Get X and Y co-ordinate of microblock
	int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

    int mb_x = thread_id % mb_cols;
    int mb_y = thread_id / mb_cols;

    // Make sure that we have a valid microblock to process
    if((mb_x >= mb_cols) || (mb_y >= mb_rows))
        return;
		
	int neighbors[8];
	neighbors[p_sw] = ((mb_x != 0)             && (mb_y != 0)        ) ? ((mb_y-1)*mb_cols + (mb_x-1)) : (NO_MB);
	neighbors[p_s ] = (                           (mb_y != 0)        ) ? ((mb_y-1)*mb_cols + (mb_x  )) : (NO_MB);
	neighbors[p_se] = ((mb_x != mb_cols-1)     && (mb_y != 0)        ) ? ((mb_y-1)*mb_cols + (mb_x+1)) : (NO_MB);

	neighbors[p_w ] = ((mb_x != 0)                                   ) ? ((mb_y  )*mb_cols + (mb_x-1)) : (NO_MB);
	neighbors[p_e ] = ((mb_x != mb_cols-1)                           ) ? ((mb_y  )*mb_cols + (mb_x+1)) : (NO_MB);

	neighbors[p_nw] = ((mb_x != 0)             && (mb_y != mb_rows-1)) ? ((mb_y+1)*mb_cols + (mb_x-1)) : (NO_MB);
	neighbors[p_n ] = (                           (mb_y != mb_rows-1)) ? ((mb_y+1)*mb_cols + (mb_x  )) : (NO_MB);
	neighbors[p_ne] = ((mb_x != mb_cols-1)     && (mb_y != mb_rows-1)) ? ((mb_y+1)*mb_cols + (mb_x+1)) : (NO_MB);
	
	for(int i = 0; i < mb_list[thread_id].n; ++i)
	{
		int cp = mb_list[thread_id].p_idx[i];
		particles[cp].ax = particles[cp].ay = 0;
		
		// Collide with particles in my block
		for(int j = 0; j < mb_list[thread_id].n; ++j)
		{
			if(i != j) // do not collide with self
			{
				apply_force_gpu(particles[cp], particles[mb_list[thread_id].p_idx[j]]);
			}
		}
		
		// Collide with particles in neighboring blocks
		for(int k = 0; k < 8; ++k)
		{
			if(neighbors[k] != NO_MB) // Make sure I have a neighbor!
			{
				for(int j = 0; j < mb_list[neighbors[k]].n; ++j)
				{
					apply_force_gpu(particles[cp], particles[mb_list[neighbors[k]].p_idx[j]]);
				}
			}
		}
	}
}

__global__ void move_gpu (particle_t * particles, int n, double size)
{
  // Get thread (particle) ID
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= n) return;

  particle_t * p = &particles[tid];
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p->vx += p->ax * dt;
    p->vy += p->ay * dt;
    p->x  += p->vx * dt;
    p->y  += p->vy * dt;

    //
    //  bounce from walls
    //
    while( p->x < 0 || p->x > size )
    {
        p->x  = p->x < 0 ? -(p->x) : 2*size-p->x;
        p->vx = -(p->vx);
    }
    while( p->y < 0 || p->y > size )
    {
        p->y  = p->y < 0 ? -(p->y) : 2*size-p->y;
        p->vy = -(p->vy);
    }

}


