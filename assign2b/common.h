#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

#include <stdio.h>

inline int mymin( int a, int b ) { return a < b ? a : b; }
inline int mymax( int a, int b ) { return a > b ? a : b; }

#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

#define NUM_THREADS 256

//
//  saving parameters
//
const int NSTEPS = 100;
const int SAVEFREQ = 10;

//
// particle data structure
//
typedef struct 
{
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;
} particle_t;

//================================================================================
//====================== Microblock structure ====================================
//================================================================================
const int max_particles_per_mb = 8;

//Possible values for valid array
#define INVALID 0
#define VALID 1
#define NEWLY_ADDED 2

typedef struct {
  int p_idx[max_particles_per_mb];
  int valid[max_particles_per_mb];
} microblock;

//================================================================================
//========================= Kernels ==============================================
//================================================================================
__global__ void distribute_gpu (microblock* mb_list, int mb_rows, int mb_cols, particle_t* particles, int n, double phys_size);
__device__ void apply_force_gpu(particle_t &particle, particle_t &neighbor);
__global__ void compute_forces_gpu (microblock* mb_list, int mb_rows, int mb_cols, particle_t* particles);
__global__ void move_gpu (particle_t * particles, int n, double size);
__global__ void migrate_particles_gpu (microblock* mb_list, int mb_rows, int mb_cols, particle_t* particles, int n, double phys_size);

#define NO_MB -1

#define p_sw 0
#define p_s  1
#define p_se 2
#define p_w  3
#define p_e  4
#define p_nw 5
#define p_n  6
#define p_ne 7

//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//

double set_size( int n );
void init_particles( int n, particle_t *p );
void apply_force( particle_t &particle, particle_t &neighbor );
void move( particle_t &p );

//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif
