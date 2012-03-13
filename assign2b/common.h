#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

inline int mymin( int a, int b ) { return a < b ? a : b; }
inline int mymax( int a, int b ) { return a > b ? a : b; }

#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

//
//  saving parameters
//
const int NSTEPS = 1000;
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
const int max_particles_per_mb = 10;

typedef struct {
  int n;
  int p_idx[max_particles_per_mb];
} microblock;

__global__ void distribute_gpu (microblock* mb_list, int mb_rows, int mb_cols, particle_t* particles, int n, double phys_size);
__device__ void apply_force_gpu(particle_t &particle, particle_t &neighbor);
__global__ void compute_forces_gpu (microblock* mb_list, int mb_rows, int mb_cols, particle_t* particles);
__global__ void move_gpu (particle_t * particles, int n, double size);

#define NO_MB -1

const int p_sw = 0;
const int p_s  = 1;
const int p_se = 2;
const int p_w  = 3;
const int p_e  = 4;
const int p_nw = 5;
const int p_n  = 6;
const int p_ne = 7;

//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//

void set_size( int n );
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
