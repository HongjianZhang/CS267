#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;

const int p_sw = 0;
const int p_s  = 1;
const int p_se = 2;
const int p_w  = 3;
const int p_e  = 4;
const int p_nw = 5;
const int p_n  = 6;
const int p_ne = 7;

const int P_SW = 0;
const int P_S  = 1;
const int P_SE = 2;
const int P_W  = 3;
const int P_E  = 4;
const int P_NW = 5;
const int P_N  = 6;
const int P_NE = 7;

const int NONE = -1;

#define VALID   1
#define INVALID 0

// moved from common.cpp
#define cutoff  0.01

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
  int globalID;
} particle_t;

//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
double set_size( int n );
void init_particles( int n, particle_t *p);
void apply_force( particle_t &particle, particle_t &neighbor );
void move( particle_t &p );

//
//  particle array management routines
//
int add_particle(particle_t &new_particle, int array_sz, particle_t *particles, char* p_valid);
void remove_particle(int index, char* p_valid);

//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p);

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif
