#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

//================================================================================
//========================== Structures ==========================================
//================================================================================

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

//Collision Token
typedef struct {
  double position;
  int type;
  int particle_id;
} token;

//Active collisions
//List of [particle_id1 particle_id2] pairs
typedef struct {
  int id1;
  int id2;
} collision;

//Holds the information for a single partition
typedef struct {
  //Particles
  particle_t* particles;
  int num_particles;
  int max_particles;

  //Active IDs
  int* is_id_active;
  int* free_ids;
  int num_used_ids;

  //Ghost flags
  int* is_ghost;

  //Collision Tokens
  token* xtokens;
  token* ytokens;

  //Active Collisions
  collision* active_collisions;
  int num_active_collisions;

  //Collision table
  char* collision_table;
} partition_t;

//================================================================================
//====================== Collision Detector Interface ============================
//================================================================================
partition_t* alloc_partition(int max_particles);
void free_partition_t(partition_t* p);
int add_particle(partition_t* part, particle_t p);
void remove_particle(partition_t* p, int id);
void set_ghost(partition_t* p, int id, int is_ghost);
void set_state(partition_t* p, int id, double x, double y, double vx, double vy);
particle_t* get_particle(partition_t* p, int id);

void update_particles(partition_t* p);
void update_particles_mpi(partition_t* p, double left_x, double right_x, double bottom_y, double top_y);

//
//  saving parameters
//
const int NSTEPS = 100;
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
#define mass    0.01
#define min_r   (cutoff/100)

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
