#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include "common.h"

//================================================================================
//======================= Simulation Parameters ==================================
//================================================================================

double simsize;
FILE* fsave;

//================================================================================
//======================== Thread Barriers =======================================
//================================================================================
typedef struct {
    int needed;
    int called;
    pthread_mutex_t mutex;
    pthread_cond_t cond;
} barrier_t;

int barrier_init(barrier_t *barrier,int needed)
{
    barrier->needed = needed;
    barrier->called = 0;
    pthread_mutex_init(&barrier->mutex,NULL);
    pthread_cond_init(&barrier->cond,NULL);
    return 0;
}

int barrier_destroy(barrier_t *barrier)
{
    pthread_mutex_destroy(&barrier->mutex);
    pthread_cond_destroy(&barrier->cond);
    return 0;
}

int barrier_wait(barrier_t *barrier)
{
    pthread_mutex_lock(&barrier->mutex);
    barrier->called++;
    if (barrier->called == barrier->needed) {
        barrier->called = 0;
        pthread_cond_broadcast(&barrier->cond);
    } else {
        pthread_cond_wait(&barrier->cond,&barrier->mutex);
    }
    pthread_mutex_unlock(&barrier->mutex);
    return 0;
}

//================================================================================
//=========================== Threading ==========================================
//================================================================================

//Thread Barrier
barrier_t barrier;
barrier_t end_barrier;

//================================================================================
//======================== Region Manager ========================================
//================================================================================

//Region Identifiers
const int N = 0;
const int NE = 1;
const int E = 2;
const int SE = 3;
const int S = 4;
const int SW = 5;
const int W = 6;
const int NW = 7;

//Region Structure

typedef struct region {
  //Partition
  partition_t* partition;
  //Bounds
  double min_x,min_y,max_x,max_y;
  //Id mapping
  int* ids;
  //Send Queue
  char* send;  
  //Neighbours;
  struct region* neighbours[8];
} region;

//Allocate a new region
region* alloc_region(int max_particles, int total_particles, double min_x, double min_y, double max_x, double max_y){
  //Allocate region and set bounds
  region* r = (region*)malloc(sizeof(region));
  r->min_x = min_x;
  r->max_x = max_x;
  r->min_y = min_y;
  r->max_y = max_y;
  
  //Allocate partition
  r->partition = alloc_partition(max_particles);

  //Create Ids
  r->ids = (int*)malloc(total_particles * sizeof(int));
  for(int i=0; i<total_particles; i++)
    r->ids[i] = -1;

  //Create Send Queue
  r->send = (char*)malloc(max_particles * sizeof(char));
  for(int i=0; i<total_particles; i++)
    r->send[i] = 0;

  //Initialize Neighbours
  for(int i=0; i<8; i++)
    r->neighbours[i] = 0;

  //Return
  return r;
}

//Add a new particle to this region
void add_region_particle(region* r, particle_t* p){
  //Error if particle already exists
  if(r->ids[p->globalID] >= 0){
    printf("Particle %d already exists in region.\n", p->globalID);
    exit(-1);
  }

  //Add to partition
  int id = add_particle(r->partition, *p);
  
  //Create ID
  r->ids[p->globalID] = id;
}

//Remove a particle from this region
void remove_region_particle(region* r, particle_t* p){
  int id = r->ids[p->globalID];
  remove_particle(r->partition, id);
  r->ids[p->globalID] = -1;
}

//================================================================================
//======================== Region Threading ======================================
//================================================================================

int inline in_sending_bounds(region* r, particle_t* p){
  return
    p->y <= r->min_y + cutoff || p->y >= r->max_y - cutoff ||
    p->x <= r->min_x + cutoff || p->x >= r->max_x - cutoff;
}

int inline out_of_bounds(region* r, particle_t* p){
  return
    p->y < r->min_y - cutoff || p->y > r->max_y + cutoff ||
    p->x < r->min_x - cutoff || p->x > r->max_x + cutoff;
}

int inline in_ghost_bounds(region* r, particle_t* p){
  return
    p->y < r->min_y || p->y > r->max_y ||
    p->x < r->min_x || p->x > r->max_x;
}

int inline in_neighbour_bounds(region* r, particle_t* p, int dir){
  switch(dir){
  case N:
    return p->y <= r->min_y + 4*cutoff;
  case E:
    return p->x >= r->max_x - 4*cutoff;
  case S:
    return p->y >= r->max_y - 4*cutoff;
  case W:
    return p->x <= r->min_x + 4*cutoff;
  case NE:
    return in_neighbour_bounds(r, p, N) && in_neighbour_bounds(r, p, E);
  case SE:
    return in_neighbour_bounds(r, p, S) && in_neighbour_bounds(r, p, E);
  case NW:
    return in_neighbour_bounds(r, p, N) && in_neighbour_bounds(r, p, W);
  case SW:
    return in_neighbour_bounds(r, p, S) && in_neighbour_bounds(r, p, W);
  }
  
  printf("Invalid Direction %d\n", dir);
  exit(-1);
}

//Send determination
void determine_sends(region* r){
  partition_t* part = r->partition;
  for(int i=0; i<part->num_particles; i++){    
    if(part->is_id_active[i] && !part->is_ghost[i]) {      
      r->send[i] = r->send[i] || in_sending_bounds(r, &(part->particles[i]));
    } else
      r->send[i] = 0;
  }
}

//Send
void send_to_neighbour(region* r, int dir){
  partition_t* part = r->partition;
  for(int i=0; i<part->num_particles; i++)
    if(r->send[i]){
      //Destination Region
      region* dest = r->neighbours[dir];

      //Get particle
      particle_t* p = &(part->particles[i]);
      if(dest && in_neighbour_bounds(r, p, dir)){
        //Add/send
        int id = dest->ids[p->globalID];
        if(id < 0)
          add_region_particle(dest, p);
        else
          set_state(dest->partition, id, p->x, p->y, p->vx, p->vy);
      }
    }
}

//Remove/Ghost
void remove_and_ghost(region* r){
  partition_t* part = r->partition;
  for(int i=0; i<part->num_particles; i++){      
    if(part->is_id_active[i]){
      //Compute which particles are currently in bounds
      r->send[i] = r->send[i] && in_sending_bounds(r, &(part->particles[i]));

      //Either remove or set as ghost
      particle_t* p = &(part->particles[i]);      
      if(out_of_bounds(r, p))
        remove_region_particle(r, p);
      else
        set_ghost(part, i, in_ghost_bounds(r, p));
    }
  }
}

//================================================================================
//==================== Main Driver ===============================================
//================================================================================

//Thread Entry
void* thread_entry(void* arg){
  region* r = (region*)arg;

  for(int step=0; step<NSTEPS; step++){
    //Determine which particles to send
    determine_sends(r);
    barrier_wait(&barrier);

    //Send particles to neighbours
    for(int dir=0; dir<8; dir++){
      send_to_neighbour(r, dir);
      barrier_wait(&barrier);
    }

    //Update particles
    remove_and_ghost(r);
    update_particles(r->partition);
    barrier_wait(&end_barrier);
  }

  return 0;
}

double run_simulation(particle_t* ps, int n, int rows, int cols, FILE* fsave){
  //Regions
  int num_regions = rows*cols;
  region* regions[num_regions];

  //Allocate Regions
  for(int i=0; i<rows; i++)
    for(int j=0; j<cols; j++){
      //Find Boundaries
      double minx = (simsize / cols)*j;
      double maxx = (simsize / cols)*(j + 1);
      double miny = (simsize / rows)*i;
      double maxy = (simsize / rows)*(i + 1);

      //Allocate Region
      int max_particles = max(100, (int)(4 * n / num_regions));      
      regions[i+j*rows] = alloc_region(max_particles, n, minx, miny, maxx, maxy);

      //Add Particles to Region
      for(int k=0; k<n; k++) {
        particle_t* p = &ps[k];
        if(p->x >= minx && p->x <= maxx && p->y >= miny && p->y <= maxy)
          add_region_particle(regions[i+j*rows], p);
        remove_and_ghost(regions[i+j*rows]);
      }
    }

  //Connect Regions to Neighbours
  for(int i=0; i<rows; i++)
    for(int j=0; j<cols; j++){
      if(i>0 && j>0)                  regions[i+j*rows]->neighbours[NW] = regions[i-1+(j-1)*rows];
      if(i>0)                         regions[i+j*rows]->neighbours[N] = regions[i-1+j*rows];
      if(i>0 && j<cols-1)             regions[i+j*rows]->neighbours[NE] = regions[i-1+(j+1)*rows];

      if(i<rows-1 && j>0)             regions[i+j*rows]->neighbours[SW] = regions[(i+1)+(j-1)*rows];
      if(i<rows-1)                    regions[i+j*rows]->neighbours[S] = regions[(i+1)+j*rows];
      if(i<rows-1 && j<cols-1)        regions[i+j*rows]->neighbours[SE] = regions[(i+1)+(j+1)*rows];
            
      if(j>0)                         regions[i+j*rows]->neighbours[W] = regions[i+(j-1)*rows];
      if(j<cols-1)                    regions[i+j*rows]->neighbours[E] = regions[i+(j+1)*rows];
    }

  //Initialize Partitions
  for(int i=0; i<num_regions; i++)
    sweep_and_prune(regions[i]->partition);

  //Initialize Timer
  double simulation_time = read_timer();

  //Create Threads
  pthread_t threads[num_regions];
  barrier_init(&barrier, num_regions);
  barrier_init(&end_barrier, num_regions+1);
  for(int i=0; i<num_regions; i++)
    pthread_create(&threads[i], NULL, &thread_entry, regions[i]);

  //Start Processing
  for(int step=0; step<NSTEPS; step++){
    //Synchronize
    barrier_wait(&end_barrier);

    //Save
    if(fsave && (step%SAVEFREQ) == 0){
      for(int i=0; i<num_regions; i++){
        partition_t* p1 = regions[i]->partition;
        for(int i=0; i<p1->num_particles; i++)
          if(p1->is_id_active[i] && !p1->is_ghost[i])
            ps[p1->particles[i].globalID] = p1->particles[i];
      }
      save(fsave, n, ps);
    }    
  }

  //Calculate Time Elapsed
  simulation_time = read_timer() - simulation_time;
  return simulation_time;
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
    int num_threads = read_int(argc, argv, "-p", 4);

    // Determine how many rows and columns
    int cols, rows;
    for(int i = (int)sqrt(num_threads); i >= 1; i--){
      if(num_threads % i == 0){
        cols = i;
        rows = num_threads/i;
        break;
      }
    }    

    char *savename = read_string( argc, argv, "-o", NULL );

    //Open file
    fsave = savename ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    simsize = set_size( n );
    init_particles( n, particles );
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = run_simulation(particles, n, rows, cols, fsave);
    
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time);
    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
