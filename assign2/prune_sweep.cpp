#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"

//================================================================================
//=========================== Global Constants ===================================
//================================================================================

//Particle Radius
double radius = cutoff/2;

//Token types
const int R = 1;
const int L = -1;

//================================================================================
//========================== Triangular Matrix Utilities =========================
//================================================================================

//Compute the number of elements in an upper triangular matrix excluding the diagonal
inline int num_triangle_elements(int n){
  return (n-1)*n/2;
}

//Compute the idx of the element at (row,col) in an upper triangular matrix excluding diagonal
inline int triangle_idx(int n, int row, int col){
  return n*row - (row+1)*(row+2)/2 + col;
}

//================================================================================
//========================= Partition Manager ====================================
//================================================================================

//Create a new partition
partition* alloc_partition(int max_particles){
  //Allocate partition structure
  partition* p = (partition*)malloc(sizeof(partition));

  //Allocate space for particles
  p->max_particles = max_particles;
  p->num_particles = 0;
  p->particles = (particle_t*)malloc(max_particles * sizeof(particle_t));

  //Allocate buffers
  p->xtokens = (token*)malloc(2 * max_particles * sizeof(token));
  p->ytokens = (token*)malloc(2 * max_particles * sizeof(token));

  //Allocate active collision list  
  p->num_active_collisions = 0;
  p->active_collisions = (collision*)malloc(100 * max_particles * sizeof(collision));

  //Allocate collision table
  int size = num_triangle_elements(max_particles);
  p->collision_table = (char*)malloc(size * sizeof(char));
  for(int i=0; i<size; i++)
    p->collision_table[i] = 0;

  //Active and Free IDs
  p->is_id_active = (int*)malloc(max_particles * sizeof(int));
  p->free_ids = (int*)malloc(max_particles * sizeof(int));
  for(int i=0; i<max_particles; i++){
    p->is_id_active[i] = 0;
    p->free_ids[i] = i;
  }
  p->num_used_ids = 0;

  //Allocate ghost flags
  p->is_ghost = (int*)malloc(max_particles * sizeof(int));
  for(int i=0; i<max_particles; i++)
    p->is_ghost[i] = 0;

  //Return new partition
  return p;
}

//Free a partition
void free_partition(partition* p){
  free(p->particles);
  free(p->is_id_active);
  free(p->free_ids);
  free(p->is_ghost);
  free(p->xtokens);
  free(p->ytokens);
  free(p->active_collisions);
  free(p->collision_table);
  free(p);
}

//================================================================================
//======================= Collision Detector =====================================
//================================================================================

//Register an active collision
void register_active_collision(partition* p, int id1, int id2){
  collision* c = &((p->active_collisions)[p->num_active_collisions]);
  p->num_active_collisions++;
  c->id1 = id1;
  c->id2 = id2;
}

//Mark id1 and id2 as intersecting
void mark_intersection(partition* p, int id1, int id2){
  int min_id = min(id1, id2);
  int max_id = max(id1, id2);
  
  int idx = triangle_idx(p->max_particles, min_id, max_id);
  char num_intersections = p->collision_table[idx];
  p->collision_table[idx] = num_intersections + 1;

  if(num_intersections == 1)
    register_active_collision(p, min_id, max_id);
}

//Unmark id1 and id2 as intersecting
void unmark_intersection(partition* p, int id1, int id2){
  int min_id = min(id1, id2);
  int max_id = max(id1, id2);
  
  int idx = triangle_idx(p->max_particles, min_id, max_id);
  p->collision_table[idx]--;
}

//Notify collision detector that token t1 and t2 has been swapped
void swap(partition* p, token t1, token t2){
  if(t1.type == R && t2.type == L)
    mark_intersection(p, t1.particle_id, t2.particle_id);
  else if(t1.type == L && t2.type == R)
    unmark_intersection(p, t1.particle_id, t2.particle_id);
}

//Sinking sort with collision detection logic
//Sinks element i down to its rightful position.
//Assumes list up to i is sorted.
inline void sweep_down(partition* p, token* tokens, int i){
  while(i>0) {
    token t1 = tokens[i-1];
    token t2 = tokens[i];
    
    if(t2.position < t1.position){
      //Swap
      tokens[i] = t1;
      tokens[i-1] = t2;
      //Notify swap
      swap(p, t1, t2);
      i--;
    } else{
      //Found rightful position
      break;
    }
  }
}

//Sweep through tokens from i = 1 ... num_tokens
//Sinking sort elements
void sweep(partition* p, token* tokens){
  for(int i=1; i<2 * p->num_particles; i++)
    sweep_down(p, tokens, i);
}

//Prune the active collisions
//Keep only the ones where num intersections is >= 2
void prune(partition* p){
  collision* dest = p->active_collisions;

  int collision_length = p->num_active_collisions;  
  for(int i=0; i<collision_length; i++){
    collision* c = &(p->active_collisions[i]);

    //Check whether ids are active
    int is_active = p->is_id_active[c->id1] && p->is_id_active[c->id2];

    //Check whether both ids are ghosts
    int is_ghost = p->is_ghost[c->id1] && p->is_ghost[c->id2];
    
    //Check whether collision is active
    int idx = triangle_idx(p->max_particles, c->id1, c->id2);
    char num_intersections = p->collision_table[idx];

    //If the collision is active then copy to dest
    if(is_active && !is_ghost && num_intersections == 2){
      if(dest != c)
        *dest = *c;
      dest++;
    }else{
      //Otherwise, discount this collision
      p->num_active_collisions--;
    }
  }
}

//Update the token information in the collision detector
void update_tokens(partition* part){
  //Update x tokens
  for(int i=0; i<2 * part->num_particles; i++){
    token* t = &(part->xtokens[i]);
    if(part->is_id_active[t->particle_id]){
      particle_t p = part->particles[t->particle_id];
      t->position = p.x + (t->type * radius);
    }
  }

  //Update y tokens
  for(int i=0; i<2 * part->num_particles; i++){
    token* t = &(part->ytokens[i]);
    if(part->is_id_active[t->particle_id]){
      particle_t p = part->particles[t->particle_id];
      t->position = p.y + (t->type * radius);
    }
  }
}

//Full Collision Detection Sweep
void sweep_and_prune(partition* p){
  update_tokens(p);
  sweep(p,p->xtokens);
  sweep(p,p->ytokens);
  prune(p);
}

//================================================================================
//====================== Collision Detector Interface ============================
//================================================================================

int add_particle(partition* part, particle_t p){
  //Safety check
  if(part->num_used_ids >= part->max_particles){
    printf("Can not add another particle. Maximum number of particles reached.\n");
    exit(-1);
  }
  
  //Get id
  int id = part->free_ids[part->num_used_ids];
  part->num_used_ids++;
  //Set active
  part->is_id_active[id] = 1;
  //Copy data
  part->particles[id] = p;

  //If we are using more ids than particles, then create more
  if(part->num_used_ids > part->num_particles){    
    //Create x tokens
    token* L_xtoken = &(part->xtokens[2 * part->num_particles]);
    L_xtoken->type = L;
    L_xtoken->particle_id = id;
  
    token* R_xtoken = &(part->xtokens[2 * part->num_particles + 1]);
    R_xtoken->type = R;
    R_xtoken->particle_id = id;

    //Create y tokens
    token* L_ytoken = &(part->ytokens[2 * part->num_particles]);
    L_ytoken->type = L;
    L_ytoken->particle_id = id;
  
    token* R_ytoken = &(part->ytokens[2 * part->num_particles + 1]);
    R_ytoken->type = R;
    R_ytoken->particle_id = id;

    //Increment num_particles
    part->num_particles++;

    //Clear Collision Table
    for(int i=0; i<part->num_particles; i++){
      if(i != id){
        int idx = triangle_idx(part->max_particles, min(i,id), max(i,id));
        part->collision_table[idx] = 0;
      }
    }
  }

  //Return id
  return id;
}

void ensure_active(partition* p, int id){
  if(!p->is_id_active[id]){
    printf("Particle %d is not active.\n", id);
    exit(-1);
  }
}

void remove_particle(partition* p, int id){
  ensure_active(p, id);

  //Set inactive
  p->is_id_active[id] = 0;
  //Return id
  p->num_used_ids--;
  p->free_ids[p->num_used_ids] = id;
}

void set_ghost(partition* p, int id, int is_ghost){
  ensure_active(p, id);
  p->is_ghost[id] = is_ghost;
}

void set_state(partition* p, int id, double x, double y, double vx, double vy){
  ensure_active(p, id);
  p->particles[id].x = x;
  p->particles[id].y = y;
  p->particles[id].vx = vx;
  p->particles[id].vy = vy;
}

particle_t* get_particle(partition* p, int id){
  ensure_active(p, id);
  return &(p->particles[id]);
}

//================================================================================
//================== Physics Calculations ========================================
//================================================================================

//Apply a force to p1, and an equal and opposing force (ala Newton's 3rd law) to p2
void apply_pairwise_force(particle_t* p1, particle_t* p2) {  
  double dx = p2->x - p1->x;
  double dy = p2->y - p1->y;
  double r2 = dx * dx + dy * dy;
  if(r2 < cutoff*cutoff) {
      
    //Limit the maximum force
    r2 = fmax( r2, min_r*min_r );
    double r = sqrt( r2 );

    //Repulsive Force
    double coef = ( 1 - cutoff / r ) / r2 / mass;

    //Propel both particles
    p1->ax += coef * dx;
    p1->ay += coef * dy;
    p2->ax -= coef * dx;
    p2->ay -= coef * dy;
  }
}

void update_particles(partition* p){
  //Calculate active collisions
  sweep_and_prune(p);

  //Reset acceleration
  for(int i=0; i<p->num_particles; i++)
    if(p->is_id_active[i])
      p->particles[i].ax = p->particles[i].ay = 0;

  //Accumulate acceleration
  for(int i=0; i<p->num_active_collisions; i++){
    collision c = p->active_collisions[i];
    apply_pairwise_force(&(p->particles[c.id1]), &(p->particles[c.id2]));
  }

  //Move Particles
  for(int i=0; i<p->num_particles; i++)
    if(p->is_id_active[i])
      move(p->particles[i]);
}
