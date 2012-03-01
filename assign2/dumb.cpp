#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"

//================================================================================
//========================= Partition Manager ====================================
//================================================================================

//Create a new partition
partition_t* alloc_partition(int max_particles){
  //Allocate partition structure
  partition_t* p = (partition_t*)malloc(sizeof(partition_t));

  //Allocate space for particles
  p->max_particles = max_particles;
  p->num_particles = 0;
  p->particles = (particle_t*)malloc(max_particles * sizeof(particle_t));

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
void free_partition(partition_t* p){
  free(p->particles);
  free(p->is_id_active);
  free(p->free_ids);
  free(p->is_ghost);
  free(p);
}

//================================================================================
//====================== Collision Detector Interface ============================
//================================================================================

int add_particle(partition_t* part, particle_t p){
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
  part->is_ghost[id] = 0;
  //Copy data
  part->particles[id] = p;

  //If we are using more ids than particles, then create more
  if(part->num_used_ids > part->num_particles){    
    //Increment num_particles
    part->num_particles++;
  }

  //Return id
  return id;
}

void ensure_active(partition_t* p, int id){
  if(!p->is_id_active[id]){
    printf("Particle %d is not active.\n", id);
    exit(-1);
  }
}

void remove_particle(partition_t* p, int id){
  ensure_active(p, id);

  //Set inactive
  p->is_id_active[id] = 0;
  //Return id
  p->num_used_ids--;
  p->free_ids[p->num_used_ids] = id;
}

void set_ghost(partition_t* p, int id, int is_ghost){
  ensure_active(p, id);
  p->is_ghost[id] = is_ghost;
}

void set_state(partition_t* p, int id, double x, double y, double vx, double vy){
  ensure_active(p, id);
  p->particles[id].x = x;
  p->particles[id].y = y;
  p->particles[id].vx = vx;
  p->particles[id].vy = vy;
}

particle_t* get_particle(partition_t* p, int id){
  ensure_active(p, id);
  return &(p->particles[id]);
}

//================================================================================
//================== Physics Calculations ========================================
//================================================================================

void update_particles(partition_t* p){

  //Calculate forces
  for(int i=0; i<p->num_particles; i++)
  {
    if(p->is_id_active[i] && !p->is_ghost[i])
	{
      p->particles[i].ax = p->particles[i].ay = 0;
	  
	  for(int j=0; j<p->num_particles; j++)
	  {
        if(p->is_id_active[j])
		{
	      apply_force(p->particles[i],p->particles[j]);
		}
	  }
    }
  }

  //Move Particles
  for(int i=0; i<p->num_particles; i++)
  {
    if(p->is_id_active[i] && !p->is_ghost[i])
	{
      move(p->particles[i]);
	}
  }
}
