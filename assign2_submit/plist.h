#ifndef __PLIST_H
#define __PLIST_H

#include "common.h"

typedef struct {
  //Particles
  particle_t* particles;
  int end_particle;
  int max_length;

  //Active IDs
  int* is_id_active;
  int* free_ids;
  int num_used_ids; // also doubles as number of particles in system!
  
} plist;

void free_plist(plist* target);
plist* alloc_plist(int max_particles);
particle_t* add_particle(plist* basket, particle_t ball);
void rm_particle(plist* basket, particle_t* ball);

#endif

