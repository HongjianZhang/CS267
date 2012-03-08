#ifndef __PPILE_H
#define __PPILE_H

#include "common.h"

typedef struct {
  //Particles
  particle_t* particles;
  int num_particles;
  int max_particles;
  
} ppile;

void free_ppile(ppile* target);
ppile* alloc_ppile(int max_particles);
particle_t* add_particle(ppile* basket, particle_t ball);
void clear_particles(ppile* basket);

#endif

