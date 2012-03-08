#ifndef __MICROBLOCK_H
#define __MICROBLOCK_H

#include <omp.h>
#include "common.h"

const int default_mbuf_depth = 6;
#define NO_MB 0

typedef struct microblock microblock;

struct microblock
{
	// Particles
	particle_t** particles; // store addresses to particles in some other master list
	int num_particles;
	int max_particles;
	
	// Border Data
	microblock* neighbors[8];
	double left_x, right_x, bottom_y, top_y;
	
	// Sync Data
	omp_lock_t lock;
};

void mb_init(microblock* target);
void mb_add_particle(microblock* microblock, particle_t* particle_addr);
void mb_expand_particle(microblock* microblock, int new_max);
void mb_rm_particle(microblock* microblock, int pos);
void mb_clear(microblock* microblock);

#endif

