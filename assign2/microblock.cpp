#include "microblock.h"

void mb_init(microblock* target)
{
	target->particles = (particle_t**) malloc(default_mbuf_depth * sizeof(particle_t*));
	target->max_particles = default_mbuf_depth;
	target->num_particles = 0;
}

void mb_expand_particle(microblock* microblock, int new_max)
{
	// Create the new particle queue and copy old addresses into it
	particle_t** new_storage = (particle_t**) malloc(new_max * sizeof(particle_t*));
	for(int i = 0; i < microblock->num_particles; ++i)
	{
		new_storage[i] = microblock->particles[i];
	}
	
	// Free old storage
	free(microblock->particles);
	
	// Use the new particle queue
	microblock->particles = new_storage;
	microblock->max_particles = new_max;
}

void mb_add_particle(microblock* microblock, particle_t* particle_addr)
{
	// Expand queue if needed
	if(microblock->num_particles == microblock->max_particles)
	{
		mb_expand_particle(microblock, microblock->max_particles + 2);
	}
	
	// Add the new particle
	microblock->particles[microblock->num_particles] = particle_addr;
	microblock->num_particles += 1;
}

void mb_rm_particle(microblock* microblock, int pos)
{
	// Remove by overwriting target with last array value, then decrementing
	microblock->particles[pos] = microblock->particles[microblock->num_particles-1];
	microblock->num_particles -= 1;
}

void mb_clear(microblock* microblock)
{
	microblock->num_particles = 0;
}