#include "ppile.h"
#include <stdlib.h>

ppile* alloc_ppile(int max_particles)
{
	ppile* pl = (ppile*) malloc(sizeof(ppile));

	pl->max_particles = max_particles;
	pl->num_particles = 0;
	pl->particles = (particle_t*)malloc(max_particles * sizeof(particle_t));
	
	return pl;
}

void free_ppile(ppile* target)
{
	free(target->particles);
	free(target);
}

particle_t* add_particle(ppile* basket, particle_t ball)
{
	if(basket->num_particles >= basket->max_particles)
	{
		printf("Can not add another particle. Maximum number of particles reached.\n");
		exit(-1);
	}
	
	int new_id = basket->num_particles;
	
	basket->particles[new_id] = ball;
	basket->num_particles += 1;

	//Return id
	return &(basket->particles[new_id]);
}

void clear_particles(ppile* basket)
{
	basket->num_particles = 0;
}
