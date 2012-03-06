#include plist.h


plist* alloc_plist(int max_particles)
{
	plist* pl = (plist*) malloc(sizeof(plist));

	pl->max_length = max_particles;
	pl->end_particle = 0;
	pl->particles = (particle_t*)malloc(max_particles * sizeof(particle_t));
	
	  //Active and Free IDs
	pl->is_id_active = (int*)malloc(max_particles * sizeof(int));
	pl->free_ids = (int*)malloc(max_particles * sizeof(int));
	for(int i = 0; i < max_particles; i++){
		pl->is_id_active[i] = 0;
		pl->free_ids[i] = i;
	}
	pl->num_used_ids = 0;
	
	return pl;
}

void free_plist(plist* target)
{
	free(target->particles);
	free(target->is_id_active);
	free(target->free_ids);
	
	free(target);
}

particle_t* add_particle(plist* basket, particle_t ball)
{
	if(basket->num_used_ids >= basket->max_particles)
	{
		printf("Can not add another particle. Maximum number of particles reached.\n");
		exit(-1);
	}

	//Get id
	int id = basket->free_ids[basket->num_used_ids];
	basket->num_used_ids++;
	//Set active
	basket->is_id_active[id] = 1;
	//Copy data
	basket->particles[id] = ball;
	
	if(part->num_used_ids > part->num_particles){    
		//Increment num_particles
		part->num_particles++;
	}

	//Return id
	return &(basket->particles[id]);
}

void rm_particle(plist* basket, particle_t* ball)
{
	// Get the particle id of interest
	int id = ball - basket->particles;
	
	// Sanity check the id
	if((id >= max_length) || (basket->is_id_active[id] == 0))
	{
		printf("RP Fault: Particle id invalid");
		exit(-1);
	}
	
	//Set inactive
	p->is_id_active[id] = 0;
	//Return id
	p->num_used_ids--;
	p->free_ids[p->num_used_ids] = id;
}