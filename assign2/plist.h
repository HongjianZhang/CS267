typedef struct {
  //Particles
  particle_t* particles;
  int end_particle;
  int max_length;

  //Active IDs
  int* is_id_active;
  int* free_ids;
  int num_used_ids;
  
} plist;

plist* alloc_plist(int max_particles);
void free_plist(plist* target);
particle_t* add_particle(plist* basket, particle_t ball);
void rm_particle(plist* basket, particle_t* ball);