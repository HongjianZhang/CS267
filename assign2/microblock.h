#include omp.h

typedef struct microblock microblock

struct microblock
{
	// Particles
	particle_t* particles;
	int num_particles;
	int max_particles;
	
	// Border Data
	int neighbors[8];
	double left_x, right_x, bottom_y, top_y;
	
	// Sync Data
	omp_lock_t lock;
};
