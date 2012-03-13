//
// Responsible for scanning through the particle list and adding the appropriate particle indices to each microblock.
// Don't forget to set microblock.n also.
//
__global__ void distribute_gpu (microblock* mb_list, int mb_rows, int mb_cols, particle_t* particles, int n, double phys_size)
{
    // Get X and Y co-ordinate of microblock
    int mb_x = threadIdx.x + blockIdx.x * blockDim.x;
    int mb_y = threadIdx.y + blockIdx.y * blockDim.y;
    int thread_id = mb_y * mb_cols + mb_x;

    // Make sure that we have a valid microblock to process
    if((mb_x >= mb_cols) || (mb_y >= mb_rows))
        return;

    // Determine boundaries of microblock
    double mb_width = phys_size / mb_cols;
    double mb_height = phys_size / mb_rows;
    double left_x = mb_x * mb_width;
    double right_x = (mb_x + 1) * mb_width;
    double bottom_y = mb_y * mb_height;
    double top_y = (mb_y + 1) * mb_height;

    // Iterate through array of particles
    for(int i = 0; i < n; i++)
    {
        // If particle lies within this microblock, add it
        if( (particles[i].x >= left_x) && (particles[i].x < right_x) &&
            (particles[i].y >= bottom_y) && (particles[i].y < top_y) )
        {
            mb_list[thread_id].p_idx[mb_list[thread_id].n] = i;
            mb_list[thread_id].n += 1;
        }
    }
}


