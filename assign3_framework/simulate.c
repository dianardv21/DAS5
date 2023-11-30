#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "simulate.h"


/* Add any global variables you may need. */

/*
 * i_max: how many data points are on a single wave
 * t_max: how many iterations the simulation should run
 * old_array: array of size i_max filled with data for t-1
 * current_array: array of size i_max filled with data for t
 * next_array: array of size i_max. You should fill this with t+1
 */
double *simulate(const int i_max, const int t_max, double *old_array,
        double *current_array, double *next_array)
{   
    int size, comm_rank;
    int num_procs = 1;
    
    MPI_Init();
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    //partition the data in chunks


    // sync block: MPI_Ssend - nonblock: MPI_ISsend
    // buffered block: MPI_Bsend - nonblock: MPI_IBsend
    // ready mode block: MPI_Rsend - nonblock: MPI_IRsend
    // MPI_Recv receives any send mode messages
    // MPI_Barrier(comm) for all processes

    MPI_Finalize();
    return current_array;
}
