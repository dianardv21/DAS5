#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "mpi.h"
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
    int numprocs, rank;
    double c = 0.15;

    typedef struct MPI_STATUS {
        int MPI_TAG;
        int MPI_ERROR;
        int MPI_SOURCE;
    }

    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // TODO: determine data partitioning
    int start = 1, end = 998;

    // start iterations
    for(int t = 0; t < t_max; t++) {

        // send/recv halo cells, 
        double left, right;
        MPI_STATUS status;

        if (rank != numprocs-1) {
            MPI_Isend((void*)current_array[end], 1, MPI_DOUBLE, rank+1,  rank, MPI_COMM_WORLD, *status); // send end to next as start-1
            MPI_Irecv(&right, 1, MPI_DOUBLE, rank+1, rank+1, MPI_COMM_WORLD, &status); // get start from next as end+1
        } else {right = 0;} // edge of array is always 0
        if(rank != 0) {
            MPI_Isend((void*)current_array[start], 1, MPI_DOUBLE, rank-1,  rank, MPI_COMM_WORLD, *status); // send start to previous as end+1
            MPI_Irecv(&left, 1, MPI_DOUBLE, rank-1, rank-1, MPI_COMM_WORLD, &status); // get end from previous as start-1
        } else {left = 0;} // edge of array is always 0
        
        // let computation run during communication
        for(int i = start+1; i < end; i++) {
            
            next_array[i] = 2*current_array[i]-old_array[i]+c*(current_array[i-1]-(2*current_array[i]-current_array[i+1]));

        }
        // TODO: Make this MPI_Wait with request instead of barrier
        MPI_Barrier(MPI_COMM_WORLD);
        // handle halo cells
        next_array[start] = 2*current_array[start]-old_array[start]+c*(left-(2*current_array[start]-current_array[start+1]));
        next_array[end] = 2*current_array[end]-old_array[end]+c*(current_array[end-1]-(2*current_array[end]-right));
        
        // swap locally
        double *temp = old_array;
        old_array = current_array;
        current_array = next_array;
        next_array = temp;

    }

    


    // MPI_ANY_TAG
    // sync block: MPI_Ssend - nonblock: MPI_ISsend
    // buffered block: MPI_Bsend - nonblock: MPI_IBsend
    // ready mode block: MPI_Rsend - nonblock: MPI_IRsend
    // MPI_Recv receives any send mode messages
    // MPI_Barrier(comm) for all processes

    MPI_Finalize();
    return current_array;
}
