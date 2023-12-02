#include <stdio.h>
#include <stdlib.h>
#include "simulate.h"
#include <mpi.h>
#include <string.h>

/* Add any global variables you may need. */

/*
 * i_max: how many data points are on a single wave
 * t_max: how many iterations the simulation should run
 * old_array: array of size i_max filled with data for t-1
 * current_array: array of size i_max filled with data for t
 * next_array: array of size i_max. You should fill this with t+1
 */
double *simulate1(const int i_max, const int t_max, double *old_array,
        double *current_array, double *next_array)
{    

    int numprocs, rank;
    double c = 0.15;
    double left, right; // halos
    int req_count = 0;

    // handles for comms
    MPI_Request reqs[6];
    MPI_Status stats[6];

    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // partition for start-end indices
    int start = 1, end;
    int jump = i_max / numprocs;
    int mod = i_max % numprocs;
    int edges[numprocs][2];
    for (int k = 0; k < numprocs; k++) {
        end = start + (jump - 1) + mod;
        edges[k][0] = start;
        edges[k][1] = end;
        start = end + 1;
        mod = 0;
    }
    
    // determine process domain
    edges[numprocs-1][1] -= 2;
    start = edges[rank][0];
    end = edges[rank][1];
    //if(rank == numprocs-1) end = end - 2;
    printf("\n%i %i\n", start, end);

    // start iterations
    for(int t = 0; t < t_max; t++) {
        
        // send/recv halo cells, 
        if (rank != numprocs-1) {
            MPI_Isend(&current_array[end], 1, MPI_DOUBLE, rank+1,  rank, MPI_COMM_WORLD, &reqs[0]); // send end to next as start-1
            MPI_Irecv(&right, 1, MPI_DOUBLE, rank+1, rank+1, MPI_COMM_WORLD, &reqs[1]); // get start from next as end+1
            req_count += 2*(numprocs-2);
        } else {right = 0;} // edge of array is always 0
        if(rank != 0) {
            MPI_Isend(&current_array[start], 1, MPI_DOUBLE, rank-1,  rank, MPI_COMM_WORLD, &reqs[2]); // send start to previous as end+1
            MPI_Irecv(&left, 1, MPI_DOUBLE, rank-1, rank-1, MPI_COMM_WORLD, &reqs[3]); // get end from previous as start-1
            req_count += 2*(numprocs-2);
        } else {left = 0;} // edge of array is always 0
        
        // let computation run during communication
        for(int i = start+1; i < end; i++) {
            
            next_array[i] = 2*current_array[i]-old_array[i]+c*(current_array[i-1]-(2*current_array[i]-current_array[i+1]));

        }
        
        // wait for comms and compute halo cells
        MPI_Waitall(req_count, reqs, MPI_STATUS_IGNORE);
        next_array[start] = 2*current_array[start]-old_array[start]+c*(left-(2*current_array[start]-current_array[start+1]));
        next_array[end] = 2*current_array[end]-old_array[end]+c*(current_array[end-1]-(2*current_array[end]-right));
        
        if(rank == 1){
        for (int i=0;i<i_max;i++){
        printf("CURSADR: %f  r: %i  i: %i  \n", current_array[i], rank, i);
    }
        }

        // swap locally
        double *temp = old_array;
        old_array = current_array;
        current_array = next_array;
        next_array = temp;

    }
    
    printf("\n\n Proc: %i \n", rank);
    for (int i=0;i<i_max;i++){
    printf("CURR: %f  r: %i  i: %i  \n", current_array[i], rank, i);
    }
   



    if (numprocs > 1) { // no comms necessary if only one process
        if(rank != 0) {
            // send all current_arrays to master process
            double send_array[i_max];
            memcpy(send_array, current_array, i_max*sizeof(double));
            for (int j = 0; j<i_max;j++){
                    printf("send: %f, curr: %f\n", send_array[j], current_array[j]);
            }
            MPI_Isend(send_array, i_max, MPI_DOUBLE, 0,  rank, MPI_COMM_WORLD, &reqs[4]);
        }
        else {
            double buffer_array[i_max]; // buffer to store received array domains
            for (int i = 1; i < numprocs; i++) {
                // for each non-master process get domain and copy only its domain to current_array
                start = edges[i][0];
                end = edges[i][1];
                printf("\nstart: %i,  end: %i, rank: %i\n", start, end, rank);
                
                // receive current_array from other processes
                MPI_Recv(buffer_array, i_max, MPI_DOUBLE, i,  i, MPI_COMM_WORLD, &stats[5]);
                for (int j = 0; j<i_max;j++){
                    printf("curr: %f, buff: %f\n", current_array[j],buffer_array[j]);
                }
                // copy relevant part of buffer to relevant part of current_array
                memcpy(current_array + start, buffer_array + start, (end-start+1)*sizeof(double));
            }
        }
        
    }

return current_array;
MPI_Finalize();
    
}











double *simulate(const int i_max, const int t_max, double *old_array,
        double *current_array, double *next_array)
{    

    int numprocs, rank;
    double c = 0.15;
    double left, right; // halos
    int req_count = 0;

    // handles for comms
    MPI_Request reqs[6];
    MPI_Status stats[6];

    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // partition for start-end indices
    int start = 1, end;
    int jump = i_max / numprocs;
    int mod = i_max % numprocs;
    int edges[numprocs][2];
    for (int k = 0; k < numprocs; k++) {
        end = start + (jump - 1) + mod;
        edges[k][0] = start;
        edges[k][1] = end;
        start = end + 1;
        mod = 0;
    }

    // determine process domain
    edges[numprocs-1][1] -= 2;
    start = edges[rank][0];
    end = edges[rank][1];

    for (int i=0;i<i_max;i++){
        printf("BEFORECurr: %f  r: %i  i: %i  \n", current_array[i], rank, i);
    }

    printf("\ns: %i -- e: %i\n", start, end);
    current_array[start] = 15;
    current_array[end] = 9;
    // send/recv halo cells, 
    if (rank != numprocs-1) {
        MPI_Isend(&current_array[end], 1, MPI_DOUBLE, rank+1,  rank, MPI_COMM_WORLD, &reqs[0]); // send end to next as start-1
        MPI_Irecv(&right, 1, MPI_DOUBLE, rank+1, rank+1, MPI_COMM_WORLD, &reqs[1]); // get start from next as end+1
        req_count += 2*(numprocs-2);
        printf("\nr: %f with rank %i\n", right, rank);
    } else {right = 0;} // edge of array is always 0
    if(rank != 0) {
        MPI_Isend(&current_array[start], 1, MPI_DOUBLE, rank-1,  rank, MPI_COMM_WORLD, &reqs[2]); // send start to previous as end+1
        MPI_Irecv(&left, 1, MPI_DOUBLE, rank-1, rank-1, MPI_COMM_WORLD, &reqs[3]); // get end from previous as start-1
        req_count += 2*(numprocs-2);
        printf("\nl :%f with rank %i\n", left, rank);
    } else {left = 0;} // edge of array is always 0
    
    // let computation run during communication
    for(int i = start+1; i < end; i++) {
        
        next_array[i] = 2*current_array[i]-old_array[i]+c*(current_array[i-1]-(2*current_array[i]-current_array[i+1]));

    }
    
    // wait for comms and compute halo cells
    MPI_Waitall(req_count, reqs, MPI_STATUS_IGNORE);
    next_array[start] = 2*current_array[start]-old_array[start]+c*(left-(2*current_array[start]-current_array[start+1]));
    next_array[end] = 2*current_array[end]-old_array[end]+c*(current_array[end-1]-(2*current_array[end]-right));

    double *temp = old_array;
    old_array = current_array;
    current_array = next_array;
    next_array = temp;

    
    if(rank == 1){
        for (int i=0;i<i_max;i++){
        printf("CURSADR: %f  r: %i  i: %i  \n", current_array[i], rank, i);
    }}
    //if(rank == 0){
    //    for (int i=0;i<i_max;i++){
    //    printf("CURR: %f  r: %i  i: %i  \n", current_array[i], rank, i);
    //}}

    

MPI_Finalize();
return current_array;
    
}
