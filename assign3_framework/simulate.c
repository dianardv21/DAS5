#include <stdio.h>
#include <stdlib.h>
#include "simulate.h"
#include <mpi.h>
#include <string.h>

#include <unistd.h> // for sleep()

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

    int numprocs, rank, req_count;
    double c = 0.15;
    double left, right; // halos

    // handles for comms
    MPI_Request reqs[6];
    MPI_Status stats[6];

    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // partition for start-end indices
    int start = 1, end;
    int jump = (i_max-2) / numprocs;
    int mod = (i_max-2) % numprocs;
    int edges[numprocs][2];

    for (int k = 0; k < numprocs; k++) {
        end = start + (jump - 1);
        if (mod) {end++; mod--;} // distribute remainder over chunks
        edges[k][0] = start;
        edges[k][1] = end;
        start = end + 1;
    }
    
    // determine process domain
    start = edges[rank][0];
    end = edges[rank][1];

    // start iterations
    for(int t = 0; t < t_max; t++) {
        
        req_count = 0; // reset counting

        // send/recv halo cells
        if (rank != 0) { // exclude leftmost process
            MPI_Isend(&current_array[start], 1, MPI_DOUBLE, rank-1,  rank, MPI_COMM_WORLD, &reqs[req_count++]); // send start to previous as end+1
            MPI_Irecv(&left, 1, MPI_DOUBLE, rank-1, rank-1, MPI_COMM_WORLD, &reqs[req_count++]); // get end from previous as start-1
        } else {left = 0;} // edge of array is always 0

        if (rank != numprocs-1) { // exclude rightmost process
            MPI_Isend(&current_array[end], 1, MPI_DOUBLE, rank+1,  rank, MPI_COMM_WORLD, &reqs[req_count++]); // send end to next as start-1
            MPI_Irecv(&right, 1, MPI_DOUBLE, rank+1, rank+1, MPI_COMM_WORLD, &reqs[req_count++]); // get start from next as end+1
        } else {right = 0;} // edge of array is always 0

        
        
        // let computation run during communication
        for(int i = start+1; i < end; i++) {
            
            next_array[i] = 2*current_array[i]-old_array[i]+c*(current_array[i-1]-(2*current_array[i]-current_array[i+1]));

        }
        
        // wait for comms and compute halo cells
        MPI_Waitall(req_count, reqs, stats);
        next_array[start] = 2*current_array[start]-old_array[start]+c*(left-(2*current_array[start]-current_array[start+1]));
        next_array[end] = 2*current_array[end]-old_array[end]+c*(current_array[end-1]-(2*current_array[end]-right));

        // swap locally
        double *temp = old_array;
        old_array = current_array;
        current_array = next_array;
        next_array = temp;

    }
   

    // collect results from other processes
    if (numprocs > 1) { // no comms necessary if only one process
        if(rank != 0) {
            // send current to master no need for non-blocking here
            MPI_Send(current_array, i_max, MPI_DOUBLE, 0,  rank, MPI_COMM_WORLD);
        }
        else {
            double buffer_array[i_max]; // buffer to store received array domains
            for (int i = 1; i < numprocs; i++) {
                // blocking receive data chunk, otherwise buffer_array gets overwritten
                MPI_Recv(&buffer_array, i_max, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &stats[1]);
                
                // for each non-master process get domain and copy only its computed domain to current_array
                start = edges[i][0];
                end = edges[i][1];
                
                // copy relevant part of buffer to relevant part of current_array
                memcpy(current_array + start, buffer_array + start, (end-start+1)*sizeof(double));
            }
        }
    }

MPI_Finalize();
return current_array;
    
}

// Hoeveel MPI request en stats moet ik aanmaken?
// Wat is de bedoeling van statuses en wat kan er misgaan met MPI_STATUS_IGNORE
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%% Blocking receive %%%%%%%%%%%%%%%%%%%%%%%%%% //
// %%%%%%%%%%%%%%%%%%%%%%%%%%% Blocking receive %%%%%%%%%%%%%%%%%%%%%%%%%% //
// %%%%%%%%%%%%%%%%%%%%%%%%%%% Blocking receive %%%%%%%%%%%%%%%%%%%%%%%%%% //
// %%%%%%%%%%%%%%%%%%%%%%%%%%% Blocking receive %%%%%%%%%%%%%%%%%%%%%%%%%% //

double *simulate(const int i_max, const int t_max, double *old_array,
        double *current_array, double *next_array)
{    

    int numprocs, rank, req_count;
    double c = 0.15;
    double left, right; // halos

    // handles for comms
    MPI_Request reqs[6];
    MPI_Status stats[6];

    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // partition for start-end indices
    int start = 1, end;
    int jump = (i_max-2) / numprocs;
    int mod = (i_max-2) % numprocs;
    int edges[numprocs][2];

    for (int k = 0; k < numprocs; k++) {
        end = start + (jump - 1);
        if (mod) {end++; mod--;} // distribute remainder over chunks
        edges[k][0] = start;
        edges[k][1] = end;
        start = end + 1;
    }
    
    // determine process domain
    start = edges[rank][0];
    end = edges[rank][1];

    // start iterations
    for(int t = 0; t < t_max; t++) {
        
        req_count = 0; // reset counting

        // send/recv halo cells
        if (rank != 0) { // exclude leftmost process
            MPI_Isend(&current_array[start], 1, MPI_DOUBLE, rank-1,  rank, MPI_COMM_WORLD, &reqs[req_count++]); // send start to previous as end+1
        } else {left = 0;} // edge of array is always 0

        if (rank != numprocs-1) { // exclude rightmost process
            MPI_Isend(&current_array[end], 1, MPI_DOUBLE, rank+1,  rank, MPI_COMM_WORLD, &reqs[req_count++]); // send end to next as start-1
        } else {right = 0;} // edge of array is always 0

        
        
        // let computation run during communication
        for(int i = start+1; i < end; i++) {
            
            next_array[i] = 2*current_array[i]-old_array[i]+c*(current_array[i-1]-(2*current_array[i]-current_array[i+1]));

        }
        
        // wait for comms and compute halo cells
        if (rank != 0) {
            MPI_Recv(&left, 1, MPI_DOUBLE, rank-1, rank-1, MPI_COMM_WORLD, &stats[1]);
        } else {left = 0;}
        if (rank != numprocs-1){
            MPI_Recv(&right, 1, MPI_DOUBLE, rank+1, rank+1, MPI_COMM_WORLD, &stats[2]);
        } else {right = 0;} 
        
        next_array[start] = 2*current_array[start]-old_array[start]+c*(left-(2*current_array[start]-current_array[start+1]));
        next_array[end] = 2*current_array[end]-old_array[end]+c*(current_array[end-1]-(2*current_array[end]-right));

        // swap locally
        double *temp = old_array;
        old_array = current_array;
        current_array = next_array;
        next_array = temp;

    }
   

    // collect results from other processes
    if (numprocs > 1) { // no comms necessary if only one process
        if(rank != 0) {
            // send current to master no need for non-blocking here
            MPI_Ssend(current_array, i_max, MPI_DOUBLE, 0,  rank, MPI_COMM_WORLD);
        }
        else {
            double buffer_array[i_max]; // buffer to store received array domains
            for (int i = 1; i < numprocs; i++) {
                // blocking receive data chunk, otherwise buffer_array gets overwritten
                MPI_Recv(&buffer_array, i_max, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &stats[1]);
                
                // for each non-master process get domain and copy only its computed domain to current_array
                start = edges[i][0];
                end = edges[i][1];
                
                // copy relevant part of buffer to relevant part of current_array
                printf("s:%i, e:%i\n", start, end);
                memcpy(current_array + start, buffer_array + start, (end-start+1)*sizeof(double));
            }
        }
    }

MPI_Finalize();
return current_array;
    
}