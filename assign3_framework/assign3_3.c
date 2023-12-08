#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

// MYMPI_Bcast function parameters:
//
// void *buffer  INOUT : buffer address
// int count IN : buffer size
// MPI_Datatype datatype , // IN : datatype of entry
// int root , // IN : root process ( sender )
// MPI_Comm communicator // IN : commuicator

// Used lecture slides and https://hpc-tutorials.llnl.gov/mpi/non_blocking/ link from Canvas
// => do clockwise and counter-clockwise implementation

int MYMPI_Bcast (void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm communicator){
    int size, rank, next, prev, last;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status stat; 

    //Finding the left and right neighbours:
    prev = rank - 1;
    next = rank + 1;
    if (rank == 0){
        prev = size - 1;
    }
    if (rank ==  (size - 1)){
        next = 0;
    }
    //Finding last node that has to receive:
    if (root == 0)
        last = size - 1;
    else last = (root-1);

    //double start = MPI_Wtime();
    // Sending the messages:
    if (rank == root){ // Root just sends
        MPI_Send(buffer, count, datatype, next, 1, communicator);
        printf("Node %d (root) started broadcasting. Message: %d\n", rank, *(int*)buffer);
    }
    else if(rank == last){ //the last node just receives
        MPI_Recv(buffer, count, datatype, prev, 1, communicator, &stat);     
        printf("Node %d received message %d\n", rank, *(int*)buffer);

    }
    else{ // Every other node receives then sends to its right neighbour:
        MPI_Recv(buffer, count, datatype, prev, 1, communicator, &stat);     
        printf("Node %d received message %d\n", rank, *(int*)buffer);

        MPI_Send(buffer, count, datatype, next, 1, communicator);
        printf("Node %d sent message\n", rank);
    }

    //double end = MPI_Wtime();
    //printf("The process %d took %f seconds to run.\n", rank, (end - start));
    return MPI_SUCCESS;
}



int main(int argc, char *argv[]){
    int rc;
    rc = MPI_Init(&argc, &argv);
    if (rc != MPI_SUCCESS) {
        printf("Unable to set up MPI \n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    int test = 767;
    MYMPI_Bcast(&test, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
