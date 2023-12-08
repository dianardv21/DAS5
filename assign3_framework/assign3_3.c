#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

// MYMPI_Bcast function parameters:
//
// void *buffer  INOUT : buffer address
// int count IN : buffer size
// MPI_Datatype datatype , // IN : datatype of entry
// int root , // IN : root process ( sender )
// MPI_Comm communicator // IN : commuicator

// Used lecture slides and https://hpc-tutorials.llnl.gov/mpi/non_blocking/


int MYMPI_Bcast (void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm communicator){
    int size, rank, next, prev;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Request reqsend, reqrec;  
    MPI_Status stat; 

    //Finding the left and right neighbours:
    prev = rank - 1;
    next = rank + 1;
    if (rank == 0)
        prev = size - 1;
    if (rank ==  (size - 1))
        next = 0;
    // Sending the messages:
    if (rank == root){ // Root just sends
        MPI_Isend(buffer, count, datatype, next, 1, communicator, &reqsend);

    }
    else{ // Every other node receives then sends to its next neighbour:
        MPI_Irecv(buffer, count, datatype, prev, 1, communicator, &reqrec);
        MPI_Wait(&reqrec, &stat); // Make sure the data is received
        MPI_Isend(buffer, count, datatype, next, 1, communicator, &reqsend);
    }

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
