#include “mpi.h”

// MYMPI_Bcast function parameters:
//
// void *buffer  INOUT : buffer address
// int count IN : buffer size
// MPI_Datatype datatype , // IN : datatype of entry
// int root , // IN : root process ( sender )
// MPI_Comm communicator // IN : commuicator

// Taking inspiration from https://hpc-tutorials.llnl.gov/mpi/non_blocking/

// What it gotta do: 
// The function should implement broadcast communication by means of point-to-point communication.
// Each MPI process of the given communicator is assumed to execute the broadcast function with the
// same message parameters and root process argument. The root process is supposed to send the
// content of its own buffer to all other processes of the communicator. Any non-root process uses the
// given buffer for receiving the corresponding data from the root process


// For your implementation of broadcast assume a 1-dimensional
// ring topology, where each node only has communication links with its two direct neighbours with
// While any MPI process may, nonetheless, send
// messages to any other MPI process, messages need to be routed through a number of intermediate
// nodes. Communication cost can be assumed to be linear in the number of nodes involved. Aim for
// an efficient implementation of your function on an (imaginary) ring network topology
// (circularly) increasing and decreasing MPI process ids. 

// -> prbably non roots have to wait first until the receive and then send
// -> smth to check if msg was sent? if there was space for the buffer? or maybe not cuz we dont create
// the buffer inside the function, we just pass it from main
// -> IDEA!!: root passes msg to both prev and next nodes => find index where we stop sending

int MYMPI_Bcast (void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm communicator){
    int size, rank, next, prev;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Request reqs[size];   // required variable for non-blocking calls
    MPI_Status stats[size]; // i dont have size as a parameter

    //Finding the left and right neighbours:
    prev = rank - 1;
    next = rank + 1;
    if (rank == 0)
        prev = size - 1;
    if (rank ==  (size - 1))
        next = 0;
    // Sending the messages:
    if (rank == root){
        MPI_Isend(buffer, count, datatype, next, 1, communicator, &reqs[rank]);

    }
    else{
        MPI_Recv(buffer, count, datatype, rank, 1, communicator, &stats[rank]);
        if(next != root)
            MPI_Isend(&buffer, count, datatype, next, 1, communicator, &reqs[rank]);
    }
    return MPI_SUCCESS;
}



main(int argc, char *argv[]){
    rc = MPI_Init(&argc, &argv);
    if(rc != MPI_SUCCESS){
        printf(“Unable to set up MPI \n”);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }



    MPI_Finalize();
}
