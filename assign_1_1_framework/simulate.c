#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include "simulate.h"

// pre-declaration of thread,
void *thread(void *vargp);

double c = 0.15;

typedef struct {// data to send to each thread
    double *old_array, *current_array, *next_array;
    int id, start, end, i_max, t_max;
} t_data;

pthread_mutex_t lock;
pthread_barrier_t barrier;

double *simulate1(const int i_max, const int t_max, const int num_threads,
        		 double *old_array, double *current_array, double *next_array) 
{
	pthread_t *thd;
	t_data *data = malloc(sizeof(t_data)*num_threads);
	thd = (pthread_t*) malloc (sizeof(pthread_t)*num_threads);

	// initialize mutexes and thread barriers
	pthread_mutex_init(&lock, NULL);
	pthread_barrier_init(&barrier, NULL, num_threads);

	// account for uneven distribution of datapoints
	int mod = i_max % num_threads; // 0 if data is evenly distributed, else this gets added to first thread workload.
	int jump = i_max / num_threads; // amount of data per thread
	int start = 1;
	int end;
	for (int i=0;i<num_threads;i++) { // calc and pass start-end indices to each thread
		// get start-end indices, remove last array index from end
		end = start + (jump - 1) + mod;
		if (i == num_threads-1) end = end - 2;
		// build data for thread
		data[i].id = i;
		data[i].i_max = i_max;
		data[i].t_max = t_max;
		data[i].start = start;
		data[i].end = end;
		data[i].old_array = old_array;
		data[i].current_array = current_array; 
		data[i].next_array = next_array; 
		// prepare for start-end calculation of next thread
		start = end + 1;
		
		mod = 0; // only add mod for first thread.
		
		// create the thread
		pthread_create(&thd[i], NULL, *thread, (void*)&data[i]);
		
	}

	// join all threads back up again
	for (int i = 0; i < num_threads; i++) {
		pthread_join(thd[i], NULL);
	}
	
	// free all allocated memory 
	free(thd);
	free(data);

    // return the result
    return current_array;
}



void *thread(void *vargp)
{	
	// get data from passed struct
	int id = ((t_data*) vargp)->id;
	int i_max = ((t_data*) vargp)->i_max; //using this for memcpy
	int t_max = ((t_data*) vargp)->t_max;
    int start = ((t_data*) vargp)->start;
    int end = ((t_data*) vargp)->end;
	double *old_array = ((t_data*) vargp)->old_array;
	double *current_array = ((t_data*) vargp)->current_array;
	double *next_array = ((t_data*) vargp)->next_array;

	// print start-end to inspect indexes
	printf("start: %i - end: %i \n", start, end); 

	//compute all iterations
  	for (int t = 0; t<t_max; t++){
		//end +1 to include end index
      	for(int i = start; i<end+1; i++)
        {	
        	next_array[i] = 2*current_array[i]-old_array[i]+c*(current_array[i-1]-(2*current_array[i]-current_array[i+1]));
        }

		// race cond: wait for threads before buffer swap
  		pthread_barrier_wait(&barrier);

		// only let first thread swap buffers to avoid conflict
		if(id == 0) {
			memcpy(old_array, current_array, i_max * sizeof(double));
      		memcpy(current_array, next_array, i_max * sizeof(double));
		}

		// let threads wait again for thread_0 to finish swapping
		pthread_barrier_wait(&barrier);
	}
	
    return NULL;
}


double *simulate(const int i_max, const int t_max, const int num_threads,
        		 double *old_array, double *current_array, double *next_array)  {
    double c = 0.15;
    for(int t = 0; t < t_max; t++) {
        
        for(int i = 1; i < i_max-1; i++) {
            next_array[i] = 2*current_array[i]-old_array[i]+c*(current_array[i-1]-(2*current_array[i]-current_array[i+1]));
        }
        
        double *temp = old_array;
        old_array = current_array;
        current_array = next_array;
        next_array = temp;
    }
    
    return current_array;
}
