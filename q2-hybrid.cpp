#include <stdlib.h>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <float.h>
#include <cmath>
#include <mpi.h>
#include <omp.h>

// Number of iterations in Monte Carlo
#define NUM_ITER 1000000

// Floats for single precision decimal values
float end_x, start_x, min_x;

// Number of MPI processes and OMP threads
int nprocs, nthreads;
float *minimum;
float local_min;
void *eval_func(void*);
double myclock(void);
float f(float);

// Function to minimize via Monte Carlo.
float f(float x) {
	return cos(x) + (pow(fabs(7.0 - x), 2.0/15.0)) + 2*(pow(fabs(5.0 - x), 4.0/35.0));
}

// Evaluate the function at points using Monte Carlo analysis.
void *eval_func(void* rank) {
	// Segment the range based on rank.
	int ID = *((int *) rank);
	float start = (ID * ((end_x - start_x) / nprocs)) + start_x; 
	float end = start + ((end_x - start_x) / nprocs); 
	local_min = FLT_MAX;
	float f_min = FLT_MAX;
	float random, min;
	int x;

	omp_set_dynamic(0);
	omp_set_num_threads(nthreads);
	#pragma omp parallel for shared(f_min, local_min) private(x, min, random)
	// Generate a random floating point number in the subdivided range.
	for (x = 0; x < NUM_ITER; x++) {
		random = (rand() * (end - start) / RAND_MAX) + start;
		min = f(random);
		#pragma omp critical // Getting wrong values 2+ threads
		{
			if (min < f_min) {
				f_min = min;
				local_min = random;
			}
		}
	}

	// Send the value that minimizes the function back to the main process.
	if (ID != 0) MPI_Send(&local_min, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
	else minimum[0] = local_min;
}

// Used to compute the amount of time the program takes to run.
double myclock() {
   static time_t t_start = 0;  // Save and subtract off each time

   struct timespec ts;
   clock_gettime(CLOCK_REALTIME, &ts);
   if( t_start == 0 ) t_start = ts.tv_sec;

   return (double) (ts.tv_sec - t_start) + ts.tv_nsec * 1.0e-9;
}

int main(int argc, char *argv[]) {
	int rank;
	double tstart, ttotal;
	nprocs = atoi(getenv("NSLOTS"));

	// Data validation.
	if (argc < 3) {
		std::cout << "Please enter the start and end points." << std::endl;
		return -1;
	}
	
	// Initialize the MPI processes.
	int rc = MPI_Init(&argc, &argv);
	if (rc != MPI_SUCCESS) {
		std::cout << "Error with MPI." << std::endl;
		MPI_Abort(MPI_COMM_WORLD, rc);
	}
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		
	if (rank == 0) {
		minimum = new float[nprocs];
		// Grab the starting and ending values of x.
		start_x = atof(argv[1]);
		end_x = atof(argv[2]);
		std::cout << "start, end: " << start_x << ", " << end_x << std::endl;
		tstart = myclock(); tstart = myclock();
	}

	// Grab the number of threads per process from the command line.
	nthreads = atoi(argv[3]);

	// Broadcast the starting and ending points to processes.
	MPI_Bcast(&start_x, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&end_x, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	
	// Evaulate the function to find minimum values.
	eval_func(&rank);

	if (rank == 0) {
		float min = 0;
		MPI_Status stat;
		// Pass back the values to prevent idle processes. 
		for (int i = 1; i < nprocs; i++) {
			MPI_Recv(&min, 1, MPI_FLOAT, i, 0, MPI_COMM_WORLD, &stat);
			minimum[i] = min; 
		}
		
		min = FLT_MAX;
		min_x = FLT_MAX;
		// Find the global minimum.
		for (int i = 0; i < nprocs; i++) {
			float x = minimum[i];
			float y = f(x); 
			if (y < min) { min = y; min_x = x;}	
		}
		ttotal = myclock() - tstart;
		printf("DATA, %f, %d, %d, %lf\n", min_x, nprocs, nthreads, ttotal);	
	}

	// Free memory and clean up.
	if (rank == 0) delete[] minimum;
	MPI_Finalize();
	return 0;
}

