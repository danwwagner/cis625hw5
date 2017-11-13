#include <stdlib.h>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <float.h>
#include <cmath>
#include <mpi.h>
#include <omp.h>

// Floats for single precision decimal values
float end_x, start_x, min_x;

// Number of MPI processes and OMP threads
int nprocs, nthreads;
float *minimum;
float local_min;
void *eval_func(void*);
double myclock(void);

// Evaluate the function at points using Monte Carlo analysis
void *eval_func(void* rank) {
	// Segment the range based on rank
	int ID = *((int *) rank);
	float start;// = ((long) ID) * (start_x * 
	float end;// = 
	local_min = FLT_MAX;
	// Generate a random floating point number in the subdivided range.

}

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
	nprocs = 2; //atoi(getenv("NSLOTS"));

	// Data validation
	if (argc != 3) {
		std::cout << "Please enter the start and end points." << std::endl;
		return -1;
	}
	
	// Initialize the MPI process
	int rc = MPI_Init(&argc, &argv);
	if (rc != MPI_SUCCESS) {
		std::cout << "Error with MPI." << std::endl;
		MPI_Abort(MPI_COMM_WORLD, rc);
	}
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		
	if (rank == 0) {
		minimum = new float[nprocs];
		// Grab the starting and ending values of x	
		start_x = atof(argv[1]);
		end_x = atof(argv[2]);
		std::cout << "start, end: " << start_x << ", " << end_x << std::endl;
		tstart = myclock(); tstart = myclock();
	}


	// Broadcast the starting and ending points to processes
	MPI_Bcast(&start_x, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&end_x, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	
	// Evaulate the function to find minimum values
	eval_func(&rank);

	if (rank == 0) {
		float min = 0;
		MPI_Status stat;
		// Pass back the values to prevent idle processes 
		for (int i = 0; i < nprocs; i++) {
			MPI_Recv(&min, 1, MPI_FLOAT, i, 0, MPI_COMM_WORLD, &stat);
			minimum[i] = min; 
		}
		
		min = FLT_MAX;
		min_x = 0;
		// Find the global minimum
		for (int i = 0; i < nprocs; i++) {
			float x = minimum[i];
			float y = cos(x) + (pow(fabs(7.0 - x), 2.0/15.0)) + 2*(pow(fabs(5.0 - x), 4.0/35.0));
			if (y < min) { min = y; min_x = x;}	
		}
		ttotal = myclock() - tstart;
		printf("DATA, %f, %d, %d, %lf\n", min_x, nprocs, nthreads, ttotal);	
	}
	else {
		MPI_Send(&local_min, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);	
	}
	// Free memory
	if (rank == 0) delete[] minimum;
	MPI_Finalize();
	return 0;
}

