#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <mpi.h>

#define R_MIN 0
#define R_MAX 1000
#define N 1000000000

int num_threads = 1;

void do_work(int rank)
{
	int n_per_proc = N / num_threads;

	if (rank == num_threads-1) n_per_proc += N - (n_per_proc*(num_threads-1));

	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_int_distribution<int> dist(R_MIN, R_MAX);

	// Generate sample space.
	std::vector<int> sample_space;
	for (int i = 0; i < n_per_proc; ++i)
	{
		auto r = dist(mt);
		sample_space.push_back(r);
	}	

	auto mean = std::accumulate(sample_space.begin(), sample_space.end(), 0LL) / sample_space.size();

	// send mean back to server so it can be incorporated with other nodes' results.
	MPI_Send(&mean, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

	// receive new mean from server.
	MPI_Recv(&mean, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	int64_t sum_squares = 0;
	for (int i = 0; i < n_per_proc; ++i)
	{
		sum_squares += (sample_space[i] - mean)*(sample_space[i] - mean);
	}

	MPI_Send(&sum_squares, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD);
}

double myclock() {
   static time_t t_start = 0;  // Save and subtract off each time

   struct timespec ts;
   clock_gettime(CLOCK_REALTIME, &ts);
   if( t_start == 0 ) t_start = ts.tv_sec;

   return (double) (ts.tv_sec - t_start) + ts.tv_nsec * 1.0e-9;
}

int main(int argc, char** argv)
{
    int rc;
	int rank;
	double tstart, ttotal;
    if ((rc = MPI_Init(&argc, &argv)) != MPI_SUCCESS)
   	{
		std::cout << "error: cannot start mpi" << std::endl;
		MPI_Abort(MPI_COMM_WORLD, rc);	  
    } 

	MPI_Comm_size(MPI_COMM_WORLD,&num_threads);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
	if (rank == 0)
	{
		tstart = myclock();	
		tstart = myclock();	

		uint64_t sum_squares = 0;
		std::vector<int> means;
		for (int i = 1; i < num_threads; ++i)
		{
			int mean;
			MPI_Recv(&mean, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			means.push_back(mean);
		}

		auto new_mean = std::accumulate(means.begin(), means.end(), 0LL) / means.size();

		for (int i = 1; i < num_threads; ++i)
		{
			MPI_Send(&new_mean, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		}

		for (int i = 1; i < num_threads; ++i)
		{
			int64_t result;
			MPI_Recv(&result, 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			sum_squares += result;	
		}

		// Compute standard deviation.
		auto sd = std::sqrt(sum_squares / (N-1));
   		ttotal = myclock() - tstart;
		std::cout << "Standard deviation: " << sd << std::endl;
		std::cout << "Time: " << ttotal << " seconds" << std::endl;
	}
	else
	{
		do_work(rank);
	}

	MPI_Finalize();
	return 0;
}
