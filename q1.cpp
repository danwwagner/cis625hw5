#include <iostream>
#include <cmath>
#include <mpi.h>

#define R_SZ 0.00001

int num_threads = 1;

double f(int x, int y)
{
	return -std::cos(x) * std::sin(y) * std::exp(-((x - M_PI)*(x - M_PI) + (y - M_PI)*(y - M_PI)));	
}

double do_work(int rank, int x_start, int x_end, int y_start, int y_end)
{
	int iters = x_end - x_start;
	int loop_start = x_start + (rank * (iters / num_threads) + (iters % num_threads));
	int loop_end = x_start + (loop_start + ((iters / num_threads) + (iters % num_threads)));

	if (rank == num_threads-1) loop_end = x_end;

	std::cout << "rank: " << rank << " start: " << loop_start << " end: " << loop_end << std::endl;	

	double result = 0.0;
	for (int i = loop_start; i < loop_end; ++i)
	{
		for (int j = y_start; j < y_end; ++j)
		{
			result += R_SZ*R_SZ*f(x_start+i*R_SZ, y_start+j*R_SZ);	
		}
	}

	return result;
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
	if (argc != 5)
	{
		std::cout << "usage: " << argv[0] << " [x_start] [x_end] [y_start] [y_end]" << std::endl;
		return 1;
	}

	int x_start = std::atoi(argv[1]);
	int x_end = std::atoi(argv[2]);
	int y_start = std::atoi(argv[3]);
	int y_end = std::atoi(argv[4]);
    int rc;
	int rank = 0;
	double tstart, ttotal;

    if ((rc = MPI_Init(&argc, &argv)) != MPI_SUCCESS)
   	{
		std::cout << "error: cannot start mpi" << std::endl;
		MPI_Abort(MPI_COMM_WORLD, rc);	  
    } 

	
	MPI_Comm_size(MPI_COMM_WORLD,&num_threads);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	tstart = myclock();	
	tstart = myclock();	

	double result = do_work(rank, x_start, x_end, y_start, y_end);

	double global_result;
	MPI_Reduce(&result, &global_result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

   	ttotal = myclock() - tstart;

	if (rank == 0)
	{
		std::cout << "Result: " << global_result << std::endl;
		std::cout << "Time: " << ttotal << " seconds" << std::endl;
	}

	MPI_Finalize();
	return 0;
}
