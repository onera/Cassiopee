#include <cmath>
#include <iostream>

#include "test/xmpi_t1.hpp"
#include "xmpi/context.hpp"

namespace xcore
{
	namespace
	{
		int basic_test()
		{
			xcore::communicator comm; // Communicateur global par defaut
			int rk = comm.rank;
			int np = comm.size;
			int rank, nbp;
#if defined(_MPI)
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			MPI_Comm_size(MPI_COMM_WORLD, &nbp );
#else
			rank = 0;
			nbp  = 1;
#endif
			if ( (rank != rk) or (nbp != np) ) return 0;
			return 1;
		}

		int ring_test()
		{
			const int tokenref = 31415;
			xcore::communicator comm; // Communicateur global par defaut
			int rk = comm.rank;
			int np = comm.size;

			int token;
			if (np > 1) {
			if ( rk == 0 ) { 
				token = tokenref;
				comm.send(token, 1);
				//xcore::status st = 
                comm.recv(token, np-1);
			} else
			{
				comm.recv(token,rk-1);
				comm.send(token,(rk+1)%np);
			}
			}
			else token = tokenref;
			if ( token != tokenref ) return 0;
			else return 1;
		}

		int ring_test_async()
		{
			const int tokenbase = 31415;
			xcore::communicator globCom;
			int rk = globCom.rank;
			int np = globCom.size;

			int token = tokenbase + rk;
			int tokenrcv;
			if (np > 1) {
			xcore::request rq = globCom.irecv(tokenrcv, (rk+np-1)%np);
			globCom.send(token, (rk+1)%np);
			rq.wait();
			} else tokenrcv = token;
			int tokensol = (rk == 0 ? token + np - 1 : token - 1);
			if ( tokenrcv != tokensol ) return 0;
			return 1;
		}

		static double dx_arctan(double x)
		{    
	    	return (1.0 / (1.0 + x*x));
		}


		int compute_pi()
		{
		   int n;
		   double PI25DT = 3.141592653589793238462643;
		   double mypi=0., h, pi=0., i, sum;
	       int myid, numprocs;

 	       xcore::communicator globCom;
	       numprocs = globCom.size;
	       myid     = globCom.rank;

	       globCom.barrier();
	       if ( myid == 0 ) n = 10000;
		   /* Share intervals with other processors */
    	   globCom.bcast(n, n, 0);

	       sum = 0.0;
           h   = 1.0/n;
    	   /* Compute and Sum the "Heights" of each bar of the integration */
	       for (i = myid+0.5; i < n; i += numprocs)
    	   {
           		sum += dx_arctan(i*h);
    	   }
    	   /* Multiply by the "Widths" of each bar and 4.0 (arctan(1)=Pi/4) */
    	   mypi = 4.0*h*sum;

	       /* Consolidate and Sum Results on rank 0 (default) */
    	   globCom.reduce(mypi, pi, xcore::sum);

	       if (myid == 0)
    	   {
    	   		double error = fabs(pi - PI25DT);
    	   		if ( error > h ) return 0;
    	   }
    	   return 1;
		}
    }
	PyObject* test_all( PyObject *self, PyObject *args )
    {
	    xcore::context context;
   		int results[4];
    	results[0] = basic_test();
	    results[1] = ring_test();
    	results[2] = ring_test_async();
    	results[3] = compute_pi();
    	if (not results[0]) std::cerr << "basic test for xmpi failed : " << results[0] << std::endl;
    	if (not results[1]) std::cerr << "ring test for xmpi failed : " << results[1] << std::endl;
    	if (not results[2]) std::cerr << "asynchronous ring test for xmpi failed : " << results[2] << std::endl;
    	if (not results[3]) std::cerr << "parallel pi computing test for xmpi failed : " << results[3] << std::endl;
    	long flag =  (results[0] & results[1] & results[2] & results[3]);
    	return PyBool_FromLong(long(flag));
	}
}
 
