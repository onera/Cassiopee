# - test de base sur mpi4py -
import KCore.test as test
try: from mpi4py import MPI
except: raise ImportError('FAILED to import mpi4py')

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
print(rank, size)
if rank == 0: test.testO([rank, size], 1)

try: import Converter.Mpi as Cmpi
except: raise ImportError('FAILED to import Cmpi')

rank = Cmpi.rank
size = Cmpi.size
print(rank, size)
if rank == 0: test.testO([rank, size], 2)
