from mpi4py import MPI
COMM_WORLD = MPI.COMM_WORLD
KCOMM = COMM_WORLD
rank = KCOMM.rank
size = KCOMM.size
print rank, size
