# - convertFile2PartialPyTreeFromPath (pyTree) -

# ------------------------------------------------------------------------
# 1/ Import section
import Converter.PyTree as C
from mpi4py import MPI

# ------------------------------------------------------------------------
# 2/ Initilise MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# ------------------------------------------------------------------------
# 3/ Prepare and set parameter
filenamein  = 'CaseU_C2_AxiTransBump.hdf'

path        = ['/Base/Zone1/GridCoordinates/CoordinateX',
               '/Base/Zone1/GridCoordinates/CoordinateY']

BegR = 128
NbE  = 12
# ------------------------------------------------------------------------
# 4/ Define filter
DataSpaceMMRY = [[0   ], [1], [NbE], [1]]
DataSpaceFILE = [[BegR], [1], [NbE], [1]]
DataSpaceGLOB = [[0]]

Filter = dict()
Filter[path[0]] = DataSpaceMMRY+DataSpaceFILE+DataSpaceGLOB
Filter[path[1]] = DataSpaceMMRY+DataSpaceFILE+DataSpaceGLOB
# ------------------------------------------------------------------------
# 5/ Partial read
t = C.convertFile2PartialPyTreeFromPath(filenamein, Filter, comm)

# ------------------------------------------------------------------------
# 6/ Print CGNS Tree
print('*'*23, t[0][1])
print('*'*23, t[1][1])
C.printTree(t)
