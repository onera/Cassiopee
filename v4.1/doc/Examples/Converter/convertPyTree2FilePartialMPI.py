# - convertPyTree2FilePartial (pyTree) -

# ------------------------------------------------------------------------
# 1/ Import section
import Converter.PyTree   as C
import Converter.Internal as I

from   mpi4py             import MPI
import numpy
import Distribution       as     DIST

# ------------------------------------------------------------------------
# 2/ Initialise MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# ------------------------------------------------------------------------
# 3/ Prepare and set parameter
# filenamein  = 'CaseU_C2_AxiTransBump.hdf'
filenameout = 'Test_T3.hdf'

path        = ['/SuperName/', '/SuperName/SuperGirl/']
NbVtx       = 1000

# > Prepare
distribVtx     = numpy.empty( (size + 1), order='F', dtype='int32')
sVtx, rVtx     = DIST.computeStepAndReminder(NbVtx, size)

# > Compute Distribution
DIST.computeDistribution(distribVtx, sVtx, rVtx)

# > Compute NbEntry
NbE = distribVtx[rank+1] - distribVtx[rank]

# > Create numpy array
CoordX = numpy.ones(NbE, order='F', dtype=numpy.float64)*rank+1
CoordY = numpy.ones(NbE, order='F', dtype=numpy.int32)*rank

# ------------------------------------------------------------------------
# 4/ Define filter
DataSpaceMMRY = [[0               ], [1], [NbE], [1]]
DataSpaceFILE = [[distribVtx[rank]], [1], [NbE], [1]]
DataSpaceGLOB = [[NbVtx           ]]

Filter = dict()
Filter[path[0]] = DataSpaceMMRY+DataSpaceFILE+DataSpaceGLOB
Filter[path[1]] = DataSpaceMMRY+DataSpaceFILE+DataSpaceGLOB
if rank == 0:
    print(rank, Filter)
comm.barrier()
if rank == 1:
    print(rank, Filter)
comm.barrier()

# ------------------------------------------------------------------------
# 5/ Partial write
dt = I.newCGNSTree()
dn = I.addChild(dt, ['SuperName', CoordX, [], 'DataArray_t'])
I._addChild(    dn, ['SuperGirl', CoordY, [], 'DataArray_t'])
C.printTree(dt)

st = I.newCGNSTree()
sn = I.addChild(st, ['SuperName', None, [], 'DataArray_t'])
I._addChild(sn, ['SuperGirl', None, [], 'DataArray_t'])
C.printTree(st)

if rank == 0: C.convertPyTree2File(st, filenameout)
comm.barrier()

C.convertPyTree2FilePartial(dt, filenameout, comm, Filter)

# ------------------------------------------------------------------------
# 6/ Print CGNS Tree
# print '*'*23, t[1][1]
# C.printTree(t)
