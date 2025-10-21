# - adaptMesh_Exit (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import XCore.PyTree as XC
import Converter.Internal as Internal
import KCore.test as test
import numpy
from mpi4py import MPI # for MPI_Init

# 2D QUAD
a = G.cartHexa((0,0,0),(0.1,0.1,0.1),(11,11,2))
a = C.convertArray2NGon(a)
C._fillEmptyBCWith(a, 'nref', 'BCFarfield', dim=2)
Internal._adaptNGon32NGon4(a)

normal2D = numpy.array([0.0,0.0,1.0], dtype=numpy.float64)
ngonelts = Internal.getNGonNode(a)
ER = Internal.getNodeFromName(ngonelts, 'ElementRange')[1]
nfaces = ER[1]
nfaceselts = Internal.getNFaceNode(a)
ER = Internal.getNodeFromName(nfaceselts, 'ElementRange')[1]
ncells = ER[1]-ER[0]+1
gcells = numpy.arange(0, ncells)
gfaces = numpy.arange(1, nfaces+1)

for i in range(2):
    test.printMem("Creating hook on mesh adaptation")
    AM = XC.AdaptMesh_Init(a, normal2D, comm=[], gcells=gcells, gfaces=gfaces)
    XC.AdaptMesh_Exit(AM)
    test.printMem("Releasing hook on mesh adaptation.")
print("done.")
