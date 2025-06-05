# - adaptMesh (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import XCore.PyTree as XC
import Converter.Internal as Internal
import numpy

# HEXA but no adaptation in the 3rd direction (2D adaptation)
a = G.cartHexa((0,0,0),(0.1,0.1,0.1),(11,11,2))
a = C.convertArray2NGon(a)
C._fillEmptyBCWith(a, 'nref', 'BCFarfield', dim=2)
Internal._adaptNGon32NGon4(a)

normal2D = numpy.array([0.0,0.0,1.0])
ngonelts = Internal.getNodeFromName(a, 'NGonElements')
ER = Internal.getNodeFromName(ngonelts, 'ElementRange')[1]
nfaces = ER[1]
nfaceselts = Internal.getNodeFromName(a, 'NFaceElements')
ER = Internal.getNodeFromName(nfaceselts, 'ElementRange')[1]
ncells = ER[1]-ER[0]+1
gcells = numpy.arange(0, ncells)
gfaces = numpy.arange(1,nfaces+1)

AM = XC.AdaptMesh_Init(a, normal2D, comm=[], gcells=gcells, gfaces=gfaces)
C._initVars(a, '{centers:indicator}=({centers:CoordinateX}>0.5)')
f = Internal.getNodeFromName(a, 'indicator')[1]
REF = f.astype(dtype=Internal.E_NpyInt)
XC.AdaptMesh_AssignRefData(AM, REF)
XC.AdaptMesh_Adapt(AM)
a = XC.AdaptMesh_ExtractMesh(AM, conformize=0)
C.convertPyTree2File(a, "out.cgns")
