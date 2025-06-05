# - adaptMesh (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import XCore.PyTree as XC
import Converter.Mpi as Cmpi
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import KCore.test as test
import numpy

LOCAL = test.getLocal()

# 2D QUAD
no_test = 0
a = G.cartHexa((0,0,0),(0.1,0.1,0.1),(11,11,2))
a = C.convertArray2NGon(a)
C._fillEmptyBCWith(a, 'nref', 'BCFarfield', dim=2)
Internal._adaptNGon32NGon4(a)
if Cmpi.rank == 0: C.convertPyTree2File(a, LOCAL+'/tmp.cgns')
Cmpi.barrier()

a, res = XC.loadAndSplitNGon(LOCAL+'/tmp.cgns')
gcells = res[5]
gfaces = res[6]
comm = res[1]

normal2D = numpy.array([0.0,0.0,1.0])
AM = XC.AdaptMesh_Init(a, normal2D, comm=[], gcells=gcells, gfaces=gfaces)
C._initVars(a, 'centers:indicator', 1.)
f = Internal.getNodeFromName(a, 'indicator')[1]
REF = f.astype(dtype=Internal.E_NpyInt)
XC.AdaptMesh_AssignRefData(AM, REF)
XC.AdaptMesh_Adapt(AM)
a = XC.AdaptMesh_ExtractMesh(AM, conformize=0)
test.testT(a, no_test*Cmpi.size+Cmpi.rank)
no_test += 1

C._initVars(a, '{centers:indicator}=({centers:CoordinateX}<0.5)')
f = Internal.getNodeFromName(a, 'indicator')[1]
REF = f.astype(dtype=Internal.E_NpyInt)
XC.AdaptMesh_AssignRefData(AM, REF)
XC.AdaptMesh_Adapt(AM)
a = XC.AdaptMesh_ExtractMesh(AM, conformize=0)
test.testT(a, no_test*Cmpi.size+Cmpi.rank)
no_test += 1

a = XC.AdaptMesh_ExtractMesh(AM, conformize=1)
test.testT(a, no_test*Cmpi.size+Cmpi.rank)
no_test += 1

# 3D HEXA
a = G.cartHexa((0,0,0),(0.1,0.1,0.1),(11,11,11))
a = C.convertArray2NGon(a)
C._fillEmptyBCWith(a, 'nref', 'BCFarfield', dim=3)
Internal._adaptNGon32NGon4(a)
if Cmpi.rank == 0: C.convertPyTree2File(a, LOCAL+'/tmp.cgns')
Cmpi.barrier()

a, res = XC.loadAndSplitNGon(LOCAL+'/tmp.cgns')
gcells = res[5]
gfaces = res[6]
comm = res[1]

normal2D = None

AM = XC.AdaptMesh_Init(a, normal2D, comm=[], gcells=gcells, gfaces=gfaces)
C._initVars(a, 'centers:indicator', 1.)
f = Internal.getNodeFromName(a, 'indicator')[1]
REF = f.astype(dtype=Internal.E_NpyInt)
XC.AdaptMesh_AssignRefData(AM, REF)
XC.AdaptMesh_Adapt(AM)
a = XC.AdaptMesh_ExtractMesh(AM, conformize=0)
test.testT(a,no_test*Cmpi.size+Cmpi.rank)
no_test += 1

C._initVars(a, '{centers:indicator}=({centers:CoordinateX}<0.5)')
f = Internal.getNodeFromName(a, 'indicator')[1]
REF = f.astype(dtype=Internal.E_NpyInt)
XC.AdaptMesh_AssignRefData(AM, REF)
XC.AdaptMesh_Adapt(AM)
a = XC.AdaptMesh_ExtractMesh(AM, conformize=0)
test.testT(a, no_test*Cmpi.size+Cmpi.rank)
no_test += 1

a = XC.AdaptMesh_ExtractMesh(AM, conformize=1)
test.testT(a, no_test*Cmpi.size+Cmpi.rank)
no_test += 1
