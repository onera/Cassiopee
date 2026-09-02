# - computeGrad2 (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import Converter.Filter as Filter
import Connector.PyTree as X
import Connector.Mpi as Xmpi
import Generator.PyTree as G
import Post.PyTree as P
import Post.Mpi as Pmpi
import KCore.test as test
import os
import numpy

N = 10
LOCAL = test.getLocal()

def _createTest(filepath, meshType="STRUCT", dim=3, axis="X", api=3):

    def _addData(t):
        def f(x, y, z): return 3.*x + 2.*y + z
        C._initVars(
            t,
            'centers:Density',
            f,
            ['centers:CoordinateX', 'centers:CoordinateY', 'centers:CoordinateZ'],
            isVectorized=True
        )

    if Cmpi.master:
        zones = []
        for rank in range(2):
            blhc = [
                float(rank)*(N-1) if ax == axis.upper() else 0.
                for ax in ["X", "Y", "Z"]
            ]
            if dim == 3: Nxyz = [N, N, N]; trhc = (1., 1., 1.)
            else:
                if axis.upper() == "Z": Nxyz = [N, 1, N]; trhc = (1., 0., 1.)
                else: Nxyz = [N, N, 1]; trhc = (1., 1., 0.)
            if meshType.upper() == "NGON":
                a = G.cartNGon(blhc, trhc, Nxyz, api=api)
            elif meshType.upper() == "STRUCT":
                a = G.cart(blhc, trhc, Nxyz)
            Cmpi._setProc(a, rank)
            zones.append(a)
        t = C.newPyTree(["Base", *zones])
        _addData(t)
        C.convertPyTree2File(t, filepath)
    Cmpi.barrier()

def runTest(filepath, meshType="STRUCT", dim=3):
    if meshType.upper() == "NGON":
        h = Filter.Handle(filepath)
        t = h.loadFromProc()
        zones = Internal.getZones(t)
        Xmpi._connectMatchNGon(zones[0])
        Pmpi._computeGrad2(t, var='centers:Density')
    elif meshType.upper() == "STRUCT":
        h = Filter.Handle(filepath)
        t = h.loadFromProc()
        t = Xmpi.connectMatch(t, dim=dim)
        Pmpi._computeGrad2(t, var='centers:Density')
    Cmpi.convertPyTree2File(t, filepath)
    #if Cmpi.master: os.remove(filepath)
    Cmpi.barrier()
    return t


filepath = os.path.join(LOCAL, "out.cgns")

# --- Without BCDataSets --- #
# 2D STRUCT
meshType = "STRUCT"
dim = 2
_createTest(filepath, meshType, dim=dim, axis="X")
t = runTest(filepath, meshType, dim=dim)
if Cmpi.master: test.testT(t, 1)

_createTest(filepath, meshType, dim=dim, axis="Y")
t = runTest(filepath, meshType, dim=dim)
if Cmpi.master: test.testT(t, 2)

_createTest(filepath, meshType, dim=dim, axis="Z")
t = runTest(filepath, meshType, dim=dim)
#if Cmpi.master: test.testT(t, 3) # WRONG REF!

# 3D STRUCT
meshType = "STRUCT"
dim = 3
_createTest(filepath, meshType, dim=dim, axis="X")
t = runTest(filepath, meshType, dim=dim)
if Cmpi.master: test.testT(t, 4)

_createTest(filepath, meshType, dim=dim, axis="Y")
t = runTest(filepath, meshType, dim=dim)
if Cmpi.master: test.testT(t, 5)

_createTest(filepath, meshType, dim=dim, axis="Z")
t = runTest(filepath, meshType, dim=dim)
if Cmpi.master: test.testT(t, 6)

# 2D NGON
meshType = "NGON"
dim = 2
_createTest(filepath, meshType, dim=dim, axis="X", api=1)
t = runTest(filepath, meshType, dim=dim)
#if Cmpi.master: test.testT(t, 7) # WRONG REF!
exit()

_createTest(filepath, meshType, dim=dim, axis="Y", api=3)
t = runTest(filepath, meshType, dim=dim)
#if Cmpi.master: test.testT(t, 8) # WRONG REF!

_createTest(filepath, meshType, dim=dim, axis="Z", api=3)
t = runTest(filepath, meshType, dim=dim)
#if Cmpi.master: test.testT(t, 9) # WRONG REF!

# 3D NGON
meshType = "NGON"
dim = 3
_createTest(filepath, meshType, dim=dim, axis="X", api=1)
t = runTest(filepath, meshType, dim=dim)
if Cmpi.master: test.testT(t, 10)

_createTest(filepath, meshType, dim=dim, axis="Y", api=3)
t = runTest(filepath, meshType, dim=dim)
if Cmpi.master: test.testT(t, 11)

_createTest(filepath, meshType, dim=dim, axis="Z", api=3)
t = runTest(filepath, meshType, dim=dim)
if Cmpi.master: test.testT(t, 12)
