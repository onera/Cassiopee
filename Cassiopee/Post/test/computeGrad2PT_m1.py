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
            if dim == 3: Nxyz = [N, N, N]
            else:
                if axis.upper() == "Z": Nxyz = [N, 1, N]
                else: Nxyz = [N, N, 1]
            if meshType.upper() == "NGON":
                a = G.cartNGon(blhc, (1., 1., 1.), Nxyz, api=api)
            elif meshType.upper() == "STRUCT":
                a = G.cart(blhc, (1., 1., 1.), Nxyz)
            Cmpi._setProc(a, rank)
            zones.append(a)
        t = C.newPyTree(["Base", *zones])
        _addData(t)
        C.convertPyTree2File(t, filepath)
    Cmpi.barrier()

def runTest(filepath, meshType="STRUCT"):
    if meshType.upper() == "NGON":
        h = Filter.Handle(filepath)
        t = h.loadFromProc()
        zones = Internal.getZones(t)
        Xmpi._connectMatchNGon(zones[0])
        Pmpi._computeGrad2(t, var='centers:Density')
    elif meshType.upper() == "STRUCT":
        h = Filter.Handle(filepath)
        t = h.loadFromProc()
        t = Xmpi.connectMatch(t)
        Pmpi._computeGrad2(t, var='centers:Density')
    #Cmpi.convertPyTree2File(t, filepath)
    if Cmpi.master: os.remove(filepath)
    Cmpi.barrier()
    return t


filepath = os.path.join(LOCAL, "out.cgns")

# --- Without BCDataSets --- #
# 2D STRUCT
meshType = "STRUCT"
_createTest(filepath, meshType, dim=2, axis="X")
t = runTest(filepath, meshType)
#if Cmpi.master: test.testT(t, 1) # WRONG REF!

_createTest(filepath, meshType, dim=2, axis="Y")
t = runTest(filepath, meshType)
#if Cmpi.master: test.testT(t, 2) # WRONG REF!

_createTest(filepath, meshType, dim=2, axis="Z")
t = runTest(filepath, meshType)
#if Cmpi.master: test.testT(t, 3) # WRONG REF!

# 3D STRUCT
_createTest(filepath, meshType, dim=3, axis="X")
t = runTest(filepath, meshType)
#if Cmpi.master: test.testT(t, 4) # WRONG REF!

_createTest(filepath, meshType, dim=3, axis="Y")
t = runTest(filepath, meshType)
#if Cmpi.master: test.testT(t, 5) # WRONG REF!

_createTest(filepath, meshType, dim=3, axis="Z")
t = runTest(filepath, meshType)
#if Cmpi.master: test.testT(t, 6) # WRONG REF!

# 2D NGON
meshType = "NGON"
_createTest(filepath, meshType, dim=2, axis="X", api=1)
t = runTest(filepath, meshType)
#if Cmpi.master: test.testT(t, 7) # WRONG REF!

_createTest(filepath, meshType, dim=2, axis="Y", api=3)
t = runTest(filepath, meshType)
#if Cmpi.master: test.testT(t, 8) # WRONG REF!

_createTest(filepath, meshType, dim=2, axis="Z", api=3)
t = runTest(filepath, meshType)
#if Cmpi.master: test.testT(t, 9) # WRONG REF!

# 3D NGON
_createTest(filepath, meshType, dim=3, axis="X", api=1)
t = runTest(filepath, meshType)
if Cmpi.master: test.testT(t, 10)

_createTest(filepath, meshType, dim=3, axis="Y", api=3)
t = runTest(filepath, meshType)
if Cmpi.master: test.testT(t, 11)

_createTest(filepath, meshType, dim=3, axis="Z", api=3)
t = runTest(filepath, meshType)
if Cmpi.master: test.testT(t, 12)

