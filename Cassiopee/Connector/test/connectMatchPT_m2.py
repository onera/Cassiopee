# - connectMatchNGon (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import Generator.PyTree as G
import Connector.Mpi as Xmpi
import Converter.Filter as Filter
import KCore.test as test
import os

LOCAL = test.getLocal()

N = 5

def _createTest(filepath, dim=3, api=3):
    if Cmpi.master:
        zones = []
        for rank in range(2):
            blhc = [float(rank)*(N-1), 0., 0.]
            if dim == 3: Nxyz = [N, N, N]; trhc = (1., 1., 1.)
            else: Nxyz = [N, N, 1]; trhc = (1., 1., 0.)
            a = G.cartNGon(blhc, trhc, Nxyz, api=api)
            def f(x, y, z): return 3.*x + 2.*y + z
            C._initVars(
                a,
                'centers:Density',
                f,
                ['centers:CoordinateX', 'centers:CoordinateY', 'centers:CoordinateZ'],
                isVectorized=True
            )
            C._initVars(a, 'F', 0.)
            Cmpi._setProc(a, rank)
            zones.append(a)
        t = C.newPyTree(["Base", *zones])
        C.convertPyTree2File(t, filepath)
    Cmpi.barrier()

def runTest(filepath):
    h = Filter.Handle(filepath)
    t = h.loadFromProc()
    zones = Internal.getZones(t)
    Xmpi._connectMatchNGon(zones[0])
    if Cmpi.master: os.remove(filepath)
    Cmpi.barrier()
    return t


filepath = os.path.join(LOCAL, "out.cgns")

# 2D NGON
dim = 2
_createTest(filepath, dim=dim, api=1)
t = runTest(filepath)
if Cmpi.master: test.testT(t, 1)

_createTest(filepath, dim=dim, api=3)
t = runTest(filepath)
if Cmpi.master: test.testT(t, 2)

# 3D NGON
dim = 3
_createTest(filepath, dim=dim, api=1)
t = runTest(filepath)
if Cmpi.master: test.testT(t, 3)

_createTest(filepath, dim=dim, api=3)
t = runTest(filepath)
if Cmpi.master: test.testT(t, 4)
