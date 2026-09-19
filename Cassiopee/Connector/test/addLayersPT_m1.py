# - setHoleInterpolatedPts (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import Generator.PyTree as G
import Connector.Mpi as Xmpi
import Converter.Filter as Filter
import KCore.test as test

LOCAL = test.getLocal()
rank = Cmpi.rank
cellNName = "cellN"

def sphere(x,y,z):
    if x*x+y*y+z*z < 0.48**2: return 0.
    else: return 1.

def sphere2(x,y,z):
    if x*x+y*y+z*z < 0.48**2: return 1.
    else: return 0.

def createTest(func):
    # Field located at cell centers - NGON
    if Cmpi.master:
        a = G.cartNGon((-2.5,-1.,-1.), (0.1,0.1,0.1), (21,21,21))
        b = G.cartNGon((-0.5,-1.,-1.), (0.1,0.1,0.1), (21,21,21))
        t = C.newPyTree(['Cart', a, b])
        C._initVars(
            t,
            f'centers:{cellNName}',
            func,
            ['centers:CoordinateX', 'centers:CoordinateY', 'centers:CoordinateZ']
        )
        zones = Internal.getZones(t)
        for i, z in enumerate(zones):
            z[0] = f'zone.{i}'
            Cmpi._setProc(z, i)
        C.convertPyTree2File(t, LOCAL+'/out.cgns')
    Cmpi.barrier()

# Sphere filled with zeros
createTest(sphere)
h = Filter.Handle(LOCAL+'/out.cgns')
t = h.loadFromProc()
zones = Internal.getZones(t)
Xmpi._connectMatchNGon(zones[0])

nod = 1
for d in [-2, -1, 0, 1, 2, 5]:
    t2 = Internal.copyTree(t)
    Xmpi._addLayers(t2, depth=d, cellNName=cellNName)
    if Cmpi.master: test.testT(t2, nod)
    nod += 1

# Sphere filled with ones
createTest(sphere2)
h = Filter.Handle(LOCAL+'/out.cgns')
t = h.loadFromProc()
zones = Internal.getZones(t)
Xmpi._connectMatchNGon(zones[0])

for d in [-2, -1, 0, 1, 2, 5]:
    t2 = Internal.copyTree(t)
    Xmpi._addLayers(t2, depth=-d, cellNName=cellNName)
    if Cmpi.master: test.testT(t2, nod)
    nod += 1