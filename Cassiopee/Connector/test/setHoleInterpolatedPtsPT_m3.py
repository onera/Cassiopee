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
cellNName = "cellN_shock"

def sphere(x,y,z):
    if x*x+y*y+z*z < 0.48**2: return 0.
    else: return 1.

# Field located at cell centers - NGON
if Cmpi.master:
    a = G.cartNGon((-2.,-1.,-1.), (0.1,0.1,0.1), (21,21,21))
    b = G.cartNGon((0.,-1.,-1.), (0.1,0.1,0.1), (21,21,21))
    t = C.newPyTree(['Cart', a, b])
    C._initVars(
        t,
        f'centers:{cellNName}',
        sphere,
        ['centers:CoordinateX', 'centers:CoordinateY', 'centers:CoordinateZ']
    )
    zones = Internal.getZones(t)
    for i, z in enumerate(zones):
        z[0] = f'zone.{i}'
        Cmpi._setProc(z, i)
    C.convertPyTree2File(t, LOCAL+'/out.cgns')
Cmpi.barrier()

h = Filter.Handle(LOCAL+'/out.cgns')
t = h.loadFromProc()
zones = Internal.getZones(t)
Xmpi._connectMatchNGon(zones[0])

depth = 2
Xmpi._setHoleInterpolatedPoints(t, depth, cellNName)
#Cmpi.convertPyTree2File(t, 'out1.cgns')
if Cmpi.master: test.testT(t, 1)
