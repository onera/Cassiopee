# - getMaxValue (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import Converter.Mpi as Cmpi
import KCore.test as test

def F(x,y,z): return 3*x + 2*y + z

N = 5
a = G.cartHexa((0.,0.,0.), (1.,1.,1.), (N,N,N))
b = G.cartHexa((N-1,0.,0.), (1.,1.,1.), (N,N,N))
t = C.newPyTree(['Base', [a, b]])

zones = Internal.getZones(t)
for i, z in enumerate(zones):
    Cmpi._setProc(z, i)

C._initVars(t, 'Density', F, ['CoordinateX', 'CoordinateY', 'CoordinateZ'], isVectorized=True)

if Cmpi.size > 1:
    Cmpi._convert2PartialTree(t, Cmpi.rank)

maxv = Cmpi.getMaxValue(t, 'Density')
if Cmpi.master:
    test.testO(maxv, 1)
