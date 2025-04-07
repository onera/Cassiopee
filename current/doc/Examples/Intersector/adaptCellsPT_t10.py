# - adapts a cells with respect to b points (array) -
import Intersector.PyTree as XOR
import Converter.PyTree as C
import Converter.Internal as I
import Generator.PyTree as G
import Converter.Internal as I
import Transform.PyTree as T
import Connector.PyTree as X
import numpy
import KCore.test as test


N = 5
t = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (N,2,2))

P=2
t = T.splitNParts(t, P, multigrid=0, dirs=[1])

t = C.convertArray2NGon(t); t = G.close(t)

t = X.connectMatch(t)
t = C.fillEmptyBCWith(t, 'wall', 'BCWall')
t = C.initVars(t, '{centers:Density} = {centers:CoordinateX} + {centers:CoordinateY}')

zs = I.getZones(t)
n = C.getNPts(zs[0])

nv0 = numpy.empty((n,), dtype=I.E_NpyInt)
nv0[:] = 2

nv = []
nv.append(nv0)
for i in range(P-1):
	nv.append(None)
#print(nv)

XOR._setZonesAndJoinsUId(t)

ta = XOR.adaptCells(t,nv, sensor_type=2)
ta = XOR.closeCells(ta)

I._rmNodesByName(ta, 'zid')
I._rmNodesByName(ta, 'rid')

test.testT(ta,1)
#C.convertPyTree2File(ta, 'PT_t10_1.cgns')
