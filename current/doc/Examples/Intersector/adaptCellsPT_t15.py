# ISO_PHQ4 test
import Geom.PyTree as D
import Converter.PyTree as C
import Intersector.PyTree as XOR
import Generator.PyTree as G
import Transform.PyTree as T
import Post.PyTree as P
import Converter.Internal as I
import time
import KCore.test as test
import numpy

a = D.triangle((0,0,0), (0,1,0), (1, 0, 0))
d = G.cart((1,0.,0.), (1,1,1),(2,1,1))
a = G.addNormalLayers(a, d)
a = C.convertArray2NGon(a)
a = G.close(a)
#C.convertPyTree2File(a,'prism.tp')

cv0 = numpy.empty((1,), dtype=I.E_NpyInt)
cv0[0]=1

cv = []
cv.append(cv0)


a = C.fillEmptyBCWith(a, 'wall', 'BCWall')
a = C.initVars(a, '{centers:Density} = {centers:CoordinateX} + {centers:CoordinateY}')
#C.convertPyTree2File(a, 'out0.cgns')

XOR._setZonesAndJoinsUId(a)

m1 = XOR.adaptCells(a, cv, sensor_type=3, subdiv_type=1) # ISO_HEX
m1 = XOR.closeCells(m1)

I._rmNodesByName(m1, 'zid')
I._rmNodesByName(m1, 'rid')

#C.convertPyTree2File(m1, 'PT_t15_1.cgns')
test.testT(m1,1)

m2 = XOR.adaptCells(a, cv, sensor_type=3, subdiv_type=1, smoothing_type=1) # ISO_HEX
m2 = XOR.closeCells(m2)

I._rmNodesByName(m2, 'zid')
I._rmNodesByName(m2, 'rid')

#C.convertPyTree2File(m2, 'PT_t15_2.cgns')
test.testT(m2,2)
