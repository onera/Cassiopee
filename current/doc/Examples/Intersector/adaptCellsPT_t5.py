import Geom.PyTree as D
import Converter.PyTree as C
import Intersector.PyTree as XOR
import Generator.PyTree as G
import Transform.PyTree as T
import Post.PyTree as P
import Converter.Internal as I
import time
import KCore.test as test

a = G.cartHexa((0.,0.,0.), (1,1,0.2), (2,2,1))
a = G.quad2Pyra(a)
a = C.convertArray2NGon(a)
a = G.close(a)

b = G.cart((0.2,0.2,0.2), (0.01,0.01,0.01),(5,5,5))
b = C.convertArray2NGon(b)
b = G.close(b)

a = C.fillEmptyBCWith(a, 'wall', 'BCWall')
a = C.initVars(a, '{centers:Density} = {centers:CoordinateX} + {centers:CoordinateY}')

#C.convertPyTree2File(a, 'out5.cgns')

XOR._setZonesAndJoinsUId(a)

m1 = XOR.adaptCells(a,b, sensor_type=0)
m1 = XOR.closeCells(m1)

I._rmNodesByName(m1, 'zid')
I._rmNodesByName(m1, 'rid')

#C.convertPyTree2File(m1, 'PT_t5_1.cgns')
test.testT(m1,1)

m2 = XOR.adaptCells(a,b, sensor_type=0, smoothing_type=1)
m2 = XOR.closeCells(m2)

I._rmNodesByName(m2, 'zid')
I._rmNodesByName(m2, 'rid')

#C.convertPyTree2File(m2, 'PT_t5_2.cgns')
test.testT(m2,2)
