# - adapts a cells with respect to b points (PyTree) -
import Intersector.PyTree as XOR
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test
import Converter.Internal as I

a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.convertArray2NGon(a); a = G.close(a)
#C.convertPyTree2File(a, 'a.cgns')
b = G.cartHexa((0.,0.,0.), (0.005,0.005,0.02), (5,5,2))
#C.convertPyTree2File(b, 'b.cgns')

a = C.fillEmptyBCWith(a, 'wall', 'BCWall')
a = C.initVars(a, '{centers:Density} = {centers:CoordinateX} + {centers:CoordinateY}')
#C.convertPyTree2File(a, 'm0.cgns')

XOR._setZonesAndJoinsUId(a)

# static adaptation : ISO_MIN for xsensor2
m = XOR.adaptCells(a,b, sensor_type=4, sensor_metric_policy=0)
m = XOR.closeCells(m)

I._rmNodesByName(m, 'zid')
I._rmNodesByName(m, 'rid')

test.testT(m,1)
#C.convertPyTree2File(m, 'PT_t16_1.cgns')


#static adaptation : ISO_MAX for xsensor2
m = XOR.adaptCells(a,b, sensor_type=4, sensor_metric_policy=2)
m = XOR.closeCells(m)

I._rmNodesByName(m, 'zid')
I._rmNodesByName(m, 'rid')

test.testT(m,2)
#C.convertPyTree2File(m, 'PT_t16_2.cgns')
