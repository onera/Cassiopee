# - adapts a cells with respect to b points (PyTree) -
import Intersector.PyTree as XOR
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test
import Converter.Internal as I

a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.convertArray2NGon(a); a = G.close(a)
b = G.cartHexa((0.,0.,0.), (0.005,0.005,0.02), (5,5,2))

m = XOR.adaptCells(a,b, sensor_type=4, sensor_metric_policy=0, conformize=0)

test.testT(m,1)
#C.convertPyTree2File(m, 'PT_t17.cgns')
