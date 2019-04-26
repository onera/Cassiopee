# - adapts a cells with respect to b points (PyTree) -
import Intersector.PyTree as XOR
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.convertArray2Tetra(a, split='withBarycenters')
a = C.convertArray2NGon(a); a = G.close(a)
#C.convertArrays2File([a], 'a.plt')
b = G.cartHexa((0.,0.,0.), (0.005,0.005,0.005), (5,5,5))
#C.convertArrays2File([b], 'b.plt')

m = XOR.adaptCells(a,b, sensor_type=0)
m = XOR.closeOctalCells(m)
test.testT(m,1)

m = XOR.adaptCells(a,b, sensor_type=1)
m = XOR.closeOctalCells(m)
C.convertPyTree2File(m, 'out.cgns')
test.testT(m,2)

