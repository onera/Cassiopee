# - adapts a cells with respect to b points (array) -
import Intersector as XOR
import Converter as C
import Generator as G
import KCore.test as test

a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.convertArray2Tetra(a, split='withBarycenters')
a = C.convertArray2NGon(a); a = G.close(a)
#C.convertArrays2File([a], 'a.plt')
b = G.cartHexa((0.,0.,0.), (0.005,0.005,0.005), (5,5,5))
#C.convertArrays2File([b], 'b.plt')

m = XOR.adaptCells(a,b, sensor_type=0)
m = XOR.closeOctalCells(m)
test.testA(m,1)

m = XOR.adaptCells(a,b, sensor_type=1)
m = XOR.closeOctalCells(m)
test.testA(m,2)

