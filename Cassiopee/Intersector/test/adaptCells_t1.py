# - adapts a cells with respect to b points (array) -
import Intersector as XOR
import Converter as C
import Generator as G
import KCore.test as test

a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.convertArray2NGon(a); a = G.close(a)
#C.convertArrays2File([a], 'a.plt')
b = G.cartHexa((0.,0.,0.), (0.005,0.005,0.005), (5,5,5))
#C.convertArrays2File([b], 'b.plt')

m = XOR.adaptCells(a,b, sensor_type=0)
m = XOR.closeCells(m[0])
test.testA(m,1)
#C.convertArrays2File([m], 't1_1.plt')

m = XOR.adaptCells(a,b, sensor_type=1)
m = XOR.closeCells(m[0])
test.testA(m,2)
#C.convertArrays2File([m], 't1_2.plt')
