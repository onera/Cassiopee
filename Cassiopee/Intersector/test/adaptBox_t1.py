# - adapt the bounding box of a point cloud (array) -
import Intersector as XOR
import Converter as C
import Generator as G
import KCore.test as test

a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.convertArray2NGon(a); a = G.close(a)

m = XOR.adaptBox(a, box_ratio=10.)
m = XOR.closeCells(m)
test.testA(m,1)
