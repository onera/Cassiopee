# - adapt the bounding box of a point cloud (array) -
import Intersector.PyTree as XOR
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.convertArray2NGon(a); a = G.close(a)

m = XOR.adaptBox(a, box_ratio=10.)

m = XOR.closeCells(m)
test.testT(m,1)
