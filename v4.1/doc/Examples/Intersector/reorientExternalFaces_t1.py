# - boolean reorientExternalFaces (array) -
import Generator as G
import Converter as C
import Intersector as XOR
import KCore.test as test


a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
a = C.convertArray2NGon(a)
a = XOR.reorient(a)

test.testA(a,1)