# - boolean reorientExternalFaces (array) -
import Generator.PyTree as G
import Converter.PyTree as C
import Intersector.PyTree as XOR
import KCore.test as test


a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
a = C.convertArray2NGon(a)
a = XOR.reorient(a)

t = C.newPyTree(['Base',2]); t[2][1][2].append(a)
test.testT(t,1)
