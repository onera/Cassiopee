# - boolean reorientExternalFaces (array) -
import Generator as G
import Converter as C
import Intersector as XOR


a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
a = C.convertArray2NGon(a)

err = XOR.checkCellsClosure(a)
