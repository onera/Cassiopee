# - boolean reorientExternalFaces (array) -
import Generator.PyTree as G
import Converter.PyTree as C
import Intersector.PyTree as XOR


a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
a = C.convertArray2NGon(a)
a = XOR.reorient(a)

C.convertPyTree2File(a, 'out.cgns')