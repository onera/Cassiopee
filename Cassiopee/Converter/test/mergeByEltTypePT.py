# - mergeByEltType (pyTree)
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cartHexa((0.,0.4,0.), (0.1,0.1,0.1), (5,5,5))
b = G.cartHexa((0.4,0.4,0.), (0.1,0.1,0.1), (5,5,5))
a = C.mergeConnectivity([a, b], None)
a = C.mergeByEltType(a)
C.convertPyTree2File(a, 'out.cgns')
