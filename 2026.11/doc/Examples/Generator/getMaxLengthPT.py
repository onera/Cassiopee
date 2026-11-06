# - getMaxLength(pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C

a = G.cart((0.,0.,0.), (0.1,0.1,0.1), (11,11,11))
a = G.enforcePlusX(a,1e-6,(5,50))
a = G.getMaxLength(a)
C.convertPyTree2File(a, "out.cgns")
