# - symmetrize (PyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C

a = G.cart((0,0,0), (1,1,1), (10,10,2))
# Symmetrize regarding plane (x,z)
b = T.symmetrize(a, (0.,0.,0.), (1,0,0), (0,0,1)); b[0]='cart2'
C.convertPyTree2File([a,b], "out.cgns")
