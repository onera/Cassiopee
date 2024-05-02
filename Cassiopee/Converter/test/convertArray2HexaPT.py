# - convertArray2Hexa (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

# 2D: quad
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
b = C.convertArray2Hexa(a)

# 3D: hexa
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
b = C.convertArray2Hexa(a)

C.convertPyTree2File([a,b], 'out.cgns')
