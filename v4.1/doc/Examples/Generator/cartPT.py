# - cart (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,11,12))
C.convertPyTree2File(a, 'out.cgns')
