# - cartNGon (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C

a = G.cartNGon((0.,0.,0.), (0.1,0.1,0.2), (2,2,2))
C.convertPyTree2File(a, 'out.cgns')
