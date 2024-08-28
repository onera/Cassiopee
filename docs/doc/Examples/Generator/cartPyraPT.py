# - cartPyra (pyTree)-
import Generator.PyTree as G
import Converter.PyTree as C

a = G.cartPyra((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
C.convertPyTree2File(a, 'out.cgns')
