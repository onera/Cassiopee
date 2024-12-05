# - getRegularityMapPT (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C

a = G.cart((0,0,0), (1,1,1), (50,50,1))
a = G.enforceX(a, 25, 0.1, 10, 10)
a = G.getRegularityMap(a)
C.convertPyTree2File(a, 'out.cgns')
