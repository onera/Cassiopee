# - createGlobalIndex (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T

a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._createGlobalIndex(a)
C.convertPyTree2File(a, 'out.cgns')

