# - convertPyTree2FFD (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = C.convertArray2NGon(a)
C.convertPyTree2FFD(a)
