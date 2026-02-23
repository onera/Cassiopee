# - reorder (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C

a = G.cart((0,0,0), (1,1,1), (8,9,20))
a = T.reorder(a, (2,1,-3))
C.convertPyTree2File(a, "out.cgns")
