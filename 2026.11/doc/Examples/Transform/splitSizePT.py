# - splitSize (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C

a = G.cart((0,0,0),(1,1,1),(50,20,10))
t = C.newPyTree(['Base',a])
t = T.splitSize(t, 300, type=0)
C.convertPyTree2File(t, 'out.cgns')
