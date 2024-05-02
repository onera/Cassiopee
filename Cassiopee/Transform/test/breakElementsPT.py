# - breakElements (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T

a = G.cartTetra((0,0,0),(1,1,1),(3,3,2))
a = C.convertArray2NGon(a)
a = G.close(a)
b = G.cartNGon((2,0,0),(1,1,1),(3,3,1))
res = T.join(a,b)
res = T.breakElements(res)
C.convertPyTree2File(res, 'out.cgns')
