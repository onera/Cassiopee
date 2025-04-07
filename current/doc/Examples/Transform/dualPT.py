# - dual (pyTree)
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T

ni = 5; nj = 5
a = G.cart((0,0,0),(1,1,1),(ni,nj,1))
a = C.convertArray2NGon(a); a = G.close(a)
res = T.dual(a)
C.convertPyTree2File(res, 'out.cgns')
