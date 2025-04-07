# - breakElements (array) -
import Converter as C
import Generator as G
import Transform as T

a = G.cartTetra((0,0,0),(1,1,1),(3,3,2))
a = C.convertArray2NGon(a)
a = G.close(a)
b = G.cartNGon((2,0,0),(1,1,1),(3,2,2))
res = T.join(a,b)
res = T.breakElements(res)
C.convertArrays2File(res, 'out.plt')
