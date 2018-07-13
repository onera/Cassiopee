# - conformizeNGon (array) -
import Generator as G
import Converter as C
import Transform as T
a = G.cartNGon((0,0,0),(0.1,0.1,1),(11,11,1))
b = G.cartNGon((1.,0,0),(0.1,0.2,1),(11,6,1))
a = G.cartNGon((0,0,0),(1,1,1),(3,3,1))
b = G.cartNGon((2.,0,0),(2,2,1),(2,2,1))
res = T.join(a,b)
res2 = C.conformizeNGon(res)
C.convertArrays2File(res2, 'out.plt')
