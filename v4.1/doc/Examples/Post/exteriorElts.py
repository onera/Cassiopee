# - exteriorElts (array) -
import Converter as C
import Post as P
import Generator as G

a = G.cartTetra((0,0,0), (1,1,1), (10,10,10))
b = P.exteriorElts(a)
C.convertArrays2File(b, 'out.plt')
