# - exteriorEltsStructured (array) -
import Converter as C
import Post as P
import Generator as G

a = G.cart((0,0,0), (1,1,1), (10,10,1))
p = P.exteriorEltsStructured(a, 2)
C.convertArrays2File(p, 'out.plt')
