# - usurp (array) -
import Converter as C
import Generator as G

ni = 30; nj = 40
m = G.cart((0,0,0), (10,10,1), (ni,nj,1))
m = C.addVar(m, 'ro')

c = C.array('cellN', ni-1, nj-1, 1)
c = C.initVars(c, 'cellN', 1)

C.convertArrays2File([m], "out.plt", "bin_tp")

# Basic test case
import Post as P
r = P.usurp([m], [c]); print(r)
