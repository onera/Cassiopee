# - frontFaces (array) -
import Converter as C
import Generator as G
import Post as P

a = G.cart( (0,0,0), (1,1,1), (11,11,11) )
def F(x, y, z):
    if (x + 2*y + z   > 20.): return 1
    else: return 0
a = C.initVars(a, 'tag', F, ['x', 'y', 'z'])
t = C.extractVars(a, ['tag'])
f = P.frontFaces(a, t)
C.convertArrays2File([a,f], 'out.plt')
