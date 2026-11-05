# - frontFaces (array) -
import Converter as C
import Generator as G
import Post as P

a = G.cart( (0,0,0), (1,1,1), (11,11,11) )

def F(x, y, z):
    return (x + 2*y + z > 20.).astype(float)

a = C.initVars(a, 'tag', F, ['x', 'y', 'z'], isVectorized=True)
t = C.extractVars(a, ['tag'])
f = P.frontFaces(a, t)
C.convertArrays2File([a,f], 'out.plt')
