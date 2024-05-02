# - refine (array) -
import Post as P
import Converter as C
import Generator as G

# Refine using butterfly
a = G.cartTetra((0,0,0), (2,1,1), (3,3,3))
a = P.exteriorFaces(a)
for i in range(5):
    a = P.refine(a, w=1./64.)
C.convertArrays2File(a, 'out.plt')
