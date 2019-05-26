# - exteriorFaces (array) -
import Converter as C
import Post as P
import Generator as G

a = G.cartTetra((0,0,0), (1,1,1), (20,20,20))
indices = []
b = P.exteriorFaces(a, indices=indices)
print(indices)
C.convertArrays2File(b, 'out.plt')
