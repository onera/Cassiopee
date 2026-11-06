# - interiorFaces (array) -
import Converter as C
import Post as P
import Generator as G

# Get interior faces in broad sense:
# faces with 2 neighbours
a = G.cartTetra((0,0,0), (1,1.,1), (20,2,1))
b = P.interiorFaces(a)
C.convertArrays2File([a,b], 'out1.plt')

# Get interior faces in strict sense:
# faces having only interior nodes
a = G.cartTetra((0,0,0), (1,1.,1), (20,3,1))
b = P.interiorFaces(a,1)
C.convertArrays2File([a,b], 'out2.plt')
