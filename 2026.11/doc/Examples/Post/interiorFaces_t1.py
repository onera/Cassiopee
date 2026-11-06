# - interiorFaces (array) -
import Post as P
import Generator as G
import KCore.test as test

# test faces interieures au sens large
# faces ayant 2 voisins
a = G.cartTetra((0,0,0), (1,1.,1), (20,2,1))
b = P.interiorFaces(a)
test.testA([b], 1)

# test faces interieures au sens strict :
# faces n'ayant que des noeuds interieurs
a = G.cartTetra((0,0,0), (1,1.,1), (20,3,1))
b = P.interiorFaces(a,1)
test.testA([b], 2)

# test faces interieures au sens strict :
#faces n'ayant que des noeuds interieurs
# ici aucune
a = G.cartTetra((0,0,0), (1,1.,1), (20,2,1))
b = P.interiorFaces(a,1)
if b[1].shape[1] != 0: print('FAILED...')
