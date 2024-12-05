# - interiorFaces (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

# test faces interieures au sens large
# faces ayant 2 voisins
a = G.cartTetra((0,0,0), (1,1.,1), (20,2,1))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
b = P.interiorFaces(a)
test.testT(b, 1)

# test faces interieures au sens strict :
# faces n'ayant que des noeuds interieurs
a = G.cartTetra((0,0,0), (1,1.,1), (20,3,1))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
b = P.interiorFaces(a,1)
test.testT(b, 2)

# test faces interieures au sens strict :
# faces n'ayant que des noeuds interieurs
# ici aucune
a = G.cartTetra((0,0,0), (1,1.,1), (20,2,1))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
b = P.interiorFaces(a,1)
test.testT(b,3)
