# - exteriorFaces (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

# test faces exterieures au sens large
# faces ayant 2 voisins
a = G.cartTetra((0,0,0), (1,1.,1), (20,2,1))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
b = P.exteriorFaces(a)
test.testT(b, 1)

# test faces interieures au sens strict :
# faces n'ayant que des noeuds exterieurs
a = G.cartTetra((0,0,0), (1,1.,1), (20,2,1))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
b = P.exteriorFaces(a)
test.testT(b,2)

# cas 3D
a = G.cart((0,0,0), (1,1,1), (4,4,6))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
b = P.exteriorFaces(a)
test.testT(b,3)

# cas avec indices
indices = []
a = G.cartTetra((0,0,0), (1,1.,1), (20,2,1))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
b = P.exteriorFaces(a,indices=indices)
test.testO(indices, 4)
