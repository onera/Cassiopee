# - exteriorFaces (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

# tetra
a = G.cartTetra((0,0,0), (1,1.,1), (20,2,1))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
b = P.exteriorFaces(a)
test.testT(b, 1)

# struct 3D
a = G.cart((0,0,0), (1,1,1), (4,4,6))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
b = P.exteriorFaces(a)
test.testT(b, 3)

# tetra avec indices
indices = []
a = G.cartTetra((0,0,0), (1,1.,1), (20,2,1))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
b = P.exteriorFaces(a, indices=indices)
test.testO(indices, 4)
