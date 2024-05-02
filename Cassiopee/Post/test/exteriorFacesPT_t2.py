# - exteriorFaces (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

# test faces exterieures au sens large
# faces ayant 2 voisins
a = G.cartTetra((0,0,0), (1,1.,1), (20,2,1))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
t = C.newPyTree(['Base',1,a])
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = P.exteriorFaces(t)
test.testT(t,1)

# test faces interieures au sens strict :
# faces n'ayant que des noeuds exterieurs
a = G.cartTetra((0,0,0), (1,1.,1), (20,2,1))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
t = C.newPyTree(['Base',1,a])
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = P.exteriorFaces(t)
test.testT(t,2)

# cas 3D
a = G.cart((0,0,0), (1,1,1), (4,4,6))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
t = C.newPyTree(['Base',1,a])
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = P.exteriorFaces(t)
test.testT(t,3)
