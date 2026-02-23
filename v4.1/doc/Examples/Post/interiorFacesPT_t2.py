# - interiorFaces (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

# test faces interieures au sens large
# faces ayant 2 voisins
a = G.cartTetra((0,0,0), (1,1.,1), (20,2,1))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
t = C.newPyTree(['Base',2,a])
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = P.interiorFaces(t)
test.testT(t,1)

# test faces interieures au sens strict :
# faces n'ayant que des noeuds interieurs
a = G.cartTetra((0,0,0), (1,1.,1), (20,3,1))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
t = C.newPyTree(['Base',2,a])
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = P.interiorFaces(t,1)
test.testT(t,2)

# test faces interieures au sens strict :
# faces n'ayant que des noeuds interieurs
# ici aucune
a = G.cartTetra((0,0,0), (1,1.,1), (20,2,1))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
t = C.newPyTree(['Base',2,a])
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = P.interiorFaces(t,1)
test.testT(t,3)
