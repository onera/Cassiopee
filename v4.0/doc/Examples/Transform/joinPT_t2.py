# - join (pyTree) -
import Transform.PyTree as T
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# Join 2 NGON avec TETRA : 2 zones
a1 = G.cartTetra((0.,0.,0.), (1.,1.,1), (11,11,10))
a1 = C.convertArray2NGon(a1)
a1 = C.initVars(a1, 'F', 2.); a1 = C.initVars(a1, 'centers:G', 1)
a2 = G.cartTetra((10.,0.,0.), (1.,1.,1), (10,10,10))
a2 = C.convertArray2NGon(a2)
a2 = C.initVars(a2, 'F', 3.); a2 = C.initVars(a2, 'centers:G', 3)
a = T.join(a1, a2); t = C.newPyTree(['Base', a])
test.testT(t, 1)
#
# Join sur une liste de zones
#
t = C.newPyTree(['Base']); t[2][1][2].append(a1); t[2][1][2].append(a2)
t[2][1][2].append(T.join(t[2][1][2]))
test.testT(t,2)
#
# Join sur un arbre
#
t = C.newPyTree(['Base']); t[2][1][2].append(a1); t[2][1][2].append(a2)
z = T.join(t)
t = C.newPyTree(['Base', z])
test.testT(t,3)
#
# Join 2 NGON issus d'un TETRA et d'un HEXA
a1 = G.cartHexa((0.,0.,0.), (1.,1.,1), (11,11,10))
a1 = C.convertArray2NGon(a1)
a1 = C.initVars(a1, 'F', 2.); a1 = C.initVars(a1, 'centers:G', 1)
a2 = G.cartTetra((10.,0.,0.), (1.,1.,1), (10,10,10))
a2 = C.convertArray2NGon(a2)
a2 = C.initVars(a2, 'F', 3.); a2 = C.initVars(a2, 'centers:G', 3)
a = T.join(a1, a2); t = C.newPyTree(['Base', a])
test.testT(t,4)
