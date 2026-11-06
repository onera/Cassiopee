# - convertArray2Hexa (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# cas non structures
# 2D: tri->quad
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (3,3,1))
C._initVars(a, '{F}={CoordinateX}')
C._initVars(a, '{centers:G}={centers:CoordinateY}')
a = C.convertArray2Hexa(a)
t = C.newPyTree(['Base',2,a])
test.testT(t,1)

# 3D: tetra->hexa
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
C._initVars(a,'{F}={CoordinateX}')
C._initVars(a, '{centers:G}={centers:CoordinateZ}')
a = C.convertArray2Hexa(a)
t = C.newPyTree(['Base',3,a])
test.testT(t,2)

# 3D: penta->hexa
a = G.cartPenta((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
C._initVars(a,'{F}={CoordinateX}')
C._initVars(a, '{centers:G}={centers:CoordinateZ}')
a = C.convertArray2Hexa(a)
t = C.newPyTree(['Base',3,a])
test.testT(t,3)

# Sur un arbre
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
b = G.cartPenta((10.,0.,0.), (0.1,0.1,0.2), (10,10,10))
t = C.newPyTree(['Base',3]); t[2][1][2] += [a,b]
C._initVars(t, '{F}={CoordinateX}' ); C._initVars(t, '{centers:G}={centers:CoordinateZ}')
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = C.convertArray2Hexa(t)
test.testT(t,4)

# Sur une liste de zones
A = C.convertArray2Hexa([a,b])
test.testT(A,5)
