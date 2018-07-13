# - F (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

def F(x): return x

# Une zone
a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
a = C.initVars(a, '{F}={CoordinateX}+{CoordinateY}')
a = C.initVars(a, '{centers:G}={centers:CoordinateX}+{centers:CoordinateY}')
a = C.addBC2Zone(a, 'wall', 'BCWall', 'imin')
b = F(a)
test.testT(b, 1)

# Une liste de zones
a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
b = G.cartTetra( (0,0,0), (1,1,1), (10,10,10) )
B = F([a,b])
test.testT(B, 2)

# Un arbre
a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
t = C.newPyTree(['Base'])
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t[2][1][2] += [a]
t = F(t)
test.testT(t, 3)
