# - connectNearMatch 2D (pyTree)-
import Generator.PyTree as G
import Converter.PyTree as C
import Connector.PyTree as X
import Transform.PyTree as T
import KCore.test as test

a = G.cart((0.,0.,0.), (0.1, 0.1, 0.1), (11, 21, 1))
# --- champ au centre
a = C.initVars(a, 'centers:Density', 1.)
# --- CL
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
# partiellement coincident
a2 = G.cart((1., 0.4, 0.), (0.1, 0.1, 0.1), (11, 21, 1))
a2 = T.oneovern(a2,(2,2,1))
# --- champ au noeud
a2 = C.initVars(a2, 'Density',1.)
a2 = C.addBC2Zone(a2, 'overlap1', 'BCOverlap', 'imax')
t = C.newPyTree(['Base',2]); t[2][1][2] += [a, a2]
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
t = X.connectNearMatch(t, dim=2)
test.testT(t,1)

# coincident 1s2
a3 = G.cart((1.00001, -1.6, 0.), (0.1, 0.1, 0.1), (11, 21, 1))
a3 = T.oneovern(a3,(2,2,1)); a3[0] = 'cart3'
t[2][1][2] += [a3]
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
t = C.initVars(t, 'centers:Density', 1.)
t = C.initVars(t, 'Density', 1.)
t = X.connectNearMatch(t, dim=2, tol=2.e-5)
test.testT(t,2)
