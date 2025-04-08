# - connectMatch 2D (pyTree)-
import Generator.PyTree as G
import Converter.PyTree as C
import Connector.PyTree as X
import KCore.test as test

a = G.cart((0.,0.,0.), (0.1, 0.1, 0.1), (11, 21, 1))
# --- champ au centre
C._initVars(a, 'centers:Density', 1.)
# --- CL
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
# partiellement coincident
a2 = G.cart((1., 0.5, 0.), (0.1, 0.1, 0.1), (11, 21, 1))
# --- champ au noeud
C._initVars(a2, 'Density',1.)
a2 = C.addBC2Zone(a2, 'wall1', 'BCOverlap', 'imax')
t = C.newPyTree(['Base',2]); t[2][1][2] += [a, a2]
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
t = X.connectMatch(t, dim=2)
test.testT(t,1)

# coincident
a3 = G.cart((1.00001, 0., 0.), (0.1, 0.1, 0.1), (11, 21, 1))
t[2][1][2].append(a3)
t = X.connectMatch(t, dim=2)
test.testT(t,2)
