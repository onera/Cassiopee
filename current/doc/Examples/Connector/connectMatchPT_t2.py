# - connectMatch 3D (pyTree)-
import Generator.PyTree as G
import Converter.PyTree as C
import Connector.PyTree as X
import KCore.test as test

a = G.cart((0.,0.,0.), (0.1, 0.1, 0.1), (11, 21, 2))

# partiellement coincident
# --- CL
a2 = G.cart((1., 0.5, 0.), (0.1, 0.1, 0.1), (11, 21, 2))
a2 = C.addBC2Zone(a2, 'wall1', 'BCWall', 'jmax')
a2 = C.addBC2Zone(a2, 'overlap1', 'BCOverlap', 'jmin')
# --- champ aux centres
C._initVars(a2, 'centers:Density', 1.)
# --- champ aux noeuds
C._initVars(a2, 'cellN', 2.)
t = C.newPyTree(['Base',a,a2])
# --- Equation state
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t = X.connectMatch(t)
test.testT(t,1)

# coincident
a3 = G.cart((1.00001, 0., 0.), (0.1, 0.1, 0.1), (11, 21, 2))
t[2][1][2].append(a3)
t = X.connectMatch(t)
test.testT(t,2)
