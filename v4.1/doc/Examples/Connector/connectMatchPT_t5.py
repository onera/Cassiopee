# - connectMatch 3D (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import KCore.test as test

# liste de zones
a = G.cart((0.,0.,0.), (0.1, 0.1, 0.1), (11, 21, 2)); a[0] = 'cart1'

# partiellement coincident
# --- CL
a2 = G.cart((1., 0.5, 0.), (0.1, 0.1, 0.1), (11, 21, 2)); a2[0] = 'cart2'
a2 = C.addBC2Zone(a2, 'wall1', 'BCWall', 'jmax')
a2 = C.addBC2Zone(a2, 'overlap1', 'BCOverlap', 'jmin')
# --- champ aux centres
a2 = C.initVars(a2, 'centers:Density', 1.)
# --- champ aux noeuds
a2 = C.initVars(a2, 'cellN', 2.)
t = C.newPyTree(['Base',a,a2])
# --- Equation state
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t[2][1][2] = X.connectMatch(t[2][1][2])
test.testT(t,1)

# coincident
a3 = G.cart((1.00001, 0., 0.), (0.1, 0.1, 0.1), (11, 21, 2)); a3[0] = 'cart3'
t = C.newPyTree(['Base']); t[2][1][2] += [a, a3]
C._addState(t[2][1], 'EquationDimension', 3)
t[2][1][2] = X.connectMatch(t[2][1][2], tol=2.e-5)
test.testT(t,2)

# CAS NGON
a1 = G.cartNGon((0.,0.,0.), (0.1, 0.1, 0.1), (11, 21, 11)); a1[0] = 'cart1'
a2 = G.cartNGon((1.,0.,0.), (0.1, 0.1, 0.1), (11, 21, 11)); a2[0] = 'cart2'
t = C.newPyTree(['Base']); t[2][1][2] += [a1,a2]
C._addState(t[2][1], 'EquationDimension', 3)
C._initVars(t,'centers:Density', 1.); C._initVars(t, 'F', 0.)
t = X.connectMatch(t)
test.testT(t,3)

# Cas hybride
a1 = G.cartNGon((0.,0.,0.), (0.1, 0.1, 0.1), (11, 21, 11)); a1[0] = 'cartNGon'
a2 = G.cart((1.,0.,0.), (0.1, 0.1, 0.1), (11, 21, 11)); a2[0] = 'cart2'
t = C.newPyTree(['Base',a1,a2])
C._addState(t[2][1], 'EquationDimension', 3)
C._initVars(t,'centers:Density', 1.); C._initVars(t, 'F', 0.)
t = X.connectMatch(t, type='hybrid')
test.testT(t,4)
