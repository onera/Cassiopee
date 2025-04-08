# - getCEBBIntersectingDomains 2D (pyTree) -
import Connector.PyTree as X
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

a = G.cart((0.,0.,0.),(0.1,0.1,1),(10,10,1)); a[0] = 'cart1'
b = G.cart((0.5,0.,0.),(0.1,0.1,1),(10,10,1)); b[0] = 'cart2'
c = G.cart((0.75,0.,0.),(0.1,0.1,1),(10,10,1)); c[0] = 'cart3'
# --- CL
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
b = C.addBC2Zone(b, 'wall1', 'BCWall', 'imin')
b = C.addBC2Zone(b, 'overlap1', 'BCOverlap', 'jmax')
c = C.addBC2Zone(c, 'wall1', 'BCWall', 'imin')

t = C.newPyTree(['Cart'])
# --- Equation state
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
t[2][1][2] += [a, b, c]
# --- champ aux centres
C._initVars(t, 'centers:cellN', 1.)
# --- champ aux noeuds
C._initVars(t, 'F', 2.)
bases = Internal.getNodesFromType(t,'CGNSBase_t')
base = bases[0]
doms = X.getCEBBIntersectingDomains(base, bases, sameBase=1)
test.testO(doms)
doms = X.getCEBBIntersectingDomains(base, bases, sameBase=0)
test.testO(doms,2)
