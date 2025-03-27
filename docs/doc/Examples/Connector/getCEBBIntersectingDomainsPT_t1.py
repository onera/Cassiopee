# - getCEBBIntersectingDomains 3D (pyTree) -
import Connector.PyTree as X
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(10,10,10)); a[0] = 'cart1'
b = G.cart((0.5,0.,0.),(0.1,0.1,0.1),(10,10,10)); b[0] = 'cart2'
c = G.cart((0.75,0.,0.),(0.1,0.1,0.1),(10,10,10)); c[0] = 'cart3'
# --- CL
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmax')
# --- champ aux centres
a = C.initVars(a, 'centers:Density', 1.)
a = C.initVars(a, 'cellN', 2.)
# --- champ aux noeuds
b = C.initVars(b, 'cellN', 1.)

t = C.newPyTree(['Cart'])
# --- Equation state
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t[2][1][2] += [a, b, c]
bases = Internal.getNodesFromType(t,'CGNSBase_t')
base = bases[0]
doms = X.getCEBBIntersectingDomains(base, bases,sameBase=1)
test.testO(doms)
doms = X.getCEBBIntersectingDomains(base, bases,sameBase=0)
test.testO(doms,2)
