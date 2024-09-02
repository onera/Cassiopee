# - applyBCOverlaps (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

# NGON
a = G.cylinder((0,0,0), 1., 1.5, 360., 0., 1., (30,30,10))
s = T.subzone(a,(1,30,1),(-1,-1,-1))
s = C.convertArray2NGon(s)
a = C.convertArray2NGon(a)
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', subzone=s)
C._initVars(a, 'Density', 1.)
t = C.newPyTree(['Base',a])
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t1 = X.applyBCOverlaps(t, depth=1, loc='nodes')
test.testT(t1, 1)

t1 = X.applyBCOverlaps(t, depth=1, loc='centers')
test.testT(t1, 2)

t1 = X.applyBCOverlaps(t, depth=2, loc='nodes')
test.testT(t1, 3)

t1 = X.applyBCOverlaps(t, depth=2, loc='centers')
test.testT(t1, 4)
