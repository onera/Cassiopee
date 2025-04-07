# - applyBCOverlaps (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import KCore.test as test

a = G.cylinder((0,0,0), 1., 1.5, 360., 0., 1., (30,30,10))
# --- CL
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmin')
# --- champ aux noeuds
C._initVars(a, 'Density', 1.)
t = C.newPyTree(['Base',a])
# --- Equation state
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
# --- Apply on PyTree
t1 = X.applyBCOverlaps(t, depth=1,loc='centers')
test.testT(t1, 1)

# in place
C._initVars(t,"centers:cellN=1")
X._applyBCOverlaps(t, depth=1,loc='centers',checkCellN=False)
test.testT(t, 1)


# --- Apply on a zone
a2 = X.applyBCOverlaps(a, depth=1,loc='centers')
t2 = C.newPyTree(['Base',a2])
test.testT(t2, 2)

# depth = 2
t = C.newPyTree(['Base',a])
# --- Equation state
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
# --- CL
t2 = X.applyBCOverlaps(t, depth=2,loc='centers')
test.testT(t2, 3)

# Apply on a zone
a2 = X.applyBCOverlaps(a, depth=2,loc='centers')
t2 = C.newPyTree(['Base',a2])
test.testT(t2,4)

# Apply at nodes
a2 = X.applyBCOverlaps(a, depth=2, loc='nodes')
t2 = C.newPyTree(['Base',a2])
test.testT(t2,5)
