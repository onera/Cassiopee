# - cellN2OversetHoles (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import KCore.test as test

# 3D non structure cellN en noeuds
a = G.cartHexa((0,0,0),(1,1,1),(10,10,10))
C._initVars(a, 'cellN', 0)
t1 = C.newPyTree(['Base',a])
# --- Equation state
t1[2][1] = C.addState(t1[2][1], 'EquationDimension', 3)
# application on PyTree
t1 = X.cellN2OversetHoles(t1)
# application on zone
a = X.cellN2OversetHoles(a)
t2 = C.newPyTree(['Base',a])
# --- Equation state
t2[2][1] = C.addState(t2[2][1], 'EquationDimension', 3)
test.testT(t1,1)
test.testT(t2,2)
#
# 2D non structure cellN en noeuds
#
a = G.cartHexa((0,0,0),(1,1,1),(10,10,1))
C._initVars(a, 'cellN', 0)
a = X.cellN2OversetHoles(a)
t = C.newPyTree(['Base',2,a])
# --- Equation state
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
test.testT(t,3)
