# - distance2Walls (pyTree) -
import KCore.test as test
import Dist2Walls.PyTree as DW
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D

# Body structured, celln aux centres
sphere = D.sphere((0.5,0.5,0.5), 0.2, 50)
sphere = C.initVars(sphere, 'centers:cellnf', 1.)
bodies = C.newPyTree(['Bodies',2,sphere])

a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(11,11,11))
# --- champ aux noeuds
a = C.initVars(a, 'Density', 0.5)
# --- champ aux centres
a = C.initVars(a, 'Pressure', 0.7)
t = C.newPyTree(['Base', a])
# --- Equation state
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t2 = DW.distance2Walls(t, bodies, type='mininterf', loc='centers')
test.testT(t2,1)

# Body structured, celln aux noeuds
sphere = D.sphere((0.5,0.5,0.5), 0.2, 50)
sphere = C.initVars(sphere,'cellnf', 1.)
bodies = C.newPyTree(['Bodies',2,sphere])

a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(11,11,11))
# --- champ aux noeuds
a = C.initVars(a, 'Density', 0.5)
# --- champ aux centres
a = C.initVars(a, 'Pressure', 0.7)
t = C.newPyTree(['Base',a])
# --- Equation state
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t2 = DW.distance2Walls(t, bodies, type='mininterf', loc='centers')
test.testT(t2,2)
