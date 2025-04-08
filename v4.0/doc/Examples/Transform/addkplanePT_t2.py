# - addkplane (pyTree) -
# traitement par arbre
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import KCore.test as test

# 3D structure + champs en noeuds + [champs en centres] + CL
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(11,10,2))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin', [1,2,3])
a = C.addBC2Zone(a, 'nearmatch1', 'BCNearMatch', 'jmax', a, 'jmin', [1,2,3])
a = C.initVars(a, 'F',2.); a = C.initVars(a, 'centers:G',1.)
t = C.newPyTree(['Base',a])
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
t = T.addkplane(t)
test.testT(t, 1)

# 2D structure + champs en noeuds + [champs en centres] + CL
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(11,10,1))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin', [1,2])
a = C.addBC2Zone(a, 'nearmatch1', 'BCNearMatch', 'jmax', a, 'jmin', [1,2,3])
a = C.initVars(a, 'F',2.); a = C.initVars(a, 'centers:G',1.)
t = C.newPyTree(['Base',2,a])
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
t = T.addkplane(t)
test.testT(t, 2)

# BAR
a = G.cart((0.,0.,0.),(0.1,1,1),(11,1,1))
a = C.convertArray2Tetra(a)
a = T.addkplane(a)
a = C.initVars(a, 'F',2.); a = C.initVars(a, 'centers:G',1.)
t = C.newPyTree(['Base',a])
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t = T.addkplane(t)
test.testT(t, 3)

# QUAD
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(11,10,1))
a = C.convertArray2Hexa(a)
a = C.initVars(a, 'F',2.); a = C.initVars(a, 'centers:G',1.)
t = C.newPyTree(['Base',a])
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t = T.addkplane(t)
test.testT(t, 4)

# TRI
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(11,10,1))
a = C.convertArray2Tetra(a)
a = C.initVars(a, 'F',2.); a = C.initVars(a, 'centers:G',1.)
t = C.newPyTree(['Base',a])
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t = T.addkplane(t)
test.testT(t, 5)
