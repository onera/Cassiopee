# - getVarNames (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# Sur une zone
a = G.cart((0,0,0),(1,1,1),(10,10,10))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
C._addVars(a, ['Density', 'centers:cellN'])
val = C.getVarNames(a)
test.testO(val, 1)

# Sur un arbre
b = C.addVars(a, ['Density', 'centers:cellN']); b[0] = 'cart2'
t = C.newPyTree(['Base',a,b])
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
val = C.getVarNames(t)
test.testO(val, 2)

# Sur une liste de zones
val = C.getVarNames([a,b])
test.testO(val, 3)
