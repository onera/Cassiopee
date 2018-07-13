# - rmVars (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# Sur une zone
a = G.cart((0,0,0),(1,1,1),(10,10,10))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
a = C.addVars(a, ['Density', 'centers:cellN', 'rou', 'rov', 'Hx', 'centers:Hy'])
a = C.rmVars(a, 'Density')
a = C.rmVars(a, ['Hx', 'centers:Hy'])
a = C.rmVars(a, 'FlowSolution#Centers')
t = C.newPyTree(['Base',a])
test.testT(t, 1)

# Sur un arbre
a = G.cart((0,0,0),(1,1,1),(10,10,10))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
a = C.addVars(a, ['Density', 'centers:cellN'])
b = G.cart((10,0,0),(1,1,1),(10,10,10))
b = C.addVars(b, ['Density', 'centers:cellN'])
t = C.newPyTree(['Base',a,b])
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = C.rmVars(t, 'Density')
test.testT(t, 2)

# Sur une liste de zones
A = C.rmVars([a,b], 'Density')
t = C.newPyTree(['Base']); t[2][1][2] += A
test.testT(t, 3)
