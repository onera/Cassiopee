# - rmVars (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# Sur une zone
a = G.cart((0,0,0),(1,1,1),(10,10,10))
C._addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
C._addVars(a, ['Density', 'centers:cellN', 'rou', 'rov', 'Hx', 'centers:Hy'])
C._rmVars(a, 'Density')
C._rmVars(a, ['Hx', 'centers:Hy'])
C._rmVars(a, 'FlowSolution#Centers')
t = C.newPyTree(['Base',a])
test.testT(t, 1)

# Sur un arbre
a = G.cart((0,0,0),(1,1,1),(10,10,10))
C._addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
C._addVars(a, ['Density', 'centers:cellN'])
b = G.cart((10,0,0),(1,1,1),(10,10,10))
C._addVars(b, ['Density', 'centers:cellN'])
t = C.newPyTree(['Base',a,b])
C._addState(t[2][1], 'Mach', 0.6)
C._rmVars(t, 'Density')
test.testT(t, 2)

# Sur une liste de zones
A = C.rmVars([a,b], 'Density')
t = C.newPyTree(['Base']); t[2][1][2] += A
test.testT(t, 3)
