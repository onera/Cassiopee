# - rmNodes (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# Sur un zone
a = G.cart((0,0,0),(1,1,1),(10,10,10))
a = C.addVars(a, ['Density', 'centers:cellN', 'rou', 'rov', 'Hx', 'Hy'])
a = C.rmNodes(a, 'FlowSolution#Centers')
t = C.newPyTree(['Base']); t[2][1][2].append(a)
test.testT(t, 1)

# Sur un arbre
a = G.cart((0,0,0),(1,1,1),(10,10,10))
a = C.addVars(a, ['Density', 'centers:cellN', 'rou', 'rov', 'Hx', 'Hy'])
b = G.cart((10,0,0),(1,1,1),(10,10,10))
b = C.addVars(b, ['Density', 'centers:cellN', 'rou', 'rov', 'Hx', 'Hy'])
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = C.rmNodes(t, 'FlowSolution#Centers')
test.testT(t, 2)
# Sur une liste de zones
A = C.rmNodes([a,b], 'FlowSolution#Centers')
t = C.newPyTree(['Base']); t[2][1][2] += A
test.testT(t, 3)
