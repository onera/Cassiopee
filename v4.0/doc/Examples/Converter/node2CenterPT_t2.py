# - node2Center (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# CAS 1D
def F(x,y): return 2*x+y

a = G.cart((0,0,0), (1,1,1), (10,1,1))
a = C.addVars(a, ['Density', 'H', 'Helio'])
a = C.initVars(a, 'Density', F, ['CoordinateX','CoordinateY'])
a = C.initVars(a, 'Helio', F, ['CoordinateX','CoordinateY'])
a = C.initVars(a, 'centers:cellN', 1.)
a = C.node2Center(a, ['Density', 'Helio'])
t = C.newPyTree(['Base',1]); t[2][1][2].append(a)
test.testT(t)
