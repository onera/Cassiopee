# - exteriorFacesStructured (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

# 1D
a = G.cart((0,0,0), (1,1,1), (10,1,1))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
t = C.newPyTree(['Base',1,a])
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
zones = P.exteriorFacesStructured(t)
test.testT(zones,1)

# 2D
a = G.cart((0,0,0), (1,1,1), (1,6,10))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
t = C.newPyTree(['Base',2,a])
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
zones = P.exteriorFacesStructured(t)
test.testT(zones,2)

# 3D
a = G.cart((0,0,0), (1,1,1), (4,4,6))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
t = C.newPyTree(['Base',3,a])
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
zones = P.exteriorFacesStructured(t)
test.testT(zones,3)
