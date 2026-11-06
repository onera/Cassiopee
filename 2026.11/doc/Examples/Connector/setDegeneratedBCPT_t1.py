# - setDegeneratedBC (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import KCore.test as test

# 3D
a = G.cylinder((0,0,0), 0., 1., 360., 0., 1, (21,21,21))
C._initVars(a,'F',1.); C._initVars(a,'centers:G', 2.)
a = C.addBC2Zone(a,'nref','BCFarfield','jmax')
t = C.newPyTree(['Base', a])
t = X.setDegeneratedBC(t)
test.testT(t,1)

# 2D DegeneratePoint
a = G.cylinder((0,0,0), 0., 1., 360., 0., 1, (21,21,1))
C._initVars(a,'F', 1.); C._initVars(a,'centers:G', 2.)
a = C.addBC2Zone(a,'nref','BCFarfield','jmax')
t = C.newPyTree(['Base', 2, a])
t = X.setDegeneratedBC(t, dim=2)
test.testT(t,2)
