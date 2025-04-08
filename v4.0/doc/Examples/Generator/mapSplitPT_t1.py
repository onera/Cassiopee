# - mapSplit (pyTree) -
import Generator.PyTree as G
import Geom.PyTree as D
import KCore.test as test
import Converter.PyTree as C

# polyline
a = D.polyline([(0,0,0),(1,0,0),(1,1,0),(2,3,0),(1.5,3,0),(1,1.5,0),(0,0,0)])
C._addVars(a,'Density'); C._initVars(a,'centers:cellN',1.)
# distribution
Ni = 41
dist = G.cart((0,0,0),(1./(Ni-1),1,1),(Ni,1,1))
dist = G.enforceX(dist, 15.5/(Ni-1), 0.005, 2,5)
dist = G.enforceX(dist, 27.5/(Ni-1), 0.005, 2,5)

C._initVars(a,'F',1.); C._initVars(a,'centers:F',2.)
zones = G.mapSplit(a, dist, 0.25)
t = C.newPyTree(['Base',1]); t[2][1][2] += zones
test.testT(t)
