# - getConnectedZones (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Connector.PyTree as X
import KCore.test as test

a = G.cart( (0,0,0), (1,1,1), (11,11,1) )
b = G.cart( (10,0,0), (1,1,1), (11,11,1) )

t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
t = X.connectMatch(t, dim=2)

zones = C.getConnectedZones(t[2][1][2][1], topTree=t)
test.testT(zones, 1)
