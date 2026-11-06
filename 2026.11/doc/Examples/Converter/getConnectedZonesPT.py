# - getConnectedZones (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Connector.PyTree as X

a = G.cart((0,0,0), (1,1,1), (11,11,1))
b = G.cart((10,0,0), (1,1,1), (11,11,1))

t = C.newPyTree(['Base',a,b])
t = X.connectMatch(t, dim=2)

zones = C.getConnectedZones(t[2][1][2][1], topTree=t)
for z in zones: print(z[0])
#>> cart
