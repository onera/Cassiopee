# - adaptNearMatch (pyTree) -
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.PyTree as C
import Transform.PyTree as T
import Converter.elsAProfile as elsAProfile
import KCore.test as test

a1 = G.cart((0.,0.,0.), (0.1, 0.1, 0.1), (11, 21, 3)); a1[0] = 'cart1'
a2 = G.cart((1., 0.2, 0.), (0.1, 0.1, 0.1), (11, 21, 3)); a2[0] = 'cart2'
a2 = T.oneovern(a2,(1,2,1))
t = C.newPyTree(['Base',a1,a2])
t = X.connectNearMatch(t)
tp = elsAProfile.adaptNearMatch(t)
test.testT(tp,1)
elsAProfile._adaptNearMatch(t)
test.testT(t,1)
