# - connectNearMatch (pyTree) -
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.PyTree as C
import Transform.PyTree as T

a1 = G.cart((0.,0.,0.), (0.1, 0.1, 0.1), (11, 21, 3)); a1[0] = 'cart1'
a2 = G.cart((1., 0.2, 0.), (0.1, 0.1, 0.1), (11, 21, 3)); a2[0] = 'cart2'
a2 = T.oneovern(a2,(1,2,1))
t = C.newPyTree(['Base',a1,a2])
t = X.connectNearMatch(t)
C.convertPyTree2File(t, 'out.cgns')
