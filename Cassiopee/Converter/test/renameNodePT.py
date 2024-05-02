# - renameNode (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import Connector.PyTree as X

a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = G.cart((9,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a,b])
t = X.connectMatch(t)

# Change the zone named 'cart' and its reference in BCMatch (GridConnectivity)
t = Internal.renameNode(t, 'cart', 'myCart')
C.convertPyTree2File(t, 'out.cgns')
