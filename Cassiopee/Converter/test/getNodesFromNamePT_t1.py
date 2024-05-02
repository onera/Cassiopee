# - getNodesFromName (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
t = C.newPyTree(['Base']); t[2][1][2].append(a)

nodes = Internal.getNodesFromName(t, 'cart')
test.testO(nodes, 1)

nodes = Internal.getNodesFromName(t, 'Coordinate*')
test.testO(nodes, 2)
