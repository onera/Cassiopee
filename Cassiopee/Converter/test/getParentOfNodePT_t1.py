# - getParentOfNode (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
t = C.newPyTree(['Base',a])

# from tree
(p, c) = Internal.getParentOfNode(t, a)
test.testO([p,c], 1)

# from list
(p, c) = Internal.getParentOfNode([a], a)
test.testO(p, 2)

# with limited search
(p, c) = Internal.getParentOfNode2(t, a)
test.testO([p,c], 3)
