# - getNodeFromName (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
t = C.newPyTree(['Base', a])

# Retourne le noeud de nom cart
node = Internal.getNodeFromName(t, 'cart')
test.testT(node, 1)
