# - getNodeFromPath (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
t = C.newPyTree(['Base']); t[2][1][2].append(a)

# Access au noeud root (les 3 sont equivalents)
root = Internal.getNodeFromPath(t, '')
root = Internal.getNodeFromPath(t, '/')
root = Internal.getNodeFromPath(t, 'CGNSTree')
test.testT(root, 3)

# Acces par le noeud root (les 3 sont equivalents)
coords = Internal.getNodeFromPath(t, 'CGNSTree/Base/cart/GridCoordinates')
coords = Internal.getNodeFromPath(t, '/Base/cart/GridCoordinates')
coords = Internal.getNodeFromPath(t, 'Base/cart/GridCoordinates')
test.testT(coords, 1)

# Acces relatif (les 2 sont equivalents)
base = Internal.getNodeFromPath(t, 'Base')
coords = Internal.getNodeFromPath(base, 'cart')
coords = Internal.getNodeFromPath(base, '/cart')
test.testT(coords, 2)
