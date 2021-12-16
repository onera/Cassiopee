# - adaptNGon12NGon2 (pyTree) -
import Generator.PyTree as G
import Converter.Internal as Internal

# NGon1
a = G.cartNGon((0,0,0), (1,1,1), (3,3,3))
Internal.printTree(a)

# NGon2
b = Internal.adaptNGon12NGon2(a)
Internal.printTree(b)

# Back to NGon1
c = Internal.adaptNGon22NGon1(b)
Internal.printTree(c)
