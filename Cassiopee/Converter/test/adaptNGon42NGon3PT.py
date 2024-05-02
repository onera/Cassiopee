# - adaptNGon42NGon3 (pyTree) -
import Generator.PyTree as G
import Converter.Internal as Internal

# NGon (CGNSv4)
a = G.cartNGon((0,0,0), (1,1,1), (3,3,3), api=3)
Internal.printTree(a)

# NGon (CGNSv3)
b = Internal.adaptNGon42NGon3(a)
Internal.printTree(b)