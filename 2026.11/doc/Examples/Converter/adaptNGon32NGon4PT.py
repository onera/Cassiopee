# - adaptNGon32NGon4 (pyTree) -
import Generator.PyTree as G
import Converter.Internal as Internal

# NGon (CGNSv3)
a = G.cartNGon((0,0,0), (1,1,1), (3,3,3))
Internal.printTree(a)

# NGon (CGNSv4)
b = Internal.adaptNGon32NGon4(a)
Internal.printTree(b)