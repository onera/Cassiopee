# - adaptNGon3NGon4 (pyTree) -
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

# NGon (CGNSv3)
a = G.cartNGon((0,0,0), (1,1,1), (3,3,3))

# NGon (CGNSv4)
b = Internal.adaptNGon32NGon4(a)
test.testT(b,1)

# Back to CGNSv3
c = Internal.adaptNGon42NGon3(b)
test.testT(c,2)
