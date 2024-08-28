# - getPath (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base', a])
b = Internal.getNodeFromName(t, 'Base')
c = Internal.newCGNSTree()

s = False
p1 = Internal.getPath(t, a, pyCGNSLike=s)
p2= Internal.getPath(t, b, pyCGNSLike=s)
p3 = Internal.getPath(b, a, pyCGNSLike=s)
p4 = Internal.getPath(b, b, pyCGNSLike=s)
p5 = Internal.getPath(t, c, pyCGNSLike=s)
test.testO([p1,p2,p3,p4,p5], 1)

s = True
p1 = Internal.getPath(t, a, pyCGNSLike=s)
p2= Internal.getPath(t, b, pyCGNSLike=s)
p3 = Internal.getPath(b, a, pyCGNSLike=s)
p4 = Internal.getPath(b, b, pyCGNSLike=s)
p5 = Internal.getPath(t, c, pyCGNSLike=s)
test.testO([p1,p2,p3,p4,p5], 2)
