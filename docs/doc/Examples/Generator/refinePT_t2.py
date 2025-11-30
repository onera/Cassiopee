# - refine (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Connector.PyTree as X
import KCore.test as test

# 3D
a = G.cart( (0,0,0), (0.1,0.1,0.1), (11,11,11))
b = G.cart( (0,1,0), (0.1,0.1,0.1), (11,11,11))
c = G.cart( (0,0,1), (0.1,0.1,0.1), (11,11,11))
d = G.cart( (0,1,1), (0.2,0.1,0.1), (6,11,11))

t = C.newPyTree(["BASE", [a,b,c,d]])
t = X.connectMatch(t)
t = X.connectNearMatch(t, ratio=[2,1,1])
C._fillEmptyBCWith(t, "wall", "BCWall")

t2 = G.refine(t, power=2., dir=0)
test.testT(t2, 1)
#
a = G.cart( (0,0,0), (0.05,0.05,0.05), (21,21,21))
b = G.cart( (0,1,0), (0.05,0.05,0.05), (21,21,21))
c = G.cart( (0,0,1), (0.05,0.05,0.05), (21,21,21))
d = G.cart( (0,1,1), (0.1,0.05,0.05), (11,21,21))

t = C.newPyTree(["BASE", [a,b,c,d]])
t = X.connectMatch(t)
t = X.connectNearMatch(t, ratio=[2,1,1])
C._fillEmptyBCWith(t, "wall", "BCWall")
t2 = G.refine(t, power=0.5, dir=0)
test.testT(t2, 2)
