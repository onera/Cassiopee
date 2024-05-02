# - copyTree (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a])

# t2 is a full copy of t (values/numpys are also copied)
t2 = Internal.copyTree(t)
test.testT(t2, 1)
