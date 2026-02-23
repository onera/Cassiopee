# - merge (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test

t1 = C.newPyTree(['Base1', 2, 'Base2', 3])
t2 = C.newPyTree(['Other1', 'Other2', 'Base2'])
a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
t1[2][1][2] += [a]
t = Internal.merge([t1, t2])
test.testT(t, 1)
