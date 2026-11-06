# - append (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

t = C.newPyTree(['Base', 'Base2'])

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
t = Internal.append(t, a, '/Base')
test.testT(t, 1)
