# - appendBaseName2ZoneName (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,10)); a[0] = 'b'
b = G.cart((11,0,0), (1,1,1), (10,10,10)); b[0] = 'a'
t = C.newPyTree(['Base',a,b])

t = Internal.appendBaseName2ZoneName(t)
test.testT(t, 1)
