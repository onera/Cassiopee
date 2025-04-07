# - rmNodesByName (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

t = C.newPyTree(['Base', 'Base2'])
for i in range(10):
    a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
    a[0] = 'Cart'+str(i)
    t[2][1][2].append(a)

t = Internal.rmNodesByName(t, 'Cart0')
t = Internal.rmNodesByName(t, 'Cart*')
C.convertPyTree2File(t, 'out.cgns')
