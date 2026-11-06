# - rmNodeByPath (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

t = C.newPyTree(['Base', 'Base2'])
for i in range(10):
    a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
    a[0] = 'Cart'+str(i)
    t[2][1][2].append(a)

t = Internal.rmNodeByPath(t, 'Base/Cart2')
t = Internal.rmNodeByPath(t, 'Base/Cart4')
C.convertPyTree2File(t, 'out.cgns')
