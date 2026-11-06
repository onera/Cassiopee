# - copySort (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart((0,0,0), (1,1,1), (10,10,10)); a[0] = 'a'
b = G.cart((12,0,0), (1,1,1), (10,10,10)); b[0] = 'b'
t1 = C.newPyTree(['Base',a,b])
t2 = C.newPyTree(['Base',b,a])
t2 = Internal.copySort(t1,t2)
C.convertPyTree2File(t2, 'out.cgns')
