# - rmNodesByValue (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
a = C.fillEmptyBCWith(a, 'far', 'BCFarfield')

a = Internal.rmNodesByValue(a, 'BCFarfield')
t = C.newPyTree(['Base',a])
C.convertPyTree2File(t, 'out.cgns')
