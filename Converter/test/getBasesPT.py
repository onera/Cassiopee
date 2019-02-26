# - getBases (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart((0,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base', a])

# Return nodes of type 'Zone_t'
bases = Internal.getBases(t); print(bases)
#>> [['Base', array(..), [..], 'CGNSBase_t']]
