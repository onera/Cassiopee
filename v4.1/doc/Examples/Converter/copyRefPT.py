# - copyRef (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a])

# t2 and t are sharing values (also data numpys)
t2 = Internal.copyRef(t)
