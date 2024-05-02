# - getSizeOf (pyTree) -
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart((0,0,0), (1,1,1), (10,10,10))
print(Internal.getSizeOf(a), 'octets')
#>> 24242 octets