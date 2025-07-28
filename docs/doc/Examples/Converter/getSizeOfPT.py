# - getSizeOf (pyTree) -
import Generator.PyTree as G
import Converter.Internal as Internal
import numpy

a = G.cart((0,0,0), (1,1,1), (10,10,10))
print(Internal.getSizeOf(a), 'octets')
#>> 24242 octets

a = [numpy.zeros((100), dtype=numpy.float64)]
print(Internal.getSizeOf(a), 'octets')
#>> 800 octets
