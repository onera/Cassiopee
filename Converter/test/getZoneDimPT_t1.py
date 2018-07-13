# - getZoneDim (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
dim = Internal.getZoneDim(a)
test.testO(dim, 1)

a = G.cartTetra( (0,0,0), (1,1,1), (10,10,10) )
dim = Internal.getZoneDim(a)
test.testO(dim, 2)
