# - getZoneType (pyTree) -
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
type = Internal.getZoneType(a)
test.testO(type, 1)

a = G.cartTetra( (0,0,0), (1,1,1), (10,10,10) )
type = Internal.getZoneType(a)
test.testO(type, 2)
