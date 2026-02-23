# - getSizeOf (pyTree) -
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,10))
mem = Internal.getSizeOf(a)
test.testO(mem)
