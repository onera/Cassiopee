# - getElementRange (pyTree) -
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

a = G.cartTetra( (0,0,0), (1,1,1), (10,10,10) )
r1 = Internal.getElementRange(a, number=0)
r2 = Internal.getElementRange(a, type='TETRA')
r3 = Internal.getElementRange(a, name='GridElements')
test.testO(r1+r2+r3, 1)

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
r1 = Internal.getElementRange(a, number=0)
test.testO(r1, 2)

a = G.cartNGon( (0,0,0), (1,1,1), (10,10,10) )
r1 = Internal.getElementRange(a, number=0)
r2 = Internal.getElementRange(a, type='NGON')
r3 = Internal.getElementRange(a, type='NFACE')
test.testO(r1+r2+r3, 3)
