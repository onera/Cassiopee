# - splitManifold (array) -
# Conforming 1 or 2 TRI/BAR together (same type for both operands)
import Converter as C
import Generator as G
import Intersector as XOR
import Geom as D
from Geom.Parametrics import base
import Transform as T
import KCore.test as test

s1 = D.sphere( (0,0,0), 1, N=20 )

s2 = D.surface(base['plane'], N=30)
s2 = T.translate(s2, (0.2,0.2,0.2))

s1 = C.convertArray2Tetra(s1); s1 = G.close(s1)
s2 = C.convertArray2Tetra(s2); s2 = G.close(s2)

x = XOR.conformUnstr(s1, s2, 0., 2)
x = T.splitManifold(x)
test.testA(x,1)


a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,1,50))
c1 = T.subzone(a,(1,1,1),(50,1,1))
c2 = T.subzone(a,(1,1,50),(50,1,50))
c3 = T.subzone(a,(1,1,1),(1,1,50))
c = [c1,c2,c3]; c = C.convertArray2Hexa(c)
c = T.join(c)
x = T.splitManifold(c)
test.testA(x,2)
