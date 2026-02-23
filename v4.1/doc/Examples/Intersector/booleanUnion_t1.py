# - boolean union (array) -
import Intersector as XOR
import Generator as G
import Converter as C
import Geom as D
import KCore.test as test

s1 = D.sphere((0,0,0), 1, N=20)
s2 = D.sphere((0.,1.,0.), 1, N=30)
s1 = C.convertArray2Tetra(s1); s1 = G.close(s1)
s2 = C.convertArray2Tetra(s2); s2 = G.close(s2)
x = XOR.booleanUnion(s1, s2, tol=0.)
test.testA([x],1)

# Test avec resserements
s1 = D.sphere((0,0,0), 1, N=20)
s2 = D.sphere((0.,1.,0.), 1, N=30)
d = G.cart((0,0,0), (1./30.,1,1), (31,1,1))
d1 = G.enforceX(d, 0.5, 0.001, (10,20))
d2 = G.enforceX(d, 0.4, 0.001, (10,20))
s2 = G.map(s2, d2, dir=1)
s1 = G.map(s1, d1, dir=1)
s1 = C.convertArray2Tetra(s1); s1 = G.close(s1)
s2 = C.convertArray2Tetra(s2); s2 = G.close(s2)
x = XOR.booleanUnion(s1, s2, tol=0.)
test.testA([x],2)
