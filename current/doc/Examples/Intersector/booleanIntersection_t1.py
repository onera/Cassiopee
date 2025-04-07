# - boolean intersection (array) -
import Intersector as XOR
import Generator as G
import Converter as C
import Geom as D
import KCore.test as test

# Intersection de deux spheres
s1 = D.sphere((0,0,0), 1, N=20)
s2 = D.sphere((0.,1.,0.), 1, N=30)
s1 = C.convertArray2Tetra(s1); s1 = G.close(s1)
s2 = C.convertArray2Tetra(s2); s2 = G.close(s2)

x = XOR.booleanIntersection(s1, s2, tol=0.)
test.testA([x],1)

# Intersection avec resserements locaux
s1 = D.sphere((0,0,0), 1, N=20)
s2 = D.sphere((0.,1.,0.), 1, N=30)
d = G.cart((0,0,0), (1./30.,1,1), (31,1,1))
d = G.enforceX(d, 0.5, 0.001, (10,20))
s2 = G.map(s2, d, dir=1)
s1 = G.map(s1, d, dir=1)
s1 = C.convertArray2Tetra(s1); s1 = G.close(s1)
s2 = C.convertArray2Tetra(s2); s2 = G.close(s2)

x = XOR.booleanIntersection(s1, s2, tol=0.)
test.testA([x], 2)
