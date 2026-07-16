# - consSmooth (array) -
import Transform as T
import Geom as D
import Generator as G
import Converter as C
import Post as P
import KCore.test as test

l1 = D.line((0.,0.,0.), (0.,1.,0.), N=5)
l2 = D.line((0.,1.,0.), (1.,1.,0.), N=5)
l3 = D.line((1.,1.,0.), (1.,0.,0.), N=5)
l4 = D.line((1.,0.,0.), (0.,0.,0.), N=5)

# struct 1D closed curve
a = T.join([l1,l2,l3,l4])
a = T.consSmooth(a, sweeps=5)
test.testA([a], 1)

# struct 1D open curve
b = T.join([l1,l2,l3])
b = T.consSmooth(b, sweeps=5)
test.testA([b], 2)

# bar1D closed curve
c = C.convertArray2Tetra(a)
c = T.consSmooth(c, sweeps=5)
test.testA([c], 3)

# bar1D open curve
d = C.convertArray2Tetra(b)
d = T.consSmooth(d, sweeps=5)
test.testA([d], 4)

# 3D : cube
e = G.cart((0.0, 0.0, 0.0), (0.1, 0.1, 0.1), (11, 11, 11))
e = C.convertArray2Hexa(e)
e = P.exteriorFaces(e)
e = C.convertArray2Tetra(e, split='withBarycenters')
e = G.close(e)
e = T.consSmooth(e, sweeps=5, omega=0.1)
test.testA([e], 5)

# 3D : open cube
N = 11
h = 0.1
a = G.cart((0.0, 0.0, 0.0), (h, h, 1.), (N, N, 1))
b = G.cart((0.0, 0.0, 0.0), (1., h, h), (1, N, N))
c = G.cart((1.0, 0.0, 0.0), (1., h, h), (1, N, N))
d = G.cart((0.0, 0.0, 0.0), (h, 1., h), (N, 1, N))
e = G.cart((0.0, 1.0, 0.0), (h, 1., h), (N, 1, N))
f = [a, b, c, d, e]
f = C.convertArray2Hexa(f)
f = C.convertArray2Tetra(f, split='withBarycenters')
f = T.join(f)
f = G.close(f, 1e-6)
f = T.reorder(f, (1,))
f = T.consSmooth(f, sweeps=5, omega=0.1)
test.testA([f], 6)
