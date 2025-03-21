# - extractPlane (array) -
import Converter as C
import Post as P
import Transform as T
import Generator as G
import KCore.test as test

# Cree un cylindre
m = G.cylinder((0,0,0), 1, 5, 0., 360., 10., (30,30,30))
m = T.rotate(m , (0,0,0), (1,0,0), 35.)
for i in [2,3]:
    a = P.extractPlane([m], (0.5, 1., 0., 1),i)
    test.testA([a], i)

# cas non structure tetra
m = G.cylinder((0,0,0), 1, 5, 0., 360., 10., (30,30,30))
m = T.rotate(m , (0,0,0), (1,0,0), 35.)
m = C.convertArray2Tetra(m)
a = P.extractPlane([m], (0.5, 1., 0., 1))
test.testA([a], 1)
