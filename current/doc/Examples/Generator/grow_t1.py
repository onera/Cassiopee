# - grow (array) -
import Converter as C
import Generator as G
import Geom as D
import Transform as T
import KCore.test as test

# Structure
a = D.sphere( (0,0,0), 1., 50 )
a = C.initVars(a, 'F', 1.)
n = G.getNormalMap(a)
n = C.center2Node(n)
n[1] = n[1]*100.
b = G.grow(a, n)
test.testA([b], 1)

# Non structure
a = D.sphere( (0,0,0), 1., 50 )
a = C.convertArray2Hexa(a)
a = C.initVars(a, 'F', 1.)
n = G.getNormalMap(a)
n = C.center2Node(n)
n[1] = n[1]*100.
b = G.grow(a, n)
test.testA([b], 2)

# Structure 1D
a = D.line( (0,0,0), (1,0,0), 20 )
c = T.addkplane(a)
n = G.getNormalMap(c)
n = C.center2Node(n)
n[1] = n[1]*100.
b = G.grow(a, n)
test.testA([b], 3)

# BAR
a = D.line( (0,0,0), (1,0,0), 20 )
a = C.convertArray2Tetra(a); a = G.close(a)
a = T.addkplane(a)
n = G.getNormalMap(a)
n = C.center2Node(n)
n[1] = n[1]*100.
b = G.grow(a, n)
test.testA([b], 4)
