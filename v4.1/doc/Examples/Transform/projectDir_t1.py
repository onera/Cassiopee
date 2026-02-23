# - projectDir (array) -
import Geom as D
import Converter as C
import Generator as G
import Transform as T
import KCore.test as test

# Structure
a = D.sphere((0,0,0), 1., 20)
a = C.initVars(a, 'F', 1)
b = G.cart((1.1,-0.1,-0.1),(0.1,0.1,0.1), (1,5,5))
b = C.initVars(b, 'F', 1)
c = T.projectDir([b], [a], (1.,0,0))
test.testA([a,b]+c, 1)

# Non structure
a = D.sphere((0,0,0), 1., 20)
a = C.initVars(a, 'F', 1)
b = G.cartTetra((1.1,-0.1,-0.1),(0.1,0.1,0.1), (1,5,5))
b = C.initVars(b, 'F', 1)
c = T.projectDir([b], [a], (1.,0,0))
test.testA([a,b]+c, 2)

a = C.convertArray2Tetra(a)
b = G.cartTetra((1.1,-0.1,-0.1),(0.1,0.1,0.1), (1,5,5))
b = C.initVars(b, 'F', 1)
c = T.projectDir([b], [a], (1.,0,0))
test.testA([a,b]+c, 3)

# Structure smoothe
a = D.sphere((0,0,0), 1., 20)
a = C.initVars(a, 'F', 1)
b = G.cart((1.1,-0.1,-0.1),(0.1,0.1,0.1), (1,5,5))
b = C.initVars(b, 'F', 1)
c = T.projectDir([b], [a], (1.,0,0),smooth=1)
test.testA([a,b]+c, 4)
