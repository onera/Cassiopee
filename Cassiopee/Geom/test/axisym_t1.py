# - axisym -
import Generator as G
import Converter as C
import Geom as D
import KCore.test as test

def F(x):
    return x

# test 1D structure (i-array) + variable
a = D.circle( (0,0,0), 1., 20., 60.)
a = C.addVars(a, 'F'); a = C.initVars(a, 'F', F, ['y'])
a = D.axisym(a, (0,0,0), (0,1,0), 360., 50)
test.testA([a], 1)

# test 2D structure (i,j-array)
a = G.cart((0.,0.,0.), (0.1,0.1,0.2),(10,10,1))
a = D.axisym(a,(1.,0.,0.),(0.,1.,0.),30.,20)
test.testA([a], 2)

# test BAR-array + variable
a = D.circle((0,0,0), 1., 20., 60., 10)
a = C.convertArray2Tetra(a)
a = C.addVars(a, 'F'); a = C.initVars(a, 'F', F, ['y'])
a = D.axisym(a, (0,0,0), (0,1,0), 360., 50)
test.testA([a], 3)

# test TRI-array
a = G.cart((0.,0.,0.), (0.1,0.1,0.2),(10,10,1))
a = C.convertArray2Tetra(a)
a = D.axisym(a,(1.,0.,0.),(0.,1.,0.),30.,30)
test.testA([a], 4)

# test QUAD-array
a = G.cart((0.,0.,0.), (0.1,0.1,0.2),(10,10,1))
a = C.convertArray2Hexa(a)
a = D.axisym(a, (1.,0.,0.),(0.,1.,0.),30.,20)
test.testA([a], 5)
