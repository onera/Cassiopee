# - axisym (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test

def F(x): return x

# test 1D structure (i-array) + variable en noeuds
a = D.circle((0,0,0), 1., 20., 60.)
C._addVars(a, 'F'); C._initVars(a, 'F', F, ['CoordinateY'])
a = D.axisym(a, (0,0,0), (0,1,0), 360., 50)
t = C.newPyTree(['Base', 2, a])
test.testT(t, 1)

# test 2D structure (i,j-array) + variable en centres
a = G.cart((0.,0.,0.), (0.1,0.1,0.2),(10,10,1))
a = C.addVars(a, 'centers:F')
a = D.axisym(a,(1.,0.,0.),(0.,1.,0.),30.,20)
t = C.newPyTree(['Base', 3, a])
test.testT(t, 2)

# test BAR-array + variable
a = D.circle((0,0,0), 1., 20., 60., 10)
a = C.convertArray2Tetra(a)
a = C.addVars(a, 'F'); a = C.initVars(a, 'F', F, ['CoordinateY'])
a = D.axisym(a, (0,0,0), (0,1,0), 360., 50)
t = C.newPyTree(['Base', 2, a])
test.testT(t, 3)

# test TRI-array
a = G.cart((0.,0.,0.), (0.1,0.1,0.2),(10,10,1))
a = C.convertArray2Tetra(a)
a = D.axisym(a,(1.,0.,0.),(0.,1.,0.),30.,30)
t = C.newPyTree(['Base', 3, a])
test.testT(t, 4)

# test QUAD-array
a = G.cart((0.,0.,0.), (0.1,0.1,0.2),(10,10,1))
a = C.convertArray2Hexa(a)
a = D.axisym(a, (1.,0.,0.),(0.,1.,0.),30.,20)
t = C.newPyTree(['Base', 3, a])
test.testT(t, 5)
