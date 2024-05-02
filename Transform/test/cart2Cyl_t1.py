# - cart2Cyl (array) -
import Transform as T
import Generator as G
import Converter as C
import Geom as D
import KCore.test as test

# structure 3D
a = G.cylinder((0.,0.,0.), 0.5, 1., 0., 359, 1., (360,20,10))
a = T.cart2Cyl(a, (0.,0.,0.),(0,0,1), thetaShift=0.)
test.testA([a], 1)

# structure 2D
a = D.cone((0,0,0), 1. , 0.5, 1.)
axis = (0,0,1)
a = T.cart2Cyl(a, (0.,0.,0.),axis, thetaShift=0.)
test.testA([a],2)

# HEXA
a = G.cylinder((0.,0.,0.), 0.5, 1., 0., 359, 1., (360,20,10))
a = C.convertArray2Hexa(a)
a = T.cart2Cyl(a, (0.,0.,0.),(0,0,1), thetaShift=0.)
test.testA([a],3)

# TETRA
a = G.cylinder((0.,0.,0.), 0.5, 1., 0., 359, 1., (360,20,10))
a = C.convertArray2Tetra(a)
a = T.cart2Cyl(a, (0.,0.,0.),(0,0,1), thetaShift=0.)
test.testA([a],4)
