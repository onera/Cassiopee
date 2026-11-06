# - cart2Cyl (array) -
import Transform as T
import Generator as G
import KCore.test as test

# test (X,Y,Z) -> (R,Theta,Z)
a = G.cylinder((0.,0.,0.), 0.5, 1., 0., 360, 1., (360,20,10))
a = T.cart2Cyl(a, (0.,0.,0.),(0,0,1))
test.testA(a,1)

# (X,Y,Z)-> (X,R,Theta)
a = G.cylinder((0.,0.,0.), 0.5, 1., 0., 360, 1., (360,20,10))
a = T.rotate(a,(0,0,0),(0,1,0),90.)
a = T.cart2Cyl(a, (0.,0.,0.),(1,0,0))
test.testA(a,2)

# (X,Y,Z)-> (Theta,Y,R)
a = G.cylinder((0.,0.,0.), 0.5, 1., 0., 360, 1., (360,20,10))
a = T.rotate(a,(0,0,0),(1,0,0),90.)
a = T.cart2Cyl(a, (0.,0.,0.),(0,1,0))
test.testA(a,3)
