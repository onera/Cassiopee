# - distance2Walls (array) -
# test : 2D
import Dist2Walls
import Generator as G
import KCore.test as test
import Geom as D

a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(21,21,1))
cyl = D.circle((0.5,0.5,0),1.,N=10)
dist = Dist2Walls.distance2Walls(a, [cyl], loc='nodes', type='ortho',
                                 signed=0, dim=2)
test.testA([dist],1)
dist = Dist2Walls.distance2Walls(a, [cyl], loc='nodes', type='ortho',
                                 signed=1, dim=2)
test.testA([dist],2)
