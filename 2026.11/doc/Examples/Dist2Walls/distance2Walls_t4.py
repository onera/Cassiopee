# - distance2Walls (array) -
# test : array structure, bodies non structures
import Dist2Walls
import Generator as G
import KCore.test as test
import Geom as D
import Converter as C

a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(21,21,11))
sphere = D.sphere((0.5,0.5,0.5),0.2,N=20)
sphere1 = C.convertArray2Hexa(sphere)
sphere = D.sphere((1.5,0.5,0.5),0.2,N=20)
sphere2 = C.convertArray2Tetra(sphere)
dist = Dist2Walls.distance2Walls(a, [sphere1,sphere2])
test.testA([dist],1)
