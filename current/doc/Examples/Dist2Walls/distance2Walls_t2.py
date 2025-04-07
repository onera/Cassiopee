# - distance2Walls (array) -
# test des types d array
import Dist2Walls
import Generator as G
import KCore.test as test
import Geom as D

sphere = D.sphere((0.5,0.5,0.5),0.2,20)
# sur une liste de zones
a1 = G.cart((0.,0.,0.),(0.1,0.1,0.1),(11,11,11))
dist = Dist2Walls.distance2Walls([a1], [sphere])
test.testA(dist,1)

# sur un array HEXA
a1 = G.cartHexa((0.,0.,0.),(0.1,0.1,0.1),(11,11,11))
dist = Dist2Walls.distance2Walls([a1], [sphere])
test.testA(dist,2)

# sur un array TETRA
a1 = G.cartTetra((0.,0.,0.),(0.1,0.1,0.1),(11,11,11))
dist = Dist2Walls.distance2Walls([a1], [sphere])
test.testA(dist,3)

# MIXTE
a1 = G.cart((0.,0.,0.),(0.1,0.1,0.1),(11,11,11))
a2 = G.cartHexa((1.,0.,0.),(0.1,0.1,0.1),(11,11,11))
a3 = G.cartTetra((-1.,0.,0.),(0.1,0.1,0.1),(11,11,11))
dist = Dist2Walls.distance2Walls([a1], [sphere])
test.testA(dist,4)
