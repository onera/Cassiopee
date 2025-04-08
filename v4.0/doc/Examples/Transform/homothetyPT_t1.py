# - homothety (pyTree) -
# traitement par zone
import Transform.PyTree as T
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# Structure 3D + champs + CL
a = G.cart((0,0,0), (1,1,1), (10,10,3))
C._addVars(a, 'Density'); C._initVars(a, 'centers:cellN', 1)
a = C.addBC2Zone(a, 'wall','BCWall','imin')
a = C.addBC2Zone(a, 'overlap','BCOverlap','jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin', [1,2,3])
a = T.homothety(a, (0.,0.,0.), 2.)
test.testT(a, 1)

# Structure 2D
a = G.cart((0,0,0), (1,1,1), (10,10,1))
C._addVars(a,'Density'); C._initVars(a, 'centers:cellN', 1)
a = C.addBC2Zone(a, 'wall','BCWall','imin')
a = C.addBC2Zone(a, 'overlap','BCOverlap','jmin')
a = T.homothety(a, (0.,0.,0.), 2.)
test.testT(a, 2)

# TETRA
a = G.cartTetra( (0,0,0), (1,1,1), (10,10,3))
C._addVars(a,'Density'); C._initVars(a, 'centers:cellN', 1)
a = T.homothety(a, (0.,0.,0.), 2.)
test.testT(a, 3)

# HEXA
a = G.cartHexa( (0,0,0), (1,1,1), (10,10,3))
C._addVars(a,'Density'); C._initVars(a, 'centers:cellN', 1)
a = T.homothety(a, (0.,0.,0.), 2.)
test.testT(a,4)

# TRI
a = G.cartTetra( (0,0,0), (1,1,1), (10,10,1))
C._addVars(a,'Density'); C._initVars(a, 'centers:cellN', 1)
a = T.homothety(a, (0.,0.,0.), 2.)
test.testT(a,5)

# QUAD
a = G.cartHexa( (0,0,0), (1,1,1), (10,10,1))
C._addVars(a,'Density'); C._initVars(a, 'centers:cellN', 1)
a = T.homothety(a, (0.,0.,0.), 2.)
test.testT(a,6)
