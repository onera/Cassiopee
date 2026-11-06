# - zip (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a1 = G.cart((0,0,0), (1,1,1), (10,10,1)); a1[0] = 'cart1'
a2 = G.cart((9+1.e-2,0,0), (1,1,1), (10,10,1)); a2[0] = 'cart2'
a3 = G.cart((0,-5.01,0),(1,1,1),(19,6,1)); a3[0] = 'cart3'
a4 = G.cart((0,9.0001,0),(1,1,1),(10,6,1)); a4[0] = 'cart4'
a5 = G.cart((9.01,9.0002,0),(1,1,1),(10,6,1)); a5[0] = 'cart5'

zones = [a1,a2,a3,a4,a5]
for z in zones:
    C._addBC2Zone(z, 'match1','BCMatch','imin', z,'imax',[1,2])
    C._addBC2Zone(z, 'match2','BCMatch','imax', z,'imin',[1,2])
    C._addBC2Zone(z, 'wall1','BCWall','jmin')
    C._addBC2Zone(z, 'wall2','BCWall','jmax')

zones = C.addVars(zones, 'Density'); zones = C.addVars(zones, 'centers:cellN')

# Close une liste de maillages sur leurs frontieres
zones2 = G.zip(zones, 1e-1)
test.testT(zones2, 1)

zones2 = C.convertArray2Hexa(zones)
zones2 = G.zip(zones2, 1e-1)
test.testT(zones2, 2)

zones2 = C.convertArray2Tetra(zones)
zones2 = G.zip(zones2, 1e-1)
test.testT(zones2, 3)

zones2 = C.convertArray2NGon(zones)
zones2 = G.zip(zones2, 1e-1)
test.testT(zones2, 4)
