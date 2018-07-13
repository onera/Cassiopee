# - close (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

a1 = G.cart((0,0,0), (1,1,1), (10,10,1)); a1[0] = 'cart1' 
a2 = G.cart((9+1.e-2,0,0), (1,1,1), (10,10,1)); a2[0] = 'cart2'
a3 = G.cart((0,-5.01,0),(1,1,1),(19,6,1)); a3[0] = 'cart3'
a4 = G.cart((0,9.0001,0),(1,1,1),(10,6,1)); a4[0] = 'cart4'
a5 = G.cart((9.01,9.0002,0),(1,1,1),(10,6,1)); a5[0] = 'cart5'

zones = [a1,a2,a3,a4,a5]
for i in xrange(len(zones)):
    zones[i] = C.addBC2Zone(zones[i], 'match1','BCMatch','imin',\
                            zones[i],'imax',[1,2])
    zones[i] = C.addBC2Zone(zones[i], 'match2','BCMatch','imax',\
                            zones[i],'imin',[1,2])
    zones[i] = C.addBC2Zone(zones[i], 'wall1','BCWall','jmin')
    zones[i] = C.addBC2Zone(zones[i], 'wall2','BCWall','jmax')

zones = C.addVars(zones,'Density'); zones = C.addVars(zones,'centers:cellN')

# Close une liste de maillages structures
zones = G.close(zones, 1e-1)
test.testT(zones, 1)
