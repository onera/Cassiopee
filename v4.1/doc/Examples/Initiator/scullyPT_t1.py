# - initScully (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Initiator.PyTree as I
import KCore.test as test

# Pas de variables conservatives dans l'arbre
NI = 200; NJ = 200
HI = 1./(NI-1); HJ = 1./(NJ-1)
a = G.cart((0.,0.,0.), (HI,HJ,1.), (NI,NJ,2))
a = C.addBC2Zone(a, 'wall', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'overlap', 'BCOverlap', 'jmin')
t = C.newPyTree(['Base']); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'Mach', 0.8)
# constant entropy outward + constant temperature in the vortex
t = I.initScully(t,(0.5,0.5), -0.2, 0.05, 0.8, 0, loc='centers')
test.testT(t,1)

# constant entropy
t = C.newPyTree(['Base']); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'Mach', 0.8)
t = I.initScully(a,(0.5,0.5), -0.2, 0.05, 0.8, 1, loc='centers')
test.testT(t,2)
