# - initVisbal (pyTree) -
import Generator.PyTree as G
import Initiator.PyTree as I
import KCore.test as test
import Converter.PyTree as C

# Pas de variables conservatives dans l arbre
NI = 100; NJ = 100
HI = 50./(NI-1); HJ = 50./(NJ-1)
a = G.cart( (0.,0.,0.), (HI,HJ,1.), (NI,NJ,2))
a = C.addBC2Zone(a, 'wall', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'overlap', 'BCOverlap', 'jmin')
t = C.newPyTree(['Base']); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'Mach', 0.8)
t = C.initVars(t,'F', 1.); t = C.initVars(t,'centers:G', 2.)
t = I.initVisbal(t, position=(7.,7.), Gamma=2., MInf=0.8, loc='centers')
test.testT(t,1)
