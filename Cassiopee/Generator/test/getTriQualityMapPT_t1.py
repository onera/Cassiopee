# - getTriQualityMap (pyTree) -
import Geom.PyTree as D
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

a = D.sphere((0,0,0), 1, 50); a = C.convertArray2Tetra(a); a = G.close(a)
t = C.newPyTree(['Base',2]); t[2][1][2].append(a)
t = C.initVars(t,'Density',2.); t = C.initVars(t,'centers:cellN',1.)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = G.getTriQualityMap(t)
test.testT(t)
