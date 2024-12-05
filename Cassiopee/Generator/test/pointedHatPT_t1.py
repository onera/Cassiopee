# - pointedHat (pyTree) -
import Geom.PyTree as D
import Generator.PyTree as G
import KCore.test as test
import Converter.PyTree as C

# i-array
c = D.circle( (0,0,0), 1., 360., 0., 100)
c = C.addVars(c,'centers:cellN'); c = C.addVars(c,'Density')
surf = G.pointedHat(c,(0.,0.,1.))
test.testT(surf,1)

# ij-array
a = G.cylinder((0,0,0), 1., 1.5, 360., 0., 1., (100,30,1))
a = C.addVars(a,'centers:cellN'); a = C.addVars(a,'Density')
surf = G.pointedHat(a,(0.,0.,1.))
test.testT(surf,2)
