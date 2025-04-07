# - addRender2Zone (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import CPlot.PyTree as CPlot
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,1))
a = C.addBC2Zone(a,'wall','BCWall','imin')
C._initVars(a,'F',1.); C._initVars(a,'centers:G',2.)
a = CPlot.addRender2Zone(a, material='Glass', color='Blue')
t = C.newPyTree(['Base', 2, a])
test.testT(t, 1)
