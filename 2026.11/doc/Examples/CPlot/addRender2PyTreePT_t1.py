# - addRender2PyTree (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import CPlot.PyTree as CPlot
import KCore.test as test

a = G.cart( (0,0,0), (1,1,1), (10,10,1) )
a = C.addBC2Zone(a,'wall','BCWall','imin')
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
t = C.newPyTree(['Base', 2, a])
t = CPlot.addRender2PyTree(t, slot=0, posCam=(1,1,1), posEye=(10,1,1),
                           mode='Solid', niso=10, isoScales=[0, 10, 0., 1.],
                           isoEdges=0.1, isoLight=1, colormap='Blue2Red')
test.testT(t, 1)
