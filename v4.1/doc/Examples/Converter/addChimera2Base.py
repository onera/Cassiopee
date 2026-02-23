# - addChimera2Base (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cylinder((0,0,0), 1., 1.5, 0., 360., 1., (80,30,2))
t = C.newPyTree(['Base', a])
t[2][1] = C.addChimera2Base(t[2][1], 'XRayTol', 1.e-6)
t[2][1] = C.addChimera2Base(t[2][1], 'XRayDelta', 0.1)
t[2][1] = C.addChimera2Base(t[2][1], 'DoubleWallTol', 100.)
t[2][1] = C.addChimera2Base(t[2][1], 'Priority', 1)
C.convertPyTree2File(t, 'out.cgns')
