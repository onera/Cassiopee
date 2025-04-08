# - grow (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D
import KCore.test as test
import Converter.Internal as Internal

a = D.sphere( (0,0,0), 1., 50 )
a = G.getNormalMap(a)
a = C.center2Node(a, Internal.__FlowSolutionCenters__)
a = C.rmVars(a, Internal.__FlowSolutionCenters__)
a = C.initVars(a, 'Density',2.); a = C.initVars(a, 'centers:cellN',1.)
a = C.addBC2Zone(a, 'overlap', 'BCOverlap', 'jmax')
a = C.addBC2Zone(a, 'wall', 'BCWall', 'jmin')
b = G.grow(a, ['sx','sy','sz'])
t = C.newPyTree(['Base1',2,'Base2',3])
t[2][1][2].append(a); t[2][2][2].append(b)
test.testT(b)
