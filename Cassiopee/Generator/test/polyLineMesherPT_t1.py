# - polyLineMesher (pyTree)-
import Generator.PyTree as G
import Geom.PyTree as D
import KCore.test as test
import Converter.PyTree as C

a = D.polyline([(0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.),(0.,0.,0.)])
a = C.convertArray2Tetra(a)
a = G.close(a, 1.e-2)

# Donnees
h = 0.05; hf = 0.0001; density = 200
res = G.polyLineMesher(a, h, hf, density)
zones = res[0]; h = res[1]; density = res[2]
tb = C.newPyTree(['Base']); tb[2][1][2] += zones
test.testT(tb)
