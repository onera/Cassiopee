# - polyline (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import KCore.test as test

a = D.polyline([(0.,0.,0.),(1.,1.,0.),(2.,0.,0.)])
t = C.newPyTree(['Base',1,a])
test.testT(t, 1)
