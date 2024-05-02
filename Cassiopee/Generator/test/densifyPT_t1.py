# - densify (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test

a = D.polyline([(0.,0.,0.),(1.,1.,0.),(2.,0.,0.),(2.5,0.5,0.),(3.5,0.,0.)])
a = C.convertArray2Tetra(a)
a = C.addVars(a,'Density'); a = C.initVars(a,'centers:cellN',1)
c = G.densify(a, 0.1)
test.testT(c,1)

a = D.polyline([(0.,0.,0.),(1.,1.,0.),(2.,0.,0.),(2.5,0.5,0.),(3.5,0.,0.)])
a = C.addVars(a,'Density'); a = C.initVars(a,'centers:cellN',1)
c = G.densify(a, 1.05)
c = C.initVars(a,'F',1.)
test.testT([c], 2)
