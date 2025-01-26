# - densify (array) -
import Generator as G
import Converter as C
import Geom as D
import KCore.test as test

a = D.polyline([(0.,0.,0.),(1.,1.,0.),(2.,0.,0.),(2.5,0.5,0.),(3.5,0.,0.)])
a = C.convertArray2Tetra(a)
c = C.initVars(a,'F',1.)
c = G.densify(a, 1.05)
test.testA([c], 1)

a = D.polyline([(0.,0.,0.),(1.,1.,0.),(2.,0.,0.),(2.5,0.5,0.),(3.5,0.,0.)])
c = C.initVars(a,'F',1.)
c = G.densify(a, 1.05)
test.testA([c], 2)
