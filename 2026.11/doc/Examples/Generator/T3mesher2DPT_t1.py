# - T3mesher2D (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test

a = D.circle((0,0,0), 1, N=50)
a = C.convertArray2Tetra(a); a = G.close(a)
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
b = G.T3mesher2D(a, triangulateOnly=0, grading=1.2, metricInterpType=1) # geometric metric interpolation
test.testT(b)
