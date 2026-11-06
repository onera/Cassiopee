# - smooth (pyTree) -
import Transform.PyTree as T
import Geom.PyTree as D
import Converter.PyTree as C
import KCore.test as test

# TRI
a = D.sphere((0.,0.,0.), 1)
b = C.convertArray2Tetra(a)
C._addVars(b, 'Density'); C._initVars(b,'centers:cellN',1.)
b = T.smooth(b, eps=0.5, niter=10)
test.testT(b,1)

# QUAD
b = C.convertArray2Hexa(a)
C._addVars(b,'Density'); C._initVars(b,'centers:cellN',1.)
b = T.smooth(b, eps=0.5, niter=10)
test.testT(b,2)
