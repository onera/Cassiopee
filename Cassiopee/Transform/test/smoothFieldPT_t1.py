# - smoothField (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Transform.PyTree as T
import KCore.test as test
import numpy

a = G.cartTetra((0,0,0), (1,1,1), (10,10,1))
a = C.initVars(a, '{Density}={CoordinateX}')
a = C.initVars(a, '{centers:Density2}={centers:CoordinateX}')

eps = 0.1; niter = 1; type = 0
T._smoothField(a, eps, niter, type, ['Density'])
test.testT(a, 1)

#eps = 0.1; niter = 1; type = 0
#T._smoothField(a, eps, niter, type, ['centers:Density2'])
#C.convertPyTree2File(a, 'out.cgns')
#test.testT(a, 2)

eps = numpy.empty((C.getNPts(a)), dtype=numpy.float64)
eps[:] = 0.1
T._smoothField(a, eps, niter, type, ['Density'])
test.testT(a, 3)
