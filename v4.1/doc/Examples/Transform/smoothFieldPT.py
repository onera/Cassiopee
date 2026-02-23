# - smoothField (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Transform.PyTree as T
import numpy

a = G.cartTetra((0,0,0), (1,1,1), (10,10,1))
a = C.initVars(a, '{Density}={CoordinateX}')

eps = 0.1; niter = 1; type = 0
#b = T.smoothField(a, eps, niter, type, ['ro'])
T._smoothField(a, eps, niter, type, ['Density'])

eps = numpy.empty( (C.getNPts(a)), dtype=numpy.float64)
eps[:] = 0.1
T._smoothField(a, eps, niter, type, ['Density'])

C.convertPyTree2File(a, 'out.cgns')
