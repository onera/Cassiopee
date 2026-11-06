# - smoothField (array) -
import Generator as G
import Converter as C
import Transform as T
import numpy
a = G.cartTetra((0,0,0), (1,1,1), (10,10,1))
a = C.initVars(a, '{ro}={x}')

eps = 0.1; niter = 1; type = 0
#Transform.transform._smoothField(a, eps, None, niter, type, ['ro'])
b = T.smoothField(a, eps, niter, type, ['ro'])
T._smoothField(a, eps, niter, type, ['ro'])

eps = numpy.empty( (C.getNPts(a)), dtype=numpy.float64)
eps[:] = 0.1
T._smoothField(a, eps, niter, type, ['ro'])
#Transform.transform._smoothField(a, 0., eps, niter, type, ['ro'])

C.convertArrays2File(a, 'out.plt')
