# - setPartialFields1 (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import KCore.test as test
import numpy

a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = C.initVars(a, 'F', 2.)
a = C.initVars(a, 'G', 2.)
f1 = numpy.array([0.,0.,0.,0.,0.,0.,0.,0.,0.,0.], dtype=numpy.float64)
inds = numpy.array([0,1,2,3,4,5,6,7,8,9], dtype=Internal.E_NpyInt)
t = C.newPyTree(['Base', a])
t = C.setPartialFields1(t, [ [[],f1] ] , [inds], loc='nodes')
test.testT(t, 1)
