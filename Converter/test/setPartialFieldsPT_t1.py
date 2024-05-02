# - setPartialFields (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter
import KCore.test as test
import numpy

a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = C.initVars(a, 'F', 2.)
f1 = Converter.array('F', 5,1,1)
f1[1][:] = 1.
inds = numpy.array([1,2,3], dtype=Internal.E_NpyInt)
b = C.setPartialFields(a, [f1], [inds], loc='nodes',startFrom=1)
t = C.newPyTree(['Base',b])
test.testT(t, 1)
