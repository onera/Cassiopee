# - setPartialFields (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter
import numpy

a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = C.initVars(a, 'F', 2.)
f1 = Converter.array('F', 5,1,1)
f1[1][:] = 1.
inds = numpy.array([0,1,2], dtype=numpy.int32)
b = C.setPartialFields(a, [f1], [inds], loc='nodes')
t = C.newPyTree(['Base',b])
C.convertPyTree2File(t, 'out.cgns')
