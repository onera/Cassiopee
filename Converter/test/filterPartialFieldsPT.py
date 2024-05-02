# - filterPartialFields (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter
import numpy

a = G.cart((0,0,0), (1,1,1), (10,10,1))
C._initVars(a, 'F', 2.)
f1 = Converter.array('F,G', 5,1,1)
f1 = Converter.initVars(f1,'F=10.')
f1 = Converter.initVars(f1,'G=1.')
f2 = Converter.array('F,G', 5,1,1)
f2 = Converter.initVars(f2,'F=-10.')
f2 = Converter.initVars(f2,'G=0.1')

inds = numpy.array([0,1,2], dtype=numpy.int32)
print(inds)
C._filterPartialFields(a, [f1,f2], inds, loc='nodes',filterName='G')
t = C.newPyTree(['Base',a])
C.convertPyTree2File(t, 'out.cgns')
