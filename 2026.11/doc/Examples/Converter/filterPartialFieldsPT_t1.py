# - filterPartialFields (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter
import numpy
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,1))
C._initVars(a, 'F', 2.)
f1 = Converter.array('F,G', 5,1,1)
f1 = Converter.initVars(f1,'F=10.')
f1 = Converter.initVars(f1,'G=1.')
f2 = Converter.array('F,G', 5,1,1)
f2 = Converter.initVars(f2,'F=-10.')
f2 = Converter.initVars(f2,'G=0.1')
inds = numpy.array([0,1,2], dtype=Internal.E_NpyInt)
C._filterPartialFields(a, [f1,f2], inds, loc='nodes',filterName='G')
test.testT(a)

a = G.cart((0,0,0), (1,1,1), (10,10,1))
C._initVars(a, 'F', 2.)
C._initVars(a, 'G', 1.)
f1 = Converter.array('F,G', 5,1,1)
f1 = Converter.initVars(f1,'F=10.')
f1 = Converter.initVars(f1,'G=1.')
f2 = Converter.array('F,G', 5,1,1)
f2 = Converter.initVars(f2,'F=-10.')
f2 = Converter.initVars(f2,'G=0.1')
inds = numpy.array([0,1,2], dtype=Internal.E_NpyInt)
C._filterPartialFields(a, [f1,f2], inds, loc='nodes',filterName='G')
test.testT(a,2)
