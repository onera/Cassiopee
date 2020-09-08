# - compressCellN (pyTree) -
import Compressor.PyTree as Compressor
import Compressor.compressor as Co
import Generator.PyTree as G
import Converter.PyTree as C
import numpy

a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._initVars(a, '{centers:cellN}=1.')
Compressor._compressCellN(a)
C.convertPyTree2File(a, 'out.cgns')

# Test 2
print("fortran order")
a = numpy.zeros( (5,3,3), dtype=numpy.float64, order='F')
a[:,0,0] = 1.
print('original a')
print(a)

za = Co.compressCellN(a)
b = Co.uncompressCellN(za)
print("deflated a")
print(b)

# Test 3
print("c order")
a = numpy.zeros( (5,3,3), dtype=numpy.float64)
a[:,0,0] = 1.
print("original a")
print(a)

za = Co.compressCellN(a)
b = Co.uncompressCellN(za)
print("deflated a")
print(b)
