# - compress/uncompress using zfp -
import numpy
import Compressor.zfp as zfp

# Test 1
a = numpy.arange(0.,1.,0.01, dtype=numpy.float64)
print(64*a.size)
print(a)

za = zfp.pack(a, accuracy=1.e-10)
#print(za[1].size)
#print(za)

b = zfp.unpack(za, accuracy=1.e-10)
print(b)

# Test 2
print("fortran order")
a = numpy.zeros( (5,3,3), dtype=numpy.float64, order='F')
a[:,0,0] = 1.
print('original a')
print(a)

za = zfp.pack(a, accuracy=1.e-8)
b = zfp.unpack(za, accuracy=1.e-8)
print("deflated a")
print(b)

# Test 3
print("c order")
a = numpy.zeros( (5,3,3), dtype=numpy.float64)
a[:,0,0] = 1.
print("original a")
print(a)

za = zfp.pack(a, accuracy=1.e-8)
b = zfp.unpack(za, accuracy=1.e-8)
print("deflated a")
print(b)
