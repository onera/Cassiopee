# - compress/uncompress using sz -
import numpy
import Compressor.sz as sz

# Test 1 - flat C ordered
a = numpy.arange(0.,1.,0.01, dtype=numpy.float64)
#print(64*a.size)
#print(a)
za = sz.pack(a, {'relBoundRatio':1.e-8})
#print(za[1].size)
#print(za)
b = sz.unpack(za, {'relBoundRatio':1.e-8})
#print(b)

# Test 2 - shaped fortran ordered
a = numpy.zeros( (5,3,3), dtype=numpy.float64, order='F')
a[:,0,0] = 1.
print(a)

za = sz.pack(a, {'relBoundRatio':1.e-8})
#print(za)
b = sz.unpack(za, {'relBoundRatio':1.e-8})
print(b)

# Test 3 - shaped C ordered
a = numpy.zeros( (5,3,3), dtype=numpy.float64)
a[:,0,0] = 1.
#print(a)

za = sz.pack(a, {'relBoundRatio':1.e-8})
#print(za)
b = sz.unpack(za, {'relBoundRatio':1.e-8})
#print(b)
