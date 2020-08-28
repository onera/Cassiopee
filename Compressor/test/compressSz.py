# - compress/uncompress using sz -
import numpy
import Compressor.sz as sz

a = numpy.arange(0.,1.,0.01, dtype=numpy.float64)
print(64*a.size)
print(a)

za = sz.pack(a, {'relBoundRatio':1.e-8})
#print(za[1].size)
print(za)

b = sz.unpack(za, {'relBoundRatio':1.e-8})
print(b)
