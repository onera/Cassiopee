# - compress/uncompress using zfp -
import numpy
import Compressor.zfp as zfp

a = numpy.arange(0.,1.,0.01, dtype=numpy.float64)
print(64*a.size)
print(a)

za = zfp.pack(a, accuracy=1.e-10)
print(za[1].size)
#print(za)

b = zfp.unpack(za, accuracy=1.e-10)
print(b)
