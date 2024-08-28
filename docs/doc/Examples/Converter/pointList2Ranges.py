# - pointList2Ranges -
import Converter as C
import numpy

a = numpy.array([0,3,6,9], dtype=numpy.int32)

r = C.converter.pointList2Ranges(a, 3,3,3)
print(r)
