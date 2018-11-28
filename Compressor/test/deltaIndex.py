# - deltaIndex -
import numpy
import Compressor

# Liste des indexes de reference
indRef = numpy.array([1,2,3,4,5], dtype='int32')

# Liste des indexes a comparer a la reference
index = numpy.array([1,2,3,4], dtype='int32')

delta = Compressor.deltaIndex(index, indRef)
print(delta)
