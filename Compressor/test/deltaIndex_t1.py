# - deltaIndex -
import numpy
import Compressor
import Converter.Internal as Internal
import KCore.test as test

# Liste des indexes de reference
indRef = numpy.array([1,2,3,4,5], dtype=Internal.__E_NPY_INT__)

# Liste des indexes a comparer a la reference
index = numpy.array([1,2,3,4], dtype=Internal.__E_NPY_INT__)
delta = Compressor.deltaIndex(index, indRef)
test.testO(delta)
