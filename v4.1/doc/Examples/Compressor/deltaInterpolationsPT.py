# - deltaInterpolations -
import numpy
import Compressor.PyTree as Compressor

# Liste des donnees d interpolations de reference
rcvIndices = numpy.array([1,2,3], dtype='int32')
donorIndices = numpy.array([1,5,6], dtype='int32')
periodicity = numpy.array([100,100,100], dtype='int32')
coefs1 = numpy.array([0.,0.,0.,0.,0.5,0.5,0.5], dtype='float')
coefs2 = numpy.array([0.,0.,0.,0.,0.5,0.5,0.5], dtype='float')
coefs3 = numpy.array([0.,0.,0.,0.,0.5,0.5,0.5], dtype='float')
coefficients = [coefs1,coefs2,coefs3]
indRef = [rcvIndices,donorIndices,periodicity,coefficients]

# Liste des donnees d interpolations a comparer a la reference
rcvIndices = numpy.array([2,3,4], dtype='int32')
donorIndices = numpy.array([5,6,7], dtype='int32')
periodicity = numpy.array([101,100,100], dtype='int32')
coefs2 = numpy.array([0.,0.,0.,0.,0.5,0.5,0.5], dtype='float')
coefs3 = numpy.array([0.,0.,0.,0.,0.5,0.5,0.5], dtype='float')
coefs4 = numpy.array([0.,0.,0.,0.,0.5,0.5,0.5], dtype='float')
coefficients = [coefs2,coefs3,coefs4]
index = [rcvIndices,donorIndices,periodicity,coefficients]

# Liste des indexes a comparer a la reference
delta = Compressor.deltaInterpolations(index, indRef, loc='cell')
print(delta)
