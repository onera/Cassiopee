# - newDataArray (pyTree) -
import Converter.Internal as Internal
import numpy

# Create a DataArray node
n = Internal.newDataArray('CoordinateX', numpy.zeros(10)); print(n)
#>> ['CoordinateX', array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.]), [], 'DataArray_t']

# Attach it to a parent node
g = Internal.newGridCoordinates()
Internal.newDataArray('CoordinateX', value=numpy.arange(0,10), parent=g)
Internal.newDataArray('CoordinateY', value=numpy.zeros(10), parent=g)
Internal.newDataArray('CoordinateZ', value=numpy.zeros(10), parent=g); Internal.printTree(g)
#>> ['GridCoordinates',None,[3 sons],'GridCoordinates_t']
#>>   |_['CoordinateX',array(shape=(10,),dtype='int64',order='F'),[0 son],'DataArray_t']
#>>   |_['CoordinateY',array(shape=(10,),dtype='float64',order='F'),[0 son],'DataArray_t']
#>>   |_['CoordinateZ',array(shape=(10,),dtype='float64',order='F'),[0 son],'DataArray_t']
