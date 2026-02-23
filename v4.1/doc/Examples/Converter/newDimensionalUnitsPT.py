# - newDimensionalUnits (pyTree) -
import Converter.Internal as Internal
import numpy

# Create a DimensionalUnits node
n = Internal.newDimensionalUnits(massUnit='Kilogram')
Internal.printTree(n)
#>> ['DimensionalUnits',array(shape=(32, 5),dtype='|S1',order='F'),[0 son],'DimensionalUnit_t']

# Attach it to a parent node
d = Internal.newDataArray('CoordinateX', numpy.zeros(10))
Internal.newDimensionalUnits(lengthUnit='Meter', parent=d)
