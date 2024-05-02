# - newDataArray (pyTree) -
import Converter.Internal as Internal
import KCore.test as test
import numpy

# Create a DataArray node
n = Internal.newDataArray('CoordinateX', numpy.zeros(10))
test.testT(n, 1)

# Attach it to a parent node
g = Internal.newGridCoordinates()
Internal.newDataArray('CoordinateX', value=numpy.arange(0,10), parent=g)
Internal.newDataArray('CoordinateY', value=numpy.zeros(10), parent=g)
Internal.newDataArray('CoordinateZ', value=numpy.zeros(10), parent=g)
test.testT(n, 2)

# Iterative nodes
n = Internal.newDataArray('FlowSolutionPointers', ['First', 'Second'])
test.testT(n, 3)
