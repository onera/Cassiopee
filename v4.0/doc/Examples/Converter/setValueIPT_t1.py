# - setValue (pyTree) -
import Converter.Internal as Internal
import KCore.test as test
import numpy

node = Internal.createNode('node1', 'DataArray_t', value=12.)

# Set a scalar value in node
Internal.setValue(node, 1.)
test.testO(node, 1)

# Set a numpy array in node
Internal.setValue(node, numpy.zeros(10))
test.testO(node, 2)

# Set an array as a list
Internal.setValue(node, [1.,12.,13.])
test.testO(node, 3)

# Set a string
Internal.setValue(node, 'toto')
test.testO(node, 4)

# Set a list of strings
Internal.setValue(node, ['toto', 'tata'])
test.testO(node, 5)
