# - setValue (pyTree) -
import Converter.Internal as Internal
import numpy

node = Internal.createNode('node1', 'DataArray_t', value=12.)

# Set a scalar value in node
Internal.setValue(node, 1.); print(node)
#>> ['node1', array([ 1.]), [], 'DataArray_t']

# Set a numpy array in node
Internal.setValue(node, numpy.zeros(10)); print(node)
#>> ['node1', array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.]), [], 'DataArray_t']

# Set an array as a list
Internal.setValue(node, [1.,12.,13.]); print(node)
#>> ['node1', array([  1.,  12.,  13.]), [], 'DataArray_t']
