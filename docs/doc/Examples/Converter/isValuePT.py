# - isValue (pyTree) -
import Converter.Internal as Internal
import numpy

# Check a scalar value
node = Internal.createNode('node1', 'DataArray_t', value=1.)
print(Internal.isValue(node, 1.))
#>> True

# Check a numpy array values
node = Internal.createNode('node1', 'DataArray_t', value=numpy.zeros(10))
print(Internal.isValue(node, numpy.zeros(10)))
#>> True

# Check a string value
node = Internal.createNode('node1', 'DataArray_t', value='toto')
print(Internal.isValue(node, 'toto'))
#>> True
