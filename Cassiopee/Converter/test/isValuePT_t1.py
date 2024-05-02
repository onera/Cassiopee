# - isValue (pyTree) -
import Converter.Internal as Internal
import numpy
import KCore.test as test

# Check a float value
node = Internal.createNode('node1', 'DataArray_t', value=1.)
ret = Internal.isValue(node, 1.)
test.testO(ret,1)
#>> True

# Check a int value
node = Internal.createNode('node1', 'DataArray_t', value=1)
ret = Internal.isValue(node, 1)
test.testO(ret,2)
#>> True

# Check a numpy array float values
node = Internal.createNode('node1', 'DataArray_t', value=numpy.zeros(10))
ret = Internal.isValue(node, numpy.zeros(10))
test.testO(ret,3)
#>> True

# Check a string value
node = Internal.createNode('node1', 'DataArray_t', value='toto')
ret = Internal.isValue(node, 'toto')
test.testO(ret,4)
#>> True

# Check a string value
node = Internal.createNode('node1', 'DataArray_t', value='toto')
ret = Internal.isValue(node, 'toto')
test.testO(ret,5)
#>> True
