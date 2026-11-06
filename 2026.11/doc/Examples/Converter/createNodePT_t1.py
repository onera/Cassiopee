# - createNode (pyTree) -
import Converter.Internal as Internal
import numpy
import KCore.test as test

# Create a node named myNode, of type DataArray_t, with one value 1.
node1 = Internal.createNode('myNode', 'DataArray_t', value=1., children=[])

# Create a node named myNode, of type DataArray_t, with an array valued in a list
node2 = Internal.createNode('myNode', 'DataArray_t', value=[12.,14.,15.], children=[])

# Create a node named myNode, of type DataArray_t, with a given numpy
a = numpy.zeros( (10) )
a[1] = 1.; a[8] = 2.
node3 = Internal.createNode('myNode', 'DataArray_t', value=a, children=[])

# Create a node named GoverningEquation, of type 'GoverningEquation_t' and value 'Euler'
node4 = Internal.createNode('GoverningEquation', 'GoverningEquation_t', value='Euler')

node = Internal.createNode('Father', 'None_t', value=None, children=[node1,node2,node3,node4])
test.testT(node, 1)
