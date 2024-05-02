# - createChild (pyTree) -
import Converter.Internal as Internal
import KCore.test as test

node = Internal.createNode('myNode', 'DataArray_t', value=1., children=[])
Internal.createChild(node, 'childName', 'DataArray_t', value=2., children=[])
test.testT(node, 1)
