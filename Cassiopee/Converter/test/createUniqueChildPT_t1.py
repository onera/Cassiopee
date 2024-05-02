# - createUniqueChild (pyTree) -
import Converter.Internal as Internal
import KCore.test as test

node = Internal.createNode('myNode', 'DataArray_t', value=1.)
Internal.createUniqueChild(node, 'childName', 'DataArray_t', value=2.)
# Since childName node already exists. Only the value will be set.
Internal.createUniqueChild(node, 'childName', 'DataArray_t', value=3.)
test.testT(node, 1)
