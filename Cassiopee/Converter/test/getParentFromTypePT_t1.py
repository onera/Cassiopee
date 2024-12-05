# - getParentFromType (pyTree) -
import Converter.Internal as Internal
import KCore.test as test

a = Internal.createNode('level0', 'DataArray_t', 0)
b = Internal.createChild(a, 'level1', 'DataArray_t', 1)
c = Internal.createChild(b, 'level2', 'DataArray_t', 2)

p = Internal.getParentFromType(a, c, 'DataArray_t')
test.testO(p[0],1)
