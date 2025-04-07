# - getParentsFromType (pyTree) -
import Converter.Internal as Internal

a = Internal.createNode('level0', 'DataArray_t', 0)
b = Internal.createChild(a, 'level1', 'DataArray_t', 1)
c = Internal.createChild(b, 'level2', 'DataArray_t', 2)

p = Internal.getParentsFromType(a, c, 'DataArray_t')
for i in p: print(i[0])
#>> level0
#>> level1