# - writeNodesFromPaths (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Filter as Filter

t = C.newPyTree(['Base'])
C.convertPyTree2File(t, 'out.hdf')

# Append
a = G.cart((0,0,0), (1,1,1), (10,10,10))
Filter.writeNodesFromPaths('out.hdf', 'CGNSTree/Base', a)

# Append and replace
a = G.cart((1,1,1), (1,1,1), (10,10,10)); a[0] = 'cart'
Filter.writeNodesFromPaths('out.hdf', 'CGNSTree/Base', a)
