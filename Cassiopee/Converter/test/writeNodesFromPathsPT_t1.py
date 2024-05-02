# - writeNodesFromPaths (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Filter as Filter
import KCore.test as test

LOCAL = test.getLocal()

t = C.newPyTree(['Base'])
C.convertPyTree2File(t, LOCAL+'/out.adf')

# Append
a = G.cart((0,0,0), (1,1,1), (10,10,10))
Filter.writeNodesFromPaths(LOCAL+'/out.adf', 'CGNSTree/Base', a)

# Append and Replace
a = G.cart((1,1,1), (1,1,1), (10,10,10)); a[0] = 'cart'
Filter.writeNodesFromPaths(LOCAL+'/out.adf', 'CGNSTree/Base', a)

# Replace
a = G.cart((1,1,1), (0.1,1,1), (10,10,10)); a[0] = 'cart'
Filter.writeNodesFromPaths(LOCAL+'/out.adf', 'CGNSTree/Base/cart', a, mode=1)

#test.testF('out.adf', 1)


t = C.newPyTree(['Base'])
C.convertPyTree2File(t, LOCAL+'/out.hdf')

# Append
a = G.cart((0,0,0), (1,1,1), (10,10,10))
Filter.writeNodesFromPaths(LOCAL+'/out.hdf', 'CGNSTree/Base', a)

# Append and Replace
a = G.cart((1,1,1), (1,1,1), (10,10,10)); a[0] = 'cart'
Filter.writeNodesFromPaths(LOCAL+'/out.hdf', 'CGNSTree/Base', a)

# Replace
a = G.cart((1,1,1), (0.1,1,1), (10,10,10)); a[0] = 'cart'
Filter.writeNodesFromPaths(LOCAL+'/out.hdf', 'CGNSTree/Base/cart', a, mode=1)

#test.testF('out.hdf', 2)
