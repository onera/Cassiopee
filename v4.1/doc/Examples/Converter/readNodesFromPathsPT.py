# - readNodesFromPaths (pyTree) -
import Converter.PyTree as C
import Converter.Filter as Filter
import Converter.Internal as Internal
import Generator.PyTree as G

# Cree le fichier test
a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = G.cart((12,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a,b])
C.convertPyTree2File(t, 'test.hdf')

# Relit les noeuds par leur paths
nodes = Filter.readNodesFromPaths('test.hdf', ['CGNSTree/Base/cart/GridCoordinates'])
Internal.printTree(nodes)
#>> ['GridCoordinates',None,[3 sons],'GridCoordinates_t']
#>>   |_['CoordinateX',array(shape=(10, 10, 10),dtype='float64',order='F'),[0 son],'DataArray_t']
#>>   |_['CoordinateY',array(shape=(10, 10, 10),dtype='float64',order='F'),[0 son],'DataArray_t']
#>>   |_['CoordinateZ',array(shape=(10, 10, 10),dtype='float64',order='F'),[0 son],'DataArray_t']
