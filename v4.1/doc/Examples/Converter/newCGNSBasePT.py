# - newCGNSBase (pyTree) -
import Converter.Internal as Internal

# Create a base node
b = Internal.newCGNSBase('Base'); print(b)
#>> ['Base', array([3, 3], dtype=int32), [], 'CGNSBase_t']

# Create a base node and attach it to tree
t = Internal.newCGNSTree()
Internal.newCGNSBase('Base', 3, 3, parent=t); Internal.printTree(t)
#>> ['CGNSTree',None,[2 sons],'CGNSTree_t']
#>>   |_['CGNSLibraryVersion',array([3.1],dtype='float64'),[0 son],'CGNSLibraryVersion_t']
#>>   |_['Base',array(shape=(2,),dtype='int32',order='F'),[0 son],'CGNSBase_t']
