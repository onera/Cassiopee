# - addBase2PyTree (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal

t = C.newPyTree(['Base', 3]) # must contain volume zones
t = C.addBase2PyTree(t, 'Base2', 2) # must contain surface zones
Internal.printTree(t)
#>> ['CGNSTree',None,[3 sons],'CGNSTree_t']
#>>   |_['CGNSLibraryVersion',array([3.1],dtype='float64'),[0 son],'CGNSLibraryVersion_t']
#>>   |_['Base',array(shape=(2,),dtype='int32',order='F'),[0 son],'CGNSBase_t']
#>>   |_['Base2',array(shape=(2,),dtype='int32',order='F'),[0 son],'CGNSBase_t']
