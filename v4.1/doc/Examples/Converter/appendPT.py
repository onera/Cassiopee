# - append (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

t = C.newPyTree(['Base', 'Base2'])
a = G.cart((0,0,0), (1,1,1), (10,10,10))

# Append a to 'Base' node of t
t = Internal.append(t, a, 'Base')

#>> ['CGNSTree',None,[3 sons],'CGNSTree_t']
#>>    |_['CGNSLibraryVersion',array([3.1],dtype='float64'),[0 son],'CGNSLibraryVersion_t']
#>>    |_['Base',array(shape=(2,),dtype='int32',order='F'),[1 son],'CGNSBase_t']
#>>    |   |_['cart',array(shape=(3, 3),dtype='int32',order='F'),[2 sons],'Zone_t']
#>>    |       ...
#>>    |_['Base2',array(shape=(2,),dtype='int32',order='F'),[0 son],'CGNSBase_t']
C.convertPyTree2File(t, 'out.cgns')
