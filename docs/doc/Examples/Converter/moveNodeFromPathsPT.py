# - moveNodeFromPaths (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal

a = G.cart((0,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base1', a, 'Base2'])
Internal.printTree(t)
#>> ['CGNSTree',None,[3 sons],'CGNSTree_t']
#>>   |_['Base1',array(shape=(2,),dtype='int32',order='F'),[1 son],'CGNSBase_t']
#>>   |   |_['cart',array(shape=(3, 3),dtype='int32',order='F'),[2 sons],'Zone_t']
#>>   |_['Base2',array(shape=(2,),dtype='int32',order='F'),[0 son],'CGNSBase_t']
Internal._moveNodeFromPaths(t, 'Base1/cart', 'Base2')
Internal.printTree(t)
#>> ['CGNSTree',None,[3 sons],'CGNSTree_t']
#>>   |_['Base1',array(shape=(2,),dtype='int32',order='F'),[0 son],'CGNSBase_t']
#>>   |_['Base2',array(shape=(2,),dtype='int32',order='F'),[1 son],'CGNSBase_t']
#>>       |_['cart',array(shape=(3, 3),dtype='int32',order='F'),[2 sons],'Zone_t']
