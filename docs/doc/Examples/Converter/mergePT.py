# - merge (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G

t1 = C.newPyTree(['Base1', 2, 'Base2', 3])
t2 = C.newPyTree(['Other1', 'Other2', 'Base2'])
a = G.cart((0,0,0), (1,1,1), (10,10,10)); t1[2][1][2] += [a]
t = Internal.merge([t1, t2])

#>> ['CGNSTree',None,[5 sons],'CGNSTree_t']
#>>   |_['CGNSLibraryVersion',array([3.1],dtype='float64'),[0 son],'CGNSLibraryVersion_t']
#>>   |_['Base1',array(shape=(2,),dtype='int32',order='F'),[1 son],'CGNSBase_t']
#>>   |   |_['cart',array(shape=(3, 3),dtype='int32',order='F'),[2 sons],'Zone_t']
#>>   |       ...
#>>   |_['Base2',array(shape=(2,),dtype='int32',order='F'),[0 son],'CGNSBase_t']
#>>   |_['Other1',array(shape=(2,),dtype='int32',order='F'),[0 son],'CGNSBase_t']
#>>   |_['Other2',array(shape=(2,),dtype='int32',order='F'),[0 son],'CGNSBase_t']
C.convertPyTree2File(t, 'out.cgns')
