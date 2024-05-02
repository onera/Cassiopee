# - appendBaseName2ZoneName (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart((0,0,0), (1,1,1), (10,10,10)); a[0] = 'a'
b = G.cart((11,0,0), (1,1,1), (10,10,10)); b[0] = 'b'
t = C.newPyTree(['Base',a,b])

t = Internal.appendBaseName2ZoneName(t)

#>> ['CGNSTree',None,[2 sons],'CGNSTree_t']
#>>    |_['CGNSLibraryVersion',array([3.1],dtype='float64'),[0 son],'CGNSLibraryVersion_t']
#>>    |_['Base',array(shape=(2,),dtype='int32',order='F'),[2 sons],'CGNSBase_t']
#>>        |_['Base_a',array(shape=(3, 3),dtype='int32',order='F'),[2 sons],'Zone_t']
#>>        |   |_['ZoneType',array('Structured',dtype='|S1'),[0 son],'ZoneType_t']
#>>        |   ...
#>>        |_['Base_b',array(shape=(3, 3),dtype='int32',order='F'),[2 sons],'Zone_t']
#>>            |_['ZoneType',array('Structured',dtype='|S1'),[0 son],'ZoneType_t']
#>>            ...
C.convertPyTree2File(t, 'out.cgns')
