# - addProcNode (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import Distributor2.PyTree as D2

a = G.cart((0,0,0), (1,1,1), (10,10,10))
D2._addProcNode(a, 12)
Internal.printTree(a)
#>> ['cart',array(shape=(3, 3),dtype='int32',order='F'),[3 sons],'Zone_t']
#>>   |_['ZoneType',array('b'Structured'',dtype='|S1'),[0 son],'ZoneType_t']
#>>   | ...
#>>   |_['.Solver#Param',None,[1 son],'UserDefinedData_t']
#>>       |_['proc',array([12],dtype='int32'),[0 son],'DataArray_t']
C.convertPyTree2File(a, 'out.cgns')
