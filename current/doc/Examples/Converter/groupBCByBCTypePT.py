# - groupBCByType (PyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G

a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(10,10,10))
a = C.fillEmptyBCWith(a, 'wall', 'BCWall')

t = C.newPyTree(['Base',a])
Internal._groupBCByBCType(t, 'BCWall', 'FamWall')

#>> ['CGNSTree',None,[2 sons],'CGNSTree_t']
#>>   |_['CGNSLibraryVersion',array([3.1],dtype='float64'),[0 son],'CGNSLibraryVersion_t']
#>>   |_['Base',array(shape=(2,),dtype='int32',order='F'),[2 sons],'CGNSBase_t']
#>>       |_['cart',array(shape=(3, 3),dtype='int32',order='F'),[3 sons],'Zone_t']
#>>       |   ...
#>>       |   |_['ZoneBC',None,[6 sons],'ZoneBC_t']
#>>       |       |_['wall1',array('FamilySpecified',dtype='|S1'),[2 sons],'BC_t']
#>>       |       |   |_['PointRange',array(shape=(3, 2),dtype='int32',order='F'),[0 son],'IndexRange_t']
#>>       |       |   |_['FamilyName',array('FamWall',dtype='|S1'),[0 son],'FamilyName_t']
#>>       |       |_['wall2',array('FamilySpecified',dtype='|S1'),[2 sons],'BC_t']
#>>       |       |   |_['PointRange',array(shape=(3, 2),dtype='int32',order='F'),[0 son],'IndexRange_t']
#>>       |  ...
#>>       |_['FamWall',None,[1 son],'Family_t']
#>>           |_['FamilyBC',array('BCWall',dtype='|S1'),[0 son],'FamilyBC_t']
C.convertPyTree2File(t, 'out.cgns')
