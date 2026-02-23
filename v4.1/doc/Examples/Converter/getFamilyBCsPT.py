# - getFamilyBCs (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart((0.,0.,0), (0.01,0.01,1.), (20,20,2))
b = G.cart((1.,0.,0), (0.01,0.01,1.), (20,20,2))

a = C.addBC2Zone(a, 'walla', 'FamilySpecified:CARTER', 'imin')
b = C.addBC2Zone(b, 'wallb', 'FamilySpecified:CARTER', 'jmin')

t = C.newPyTree(['Base',a,b])
C._addFamily2Base(t[2][1], 'CARTER', bndType='BCWall')
B1 = C.getFamilyBCs(t, 'CARTER'); Internal.printTree(B1)
#>> ['walla',array('FamilySpecified',dtype='|S1'),[2 sons],'BC_t']
#>>   |_['PointRange',array(shape=(3, 2),dtype='int32',order='F'),[0 son],'IndexRange_t']
#>>   |_['FamilyName',array('CARTER',dtype='|S1'),[0 son],'FamilyName_t']
#>> ['wallb',array('FamilySpecified',dtype='|S1'),[2 sons],'BC_t']
#>>   |_['PointRange',array(shape=(3, 2),dtype='int32',order='F'),[0 son],'IndexRange_t']
#>>   |_['FamilyName',array('CARTER',dtype='|S1'),[0 son],'FamilyName_t']
