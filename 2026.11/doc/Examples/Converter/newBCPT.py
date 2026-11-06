# - newBC (pyTree) -
import Converter.Internal as Internal

# Create a BC node
n = Internal.newBC(name='BC', pointList=[22,1036,101,43], btype='BCFarfield')
Internal.printTree(n)
#>> ['BC',array('BCFarfield',dtype='|S1'),[1 son],'BC_t']
#>>    |_['PointList',array(shape=(4,),dtype='int32',order='F'),[0 son],'IndexArray_t']

# Attach it to a parent node
d = Internal.newZoneBC()
Internal.newBC('BC', [1,45,1,21,1,1], 'BCWall', parent=d)
