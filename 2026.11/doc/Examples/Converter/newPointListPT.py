# - newPointList (pyTree) -
import Converter.Internal as Internal

# Create a Descriptor node
n = Internal.newPointList(name='PointList', value=[101,51,22,1036,2]); Internal.printTree(n)
#>> ['PointList',array(shape=(5,),dtype='int32',order='F'),[0 son],'IndexArray_t']

# Attach it to a parent node
d = Internal.newBC('wall', [1,80,30,30,1,2], 'BCWall')
Internal.newPointList('PointList', [101,51,22,1036,2], parent=d)
