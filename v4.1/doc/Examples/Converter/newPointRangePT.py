# - newPointRange (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newPointRange(name='PointRange', value=[1,10,1,10,5,5]); Internal.printTree(n)
#>> ['PointRange',array(shape=(3, 2),dtype='int32',order='F'),[0 son],'IndexRange_t']

# Attach it to a parent node
d = Internal.newBC('wall', [1,80,30,30,1,2], 'BCWall')
Internal.newPointRange('PointRange', [1,71,29,29,1,5], parent=d)
