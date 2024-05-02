# - newBCDataSet (pyTree) -
import Converter.Internal as Internal

# Create a BC data set node
n = Internal.newBCDataSet(name='BCDataSet', value='BCWall'); Internal.printTree(n)
#>> ['BCDataSet',array('BCWall',dtype='|S1'),[0 son],'BCDataSet_t']

# Attach it to a parent node
d = Internal.newBC(name='BC', pointList=[22,1036,101,43], btype='BCFarfield')
Internal.newBCDataSet(name='BCDataSet', value='UserDefined', parent=d)

# Complete BC + BCDataSet + BCData
b = Internal.newBC(name='BC', pointList=[22,1036,101,43], btype='BCFarfield')
d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined', gridLocation='FaceCenter', parent=b)
d = Internal.newBCData('BCNeumann', parent=d)
d = Internal.newDataArray('Density', value=[1.,1.,1.,1.], parent=d)
