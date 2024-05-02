# - newBCProperty (pyTree) -
import Converter.Internal as Internal

# Create a BC property node
n = Internal.newBCProperty(wallFunction='Null', area='Null'); Internal.printTree(n)
#>> ['BCProperty',None,[2 sons],'BCProperty_t']
#>>    |_['WallFunctionType',array('Null',dtype='|S1'),[0 son],'WallFunctionType_t']
#>>    |_['Area',array('Null',dtype='|S1'),[0 son],'Area_t']

# Attach it to a parent node
d = Internal.newBC(name='BC', pointList=[22,1036,101,43], btype='BCWall')
Internal.newBCProperty(wallFunction='Null', area='Null', parent=d)
