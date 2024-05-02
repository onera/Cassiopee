# - newZone (pyTree) -
import Converter.Internal as Internal

# Create a zone node
z = Internal.newZone('Zone', zsize=[[10, 2, 0]], ztype='Structured'); Internal.printTree(z)
#>> ['Zone',array(shape=(3, 1),dtype='int32',order='F'),[1 son],'Zone_t']
#>>   |_['ZoneType',array('Structured',dtype='|S1'),[0 son],'ZoneType_t']

# Create a zone node and attach it to tree
t = Internal.newCGNSTree()
b = Internal.newCGNSBase('Base', 3, 3, parent=t)
z = Internal.newZone('Zone', [[10],[2],[0]], 'Structured', parent=b)
