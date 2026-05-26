# - newGridCoordinates (pyTree) -
import Converter.Internal as Internal

# Create a GridCoordinates node
n = Internal.newGridCoordinates(); print(n)
#>> ['GridCoordinates', None, [], 'GridCoordinates_t']

# Create a zone node
z = Internal.newZone('Zone', zsize=[[10,9,0]], ztype='Unstructured')
n = Internal.newGridCoordinates(parent=z); Internal.printTree(z)
#>> ['Zone',array(shape=(1, 3),dtype='int32',order='F'),[2 sons],'Zone_t']
#>>   |_['ZoneType',array('Structured',dtype='|S1'),[0 son],'ZoneType_t']
#>>   |_['GridCoordinates',None,[0 son],'GridCoordinates_t']
