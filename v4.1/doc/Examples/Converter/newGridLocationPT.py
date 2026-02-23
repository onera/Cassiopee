# - newGridLocation (pyTree) -
import Converter.Internal as Internal

# Create a GridLocation node
n = Internal.newGridLocation(value='CellCenter'); Internal.printTree(n)
#>> ['GridLocation',array('CellCenter',dtype='|S1'),[0 son],'GridLocation_t']

# Attach it to a parent node
d = Internal.newBC('wall', [1,80,30,30,1,2], 'BCWall')
Internal.newGridLocation('Vertex', parent=d)
