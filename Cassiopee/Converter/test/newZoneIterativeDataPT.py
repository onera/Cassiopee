# - newZoneIterativeData (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newZoneIterativeData(name='ZoneIterativeData'); Internal.printTree(n)
#>> ['ZoneIterativeData',None,[0 son],'ZoneIterativeData_t']

# Attach it to a parent node
z = Internal.newZone('Zone', zsize=[[10],[2],[0]], ztype='Structured')
Internal.newZoneIterativeData(name='ZoneIterativeData', parent=z)
