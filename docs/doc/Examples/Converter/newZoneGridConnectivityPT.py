# - newZoneGridConnectivity (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newZoneGridConnectivity(name='ZoneGridConnectivity'); Internal.printTree(n)
#>> ['ZoneGridConnectivity',None,[0 son],'ZoneGridConnectivity_t']

# Attach it to a parent node
z = Internal.newZone('Zone', zsize=[[10],[2],[0]], ztype='Structured')
Internal.newZoneGridConnectivity('ZoneGridConnectivity', parent=z)
