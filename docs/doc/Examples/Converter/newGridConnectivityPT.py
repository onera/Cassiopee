# - newGridConnectivity (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newGridConnectivity(name='Match', donorName='blk1', ctype='Abutting1to1'); Internal.printTree(n)
#>> ['Match',array('blk1',dtype='|S1'),[1 son],'GridConnectivity_t']
#>>    |_['GridConnectivityType',array('Abutting1to1',dtype='|S1'),[0 son],'GridConnectivityType_t']

# Attach it to a parent node
d = Internal.newZoneGridConnectivity(name='ZoneGridConnectivity')
Internal.newGridConnectivity(name='Match', donorName='blk1', ctype='Abutting1to1', parent=d)
