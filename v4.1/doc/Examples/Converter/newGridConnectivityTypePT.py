# - newGridConnectivityType (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newGridConnectivityType(ctype='Abutting1to1'); Internal.printTree(n)
#>> ['GridConnectivityType',array('Abutting1to1',dtype='|S1'),[0 son],'GridConnectivityType_t']

# Attach it to a parent node
d = Internal.newGridConnectivity(name='Match', donorName='blk1', ctype=None)
Internal.newGridConnectivityType(ctype='Abutting1to1', parent=d)
