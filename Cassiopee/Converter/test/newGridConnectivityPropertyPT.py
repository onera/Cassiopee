# - newGridConnectivityProperty (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newGridConnectivityProperty(); Internal.printTree(n)
#>> ['GridConnectivityProperty',None,[0 son],'GridConnectivityProperty_t']

# Attach it to a parent node
d = Internal.newGridConnectivity(name='Match', donorName='blk1', ctype='Abutting1to1')
Internal.newGridConnectivityProperty(parent=d)
