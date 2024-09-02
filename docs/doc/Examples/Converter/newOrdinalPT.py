# - newOrdinal (pyTree) -
import Converter.Internal as Internal

# Create a ordinal node
b = Internal.newOrdinal(value=1); Internal.printTree(b)
#>> ['Ordinal',array([1],dtype='int32'),[0 son],'Ordinal_t']

# Attach it to zone
z = Internal.newZone('Zone', [[10],[2],[0]], 'Structured')
Internal.newOrdinal(value=1, parent=z)
