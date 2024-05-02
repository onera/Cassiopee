# - newDiscreteData (pyTree) -
import Converter.Internal as Internal

# Create a discrete data node
b = Internal.newDiscreteData(name='DiscreteData'); Internal.printTree(b)
#>> ['DiscreteData',None,[0 son],'DiscreteData_t']

# Attach it to zone
z = Internal.newZone('Zone', [[10],[2],[0]], 'Structured')
Internal.newDiscreteData('DiscreteData', parent=z)
