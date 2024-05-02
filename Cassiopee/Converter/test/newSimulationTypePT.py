# - newSimulationType (pyTree) -
import Converter.Internal as Internal

# Create a simulation type node
b = Internal.newSimulationType(value='NonTimeAccurate'); Internal.printTree(b)
#>> ['SimulationType',array('NonTimeAccurate',dtype='|S1'),[0 son],'SimulationType_t']

# Attach it to parent
b = Internal.newCGNSBase('Base', 3, 3)
Internal.newSimulationType('TimeAccurate', parent=b)
