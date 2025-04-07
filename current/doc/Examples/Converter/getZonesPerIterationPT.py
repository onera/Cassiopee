# - getZonesPerIterations (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart((0,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a])
b = Internal.getNodeFromName(t, 'Base')
n = Internal.newBaseIterativeData(name='BaseIterativeData', nsteps=1, itype='IterationValues', parent=b)
Internal.newDataArray('TimeValues', value=[0.], parent=n)
Internal.newDataArray('NumberOfZones', value=[1], parent=n)
Internal.newDataArray('ZonePointers', value=[['cart']], parent=n)

# List of zones of it=0
zones = Internal.getZonesPerIteration(t, iteration=0)
# List of zones of time 0.
zones = Internal.getZonesPerIteration(t, time=0.)
# All zones of all iterations
zones = Internal.getZonesPerIteration(t)
