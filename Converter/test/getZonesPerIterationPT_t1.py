# - getZonesPerIterations (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (3,3,3))
t = C.newPyTree(['Base',a])
b = Internal.getNodeFromName(t, 'Base')
n = Internal.newBaseIterativeData(name='BaseIterativeData', nsteps=1, itype='IterationValues', parent=b)
Internal.newDataArray('TimeValues', value=[0.], parent=n)
Internal.newDataArray('NumberOfZones', value=[1], parent=n)
Internal.newDataArray('ZonePointers', value=[['cart']], parent=n)

# List of zones of it=0
zones = Internal.getZonesPerIteration(t, iteration=0)
test.testT(C.newPyTree(['Base',zones]),1)
# List of zones of time 0.
zones = Internal.getZonesPerIteration(t, time=0.)
test.testT(C.newPyTree(['Base',zones]),2)
# All zones of all iterations
zones = Internal.getZonesPerIteration(t)
test.testT(C.newPyTree(['Base',zones]),3)
