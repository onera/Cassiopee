# - loadZonesVars (pyTree) -
import Converter.Filter as Filter
import Converter.PyTree as C
import Generator.PyTree as G

# Create test case
a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = G.cart((10,0,0), (1,1,1), (10,10,10))
C.convertPyTree2File([a,b], 'file.hdf')

# Create a handle on a CGNS file
h = Filter.Handle('file.hdf')

# Load skeleton
a = h.loadSkeleton()

# Fully load all zones
h._loadZones(a)

# Fully load two zones
h._loadZones(a, znp=['Base/cart', 'Base/cart.0'])
