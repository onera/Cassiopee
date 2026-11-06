# - loadVariables (pyTree) -
import Converter.Filter as Filter
import Converter.PyTree as C
import Generator.PyTree as G

# Create test case
a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._initVars(a, 'F', 0)
C._initVars(a, 'centers:G', 1)
C.convertPyTree2File(a, 'file.hdf')

# Create a handle on a CGNS file
h = Filter.Handle('file.hdf')

# Load skeleton
a = h.loadSkeleton()

# Load all zones without variable
h._loadZonesWoVars(a)

# Load given variables for all zones
h._loadVariables(a, var=['F', 'centers:G'])

# Load given variables for given zone
h._loadVariables(a, var=['F', 'centers:G'], znp=['Base/cart'])
