# - loadSkeleton (pyTree) -
import Converter.Filter as Filter
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# Create test case
a = G.cart((0,0,0), (1,1,1), (10,10,10))
C.convertPyTree2File(a, 'file.hdf')

# Create a handle on a CGNS file
h = Filter.Handle('file.hdf')

# Load skeleton
a = h.loadSkeleton()

test.testT(a, 1)
