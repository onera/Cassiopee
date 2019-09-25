# - loadAndSplit (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Filter as Filter

# Case
a = G.cart((0,0,0), (1,1,1), (10,10,10))
C.convertPyTree2File(a, 'out.hdf')

# With skeleton
h = Filter.Handle('out.hdf')
t = h.loadAndSplitSkeleton(NParts=3)
h._loadContainerPartial(t, variablesN=['GridCoordinates/CoordinateX','GridCoordinates/CoordinateY','GridCoordinates/CoordinateZ'])
C.convertPyTree2File(t, 'out.cgns')

# In one go
h = Filter.Handle('out.hdf')
t = h.loadAndSplit(NParts=3)
C.convertPyTree2File(t, 'out.cgns')
