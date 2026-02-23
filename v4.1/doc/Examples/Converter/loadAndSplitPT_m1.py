# - loadAndSplit (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Filter as Filter
import Converter.Mpi as Cmpi
import KCore.test as test

LOCAL = test.getLocal()

# Build case
if Cmpi.rank == 0:
    a = G.cart((0,0,0), (1,1,1), (100,100,100))
    C.convertPyTree2File(a, LOCAL+'/out.hdf')
Cmpi.barrier()

# With skeleton
h = Filter.Handle(LOCAL+'/out.hdf')
t = h.loadAndSplitSkeleton(NParts=3)
h._loadContainerPartial(t, variablesN=['GridCoordinates/CoordinateX','GridCoordinates/CoordinateY','GridCoordinates/CoordinateZ'])
Cmpi.convertPyTree2File(t, LOCAL+'/out.cgns')

# In one go
h = Filter.Handle(LOCAL+'/out.hdf')
t = h.loadAndSplit(NParts=3)
Cmpi.convertPyTree2File(t, LOCAL+'/out.cgns')
