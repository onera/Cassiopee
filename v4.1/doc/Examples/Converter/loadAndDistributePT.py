# - loadAndDistribute (pyTree) -
import Converter.Filter as Filter
import Converter.Internal as Internal
import Converter.Mpi as Cmpi

# Build case
if Cmpi.rank == 0:
    import Converter.PyTree as C
    import Generator.PyTree as G
    import Transform.PyTree as T
    a = G.cart((0,0,0), (1,1,1), (100,50,50))
    a = T.splitNParts(a, Cmpi.size)
    C.convertPyTree2File(a, 'case1.cgns')
Cmpi.barrier()

h = Filter.Handle('case1.cgns')
a = h.loadAndDistribute()
Internal.printTree(a)
Cmpi.convertPyTree2File(a, 'out.cgns')
