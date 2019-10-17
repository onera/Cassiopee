import Converter.Filter as Filter
import Converter.Internal as Internal
import Converter.Mpi as Cmpi

h = Filter.Handle('case1.cgns')
a = h.loadFromProc()
Internal.printTree(a)
Cmpi.convertPyTree2File(a, 'out.cgns')
