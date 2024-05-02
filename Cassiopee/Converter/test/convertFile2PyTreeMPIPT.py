# - convertFile2PyTree (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Mpi as Cmpi
import Converter.Internal as Internal

if Cmpi.rank == 0:
    a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(11,11,11))
    t = C.newPyTree(['Base',a])
    C.convertPyTree2File(t, 'in.cgns')
Cmpi.barrier()
# Identique sur tous les procs
t1 = Cmpi.convertFile2PyTree('in.cgns')
if Cmpi.rank == 1 or Cmpi.rank == 0:
    Internal.printTree(t1)
