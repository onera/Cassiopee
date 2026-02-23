# - getProcDict (pyTree) -
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Distributor2.PyTree as Distributor2
import Generator.PyTree as G

# Cree le fichier test
if Cmpi.rank == 0:
    a = G.cart((0,0,0), (1,1,1), (10,10,10))
    b = G.cart((12,0,0), (1,1,1), (10,10,10))
    t = C.newPyTree(['Base',a,b])
    C.convertPyTree2File(t, 'in.cgns')
Cmpi.barrier()

# Relit des zones par procs
t = Cmpi.convertFile2SkeletonTree('in.cgns')
(t, dic) = Distributor2.distribute(t, NProc=2, algorithm='fast')

procDict = Cmpi.getProcDict(t)
if Cmpi.rank == 0: print(procDict)
#>> {'cart.0': 1, 'cart': 0}
