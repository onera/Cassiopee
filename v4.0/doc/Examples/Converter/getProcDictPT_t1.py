# - getProcDict (pyTree) -
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Distributor2.PyTree as Distributor2
import Generator.PyTree as G
import KCore.test as test

LOCAL = test.getLocal()

# Cree le fichier test
if Cmpi.rank == 0:
    a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
    b = G.cart( (12,0,0), (1,1,1), (10,10,10) )
    t = C.newPyTree(['Base',a,b])
    C.convertPyTree2File(t, LOCAL+'/in.cgns')
Cmpi.barrier()

# Relit des zones par procs
t = Cmpi.convertFile2SkeletonTree(LOCAL+'/in.cgns')
(t, dic) = Distributor2.distribute(t, NProc=2, algorithm='fast')

# Squelette
d = Cmpi.getProcDict(t)
if Cmpi.rank == 0: test.testO(d, 1)

t = Cmpi.readZones(t, LOCAL+'/in.cgns', rank=Cmpi.rank)
# Squelette charge
d = Cmpi.getProcDict(t)
if Cmpi.rank == 0: test.testO(d, 2)

t = Cmpi.convert2PartialTree(t)
# Partiel charge
d = Cmpi.getProcDict(t)
if Cmpi.rank == 0: test.testO(d, 3)
