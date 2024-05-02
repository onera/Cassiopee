# - readZones (pyTree) -
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Distributor2.PyTree as Distributor2
import Generator.PyTree as G

# Cree le fichier test
if Cmpi.rank == 0:
    a = G.cart((0,0,0), (1,1,1), (10,10,10))
    b = G.cart((12,0,0), (1,1,1), (10,10,10))
    t = C.newPyTree(['Base',a,b])
    C.convertPyTree2File(t, 'test.cgns')
Cmpi.barrier()

# Relit les zones par paths
t = Cmpi.convertFile2SkeletonTree('test.cgns')
Cmpi._readZones(t, 'test.cgns', zoneNames=['Base/cart'])

# Relit des zones par procs
t = Cmpi.convertFile2SkeletonTree('test.cgns')
(t, dic) = Distributor2.distribute(t, NProc=Cmpi.size, algorithm='fast')
t = Cmpi.readZones(t, 'test.cgns', rank=Cmpi.rank)
if Cmpi.rank == 0: C.convertPyTree2File(t, 'out.cgns')
