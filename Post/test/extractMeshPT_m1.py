# - extractMesh (pyTree) - 
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Post.PyTree as P
import Post.Mpi as Pmpi
import Distributor2.PyTree as Distributor2
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

# Create case
if Cmpi.rank == 0:
    ni = 100; nj = 100; nk = 100
    m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
    m = C.initVars(m, 'Density', 2.)
    m = C.initVars(m, 'centers:cellN', 1)
    m = T.splitNParts(m, 4)
    C.convertPyTree2File(m, 'in.cgns')

    # Extraction mesh
    a = G.cart((0.,0.,0.5), (1., 0.1, 1.), (20, 20, 1)); a[0] = 'extraction'
    a = T.splitNParts(a, 2)
    C.convertPyTree2File(a, 'in2.cgns')
Cmpi.barrier()

# Extract solution on extraction mesh
m = Cmpi.convertFile2SkeletonTree('in.cgns')
(m, dic) = Distributor2.distribute(m, NProc=Cmpi.size, algorithm='fast')
m = Cmpi.readZones(m, 'in.cgns', rank=Cmpi.rank)

a = Cmpi.convertFile2SkeletonTree('in2.cgns')
(a, dic) = Distributor2.distribute(a, NProc=Cmpi.size, algorithm='fast')
a = Cmpi.readZones(a, 'in2.cgns', rank=Cmpi.rank)

a = Pmpi.extractMesh(m, a)
if Cmpi.rank == 0: test.testT(a, 1)

