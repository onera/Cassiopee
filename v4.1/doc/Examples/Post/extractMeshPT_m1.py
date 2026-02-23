# - extractMesh (pyTree) -
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Post.Mpi as Pmpi
import Distributor2.PyTree as Distributor2
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

LOCAL = test.getLocal()

# Create case
if Cmpi.rank == 0:
    ni = 100; nj = 100; nk = 100
    m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
    C._initVars(m, 'Density', 2.)
    C._initVars(m, 'centers:cellN', 1)
    m = T.splitNParts(m, 4)
    C.convertPyTree2File(m, LOCAL+'/in.cgns')

    # Extraction mesh
    a = G.cart((0.,0.,0.5), (1., 0.1, 1.), (20, 20, 1)); a[0] = 'extraction'
    a = T.splitNParts(a, 2)
    C.convertPyTree2File(a, LOCAL+'/in2.cgns')
Cmpi.barrier()

# Extract solution on extraction mesh
m = Cmpi.convertFile2SkeletonTree(LOCAL+'/in.cgns')
dic = Distributor2._distribute(m, NProc=Cmpi.size, algorithm='fast')
m = Cmpi.readZones(m, LOCAL+'/in.cgns', rank=Cmpi.rank)

a = Cmpi.convertFile2SkeletonTree(LOCAL+'/in2.cgns')
dic = Distributor2._distribute(a, NProc=Cmpi.size, algorithm='fast')
a = Cmpi.readZones(a, LOCAL+'/in2.cgns', rank=Cmpi.rank)

a = Pmpi.extractMesh(m, a)
if Cmpi.rank == 0: test.testT(a, 1)
