# - addBXZones (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.Filter as Filter
import Converter.Mpi    as Cmpi
import KCore.test       as test

LOCAL = test.getLocal()

# Create case
if Cmpi.rank == 0:
    ni = 50 ; nj = 50 ; nk = 50
    a = G.cart(( 0,0,0), (10./(ni-1),10./(nj-1),10./(nk-1)), (ni,nj,nk))
    b = G.cart((10,0,0), (10./(ni-1),10./(nj-1),10./(nk-1)), (ni,nj,nk))

    b = T.reorder(b,(-3,2,1))

    t = C.newPyTree(['Base',a,b])
    t = T.splitNParts(t, 5, multigrid=0, dirs=[1,2,3])

    C.convertPyTree2File(t, LOCAL+'/case.cgns')

Cmpi.barrier()

# Load
h = Filter.Handle(LOCAL+'/case.cgns')
a = h.loadAndDistribute()

# Ajout XZones
Cmpi._addBXZones(a)

# Test
if Cmpi.rank == 0: test.testT(a, 1)
