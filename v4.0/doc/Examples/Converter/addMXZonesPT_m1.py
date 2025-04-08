# - addMXZones (pyTree) -
import Converter.PyTree as CP
import Generator.PyTree as GP
import Transform.PyTree as TP
import Converter.Filter as Filter
import Converter.Mpi    as Cmpi
import KCore.test       as test

LOCAL = test.getLocal()

# Create case
if Cmpi.rank == 0:
    ni = 50 ; nj = 50 ; nk = 50
    a = GP.cart(( 0,0,0), (10./(ni-1),10./(nj-1),10./(nk-1)), (ni,nj,nk))
    b = GP.cart((10,0,0), (10./(ni-1),10./(nj-1),10./(nk-1)), (ni,nj,nk))
    b = TP.reorder(b,(-3,2,1))
    t = CP.newPyTree(['Base',a,b])
    t = TP.splitNParts(t, 5, multigrid=0, dirs=[1,2,3])
    CP.convertPyTree2File(t, LOCAL+'/case.cgns')
Cmpi.barrier()

# Load
h = Filter.Handle(LOCAL+'/case.cgns')
a = h.loadAndDistribute()

# # Ajout XZones
Cmpi._addMXZones(a)

# Test
if Cmpi.rank == 0: test.testT(a, 1)
