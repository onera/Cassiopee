# - connectMatchPeriodic 3D MPI (pyTree)-
import Generator.PyTree    as G
import Converter.PyTree    as C
import Transform.PyTree    as T
import Converter.Mpi       as Cmpi
import Connector.Mpi       as Xmpi
import Distributor2.PyTree as Distributor2
import Converter.Filter    as Filter
import KCore.test          as test
LOCAL = test.getLocal()

# Cree le fichier test
if Cmpi.rank == 0:
    a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 90., 5., (11,11,11))
    C._initVars(a,'F',1.); C._initVars(a,'centers:G',2.)
    C._addBC2Zone(a,'wall','BCWall','jmin')
    C._addBC2Zone(a,'overlap','BCOverlap','jmax')
    t = C.newPyTree(['Base']); t[2][1][2] += [a]
    C._addState(t[2][1], 'EquationDimension', 3)
    t = T.splitNParts(t, 4, multigrid=0, dirs=[1,2,3])
    C.convertPyTree2File(t, LOCAL+'/in.cgns')

Cmpi.barrier()

h  = Filter.Handle(LOCAL+'/in.cgns')
a  = h.loadAndDistribute()
a  = Xmpi.connectMatchPeriodic(a,rotationCenter=[0.,0.,0.],rotationAngle=[0.,0.,90.])

if Cmpi.rank == 0:
    test.testT(t,1)
