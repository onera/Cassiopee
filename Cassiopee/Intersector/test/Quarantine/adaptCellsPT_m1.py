# - adaptCells (pyTree) -
import Intersector.Mpi as XORMPI
import Intersector.PyTree as XOR

import Converter.Internal as I
import Converter.Mpi as Cmpi
import Converter.PyTree as C
import Connector.PyTree as X
import Transform.PyTree as T
import Generator.PyTree as G
import Distributor2.PyTree as D2
import Converter.Distributed as D

import numpy as numpy
import KCore.test as test

LOCAL = test.getLocal()

Nprocs = Cmpi.size
ifname = LOCAL + '/case_' + str(Nprocs) + '.cgns'
ofname = LOCAL + '/out_' + str(Nprocs) + '.cgns'

if Cmpi.rank == 0:
    # Build case
    N = 7
    a = G.cart((0,0,0), (1,1,1), (N,N,N))

    a = T.splitNParts(a, Nprocs)

    a = C.convertArray2NGon(a)

    a = X.connectMatch(a)
    dz = float(N-1)
    a = X.connectMatchPeriodic(a, translation=[0.,0.,dz], tol=1.e-6, dim=3, unitAngle='Degree')

    D2._distribute(a, Nprocs)
    XOR._setZonesAndJoinsUId(a)

    C.convertPyTree2File(a, ifname)

Cmpi.barrier()

t = Cmpi.convertFile2SkeletonTree(ifname)        # read the case as a skeleton
t = Cmpi.readZones(t, ifname, rank=Cmpi.rank)    # load proc zones => ti is a 'LS' : loaded skeleton
procDico = Cmpi.getProcDict(t, reduction=False) # t is a LS => reduction=False : no need for MPI com
zidDico = Cmpi.getPropertyDict(t, 'zid', reduction=False)

Cmpi._convert2PartialTree(t)                    # now t is a partial tree (contains only loaded zones)

# build subdivision requirement
CVmax=3
zs = I.getZones(t)
cell_vals = []
for z in zs:
    n = C.getNCells(I.getZones(z))
    cv = numpy.empty((n,), dtype=I.E_NpyInt)
    cv[:] = 0
    if Cmpi.rank%2 == 0: cv[0] = CVmax
    cell_vals.append(cv)

# add dummy BC and fields
for z in zs:
    C._fillEmptyBCWith(z, 'wall', 'BCWall')
    C._initVars(z, '{centers:Density} = {centers:CoordinateX} + {centers:CoordinateY}')
    C._initVars(z, '{centers:Temperature} = {centers:CoordinateX} + {centers:CoordinateY}')
    C._initVars(z, '{centers:var1} = {centers:CoordinateX} + {centers:CoordinateY}')
    C._initVars(z, '{centers:var2} = {centers:CoordinateX} + {centers:CoordinateY}')
    C._initVars(z, '{centers:var3} = {centers:CoordinateX} + {centers:CoordinateY}')

at = XORMPI.adaptCells(t, cell_vals, sensor_type=3, subdiv_type=0, procDict=procDico, zidDict=zidDico)
at = XORMPI.closeCells(at, procDico, zidDico)

Cmpi.convertPyTree2File(at, ofname)

if Cmpi.rank == 0:
    r = C.convertFile2PyTree(ofname)
    test.testT(r, 1)
