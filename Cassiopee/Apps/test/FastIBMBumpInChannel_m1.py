import FastC.PyTree as FastC
import FastS.Mpi as FastS
import Apps.Fast.FastIBM as FastIBM
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import KCore.test as test

##PREP
tb = FastIBM.bumpInChannel(snear=0.005)
t,tc = FastIBM.prepareIBMData(tb, None, None, vmin=21, expand=3, frontType=1)

if Cmpi.rank==0:
    test.testT(t , 1)
    test.testT(tc, 2)

##COMPUTE
numb = {}
numb["temporal_scheme"]    = "implicit"
numb["ss_iteration"]       = 1
numb["omp_mode"]           = 0
numb["modulo_verif"]       = 100

numz = {}
numz["time_step"]          = 3.e-5
numz["time_step_nature"]   = "local"
numz["cfl"]                = 5.
numz["scheme"]             = "roe"

it0 = 0.; time0 = 0.
FastC._setNum2Base(t, numb); FastC._setNum2Zones(t, numz)

tcskel = Cmpi.convert2SkeletonTree(tc)
tcskel = Cmpi.allgatherTree(tcskel)
graph = FastC.prepGraphs(tcskel)
del tcskel

t, tc, metrics = FastS.warmup(t, tc, graph=graph)

for it in range(200):
    if Cmpi.rank == 0 and it%100 == 0: print("it %d / 200"%it, flush=True)
    FastS._compute(t, metrics, it, tc, graph=graph)

Internal._rmNodesFromName(t, 'Parameter_int')
Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
Internal._rmNodesFromName(tc, 'Parameter_int')
Internal._rmNodesByName(tc, '.Solver#Param')
Internal._rmNodesByName(tc, '.Solver#ownData')

if Cmpi.rank==0:
    test.testT(t , 3)
    test.testT(tc, 4)

##POST
graphIBCDPost, ts = FastIBM.prepareSkinReconstruction(tb, tc, dimPb=2, ibctypes=[3])
FastIBM._computeSkinVariables(ts, tc, graphIBCDPost, ibctypes=[3], dimPb=2)

wall, aeroLoads = FastIBM.computeAerodynamicLoads(ts, dimPb=2, famZones=[], Pref=None, verbose=0)
wall, aeroLoads = FastIBM.computeAerodynamicCoefficients(wall, aeroLoads, dimPb=2, Sref=1., Lref=1., Qref=None, alpha=0., beta=0., verbose=0)

import Converter.PyTree as C
C._rmVars(wall, ['yplus', 'yplusIP'])

if Cmpi.rank==0:
    test.testT(wall, 5)