import Fast.PyTree as Fast
import FastS.PyTree as FastS
import Apps.Fast.FastIBM as FastIBM
import Converter.Internal as Internal
import Converter.PyTree as C
import Transform.PyTree as T
import KCore.test as test

import numpy
LOCAL = test.getLocal()

##PREP
def multiElementAirfoil(snear=0.001, ibctype='Musker', alpha=16.):
    """Generate an IBM case for the canonical subsonic 2DMEA test-case."""
    t = C.convertFile2PyTree('2DMEA.stp', hmax=10*snear)
    Internal._rmNodesByName(t, 'CAD')
    Internal._rmNodesByName(t, 'FACES')
    Internal._rmNodesByType(t, Internal.__FlowSolutionNodes__)

    t = T.scale(t, 1./25.4, (0.,0.,0.))
    # t = T.translate(t, (0., -0.98965, 0.))
    t = C.convertArray2Tetra(t)
    t = T.join(t)
    t = T.splitConnexity(t)
    zones = sorted([z for z in Internal.getZones(t)], key=lambda z: numpy.min(Internal.getNodeFromName(z, 'CoordinateX')[1]))
    zones[0][0] = 'slat'
    zones[1][0] = 'main'
    zones[2][0] = 'flap'

    t = C.newPyTree(['Base', zones])

    FastIBM._setSnear(t, snear)
    FastIBM._setIBCType(t, ibctype)
    FastIBM._setDfar(t, 100)

    C._addState(t, adim='adim1', MInf=0.2, alphaZ=alpha, alphaY=0., ReInf=5.e6,\
                MutSMuInf=0.2, TurbLevelInf=0.0001, EquationDimension=2, GoverningEquations='NSTurbulent')

    return t

tb = multiElementAirfoil(snear=0.005)
t,tc = FastIBM.prepareIBMData(tb, None, None, vmin=21, expand=3, frontType=1)

test.testT(t , 1)
test.testT(tc, 2)
test.testT(tb, 3)

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
Fast._setNum2Base(t, numb); Fast._setNum2Zones(t, numz)

t, tc, metrics = FastS.warmup(t, tc)

for it in range(200):
    if it%100 == 0: print("it %d / 200"%it, flush=True)
    FastS._compute(t, metrics, it, tc)

Internal._rmNodesFromName(t, 'Parameter_int')
Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
Internal._rmNodesFromName(tc, 'Parameter_int')
Internal._rmNodesByName(tc, '.Solver#Param')
Internal._rmNodesByName(tc, '.Solver#ownData')

test.testT(t , 4)
test.testT(tc, 5)

##POST
graphIBCDPost, ts = FastIBM.prepareSkinReconstruction(tb, tc, dimPb=2, ibctypes=[3])
FastIBM._computeSkinVariables(ts, tc, graphIBCDPost, ibctypes=[3], dimPb=2)
C._rmVars(ts, ['yplus', 'yplusIP'])

test.testT(ts, 6)