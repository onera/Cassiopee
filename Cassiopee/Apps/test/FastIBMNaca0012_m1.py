import FastC.PyTree as FastC
import FastS.Mpi as FastS
import Apps.Fast.FastIBM as FastIBM
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import KCore.test as test

##PREP
tb = FastIBM.naca0012(snear=0.005, alpha=0.)
t,tc = FastIBM.prepareIBMData(tb, None, None, vmin=21, expand=3, frontType=1)

if Cmpi.rank==0:
    ####
    # The following lines are to avoid regression since the bug fix for duplicate information in tc
    ####
    dictOfDoublons = {
        'CartX0'    :  ['ID_Cart.28X1', 'ID_CartX1'],
        'Cart.1X0'  :['ID_Cart.2X0', 'ID_Cart.15X0', 'ID_Cart.4X0'],
        'Cart.2X0'  :['ID_Cart.11X0', 'ID_Cart.10X0', 'ID_Cart.1X0', 'ID_Cart.14X0'],
        'Cart.4X0'  :['ID_Cart.15X0', 'ID_Cart.1X0', 'ID_Cart.9X0', 'IBCD_3_Cart.4X0', 'ID_Cart.11X0'],
        'Cart.8X0'  :['ID_Cart.8X1', 'ID_Cart.27X1', 'ID_Cart.15X0'],
        'Cart.11X0' :['ID_Cart.2X0', 'ID_Cart.10X0', 'ID_Cart.4X0'],
        'Cart.19X0' :['ID_Cart.25X1', 'ID_Cart.6X0'],
        'Cart.13X1' :['ID_Cart.10X0', 'ID_Cart.14X0', 'ID_Cart.12X1'],
        'Cart.22X1' :['ID_CartX1'],
        'Cart.26X1' :['ID_Cart.3X1', 'ID_Cart.29X1', 'ID_Cart.27X1'],
        'Cart.28X1' :['ID_CartX0', 'ID_Cart.5X0', 'ID_CartX1', 'ID_Cart.6X0'],
        'Cart.8X1'  :['ID_Cart.15X0', 'ID_Cart.8X0'],
        'Cart.11X1' :['ID_Cart.20X1', 'ID_Cart.31X1', 'ID_Cart.12X1'],
        'Cart.14X0' :['ID_Cart.2X0', 'ID_Cart.10X0', 'ID_Cart.13X1'],
        'Cart.15X0' :['ID_Cart.4X0', 'ID_Cart.8X0', 'ID_Cart.1X0', 'ID_Cart.8X1'],
        'Cart.25X1' :['ID_Cart.2X1', 'ID_Cart.29X1', 'ID_Cart.3X1', 'ID_Cart.19X0'],
        'Cart.5X0'  :['ID_Cart.6X0', 'ID_Cart.28X1'],
        'Cart.0X1'  :['ID_Cart.24X1', 'ID_Cart.1X1', 'ID_Cart.21X1'],
        'Cart.12X0' :['ID_Cart.15X1'],
        'Cart.14X1' :['ID_Cart.17X0', 'ID_Cart.30X1'],
        'Cart.15X1' :['ID_Cart.9X1', 'ID_Cart.35X1', 'ID_Cart.12X0', 'ID_Cart.33X1'],
        'Cart.17X0' :['ID_Cart.30X1', 'ID_Cart.14X1', 'ID_Cart.13X0', 'ID_Cart.18X1', 'ID_Cart.16X0'],
        'Cart.20X0' :['ID_Cart.23X1', 'ID_Cart.35X1', 'ID_Cart.17X1'],
        'Cart.21X1' :['ID_Cart.9X1', 'ID_Cart.0X1', 'ID_Cart.10X1', 'ID_Cart.1X1'],
        'Cart.30X1' :['ID_Cart.14X1', 'ID_Cart.17X0', 'ID_Cart.18X1', 'ID_Cart.13X0'],
        'Cart.33X1' :['ID_Cart.9X1', 'ID_Cart.35X1', 'ID_Cart.17X1', 'ID_Cart.15X1', 'ID_Cart.32X1'],
        'Cart.34X1' :['ID_Cart.7X0', 'ID_Cart.18X0', 'ID_Cart.10X1', 'ID_Cart.32X1'],
        'Cart.4X1'  :['ID_Cart.19X1', 'ID_Cart.5X1'],
        'Cart.5X1'  :['ID_Cart.19X1', 'ID_Cart.4X1', 'ID_Cart.0X0']
    }

    for b in Internal.getBases(tc):
        for z in Internal.getZones(b):
            if z[0] in dictOfDoublons:
                pos = 0
                z2 = Internal.copyRef(z)
                for zs in z2[2]:
                    if ('ID' in zs[0] or 'IBCD' in zs[0]) and zs[0] in dictOfDoublons[z[0]]:
                        Internal.addChild(z, zs, pos)
                        pos +=2
                    else:
                        pos += 1
    ####

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