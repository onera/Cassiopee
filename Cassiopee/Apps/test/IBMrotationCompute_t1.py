# case - cylinder - IBM rotation - compute
import Connector.IBM as X_IBM
import Connector.Mpi as Xmpi
import Connector.PyTree as X
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import Converter.PyTree as C
import Fast.PyTree as Fast
import FastS.PyTree as FastS
import RigidMotion.PyTree as R
import Transform.PyTree as T
import KCore.test as test
import math, os
import numpy

LOCAL = test.getLocal()

from IBMrotation import *

## =======================
## =======================
UInf = 0.0001

# necessaire pour les IBM
Model        = 'NSLaminar'
dimPb        = 2
modulo_verif = 5

rpm  = 100
OMG  = rpm*2*math.pi/60
time_step= 1/(OMG*180/math.pi)/600
## =======================
## =======================

NIT = 10
numb={"temporal_scheme": "implicit_local",
      "omp_mode":0,
      "modulo_verif":modulo_verif}

numz={"time_step": time_step,
      "time_step_nature": "global",
      "epsi_newton": 0.1}

numz1 = numz.copy()
numz1['scheme'] = 'ausmpred'

numz2 = numz.copy()
numz2['scheme'] = 'ausmpred'

# load t, tc
fileTIn ='/tRotComp.cgns'
fileTcIn='/tcRotComp.cgns'
t,tc,ts,graph = Fast.load(LOCAL+fileTIn,LOCAL+fileTcIn)
baseNameIBM   ='CARTESIAN_NEARBODY'
baseNameBKGD  ='CARTESIAN_OFFBODY'


Fast._setNum2Base(t, numb)
for n in [baseNameIBM]:
    b = Internal.getNodeFromName1(t, n)
    Fast._setNum2Zones(b, numz1)
for n in [baseNameBKGD]:
    b = Internal.getNodeFromName1(t, n)
    Fast._setNum2Zones(b, numz2)


# load bodies
tb = C.convertFile2PyTree(LOCAL+'/bodiesBlankRotComp.cgns')

# Warmup
it0 = 0; time0 = 0.
first = Internal.getNodeFromName1(t, 'Iteration')
if first is not None: it0 = Internal.getValue(first)
first = Internal.getNodeFromName1(t, 'Time')
if first is not None: time0 = Internal.getValue(first)
time_step = Internal.getNodeFromName(t, 'time_step')
time_step = Internal.getValue(time_step)

# Tag les zones pour le motion
zones = Internal.getZones(t)
for z in zones:
    timeMotion = Internal.getNodeFromName(z, 'TimeMotion')
    if timeMotion is not None:
        define = Internal.getNodeFromName1(z, '.Solver#define')
        Internal._createUniqueChild(define, 'motion', 'DataArray_t', value='deformation')

# Set interpolated points sur les bases IBM
C._initVars(t, "{centers:cellN#Motion}=1.")

##IBM CART
X._applyBCOverlaps(Internal.getNodeFromName1(t,baseNameIBM), depth=2, val=2, cellNName='cellN#Motion')
C._initVars(t, "{centers:cellN#MotionInit}={centers:cellN#Motion}")

# Data for chimera motion
nbases  = len(Internal.getBases(t))
nbodies = len(Internal.getBases(tb))
BM      = numpy.zeros((nbases, nbodies),dtype=numpy.int32)
BM[1,0] = 1.

# Creation BB tree commun
tBB = Cmpi.createBBoxTree(t)

# Attachement des donneurs en dur (id on all procs)
intersectionDict={}
for z1 in Internal.getZones(Internal.getNodeFromName1(tBB,baseNameIBM)):
    for z2 in Internal.getZones(Internal.getNodeFromName1(tBB,baseNameBKGD)):
        Fast._addPair(intersectionDict, z1[0], z2[0])

procDict = None;
graphX   = {}
procDict = Cmpi.getProcDict(tBB)
graphX   = Cmpi.computeGraph(tBB, type='bbox2', t2=None, procDict=procDict, reduction=False,intersectionsDict=intersectionDict)

R._evalPositionIBC(tc,time0)

(t, tc, metrics) = FastS.warmup(t, tc, graph)

# Doivent etre apres le warmup
dictOfADT={}
(dictOfNobOfRcvZones,dictOfNozOfRcvZones)   = Fast.getDictOfNobNozOfRcvZones(t, intersectionDict)
(dictOfNobOfRcvZonesC,dictOfNozOfRcvZonesC) = Fast.getDictOfNobNozOfRcvZones(tc, intersectionDict)
(dictOfNobOfDnrZones,dictOfNozOfDnrZones)   = Fast.getDictOfNobNozOfDnrZones(tc, intersectionDict, dictOfADT,cartFilter='CARTESIAN',isIbmAle=True)

# modifie la vitesse de reference (si M=0)
cont  = Internal.getNodeFromType(t, 'ReferenceState_t')
Minf  = Internal.getNodeFromName(cont, 'Mach')
zones = Internal.getZones(t)
for z in zones:
    n = Internal.getNodeFromName2(z, 'Parameter_real')[1]
    n[5] = max(30, UInf)

timeiter  = time0
varType   = 0

# Get models
eqs    = Internal.getNodeFromType(Internal.getNodeFromName1(t,baseNameIBM), 'GoverningEquations_t')
Model1 = Internal.getValue(eqs)
eqs    = Internal.getNodeFromType(Internal.getNodeFromName1(t, baseNameBKGD), 'GoverningEquations_t')
Model2 = Internal.getValue(eqs)

if Model1 == 'NSTurbulent' and Model2 == 'NSTurbulent': varType = 1

RefState = Internal.getNodeFromType(t,'ReferenceState_t')
VInf     = Internal.getNodeFromName1(RefState,'VelocityX')[1][0]
PInf     = Internal.getNodeFromName1(RefState,'Pressure')[1][0]
rhoInf   = Internal.getNodeFromName1(RefState,'Density')[1][0]

for it in range(it0,NIT+it0):
    if Cmpi.rank == 0: print("it=%d, time=%f, angle=%f, N_rot=%d"%(it, timeiter, (timeiter*OMG*180./numpy.pi) % 360 , (timeiter*OMG*180./numpy.pi)//360 ),flush=True)
    timeiter += time_step

    # bouge tout
    R._evalPosition(tb, timeiter)
    R._evalPosition(t, timeiter)
    R._evalPosition(tc, timeiter)
    R._evalGridSpeed(t, timeiter)

    R._evalPositionIBC(tc,timeiter)

    FastS.copy_velocity_ale(t, metrics, it=it) # get motion velocities @ face centers

    # Masquage
    C._initVars(t, "{centers:cellN#Motion}={centers:cellN#MotionInit}")
    bodies=[]
    for base in Internal.getBases(tb):
        bodies.append(Internal.getZones(base))

    XRAYDIM1 = 2000;RAYDIM2 = XRAYDIM1
    if dimPb == 2: XRAYDIM2 = 2
    X._blankCells(t, bodies, BM, cellNName='cellN#Motion',XRaydim1=XRAYDIM1, XRaydim2=XRAYDIM2, blankingType='cell_intersect',dim=dimPb)
    X._setHoleInterpolatedPoints(Internal.getNodeFromName1(t,baseNameBKGD), depth=2, loc='centers',addGC=False, cellNName='cellN#Motion', dir=2)
    C._initVars(t, "{centers:cellN#Motion}=({centers:cellN#Static}<2)*{centers:cellN#Motion}+({centers:cellN#Static}>1)*minimum(1,{centers:cellN#Motion})")
    C._initVars(t,"{centers:cellN}=minimum({centers:cellN#Motion}*{centers:cellN#Static},2.)")
    C._cpVars(t, 'centers:cellN',tc, 'cellN')
    ucData = (graphX, intersectionDict, dictOfADT,
              dictOfNobOfRcvZones, dictOfNozOfRcvZones,
              dictOfNobOfDnrZones, dictOfNozOfDnrZones,
              dictOfNobOfRcvZonesC, dictOfNozOfRcvZonesC,
              timeiter, procDict, True, varType, 1, 2, 1)

    FastS._compute(t, metrics, it, tc, graph, layer="Python", ucData=ucData)

    if it%modulo_verif==0:
        FastS.display_temporal_criteria(t, metrics, it, format='double')

    Cmpi.barrier()

for adt in dictOfADT.values():
    if adt is not None: C.freeHook(adt)

Internal.createUniqueChild(t, 'Iteration', 'DataArray_t', value=NIT+it0)
Internal.createUniqueChild(t, 'Time', 'DataArray_t', value=timeiter)

Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
Internal._rmNodesByName(tc, '.Solver#Param')
Internal._rmNodesByName(tc, '.Solver#ownData')

test.testT(t,1)
test.testT(tc,2)

#Cmpi.convertPyTree2File(t , LOCAL+'/t_restart.cgns')
#Cmpi.convertPyTree2File(tc, LOCAL+'/tc_restart.cgns')

os.remove(LOCAL+'/tRotComp.cgns')
os.remove(LOCAL+'/tcRotComp.cgns')
os.remove(LOCAL+'/bodiesBlankRotComp.cgns')
