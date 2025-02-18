# - compute (pyTree) -
import Apps.Fast.IBM as Apps
import Apps.Fast.FastIBM as FastIBM
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import Converter.PyTree as C
import Fast.PyTree as Fast
import FastC.PyTree as FastC
import FastS.PyTree as FastS
import Geom.IBM as D_IBM
import Geom.PyTree as D
import KCore.test as test
import Connector.IBM as X_IBM
import os

test.TOLERANCE = 5.e-7
LOCAL           = test.getLocal()
FastC.MX_SYNCHRO= 1761

NP   = Cmpi.size
rank = Cmpi.rank

filename_tcase   = 'windTunnel.cgns'
DIRECTORY_PROBES = LOCAL+'/'
## =========================
##       PREP
## =========================
Diam       = 1.83
Ltube      = 15.
Cp         = 1004.0
Rgaz       = 287.06
Mu0        = 1.791e-5
Psi        = 0.02
Mach       = 0.082
Pgen0      = 94450.0
Tgen0      = 297.2
Pout       = 93970.6
Psta0      = Pgen0/(1+0.2*Mach*Mach)**(3.5)
Tsta0      = Tgen0/(1+0.2*Mach*Mach)
Cson       = (1.4*Rgaz*Tsta0)**(0.5)
Vref       = Mach*Cson
Ro0        = Psta0/Tsta0/Rgaz
Stag_Enthal= Cp*Tgen0
Pout       = Psta0-0.5*Ro0*Vref*Vref*Psi*Ltube/Diam
Red        = Ro0*Vref*Diam/Mu0
ReL        = Ro0*Vref*Ltube/Mu0

tb = C.convertFile2PyTree(filename_tcase)

for base in Internal.getBases(tb):
    fl = Internal.newFlowEquationSet(parent=base)
    gov = Internal.newGoverningEquations(value='NSTurbulent', parent=fl)

for node in Internal.getNodesByName(tb,'FlowEquationSet'):
    Internal._createUniqueChild(node, 'EquationDimension', '"int"', value=2, children=[])
    Internal._createUniqueChild(node,'TurbulenceModel', 'TurbulenceModel_t', value='OneEquation_SpalartAllmaras', children=[])

C._addState(tb, adim='dim1', UInf=Vref, TInf=Tsta0, PInf=Psta0, LInf=Diam, alphaZ=0., alphaY=0., MutSMuInf=0.01, TurbLevelInf=0.0015)
D_IBM._snearFactor(tb, 3)
itExtrctPrb = 2
D_IBM._setOutPressControlParam(tb,probeName='point', AtestSection=0.83721966959, AOutPress=0.83721966959,
                               machTarget=0.082, pStatTarget=99445, tStatTarget=297.2,lmbd=0.1,
                               cxSupport=0.6, sSupport=0.1, itExtrctPrb=itExtrctPrb)
test.testT(tb,1)

t,tc=X_IBM.prepareIBMData(tb         , None   , None   ,
                          snears=0.01, dfars=0, vmin=11, cartesian=False)

[RoInf, RouInf, RovInf, RowInf, RoeInf, PInf, TInf, cvInf, MInf,
 ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
 Mus, Cs, Ts, Pr] = C.getState(tc)
Tsta    = Tgen0/(1.+0.5*(Gamma-1)*Mach*Mach);
Psta    = Pgen0/(1.+0.5*(Gamma-1)*Mach*Mach)**(Gamma/(Gamma-1));
Ro      = Psta/Tsta/Rgaz
Cson    = (1.4*Rgaz*Tsta)**(0.5)
Vref    = Mach*Cson
Pout    = Psta-0.5*Ro*Vref*Vref*Psi*Ltube/Diam
stagEnt = Gamma*cvInf*Tgen0
Apps._initOutflow(tc, 'outlet', Pout, isDensityConstant=False)
Apps._initInj(tc, 'inlet', Pgen0, stagEnt, injDir=[1.,0.,0.])

test.testT(tc,2)
test.testT(t ,3)

a = D.point((-1,0,0))
C.convertPyTree2File(a,DIRECTORY_PROBES+'probes.cgns')

## =========================
##       COMPUTE
## =========================

NIT   = 20
graph = None

it0 = 0; time0 = 0.
first = Internal.getNodeFromName1(t, 'Iteration')
if first is not None: it0 = Internal.getValue(first)
first = Internal.getNodeFromName1(t, 'Time')
if first is not None: time0 = Internal.getValue(first)

# Numerics
modulo_verif       =  5
buffer_size        =  100
numb = {}
numb["temporal_scheme"]    = "implicit"
numb["ss_iteration"]       = 2
numb["modulo_verif"]       = modulo_verif
numb["omp_mode"]           = 0
numz = {}
numz["time_step_nature"]   = "local"
numz["cfl"]                = 5
numz["scheme"]             = "roe_min"
numz["psiroe"]             = 0.01

Fast._setNum2Zones(t, numz); Fast._setNum2Base(t, numb)

[RoInf, RouInf, RovInf, RowInf, RoeInf, PInf, TInf, cvInf, MInf,
 ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
 Mus, Cs, Ts, Pr] = C.getState(tc)

(t, tc, metrics) = FastS.warmup(t, tc, graph=graph)


#0: Start of P mod algo.     [# iter] (multiple of probe recording value)
#1: Frequency of P mod algo. [# iter] (same as the probe recording value)
#2: WT response time         [# iter]
##Note: Frequency of P mod algo.= 0.1* WT response time [good default values]
##E.g.  WT response time = 2000 & Frequency of P mod algo. = 200
itValues4gain=[5,2,100]

isRestartProbe = False
values4gain,controlProbeName,itExtrctPrb=FastIBM.getInfoOutletPressure(tb, familyName='outlet')
FastIBM._setUpOutletPressure(values4gain, itValues4gain)

probe_in = os.path.join(DIRECTORY_PROBES, "probes.cgns")
dictOfPointProbes = FastIBM.initPointProbes(t, probe_in, fields=['centers:Mach'], bufferSize=buffer_size, append=isRestartProbe, historyDirectory=DIRECTORY_PROBES)

for it in range(NIT):
    FastS._compute(t, metrics, it, tc, graph=graph, layer='Python')
    if it%modulo_verif == 0:
        FastS.display_temporal_criteria(t, metrics, it)
    if it%itExtrctPrb == 0:
        FastIBM._updatePointProbes(t, dictOfPointProbes, it, ['centers:Mach'])
        if it>itValues4gain[0] and it%itValues4gain[1] == 0:
            FastIBM._controlOutletPressureMachProbe(tc,dictOfPointProbes,controlProbeName,values4gain,it,familyName='outlet')

#C.convertPyTree2File(tc,LOCAL+'/tc_restart.cgns')
#C.convertPyTree2File(t ,LOCAL+'/t_restart.cgns' )

for name, probe in dictOfPointProbes.items(): probe.flush()
os.remove(DIRECTORY_PROBES+'/probe_point.cgns')
os.remove(DIRECTORY_PROBES+'/probes.cgns')

Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')

Internal._rmNodesByName(tc, '.Solver#Param')
Internal._rmNodesByName(tc, '.Solver#ownData')
Internal._rmNodesFromName(tc, 'Parameter_int')
Internal._rmNodesFromName(tc, 'Parameter_real')

test.testT(tc,4)
#test.testT(t ,5)
