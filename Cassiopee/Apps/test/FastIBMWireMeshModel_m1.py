## test case - Wire Mesh Model 
import Connector.IBM as X_IBM
import Apps.Fast.IBM as App
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Converter.Internal as Internal
import KCore.test as test
import Fast.PyTree as Fast
import FastC.PyTree as FastC
import FastS.Mpi as FastS
import Post.PyTree as P
import Geom.PyTree as D
import Geom.IBM as D_IBM
import Geom.Offset as D_Offset
import numpy
import os

test.TOLERANCE = 5e-7
LOCAL = test.getLocal()

tFile   = LOCAL+'/t_WMM.cgns'
tcFile  = LOCAL+'/tc_WMM.cgns'

##Geometry - Vertical Line
tb = Internal.newCGNSTree()
Internal.newCGNSBase('IBCFil_Base1', 3, 3, parent=tb);
C._addState(tb, 'EquationDimension', 2)
C._addState(tb, 'GoverningEquations', 'NSTurbulent')
C._addState(tb, 'TurbulenceModel', 'OneEquation_SpalartAllmaras')

base = Internal.getNodeByName(tb, 'IBCFil_Base1')
Internal.addChild(base, D.line((-0.09128554453599108,-0.19576248199991644,0), (0.09128554453599105,0.19576248199991644,0),N=800))

uinf         = 69.22970250694424*numpy.cos(4* numpy.pi/180)
Lcharac      = 0.03362355
C._addState(tb, adim='dim4', UInf=uinf, TInf=298.15, PInf=101325., LInf=Lcharac,Mus=1.78938e-5)
D_IBM._setSnear(tb, 0.0025)
D_IBM._setDfar(tb, 0.75)
D_IBM._setIBCType(tb, 'wiremodel')

tboffset = D_Offset.offsetSurface(tb, offset=0.025*3, pointsPerUnitLength=400, dim=2)
D_IBM._setSnear(tboffset, 0.0025)
tboffset = C.newPyTree(['Base', tboffset])

##PREP
dfars     = 5
snears    = 1
vmin      = 11

X_IBM.prepareIBMData(tb               , tFile        , tcFile   , tbox=tboffset,      
                     snears=snears    , dfars=dfars  , vmin=vmin, 
                     check=False       , frontType=1  , cartesian=False)
App._distribute(tFile, tcFile, NP=Cmpi.size)
t       = Fast.loadTree(os.path.basename(tFile), directory=LOCAL, split='single',  mpirun=True)
tc,graph= Fast.loadFile(tcFile, split='single',  mpirun=True, graph=True)

if Cmpi.rank == 0:
    test.testT(t , 1)
    test.testT(tc, 2)

##COMPUTE
NIT                        = 25   # number of iterations 
display_probe_freq         = 5    # iteration frequency to display modulo_verif
numb = {}
numb["temporal_scheme"]    = "implicit"
numb["ss_iteration"]       = 30
numb["modulo_verif"]       = display_probe_freq
numz = {}
numz["time_step"]          = 5.0e-5 
numz["time_step_nature"]   = "global"
numz["epsi_newton"]        = 0.1
Red                        = 360
beta                       = 0.64
kwire_local                = (0.5+26/Red)*(1.-beta*beta)/(beta*beta)
dv                         = numpy.sqrt((0.25*kwire_local)**2+1)-0.25*kwire_local
numz["DiameterWire_p"]     = 0.08e-03
numz["CtWire_p"]           = 0.065
numz["KWire_p"]            = kwire_local

dim=Internal.getValue(Internal.getNodeFromName(t, 'EquationDimension'))

it0 = 0; time0 = 0.
first = Internal.getNodeFromName1(t, 'Iteration')
if first is not None: it0 = Internal.getValue(first)
first = Internal.getNodeFromName1(t, 'Time')
if first is not None: time0 = Internal.getValue(first)

# Numerics
FastC._setNum2Base(t, numb)
FastC._setNum2Zones(t, numz)
ts = None
(t, tc, metrics) = FastS.warmup(t, tc, graph,tmy=ts)
graphInvIBCD_WM  = Cmpi.computeGraph(tc, type='INV_IBCD', procDict=graph['procDict'])

file_select = 1
time_step   = Internal.getNodeFromName(t, 'time_step')
time_step   = Internal.getValue(time_step)

for it in range(NIT):
    FastS._compute(t, metrics, it, tc, graph, layer="Python", graphInvIBCD_WM=graphInvIBCD_WM)
    time0 += time_step

    if it%display_probe_freq == 0:
        if Cmpi.rank==0: print('- %d / %d - %f'%(it+it0, NIT+it0, time0))
        FastS.display_temporal_criteria(t, metrics, it, format='double')

##TO VISUALIZE
#Cmpi.convertPyTree2File(t,LOCAL+'/t_restart_WMM_checkMPI.cgns')

if Cmpi.rank == 0:
    Internal._rmNodesFromType(t, 'Rind_t')
    Internal._rmNodesByName(t, '.Solver#Param')
    Internal._rmNodesByName(t, '.Solver#ownData')
    Internal._rmNodesByName(t, '.Solver#dtloc')
    Internal._rmNodesFromType(tc, 'Rind_t')
    Internal._rmNodesByName(tc, '.Solver#Param')
    Internal._rmNodesByName(tc, '.Solver#ownData')
    Internal._rmNodesByName(tc, '.Solver#dtloc')

    test.testT(t , 3)
    test.testT(tc, 4)

    os.remove(tFile)
    os.remove(tcFile)
