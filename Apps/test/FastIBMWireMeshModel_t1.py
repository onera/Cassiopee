import Apps.Fast.IBM as App
import Converter.PyTree as C
import Converter.Internal as Internal
import KCore.test as test
import Fast.PyTree as Fast
import FastC.PyTree as FastC
import FastS.PyTree as FastS
import Post.PyTree as P
import Geom.PyTree as D
import numpy

LOCAL = test.getLocal()

##Vertical Line
t = Internal.newCGNSTree()
Internal.newCGNSBase('Base1', 3, 3, parent=t);
C._addState(t, 'EquationDimension', 2)
C._addState(t, 'GoverningEquations', 'NSTurbulent')
C._addState(t, 'TurbulenceModel', 'OneEquation_SpalartAllmaras')

base=Internal.getNodeByName(t,'Base1')
Internal.addChild(base, D.line((-0.09128554453599108,-0.19576248199991644,0), (0.09128554453599105,0.19576248199991644,0),N=800))

uinf            = 69.22970250694424*numpy.cos(4* numpy.pi/180)
Lcharac         = 0.03362355
C._addState(t, adim='dim4', UInf=uinf, TInf=298.15, PInf=101325,LInf=Lcharac,Mus=1.78938e-5)
App._setSnear(t, 0.0025)
App._setDfar(t, 0.75)
App._setIBCType(t, 'wiremodel')

tFile                 = LOCAL+'/t.cgns'
tcFile                = LOCAL+'/tc.cgns'

##PREP
dfar      = 5
snears    = 1
vmin      = 21;

t,tc=App.prepare1(t                , tFile      , tcFile             , 
                  snears=snears    , dfar=dfar  , vmin=vmin          ,
                  check=False      , frontType=1, isFilamentOnly=True, isWireModel=True )

test.testT(t , 1)
test.testT(tc, 2)

##NEEDED FOR THE MPI version
C.convertPyTree2File(t ,LOCAL+'/t_WMM.cgns')
C.convertPyTree2File(tc,LOCAL+'/tc_WMM.cgns')


##COMPUTE
NIT                        = 25  # number of iterations 
display_probe_freq         = 5   # iteration frequency to display modulo_verif

numb = {}
numb["temporal_scheme"]    = "implicit"
numb["ss_iteration"]       = 30
numb["modulo_verif"]       = display_probe_freq
numz = {}
numz["time_step"]          = 1.0e-4 # CFL 1
numz["time_step_nature"]   = "global"
numz["epsi_newton"]        = 0.1
Red                        = 360
beta                       = 0.64
kwire_local                = (0.5+26/Red)*(1.-beta*beta)/(beta*beta)
dv                         = numpy.sqrt((0.25*kwire_local)**2+1)-0.25*kwire_local
numz["DiameterWire_p"]     =0.08e-03
numz["CtWire_p"]           =0.065

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
graph=None
ts   =None
(t, tc, metrics) = FastS.warmup(t, tc, graph,tmy=ts)

file_select = 1
time_step   = Internal.getNodeFromName(t, 'time_step')
time_step   = Internal.getValue(time_step)


for it in range(NIT):
    FastS._compute(t, metrics, it, tc, graph,layer="Python")
    time0 += time_step
        
    if it%display_probe_freq == 0:
        print('- %d / %d - %f'%(it+it0, NIT+it0, time0))
        FastS.display_temporal_criteria(t, metrics, it, format='double')

# time stamp
Internal.createUniqueChild(t, 'Iteration', 'DataArray_t', value=it0+NIT)
Internal.createUniqueChild(t, 'Time', 'DataArray_t', value=time0)
Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
Internal._rmNodesByName(tc, '.Solver#Param')
Internal._rmNodesByName(tc, '.Solver#ownData')

test.testT(t , 3)
test.testT(tc, 4)

##TO VISUALIZE
#C.convertPyTree2File(t ,LOCAL+'/restart_WMM.cgns')
#C.convertPyTree2File(tc,LOCAL+'/tc_restart_WMM.cgns')
