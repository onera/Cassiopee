## test case - Wire Mesh Model - OpenMP only
## using a slanted line to make the test the MPI with the same geometry
import Connector.IBM as X_IBM
import Converter.PyTree as C
import Converter.Internal as Internal
import KCore.test as test
import Fast.PyTree as Fast
import FastC.PyTree as FastC
import FastS.PyTree as FastS
import Post.PyTree as P
import Geom.PyTree as D
import Geom.IBM as D_IBM
import Geom.Offset as D_Offset
import numpy

LOCAL = test.getLocal()

##Geometry - Vertical Line
tb = Internal.newCGNSTree()
Internal.newCGNSBase('IBCFil_Base1', 3, 3, parent=tb)
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
#C.convertPyTree2File(tb,LOCAL+'/tb_check.cgns')

tboffset = D_Offset.offsetSurface(tb, offset=0.025*3, pointsPerUnitLength=400, dim=2)
D_IBM._setSnear(tboffset, 0.0025)
tboffset = C.newPyTree(['Base', tboffset])

##PREP
dfars = 5
snears = 1
vmin = 11

t,tc = X_IBM.prepareIBMData(tb               , None         , None     , tbox=tboffset,      
                            snears=snears    , dfars=dfars  , vmin=vmin, 
                            check=False      , frontType=1  , cartesian=False)
test.testT(t , 1)
test.testT(tc, 2)

##COMPUTE
NIT           = 25  # number of iterations 
modulo_verif  = 5   # iteration frequency to display modulo_verif

numb = {}
numb["temporal_scheme"]    = "implicit"
numb["ss_iteration"]       = 30
numb["modulo_verif"]       = modulo_verif
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
(t, tc, metrics) = FastS.warmup(t, tc, graph, tmy=ts)

time_step   = Internal.getNodeFromName(t, 'time_step')
time_step   = Internal.getValue(time_step)

for it in range(NIT):
    FastS._compute(t, metrics, it, tc, graph, layer="Python")
    time0 += time_step
        
    if it%modulo_verif == 0:
        print('- %d / %d - %f'%(it+it0, NIT+it0, time0))
        FastS.display_temporal_criteria(t, metrics, it)

# time stamp
Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
Internal._rmNodesByName(tc, '.Solver#Param')
Internal._rmNodesByName(tc, '.Solver#ownData')

test.testT(t , 3)
test.testT(tc, 4)

##TO VISUALIZE
#C.convertPyTree2File(t,LOCAL+'/t_restart_WMM_check.cgns')
