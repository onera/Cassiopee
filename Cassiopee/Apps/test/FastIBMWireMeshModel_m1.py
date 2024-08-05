# ...
import Apps.Fast.IBM as App
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Converter.Internal as Internal
import KCore.test as test
import Fast.PyTree as Fast
import FastC.PyTree as FastC
import FastS.Mpi as FastS
import Post.PyTree as P
import numpy

LOCAL = test.getLocal()

tFile   = LOCAL+'/t_WMM.cgns'
tcFile  = LOCAL+'/tc_WMM.cgns'

##READING SERIAL & DISTRIBUTING
App._distribute(tFile, tcFile, NP=Cmpi.size)

t       = Fast.loadTree(tFile , split='single', directory=LOCAL, mpirun=True)
tc,graph= Fast.loadFile(tcFile, split='single',  mpirun=True, graph=True)

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

file_select = 1
time_step   = Internal.getNodeFromName(t, 'time_step')
time_step   = Internal.getValue(time_step)

for it in range(NIT):
    FastS._compute(t, metrics, it, tc, graph, layer="Python")
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
    test.testT(t , 1)
    test.testT(tc, 2)
    

