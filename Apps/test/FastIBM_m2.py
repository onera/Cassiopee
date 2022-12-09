# - Fast.IBM -
# Euler, para, frontType=1
import Apps.Fast.IBM as App
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import KCore.test as test
import Converter.Internal as Internal
import Geom.PyTree as D
import FastC.PyTree as FastC

test.TOLERANCE = 1.e-8
LOCAL = test.getLocal()

FILEB = LOCAL+"/case.cgns"
FILEC = LOCAL+"/tc_restart.cgns"
FILED = LOCAL+"/tcw.cgns"

myApp = App.IBM(format='single')
myApp.set(numb={"temporal_scheme": "implicit",
                "ss_iteration":3,
                "omp_mode":0})
myApp.set(numz={"time_step": 0.0007,
                "scheme":"roe_min",
                "time_step_nature":"local",
                "cfl":4.})
# case
a = D.sphere6((0.,0.,0.),1.,N=20)
for z in a:
    App._setSnear(z,0.1)
    App._setDfar(z,10.)
    App._setIBCType(z,'Musker')
    
tb = C.newPyTree(['Body',a])
for base in Internal.getBases(tb):
    fl = Internal.newFlowEquationSet(parent=base)
    gov = Internal.newGoverningEquations(value='NSTurbulent', parent=fl)
eqdim = Internal.createNode('EquationDimension', '"int"', value=3, children=[])
turbmod = Internal.createNode('TurbulenceModel', 'TurbulenceModel_t', value='OneEquation_SpalartAllmaras', children=[])
for node in Internal.getNodesByName(tb,'FlowEquationSet'):
    Internal.addChild(node,eqdim)
    Internal.addChild(node,turbmod)
    
C._addState(tb, adim='adim1', MInf=0.1, alphaZ=0., alphaY=0., ReInf=40000, MutSMuInf=0.1, TurbLevelInf=1.e-4)
tb = C.convertArray2Tetra(tb)
noz = 0
for z in Internal.getZones(tb):
    Cmpi._setProc(z,noz%(Cmpi.size))
    noz += 1
    
if Cmpi.rank==0: C.convertPyTree2File(tb, FILEB)
Cmpi.barrier()

# Prepare
myApp.input_var.NP=Cmpi.size
t,tc = myApp.prepare(FILEB, t_out=LOCAL+'/t.cgns', tc_out=LOCAL+'/tc.cgns')
Internal._rmNodesFromType(tc,'Rind_t')
Internal._rmNodesFromName(tc,Internal.__GridCoordinates__)
if Cmpi.rank == 0: test.testT(tc, 1)

# Compute
t,tc = myApp.compute(LOCAL+'/t.cgns', LOCAL+'/tc.cgns', t_out=LOCAL+'/restart.cgns', tc_out=FILEC, nit=100)
Cmpi.barrier()

if Cmpi.rank == 0:
    t = C.convertFile2PyTree(LOCAL+'/restart.cgns')
    Internal._rmNodesFromType(t, 'Rind_t')
    Internal._rmNodesByName(t, '.Solver#Param')
    Internal._rmNodesByName(t, '.Solver#ownData')
    Internal._rmNodesByName(t, '.Solver#dtloc')
    test.testT(t,2)
    
procDictR = Cmpi.getProcDict(tb)
Cmpi._convert2PartialTree(tb,rank=Cmpi.rank)

tcw ,graphP= App._prepareSkinReconstruction(tb,tc)
App._computeSkinVariables(tb,tc, tcw, graphP)
Cmpi.convertPyTree2File(tb, LOCAL+'/wall.cgns')
Cmpi.barrier()
if Cmpi.rank==0:
    tb = C.convertFile2PyTree(LOCAL+'/wall.cgns')
    test.testT(tb,3)
