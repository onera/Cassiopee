import Apps.Fast.IBMO as App
import Converter.Mpi as Cmpi
import Transform.PyTree as T
import Converter.Internal as Internal
import KCore.test as test
test.TOLERANCE = 1.e-8

LOCAL = test.getLocal()

NP = Cmpi.size
rank = Cmpi.rank
NIT = 100
FILE = 'naca_IBMO.cgns'

myApp = App.IBMO(format='single')
myApp.set(numb={"temporal_scheme": "explicit",
                "ss_iteration":5,
                "omp_mode":0})
myApp.set(numz={"time_step": 0.002,
                "scheme":"roe_min",
                "time_step_nature":"local",
                "cfl":0.5})

t,tc = myApp.prepare(FILE, t_out=LOCAL+'/t.cgns', tc_out=LOCAL+'/tc.cgns', expand=3, vmin=11, check=False, NP=Cmpi.size, distrib=True)
cartBase = Internal.getNodeFromName(t,'CARTESIAN')
Internal._rmNodesFromType(cartBase,'Rind_t')
test.testT(t,1)

t,tc = myApp.compute(LOCAL+'/t.cgns',LOCAL+'/tc.cgns', t_out=LOCAL+'/restart.cgns', tc_out=LOCAL+'/tc_restart.cgns', nit=NIT)
Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
cartBase = Internal.getNodeFromName(t,'CARTESIAN')
Internal._rmNodesFromType(cartBase,'Rind_t')

####
# The following lines are to avoid regression since the removal of sortByName in FastS warmup
####
Internal._sortByName(t, recursive=False)
cgnslibver = Internal.getNodeByType(t, 'CGNSLibraryVersion_t')
Internal._rmNodesByType(t, 'CGNSLibraryVersion_t')
Internal.addChild(t, cgnslibver, 0)
####
test.testT(t,2)

# Suppress since it doesnt test anything
#t = T.subzone(t,(1,1,1),(-1,-1,1))
#for nob in range(len(t[2])):
#    if t[2][nob][0] != 'CARTESIAN' and t[2][nob][3] == 'CGNSBase_t':
#        Internal._rmGhostCells(t,t[2][nob],2, adaptBCs=1)
#    else:
#        Internal._rmNodesFromType(t[2][nob],'Rind_t')
