# - FastIBMO -
import Apps.Fast.IBMO as App
import Converter.Mpi as Cmpi
import Transform.PyTree as T
import Converter.Internal as Internal
import KCore.test as test

test.TOLERANCE = 1.e-6

LOCAL = test.getLocal()
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
if Cmpi.rank == 0: test.testT(t,1)
Cmpi.barrier()

t,tc = myApp.compute(LOCAL+'/t.cgns',LOCAL+'/tc.cgns', t_out=LOCAL+'/restart.cgns', tc_out=LOCAL+'/tc_restart.cgns', nit=100)
if Cmpi.rank == 0:
    Internal._rmNodesByName(t, '.Solver#Param')
    Internal._rmNodesByName(t, '.Solver#ownData')
    Internal._rmNodesByName(t, '.Solver#dtloc')
    test.testT(t,2)
Cmpi.barrier()

t = T.subzone(t,(1,1,1),(-1,-1,1))
for nob in range(len(t[2])):
    if t[2][nob][0] != 'CARTESIAN' and t[2][nob][3] == 'CGNSBase_t':
        Internal._rmGhostCells(t,t[2][nob],2, adaptBCs=1)
if Cmpi.rank == 0: test.testT(t,3)
Cmpi.barrier()
