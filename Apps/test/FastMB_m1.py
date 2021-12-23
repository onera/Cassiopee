# - Fast.MB -
import Apps.Fast.MB as App
import KCore.test as test
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.Mpi as Cmpi

LOCAL = test.getLocal()

myApp = App.MB()
myApp.set(format='single')
myApp.set(numb={"temporal_scheme": "implicit",
                "ss_iteration":3})
myApp.set(numz={"time_step": 0.0007,
                "scheme":"roe_min",
                "time_step_nature":"local",
                "cfl":4.})

if Cmpi.rank == 0: # prep en seq pour l'instant
    t, tc = myApp.prepare('naca.cgns', t_out=LOCAL+'/t.cgns', tc_out=LOCAL+'/tc.cgns', NP=Cmpi.size)
    test.testT(tc, 2)
Cmpi.barrier()

t, tc = myApp.compute(LOCAL+'/t.cgns', LOCAL+'/tc.cgns', t_out=LOCAL+'/restart.cgns', nit=300)

Cmpi.barrier()
if Cmpi.rank == 0:
    t = C.convertFile2PyTree(LOCAL+'/restart.cgns')
    Internal._rmNodesByName(t, '.Solver#Param')
    Internal._rmNodesByName(t, '.Solver#ownData')
    Internal._rmNodesByName(t, '.Solver#dtloc')
    test.testT(t, 1)
