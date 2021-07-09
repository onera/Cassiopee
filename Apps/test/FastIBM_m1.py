# - Fast.IBM -
# Euler, para, frontType=1
import Apps.Fast.IBM as App
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import KCore.test as test
import Converter.Internal as Internal
test.TOLERANCE = 1.e-6

LOCAL = test.getLocal()

myApp = App.IBM(format='single')
myApp.set(numb={"temporal_scheme": "implicit",
                "ss_iteration":3,
                "omp_mode":0})
myApp.set(numz={"time_step": 0.0007,
                "scheme":"roe_min",
                "time_step_nature":"local",
                "cfl":4.})

# Prepare
t, tc = myApp.prepare('naca1DEuler.cgns', t_out=LOCAL+'/t.cgns', tc_out=LOCAL+'/tc.cgns', NP=Cmpi.size)
Internal._rmNodesFromName(tc,Internal.__GridCoordinates__)
Internal._rmNodesFromType(tc,'Rind_t')
if Cmpi.rank == 0: test.testT(tc, 1)

# Compute
t,tc = myApp.compute(LOCAL+'/t.cgns', LOCAL+'/tc.cgns', t_out=LOCAL+'/restart.cgns', tc_out=LOCAL+'/tc_restart.cgns', nit=300)

if Cmpi.rank == 0:
    t = C.convertFile2PyTree(LOCAL+'/restart.cgns')
    Internal._rmNodesFromType(t,'Rind_t')
    test.testT(t, 2)
