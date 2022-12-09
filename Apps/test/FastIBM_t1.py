# - Fast.IBM -
# Euler, seq, frontType=1
import Apps.Fast.IBM as App
import Converter.PyTree as C
import Converter.Internal as Internal
import KCore.test as test

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
myApp.input_var.NP = 1
t, tc = myApp.prepare('naca1DEuler.cgns', t_out=LOCAL+'/t.cgns', tc_out=LOCAL+'/tc.cgns')
Internal._rmNodesFromName(tc, Internal.__GridCoordinates__)
test.testT(tc, 1)

# Compute
t,tc = myApp.compute(LOCAL+'/t.cgns', LOCAL+'/tc.cgns', t_out=LOCAL+'/restart.cgns', tc_out=LOCAL+'/tc_restart.cgns', nit=300)
t = C.convertFile2PyTree(LOCAL+'/restart.cgns')
Internal._rmNodesFromName(t, '.Solver#Param')
Internal._rmNodesFromName(t, '.Solver#ownData')
Internal._rmNodesFromName(t, '.Solver#dtloc')
Internal._rmNodesFromType(t, 'Rind_t')
test.testT(t, 2)

# Post
t, zw = myApp.post('naca1DEuler.cgns', LOCAL+'/restart.cgns', LOCAL+'/tc_restart.cgns', t_out=LOCAL+'/out.cgns', wall_out=LOCAL+'/wall.cgns')
Internal._rmNodesFromName(t, '.Solver#dtloc')
Internal._rmNodesFromName(zw, '.Solver#dtloc')
test.testT(t, 3)
test.testT(zw, 4)
