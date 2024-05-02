# - Fast.MB -
import Apps.Fast.MB as App
import Converter.PyTree as C
import Converter.Internal as Internal
import KCore.test as test

LOCAL = test.getLocal()

myApp = App.MB()
myApp.set(format='single')
myApp.set(numb={"temporal_scheme": "implicit",
                "ss_iteration":3})
myApp.set(numz={"time_step": 0.0007,
                "scheme":"roe_min",
                "time_step_nature":"local",
                "cfl":4.})

t, tc = myApp.prepare('naca.cgns', t_out=LOCAL+'/t.cgns', tc_out=LOCAL+'/tc.cgns', NP=0)
test.testT(tc, 2)

t, tc = myApp.compute(LOCAL+'/t.cgns', LOCAL+'/tc.cgns', t_out=LOCAL+'/restart.cgns', nit=300)
t = C.convertFile2PyTree(LOCAL+'/restart.cgns')
Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
test.testT(t, 1)
