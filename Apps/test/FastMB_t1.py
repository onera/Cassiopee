# - Fast.MB -
import Apps.Fast.MB as App
import KCore.test as test

myApp = App.MB()
myApp.set(NP=0, format='single')

t, tc = myApp.prepare('naca.cgns', t_out='t.cgns', tc_out='tc.cgns')

test.testT(t, 1)
test.testT(tc, 2)

