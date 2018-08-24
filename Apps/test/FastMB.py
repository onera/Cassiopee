# - Fast.MB -
import Apps.Fast.MB as App

myApp = App.MB()
myApp.set(NP=0, format='single')
myApp.set(numb={"temporal_scheme": "implicit",
                "ss_iteration":3})
myApp.set(numz={"time_step": 0.0007,
                "scheme":"roe_min",
                "time_step_nature":"local",
                "cfl":4.})

myApp.prepare('naca.cgns', t_out='t.cgns', tc_out='tc.cgns')
myApp.compute('t.cgns', 'tc.cgns', t_out='restart.cgns', nit=1000)

