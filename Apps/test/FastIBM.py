# - Fast.IBM -
import Apps.Fast.IBM as App

myApp = App.IBM(NP=0, format='single')
myApp.set(numb={"temporal_scheme": "implicit",
                "ss_iteration":3})
myApp.set(numz={"time_step": 0.0007,
                "scheme":"roe_min",
                "time_step_nature":"local",
                "cfl":4.})

# Prepare
myApp.prepare('naca1D.cgns', t_out='t.cgns', tc_out='tc.cgns')
