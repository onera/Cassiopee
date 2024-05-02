# - Fast.IBM -
# Fast IBM application
import Apps.Fast.IBM as App

myApp = App.IBM(format='single')
myApp.set(numb={"temporal_scheme": "implicit",
                "ss_iteration":3,
                "omp_mode":1})
myApp.set(numz={"time_step": 0.0007,
                "scheme":"roe_min",
                "time_step_nature":"local",
                "cfl":4.})

# Prepare
myApp.prepare('naca1DEuler.cgns', t_out='t.cgns', tc_out='tc.cgns')

# Compute
myApp.compute('t.cgns', 'tc.cgns', t_out='restart.cgns', tc_out='tc_restart.cgns', nit=300)

# Post
myApp.post('naca1DEuler.cgns', 'restart.cgns', 'tc_restart.cgns', t_out='out.cgns', wall_out='wall.cgns')
