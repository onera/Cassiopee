# - Fast.MB -
# Fast Multiblock application
import Apps.Fast.MB as App
import Converter.Mpi as Cmpi

myApp = App.MB(format='single')
myApp.set(numb={"temporal_scheme": "implicit",
                "ss_iteration":3})
myApp.set(numz={"time_step": 0.0007,
                "scheme":"roe_min",
                "time_step_nature":"local",
                "cfl":4.})

# Prepare
myApp.prepare('naca.cgns', t_out='t.cgns', tc_out='tc.cgns', NP=Cmpi.size)

# Compute
myApp.compute('t.cgns', 'tc.cgns', t_out='restart.cgns', nit=300)

# Post
myApp.post('restart.cgns', t_out='out.cgns', wall_out='wall.cgns')
