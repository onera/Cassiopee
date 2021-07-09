# - Fast.IBMO -
# Fast Overset+IBM application
import Apps.Fast.IBMO as App
import Converter.Mpi as Cmpi
import Transform.PyTree as T
import Converter.Internal as Internal

NP = Cmpi.size
rank = Cmpi.rank
NIT = 300
FILE = 'naca_IBMO.cgns'
myApp = App.IBMO(format='single')
myApp.set(numb={"temporal_scheme": "implicit",
                "ss_iteration":5,
                "omp_mode":0})
myApp.set(numz={"time_step": 0.002,
                "scheme":"roe_min",
                "time_step_nature":"local",
                "cfl":1.})

t,tc = myApp.prepare(FILE, t_out='t.cgns', tc_out='tc.cgns', expand=3, vmin=11, check=False, NP=Cmpi.size, distrib=True)
t,tc = myApp.compute('t.cgns','tc.cgns', t_out='restart.cgns', tc_out='tc_restart.cgns', nit=NIT)

t = T.subzone(t,(1,1,1),(-1,-1,1))
for nob in range(len(t[2])):
    if t[2][nob][0] != 'CARTESIAN' and t[2][nob][3]=='CGNSBase_t':
        Internal._rmGhostCells(t,t[2][nob],2, adaptBCs=1)

Cmpi.convertPyTree2File(t,"out.cgns")
