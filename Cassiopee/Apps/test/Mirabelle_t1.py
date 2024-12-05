# - Mesh.Mirabelle3 -
# propagate a new edge distribution using buildDistrib
import Converter.PyTree as C
import Transform.PyTree as T
import Apps.Mesh.Mirabelle3 as Mirabelle
import KCore.test as test
t = C.convertFile2PyTree("premirabelle.cgns")

# Getting the linelet
zdir  = 2;
block = 'zone.6' # stator
N     = 53
h1    = 8.0e-6
h2    = -1
linelet = Mirabelle.buildDistrib(t, block, zdir, h1=h1, h2=h2, N=N)

# getting match graph
g = Mirabelle.getGraph(t)

if zdir == 2: zdir = 3
elif zdir == 3: zdir = 5
stack = [(block, zdir)]

# Run
treated = []
Mirabelle._propagate(t, g, stack, treated, linelet)

t = C.rmBCOfType(t, 'BC*') # rm Wall BCs
t = C.rmBCOfType(t, 'BCMatch') # rm Wall BCs
t = C.rmBCOfType(t, 'BCOverlap') # rm Wall BCs

T._makeDirect(t)
test.testT(t, 1)
