# - Mesh.Mirabelle3 -
# propagate a new edge distribution imposing the edge values
import Converter.PyTree as C
import Transform.PyTree as T
import Apps.Mesh.Mirabelle3 as Mirabelle
import Geom.PyTree as D
import KCore.test as test
t = C.convertFile2PyTree("premirabelle.cgns")

# Getting the linelet
zdir    = 2;
block   = 'zone.6' # stator
N       = 53
h1      = 8.0e-6#-1
h2      = -1
isAvg   = True

linelet = D.line((0,0,0), (1,0,0),N=N)

# getting match graph
g = Mirabelle.getGraph(t)

if zdir == 2: zdir = 3
elif zdir == 3: zdir = 5
stack = [(block, zdir)]

# Run
treated = []
Mirabelle._propagate(t, g, stack, treated, linelet, h1, h2,isAvg=isAvg,nAvg=2)

t = C.rmBCOfType(t, 'BC*') # rm Wall BCs
t = C.rmBCOfType(t, 'BCMatch') # rm Wall BCs
t = C.rmBCOfType(t, 'BCOverlap') # rm Wall BCs

T._makeDirect(t)
test.testT(t, 1)
