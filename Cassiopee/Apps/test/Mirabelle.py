# - Mirabelle -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import Apps.Mesh.Mirabelle3 as Mirabelle

# Case
a = G.cart((0,0,0), (1,1,1), (30,30,20))
C._addBC2Zone(a, 'wall', 'BCWall', 'kmin')
C._fillEmptyBCWith(a, 'far', 'BCFarfield')
a = T.splitNParts(a, 12)
t = C.newPyTree(['Base',a])

# split full match
T._splitFullMatch(t)

C.convertPyTree2File(t, 'case.cgns')

print("== Computing linelet...")

# Trouve la premiere parois
#(block,win,zdir) = Mirabelle.getAWallWindow(t)

# Getting the linelet
zdir = 3; block = 'cart.37'
linelet = Mirabelle.buildDistrib(t, block, zdir, h1=0.01, h2=-1, N=100)

# getting match graph
g = Mirabelle.getGraph(t)

if zdir == 2: zdir = 3
elif zdir == 3: zdir = 5
stack = [(block, zdir)]

# Run
print("== Running propagate...")
treated = []
Mirabelle._propagate(t, g, stack, treated, linelet)

# Adapte les donneurs a la fin
Mirabelle._adaptDonorRanges(t)

C.convertPyTree2File(t, 'remesh.cgns')
