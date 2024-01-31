import Converter.PyTree as C
import Generator.PyTree as G
import XCore.PyTree as X
import Converter.Mpi as Cmpi
import Converter.Internal as I
import numpy as np
import Intersector.PyTree as XOR

# mpirun -np <N> python3 <script.py>

def Func(x, y, z):
    a = 0.45
    b = 0.25
    dx2 = ((x-0.5)/a)**2 
    dy2 = ((y-0.5)/b)**2
    sr = 1.
    br = 1.2
    return float(dx2 + dy2 > sr and dx2 + dy2 < br)

def make_field(cx, cy, cz, F):
    ncells = len(cx)
    a = np.empty(ncells)
    for i in range(ncells):
        a[i] = F(cx[i], cy[i], cz[i])
    return a

if Cmpi.rank == 0:
    a = G.cartHexa((0,0,0),(0.1,0.1,0.1),(11,11,2))
    a = C.convertArray2NGon(a)
    a = G.close(a)
    a = C.fillEmptyBCWith(a, 'wall', 'BCWall', dim=3)
    I._adaptNGon32NGon4(a)
    C.convertPyTree2File(a, 'case.cgns')

Cmpi.barrier()

# Parallel ongoing
a, res = X.loadAndSplitNGon('case.cgns')
Cmpi.convertPyTree2File(a, 'parted.cgns')
gcells = res[5]
gfaces = res[6]
gpoints = res[7]
comm = res[1]

hmax = XOR.edgeLengthExtrema(a) / 1.0
hmin = XOR.edgeLengthExtrema(a) / 10.0
Tr = 0.15
Tu = 0.07
eps = 0.30
itermax1 = 10
itermax2 = 0
unrefine = False
mode_2D = np.array([0.0, 0.0, 1.0])
#mode_2D = None

own, nei = X._prepareMeshForAdaptation(a)

AM = X.CreateAdaptMesh(a, own, nei, comm, Tr, Tu, eps, hmin, hmax, unrefine,
    mode_2D, gcells, gfaces, gpoints)

ncells = C.getNCells(a)

for it in range(itermax1):
    C._initVars(a, 'centers:F', Func, ['centers:CoordinateX', 'centers:CoordinateY', 'centers:CoordinateZ'])

    f = I.getNodeFromName(a, 'F')[1]

    REF = np.empty(len(f), dtype=np.int32)
    for i in range(len(f)):
        REF[i] = int(f[i])
    
    X._assignRefDataToAM(AM, REF)

    X.AdaptMesh(AM)

    m = X.ExtractLeafMesh(AM, conformize=0)

    z = I.createZoneNode('p'+str(Cmpi.rank), m)
    
    a = C.newPyTree(['Base', z])
 
Cmpi.convertPyTree2File(a, 'refined.cgns')
