import Converter.PyTree as C
import Generator.PyTree as G
import XCore.PyTree as X
import Converter.Mpi as Cmpi
import Converter.Internal as I
import numpy as np
import Intersector.PyTree as XOR
import subprocess

# mpirun -np <N> python3 <script.py>

def F1(x, y, z):
    return np.tanh(-100.*(y-0.5-0.25*np.sin(2.*np.pi*x))) + np.tanh(100.*(y-x))

def F2(x, y, z):
    e = 0.25
    u = 4.*x-2.
    v = 4.*y-2.
    return np.tanh(30.*(u*u + v*v - e)) + \
           np.tanh(30.*((u-0.75)**2 + (v-0.75)**2 - e)) + \
           np.tanh(30.*((u-0.75)**2 + (v+0.75)**2 - e)) + \
           np.tanh(30.*((u+0.75)**2 + (v-0.75)**2 - e)) + \
           np.tanh(30.*((u+0.75)**2 + (v+0.75)**2 - e))

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
itermax = 6
unrefine = False
mode_2D = np.array([0.0, 0.0, 1.0])
#mode_2D = None

own, nei = X._prepareMeshForAdaptation(a)

AM = X.CreateAdaptMesh(a, own, nei, comm, Tr, Tu, eps, hmin, hmax, unrefine,
    mode_2D, gcells, gfaces, gpoints)

ncells = C.getNCells(a)

for it in range(itermax):

    fc, fa = G.getFaceCentersAndAreas(a)
    centers = G.getCellCenters(a, fc, fa, [own], [nei])[0]
    cx = centers[0]; cy = centers[1]; cz = centers[2]
    f = make_field(cx, cy, cz, F2)
    g = X.computeGradient(AM, f, cx, cy, cz, own, nei)
    h = X.computeHessian(AM, f, g, cx, cy, cz, own, nei)
    REF = X._makeRefDataFromGradAndHess(AM, f, g, h)

    X._assignRefDataToAM(AM, REF)

    X.AdaptMesh(AM)

    m, BCs, own, nei = X.ExtractLeafMesh(AM, conformize=1)

    z = I.createZoneNode('p'+str(Cmpi.rank), m)
    
    for i in range(len(BCs)):
        cont = I.createUniqueChild(z, 'ZoneBC', 'ZoneBC_t')
        BC = BCs[i]
        bc = I.newBC(name=BC[1], pointList=BC[0], family=BC[1], parent=cont)
        I._createUniqueChild(bc, 'GridLocation', 'GridLocation_t', value='FaceCenter')

    a = C.newPyTree(['Base', z])

Cmpi.convertPyTree2File(a, 'refined.cgns')
