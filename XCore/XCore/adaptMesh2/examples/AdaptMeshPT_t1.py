import Converter.PyTree as C
import Generator.PyTree as G
import XCore.PyTree as X
import Converter.Internal as I
import numpy as np
import Intersector.PyTree as XOR
import subprocess

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

a = G.cartHexa((0,0,0),(0.02,0.02,0.01),(51,51,2))
a = C.convertArray2NGon(a)
a = G.close(a)
a = C.fillEmptyBCWith(a, 'wall', 'BCWall', dim=3)
I._adaptNGon32NGon4(a)

# Placeholders for now
hmax = XOR.edgeLengthExtrema(a)
hmin = hmax / 4.0
eps = 1e-1

# Refinement threshold
Tr = 0.1
# Unrefinement threshold
Tu = Tr/2.5
# Max iterations for first adaptation run
itermax1 = 2
# Max iterations for second adaptation run
itermax2 = 3

# Next two functions must be called before adaptation cycles
own, nei = X._prepareMeshForAdaptation(a)

AM = X.CreateAdaptMesh(a, own, nei, Tr, Tu)

ncells = C.getNCells(a)

for it in range(itermax1):
    print('\n*******************************')
    print('Iteration', it) 
   
    print("Getting face centers and areas")
    fc, fa = G.getFaceCentersAndAreas(a)
    
    print("Getting cell centers")
    centers = G.getCellCenters(a, fc, fa, [own], [nei])[0]

    cx = centers[0]; cy = centers[1]; cz = centers[2]
   
    print("Making sensor")
    field = make_field(cx, cy, cz, F1)

    print("Computing hessian")
    H = X.computeHessian(a, field, cx, cy, cz, own, nei)

    print("Making metric")
    M = X.hessianToMetric(H, hmin, hmax, eps)

    print("Making refinement data")
    X._metricToRefData(a, M, AM)

    print("Adapting")

    m, BCs, own, nei = X.AdaptMesh(AM, fc[0], cx, cy, cz)

    z = I.createZoneNode('iter%d'%it, m)

    for i in range(len(BCs)):
        cont = I.createUniqueChild(z, 'ZoneBC', 'ZoneBC_t')
        BC = BCs[i]
        bc = I.newBC(name=BC[1], pointList=BC[0], family=BC[1], parent=cont)
        I._createUniqueChild(bc, 'GridLocation', 'GridLocation_t', value='FaceCenter')
    
    a = C.newPyTree(['Base', z])

    C.convertPyTree2File(a, 'refined1.cgns')
    
    new_ncells = C.getNCells(a)
    if ncells == new_ncells: break
    ncells = new_ncells


for it in range(itermax2):
    print('\n*******************************')
    print('Iteration', it) 
   
    print("Getting face centers and areas")
    fc, fa = G.getFaceCentersAndAreas(a)
    
    print("Getting cell centers")
    centers = G.getCellCenters(a, fc, fa, [own], [nei])[0]

    cx = centers[0]; cy = centers[1]; cz = centers[2]
   
    print("Making sensor")
    field = make_field(cx, cy, cz, F2)

    print("Computing hessian")
    H = X.computeHessian(a, field, cx, cy, cz, own, nei)

    print("Making metric")
    M = X.hessianToMetric(H, hmin, hmax, eps)

    print("Making refinement data")
    X._metricToRefData(a, M, AM)

    print("Adapting")

    m, BCs, own, nei = X.AdaptMesh(AM, fc[0], cx, cy, cz)

    z = I.createZoneNode('iter%d'%it, m)

    for i in range(len(BCs)):
        cont = I.createUniqueChild(z, 'ZoneBC', 'ZoneBC_t')
        BC = BCs[i]
        bc = I.newBC(name=BC[1], pointList=BC[0], family=BC[1], parent=cont)
        I._createUniqueChild(bc, 'GridLocation', 'GridLocation_t', value='FaceCenter')
    
    a = C.newPyTree(['Base', z])

    C.convertPyTree2File(a, 'refined2.cgns')
    
    new_ncells = C.getNCells(a)
    if ncells == new_ncells: break
    ncells = new_ncells
