import Converter.PyTree as C
import Generator.PyTree as G
import XCore.PyTree as X
import Converter.Mpi as Cmpi
import Converter.Internal as I
import numpy as np
import Intersector.PyTree as XOR
import subprocess

a, res = X.loadAndSplitNGon('case_v4.cgns')
Cmpi.convertPyTree2File(a, 'parted.cgns')
gcells = res[5]
gfaces = res[6]
gpoints = res[7]
comm = res[1]

hmax = XOR.edgeLengthExtrema(a) / 1.0
hmin = XOR.edgeLengthExtrema(a) / 4.0
Tr = 0.05
Tu = 0.07
eps = 0.30
itermax = 5
unrefine = False
#mode_2D = np.array([0.0, 0.0, 1.0])
mode_2D = None

own, nei = X._prepareMeshForAdaptation(a)

AM = X.CreateAdaptMesh(a, own, nei, comm, Tr, Tu, eps, hmin, hmax, unrefine,
    mode_2D, gcells, gfaces, gpoints)

ncells = C.getNCells(a)

fc, fa = G.getFaceCentersAndAreas(a)
centers = G.getCellCenters(a, fc, fa, [own], [nei])[0]
cx = centers[0]; cy = centers[1]; cz = centers[2]

fsolc = I.getNodeFromName(a, I.__FlowSolutionCenters__)
f = I.getNodeFromName(fsolc, "Ma")[1]

g = X.computeGradient(AM, f, cx, cy, cz, own, nei)
h = X.computeHessian(AM, f, g, cx, cy, cz, own, nei)
REF = X._makeRefDataFromGradAndHess(AM, f, g, h)

X._assignRefDataToAM(AM, REF)

X.AdaptMesh(AM)

m, BCs, own, nei = X.ExtractLeafMesh(AM, conformize=1)
#m = X.ExtractLeafMesh(AM, conformize=0) # HEXA

z = I.createZoneNode('p'+str(Cmpi.rank), m)

for i in range(len(BCs)):
    cont = I.createUniqueChild(z, 'ZoneBC', 'ZoneBC_t')
    BC = BCs[i]
    bc = I.newBC(name=BC[1], pointList=BC[0], family=BC[1], parent=cont)
    I._createUniqueChild(bc, 'GridLocation', 'GridLocation_t', value='FaceCenter')

a = C.newPyTree(['Base', z])

# Create QUAD bczones
CX = I.getNodeFromName(z, 'CoordinateX')[1]
CY = I.getNodeFromName(z, 'CoordinateY')[1]
CZ = I.getNodeFromName(z, 'CoordinateZ')[1]

bcs = X.extractBoundaryMesh(AM, mode=0)

for bc in bcs:
    A = ['CoordinateX,CoordinateY,CoordinateZ', [CX,CY,CZ], [bc[0]], 'QUAD']
    zone = I.createZoneNode(bc[1].upper(), A)
    a[2][1][2].append(zone)

Cmpi.convertPyTree2File(a, 'refined.cgns')
