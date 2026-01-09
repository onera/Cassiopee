# - convertFile2PyTree -
# - tau grid format -
import Converter.PyTree as C
import Post.PyTree as P
import Geom.PyTree as D
import Transform.PyTree as T
import Generator.PyTree as G
import KCore.test as test

LOCAL = test.getLocal()

dz = 0.01
xmin, ymin, zmin, xmax, ymax, zmax = [-0.5,-0.5,0,1.5,0.5,dz]
size = 0.01
L = 1
N = 200
N_airfoil = int(L/size*2)+1

tbox = G.cartHexa((xmin,ymin,zmin),(xmax-xmin,ymax-ymin,zmin),(2,2,1))
ff = P.exteriorFaces(tbox)
ff = D.uniformize(ff, N)

airfoil = D.naca(12., N=N_airfoil)

airfoil = T.reorder(airfoil, (-1,2,3))
distrib = D.line((0,0,0), (0.1,0,0), N=10)
bl = G.addNormalLayers(airfoil, distrib, niter=100)
bl = T.reorder(bl, (-1,2,3))
bl = C.convertArray2Hexa(bl)
bl = G.close(bl)

ext = P.exteriorFaces(bl)
ext = T.splitConnexity(ext)
ext1 = T.reorder(ext[1], (-1,))
borders = T.join(ext1, ff)

m2 = G.T3mesher2D(borders, triangulateOnly=0, grading=1.1, metricInterpType=0)
m2 = T.reorder(m2, (1,))

m = [bl, m2]

T._addkplane(m)
T._contract(m, (0,0,0), (1,0,0), (0,1,0), dz)

# Add BCs
ext = []
for i in m:
    e = P.exteriorFaces(i)
    ext += T.splitSharpEdges(e, 30.) # because of prisms, ext is in ngon

m = C.mergeConnectivity(bl, m2)

for c, e  in enumerate(ext):
    res = T.breakElements(e)
    if c == 0: #1
        for r in res: C._addBC2Zone(m, 'up', 'BCSymmetryPlane', subzone=r) # sym 2D
    elif c == 1: #2
        for r in res: C._addBC2Zone(m, 'wall', 'BCWall', subzone=r)
    elif c == 2: #3
        for r in res: C._addBC2Zone(m, 'down', 'BCSymmetryPlane', subzone=r) # sym 2D
    elif c == 4: #4
        for r in res: C._addBC2Zone(m, 'down2', 'BCSymmetryPlane', subzone=r) # sym 2D
    elif c == 5: #5
        for r in res: C._addBC2Zone(m, 'up2', 'BCSymmetryPlane', subzone=r) # sym 2D
    elif c == 6: #6
        for r in res: C._addBC2Zone(m, 'far1', 'BCFarfield', subzone=r)
    elif c == 7: #7
        for r in res: C._addBC2Zone(m, 'far2', 'BCFarfield', subzone=r)
    elif c == 8: #8
        for r in res: C._addBC2Zone(m, 'far3', 'BCFarfield', subzone=r)
    elif c == 9: #9
        for r in res: C._addBC2Zone(m, 'far4', 'BCFarfield', subzone=r)

#C.convertPyTree2File(m, 'mesh.cgns')

import KCore.Dist as Dist
from KCore.config import *
(netcdf, netcdfIncDir, netcdfLibDir, netcdflibs) = Dist.checkNetcdf(additionalLibPaths,
                                                                    additionalIncludePaths)

if netcdf:
    C.convertPyTree2File(m, LOCAL+'/out.grid')
    t = C.convertFile2PyTree(LOCAL+'/out.grid')
    test.testT(t, 1)
