# - convertFile2PyTree -
# - FSDM format -
import Converter.PyTree as C
import Post.PyTree as P
import Transform.PyTree as T
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

import numpy

LOCAL = test.getLocal()

# 1 - Create ME Zone ###
a = G.cartTetra((-0.05,0.45,0.05), (0.1,0.1,0.1), (5,5,5))
b = G.cartPyra((0.4,0.4,0.), (0.1,0.1,0.1), (5,5,5))
c = G.cartPenta((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
d = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.1), (5,5,5))

ffcount = 0
for z in [a,b,c,d]:
    zNG = C.convertArray2NGon(z)
    e = P.exteriorFaces(zNG)
    e = T.splitSharpEdges(e, 75.)
    res = T.breakElements(e)
    for r in res: 
        C._addBC2Zone(z, 'farME%d'%ffcount, 'BCFarfield', subzone=r)
        ffcount += 1

BCInfo = C.getBCs([a,b,c,d])
z = C.mergeConnectivity([a,b,c,d], None)
C._recoverBCs(z, BCInfo=BCInfo)
t = C.newPyTree(['Base', [z]])
b = Internal.getBases(t)[0]

# 2 - Get ME Info ###
nvert,maxElt,maxFace = 0,0,0

cxME = Internal.getNodeFromName(t, 'CoordinateX')[1]
cyME = Internal.getNodeFromName(t, 'CoordinateY')[1]
czME = Internal.getNodeFromName(t, 'CoordinateZ')[1]
nvert = len(cxME)

for e in Internal.getElementNodes(z):
    rangeLocal = Internal.getElementRange(z, e[0])
    minR, maxR = rangeLocal
    maxElt = max(maxElt, maxR)

for e in Internal.getElementBoundaryNodes(z):
    rangeLocal = Internal.getElementRange(z, e[0])
    minR, maxR = rangeLocal
    maxFace = max(maxFace, maxR)

maxFace = max(maxElt, maxFace) # if no ME BC

# 3 - Create NGon Zone ###
zNG = G.cartNGon((0.85,0.,0.), (0.1,0.1,0.1), (5,5,5), api=3)

# add NG BCs
p = P.exteriorFaces(zNG)
C._addBC2Zone(zNG, 'farNG', 'BCFarfield', subzone=p)

# get GridCoordinates
cxNG = Internal.getNodeFromName(zNG, 'CoordinateX')[1]
cyNG = Internal.getNodeFromName(zNG, 'CoordinateY')[1]
czNG = Internal.getNodeFromName(zNG, 'CoordinateZ')[1]

# update NGonElements
e = Internal.getNodeFromName(zNG, 'NGonElements')
elementRange = Internal.getNodeFromName(e, 'ElementRange')[1]
elementConnectivity = Internal.getNodeFromName(e, 'ElementConnectivity')[1]
elementRange += maxFace
elementConnectivity += nvert

# updateNFaceElements
e = Internal.getNodeFromName(zNG, 'NFaceElements')
elementRange = Internal.getNodeFromName(e, 'ElementRange')[1]
elementRange += maxFace

# update BC
e = Internal.getNodeFromName(zNG, 'farNG')
pointList = Internal.getNodeFromName(e, 'PointList')[1]
pointList += maxFace

# 4 - merge both

# update zoneDim
z[1][0,0] += zNG[1][0,0]
z[1][0,1] += zNG[1][0,1]

# update GridCoordinates
cx = numpy.concatenate([cxME, cxNG])
cy = numpy.concatenate([cyME, cyNG])
cz = numpy.concatenate([czME, czNG])
Internal.getNodeFromName(z, 'CoordinateX')[1] = cx
Internal.getNodeFromName(z, 'CoordinateY')[1] = cy
Internal.getNodeFromName(z, 'CoordinateZ')[1] = cz

# add NGonElements, NFaceElements to zone
z[2] += [Internal.getNodeFromName(zNG, 'NGonElements'), Internal.getNodeFromName(zNG, 'NFaceElements')]

# add NGon BC to zone
ZBC = Internal.getNodeFromName(z, 'ZoneBC')
if ZBC is not None:
    ZBC[2] += Internal.getNodesByType(zNG, 'BC_t')
else:
    z[2] += [Internal.getNodeFromName(zNG, 'ZoneBC')]

# save and test out final tree (with all possible elements - ME + NG)
C.convertPyTree2File(t, LOCAL+'/out.h5')
t = C.convertFile2PyTree(LOCAL+'/out.h5')
test.testT(t, 1)