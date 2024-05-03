# - octree -
# octree: low height
# Adaptation du maillage avec des zones definies dans une base 'Refine'
# Mesh adaptation in zones in the 'Refine' base
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import Connector.PyTree as X
import Post.PyTree as P
import Transform.PyTree as T
import Dist2Walls.PyTree
import Converter.Internal as Internal
import numpy, math

# INPUT: Geometrie fermee, maillee correctement avec un pas a peu pres 
# homogene (TRI + normales exterieures!)

# INPUT: closed mesh with approximatively homogenic mesh step. The mesh should be TRI and the normals directed towards exterior
s = C.convertFile2PyTree('catR.cgns')

# taille des mailles de l'octree autour de la surface (si -1 determine
# a partir de la taille des triangles de la surface)

# Octree mesh size around the surface (if -1 : computed from surface TRI mesh size) 

snear = -1

# Taille de l'evidement en *snear (si -1 calcule automatiquement)
# Hole size (*snear) (if -1 : automatically computed)
holeFactor = -1 # Hole size (*snear)
dfar = -1 # distance des frontieres lointaine du corps , boundary distance far from body
hWall = 0.001 # first cell height ( if prism, 0. otherwise)
OUTPUT = 0 # 0: ELSA, 1: CEDRE

# Get Refine Base (if any)
refine = Internal.getNodeFromName1(s, 'Refine')
if refine is not None:
    s = Internal.rmNodesByName(s, 'Refine')

# - octree4 firewall - 
s = C.deleteFlowSolutions__(s)
s = C.convertArray2Tetra(s)
s = T.join(s); s = G.close(s)
s = T.reorder(s, (+1,))
# - end of firewall -

XMax_surf = C.getMaxValue(s,'CoordinateX')

# snear computation (if snear < 0)
sizemax = 0.
box = G.bbox(s)
sizemax = max(sizemax,abs(box[3]-box[0]),abs(box[4]-box[1]),abs(box[5]-box[2]));print('sizemax=',sizemax)
if snear < 0:
    mesh = G.getNormalMap(s)
    mesh = C.magnitude(mesh, ['centers:sx','centers:sy','centers:sz'])
    meanarea = C.getMeanValue(mesh,'centers:sMagnitude')
    minarea = C.getMinValue(mesh,'centers:sMagnitude')
    maxarea = C.getMaxValue(mesh,'centers:sMagnitude')
    print('surface area min/max/mean=',minarea,maxarea,meanarea)
    snear = 2.*meanarea**(1./2.); print('snear=', snear)
if dfar < 0: dfar = 10.*sizemax; print('dfar=', dfar)
 
# hole height
if holeFactor < 0: 
    if (hWall <= 0.): holeFactor = 2.
    else: holeFactor = 2.5
    delta = holeFactor*0.5*snear
else: delta = holeFactor*0.5*snear
print('delta=',delta)

xrayDim1 = 1500; xrayDim2 = 1500

# - Generation of the Boundary layer (if any)
if hWall > 0.:
    print('Generating boundary layers...')
    algoBL = 0; secondlayer = True
    smoothIter = 100
    if algoBL == 0: 
        raison = 1.2 # BL cell height factor = 1.2 
        nlayer = int(1+math.log(snear/(8.*hWall))/math.log(raison)); print('number of layers = ',nlayer)
        dlayer = hWall*(1.-raison**nlayer)/(1.-raison); print('first layer height = ',dlayer)
        d = G.cart((0.,0.,0.), (dlayer/nlayer,1,1), (nlayer,1,1))
        for i in range(0,nlayer): C.setValue(d, 'CoordinateX', i, hWall*(1.-raison**i)/(1.-raison))
        hLast = C.getValue(d, 'CoordinateX', nlayer-1)-C.getValue(d, 'CoordinateX', nlayer-2); print('last layer height = ',hLast)
        lay = G.addNormalLayers(s, d, check=0, niter=smoothIter)
        if secondlayer:
            # second layer with height factor = 2 
            ext = P.exteriorFaces(lay)
            ext = T.splitConnexity(ext)
            s = ext[1]
            s = T.breakElements(s)[0]
            s = T.reorderAll(s)
            hWall = hLast
            raison = 2.
            nlayer = int(1+math.log(snear/(1.*hWall))/math.log(raison)); print('nb de couches de prismes = ',nlayer)
            dlayer = hWall*(1-raison**nlayer)/(1-raison); print('epaisseur de la couche de prismes = ',dlayer)
            d = G.cart((0.,0.,0.), (dlayer/nlayer,1,1), (nlayer,1,1))
            for i in range(0,nlayer): C.setValue(d, 'CoordinateX', i, hWall*(1-raison**i)/(1-raison))
            hLast = C.getValue(d, 'CoordinateX', nlayer-1)-C.getValue(d, 'CoordinateX', nlayer-2); print('hauteur de la derniere couche de prisme = ',hLast)
            lay = C.convertArray2NGon(lay)
            lay1 = G.addNormalLayers(s, d, check=0, niter=smoothIter)
            lay = T.join(lay, lay1)
    else: # distribution with G.enforce
        nlayer = 41
        Nr = int(0.3*nlayer)
        Nr2 = int(1.*nlayer)
        d = G.cart((0,0,0),(delta/nlayer,1,1),(nlayer,1,1))
        lend = D.getLength(d)
        d = D.getCurvilinearAbscissa(d)
        distrib = C.cpVars(d,'s',d,'CoordinateX')
        C._rmVars(distrib,'s')
        C._rmVars(d,'s')
        distrib = G.enforceMoinsX(distrib,0.5*snear/lend, Nr2, Nr2)
        distrib = G.enforcePlusX(distrib,hWall/lend, Nr2, Nr2)
        # distrib = G.enforceMoinsX(distrib,0.5*snear/lend,Nr2, Nr2)
        # print(distrib)
        d = G.map(d,distrib,dir=1)
        lay = G.addNormalLayers(s, d, check=0, niter=smoothIter)
    C.convertPyTree2File(lay, 'layer.cgns')
    ext = P.exteriorFaces(lay)
    ext = T.splitConnexity(ext)
    s = ext[1]
    s = T.breakElements(s)[0]
    s = T.reorderAll(s)
    C.convertPyTree2File(s, 'externalLayer.cgns')

    # Recalcule le snear base sur la surface prismatique externe
    # Recompute snear factor based on the exterior BL layer

    if snear < 0.: # snear auto
        mesh = G.getNormalMap(s)
        mesh = C.magnitude(mesh, ['centers:sx','centers:sy','centers:sz'])
        meanarea = C.getMeanValue(mesh,'centers:sMagnitude')
        minarea = C.getMinValue(mesh,'centers:sMagnitude')
        maxarea = C.getMaxValue(mesh,'centers:sMagnitude')
        print('surface area min/max/mean=',minarea,maxarea,meanarea)
        snear = 2.*meanarea**(1./2.); print('snear=', snear)
        # hole height (upward) 
        delta = holeFactor*0.5*snear
        print('delta=',delta)

# - Octree generation
print('Generating octree...')
surfaces = [s]; snears = [snear]
o = G.octree(surfaces, snears, dfar=dfar, balancing=1)

bbox = G.bbox(o)
expand = int(delta/snear)+2
print('Expanding layers %d...'%expand)
for i in range(expand):
    o = G.expandLayer(o, level=0, balancing=1)
t = C.newPyTree(['Base']); t[2][1][2] += [o]

# - Adapt octree
if refine is not None:
    print('Adapting octree...')
    volmin = (2*snear)**3 # 2nd level
    # refinement surface
    ref = C.convertArray2Tetra(refine)
    ref = T.join(ref); ref = G.close(ref)
    BM = numpy.array([[1]])
    end = 0; count = 0
    while end == 0 and count < 4:
        t = X.blankCells(t, [[ref]], BM, blankingType='center_in')
        t = G.getVolumeMap(t)
        t = C.initVars(t, '{centers:indicator} = logical_and({centers:cellN}<0.001, {centers:vol}>%f)'%volmin)
        end = 1
        if C.getMaxValue(t,'centers:indicator') == 1.: end = 0
        # Keeping the lowest refinement level
        o = t[2][1][2][0]
        o = G.adaptOctree(o, 'centers:indicator', balancing=1)
        # Output
        t = C.newPyTree(['Base']); t[2][1][2] += [o]
        count += 1
C.convertPyTree2File(t, 'octree.cgns')

# - Extract blanking surface
print('Extract blanking surface (%s)...'%delta)
o = Dist2Walls.PyTree.distance2Walls(o, s, type='ortho', 
                                     loc='centers', signed=1)
s2 = P.isoSurfMC(o, 'centers:TurbulentDistance', delta)
C.convertPyTree2File(s2+[s], 'deltaSurface.cgns')
t = C.newPyTree(['Base']); t[2][1][2] += [o]

# Blanking with delta surface
print('Blanking octree width delta...')
BM = numpy.array([[1]])
t = X.blankCells(t, [s2], BM, blankingType='node_in', delta=1.e-12,
                 XRaydim1=xrayDim1, XRaydim2=xrayDim2)

print('Conformizing...')
t = C.convertArray2NGon(t); t = G.close(t)
t = C.conformizeNGon(t); t = G.close(t)

# Select cells with cellN=1
print('Selecting...')
t = P.selectCells2(t, 'cellN', strict=1)

# Getting exterior faces
print('Getting exterior faces...')
ext = P.exteriorFaces(t)
ext = T.splitConnexity(ext)
C.convertPyTree2File(ext, 'ext2Watch.cgns')

# Get the body exterior surface 
bbox = [bbox[0]+1.e-6,bbox[1]+1.e-6,bbox[2]+1.e-6,
        bbox[3]-1.e-6,bbox[4]-1.e-6,bbox[5]-1.e-6]
i = 0; lgt = len(ext)
for e in ext:
    bb = G.bbox(e)
    if (bb[0] < bbox[0] or bb[3] > bbox[3] or bb[1] < bbox[1] or bb[4] > bbox[4] or bb[2] < bbox[2] or bb[5] > bbox[5]): pass
    else: ext = e; break
    i += 1
if (i >= lgt): 
    print('Exterior not found.')
    import sys; sys.exit()

ext = T.reorder(ext, (+1,))
ext = C.rmVars(ext, 'cellN')

# Tetra filling between exterior faces and surface
print('Tetra filling...')
ext = C.convertArray2Tetra(ext)
ext = G.close(ext) # force
ext = T.reorder(ext, (-1,)) # warn
s = T.reorder(s, (-1,))
m = T.join(ext, s)
C.convertPyTree2File([m], 'ext.cgns')

TETRAMESHER='tetgen'
if (TETRAMESHER == 'tetgen'): algo = 1
else: algo = 0

m = G.tetraMesher(m, maxh=0.5*snear, algo=1)
C.convertPyTree2File(m, 'tetra1.cgns')

# Reconnect
print('Reconnecting...')
z = t[2][1][2][0]
z = C.rmVars(z, 'cellN')
m = G.close(m)
m = C.convertArray2NGon(m)

tp = C.newPyTree(['Base']); tp[2][1][2] += [m,z]

m = T.join(m, z)
m = G.close(m) 
m = C.conformizeNGon(m, tol=0.00001*snear)

if hWall > 0.:
    lay = C.convertArray2NGon(lay)
    m = T.join(m, lay)

print('Final closing...')
print(Internal.getZoneDim(m))
m = G.close(m)

# Create BCs
print('Creating BCs...')
ext = P.exteriorFaces(m)
(xm,ym,zm) = C.getMinValue(ext, 'GridCoordinates')
(xp,yp,zp) = C.getMaxValue(ext, 'GridCoordinates')

bc = C.newPyTree(['Base',2])
ext = T.splitSharpEdges(ext, alphaRef=90.)
for zone in Internal.getNodesFromType2(ext, 'Zone_t'):
    if C.getMaxValue(zone,'CoordinateX') != XMax_surf: bc[2][1][2].append(zone) # Delete body surface
bc2 = T.splitSharpEdges(bc, alphaRef=30.)

for zone in Internal.getNodesFromType2(bc, 'Zone_t'):
    (xpZone, ypZone) = C.getMaxValue(zone, ['CoordinateX', 'CoordinateY']) 
    if (xpZone == xp and ypZone == yp): Internal._rmNode(bc,zone)
for zone in Internal.getNodesFromType2(bc2, 'Zone_t'):
    bc[2][1][2].append(zone)
C.convertPyTree2File(bc, 'exterior.cgns')

# Add BCs
eps = 1.e-12
for zone in Internal.getByType(bc,'Zone_t')[2]:
    dimax = Internal.getZoneDim(zone)[1] - 1 ; print(zone[0],dimax,G.barycenter(zone))
    bary = G.barycenter(zone)
    if abs(bary[0] - xm) < eps : C._addBC2Zone(m,'OUT','FamilySpecified:OUT',subzone=zone); print('Xmin')
    if abs(bary[0] - xp) < eps : C._addBC2Zone(m,'IN','FamilySpecified:IN',subzone=zone); print('Xmax')
    if abs(bary[1] - ym) < eps : C._addBC2Zone(m,'RIGHT','FamilySpecified:RIGHT',subzone=zone); print('Ymin')
    if abs(bary[1] - yp) < eps : C._addBC2Zone(m,'LEFT','FamilySpecified:LEFT',subzone=zone); print('Ymax')
    if abs(bary[2] - zm) < eps : C._addBC2Zone(m,'DOWN','FamilySpecified:DOWN',subzone=zone); print('Zmin')
    if abs(bary[2] - zp) < eps : C._addBC2Zone(m,'TOP','FamilySpecified:TOP',subzone=zone); print('Zmax')

m = C.fillEmptyBCWith(m, 'WALL', 'FamilySpecified:WALL')

tp = C.newPyTree(['Base']); tp[2][1][2] += [m]

C._addFamily2Base(tp, 'TOP', bndType='BCFarfield')
C._addFamily2Base(tp, 'DOWN', bndType='BCFarfield') 
C._addFamily2Base(tp, 'IN', bndType='BCFarfield') 
C._addFamily2Base(tp, 'OUT', bndType='BCFarfield')
C._addFamily2Base(tp, 'LEFT', bndType='BCFarfield') 
C._addFamily2Base(tp, 'RIGHT', bndType='BCFarfield') 
C._addFamily2Base(tp, 'WALL', bndType='BCWall') 

C.convertPyTree2File(tp, 'catOut.cgns')

if OUTPUT == 0: # ELSA
    Internal._createElsaHybrid(tp)
    C.convertPyTree2File(tp, 'catOut2.hdf')
if OUTPUT == 1: # CEDRE
    C.convertPyTree2File(tp, 'catOut2.d')
