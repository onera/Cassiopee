# - octree -
# octree: evidemment de faible hauteur
# Adaptation du maillage avec des zones definies dans une base 'Refine'
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import Connector.PyTree as X
import Post.PyTree as P
import Transform.PyTree as T
import Dist2Walls.PyTree
import Converter.Internal as Internal
import numpy, math
import cases

# INPUT: Geometrie fermee, maillee correctement avec un pas a peu pres 
# homogene (TRI + normales exterieures!)
s = C.convertFile2PyTree('catR.cgns')

# taille des mailles de l'octree autour de la surface (si -1 determine
# a partir de la taille des triangles de la surface)
snear = -1
# Taille de l'evidement en *snear (si -1 calcule automatiquement)
holeFactor = -1 # evidement en * snear
dfar = cases.CASES[case][3] # distance des frontieres lointaine du corps
hWall = cases.CASES[case][4] # hauteur de la premiere maille (si prismes)

# Force auto
snear = -1; dfar = -1; holeFactor = -1

# INPUT: Geometrie fermee, maillee correctement avec un pas homogene snear 
# (TRI + normales exterieures!)
s = C.convertFile2PyTree(case+'R.cgns')

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

# Calcul automatique du snear (si snear < 0)
sizemax = 0.
box = G.bbox(s)
sizemax = max(sizemax,abs(box[3]-box[0]),abs(box[4]-box[1]),abs(box[5]-box[2]));print 'sizemax=',sizemax
if (snear < 0):
    mesh = G.getNormalMap(s)
    mesh = C.magnitude(mesh, ['centers:sx','centers:sy','centers:sz'])
    meanarea = C.getMeanValue(mesh,'centers:sMagnitude')
    minarea = C.getMinValue(mesh,'centers:sMagnitude')
    maxarea = C.getMaxValue(mesh,'centers:sMagnitude')
    print 'surface area min/max/mean=',minarea,maxarea,meanarea
    # snear.append(2.*(vmin-1)*meanarea**(1/2.))
    snear = 2.*meanarea**(1./2.); print 'snear=', snear
if (dfar < 0): dfar = 10.*sizemax; print 'dfar=', dfar
 
# Hauteur evidement
if (holeFactor < 0): 
    if (hWall <= 0.): holeFactor = 2.
    else: holeFactor = 2.5
    delta = holeFactor*0.5*snear
else: delta = holeFactor*0.5*snear
print 'delta=',delta

xrayDim1 = 1500; xrayDim2 = 1500

# - Genere la Boundary layer (if any)
if (hWall > 0.):
    print 'Generating boundary layers...'
    algoBL = 0; secondlayer = True
    smoothIter = 100
    if algoBL == 0: # raison geometrique 
        raison = 1.2
        nlayer = 1+math.log(snear/(8.*hWall))/math.log(raison); print 'nb de couches de prismes = ',nlayer
        nlayer = int(nlayer)
        dlayer = hWall*(1.-raison**nlayer)/(1.-raison); print 'epaisseur de la couche de prismes = ',dlayer
        d = G.cart((0.,0.,0.), (dlayer/nlayer,1,1), (nlayer,1,1))
        for i in xrange(0,nlayer): C.setValue(d, 'CoordinateX', i, hWall*(1.-raison**i)/(1.-raison))
        hLast = C.getValue(d, 'CoordinateX', nlayer-1)-C.getValue(d, 'CoordinateX', nlayer-2); print 'hauteur de la derniere couche de prisme = ',hLast
        lay = G.addNormalLayers(s, d, check=0, niter=smoothIter)
        if secondlayer:
            # 2eme couche
            ext = P.exteriorFaces(lay)
            ext = T.splitConnexity(ext)
            s = ext[1]
            s = T.breakElements(s)[0]
            s = T.reorderAll(s)
            hWall = hLast
            raison = 2.
            nlayer = 1+math.log(snear/(1.*hWall))/math.log(raison); print 'nb de couches de prismes = ',nlayer
            nlayer = int(nlayer)
            dlayer = hWall*(1-raison**nlayer)/(1-raison); print 'epaisseur de la couche de prismes = ',dlayer
            d = G.cart((0.,0.,0.), (dlayer/nlayer,1,1), (nlayer,1,1))
            for i in xrange(0,nlayer): C.setValue(d, 'CoordinateX', i, hWall*(1-raison**i)/(1-raison))
            hLast = C.getValue(d, 'CoordinateX', nlayer-1)-C.getValue(d, 'CoordinateX', nlayer-2); print 'hauteur de la derniere couche de prisme = ',hLast
            lay = C.convertArray2NGon(lay)
            lay1 = G.addNormalLayers(s, d, check=0, niter=smoothIter)
            lay = T.join(lay, lay1)
    else:
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
        # print distrib
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
    if cases.CASES[case][0] < 0.: # snear auto
        mesh = G.getNormalMap(s)
        mesh = C.magnitude(mesh, ['centers:sx','centers:sy','centers:sz'])
        meanarea = C.getMeanValue(mesh,'centers:sMagnitude')
        minarea = C.getMinValue(mesh,'centers:sMagnitude')
        maxarea = C.getMaxValue(mesh,'centers:sMagnitude')
        print 'surface area min/max/mean=',minarea,maxarea,meanarea
        # snear.append(2.*(vmin-1)*meanarea**(1/2.))
        snear = 2.*meanarea**(1./2.); print 'snear=', snear
        # Hauteur evidement (au dessus) 
        delta = holeFactor*0.5*snear
        print 'delta=',delta

# - Genere Octree
print 'Generating octree...'
surfaces = [s]; snears = [snear]
o = G.octree(surfaces, snears, dfar=dfar, balancing=1)

bbox = G.bbox(o)
expand = int(delta/snear)+2
print 'Expanding layers %d...'%expand
for i in xrange(expand):
    o = G.expandLayer(o, level=0, balancing=1)
t = C.newPyTree(['Base']); t[2][1][2] += [o]

# - Adapt octree
if refine is not None:
    print 'Adapting octree...'
    volmin = (2*snear)**3 # niveau 2
    # refinement surface
    ref = C.convertArray2Tetra(refine)
    ref = T.join(ref); ref = G.close(ref)
    BM = numpy.array([[1]])
    end = 0; count = 0
    while end == 0 and count < 4:
        #print count, end
        t = X.blankCells(t, [[ref]], BM, blankingType='center_in')
        #t = X.blankCellsTri(t, [[ref]], BM, blankingType='center_in')
        t = G.getVolumeMap(t)
        t = C.initVars(t, '{centers:indicator} = logical_and({centers:cellN}<0.001, {centers:vol}>%f)'%volmin)
        end = 1
        if C.getMaxValue(t,'centers:indicator') == 1.: end = 0
        # Maintien du niveau de raffinement le plus fin
        o = t[2][1][2][0]
        o = G.adaptOctree(o, 'centers:indicator', balancing=1)
        # Sortie
        t = C.newPyTree(['Base']); t[2][1][2] += [o]
        count += 1
C.convertPyTree2File(t, 'octree.cgns')

# - Extract blanking surface
print 'Extract blanking surface (%s)...'%delta
o = Dist2Walls.PyTree.distance2Walls(o, s, type='ortho', 
                                     loc='centers', signed=1)
s2 = P.isoSurfMC(o, 'centers:TurbulentDistance', delta)
C.convertPyTree2File(s2+[s], 'deltaSurface.cgns')
t = C.newPyTree(['Base']); t[2][1][2] += [o]

# Blanking with delta surface
print 'Blanking octree width delta...'
BM = numpy.array([[1]])
t = X.blankCells(t, [s2], BM, blankingType='node_in', delta=1.e-12,
                 XRaydim1=xrayDim1, XRaydim2=xrayDim2)
#s2 = C.convertArray2Tetra(s2)
#t = X.blankCellsTri(t, [s2], BM, blankingType='node_in')

print 'Conformizing...'
t = C.convertArray2NGon(t); t = G.close(t)
t = C.conformizeNGon(t); t = G.close(t)
#C.convertPyTree2File(t, 'blanked.cgns')

# Selectionne les cellules cellN=1
print 'Selecting...'
t = P.selectCells2(t, 'cellN', strict=1)
#C.convertPyTree2File(t, 'main.cgns')

# Recupere l'exterieur faces
print 'Getting exterior faces...'
ext = P.exteriorFaces(t)
ext = T.splitConnexity(ext)
C.convertPyTree2File(ext, 'ext2Watch.cgns')

# Selectionne la bonne frontiere (celle autour du corps)
bbox = [bbox[0]+1.e-6,bbox[1]+1.e-6,bbox[2]+1.e-6,
        bbox[3]-1.e-6,bbox[4]-1.e-6,bbox[5]-1.e-6]
i = 0; lgt = len(ext)
for e in ext:
    bb = G.bbox(e)
    if (bb[0] < bbox[0] or bb[3] > bbox[3] or bb[1] < bbox[1] or bb[4] > bbox[4] or bb[2] < bbox[2] or bb[5] > bbox[5]): pass
    else: ext = e; break
    i += 1
if (i >= lgt): 
    print 'Exterior not found.'
    import sys; sys.exit()

ext = T.reorder(ext, (+1,))
#res = P.integNormProduct(ext,['CoordinateX','CoordinateY','CoordinateZ'])
#if (res < 0): ext = T.reorder(ext, (-1,))
ext = C.rmVars(ext, 'cellN')

# Rempli l'espace entre l'exterieur faces et la surface
print 'Tetra filling...'
ext = C.convertArray2Tetra(ext)
ext = G.close(ext) # force
ext = T.reorder(ext, (-1,)) # warn
s = T.reorder(s, (-1,))
m = T.join(ext, s)
C.convertPyTree2File([m], 'ext.cgns')

# Si relecture
#t = C.convertFile2PyTree('main.cgns')
#ext = C.convertFile2PyTree('ext.cgns')
#m = ext[2][1][2][0]

TETRAMESHER='tetgen'
if (TETRAMESHER == 'tetgen'): algo = 1
else: algo = 0

m = G.tetraMesher(m, maxh=0.5*snear, algo=1)
C.convertPyTree2File(m, 'tetra1.cgns')

# Reconnect
print 'Reconnecting...'
z = t[2][1][2][0]
z = C.rmVars(z, 'cellN')
m = G.close(m)
m = C.convertArray2NGon(m)

tp = C.newPyTree(['Base']); tp[2][1][2] += [m,z]
#C.convertPyTree2File(tp, 'twoMeshs.cgns')

m = T.join(m, z)
m = G.close(m) # a forcer eventuellement
m = C.conformizeNGon(m, tol=0.00001*snear)

if hWall > 0.:
    lay = C.convertArray2NGon(lay)
    m = T.join(m, lay)

print 'Final closing...'
print Internal.getZoneDim(m)
m = G.close(m)

C.convertPyTree2File(m, case+'Out.cgns')

# Create BCs
print 'Creating BCs...'
ext = P.exteriorFaces(m)
(xm,ym,zm) = C.getMinValue(ext, 'GridCoordinates')
(xp,yp,zp) = C.getMaxValue(ext, 'GridCoordinates')

bc = C.newPyTree(['Base',2])
ext = T.splitSharpEdges(ext, alphaRef=90.)
for zone in Internal.getNodesFromType2(ext, 'Zone_t'):
    if C.getMaxValue(zone,'CoordinateX') != XMax_surf: bc[2][1][2].append(zone)
bc2 = T.splitSharpEdges(bc, alphaRef=30.)

for zone in Internal.getNodesFromType2(bc, 'Zone_t'):
    (xpZone, ypZone) = C.getMaxValue(zone, ['CoordinateX', 'CoordinateY']) 
    if (xpZone == xp and ypZone == yp): Internal._rmNode(bc,zone)
for zone in Internal.getNodesFromType2(bc2, 'Zone_t'):
    bc[2][1][2].append(zone)
C.convertPyTree2File(bc, 'exterior.cgns')

# Add BCs
for zone in Internal.getNodesFromType2(bc, 'Zone_t'):
    dimax = Internal.getZoneDim(zone)[1]-1 #print zone[0],dimax
    [x0,y0,z0] = C.getValue(zone,'GridCoordinates',0)
    [x1,y1,z1] = C.getValue(zone,'GridCoordinates',dimax)
    if (x0 == x1 and x1 == xm): 
        C._addBC2Zone(m, 'IN', 'FamilySpecified:IN', subzone=zone)
        #print 'Xmin'
    if (x0 == x1 and x1 == xp): 
        C._addBC2Zone(m, 'OUT', 'FamilySpecified:OUT', subzone=zone)
        #print 'Xmax'
    if (y0 == y1 and y1 == ym): 
        C._addBC2Zone(m, 'LEFT', 'FamilySpecified:LEFT', subzone=zone)
        #print 'Ymin'
    if (y0 == y1 and y1 == yp): 
        C._addBC2Zone(m, 'RIGHT', 'FamilySpecified:RIGHT', subzone=zone)
        #print 'Ymax'
    if (z0 == z1 and z1 == zm): 
        C._addBC2Zone(m, 'DOWN', 'FamilySpecified:DOWN', subzone=zone)
        #print 'Zmin'
    if (z0 == z1 and z1 == zp): 
        C._addBC2Zone(m, 'TOP', 'FamilySpecified:TOP', subzone=zone)
        #print 'Zmax'

m = C.fillEmptyBCWith(m, 'WALL', 'FamilySpecified:WALL')

tp = C.newPyTree(['Base']); tp[2][1][2] += [m]

C._addFamily2Base(tp, 'TOP', bndType='BCFarfield')
C._addFamily2Base(tp, 'DOWN', bndType='BCFarfield') 
C._addFamily2Base(tp, 'IN', bndType='BCFarfield') 
C._addFamily2Base(tp, 'OUT', bndType='BCFarfield')
C._addFamily2Base(tp, 'LEFT', bndType='BCFarfield') 
C._addFamily2Base(tp, 'RIGHT', bndType='BCFarfield') 
C._addFamily2Base(tp, 'WALL', bndType='BCWall') 

C.convertPyTree2File(tp, case+'Out2.cgns')
