import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import Connector.PyTree as X
import Post.PyTree as P
import Transform.PyTree as T
import Intersector.PyTree as XOR
import Dist2Walls.PyTree as DTW
import Converter.Internal as Internal
import Intersector.PyTree as XOR
import numpy, math

SYM_PLANE_TOL = 1.e-7
toler = 5.e-2
toldist = 1.e-12 # tolerance pour distance entre deux surfaces representant la meme frontiere
tol_match_perio = 1.e-6 # tolerance pour la perio entre les deux frontieres perios

# REMAILLAGE DE LA SURFACE ISSUE DE LEVELSET
def remeshLevelSet(z, X=0, dir=1):
    z = Internal.getZones(z)[0]
    zname = z[0]
    z = C.convertArray2Tetra(z)
    z = T.reorder(z,(1,))
    COORD = '{CoordinateX}'
    if dir == 2 or dir == -2:
        COORD = '{CoordinateY}'
    elif dir == 3 or dir == -3:
        COORD = '{CoordinateZ}'
    formula = '({centers:tag}>%g)'%X
    formulaopp = '({centers:tag}<%g)'%X

    C._initVars(z, '{tag}=%s'%COORD)
    z = C.node2Center(z,["tag"])
    if dir > 0:
        z2 = P.selectCells(z, formula, strict=0)
        z3 = P.selectCells(z, formulaopp, strict=0)
    else:
        z2 = P.selectCells(z, formulaopp, strict=0)
        z3 = P.selectCells(z, formula, strict=0)
    edges3 = P.exteriorFaces(z3); del z3
    G._getVolumeMap(z2)
    smax = C.getMaxValue(z2, 'centers:vol'); hmax = math.sqrt(smax)
    smin = C.getMinValue(z2, 'centers:vol'); hmin = math.sqrt(smin)
    Internal._rmNodesFromType(z2, "FlowSolution_t")
    SCALED = False
    # il faut rescaler pour la tolerance du hausdorf
    if hmax < 1e-2 or hmin < 1e-4:
        SCALED = True
        T._homothety(z, (0.,0.,0.), 1.e3)
        T._homothety(z2, (0.,0.,0.), 1.e3)
        hmin = 1.e3*hmin; hmax = 1.e3*hmax
    #z2 = G.mmgs(z2, ridgeAngle=60., hmin=hmin, hmax=hmax, hausd=0.01, grow=1.1,
    #            anisotropy=1, optim=0, fixedConstraints=[], sizeConstraints=[])
    #edges = P.exteriorFaces(z2); del z2
    z = G.mmgs(z, ridgeAngle=60., hmin=hmin, hmax=hmax, hausd=0.01, grow=1.1,
               anisotropy=1, optim=0, fixedConstraints=[], sizeConstraints=[edges3])
    if SCALED:
        T._homothety(z, (0.,0.,0.), 1.e-3)
    z[0] = zname
    return z

def createDragonMesh0(body, dictOfParams={}, check=False, directory_tmp_files='./'):
    if 'sym_plane_axis' in dictOfParams: sym = dictOfParams['sym_plane_axis'] # symmetry plane at x,y,z constant
    else: sym = None

    if 'sym_plane_values' in dictOfParams: Locsyms = dictOfParams['sym_plane_values']
    else: Locsyms = []

    if 'NParts' in dictOfParams: NParts = dictOfParams['NParts']
    else: NParts = 0

    if 'holeFactor' in dictOfParams: holeFactor = dictOfParams['holeFactor']
    else: holeFactor = 1.25

    if 'hWall' not in dictOfParams:
        hWall = 1.e-6
        print('Warning: createDragonMesh: hmin not defined. Set to %g.'%hWall)
    else: hWall = dictOfParams['hWall'] # hauteur de la premiere maille (si prismes)

    if 'h_factor' in dictOfParams: height_factor = dictOfParams['h_factor']
    else: height_factor = 4

    if 'niter' not in dictOfParams: smoothIter = 100
    else: smoothIter = dictOfParams["niter"]

    if "q" in dictOfParams: raison = dictOfParams["q"]
    else: raison = 1.2

    if 'nlayers' in dictOfParams: nblayerforced = dictOfParams['nlayers'] #0 si pas de forcage
    else: nblayerforced = 0 # auto

    if 'h_ext' in dictOfParams: snear = dictOfParams["h_ext"] # -1 pour calcul auto
    else: snear = -1.

    if 'dfar' in dictOfParams: dfar = dictOfParams["dfar"] # -1 pour calcul auto
    else: dfar = -1.

    if 'h_ext_factor' in dictOfParams: snearFactor = dictOfParams["h_ext_factor"]
    else: snearFactor = 1.


    nbsyms = len(Locsyms)
    if nbsyms > 0:
        locmin = min(Locsyms)
        locmax = max(Locsyms)

        if nbsyms == 2: print('Configuration with two symmetry planes.')
        elif nbsyms == 1: print('Configuration with one symmetry plane.')
    else: print('Configuration in free air.')

    # pour un seul plan de symetrie, on symetrise la geometrie pour assurer que l octree se decoupe sur ce plan
    if nbsyms == 1:
        if sym == 'X': syms = T.symmetrize(body, (Locsyms[0],0.,0.), (0,1,0), (0,0,1))
        elif sym == 'Y': syms = T.symmetrize(body, (0.,Locsyms[0],0.), (1,0,0), (0,0,1))
        elif sym == 'Z': syms = T.symmetrize(body, (0.,0.,Locsyms[0]), (1,0,0), (0,1,0))
        syms[0] = 'syms'
        body = T.join(body,syms); G._close(body, tol=4.e-5)

    T._reorderAll(body, dir=1)

    # on travaille sur des geometries triangulees
    body = C.deleteFlowSolutions__(body)
    body = C.convertArray2Tetra(body)
    body = T.join(body); body = G.close(body)
    T._reorderAll(body, dir=1)

    XMax_surf = C.getMaxValue(body, 'CoordinateX')

    # Calcul automatique du snear (si snear < 0), il peut etre module avec le param snearFactor
    sizemax = 0.
    box = G.bbox(body)
    sizemax = max(sizemax,abs(box[3]-box[0]),abs(box[4]-box[1]),abs(box[5]-box[2]));print('sizemax=',sizemax)
    if snear < 0.:
        G._getNormalMap(body)
        C._magnitude(body, ['centers:sx','centers:sy','centers:sz'])
        meanarea = C.getMeanValue(body,'centers:sMagnitude')
        minarea = C.getMinValue(body,'centers:sMagnitude')
        maxarea = C.getMaxValue(body,'centers:sMagnitude')
        print('surface area min/max/mean=',minarea,maxarea,meanarea)
        # snear.append(2.*(vmin-1)*meanarea**(1/2.))
        snear = snearFactor*meanarea**(1./2.); print('snear=', snear)
        if dfar < 0.: dfar = 10.*sizemax; print('dfar=', dfar)

    # Hauteur evidement, peut etre modulee avec le param holeFactor
    delta = holeFactor*snear
    print('delta=',delta)

    xrayDim1 = 1500; xrayDim2 = 1500

    lay = []
    # Generation de la Boundary layer
    print('Generating the boundary layer...')

    if nblayerforced != 0: nlayer = nblayerforced
    else: nlayer = int(math.log(snear/(height_factor*hWall))/math.log(raison));
    print('Number of prismatic layers = ',nlayer)
    dlayer = hWall*(1.-raison**nlayer)/(1.-raison)
    print('Thickness of prismatic layer: ',dlayer)
    distrib = G.cart((0.,0.,0.), (dlayer/nlayer,1,1), (nlayer,1,1))
    for i in range(0,nlayer): C.setValue(distrib, 'CoordinateX', i, hWall*(1.-raison**i)/(1.-raison))
    hLast = C.getValue(distrib, 'CoordinateX', nlayer-1)-C.getValue(distrib, 'CoordinateX', nlayer-2) ;
    print('hauteur de la derniere couche de prisme = ',hLast)
    #C.convertPyTree2File(body,"surface.cgns")
    lay1 = G.addNormalLayers(body, distrib, check=1, niter=smoothIter)
    lay.append(lay1); # layy = lay[0]  ?
    if check: C.convertPyTree2File(lay, directory_tmp_files+'layer.cgns')
    ext = P.exteriorFaces(lay[0])
    if sym is None:
        ext = T.splitConnexity(ext)
        s = ext[1]
    else:
        ext = T.splitSharpEdges(ext, alphaRef=85.)
        s = ext[-1]
    s = T.breakElements(s)[0]
    s = T.splitConnexity(s)[-1]
    T._reorderAll(s, dir=1)
    if check: C.convertPyTree2File(s, directory_tmp_files+'externalLayer.cgns')

    # Generation de l octree 3D
    print('Generating octree...')
    surfaces = [s] ; snears = [snear]
    o = G.octree(surfaces, snears, dfar=dfar, balancing=1)

    bbox = G.bbox(o)
    expand = int(delta/snear)+2
    print('Expanding layers %d...'%expand)
    for i in range(expand):
        o = G.expandLayer(o, level=0, balancing=1)

    # pour le cas de 2 plans de symetrie, on conserve un quadtree que l on va extruder en Z
    if nbsyms == 2: o = P.isoSurfMC(o,'CoordinateZ',locmin)
    if check: C.convertPyTree2File(o, directory_tmp_files+'octree.cgns')

    eps=1.e-5
    # Extrusion de l octree
    if nbsyms == 2:
        G._getVolumeMap(o)
        minsize = C.getMinValue(o,'centers:vol')**0.5;print(minsize)
        maxsize = C.getMaxValue(o,'centers:vol')**0.5;print(maxsize)
        nblevels = int(math.log(maxsize/minsize)/math.log(2))
        maxnpts = int(round(locmax/minsize))
        maxnpts = 64 # SP to TR : en dur  ?

        tsym = C.newPyTree(['Base'])
        for i in range(nblevels+1):
            area = ((2**i)*minsize)**2
            #C._initVars(t,'{centers:indic}=logical_and({centers:vol}-%g>-%g,{centers:vol}-%g<%g)'%(aera,eps,aera,eps))
            #a1 = P.selectCells2(t, 'centers:indic',strict=0)
            a1 = P.selectCells(o, 'abs({centers:vol}-%f)<%f'%(area,eps), strict=0)
            a1 = Internal.getNodeFromType(a1, 'Zone_t'); a1[0]='octree%d'%i
            npts = int(max(1,maxnpts/2**i));print(i,npts)
            T._addkplane(a1,N=npts)
            T._scale(a1, factor=(1.,1.,locmax/npts))
            #T._scale(a1, factor=(1.,1.,minsize*2**i))
            zmin = C.getMinValue(a1,'CoordinateZ')
            T._translate(a1, (0.,0.,-zmin))
            tsym[2][1][2].append(a1)

        Internal._rmNodesByName(tsym,'FlowSolution*')
        o = Internal.getZones(T.join(tsym))[0]

    # conversion en NGon et conformisation de l octree
    t = C.newPyTree(['Base',o])

    print('Conformizing...')
    t = C.convertArray2NGon(t)
    t = C.conformizeNGon(t) ; G._close(t)

    # Extract blanking surface
    print('Extract blanking surface (%s)...'%delta)
    DTW._distance2Walls(t, s, type='ortho', loc='centers', signed=1)

    s2 = P.isoSurfMC(t, 'centers:TurbulentDistance', delta)
    if nbsyms == 2:
        s2 = [T.splitConnexity(s2)[0]]
        sexts = P.exteriorFaces(s2)
        sexts = T.splitConnexity(sexts)
        for sext in sexts:
            sext = G.tetraMesher(sext, grading=0.2, maxh=0.5*snear, algo=1)
            s2.append(sext)

    # Blanking with delta surface
    print('Blanking octree width delta...')
    BM = numpy.array([[1]], dtype=Internal.E_NpyInt)
    t = X.blankCells(t, [s2], BM, blankingType='center_in', delta=1.e-12, XRaydim1=xrayDim1, XRaydim2=xrayDim2)

    if check: C.convertPyTree2File(t, directory_tmp_files+'blanked.cgns')

    # Selectionne les cellules cellN=1
    print('Selecting...')
    t = P.selectCells2(t, 'centers:cellN', strict=1)
    if check: C.convertPyTree2File(t, directory_tmp_files+'main.cgns')

    # Recupere les surfaces exterieures
    print('Getting exterior faces...')
    ext = P.exteriorFaces(t)
    if nbsyms == 2: ext = T.splitSharpEdges(ext, alphaRef=89.)
    else: ext = T.splitConnexity(ext)
    if check: C.convertPyTree2File(ext, directory_tmp_files+'ext2Watch.cgns')

    # Selectionne la bonne frontiere (celle proche corps)
    if nbsyms != 2:
        ext2 = G.getVolumeMap(ext)
        minsize = C.getMinValue(ext2,'centers:vol')**0.5

        ext = C.newPyTree(['Base'])
        for e in ext2:
            meansize = C.getMeanValue(e,'centers:vol')**0.5
            if abs(meansize-minsize) < eps: ext[2][1][2].append(e)

        ext = T.join(ext)
    else:
        ext2 = C.newPyTree(['Base'])
        for e in ext:
            if abs(C.getMinValue(e,'Coordinate%s'%sym)-C.getMaxValue(e,'Coordinate%s'%sym))>eps: ext2[2][1][2].append(e)
        ext2 = T.join(ext2)
        ext2 = T.splitConnexity(ext2)[0]
        bbox = G.bbox(ext2)
        for e in ext:
            if abs(C.getMinValue(e,'Coordinate%s'%sym)-C.getMaxValue(e,'Coordinate%s'%sym))<eps and abs(C.getMinValue(e,'Coordinate%s'%sym)-bbox[2])>eps and abs(C.getMaxValue(e,'Coordinate%s'%sym)-bbox[5])>eps: ext2[2][1][2].append(e)
        ext = T.join(ext2)

    C._addBC2Zone(t, 'toto', 'BCWall', subzone=ext)
    print("triangulateBC..")
    t = XOR.triangulateBC(t, 'BCWall')
    ext = C.extractBCOfType(t, 'BCWall')[0]
    ext = T.breakElements(ext)[0]
    T._reorderAll(ext, dir=1)

    # pour un seul plan de symetrie, on selectionne la partie du maillage qui interesse
    if nbsyms == 1:
        ext = P.selectCells(ext, '{Coordinate%s}>%f'%(sym,Locsyms[0]), strict=0)
        s = P.selectCells(s, '{Coordinate%s}>%f'%(sym,Locsyms[0]-0.0001), strict=1)
        lay = P.selectCells(lay, '{Coordinate%s}>%f'%(sym,Locsyms[0]-0.0001), strict=1)

    # Remplit l'espace entre l'exterieur faces, la surface et le remplissage tri sur le(s) plan(s) de symetrie
    print('Tetra filling...')
    ext = XOR.conformUnstr(ext, tol=0., itermax=1); ext = T.splitManifold(ext)
    if len(ext) > 1: print("Warning: Tetra/octree border has manifolds. Tetra mesher might fail.")
    ext = ext[0]
    T._reorderAll(ext, dir=1)
    T._reorderAll(s, dir=-1)

    if check: C.convertPyTree2File([s,ext], directory_tmp_files+'tetra_surfs.cgns')
    s_in = [s, ext]

    if nbsyms > 0:
        extlines = P.exteriorFaces(ext)
        extlines = T.splitConnexity(extlines)
        slines = P.exteriorFaces(s)
        slines = T.splitConnexity(slines)
        surfsmin = []; surfsmax = []
        for i in range(len(extlines)):
            extlineloc = C.getMeanValue(extlines[i],'Coordinate%s'%sym)
            if abs(extlineloc - locmin)<eps: surfsmin.append(extlines[i])
            if abs(extlineloc - locmax)<eps: surfsmax.append(extlines[i])
            slineloc = C.getMeanValue(slines[i],'Coordinate%s'%sym)
            if abs(slineloc - locmin)<eps: surfsmin.append(slines[i])
            if abs(slineloc - locmax)<eps: surfsmax.append(slines[i])

        surfsym1 = G.tetraMesher(surfsmin, grading=0.2, maxh=0.5*snear, algo=1)
        s_in += [surfsym1]
        if nbsyms == 2:
            surfsym2 = G.tetraMesher(surfsmax, grading=0.2, maxh=0.5*snear, algo=1)
            s_in += [surfsym2]

    s_in = T.join(s_in); s_in = XOR.conformUnstr(s_in, tol=0., itermax=1); T._reorderAll(s_in,dir=1)
    m = G.tetraMesher([s_in], grading=0.2, maxh=0.5*snear, algo=1)
    m = C.convertArray2NGon(m)
    if check: C.convertPyTree2File(m, directory_tmp_files+'tetra1.cgns')

    # Assemblage de la couche limite, l octree et la couche intermediaire tetra
    print('Reconnecting...')
    lay = C.convertArray2NGon(lay)
    m = T.join(m, lay)
    print('Final closing...')
    if nbsyms == 1:
        t = P.selectCells(t, '{Coordinate%s}>%f'%(sym,Locsyms[0]), strict=0)

    Internal._rmNodesByName(t,'FlowSolution*')
    if check: C.convertPyTree2File(t,directory_tmp_files+'tmp.cgns')
    m = T.join(m, t)
    G._close(m)
    #print(Internal.getZoneDim(m))

    m = XOR.closeCells(m)
    err = XOR.checkCellsClosure(m);print("CheckCellsClosure status : %d"%(err))
    if check: C.convertPyTree2File(m, directory_tmp_files+'meshOut.cgns')

    # Creation des BCs
    print('Create BCs...')
    ext = P.exteriorFaces(m)
    xm = C.getMinValue(ext,'CoordinateX'); ym = C.getMinValue(ext,'CoordinateY'); zm = C.getMinValue(ext,'CoordinateZ')
    xp = C.getMaxValue(ext,'CoordinateX'); yp = C.getMaxValue(ext,'CoordinateY'); zp = C.getMaxValue(ext,'CoordinateZ')
    #print(xm,xp,ym,yp,zm,zp)

    bc = C.newPyTree(['Base',2])

    ext = T.splitSharpEdges(ext, alphaRef=90.)
    for zone in Internal.getByType(ext,'Zone_t')[2]:
        if C.getMaxValue(zone,'CoordinateX') != XMax_surf: bc[2][1][2].append(zone)
    bc2 = T.splitSharpEdges(bc, alphaRef=30.)

    for zone in Internal.getByType(bc,'Zone_t')[2]:
        if C.getMaxValue(zone,'CoordinateX') == xp and C.getMaxValue(zone,'CoordinateY') == yp:
            Internal._rmNode(bc,zone)
    for zone in Internal.getByType(bc2,'Zone_t')[2]:
        bc[2][1][2].append(zone)
    if check: C.convertPyTree2File(bc, directory_tmp_files+'exterior.cgns')

    # Add BCs
    eps = 1.e-3
    for zone in Internal.getByType(bc,'Zone_t')[2]:
        dimax = Internal.getZoneDim(zone)[1] - 1
        bary = G.barycenter(zone)
        #print(zone[0],dimax,bary)
        if abs(bary[0] - xm) < eps: C._addBC2Zone(m,'IN','FamilySpecified:IN',subzone=zone); #print('Xmin')
        if abs(bary[0] - xp) < eps: C._addBC2Zone(m,'OUT','FamilySpecified:OUT',subzone=zone); #print('Xmax')
        if abs(bary[1] - ym) < eps: C._addBC2Zone(m,'DOWN','FamilySpecified:DOWN',subzone=zone); #print('Ymin')
        if abs(bary[1] - yp) < eps: C._addBC2Zone(m,'TOP','FamilySpecified:TOP',subzone=zone); #print('Ymax')
        if abs(bary[2] - zm) < eps: C._addBC2Zone(m,'RIGHT','FamilySpecified:RIGHT',subzone=zone); #print('Zmin')
        if abs(bary[2] - zp) < eps: C._addBC2Zone(m,'LEFT','FamilySpecified:LEFT',subzone=zone); #print('Zmax')

    C._fillEmptyBCWith(m,'WALL', 'FamilySpecified:WALL')

    tp = C.newPyTree(['Base']); tp[2][1][2] += [m]

    print('add families')
    base = Internal.getNodeFromType1(tp, 'CGNSBase_t')
    #print(base)
    if nbsyms > 0:
        C._addFamily2Base(base, 'RIGHT', bndType='BCSymmetryPlane')
        if nbsyms == 2: C._addFamily2Base(base, 'LEFT', bndType='BCSymmetryPlane')
        if nbsyms == 1: C._addFamily2Base(base, 'LEFT', bndType='BCFarfield')
    else:
        C._addFamily2Base(base, 'RIGHT', bndType='BCFarfield')
        C._addFamily2Base(base, 'LEFT', bndType='BCFarfield')
    C._addFamily2Base(base, 'IN', bndType='BCFarfield')
    C._addFamily2Base(base, 'OUT', bndType='BCFarfield')
    C._addFamily2Base(base, 'DOWN', bndType='BCFarfield')
    C._addFamily2Base(base, 'TOP', bndType='BCFarfield')
    C._addFamily2Base(base, 'WALL', bndType='BCWall')

    # Decoupage en NBlocs

    if NParts > 0:
        print('Decoupage en ',NParts,' blocs')
        tp = T.splitNParts(tp,NParts)

    return tp

def uniformizeMesh(z, N, dir):
    N = int(N)
    distrib = G.cart( (0,0,0), (1./(N-1),1,1), (N,1,1) )
    zp = G.map(z, distrib, dir)
    return zp

def createDragonMeshForBladeInChannel(ts, dictOfParams={}, check=False, directory_tmp_files='./'):
    BA2BF = -1
    Internal._rmNodesFromType(ts,'FlowSolution_t')
    bases = Internal.getZones(ts)
    if len(bases) != 3:
        raise("ValueError: createDragonMesh: 3 zones must be defined: HUB/SHROUD/BLADE")
        return None

    surf_hub = Internal.getNodeFromName2(ts, "HUB")
    if surf_hub is None or surf_hub == []:
        raise("ValueError: no base/zone of name HUB found.")
        return None

    surf_shroud = Internal.getNodeFromName2(ts,"SHROUD")
    if surf_shroud is None or surf_shroud==[]:
        raise("ValueError: no base/zone of name SHROUD found.")
        return None

    surf_blade = Internal.getNodeFromName2(ts,'BLADE')
    if surf_blade is None or surf_blade==[]:
        raise("ValueError: no base/zone of name BLADE found.")
        return None

    if 'remesh_input' in dictOfParams:
        remeshBorders = dictOfParams['remesh_input']
        # SP : preferer le remaillage de l aube en externe
        #remesh_aube = dictOfParams['remesh_input']

    if 'niter' not in dictOfParams: smoothIter = 100
    else: smoothIter = dictOfParams["niter"]

    if 'hWall' not in dictOfParams:
        hWall = 1e-6
        print('Warning: createDragonMesh: hmin not defined. Set to %g.'%hWall)
    else: hWall = dictOfParams['hWall'] # hauteur de la premiere maille (si prismes)

    if "q" in dictOfParams: raison = dictOfParams["q"]
    else: raison = 1.2

    if 'nlayers' in dictOfParams: nlayer = dictOfParams['nlayers']
    else:
        nlayer = 10 # auto
        print('Warning: createDragonMesh: nlayer not defined. Set to %d.'%nlayer)

    typePerio = -1 # rotation : 0, 1 : translation
    if "periodicity" in dictOfParams:
        if dictOfParams["periodicity"]=='rotation': typePerio = 0
        else: typePerio = 1
    else:
        print(" periodicity key not found in dictOfParams. Set to default (rotation).")
        typePerio = 0

    centre = (0.,0.,0.); axis = (0.,0.,1.); THETA = 0
    translx = 0; transly = 0; translz = 0

    if typePerio == 0:
        if 'nb_blades' in dictOfParams: nb_aubes_fan = dictOfParams["nb_blades"]
        else:
            raise ValueError("createDragonMesh: nb_blades required for periodicity by rotation.")
            return None

        if 'rotation_axis' in dictOfParams: axis = dictOfParams['rotation_axis']
        else:
            print('Warning: createDragonMesh: rotation axis not defined. Set to default :', axis)

        if 'rotation_center' in dictOfParams: centre = dictOfParams['rotation_center']
        else:
            print('Warning: createDragonMesh: rotation center not defined. Set to default :', centre)

        THETA = 360./dictOfParams["nb_blades"]
    else:
        transl_vct = dictOfParams['translation_vct']
        translx = transl_vct[0]
        transly = transl_vct[1]
        translz = transl_vct[2]

    hext = hWall*raison**(nlayer) # cell height at outer prism layer
    # projection of the blade root and tip onto the hub and shroud
    surf_aube = Internal.getZones(surf_blade)[0]
    surf_shroud = Internal.getZones(surf_shroud)[0]
    surf_hub = Internal.getZones(surf_hub)[0]
    lines_ext_aube = P.exteriorFaces(surf_aube)
    lines_ext_aube = T.splitConnexity(lines_ext_aube)
    # which one is the root/tip ?
    DTW._distance2Walls(lines_ext_aube, surf_hub,loc='nodes')
    distl0 = C.getMaxValue(lines_ext_aube[0],'TurbulentDistance')
    distl1 = C.getMaxValue(lines_ext_aube[1],'TurbulentDistance')
    if distl0 < distl1:
        lines_ext_aube[0][0] = 'line_BLADE_HUB'
        lines_ext_aube[1][0] = 'line_BLADE_SHROUD'
    else:
        lines_ext_aube[0][0] = 'line_BLADE_SHROUD'
        lines_ext_aube[1][0] = 'line_BLADE_HUB'

    line_aube_hub = Internal.getNodeFromName(lines_ext_aube, 'line_BLADE_HUB')
    line_aube_shroud = Internal.getNodeFromName(lines_ext_aube, 'line_BLADE_SHROUD')

    hook = C.createHook(surf_aube, function='nodes')
    nodesMatch = C.identifyNodes(hook, line_aube_hub)
    T._projectOrtho(line_aube_hub, surf_hub)

    for noind in range(len(nodesMatch)):
        inds = nodesMatch[noind]
        if inds > -1:
            xp = C.getValue(line_aube_hub, 'GridCoordinates', noind)
            C.setValue(surf_aube, 'CoordinateX', inds-1, xp[0])
            C.setValue(surf_aube, 'CoordinateY', inds-1, xp[1])
            C.setValue(surf_aube, 'CoordinateZ', inds-1, xp[2])

    nodesMatch = C.identifyNodes(hook, line_aube_shroud)
    T._projectOrtho(line_aube_shroud, surf_shroud)

    for noind in range(len(nodesMatch)):
        inds = nodesMatch[noind]
        if inds > -1:
            xp = C.getValue(line_aube_shroud, 'GridCoordinates', noind)
            C.setValue(surf_aube, 'CoordinateX', inds-1, xp[0])
            C.setValue(surf_aube, 'CoordinateY', inds-1, xp[1])
            C.setValue(surf_aube, 'CoordinateZ', inds-1, xp[2])

    C.freeHook(hook)

    # Creation of channel borders - HUB
    lines_hub = extractExternalLines__(surf_hub, surf_aube, BA2BF)
    lines_perios_hub = lines_hub[:2]
    line_in_h = lines_hub[2];  line_out_h = lines_hub[3]

    line_periom_h = lines_perios_hub[0]
    line_periop_h = lines_perios_hub[1]
    if typePerio == 1:
        periomdup = T.translate(line_periom_h,(translx, transly, translz))
    else:
        periomdup = T.rotate(line_periom_h, centre, axis, THETA)
    periomdup = C.diffArrays(line_periop_h,periomdup)
    C._initVars(periomdup,'{dist}=sqrt({DCoordinateX}**2+{DCoordinateY}**2+{DCoordinateZ}**2)')
    if C.getMaxValue(periomdup,'dist')<tol_match_perio:
        pass
    else:
        line_periom_h = lines_perios_hub[1]
        line_periop_h = lines_perios_hub[0]

    # Creation of channel borders -  SHROUD
    lines_shroud= extractExternalLines__(surf_shroud, surf_aube, BA2BF)
    lines_perios_shroud = lines_shroud[:2]
    line_in_s = lines_shroud[2];  line_out_s = lines_shroud[3]
    line_periom_s = lines_perios_shroud[0]
    line_periop_s = lines_perios_shroud[1]
    if typePerio == 1:
        periomdup = T.translate(line_periom_s,(translx, transly, translz))
    else:
        periomdup = T.rotate(line_periom_s, centre, axis, THETA)

    periomdup = C.diffArrays(line_periop_s,periomdup)
    C._initVars(periomdup,'{dist}=sqrt({DCoordinateX}**2+{DCoordinateY}**2+{DCoordinateZ}**2)')
    if C.getMaxValue(periomdup,'dist')<tol_match_perio:pass
    else:
        line_periom_s = lines_perios_shroud[1]
        line_periop_s = lines_perios_shroud[0]

    # retrieve orders to create inlet/outlet border from hub to shroud
    NPTS_L = C.getNPts(line_in_h)
    ptA = getMatchingPoint__(line_in_h, line_periom_h)
    ptB = getMatchingPoint__(line_in_s, line_periom_s)
    linem_in = D.line(ptA, ptB, N=NPTS_L)
    ptA = getMatchingPoint__(line_out_h, line_periom_h)
    ptB = getMatchingPoint__(line_out_s, line_periom_s)
    linem_out = D.line(ptA, ptB, N=NPTS_L)
    surf_periom = G.TFI([line_periom_h,line_periom_s, linem_in, linem_out]);
    surf_periom[0] = 'PERIOM'

    if typePerio == 0:
        linep_in  = T.rotate(linem_in, centre, axis, THETA)
        linep_out = T.rotate(linem_out, centre, axis, THETA)
    else:
        linep_in  = T.translate(linem_in,(translx,transly,translz))
        linep_out = T.translate(linem_out,(translx,transly,translz))
    linep_in[0]=linem_in[0]+'-dup'
    linep_out[0]=linem_out[0]+'-dup'

    surf_periop = T.translate(surf_periom,(translx,transly,translz))
    surf_periop[0]='PERIOP'

    surf_inlet = G.TFI([linem_in,linep_in,line_in_h,line_in_s])
    surf_inlet[0] = 'INLET'
    surf_outlet = G.TFI([linem_out,linep_out,line_out_h,line_out_s])
    surf_outlet[0] = 'OUTLET'

    lines_ext_aube = P.exteriorFaces(surf_aube)
    lines_ext_aube = T.splitConnexity(lines_ext_aube)
    #which one is hub and shroud ?
    DTW._distance2Walls(lines_ext_aube,surf_hub, loc='nodes')
    d0 = C.getMaxValue(lines_ext_aube[0], 'TurbulentDistance')
    d1 = C.getMaxValue(lines_ext_aube[1], 'TurbulentDistance')
    Internal._rmNodesFromType(lines_ext_aube,'FlowSolution_t')
    if d0 < d1:
        line_aube_hub = lines_ext_aube[0]
        line_aube_shroud = lines_ext_aube[1]
    else:
        line_aube_hub = lines_ext_aube[1]
        line_aube_shroud = lines_ext_aube[0]

    mesh_hub = generateTriMeshBetweenContours__(line_aube_hub, surf_hub,
                                                hmin=hWall,
                                                remeshBorders=remeshBorders)

    mesh_shroud = generateTriMeshBetweenContours__(line_aube_shroud, surf_shroud,
                                                   hmin=hWall,
                                                   remeshBorders=remeshBorders)

    surfs = T.join([mesh_hub,mesh_shroud, surf_aube])
    surfs = T.reorder(surfs,(-1,))

    ts = C.newPyTree(['HUB','SHROUD','BLADE','AMONT','AVAL','PERIODIC'])
    ts[2][1][2] = Internal.getZones(mesh_shroud)
    ts[2][2][2] = Internal.getZones(mesh_hub)
    ts[2][3][2] = Internal.getZones(surf_aube)
    ts[2][4][2] = Internal.getZones(surf_inlet)
    ts[2][5][2] = Internal.getZones(surf_outlet)
    ts[2][6][2] = [surf_periop, surf_periom]
    # input: surfaces bases ['WALL','AMONT','AVAL','PERIODIC']
    print('Generating boundary layers...')
    d = G.cart((0.,0.,0.), (0.1,1,1), (nlayer,1,1))
    for i in range(0,nlayer): C.setValue(d, 'CoordinateX', i, hWall*(1.-raison**i)/(1.-raison))
    h_ext = C.getValue(d,'CoordinateX',nlayer-1)-C.getValue(d,'CoordinateX',nlayer-2)

    surfs_wall=Internal.getZones(surfs)
    prismatic = True
    volmin = -1.e10
    if prismatic:
        lay = G.addNormalLayers(surfs_wall, d, check=0, niterType=0, niter=smoothIter, niterK=[],
                                smoothType=0, eps=0.4, nitLocal=3,
                                kappaType=0, kappaS=[0.2,1.6], blanking=False, algo=0)

        lay = C.convertArray2NGon(lay, recoverBC=0)
        G._getVolumeMap(lay)
        volmin = C.getMinValue(lay,'centers:vol')
        C._rmVars(lay, ["centers:vol"])
        if check:
            C.convertPyTree2File(lay,"lay.cgns")
    if volmin < -1.e-16 or not prismatic:
        print("Warning: negative cells found in prismatic layer !")
        print("Full TETRA mesh generation is applied.")
        prismatic = False
        del lay
        extBorders = [surf_inlet, surf_outlet, surf_periom, surf_periop]
        surfs_wall += extBorders
        surfs = C.convertArray2Tetra(surfs_wall,split='withBarycenters')

        surfs = T.join(surfs); surfs = T.reorderAll(surfs)
        mesh_final = G.tetraMesher(surfs)
        npts = Internal.getZoneDim(mesh_final)[1]
        if npts == 0:
            raise ValueError("createDragonMesh: tetraMesher failed. Please check input data.")
            return None
        mesh_final = C.convertArray2NGon(mesh_final)
    else:
        # Reprojection des frontieres QUAD sur les frontieres du domaine
        Internal._rmNodesByName(lay,'FlowSol*')
        bcc = P.exteriorFaces(lay)
        bcc = T.breakElements(bcc)

        ext_QUAD = []; ext_TRI = []
        for bc in bcc:
            zdim = Internal.getZoneDim(bc)
            if zdim[3] == 'QUAD': ext_QUAD.append(bc)

        ext_QUAD = T.join(ext_QUAD)
        hook = C.createHook(lay, function='nodes')
        surf_inlet[0] = 'INLET'; surf_outlet[0]='OUTLET'
        surf_periom[0]='PERIOM'; surf_periop[0]= 'PERIOP'
        extBorders = [surf_inlet, surf_outlet, surf_periom, surf_periop]
        extBorders = C.convertArray2Tetra(extBorders)

        ext_QUAD = T.splitConnexity(ext_QUAD) # hub and shroud parts
        for ext_QUAD0 in ext_QUAD:
            ext_QUAD0 = T.splitSharpEdges(ext_QUAD0,30.)
            allnodes = []
            for extq in ext_QUAD0:
                nodes = C.identifyNodes(hook, extq)
                allnodes.append(nodes)

            # force periodicity
            DTW._distance2Walls(ext_QUAD0,surf_periom,loc='nodes')
            C._initVars(ext_QUAD0,"{distM}={TurbulentDistance}")
            DTW._distance2Walls(ext_QUAD0,surf_periop, loc='nodes')
            C._initVars(ext_QUAD0,"{distP}={TurbulentDistance}")
            noperm = -1; noperp = -1
            distm = 1e10; distp = 1e10
            for nol, zl in enumerate(ext_QUAD0):
                dl = C.getMaxValue(zl,'distM')
                if dl < distm:
                    distm = dl; noperm = nol
                dl = C.getMaxValue(zl,'distP')
                if dl < distp:
                    distp = dl; noperp = nol
            ext_QUADPERP = None

            projsurf = T.projectOrtho(ext_QUAD0[noperm], extBorders)
            if typePerio == 1:
                ext_QUADPERP = T.translate(projsurf,(translx, transly, translz))
            else:
                ext_QUADPERP = T.rotate(projsurf,centre, axis, THETA)

            nodes = allnodes[noperp]
            for ind in nodes:
                coord = tuple(C.getValue(lay, 'GridCoordinates', ind-1))
                indproj = D.getNearestPointIndex(ext_QUADPERP, coord)[0]
                coordproj = C.getValue(ext_QUADPERP, 'GridCoordinates', indproj)
                C.setValue(lay, 'CoordinateX', ind-1, coordproj[0])
                C.setValue(lay, 'CoordinateY', ind-1, coordproj[1])
                C.setValue(lay, 'CoordinateZ', ind-1, coordproj[2])

            nodes = allnodes[noperm]
            for ind in nodes:
                coord = tuple(C.getValue(lay, 'GridCoordinates', ind-1))
                indproj = D.getNearestPointIndex(projsurf, coord)[0]
                coordproj = C.getValue(projsurf, 'GridCoordinates', indproj)
                C.setValue(lay, 'CoordinateX', ind-1, coordproj[0])
                C.setValue(lay, 'CoordinateY', ind-1, coordproj[1])
                C.setValue(lay, 'CoordinateZ', ind-1, coordproj[2])

        C.freeHook(hook)
        #---------------------------------------
        # on cherche la frontiere ext TRI sur le maillage reprojete
        bcc = P.exteriorFaces(lay)
        bcc = T.breakElements(bcc)

        ext_TRI = []
        for bc in bcc:
            zdim = Internal.getZoneDim(bc)
            if zdim[3]=='TRI': ext_TRI.append(bc)
        DTW._distance2Walls(ext_TRI, surfs_wall, loc='nodes')

        ext_TRI0 = T.splitConnexity(ext_TRI)
        ext_TRI=[]
        for ext in ext_TRI0:
            if C.getMaxValue(ext,'TurbulentDistance')>hWall: ext_TRI.append(ext)
        ext_TRI = T.join(ext_TRI)
        T._reorder(ext_TRI,(1,))
        Internal._rmNodesFromType(ext_TRI,'FlowSolution_t')
        tetMesh = createInternalTetraMesh__(ext_TRI, ts, hext)
        tetMesh = C.convertArray2NGon(tetMesh)
        mesh_final = T.join(tetMesh,lay)

    G._close(mesh_final, toldist)
    G._getVolumeMap(mesh_final)
    volmin = C.getMinValue(mesh_final,'centers:vol')
    print("final min volume =",volmin)
    Internal._rmNodesByName(mesh_final,'FlowSol*')

    # Add BCs
    C._addBC2Zone(mesh_final, 'HUB', 'FamilySpecified:HUB', subzone=mesh_hub)
    C._addBC2Zone(mesh_final, 'SHROUD', 'FamilySpecified:SHROUD', subzone=mesh_shroud)
    C._addBC2Zone(mesh_final, 'BLADE', 'FamilySpecified:BLADE', subzone=surf_aube)

    extFaces = P.exteriorFaces(mesh_final); extFaces = T.splitSharpEdges(extFaces, 60)
    extBCs = ['INLET','OUTLET','PERIOP','PERIOM']
    ZBC={}
    for bcname in extBCs:
        zbc = Internal.getZones(Internal.getNodeFromName(ts,bcname))[0]
        DTW._distance2Walls(extFaces, zbc, loc='nodes')
        distmin = 1.e10; efound = -1
        for noe, ze in enumerate(extFaces):
            d1 = C.getMaxValue(ze,'TurbulentDistance')
            if d1 < distmin:
                distmin = d1; efound = noe
        if efound > -1:
            C._addBC2Zone(mesh_final,bcname,'FamilySpecified:%s'%bcname,subzone=extFaces[efound])

    print('add families')
    tp = C.newPyTree(['Base']); tp[2][1][2] += [mesh_final]
    base = Internal.getNodeFromType1(tp, 'CGNSBase_t')
    C._addFamily2Base(base, 'INLET', bndType='BCInflow')
    C._addFamily2Base(base, 'OUTLET', bndType='BCOutflow')
    C._addFamily2Base(base, 'HUB', bndType='BCWall')
    C._addFamily2Base(base, 'SHROUD', bndType='BCWall')
    C._addFamily2Base(base, 'BLADE', bndType='BCWall')
    return tp

def createInternalTetraMesh__(ext_TRI, ts, hext):
    # ts  : ['HUB','SHROUD','BLADE','AMONT','AVAL','PERIODIC']
    mesh_cart = Internal.getNodesFromName1(ts, 'SHROUD')
    mesh_spin = Internal.getNodesFromName1(ts, 'HUB')
    mesh_amont = Internal.getNodesFromName1(ts, 'AMONT')
    mesh_aval = Internal.getNodesFromName1(ts, 'AVAL')
    mesh_perios =  Internal.getNodesFromName1(ts, 'PERIODIC')

    ext = P.exteriorFaces(ext_TRI)
    # ATTENTION REORDER COHERENT CAR NORMALES VERS L INTERIEUR
    # on remet les ext tous dans le meme sens
    ext = T.join(ext)
    ext = T.splitConnexity(ext)
    ext = C.convertBAR2Struct(ext)
    T._reorder(ext[1],(-1,2,3))

    lines_amont=[]; lines_aval = []; perios=[]
    for e in ext:
        DTW._distance2Walls(e, mesh_perios, loc='nodes')
        dl = D.getLength(e)/C.getNPts(e)
        toldistrel = 0.1*dl
        e1 = P.selectCells(e,'{TurbulentDistance}>%g'%toldistrel,strict=0)
        e1 = T.splitConnexity(e1) # on doit en avoir 2
        DTW._distance2Walls(e1, mesh_amont, loc='nodes')

        if C.getMinValue(e1[0],'TurbulentDistance')<C.getMinValue(e1[1],'TurbulentDistance'):
            lines_amont.append(e1[0])
            lines_aval.append(e1[1])
        else:
            lines_amont.append(e1[1])
            lines_aval.append(e1[0])
        e1 = P.selectCells(e,'{TurbulentDistance}<%g'%toldistrel,strict=1)
        e1 = T.splitConnexity(e1) # on doit en avoir 2
        perios += e1

    lines_amont = C.convertBAR2Struct(lines_amont)
    lines_aval = C.convertBAR2Struct(lines_aval)
    perios = C.convertBAR2Struct(perios)

    # TFI AMONT
    Internal._rmNodesFromType(lines_amont,'FlowSolution_t')
    Internal._rmNodesFromType(lines_aval,'FlowSolution_t')
    Internal._rmNodesFromType(perios,'FlowSolution_t')

    npts = C.getNPts(lines_amont[0])
    zperiop = Internal.getZones(mesh_perios)[0]
    DTW._distance2Walls(lines_amont, zperiop, loc='nodes')

    d0 = C.getValue(lines_amont[0],'TurbulentDistance',0)
    d1 = C.getValue(lines_amont[0],'TurbulentDistance',npts-1)
    if d0 < d1:
        pt1 = C.getValue(lines_amont[0],'GridCoordinates',0)
        pt3 = C.getValue(lines_amont[0],'GridCoordinates',npts-1)
    else:
        pt1 = C.getValue(lines_amont[0],'GridCoordinates',npts-1)
        pt3 = C.getValue(lines_amont[0],'GridCoordinates',0)

    d0 = C.getValue(lines_amont[1],'TurbulentDistance',0)
    d1 = C.getValue(lines_amont[1],'TurbulentDistance',npts-1)
    if d0 < d1:
        pt2 = C.getValue(lines_amont[1],'GridCoordinates',0)
        pt4 = C.getValue(lines_amont[1],'GridCoordinates',npts-1)
    else:
        pt2 = C.getValue(lines_amont[1],'GridCoordinates',npts-1)
        pt4 = C.getValue(lines_amont[1],'GridCoordinates',0)

    line1 = D.line(pt1, pt2,N=npts); line2 = D.line(pt3, pt4,N=npts)
    Internal._rmNodesFromType(lines_amont, 'FlowSolution_t')
    lines = T.join([line1,line2]+lines_amont)
    lines = T.splitSharpEdges(lines)
    lines = C.convertBAR2Struct(lines)
    AMONT = G.TFI(lines)

    # TFI AVAL
    npts = C.getNPts(lines_aval[0])
    DTW._distance2Walls(lines_aval, zperiop, loc='nodes')
    d0 = C.getValue(lines_aval[0],'TurbulentDistance',0)
    d1 = C.getValue(lines_aval[0],'TurbulentDistance',npts-1)
    if d0 < d1:
        pt1 = C.getValue(lines_aval[0],'GridCoordinates',0)
        pt3 = C.getValue(lines_aval[0],'GridCoordinates',npts-1)
    else:
        pt1 = C.getValue(lines_aval[0],'GridCoordinates',npts-1)
        pt3 = C.getValue(lines_aval[0],'GridCoordinates',0)

    d0 = C.getValue(lines_aval[1],'TurbulentDistance',0)
    d1 = C.getValue(lines_aval[1],'TurbulentDistance',npts-1)
    if d0 < d1:
        pt2 = C.getValue(lines_aval[1],'GridCoordinates',0)
        pt4 = C.getValue(lines_aval[1],'GridCoordinates',npts-1)
    else:
        pt2 = C.getValue(lines_aval[1],'GridCoordinates',npts-1)
        pt4 = C.getValue(lines_aval[1],'GridCoordinates',0)

    line3 = D.line(pt1, pt2,N=npts); line4 = D.line(pt3, pt4,N=npts)
    Internal._rmNodesFromType(lines_aval,'FlowSolution_t')
    lines = T.join([line3,line4]+lines_aval)
    lines = T.splitSharpEdges(lines)
    lines = C.convertBAR2Struct(lines)
    AVAL = G.TFI(lines)

    perios = C.convertArray2Tetra([line1,line2,line3,line4]+perios)
    perios = T.join(perios)
    perios = T.splitConnexity(perios)
    sides = []
    for e in perios:
        alp0 = 80.
        nsplit = 0
        while nsplit < 4:
            zsplit = T.splitSharpEdges(e, alp0)
            alp0 = alp0-10.
            nsplit = len(zsplit)
        if nsplit == 4:
            zsplit = C.convertBAR2Struct(zsplit)
            sides.append(zsplit)

    TFIPERM = G.TFI(sides[0])
    TFIPERP = G.TFI(sides[1])

    # remap
    ZONES = [AMONT,AVAL,TFIPERM,TFIPERP]
    ZONES = T.join(ZONES)
    ZONES = remapSurf__(ZONES, hmin=hext, dir=2)
    ZONES = C.convertArray2Tetra(ZONES)
    surfs = Internal.getZones(ZONES)+Internal.getZones(ext_TRI)
    surfs = T.join(surfs)
    surfs = T.reorder(surfs,(1,))
    tetMesh = G.tetraMesher(surfs)
    return tetMesh

def remapSurf__(z, hmin=1e-6, dir=1):
    diri = dir
    dimZ = Internal.getZoneDim(z)
    NIDIR = dimZ[diri]
    if diri == 1: zt = T.subzone(z, (1,1,1), (NIDIR,1,1))
    else:
        zt = T.subzone(z, (1,1,1), (1,NIDIR,1))
        zt = T.reorder(zt, (2,1,3))

    l = D.getLength(zt)
    D._getCurvilinearAbscissa(zt)
    distrib = C.cpVars(zt, 's', zt, 'CoordinateX')
    C._initVars(distrib, 'CoordinateY', 0.)
    C._initVars(distrib, 'CoordinateZ', 0.)
    C._rmVars(distrib, 's')
    Npts_orig = Internal.getZoneDim(distrib)[1]
    Nr = int(0.1*Npts_orig)
    distrib = G.enforcePlusX(distrib, hmin/l, Nr, Nr)
    Npts_fin = Internal.getZoneDim(distrib)[1]
    Nr2 = int(0.1*Npts_fin)
    distrib = G.enforceMoinsX(distrib, hmin/l, Nr2, Nr2)
    z = G.map(z, distrib, diri)
    return z


# IN : surface hub or shroud
# OUT : periodic lines, inlet and outlet
def extractExternalLines__(surf, surf_aube, BA2BF=1):
    lines_ext = P.exteriorFaces(surf)
    DTW._distance2Walls(lines_ext, surf_aube, loc='nodes')
    lines_ext = T.splitConnexity(lines_ext)
    lines_ext = C.convertBAR2Struct(lines_ext)
    distl0 = C.getMaxValue(lines_ext[0],'TurbulentDistance')
    distl1 = C.getMaxValue(lines_ext[1],'TurbulentDistance')
    Internal._rmNodesFromType(lines_ext,'FlowSolution_t')

    if distl0 < distl1: line_ext = lines_ext[1]
    else: line_ext = lines_ext[0]

    lines_ext = []; nsplit = 0; alp0 = 80.
    while nsplit < 4:
        zsplit = T.splitSharpEdges(line_ext, alp0)
        alp0 = alp0-5.
        nsplit = len(zsplit)
        if nsplit >= 4:
            lines_ext = zsplit
            break

    no_inflow=-1; no_outflow = -1
    dmin = 1e10; dmax = 1e10
    bb_all = G.bbox(lines_ext)
    if BA2BF == 1:
        imin = 0; imax = 3 # amont/aval selon axe X+
    elif BA2BF ==-1:
        imin = 3; imax = 0 # amont/aval selon axe X-
    elif BA2BF == 2:
        imin = 1; imax = 4 # amont/aval selon axe Y+
    elif BA2BF ==-2:
        imin = 4; imax = 1 # amont/aval selon axe Y-
    elif BA2BF == 3:
        imin = 2; imax = 5 # amont/aval selon axe Z+
    elif BA2BF ==-3:
        imin = 5; imax = 2 # amont/aval selon axe Z-

    for nol0 in range(len(lines_ext)):
        l0 = lines_ext[nol0]
        bbl = G.bbox(l0)
        dl = abs(bbl[imax]-bbl[imin])
        if abs(bbl[imin]-bb_all[imin])<1e-10:
            if dl < dmin:
                if BA2BF > 0: no_inflow=nol0
                else: no_outflow=nol0
                dmin = dl
        elif abs(bbl[imax]-bb_all[imax])<1e-10:
            if dl < dmax:
                if BA2BF < 0: no_inflow=nol0
                else: no_outflow=nol0
                dmax = dl
    lines_perios = []
    for noz in range(len(lines_ext)):
        if noz != no_inflow and noz != no_outflow:
            lines_perios.append(lines_ext[noz])

    lines_perios = T.join(lines_perios)
    lines_perios = T.splitConnexity(lines_perios)
    lines_perios = C.convertBAR2Struct(lines_perios)
    lines_ext = C.convertBAR2Struct(lines_ext)
    return lines_perios+[lines_ext[no_inflow], lines_ext[no_outflow]]

# lineb : contour de l aube au moyeu/carter typiquement
# surfp : surface de projection (carter/moyeu)
# line_ext : ligne exterieure
def orderExteriorEdges__(surf, lineb):
    lines_ext = P.exteriorFaces(surf)
    lines_ext = C.convertArray2Hexa(lines_ext)
    lines_ext = T.splitConnexity(lines_ext)
    DTW._distance2Walls(lines_ext, lineb,type='ortho',loc='nodes')
    dmax = -1e10
    line_ext = None; line_in = None
    if len(lines_ext) != 2:
        raise ValueError("DRAGON: more than two curves defined to generate the TRI mesh.Case not taken into account.")

    nol_ext = -1
    for l0 in range(len(lines_ext)):
        d0 = C.getMaxValue(lines_ext[l0],'TurbulentDistance')
        if d0 > dmax:
            dmax = d0; line_ext = lines_ext[l0]; nol_ext=l0
    if nol_ext == 0:
        line_in = lines_ext[1]
    else:
        line_in = lines_ext[0]
    return [line_in,line_ext]

def generateTriMeshBetweenContours__(lineb, surfp, hmin=1e-6,remeshBorders=False):
    extFaces0 = P.exteriorFaces(surfp)
    extFaces = T.splitConnexity(extFaces0)
    h0 = 0
    lmax = 0
    for ef in extFaces:
        lenloc = D.getLength(ef)
        if lenloc > lmax:
            lmax = lenloc
            npts = C.getNPts(ef)
            h0 = lmax/npts

    [line_in,line_ext] = orderExteriorEdges__(surfp, lineb)

    D0 = DTW.distance2Walls(line_ext, line_in, loc='nodes', type='ortho')
    D0 = C.getMinValue(D0,'TurbulentDistance')
    N0 = int(0.5*D0/h0)+1
    distrib = G.cart((0,0,0),(h0,1,1),(N0,1,1))
    distrib = G.enforcePlusX(distrib,hmin,N0//2,N0*3//2)
    tri_in = G.surfaceWalk(surfp, lineb, distrib, constraints=[], niter=50,
                           alphaRef=180., check=0, toldist=1.e-6)
    tri_in = C.convertArray2Tetra(tri_in); tri_in = G.close(tri_in)
    [line_in,line_ext0] = orderExteriorEdges__(tri_in, lineb)

    Internal._rmNodesFromName(line_ext, 'FlowSol*')
    Internal._rmNodesFromName(line_ext0, 'FlowSol*')

    contour = T.join([line_ext0,line_ext])
    tri_ext = G.tetraMesher(contour)

    # check qu on a bien triangule la zone entre les deux ??
    # par distance des frontieres exterieures aux lignes en entree
    ext = P.exteriorFaces(tri_ext)
    DTW._distance2Walls(line_ext, tri_ext, loc='nodes', type='ortho')
    if C.getMinValue(line_ext,'TurbulentDistance')>toldist:
        line_ext0 = T.reorder(line_ext0,(-1,))
        #line_ext = T.reorder(line_ext,(-1,))
        Internal._rmNodesFromName(line_ext, 'FlowSol*')

        contour = T.join([line_ext0,line_ext])
        tri_ext = G.tetraMesher(contour)

    # projection sur la surface de depart...
    #...sans bouger les frontieres...
    _projectOrthoWithConstraints__(tri_ext, surfp, contour)
    #
    if remeshBorders:
        extFaces0 = P.exteriorFaces(tri_ext)
        extFaces = T.splitConnexity(extFaces0)
        density = 0
        lmax = 0
        for ef in extFaces:
            lenloc = D.getLength(ef)
            if lenloc > lmax:
                lmax = lenloc
                npts = C.getNPts(ef)
                density = lmax/npts

        tri_ext = G.mmgs(tri_ext,ridgeAngle=45.,
                         hmin=0., hmax=density, hausd=0.01,
                         optim=0, fixedConstraints=extFaces0,
                         sizeConstraints=extFaces0)
    tri = T.join(tri_in,tri_ext)
    return tri

# after project Ortho : recover borders of surface mesh
def _projectOrthoWithConstraints__(trimesh, surfp, lines_ext):
    hook = C.createHook(trimesh, function='nodes')

    nodes = C.identifyNodes(hook, lines_ext)
    T._projectOrtho(trimesh, surfp)
    for indbc in range(len(nodes)):
        indv = nodes[indbc]-1
        coord = C.getValue(lines_ext, 'GridCoordinates', indbc)
        C.setValue(trimesh, 'CoordinateX', indv, coord[0])
        C.setValue(trimesh, 'CoordinateY', indv, coord[1])
        C.setValue(trimesh, 'CoordinateZ', indv, coord[2])
    C.freeHook(hook)
    return None

def getMatchingPoint__(line_in_h, line_periom_h):
    npts = C.getNPts(line_periom_h)
    nptsi = C.getNPts(line_in_h)
    # INFLOW HUB
    pt11 = C.getValue(line_in_h,'GridCoordinates',0)
    pt12 = C.getValue(line_in_h,'GridCoordinates',nptsi-1)

    pt21 = C.getValue(line_periom_h,'GridCoordinates',0)
    pt22 = C.getValue(line_periom_h,'GridCoordinates',npts-1)

    dist11_21 = math.sqrt((pt11[0]-pt21[0])**2+(pt11[1]-pt21[1])**2+(pt11[2]-pt21[2])**2)
    dist11_22 = math.sqrt((pt11[0]-pt22[0])**2+(pt11[1]-pt22[1])**2+(pt11[2]-pt22[2])**2)

    if dist11_21<toldist or dist11_22<toldist: ptA = pt11
    else:
        dist21_21 = math.sqrt((pt21[0]-pt21[0])**2+(pt21[1]-pt21[1])**2+(pt21[2]-pt21[2])**2)
        dist21_22 = math.sqrt((pt21[0]-pt22[0])**2+(pt21[1]-pt22[1])**2+(pt21[2]-pt22[2])**2)
        if dist21_21<toldist or dist21_22<toldist: ptA = pt21
        else:
            raise ValueError("createDragonMesh: getMatchingPoint__: no valid point found. Please contact the support.")
            return None
    return ptA


def _mirror(vtree, axis, CoordZero):
    GCNode = Internal.getNodeFromNameAndType(vtree,'GridCoordinates','GridCoordinates_t')
    CoordNode=Internal.getNodeFromNameAndType(GCNode,'Coordinate'+axis,'DataArray_t')
    Internal.setValue(CoordNode,-Internal.getValue(CoordNode)+2.*CoordZero)
    return None

# Geometry must be lying in the Z+ plane (Y+, X+ resp)
def createDragonMeshWithSymPlanes(s, dictOfParams={}, check=False, directory_tmp_files='./'):
    locmin = C.getMinValue(s,'Coordinate'+dictOfParams['sym_plane_axis'])
    locmax = C.getMaxValue(s,'Coordinate'+dictOfParams['sym_plane_axis'])
    print('locmin=',locmin); print('locmax=',locmax)
    nbSymPlanes = dictOfParams["nb_sym_planes"]

    # Temporarily changing the geometry to work with positive coordinates
    if nbSymPlanes==2:
        doMirror=False
        if locmin <- SYM_PLANE_TOL: doMirror=True

        if doMirror:
            print('Working on mirror geometry')
            _mirror(s,dictOfParams['sym_plane_axis'],0.)
            locmin = C.getMinValue(s,'Coordinate'+dictOfParams['sym_plane_axis'])
            locmax = C.getMaxValue(s,'Coordinate'+dictOfParams['sym_plane_axis'])

    # Last parameters of the DRAGON mesh

    if nbSymPlanes == 1: dictOfParams['sym_plane_values'] = [locmin]
    elif nbSymPlanes == 2: dictOfParams['sym_plane_values'] = [locmin,locmax]

    t = createDragonMesh0(s, dictOfParams, check=check,directory_tmp_files=directory_tmp_files)

    if nbSymPlanes == 2:
        if doMirror: _mirror(t, dictOfParams['sym_plane_axis'],0.)
    return t

def createDragonMesh(s, dictOfParams={}, check=False, directory_tmp_files='./'):
    if 'topology' not in dictOfParams:
        print("Warning: createDragonMesh: topology not defined in dictOfParams. Set to farfield.")
        dictOfParams['topology'] = 'farfield'
    topo = dictOfParams['topology']

    if topo == 'blade_in_channel':
        return createDragonMeshForBladeInChannel(s, dictOfParams, check=check, directory_tmp_files=directory_tmp_files)

    elif topo == 'farfield':
        return createDragonMesh0(s, dictOfParams, check=check,directory_tmp_files=directory_tmp_files)

    elif topo == 'sym':
        return createDragonMeshWithSymPlanes(s, dictOfParams, check=check,directory_tmp_files=directory_tmp_files)

    else:
        return createDragonMesh0(s, dictOfParams, check=check, directory_tmp_files=directory_tmp_files)
