import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import Connector.PyTree as X
import Post.PyTree as P
import Transform.PyTree as T
import Dist2Walls.PyTree
import Intersector.PyTree as XOR
import Dist2Walls.PyTree as DTW
import Converter.Internal as Internal
import Intersector.PyTree as XOR
import numpy, math, copy

SYM_PLANE_TOL = 1.e-7

def createDragonMesh0(body, dictOfParams={}, check=False, directory_tmp_files='./'):
    if 'sym_plane_axis' in dictOfParams: sym = dictOfParams['sym_plane_axis'] # symmetry plane at x,y,z constant
    else: sym = None
    
    if 'sym_plane_values' in dictOfParams: Locsyms= dictOfParams['sym_plane_values']
    else: Locsyms = []
    
    if 'NParts' in dictOfParams: NParts = dictOfParams['NParts']
    else: NParts = 0

    if 'holeFactor' in dictOfParams: holeFactor = dictOfParams['holeFactor']
    else: holeFactor = 1.25

    if 'hWall' not in dictOfParams:
        hWall = 1e-6
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

    if 'h_ext' in dictOfParams:  snear = dictOfParams["h_ext"] # -1 pour calcul auto
    else: snear = -1
    
    if 'dfar' in dictOfParams: dfar = dictOfParams["dfar"] # -1 pour calcul auto
    else: dfar = -1
    
    if 'h_ext_factor' in dictOfParams: snearFactor = dictOfParams["h_ext_factor"]
    else: snearFactor = 1


    nbsyms = len(Locsyms)
    if nbsyms>0:
        locmin = min(Locsyms)
        locmax = max(Locsyms)

        if nbsyms == 2: print('Configuration with two symmetry planes.')
        elif nbsyms == 1: print('Configuration with one symmetry plane.')
    else: print('Configuration in free air.')

    
    # pour un seul plan de symetrie, on symetrise la geometrie pour assurer que l octree se decoupe sur ce plan
    if nbsyms == 1:
        if sym == 'X': syms = T.symetrize(body, (Locsyms[0],0.,0.), (0,1,0), (0,0,1))
        if sym == 'Y': syms = T.symetrize(body, (0.,Locsyms[0],0.), (1,0,0), (0,0,1))
        if sym == 'Z': syms = T.symetrize(body, (0.,0.,Locsyms[0]), (1,0,0), (0,1,0))
        syms[0]='syms'
        body = T.join(body,syms); G._close(body,tol=4.e-5)

    T._reorderAll(body, dir=1)

    # on travaille sur des geometries triangulees 
    body = C.deleteFlowSolutions__(body)
    body = C.convertArray2Tetra(body)
    body = T.join(body); body = G.close(body)
    T._reorderAll(body, dir=1)

    XMax_surf = C.getMaxValue(body,'CoordinateX')

    # Calcul automatique du snear (si snear < 0), il peut etre module avec le param snearFactor
    sizemax = 0.
    box = G.bbox(body)
    sizemax = max(sizemax,abs(box[3]-box[0]),abs(box[4]-box[1]),abs(box[5]-box[2]));print('sizemax=',sizemax)
    if (snear < 0):
        G._getNormalMap(body)
        C._magnitude(body, ['centers:sx','centers:sy','centers:sz'])
        meanarea = C.getMeanValue(body,'centers:sMagnitude')
        minarea = C.getMinValue(body,'centers:sMagnitude')
        maxarea = C.getMaxValue(body,'centers:sMagnitude')
        print('surface area min/max/mean=',minarea,maxarea,meanarea)
        # snear.append(2.*(vmin-1)*meanarea**(1/2.))
        snear = snearFactor*meanarea**(1./2.); print('snear=', snear)
        if (dfar < 0): dfar = 10.*sizemax ; print('dfar=', dfar)
 
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
    T._reorderAll(s,dir=1)
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
            a1 = Internal.getNodeFromType(a1,'Zone_t'); a1[0]='octree%d'%i
            npts = int(max(1,maxnpts/2**i));print(i,npts)
            T._addkplane(a1,N=npts)
            T._scale(a1, factor=(1.,1.,locmax/npts))
            #T._scale(a1, factor=(1.,1.,minsize*2**i))
            zmin = C.getMinValue(a1,'CoordinateZ')
            T._translate(a1, (0.,0.,-zmin))
            tsym[2][1][2].append(a1)

        Internal._rmNodesByName(tsym,'FlowSolution*')
        if check: C.convertPyTree2File(tsym,directory_tmp_files+'toto.cgns')
        o = Internal.getZones(T.join(tsym))[0]

    # conversion en NGon et conformisation de l octree
    t = C.newPyTree(['Base',o])

    print('Conformizing...')
    t = C.convertArray2NGon(t) ; G._close(t)
    t = C.conformizeNGon(t) ; G._close(t)

    # Extract blanking surface
    print('Extract blanking surface (%s)...'%delta)
    DTW._distance2Walls(t, s, type='ortho', loc='centers', signed=1)

    s2 = P.isoSurfMC(t, 'centers:TurbulentDistance', delta)
    if nbsyms == 2:
        s2 = [T.splitConnexity(s2)[0]]
        sexts = P.exteriorFaces(s2)
        sexts = T.splitConnexity(sexts)
        if check: C.convertPyTree2File(sexts,directory_tmp_files+'tmp.cgns')
        for sext in sexts: 
            sext = G.tetraMesher(sext, grading=0.2, maxh=0.5*snear, algo=1)
            s2.append(sext)
    if check: C.convertPyTree2File(s2, directory_tmp_files+'deltaSurface.cgns')

    # Blanking with delta surface
    print('Blanking octree width delta...')
    BM = numpy.array([[1]])
    t = X.blankCells(t, [s2], BM, blankingType='center_in', delta=1.e-12, XRaydim1=xrayDim1, XRaydim2=xrayDim2)

    if check: C.convertPyTree2File(t, directory_tmp_files+'blanked.cgns')

    # Selectionne les cellules cellN=1
    print('Selecting...')
    t = P.selectCells2(t, 'centers:cellN', strict=1)
    if check: C.convertPyTree2File(t, directory_tmp_files+'main.cgns')

    # Recupere les surfaces exterieures
    print('Getting exterior faces...')
    ext = P.exteriorFaces(t)
    if nbsyms ==2: ext = T.splitSharpEdges(ext, alphaRef=89.)
    else: ext = T.splitConnexity(ext)
    if check: C.convertPyTree2File(ext, directory_tmp_files+'ext2Watch.cgns')

    # Selectionne la bonne frontiere (celle proche corps)
    if nbsyms!=2:
        ext2 = G.getVolumeMap(ext)
        minsize = C.getMinValue(ext2,'centers:vol')**0.5

        ext = C.newPyTree(['Base'])
        for e in ext2:
            meansize = C.getMeanValue(e,'centers:vol')**0.5
            if abs(meansize-minsize)<eps: ext[2][1][2].append(e)
        
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
    if (len(ext)>1): print("Warning: Tetra/octree border has manifolds. Tetra mesher might fail.")
    ext = ext[0]
    T._reorderAll(ext, dir=1)
    T._reorderAll(s, dir=-1)

    if check: C.convertPyTree2File([s,ext], directory_tmp_files+'tetra_surfs.cgns')
    s_in = [s, ext]

    if nbsyms>0:
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
        if nbsyms==2:
            surfsym2 = G.tetraMesher(surfsmax, grading=0.2, maxh=0.5*snear, algo=1)
            s_in += [surfsym2]

    s_in = T.join(s_in); s_in = XOR.conformUnstr(s_in, tol=0., itermax=1); T._reorderAll(s_in,dir=1)
    if check: C.convertPyTree2File(s_in, directory_tmp_files+'ext.cgns')
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
    bc2=T.splitSharpEdges(bc, alphaRef=30.)

    for zone in Internal.getByType(bc,'Zone_t')[2]:
        if (C.getMaxValue(zone,'CoordinateX') == xp and C.getMaxValue(zone,'CoordinateY') == yp):
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
        if abs(bary[0] - xm) < eps : C._addBC2Zone(m,'IN','FamilySpecified:IN',subzone=zone); #print('Xmin')
        if abs(bary[0] - xp) < eps : C._addBC2Zone(m,'OUT','FamilySpecified:OUT',subzone=zone); #print('Xmax')
        if abs(bary[1] - ym) < eps : C._addBC2Zone(m,'DOWN','FamilySpecified:DOWN',subzone=zone); #print('Ymin')
        if abs(bary[1] - yp) < eps : C._addBC2Zone(m,'TOP','FamilySpecified:TOP',subzone=zone); #print('Ymax')
        if abs(bary[2] - zm) < eps : C._addBC2Zone(m,'RIGHT','FamilySpecified:RIGHT',subzone=zone); #print('Zmin')
        if abs(bary[2] - zp) < eps : C._addBC2Zone(m,'LEFT','FamilySpecified:LEFT',subzone=zone); #print('Zmax')

    C._fillEmptyBCWith(m,'WALL', 'FamilySpecified:WALL')

    tp = C.newPyTree(['Base']); tp[2][1][2] += [m]

    print('add families')
    base = Internal.getNodeFromType(tp,'CGNSBase_t')
    #print(base)
    if nbsyms>0:
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


def createDragonMeshForBladeInChannel(aube, dictOfParams={}, check=False, directory_tmp_files='./'):
    if 'remesh_input' in dictOfParams:
        remesh_aube = dictOfParams['remesh_input']
        if remesh_aube: aube = G.mmgs(aube)

    if 'nb_blades' in dictOfParams:
        nb_aubes_fan = dictOfParams["nb_blades"]
    else:
        raise ValueError("createDragonMesh: nb_blades required for periodicity.")
        return None

    if 'z_in' in dictOfParams:
        zamont = dictOfParams['z_in']
    else:
        raise ValueError("createDragonMesh: value of z at inflow border is required.")
        return None
    if 'z_out' in dictOfParams:
        zaval = dictOfParams['z_out']
    else:
        raise ValueError("createDragonMesh: value of z at outflow border is required.")
        return None
    
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
        
    centre = (0.,0.,0.); axis = (0.,0.,1.)
    if 'rotation_axis' in dictOfParams: axis = dictOfParams['rotation_axis']
    else:
        print('Warning: createDragonMesh: rotation axis not defined. Set to default :', axis)
    if 'rotation_center' in dictOfParams: centre = dictOfParams['rotation_center']
    else:
        print('Warning: createDragonMesh: rotation center not defined. Set to default :', centre)   
    
    mesh_aube = C.convertArray2Tetra(aube)
    mesh_aube = T.reorder(mesh_aube, (-1,))
    mesh_aube = Internal.getZones(mesh_aube)[0]
    
    lines_ext_aube = P.exteriorFaces(mesh_aube)
    lines_ext_aube = T.splitConnexity(lines_ext_aube)
    line_aube_spin = lines_ext_aube[0];T._reorder(line_aube_spin,(-1,2,3))
    line_aube_cart = lines_ext_aube[1];T._reorder(line_aube_cart,(-1,2,3))

    
    dim = C.getNPts(line_aube_spin)
    zmin_aube_spin = C.getMinValue(line_aube_spin,'CoordinateZ')
    zmax_aube_spin = C.getMaxValue(line_aube_spin,'CoordinateZ')
    for i in range(dim):
        coord = C.getValue(line_aube_spin,'GridCoordinates',i)
        if coord[2]==zmin_aube_spin: ba_aube_spin = coord; ind_ba_spin = i
        if coord[2]==zmax_aube_spin: bf_aube_spin = coord; ind_bf_spin = i

    line2_spin = T.splitBAR(line_aube_spin,ind_ba_spin,ind_bf_spin)[1]
    line2_spin = D.uniformize(line2_spin, N=50, sharpAngle=40.)
    amont_spin = copy.deepcopy(ba_aube_spin); amont_spin[2] = zamont
    aval_spin = copy.deepcopy(bf_aube_spin); aval_spin[2] = zaval
    line1_spin = D.line(amont_spin,ba_aube_spin,N=50)
    line3_spin = D.line(bf_aube_spin,aval_spin,N=50)
    line_spin = T.join([line1_spin,line2_spin,line3_spin]); line_spin = C.convertBAR2Struct(line_spin)
    
    dim = C.getNPts(line_aube_cart)
    zmin_aube_cart = C.getMinValue(line_aube_cart,'CoordinateZ')
    zmax_aube_cart = C.getMaxValue(line_aube_cart,'CoordinateZ')
    for i in range(dim):
        coord = C.getValue(line_aube_cart,'GridCoordinates',i)
        if coord[2]==zmin_aube_cart: ba_aube_cart = coord; ind_ba_cart = i
        if coord[2]==zmax_aube_cart: bf_aube_cart = coord; ind_bf_cart = i
    line2_cart = T.splitBAR(line_aube_cart,ind_ba_cart,ind_bf_cart)[1]
    line2_cart = D.uniformize(line2_cart, N=50, sharpAngle=40.)
    amont_cart = copy.deepcopy(ba_aube_cart); amont_cart[2] = zamont
    aval_cart = copy.deepcopy(bf_aube_cart); aval_cart[2] = zaval
    line1_cart = D.line(amont_cart,ba_aube_cart,N=50)
    line3_cart = D.line(bf_aube_cart,aval_cart,N=50)
    line_cart = T.join([line1_cart,line2_cart,line3_cart]); line_cart = C.convertBAR2Struct(line_cart)

    line_amont = D.line(amont_spin,amont_cart,N=50)
    line_aval = D.line(aval_spin,aval_cart,N=50)

    perio = G.TFI([line_spin,line_aval,line_cart,line_amont])
    perio = T.rotate(perio, centre, axis, -0.5*360./nb_aubes_fan)

    line_amont = T.rotate(line_amont, centre, axis, -0.5*360./nb_aubes_fan)
    surf_amont = D.axisym(line_amont, centre, axis, angle=360./nb_aubes_fan, Ntheta=30); surf_amont = C.convertArray2Tetra(surf_amont); surf_amont[0]='AMONT'
    line_aval = T.rotate(line_aval, centre, axis, -0.5*360./nb_aubes_fan)
    surf_aval = D.axisym(line_aval, centre, axis, angle=360./nb_aubes_fan, Ntheta=30); surf_aval = C.convertArray2Tetra(surf_aval); surf_aval[0]='AVAL'
    line_spin = T.rotate(line_spin, centre, axis, -0.5*360./nb_aubes_fan)
    surf_spin = D.axisym(line_spin, centre, axis, angle=360./nb_aubes_fan, Ntheta=30); surf_spin[0]='SPIN'
    line_cart = T.rotate(line_cart, centre, axis, -0.5*360./nb_aubes_fan)
    surf_cart = D.axisym(line_cart, centre, axis, angle=360./nb_aubes_fan, Ntheta=30); surf_cart[0]='CART'
    surf_perio1 = perio; surf_perio1 = C.convertArray2Tetra(surf_perio1); surf_perio1[0]='PERIO1'
    surf_perio2 = T.rotate(perio, centre, axis, 360./nb_aubes_fan); surf_perio2 = C.convertArray2Tetra(surf_perio2); surf_perio2[0]='PERIO2'


    tmp = C.newPyTree(['Base',surf_spin,surf_cart,surf_amont,surf_aval,surf_perio1,surf_perio2])

    bc = P.exteriorFaces(mesh_aube)
    bcproj = T.projectOrtho(bc,[surf_spin,surf_cart])
    hook = C.createHook(mesh_aube, function='nodes')
    nodes = C.identifyNodes(hook, bc)#; print(nodes)
    for ind in nodes:
        coord = tuple(C.getValue(mesh_aube, 'GridCoordinates', ind-1))#;print(coord)
        indproj = D.getNearestPointIndex(bcproj, coord)[0]#;print indproj
        coordproj = C.getValue(bcproj, 'GridCoordinates', indproj)#;print(coordproj)
        C.setValue(mesh_aube, 'CoordinateX', ind-1, coordproj[0])
        C.setValue(mesh_aube, 'CoordinateY', ind-1, coordproj[1])
        C.setValue(mesh_aube, 'CoordinateZ', ind-1, coordproj[2])

    lines_ext_aube = P.exteriorFaces(mesh_aube)
    lines_ext_aube = T.splitConnexity(lines_ext_aube)
    line_aube_spin = lines_ext_aube[0]#;T._reorder(line_aube_spin,(-1,2,3))
    line_aube_cart = lines_ext_aube[1];T._reorder(line_aube_cart,(-1,2,3))

    line_ext_spin = P.exteriorFaces(surf_spin)
    line_ext_spin = C.convertArray2Hexa(line_ext_spin); line_ext_spin = C.convertBAR2Struct(line_ext_spin)
    T._projectOrtho(line_ext_spin,[surf_spin])
    line_ext_cart = P.exteriorFaces(surf_cart)
    line_ext_cart = C.convertArray2Hexa(line_ext_cart); line_ext_cart = C.convertBAR2Struct(line_ext_cart)
    T._projectOrtho(line_ext_cart,[surf_cart])

    mesh_spin = G.tetraMesher([line_ext_spin,line_aube_spin])
    T._projectOrtho(mesh_spin,[surf_spin])
    
    mesh_cart = G.tetraMesher([line_ext_cart,line_aube_cart])
    T._smooth(mesh_cart, eps=0.5, niter=60, type=0, fixedConstraints=[line_ext_cart,line_aube_cart], \
              projConstraints=[], delta=0.1, point=(0, 0, 0), radius=-1.)
    T._projectOrtho(mesh_cart,[surf_cart])

    surfs = T.join([mesh_spin,mesh_cart,mesh_aube]); G._close(surfs)
    surfs = T.reorder(surfs,(-1,))
    ts = C.newPyTree(['WALL','AMONT','AVAL','PERIODIC'])
    ts[2][1][2] = [surfs]
    ts[2][2][2] = [surf_amont]
    ts[2][3][2] = [surf_aval]
    ts[2][4][2] = [surf_perio1, surf_perio2]

    # input: surfaces bases ['WALL','AMONT','AVAL','PERIODIC']
    print('Generating boundary layers...')
    surfs = Internal.getZones(Internal.getNodeFromName(ts,'WALL'))[0]
    surf_amont = Internal.getZones(Internal.getNodeFromName(ts,'AMONT'))[0]
    surf_aval = Internal.getZones(Internal.getNodeFromName(ts,'AVAL'))[0]
    surf_perio1 = Internal.getZones(Internal.getNodeFromName(ts,'PERIODIC'))[0]
    surf_perio2 = Internal.getZones(Internal.getNodeFromName(ts,'PERIODIC'))[1]

    d = G.cart((0.,0.,0.), (0.1,1,1), (nlayer,1,1))
    for i in range(0,nlayer): C.setValue(d, 'CoordinateX', i, hWall*(1.-raison**i)/(1.-raison))
    lay = G.addNormalLayers(surfs, d, check=1, niterType=0, niter=smoothIter, niterK=[], 
                            smoothType=0, eps=0.4, nitLocal=3, 
                            kappaType=0, kappaS=[0.2,1.6], blanking=False, algo=0)
    lay = C.convertArray2NGon(lay, recoverBC=0)
    XOR._volumes(lay)
    volmin = C.getMinValue(lay,'centers:volumes')
    print("final volume min =",volmin)
    Internal._rmNodesByName(lay,'FlowSol*')

    bcc = P.exteriorFaces(lay)
    bcc = T.breakElements(bcc)
    bc = bcc[1] #on selectionne les quads
    bc = T.splitSharpEdges(bc, alphaRef=60.)
    bcnames = [bc[2][0],bc[3][0],bc[4][0],bc[9][0]] # a automatiser

    bc2 = []
    for node in bc:
        if node[0] not in bcnames: bc2.append(node)

    bc = T.join(bc2)
    bc = T.splitConnexity(bc)

    bcstruct = []
    for i in range(len(bc)):
        z = bc[i]
        cont = P.exteriorFaces(z)
        cont = T.splitSharpEdges(cont, alphaRef=80.)
        cont = C.convertBAR2Struct(cont)
        cont = G.TFI(cont)
        if i%2 != 0: cont = T.reorder(cont,(1,-2,3))
        bcstruct.append(cont)


    hook = C.createHook(lay, function='nodes')

    for i in range(len(bc)):
        b = bc[i]
        ymin = C.getMinValue(b, 'CoordinateY'); print(b[0],ymin)
        nodes = C.identifyNodes(hook, b)#; print(nodes)
        if ymin<0.:
            r = T.projectOrtho(b,[surf_perio1])
        if ymin>0.:
            b = bcstruct[i-1]
            b = T.projectOrtho(b,[surf_perio1])
            r = T.rotate(b, centre, axis, 360./nb_aubes_fan)
        for ind in nodes:
            coord = tuple(C.getValue(lay, 'GridCoordinates', ind-1))
            if ymin<0.: indproj = D.getNearestPointIndex(r, coord)[0]
            if ymin>0.: indproj = D.getNearestPointIndex(bcstruct[i], coord)[0]
            coordproj = C.getValue(r, 'GridCoordinates', indproj)
            C.setValue(lay, 'CoordinateX', ind-1, coordproj[0])
            C.setValue(lay, 'CoordinateY', ind-1, coordproj[1])
            C.setValue(lay, 'CoordinateZ', ind-1, coordproj[2])

    bcc = P.exteriorFaces(lay)
    bcc = T.breakElements(bcc)
    bc = bcc[1] #on selectionne les quads
    bc = P.exteriorFaces(bc)
    bc = T.splitConnexity(bc)

    #detection auto
    minx = C.getMinValue(bc,'CoordinateX')
    maxx = C.getMaxValue(bc,'CoordinateX')
    for bound in bc:
        minxloc = C.getMinValue(bound,'CoordinateX')
        maxxloc = C.getMaxValue(bound,'CoordinateX')
        if minxloc > minx and minxloc < (0.5*(maxx+minx)): line_ext_spin = bound
        if maxxloc < maxx and maxxloc > (0.5*(maxx+minx)): line_ext_cart = bound
    
    #
    line_ext_spin = T.splitSharpEdges(line_ext_spin, alphaRef=89.)
    line_ext_spin = C.convertBAR2Struct(line_ext_spin)

    line_ext_cart = T.splitSharpEdges(line_ext_cart, alphaRef=89.)
    line_ext_cart = C.convertBAR2Struct(line_ext_cart)

    line_spin_amont=[]
    line_spin_aval=[]
    line_spin_perio1=[]
    line_spin_perio2=[]
    for seg in line_ext_spin:
        for surfproj in [surf_amont,surf_aval,surf_perio1,surf_perio2]:
            segproj = T.projectOrtho(seg,surfproj)
            coord = numpy.array(C.getValue(seg, 'GridCoordinates', (2,1,1)))#; print(coord)
            coordproj = numpy.array(C.getValue(segproj, 'GridCoordinates', (2,1,1)))#; print(coordproj)
            if numpy.allclose(coord,coordproj, atol=1.e-3):
               if surfproj[0]=='AMONT':line_spin_amont.append(seg)
               if surfproj[0]=='AVAL':line_spin_aval.append(seg)
               if surfproj[0]=='PERIO1':line_spin_perio1.append(seg)
               if surfproj[0]=='PERIO2':line_spin_perio2.append(seg)


    line_spin_amont = T.join(line_spin_amont)
    line_spin_aval = T.join(line_spin_aval)
    line_spin_perio1 = T.join(line_spin_perio1)
    line_spin_perio2 = T.join(line_spin_perio2)

    toler = 5.e-2

    line_cart_amont=[]
    line_cart_aval=[]
    line_cart_perio1=[]
    line_cart_perio2=[]
    for seg in line_ext_cart:
        for surfproj in [surf_amont,surf_aval,surf_perio1,surf_perio2]:
            segproj = T.projectOrtho(seg,surfproj)
            coord1 = numpy.array(C.getValue(seg, 'GridCoordinates', (1,1,1)))
            coordproj1 = numpy.array(C.getValue(segproj, 'GridCoordinates', (1,1,1)))
            coord2 = numpy.array(C.getValue(seg, 'GridCoordinates', (0,1,1)))
            coordproj2 = numpy.array(C.getValue(segproj, 'GridCoordinates', (0,1,1)))
            if numpy.allclose(coord1,coordproj1,atol=toler) and numpy.allclose(coord2,coordproj2,atol=toler):
               if surfproj[0]=='AMONT':line_cart_amont.append(seg)#;print('AMONT')
               if surfproj[0]=='AVAL':line_cart_aval.append(seg)#;print('AVAL')
               if surfproj[0]=='PERIO1':line_cart_perio1.append(seg)#;print('PERIO1')
               if surfproj[0]=='PERIO2':line_cart_perio2.append(seg)#;print('PERIO2')

    line_cart_amont = T.join(line_cart_amont)
    line_cart_aval = T.join(line_cart_aval)
    line_cart_perio1 = T.join(line_cart_perio1)
    line_cart_perio2 = T.join(line_cart_perio2)

    ptA = C.getValue(line_spin_perio1, 'GridCoordinates', (0,1,1))
    ptB = C.getValue(line_spin_perio1, 'GridCoordinates', (1,1,1))
    ptC = C.getValue(line_cart_perio1, 'GridCoordinates', (0,1,1))
    ptD = C.getValue(line_cart_perio1, 'GridCoordinates', (1,1,1))
    ptE = C.getValue(line_spin_perio2, 'GridCoordinates', (0,1,1))
    ptF = C.getValue(line_spin_perio2, 'GridCoordinates', (1,1,1))
    ptG = C.getValue(line_cart_perio2, 'GridCoordinates', (0,1,1))
    ptH = C.getValue(line_cart_perio2, 'GridCoordinates', (1,1,1))
    linepts = [ptA,ptB,ptC,ptD,ptE,ptF,ptG,ptH]
    # print(linepts)
    z_linepts = [pt[2] for pt in linepts]
    linepts.sort(key=lambda x:x[2])
    linepts_amont = linepts[:4]
    linepts_amont.sort(key=lambda x:x[0])
    linepts_amont_spin = linepts_amont[:2]
    linepts_amont_cart = linepts_amont[-2:]
    linepts_amont_spin.sort(key=lambda x:x[1])
    linepts_amont_cart.sort(key=lambda x:x[1])
    linepts_aval = linepts[-4:]
    linepts_aval.sort(key=lambda x:x[0])
    linepts_aval_spin = linepts_aval[:2]
    linepts_aval_cart = linepts_aval[-2:]
    linepts_aval_spin.sort(key=lambda x:x[1])
    linepts_aval_cart.sort(key=lambda x:x[1])


    line_perio1_aval = D.line(linepts_aval_spin[0],linepts_aval_cart[0],N=50)
    line_perio1_amont = D.line(linepts_amont_spin[0],linepts_amont_cart[0],N=50)
    line_perio2_aval = D.line(linepts_aval_spin[1],linepts_aval_cart[1],N=50)
    line_perio2_amont = D.line(linepts_amont_spin[1],linepts_amont_cart[1],N=50)

    dh = 0.07
    Npts = 100 #100
    D.setH(line_perio1_amont, 0, dh); D.setH(line_perio1_amont, -1, dh)
    line_perio1_amont = D.enforceh(line_perio1_amont, N=Npts)
    D.setH(line_perio1_aval, 0, dh); D.setH(line_perio1_aval, -1, dh)
    line_perio1_aval = D.enforceh(line_perio1_aval, N=Npts)
    D.setH(line_perio2_amont, 0, dh); D.setH(line_perio2_amont, -1, dh)
    line_perio2_amont = D.enforceh(line_perio2_amont, N=Npts)
    D.setH(line_perio2_aval, 0, dh); D.setH(line_perio2_aval, -1, dh)
    line_perio2_aval = D.enforceh(line_perio2_aval, N=Npts)


    if check: C.convertPyTree2File([line_spin_amont,line_perio1_amont,line_cart_amont,line_perio2_amont,line_spin_aval,line_perio1_aval,line_cart_aval,line_perio2_aval,line_spin_perio1,line_cart_perio1,line_spin_perio2,line_cart_perio2],directory_tmp_files+'alllines.cgns')


    print('mesh_amont')
    mesh_amont = G.TFI([line_spin_amont,line_perio1_amont,line_cart_amont,line_perio2_amont])
    mesh_amont = C.convertArray2Tetra(mesh_amont)

    print('mesh_aval')
    mesh_aval = G.TFI([line_spin_aval,line_perio1_aval,line_cart_aval,line_perio2_aval])
    mesh_aval = C.convertArray2Tetra(mesh_aval)

    print('mesh_perio1')
    mesh_perio1 = G.TFI([line_spin_perio1,line_perio1_amont,line_cart_perio1,line_perio1_aval])
    mesh_perio1 = C.convertArray2Tetra(mesh_perio1)

    print('mesh_perio2')
    mesh_perio2 = G.TFI([line_spin_perio2,line_perio2_amont,line_cart_perio2,line_perio2_aval])
    mesh_perio2 = C.convertArray2Tetra(mesh_perio2)

    bcc = P.exteriorFaces(lay)
    bcc = T.breakElements(bcc)
    mesh_ext = T.splitConnexity(bcc[0])[1]
    surfs = [mesh_amont,mesh_aval,mesh_perio1,mesh_perio2]
    surfs = T.join(surfs)
    G._close(surfs, tol=toler)
    surfs = [surfs,mesh_ext]
    surfs = T.join(surfs)
    G._close(surfs)
    T._reorderAll(surfs)
    if check: C.convertPyTree2File(surfs, directory_tmp_files+'surfs_ext.cgns')

    mesh_ext = G.tetraMesher(surfs)
    mesh_ext = C.convertArray2NGon(mesh_ext)
    if check: C.convertPyTree2File(mesh_ext, directory_tmp_files+'mesh_ext.cgns')

    mesh_final = T.join(mesh_ext,lay); G._close(mesh_final)

    XOR._volumes(mesh_final)
    volmin = C.getMinValue(mesh_final,'centers:volumes')
    print("final volume min =",volmin)
    Internal._rmNodesByName(mesh_final,'FlowSol*')

    # Add BCs
    C._addBC2Zone(mesh_final, 'SPIN', 'FamilySpecified:SPIN', subzone=mesh_spin)
    C._addBC2Zone(mesh_final, 'CART', 'FamilySpecified:CART', subzone=mesh_cart)
    C._addBC2Zone(mesh_final, 'AUBE', 'FamilySpecified:AUBE', subzone=mesh_aube)
    mesh_final = X.connectMatchPeriodic(mesh_final, rotationCenter=[0.,0.,0.],rotationAngle=[0.,0.,360./nb_aubes_fan], tol=1e-5)
    faceList = C.getEmptyBC(mesh_final)
    C._addBC2Zone(mesh_final, 'IN', 'FamilySpecified:IN', faceList = faceList[0])
    C._addBC2Zone(mesh_final, 'OUT', 'FamilySpecified:OUT', faceList = faceList[1])
    # convert GC2BC
    pl = []
    for gc in Internal.getNodesFromType(mesh_final,'GridConnectivity_t'):
        pl.append(Internal.getValue(Internal.getNodeFromName(gc,'PointList')))
    C._addBC2Zone(mesh_final, 'PERIOR', 'FamilySpecified:PERIOR', faceList = pl[0])
    C._addBC2Zone(mesh_final, 'PERIOL', 'FamilySpecified:PERIOL', faceList = pl[1])
    Internal._rmNodesByName(mesh_final,'ZoneGridConnectivity')
    
    print('add families')
    tp = C.newPyTree(['Base']); tp[2][1][2] += [mesh_final]
    base = Internal.getNodeFromType(tp,'CGNSBase_t')
    C._addFamily2Base(base, 'IN', bndType='BCInflow')
    C._addFamily2Base(base, 'OUT', bndType='BCOutflow')
    C._addFamily2Base(base, 'SPIN', bndType='BCWall')
    C._addFamily2Base(base, 'CART', bndType='BCWall')
    C._addFamily2Base(base, 'AUBE', bndType='BCWall')

    return tp

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

    if nbSymPlanes==1: dictOfParams['sym_plane_values']=[locmin]
    elif nbSymPlanes==2: dictOfParams['sym_plane_values']=[locmin,locmax]
    
    t = createDragonMesh0(s,dictOfParams,check=check,directory_tmp_files=directory_tmp_files)

    if nbSymPlanes==2:
        if doMirror: _mirror(t, dictOfParams['sym_plane_axis'],0.)
    return t

def createDragonMesh(s, dictOfParams={}, check=False, directory_tmp_files='./'):
    if 'topology' not in dictOfParams:
        print("Warning: createDragonMesh: topology not defined in dictOfParams. Set to farfield.")
        dictOfParams['topology'] = 'farfield'
    topo = dictOfParams['topology']    

    if topo=='blade_in_channel':
        return createDragonMeshForBladeInChannel(s,dictOfParams,check=check,directory_tmp_files=directory_tmp_files)

    elif topo=='farfield':    
        return createDragonMesh0(s,dictOfParams,check=check,directory_tmp_files=directory_tmp_files)

    elif topo == 'sym':
        return createDragonMeshWithSymPlanes(s,dictOfParams,check=check,directory_tmp_files=directory_tmp_files)

    else:
        return createDragonMesh0(s,dictOfParams,check=check,directory_tmp_files=directory_tmp_files)
