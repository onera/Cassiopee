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
import numpy, math


def createDragonMesh(s, dictOfParams={}, check=False, directory_tmp_files='./'):
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
        syms = T.symetrize(s, (0.,0.,Locsyms[0]), (1,0,0), (0,1,0)); syms[0]='syms'
        s = T.join(s,syms); G._close(s,tol=4.e-5)

    s = T.reorderAll(s)

    # on travaille sur des geometries triangulees 
    s = C.deleteFlowSolutions__(s)
    s = C.convertArray2Tetra(s)
    s = T.join(s) ; s = G.close(s)
    s = T.reorder(s, (+1,))

    XMax_surf = C.getMaxValue(s,'CoordinateX')
    #C.convertPyTree2File(s,"tmp.cgns")
    # Calcul automatique du snear (si snear < 0), il peut etre module avec le param snearFactor
    sizemax = 0.
    box = G.bbox(s)
    sizemax = max(sizemax,abs(box[3]-box[0]),abs(box[4]-box[1]),abs(box[5]-box[2]));print('sizemax=',sizemax)
    if (snear < 0):
        G._getNormalMap(s)
        C._magnitude(s, ['centers:sx','centers:sy','centers:sz'])
        meanarea = C.getMeanValue(s,'centers:sMagnitude')
        minarea = C.getMinValue(s,'centers:sMagnitude')
        maxarea = C.getMaxValue(s,'centers:sMagnitude')
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
    print('Generating the boundary layers...')

    if nblayerforced != 0: nlayer = nblayerforced
    else: nlayer = int(math.log(snear/(height_factor*hWall))/math.log(raison));
    print('Number of prismatic layers = ',nlayer)
    dlayer = hWall*(1.-raison**nlayer)/(1.-raison) ;
    print('Thickness of prismatic layer: ',dlayer)
    d = G.cart((0.,0.,0.), (dlayer/nlayer,1,1), (nlayer,1,1))
    for i in range(0,nlayer): C.setValue(d, 'CoordinateX', i, hWall*(1.-raison**i)/(1.-raison))
    hLast = C.getValue(d, 'CoordinateX', nlayer-1)-C.getValue(d, 'CoordinateX', nlayer-2) ; print('hauteur de la derniere couche de prisme = ',hLast)
    lay1 = G.addNormalLayers(s, d, check=1, niter=smoothIter)
    lay.append(lay1);layy = lay[0]
    if check: C.convertPyTree2File(lay, directory_tmp_files+'layer.cgns')
    ext = P.exteriorFaces(lay[0])
    if sym is None:
        ext = T.splitConnexity(ext)
        s = ext[1]
    else:
        ext = T.splitSharpEdges(ext, alphaRef=80.)
        s = ext[-1]
    s = T.breakElements(s)[0]
    s = T.reorderAll(s)
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
            a1 = P.selectCells(o, 'abs({centers:vol}-%g)<%g'%(area,eps), strict=0)
            a1 = Internal.getNodeFromType(a1,'Zone_t'); a1[0]='octree%d'%i
            npts = int(max(1,maxnpts/2**i));print(i,npts)
            a1 = T.addkplane(a1,N=npts)
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
    t = X.blankCells(t, [s2], BM, blankingType='node_in', delta=1.e-12, XRaydim1=xrayDim1, XRaydim2=xrayDim2)

    if check: C.convertPyTree2File(t, directory_tmp_files+'blanked.cgns')

    # Selectionne les cellules cellN=1
    print('Selecting...')
    t = P.selectCells2(t, 'cellN', strict=1)
    if check: C.convertPyTree2File(t, directory_tmp_files+'main.cgns')

    # Recupere les surfaces exterieures
    print('Getting exterior faces...')
    ext = P.exteriorFaces(t)
    if nbsyms ==2: ext = T.splitSharpEdges(ext, alphaRef=89.)
    else: ext = T.splitConnexity(ext)
    if check: C.convertPyTree2File(ext, directory_tmp_files+'ext2Watch.cgns')

    # Selectionne la bonne frontiere (celle proche corps)
    ext2 = G.getVolumeMap(ext)
    minsize = C.getMinValue(ext2,'centers:vol')**0.5

    ext = C.newPyTree(['Base'])
    for e in ext2:
        meansize = C.getMeanValue(e,'centers:vol')**0.5
        if abs(meansize-minsize)<eps: ext[2][1][2].append(e)
    
    ext = T.join(ext)
    C._addBC2Zone(t, 'toto', 'BCWall', subzone=ext)
    print("triangulateBC..")
    t = XOR.triangulateBC(t, 'BCWall')
    ext = C.extractBCOfType(t, 'BCWall')[0]
    ext = T.breakElements(ext)[0]
    ext = T.reorder(ext, (+1,))

    # pour un seul plan de symetrie, on selectionne la partie du maillage qui interesse
    if nbsyms == 1:
       ext = P.selectCells(ext, '{Coordinate%s}>%f'%(sym,Locsyms[0]), strict=0)
       s = P.selectCells(s, '{Coordinate%s}>%f'%(sym,Locsyms[0]-0.0001), strict=1)
       lay = P.selectCells(lay, '{Coordinate%s}>%f'%(sym,Locsyms[0]-0.0001), strict=1)

    # C.convertPyTree2File(ext, 'tmp.cgns')
    # C.convertPyTree2File(s, 'tmp1.cgns')
    # C.convertPyTree2File(lay, 'tmp2.cgns')
    # import sys;sys.exit()

    # Remplit l'espace entre l'exterieur faces, la surface et le remplissage tri sur le(s) plan(s) de symetrie
    print('Tetra filling...')
    s = T.reorder(s, (-1,))
    if nbsyms>0:
        extlines = P.exteriorFaces(ext)
        extlines = T.splitConnexity(extlines)
        slines = P.exteriorFaces(s)
        slines = T.splitConnexity(slines)
        #print(len(extlines),len(slines))
        surfsmin = []; surfsmax = []
        for i in range(len(extlines)):
            extlineloc = C.getMeanValue(extlines[i],'Coordinate%s'%sym)
            if abs(extlineloc - locmin)<eps: surfsmin.append(extlines[i])
            if abs(extlineloc - locmax)<eps: surfsmax.append(extlines[i])
            slineloc = C.getMeanValue(slines[i],'Coordinate%s'%sym)
            if abs(slineloc - locmin)<eps: surfsmin.append(slines[i])
            if abs(slineloc - locmax)<eps: surfsmax.append(slines[i])
        surfsym1 = G.tetraMesher(surfsmin, grading=0.2, maxh=0.5*snear, algo=1)
        if nbsyms==2:
            surfsym2 = G.tetraMesher(surfsmax, grading=0.2, maxh=0.5*snear, algo=1)
            if check: C.convertPyTree2File([surfsym1,surfsym2,ext,s], directory_tmp_files+'ext.cgns')
            m = G.tetraMesher([surfsym1,surfsym2,ext,s], grading=0.2, maxh=0.5*snear, algo=1)
        if nbsyms==1:
            if check: C.convertPyTree2File([surfsym1,ext,s], directory_tmp_files+'ext.cgns')
            m = G.tetraMesher([surfsym1,ext,s], grading=0.2, maxh=0.5*snear, algo=1)
    else:
        m = G.tetraMesher([ext,s], grading=0.2, maxh=0.5*snear, algo=1)

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
        if (C.getMaxValue(zone,'CoordinateX') == xp and C.getMaxValue(zone,'CoordinateY') == yp): print(zone[0]);Internal._rmNode(bc,zone)
    for zone in Internal.getByType(bc2,'Zone_t')[2]:
        bc[2][1][2].append(zone)
    if check: C.convertPyTree2File(bc, directory_tmp_files+'exterior.cgns')

    # Add BCs
    eps = 1.e-3
    for zone in Internal.getByType(bc,'Zone_t')[2]:
        dimax = Internal.getZoneDim(zone)[1] - 1
        bary = G.barycenter(zone)
        print(zone[0],dimax,bary)
        if abs(bary[0] - xm) < eps : C._addBC2Zone(m,'IN','FamilySpecified:IN',subzone=zone); print('Xmin')
        if abs(bary[0] - xp) < eps : C._addBC2Zone(m,'OUT','FamilySpecified:OUT',subzone=zone); print('Xmax')
        if abs(bary[1] - ym) < eps : C._addBC2Zone(m,'DOWN','FamilySpecified:DOWN',subzone=zone); print('Ymin')
        if abs(bary[1] - yp) < eps : C._addBC2Zone(m,'TOP','FamilySpecified:TOP',subzone=zone); print('Ymax')
        if abs(bary[2] - zm) < eps : C._addBC2Zone(m,'RIGHT','FamilySpecified:RIGHT',subzone=zone); print('Zmin')
        if abs(bary[2] - zp) < eps : C._addBC2Zone(m,'LEFT','FamilySpecified:LEFT',subzone=zone); print('Zmax')

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
        C._addFamily2Base(base, 'RIGHT', bndType='BCSymmetryPlane')
        C._addFamily2Base(base, 'LEFT', bndType='BCSymmetryPlane')
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
