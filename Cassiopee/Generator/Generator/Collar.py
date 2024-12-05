"""Collar grid generation module. Extension of Generator.
"""
try: range = xrange
except: pass

from . import Generator as G
__version__ = G.__version__

#=============================================================================
# Python Interface to create collar grids
#=============================================================================
try:
    import KCore
    from . import Generator as G
    import Intersector as XOR
    import Converter as C
    import Post as P
    import Transform as T
    from . import SurfaceWalk as SW
    import Dist2Walls as DTW
    import Geom as D
    import math
except:
    raise ImportError("Collar: requires Converter, Generator, Post, Transform, Dist2Walls, Intersector modules.")

#==============================================================================
# Return the list of created collar grids and the list of BCWall ranges as
# [[z1,range11,range12,...],[z2,range21,...],...]
# =============================================================================
def createCollarMesh__(s1, s2, distribj, distribk, niterj, niterk, ext, 
                       alphaRef, btype, contour, constraints1, constraints2, 
                       toldist):
    hjmin = distribj[1][0,1] 
    remap = 1 # 1: remap l'edge ...
    vars0 = 'x,y,z'
    if contour != []: contour[0] = vars0; remap = 0
    if constraints1 != []:
        for c1 in range(len(constraints1)): constraints1[c1][0] = vars0
    if constraints2 != []:
        for c2 in range(len(constraints2)): constraints2[c2][0] = vars0

    distribk[0] = vars0; distribj[0] = vars0;
    for nos1 in range(len(s1)): s1[nos1][0] = vars0
    for nos2 in range(len(s2)): s2[nos2][0] = vars0
    s1 = C.convertArray2Tetra(s1); s2 = C.convertArray2Tetra(s2)
    if isinstance(s1[0],list): s1 = T.join(s1)
    if isinstance(s2[0],list): s2 = T.join(s2)
    s1 = G.close(s1, toldist); s2 = G.close(s2, toldist)
    s0 = booleanSurface__(s1, s2, btype, toldist)
    if s0 == []: return []
    if contour == []: edges = booleanSurfaceEdges__(s1, s2, toldist)
    else: edges = [C.convertBAR2Struct(contour)]
    infos = [] # 0: TFI, 1: extrusion
    # Doit on determiner les contraintes ? oui si definies par [] en entree
    compCons1 = 0; compCons2 = 0

    for noe in range(len(edges)):
        edge = edges[noe]              
        # Remap l'edge si construit automatiquement par booleanSurfaceEdges
        if remap == 1: edge = remapIntersectionEdge__(edge,s1,s2,hjmin,constraints1+constraints2,alphaRef,toldist)
        # Creation de la collar mesh
        if btype == 'union': infos += buildUnionCollarMesh__(edge, s0, s1, s2, distribj, distribk,
                                                             constraints1, constraints2,
                                                             niterj, niterk, alphaRef, toldist)

        else: infos += buildDifferenceCollarMesh__(edge, s0, s1, s2, distribj, distribk, constraints1, constraints2,
                                                   niterj, niterk, ext+1, alphaRef, toldist)
    return infos

#-------------------------------------------------------------------------------
# split la zone a dans la direction i aux pts ou des contraintes sont appliquees
#-------------------------------------------------------------------------------
def getSplitZones__(a, constraints,toldist):
    ni = a[2]; nj = a[3]; nk = a[4]; ninj = ni*nj
    indices = []
    # determination des pts de split : ce sont les pts avec contraintes
    for e1 in constraints:
        ptS1 = C.getValue(e1,0); ptE1 = C.getValue(e1,e1[1].shape[1]-1)
        ind1,ds1 = D.getNearestPointIndex(a, (ptS1[0],ptS1[1],ptS1[2]))
        ind2,ds2 = D.getNearestPointIndex(a, (ptE1[0],ptE1[1],ptE1[2]))
        if ds1 < toldist:
            k1 = ind1/ninj; j1 = (ind1-k1*ninj)/ni; i1 = ind1-j1*ni-k1*ninj
            indices.append(i1+1)
        elif ds2 < toldist:
            k2 = ind2/ninj; j2 = (ind2-k2*ninj)/ni; i2 = ind2-j2*ni-k2*ninj
            indices.append(i2+1)

    indices.sort()
    pool = []
    if indices == []: return [a],0
    i = 0; firstsplit = 1
    if indices[0] != 1: z0 = T.subzone(a,(1,1,1),(indices[0],nj,nk)); firstsplit = 0
    while i < len(indices)-1:
        pool.append(T.subzone(a,(indices[i],1,1),(indices[i+1],nj,nk)))
        i += 1
    # fin
    if indices[i] != ni: 
        if firstsplit == 1: pool = [z0]+pool; pool += [T.subzone(a,(indices[i],1,1),(ni,nj,nk))]
        else: 
            pool += [T.join([z0,T.subzone(a,(indices[i],1,1),(ni,nj,nk))])]
    else:
        pool = [z0]+pool
    return pool

#-----------------------------------------------------------------------------
# cas ou le contour initial est obtenu par intersection geometrique entre les
# deux surfaces s1 et s2
# IN : constraints : liste des lignes contraintes a respecter
# IN : hjmin : hauteur de la 1ere maille en j a comparer avec la longueur des mailles des edges 
#-----------------------------------------------------------------------------
def remapIntersectionEdge__(edge, s1, s2, hjmin, constraints, alphaRef, toldist):
    if constraints == []:
        net = max(10,edge[1].shape[1]//2)
        db = G.cart((0,0,0), (1./(net-1.),1,1), (net,1,1) )
        edge = G.map(edge, db)
        return edge
    else:
        edges = getSplitZones__(edge,constraints,toldist)
        for noe in range(len(edges)):
            e = edges[noe]
            net = max(10,e[1].shape[1]//2)
            le = D.getLength(e)
            if le < hjmin:
                db = G.cart((0,0,0), (1.,1.,1.), (2,1,1) )
                e = G.map(e, db)
                edges[noe] = e
            else:
                db = G.cart((0,0,0), (1./(net-1.),1,1), (net,1,1) )
                e = G.map(e, db)
                edges[noe] = e
        edge = T.join(edges)        
        return edge

#-----------------------------------------------------------------------------
# angle min et max a la jonction entre s1 et s2 partageant la meme ligne i=1
#-----------------------------------------------------------------------------
def getExtremaAnglesAtJunction__(edge, s1,s2):
    ni = s1[2]; nj = s1[3]
    s = T.join([s1,s2]); s = C.convertArray2Tetra(s) 

    alp = D.getSharpestAngle(s); #alp = C.addVars([s, alp])
    alp0 = alp[1][0]    
    alpha = C.initVars(edge,'alpha',0.); alpha = C.extractVars(alpha,['alpha']); alpt = alpha[1][0]
    xt = C.extractVars(edge,['x'])[1][0]
    yt = C.extractVars(edge,['y'])[1][0]
    zt = C.extractVars(edge,['z'])[1][0]
    for ind in range(len(xt)):
        (index, d2) = D.getNearestPointIndex(s, (xt[ind],yt[ind],zt[ind]))
        alpt[ind] = alp0[index]
    alpha[1][0] = alpt
    alphamin = C.getMinValue(alpha,'alpha')
    alphamax = C.getMaxValue(alpha,'alpha')
    return alpha, alphamin, alphamax

#-----------------------------------------------------------------------------
# return the list of zones defining the collar grids at the junction s1 and s2
# where s1 and s2 are in difference
# and list of ranges1 and ranges2: indices of BCWall of each zone
# as [[z1,win1,win2,...],[z2,win1,win2]]
# nremapPts: nb de pts de remap de l edge initial
#-----------------------------------------------------------------------------
def buildDifferenceCollarMesh__(edge, s0, s1, s2, distribj, distribk, constraints1, constraints2, 
                                niterj, niterk, ext, alphaRef, toldist):
    check = 1
    s2 = T.reorder(s2,(-1,))
    edge0 = C.convertBAR2Struct(edge); ne = edge0[2]
    distribi = G.cart((0,0,0), (1./(ne-1),1,1), (ne,1,1) )
    res1 = orderContourForDifferenceSurface__(edge0, s1, s2)
    if res1 == 0: edge1 = edge0; edge1opp = T.reorder(edge0,(-1,2,3))
    else: edge1opp = edge0; edge1 = T.reorder(edge0,(-1,2,3))

    #===================================================
    # on uniformise la distribj en grossier
    # Extraction de la distribution en i
    sourcej = D.getCurvilinearAbscissa(distribj)
    sourcej = C.initVars(sourcej, '{CoordinateX}={s}')
    sourcej = C.initVars(sourcej, '{CoordinateY}=0')
    sourcej = C.initVars(sourcej, '{CoordinateZ}=0')
    sourcej = C.rmVars(sourcej, 's')
    #
    sourcek = D.getCurvilinearAbscissa(distribk)
    sourcek = C.initVars(sourcek, '{CoordinateX}={s}')
    sourcek = C.initVars(sourcek, '{CoordinateY}=0')
    sourcek = C.initVars(sourcek, '{CoordinateZ}=0')
    sourcek = C.rmVars(sourcek, 's')
    # distribution reguliere
    nis = 21
    sourceC = G.cart((0,0,0),(1./(nis-1),1,1),(nis,1,1))
    distribj = G.map(distribj,sourceC,1)
    distribk = G.map(distribk,sourceC,1)
    #
    #===================================================
    walkin1 = SW.surfaceWalk__([s1], edge1, distribj, constraints1, niterj, alphaRef, check, toldist)
    walkin2 = SW.surfaceWalk__([s2], edge1opp, distribj, constraints2, niterj, alphaRef, check, toldist)
    ni1 = walkin1[2]; nj1 = walkin1[3]; nj2 = walkin2[3]   
    if nj1 == 1: raise ValueError("collarMesh: surfaceWalk failed on 1st surface.")      
    if nj2 == 1: raise ValueError("collarMesh: surfaceWalk failed on 2nd surface.")
    if ext > 1:
        walkout1 = SW.surfaceWalk__([s1], edge1opp, distribj, constraints1, niterj, alphaRef, check, toldist)
        walkout2 = SW.surfaceWalk__([s2], edge1, distribj, constraints2, niterj, alphaRef, check, toldist)
        njo1 = walkout1[3]; njo2 = walkout2[3]
        if njo1 == 1: raise ValueError("collarMesh: surfaceWalk for extension failed on 1st surface.")      
        if njo2 == 1: raise ValueError("collarMesh: surfaceWalk for extension failed on 2nd surface.") 
        #
        # Extrusion en k
        [walkin1,walkout1] = G.addNormalLayers([walkin1,walkout1],distribk,niter=niterk,check=1)
        [walkin2,walkout2] = G.addNormalLayers([walkin2,walkout2],distribk,niter=niterk,check=1)
        #
        # Remap walkin1 et walkin2
        nk1 = walkout1[4]; nk2 = walkout2[4]
        walkin1 = G.map(walkin1,sourcej,dir=2); walkin2 = G.map(walkin2,sourcej,dir=2)
        walkout1 = G.map(walkout1,sourcej,dir=2); walkout1 = T.subzone(walkout1,(1,1,1),(ni1,ext+1,nk1))  
        walkout2 = G.map(walkout2,sourcej,dir=2); walkout2 = T.subzone(walkout2,(1,1,1),(ni1,ext+1,nk2))
        #
        # Join
        z1 = T.join(walkin1,walkout1); z2 = T.join(walkin2,walkout2)
        z1 = T.reorder(z1,(2,-1,3)); z2 = T.reorder(z2,(2,-1,3))
        #        
        # determination des parois       
        imin1 = walkin1[2]; jmin1 = walkin1[3]; imin2 = walkin2[2]; jmin2 = walkin2[3]
        krange1 = [1,imin1,1,jmin1,1,1]; krange2 = [1,imin2,1,jmin2,1,1]        
        return [[z1,krange1],[z2,krange2]] 
    else:
        [walkin1,walkin2] = G.addNormalLayers([walkin1,walkin2],distribk,niter=niterk,check=1)
        # Remap walkin1 et walkin2
        walkin1 = G.map(walkin1,sourcej,dir=2); walkin2 = G.map(walkin2,sourcej,dir=2)
        z = T.join([walkin1,walkin2]); z = T.reorder(z,(2,-1,3))
        # Map en k
        z = G.map(z,sourcek,3)
        imin1 = z[2]; jmin1 = z[3]; krange = [1,imin1,1,jmin1,1,1]
        return [[z,krange]] 

#------------------------------------------------------------------------------
# builds the collar mesh in case of union
# return the list of zones as [[z1,win1,win2,...],[z2,win1,win2]]
#------------------------------------------------------------------------------
def buildUnionCollarMesh__(edge, s0, s1, s2, distribj, distribk,
                           constraints1, constraints2, 
                           niterj, niterk, alphaRef, toldist):
    check = 1
    #===================================================
    # on uniformise la distribj en grossier
    # Extraction de la distribution en i
    sourcej = D.getCurvilinearAbscissa(distribj)
    sourcej = C.initVars(sourcej, '{CoordinateX}={s}')
    sourcej = C.initVars(sourcej, '{CoordinateY}=0')
    sourcej = C.initVars(sourcej, '{CoordinateZ}=0')
    sourcej = C.rmVars(sourcej, 's')
    #
    sourcek = D.getCurvilinearAbscissa(distribk)
    sourcek = C.initVars(sourcek, '{CoordinateX}={s}')
    sourcek = C.initVars(sourcek, '{CoordinateY}=0')
    sourcek = C.initVars(sourcek, '{CoordinateZ}=0')
    sourcek = C.rmVars(sourcek, 's')
    # distribution reguliere
    nis = 21
    sourceC = G.cart((0,0,0),(1./(nis-1),1,1),(nis,1,1))
    distribj = G.map(distribj,sourceC,1)
    distribk = G.map(distribk,sourceC,1)
    #===================================================
    #
    nj = distribj[2]; nk = distribk[2]
    edge0 = C.convertBAR2Struct(edge); ne = edge0[1].shape[1]
    res1 = orderContourForUnionSurface__(edge0, s1,s2)
    if res1==0: edge1 = edge0; edge1opp = T.reorder(edge0,(-1,2,3))
    else: edge1opp = edge0; edge1 = T.reorder(edge0,(-1,2,3))
    walkin2 = SW.surfaceWalk__([s2], edge1opp, distribj, constraints2, niterj, alphaRef, check, toldist)
    nj2 = walkin2[3]
    if nj2 == 1: raise ValueError("collarMesh: surfaceWalk failed on 2nd surface.")     
    #
    walkin1 = SW.surfaceWalk__([s1], edge1, distribj, constraints1, niterj, alphaRef, check, toldist)        
    nj1 = walkin1[3]
    if nj1 == 1: raise ValueError("collarMesh: surfaceWalk failed on 1st surface.")     
    #
    # determine the angle at the junction located at edge between s1 and s2    
    alpha,alphamin,alphamax = getExtremaAnglesAtJunction__(edge1,walkin1,T.reorder(walkin2,(-1,2,3)))
    ldk = D.getLength(distribk)
    # Remap walkin1 et walkin2
    if alphamin < 120.  : # TFI si 180-60 deg 
        if alphamax < 120.: z1 = generateCollarVolumeMesh0__(walkin1,walkin2,sourcej,toldist)
        else: z1 = generateCollarVolumeMesh1__(walkin2,walkin1,sourcej,alpha,toldist)
        return [[z1,[1,1,1,z1[3],1,z1[4]],[1,z1[2],1,z1[3],1,1]]]

    else: # angle > 120. : extrusion
        z = T.merge([walkin1,walkin2]); z = G.close(z)
        z = G.addNormalLayers(z,distribk,niter=niterk,check=1)
        info = []
        if isinstance(z[0],list):
            for z0 in z:
                z0 = G.map(z0,sourcek,dir=3)
                info.append([z0,[1,z0[2],1,z0[3],1,1]])
        else: 
            z = G.map(z,sourcek,dir=3)
            info.append([z,[1,z[2],1,z[3],1,1]])
        return info

#-----------------------------------------------------------------------------
# construit la surface booleenne a partir des 2 corps eventuellement ouverts
# retourne la surface resultante et les edges faisant la jonction entre les 2 surfaces
#-----------------------------------------------------------------------------
def booleanSurface__(s1,s2,btype,toldist):
    closed = 1  # suppose que les surfaces s1 et s2 sont deja fermees
    if closed == 1:
        if btype == 'union': s0 = XOR.booleanUnion(s1,s2)
        else: s0 = XOR.booleanMinus(s1,s2)
        s0 = G.close(s0,toldist)
        return s0

    else:
        # fermeture des corps
        surfs1 = []; surfs2 = []
        extEdges1 = []
        try:
            extEdges1 = P.exteriorFaces(s1)
            extEdges1 = T.splitConnexity(extEdges1)
            for noe in range(len(extEdges1)):
                p = G.fittingPlaster(extEdges1[noe], bumpFactor=0.)
                b = G.gapfixer(extEdges1[noe], p)
                surfs1.append(b)
        except: pass

        extEdges2 = []
        try:
            extEdges2 = P.exteriorFaces(s2)
            extEdges2 = T.splitConnexity(extEdges2)        
            for noe in range(len(extEdges2)):
                p = G.fittingPlaster(extEdges2[noe], bumpFactor=0.)
                b = G.gapfixer(extEdges2[noe], p)
                surfs2.append(b)           
        except: pass
        surf1 = T.join([s1]+surfs1); surf2 = T.join([s2]+surfs2)
        surf1 = T.reorder(surf1,(1,)); surf2 = T.reorder(surf2,(1,))
        #
        if btype == 'union': s0 = XOR.booleanUnion(surf1,surf2)
        else: s0 = XOR.booleanMinus(surf1,surf2)
        s0 = G.close(s0,toldist)
        return s0

#=============================================================================
def joinTBranches__(edges,toldist):
    nedges = len(edges)
    tag = [0]*nedges; n = 0
    lenmax = 0.
    noestart = 0
    for noe in range(nedges):
        edges[noe] = C.convertBAR2Struct(edges[noe])
        len0 = D.getLength(edges[noe])
        if len0 > lenmax: lenmax = len0; noestart = noe
    tag[noestart] = 1; edgestart = edges[noestart]

    ptS = C.getValue(edgestart,0)
    count = 1
    while count < nedges:
        ptE = C.getValue(edgestart,edgestart[2]-1)
        if ptE == ptS:  ptE = C.getValue(edgestart,0)
        voisins = []
        for noe in range(nedges):
            if tag[noe] == 0:
                ind,ds = D.getNearestPointIndex(edges[noe], (ptE[0],ptE[1],ptE[2]))
                if ds < toldist: voisins.append(noe); tag[noe] = 1
        # chercher le voisin valide
        count += len(voisins)
        if len(voisins) == 1: noev = voisins[0]; edgestart = T.join(edgestart,edges[noev])
        else:
            alpref = 180.; noref = 0
            for noev in voisins:
                alp = D.getSharpestAngle(edges[noev])
                alpmax = C.getMaxValue(alp,'alpha')
                if abs(alpmax-180.) < alpref: alpref = abs(alpmax-180.); noref = noev
            edgestart = T.join(edgestart,edges[noref])
    return edgestart
#-----------------------------------------------------------------------------
# determine the contour(s) edgesOut at junction between s1 and s2, resulting in
# booleanSurf by union or difference boolean operations.
# Return a list of i-arrays
#-----------------------------------------------------------------------------
def booleanSurfaceEdges__(s1,s2, toldist):
    err = 0
    edge = XOR.intersection(s1,s2)
    edge = G.close(edge,toldist)    
    edges = T.splitConnexity(edge)
    for noe in range(len(edges)):
        edges0 = T.splitTBranches(edges[noe])
        if len(edges0)>1: edges[noe] = joinTBranches__(edges0,toldist)
        edges[noe] = C.convertBAR2Struct(edges[noe])
    if err == 1:
        print('collarMesh: intersection between surfaces contains T-Branches. Intersection might be wrong. Check edge.plt file')
        C.convertArrays2File(edges, "edge.plt")
    return edges
#------------------------------------------------------------------------------
# Retourne la liste des contraintes sur une des surfaces initiales
# faisant un angle fort avec l edge a partir duquel la collar grid est
# construite
# il faut que edge soit discretise de la meme maniere que boolSurf
# cad soit le contour d intersection entre les 2 surfaces
#------------------------------------------------------------------------------
def getSurfaceConstraints__(iniSurf, boolSurf, edge, toldist, alphaRef=30.):
    constraints = []
    projtol = 1.e-5
    if len(edge) == 5: edgens = C.convertArray2Tetra(edge); edgens = G.close(edgens)
    else: edgens = G.close(edge)
    ne1 = edgens[1].shape[1];
    alphaS = D.getSharpestAngle(iniSurf)[1][0,:]
    x1 = C.extractVars(edgens,['x'])[1][0,:]
    y1 = C.extractVars(edgens,['y'])[1][0,:]
    z1 = C.extractVars(edgens,['z'])[1][0,:]
    alpha = C.initVars(edgens,'alpha',180.); alpha = C.extractVars(alpha,['alpha'])
    alphaMin = 360.; alphaMax = 0.
    for i in range(ne1):
        (ind, dist2) = D.getNearestPointIndex(iniSurf, (x1[i],y1[i],z1[i]))
        alp0 = alphaS[ind]
        if alp0 < alphaMin: alphaMin = alp0
        if alp0 > alphaMax: alphaMax = alp0
    if alphaMax-180. < alphaRef and 180-alphaMin < alphaRef: return []# pas de singularite dans le contour booleen
    # determination des contraintes : tri des edges candidats
    e0 = C.convertBAR2Struct(edgens)
    edges10 = P.sharpEdges(boolSurf,alphaRef=alphaRef)    
    edges10 = T.splitSharpEdges(edges10,alphaRef=alphaRef)
    edges10 = T.splitTBranches(edges10)
    # ne prendre que les contraintes issues de la surface courante iniSurf
    edges11 = []
    for noe1 in range(len(edges10)):
        dist1  = DTW.distance2Walls([edges10[noe1]],[iniSurf], type='ortho',loc='nodes')[0]
        if C.getMaxValue(dist1,'TurbulentDistance') < toldist:
            edges11.append(C.convertBAR2Struct(edges10[noe1]))
    edges1 = []; xp = []; yp = []; zp = []
    for e1 in edges11:
        x1 = C.extractVars(e1,['x'])[1][0,:]
        y1 = C.extractVars(e1,['y'])[1][0,:]
        z1 = C.extractVars(e1,['z'])[1][0,:]

        dist1 = DTW.distance2Walls([e1],[e0], type='ortho',loc='nodes')[0]
        if C.getMinValue(dist1,'TurbulentDistance') < toldist and  C.getMaxValue(dist1,'TurbulentDistance') > toldist:
            dist1 = dist1[1][0,:]; ne1 = e1[1].shape[1]
            #recherche du 1er point de distance non nulle
            if dist1[0] < toldist:
                xp.append(x1[0]); yp.append(y1[0]); zp.append(z1[0])
                edges1.append(e1)
            else:
                istart = 1; iend = 1
                found = 0
                for i in range(ne1-1):
                    if dist1[i] < toldist and dist1[i+1] > toldist and found == 0:
                        istart = i+1; found = 1
                        xp.append(x1[i]); yp.append(y1[i]); zp.append(z1[i])
                        edges1.append(T.subzone(e1,(istart,1,1),(iend,1,1)))
                        break
                    elif dist1[i] > toldist and dist1[i+1] < toldist and found == 0:
                        iend = i+2; found = 1
                        xp.append(x1[i+1]); yp.append(y1[i+1]); zp.append(z1[i+1])
                        edges1.append(T.subzone(e1,(istart,1,1),(iend,1,1)))
                        break
                if found == 0:
                    iend = ne1
                    if dist1[ne1-1] < toldist:
                        xp.append(x1[iend-1]); yp.append(y1[iend-1]); zp.append(z1[iend-1])
                        edges1.append(e1)

    xt0 = C.extractVars(e0,['x'])[1][0,:]
    yt0 = C.extractVars(e0,['y'])[1][0,:]
    zt0 = C.extractVars(e0,['z'])[1][0,:]                        
    for noe1 in range(len(edges1)):
        e1 = edges1[noe1]
        dist1 = DTW.distance2Walls([e1],[e0], type='ortho',loc='nodes')[0]
        nmatch = 0
        for d in dist1[1][0,:]:
            if d < toldist: nmatch+= 1
        if C.getMinValue(dist1,'TurbulentDistance') < toldist and C.getMaxValue(dist1,'TurbulentDistance') > toldist and nmatch == 1:
            xA = xp[noe1]; yA = yp[noe1]; zA = zp[noe1]
            ind1,d1 = D.getNearestPointIndex(e0, (xA,yA,zA))
            ind2,d2 = D.getNearestPointIndex(e1, (xA,yA,zA))

            if d1 < toldist and d2 < toldist:
                xt1 = C.extractVars(e1,['x'])[1][0,:]
                yt1 = C.extractVars(e1,['y'])[1][0,:]
                zt1 = C.extractVars(e1,['z'])[1][0,:]
                if ind1 < len(xt0)-1: ind1p = ind1+1
                else: ind1p = ind1-1
                if ind2 < len(xt1)-1: ind2p = ind2+1
                else: ind2p = ind2-1
                xB = xt0[ind1p]; yB = yt0[ind1p]; zB = zt0[ind1p] 
                xC = xt1[ind2p]; yC = yt1[ind2p]; zC = zt1[ind2p]
                # calcul de l angle entre les deux contours au pt d intersection
                xAB = xB-xA; yAB = yB-yA; zAB = zB-zA
                xAC = xC-xA; yAC = yC-yA; zAC = zC-zA
                nAB = math.sqrt(xAB**2+yAB**2+zAB**2)
                nAC = math.sqrt(xAC**2+yAC**2+zAC**2)
                ps = (xAB*xAC+yAB*yAC+zAB*zAC)/(nAB*nAC) # cos de l angle
                if abs(ps) < math.cos(alphaRef*math.pi/180.):  constraints.append(e1)
    return constraints

#=============================================================================
# IN: c: contour a partir duquel etendre
# IN: surfaces: liste des surfaces sur lesquelles on calcule la normale
#                 pour la direction d extension
# IN: projSurf: surface de projection de l extension : autre surface de la collar
# IN: toldist: tolerance de close
# IN: niter: nb d iterations de lissage
# OUT: coords: coordonnees de l extension (i,j) array
# OUT: contour: contour correspondant a l extension projete sur la surface projSurf
#=============================================================================
def createBFExtension(c, surfaces, projSurf, dhj, toldist, niter=0):
    for nos in range(len(surfaces)):
        if len(surfaces[nos]) == 5: surfaces[nos] = C.convertArray2Hexa(surfaces[nos])
    surfaces = G.close(surfaces)
    c0 = C.convertBAR2Struct(c)
    vars = c0[0]
    d0 = DTW.distance2Walls([c0],projSurf, type='ortho',loc='centers')[0]; d0 = C.center2Node(d0)
    ind0 = SW.getMinDistIndex__(d0)
    if ind0 != 0 : c0 = T.reorder(c0,(-1,2,3)) 

    c2 = T.projectOrtho(c0,surfaces); c2 = G.close(c2,tol=toldist)
    imax = c2[1].shape[1]; jmax = dhj[1].shape[1]
    loop = 0
    dx = c2[1][0,0]-c2[1][0,imax-1]
    dy = c2[1][1,0]-c2[1][1,imax-1]
    dz = c2[1][2,0]-c2[1][2,imax-1]
    if ( abs(dx) < toldist and abs(dy) < toldist and abs(dz) <toldist ): loop = 1

    # computes normals to s1
    normales = []
    for s1 in surfaces:
        n1 = G.getSmoothNormalMap(s1,niter=0)
        n1 = C.addVars([s1,n1])
        normales.append(n1)
    c2p = C.copy(c2)
    coords = C.arrayS(vars,imax,jmax,1)
    coords[1][0,0:imax] = c2p[1][0,:]
    coords[1][1,0:imax] = c2p[1][1,:]
    coords[1][2,0:imax] = c2p[1][2,:]
    eta = P.extractMesh(normales,c2p)
    eta = C.extractVars(eta,['sx','sy','sz'])

    for j1 in range(1,jmax):
        hloc = dhj[1][0,j1]-dhj[1][0,j1-1]
        istart = j1*imax; iend = istart+imax
        #ht = G.getLocalStepFactor__(c2p)
        #ht[1][0,:] = hloc*ht[1][0,:]
        eta2 = C.copy(eta)
#         eta2[1][0,:] = ht[1][0,:]* eta[1][0,:]
#         eta2[1][1,:] = ht[1][0,:]* eta[1][1,:]
#         eta2[1][2,:] = ht[1][0,:]* eta[1][2,:]
        eta2[1][0,:] = hloc* eta[1][0,:]
        eta2[1][1,:] = hloc* eta[1][1,:]
        eta2[1][2,:] = hloc* eta[1][2,:]
        c2n = T.deform(c2p, eta2)
        coords[1][0,istart:iend] = c2n[1][0,:]
        coords[1][1,istart:iend] = c2n[1][1,:]
        coords[1][2,istart:iend] = c2n[1][2,:]
        c2p = c2n
    return projectExtensionOnSurface(coords,projSurf)

#=============================================================================
# determination de l ordonnancement des contours si necessaire
# projection de eta sur surf1, si distance faible OK, sinon reorder du contour
#=============================================================================
def reorderContourForSurfaces__(c0, surf1, niter, toldist):
    # 1. calcul de eta
    c00 = C.convertArray2Hexa(c0)
    alp = D.getSharpestAngle(c00); alp = C.center2Node(alp)
    alp = C.initVars(alp,'{alpha}=abs(180-{alpha})')
    indstart = SW.getMinDistIndex__(alp)
    if indstart == c00[1].shape[1]-1: indstart = indstart-1
    x = c00[1][0,indstart]; y = c00[1][1,indstart]; z = c00[1][2,indstart]
    ksix = c00[1][0,indstart+1]-x; ksiy = c00[1][1,indstart+1]-y; ksiz = c00[1][2,indstart+1]-z
    #
    n1 = []; n2 = []
    for s1 in surf1: n1.append(G.getSmoothNormalMap(s1,niter))
    n1 = C.addVars([surf1,n1])
    pt0 = C.array('x,y,z',1,1,1); pt0[1][0,0]=x; pt0[1][1,0]=y; pt0[1][2,0]=z
    pt1 = T.projectOrtho(pt0,surf1)
    x1 = pt1[1][0,0]; y1 = pt1[1][1,0]; z1 = pt1[1][2,0]
    [n1x,n1y,n1z] = P.extractPoint(n1,(x1,y1,z1))
    etax = n1y*ksiz-n1z*ksiy; etay = n1z*ksix-n1x*ksiz; etaz = n1x*ksiy-n1y*ksix
    x1p=x1+etax; y1p=y1+etay; z1p=z1+etaz
    pt0 = C.array('x,y,z',1,1,1); pt0[1][0,0]=x1p; pt0[1][1,0]=y1p; pt0[1][2,0]=z1p
    pt1 = T.projectOrtho(pt0,surf1)
    dx = (pt1[1][0,0]-pt0[1][0,0]); dy = (pt1[1][1,0]-pt0[1][1,0]); dz = (pt1[1][2,0]-pt0[1][2,0])
    d1 =  dx*dx+dy*dy+dz*dz
    # comparaison avec l oppose de c0
    c0r = T.reorder(c0,(-1,2,3)); c00 = C.convertArray2Hexa(c0r)
    alp = D.getSharpestAngle(c00); alp = C.center2Node(alp)
    alp = C.initVars(alp,'alpha',initAlpha__,['alpha'])
    indstart = SW.getMinDistIndex__(alp)
    if indstart == c00[1].shape[1]-1: indstart = indstart-1
    x = c00[1][0,indstart]; y = c00[1][1,indstart]; z = c00[1][2,indstart]
    ksix = c00[1][0,indstart+1]-x; ksiy = c00[1][1,indstart+1]-y; ksiz = c00[1][2,indstart+1]-z
    #
    n1 = []; n2 = []
    for s1 in surf1: n1.append(G.getSmoothNormalMap(s1,niter))
    n1 = C.addVars([surf1,n1])
    pt0 = C.array('x,y,z',1,1,1); pt0[1][0,0]=x; pt0[1][1,0]=y; pt0[1][2,0]=z
    pt1 = T.projectOrtho(pt0,surf1)
    x1 = pt1[1][0,0]; y1 = pt1[1][1,0]; z1 = pt1[1][2,0]
    [n1x,n1y,n1z] = P.extractPoint(n1,(x1,y1,z1))
    etax = n1y*ksiz-n1z*ksiy; etay = n1z*ksix-n1x*ksiz; etaz = n1x*ksiy-n1y*ksix
    x1p=x1+etax; y1p=y1+etay; z1p=z1+etaz
    pt0 = C.array('x,y,z',1,1,1); pt0[1][0,0]=x1p; pt0[1][1,0]=y1p; pt0[1][2,0]=z1p
    pt1 = T.projectOrtho(pt0,surf1)
    dx = (pt1[1][0,0]-pt0[1][0,0]); dy = (pt1[1][1,0]-pt0[1][1,0]); dz = (pt1[1][2,0]-pt0[1][2,0])
    d2 =  dx*dx+dy*dy+dz*dz
    if d1 < d2: return c0
    else: return c0r

#===============================================================================
# Generation de la grille collar en TFI a partir des surfaces
# IN : calpha : contour avec angle alpha entre les deux surfaces
# IN : sourcej : distribution dans la direction d extrusion surfacique a mapper
#===============================================================================
def generateCollarVolumeMesh0__(surf1, surf2, sourcej, toldist):
    surf1 = G.map(surf1,sourcej,dir=2)
    surf2 = G.map(surf2,sourcej,dir=2)
    ni1 = surf1[2]; nj1 = surf1[3]
    ni2 = surf2[2]; nj2 = surf2[3]
    # boucle : couper en 2 les 2 surfaces
    xmin1 = surf1[1][0,0]; ymin1 = surf1[1][1,0]; zmin1 = surf1[1][2,0]
    xmax1 = surf1[1][0,ni1-1]; ymax1 = surf1[1][1,ni1-1]; zmax1 = surf1[1][2,ni1-1]
    cas1 = 0 
    mesh = []
    if abs(xmin1-xmax1) < toldist and abs(ymin1-ymax1)< toldist and abs(zmin1-zmax1)< toldist:
        surf1o = []; surf2o = []
        # verifier que ds le meme ordre
        cas1 = 1
        x1 = surf1[1][0,1]; y1 = surf1[1][1,1]; z1 = surf1[1][2,1]
        x2 = surf2[1][0,1]; y2 = surf2[1][1,1]; z2 = surf2[1][2,1]
        if abs(x1-x2)>toldist or abs(y1-y2)>toldist or abs(z1-z2)>toldist:
            surf2 = T.reorder(surf2,(-1,-2,3))# pour que les 2 coincident en i = 1
        isplit = ni1/2
        surf1_1 = T.subzone(surf1,(1,1,1),(isplit,nj1,1))
        surf1_2 = T.subzone(surf1,(isplit,1,1),(ni1,nj1,1))
        surf2_1 = T.subzone(surf2,(1,1,1),(isplit,nj2,1))
        surf2_2 = T.subzone(surf2,(isplit,1,1),(ni2,nj2,1))
        surf1o = [surf1_1,surf1_2]; surf2o = [surf2_1,surf2_2]
        m = buildTFIMeshType0__(surf1o[0],surf2o[0]); mesh.append(m)
        m = buildTFIMeshType0__(surf1o[1],surf2o[1]); mesh.append(m)
    else:
        surf1o = []; surf2o = []
        # verifier que ds le meme ordre
        x1 = surf1[1][0,1]; y1 = surf1[1][1,1]; z1 = surf1[1][2,1]
        x2 = surf2[1][0,1]; y2 = surf2[1][1,1]; z2 = surf2[1][2,1]
        surf1o = [surf1]
        if abs(x1-x2)>toldist or abs(y1-y2)>toldist or abs(z1-z2)>toldist:
            surf2o = [T.reorder(surf2,(-1,-2,3))]# pour que les 2 coincident en i = 1
        else: surf2o = [surf2]
        m = buildTFIMeshType0__(surf1o[0],surf2o[0]); mesh.append(m)

    mesh = T.join(mesh)
    if cas1 == 1: mesh = T.reorder(mesh,(2,-1,3))
    return mesh
#=============================================================================
# Generation de la grille collar en TFI a partir des surfaces
# IN : calpha : contour avec angle alpha entre les deux surfaces
#=============================================================================
def generateCollarVolumeMesh1__(surf1, surf2, distribj, calpha, toldist):
    mesh = []
    ni1 = surf1[2]; nj1 = surf1[3]
    ni2 = surf2[2]; nj2 = surf2[3]
    # calcul du vecteur directeur pour decaler r1 par rapport a r2
    if ni1 != ni2: raise ValueError("generateCollarVolumeMesh1__: i direction must be common.")
    xmin1 = surf1[1][0,0]; ymin1 = surf1[1][1,0]; zmin1 = surf1[1][2,0]
    xmax1 = surf1[1][0,ni1-1]; ymax1 = surf1[1][1,ni1-1]; zmax1 = surf1[1][2,ni1-1]
    isloop = 0
    if abs(xmin1-xmax1) < toldist and abs(ymin1-ymax1)< toldist and abs(zmin1-zmax1)< toldist: 
        isloop = 1
        # verifier que ds le meme ordre
        x1 = surf1[1][0,1]; y1 = surf1[1][1,1]; z1 = surf1[1][2,1]
        x2 = surf2[1][0,1]; y2 = surf2[1][1,1]; z2 = surf2[1][2,1]
        if abs(x1-x2)>toldist or abs(y1-y2)>toldist or abs(z1-z2)>toldist:
            surf2 = T.reorder(surf2,(-1,-2,3))# pour que les 2 coincident en i = 1
    #
    #================================================
    # 1er vecteur : tangente a surf1
    #================================================
    n1 = C.array('sx,sy,sz',ni1,1,1); vect1 = n1[1]  
    coords1 = surf1[1]
    for i1 in range(ni1):
        istart = i1; iend = i1+ni1#(nj1-1)*ni1 
        dx = coords1[0,iend]-coords1[0,istart]
        dy = coords1[1,iend]-coords1[1,istart]
        dz = coords1[2,iend]-coords1[2,istart]
        l1inv = 1./math.sqrt(dx*dx+dy*dy+dz*dz)
        vect1[0,i1] = dx*l1inv
        vect1[1,i1] = dy*l1inv
        vect1[2,i1] = dz*l1inv

    #=============================================
    # 2nd vecteur : normale a surf2
    #=============================================
    surfLoc = T.subzone(surf2,(1,1,1),(ni2,2,1))
    surfu = C.convertArray2Hexa(surfLoc); surfu = G.close(surfu)
    surfu = T.reorder(surfu, (1,))
    indicesU = KCore.indiceStruct2Unstr2([surf2], surfu, 1.e-14)[0]
    nu = G.getSmoothNormalMap(surfu)  
    c2 = T.subzone(surfLoc,(1,1,1),(ni2,1,1))
    n2 = C.array('sx,sy,sz',ni2,1,1)
    for ind in range(ni2):
        indU  = indicesU[ind]
        n2[1][:,ind] = nu[1][:,indU]

    #=============================================
    # Vecteur fonction de alpha
    #=============================================
    alp1 = C.extractVars(calpha,['alpha'])[1]
    n = C.array('sx,sy,sz',ni1,1,1);
    vn = n[1]; vn1 = n1[1]; vn2 = n2[1]
    dh1 = T.subzone(surf1,(1,1,1),(1,surf1[3],1))  
    dh1 = T.reorder(dh1,(-2,1,3))
    l1 = D.getLength(dh1)
    alphamax = 90.
    for i in range(ni1):
        alp0 = alp1[0,i]
        if alp0 < alphamax: 
            for j in range(3): vn[j,i]=l1*vn1[j,i]       
        else:
            for j in range(3): vn[j,i]=l1*vn2[j,i]
    # Petit lissage ...
    surfLoc = C.addVars([n,c2])
    surfu = C.convertArray2Hexa(surfLoc); surfu = G.close(surfu)
    indicesU = KCore.indiceStruct2Unstr2([surfLoc], surfu, 1.e-14)[0]
    nitemax = 20
    for ite in range(nitemax):
        surfu = C.node2Center(surfu)
        surfu = C.center2Node(surfu)
    surfu = C.extractVars(surfu,['sx','sy','sz'])
    for ind in range(n[1].shape[1]):
        indU  = indicesU[ind]
        n[1][:,ind] = surfu[1][:,indU]

    c2p = T.deform(c2,n)
    c2p = C.extractVars(c2p,['x','y','z'])
    #=============================================            
    # Construction de surf1opp
    #=============================================
    surf1opp = C.array('x,y,z',ni1,2,1)
    for ind in range(ni1):
        surf1opp[1][:,ind]= c2[1][:,ind]

    for ind in range(ni1,2*ni1):
        indo = ind-ni1
        surf1opp[1][:,ind]= c2p[1][:,indo]

    dh1 = T.subzone(surf1,(1,1,1),(1,nj1,1))   
    dh1 = T.reorder(dh1,(-2,1,3)) 
    sourcej = D.getCurvilinearAbscissa(dh1)
    sourcej = C.initVars(sourcej, '{CoordinateX}={s}')
    sourcej = C.initVars(sourcej, '{CoordinateY}=0')
    sourcej = C.initVars(sourcej, '{CoordinateZ}=0')
    sourcej = C.rmVars(sourcej, 's')
    surf1opp = G.map(surf1opp,sourcej,dir=2)

    #================================================
    # Creation de surf2opp
    #================================================
    surf2opp = C.array('x,y,z',ni2,2,1)
    for i in range(ni1):
        ind1 = i + (nj1-1)*ni1
        surf2opp[1][:,i] = surf1[1][:,ind1]
        surf2opp[1][:,i+ni1] = surf1opp[1][:,ind1]
    dh2 = T.subzone(surf2,(1,1,1),(1,surf2[3],1))   
    dh2 = T.reorder(dh2,(-2,1,3)) 
    sourcej2 = D.getCurvilinearAbscissa(dh2)
    sourcej2 = C.initVars(sourcej2,'{CoordinateX}={s}')
    sourcej2 = C.initVars(sourcej2,'{CoordinateY}=0')
    sourcej2 = C.initVars(sourcej2,'{CoordinateZ}=0')
    sourcej2 = C.rmVars(sourcej2,'s')
    surf2opp = G.map(surf2opp,sourcej2,dir=2)

    #================================================
    # boucle : couper en 2 les 2 surfaces pour la TFI
    #================================================
    if isloop == 1:
        surf1o = []; surf2o = []
        ni1 = surf1[2]; nj1 = surf1[3]; isplit = ni1/2
        ni2 = surf2[2]; nj2 = surf2[3]
        surf1_1 = T.subzone(surf1,(1,1,1),(isplit,surf1[3],1))
        surf1_2 = T.subzone(surf1,(isplit,1,1),(ni1,nj1,1))
        surf2_1 = T.subzone(surf2,(1,1,1),(isplit,nj2,1))
        surf2_2 = T.subzone(surf2,(isplit,1,1),(ni2,nj2,1))
        surf1o = [surf1_1,surf1_2]; surf2o = [surf2_1,surf2_2]
        nj1opp = surf1opp[3]; nj2opp = surf2opp[3]
        surf1opp_1 = T.subzone(surf1opp,(1,1,1),(isplit,nj1opp,1))
        surf1opp_2 = T.subzone(surf1opp,(isplit,1,1),(ni1,nj1opp,1))
        surf2opp_1 = T.subzone(surf2opp,(1,1,1),(isplit,nj2opp,1))
        surf2opp_2 = T.subzone(surf2opp,(isplit,1,1),(ni2,nj2opp,1))
        #
        surf1o = [surf1_1,surf1_2]; surf1oppo=[surf1opp_1,surf1opp_2]
        surf2o = [surf2_1,surf2_2]; surf2oppo=[surf2opp_1,surf2opp_2]
    else:
        surf1o = [surf1]; surf1oppo = [surf1opp]
        surf2o = [surf2]; surf2oppo = [surf2opp]        
    #
    #=============================================
    # creation des fenetres i = 1 et i = imax
    #=============================================
    for nos1 in range(len(surf1o)):
        r1 = surf1o[nos1]; r1opp = surf1oppo[nos1]
        r2 = surf2o[nos1]; r2opp = surf2oppo[nos1]

        # Calcul de r5 
        a11 = T.subzone(r1,(1,1,1),(1,r1[3],r1[4]))
        a21 = T.subzone(r1opp,(1,1,1),(1,r1opp[3],r1opp[4]))
        a31 = T.subzone(r2,(1,1,1),(1,r2[3],r2[4]))
        a41 = T.subzone(r2opp,(1,1,1),(1,r2opp[3],r2opp[4]))
        A = [a11,a21,a31,a41]
        A = T.reorder(A,(3,1,2))
        r5 = G.TFI(A)
        # Calcul de r6
        a11 = T.subzone(r1,(r1[2],1,1),(r1[2],r1[3],r1[4]))
        a21 = T.subzone(r1opp,(r1opp[2],1,1),(r1opp[2],r1opp[3],r1opp[4]))
        a31 = T.subzone(r2,(r2[2],1,1),(r2[2],r2[3],r2[4]))
        a41 = T.subzone(r2opp,(r2opp[2],1,1),(r2opp[2],r2opp[3],r2opp[4]))
        A = [a11,a21,a31,a41]; A = T.reorder(A,(3,1,2))
        r6 = G.TFI(A)
        out = [r1,r1opp,r2,r2opp,r5,r6]       
        m = G.TFI(out)
        mesh.append(m)
    mesh = T.join(mesh)
    if isloop == 1: mesh = T.reorder(mesh,(2,-1,3))
    mesh = G.map(mesh,distribj,dir=1)
    mesh = G.map(mesh,distribj,dir=3)
    return mesh

def buildTFIMeshType0__(r1,r2):
    if len(r1) != 5 or len(r2) != 5: raise ValueError("buildTFIMesh0__: only for structured arrays.")
    # calcul du vecteur directeur pour decaler r1 par rapport a r2
    ni1 = r1[2]; nj1 = r1[3]; ni2 = r2[2]; nj2  = r2[3]
    if ni1 != ni2: raise ValueError("buildTFIMesh0__: i direction must be common.")
    #=============================================
    # translation de r1 par rapport a r2 -> r1opp1
    # idem pour r2
    #=============================================
    vect = C.array('sx,sy,sz',ni1,nj1,1); vect1 = vect[1]  
    coords2 = r2[1]; r1opp = C.copy(r1)

    for j2 in range(nj2-1):
        for i2 in range(ni2):
            iend = i2+j2*ni2; istart = iend + ni2
            dx = coords2[0,iend]-coords2[0,istart]
            dy = coords2[1,iend]-coords2[1,istart]
            dz = coords2[2,iend]-coords2[2,istart]
            for j1 in range(nj1):
                vect1[0,i2+j1*ni1] = dx
                vect1[1,i2+j1*ni1] = dy
                vect1[2,i2+j1*ni1] = dz        
        r1opp = T.deform(r1opp,vect)

    vect = C.array('sx,sy,sz',ni2,nj2,1); vect2 = vect[1]  
    coords1 = r1[1]; r2opp = C.copy(r2)
    for j1 in range(nj1-1):
        for i1 in range(ni1):
            istart = i1+j1*ni1; iend = istart + ni1
            dx = coords1[0,iend]-coords1[0,istart]
            dy = coords1[1,iend]-coords1[1,istart]
            dz = coords1[2,iend]-coords1[2,istart]
            for j2 in range(nj2):
                vect2[0,i1+j2*ni2] = dx
                vect2[1,i1+j2*ni2] = dy
                vect2[2,i1+j2*ni2] = dz
        r2opp = T.deform(r2opp,vect)

    #=============================================
    # creation des fenetres i = 1 et i = imax
    # calcul de r4 : 
    #=============================================
    r1opp = C.extractVars(r1opp,['x','y','z'])
    r2opp = C.extractVars(r2opp,['x','y','z'])
    a11 = T.subzone(r1,(1,1,1),(1,r1[3],r1[4]))
    a21 = T.subzone(r1opp,(1,1,1),(1,r1opp[3],r1opp[4]))
    a31 = T.subzone(r2,(1,1,1),(1,r2[3],r2[4]))
    a41 = T.subzone(r2opp,(1,1,1),(1,r2opp[3],r2opp[4]))
    A = [a11,a21,a31,a41]; A = T.reorder(A,(3,1,2))
    r5 = G.TFI(A)
    a11 = T.subzone(r1,(ni1,1,1),(ni1,r1[3],r1[4]))
    a21 = T.subzone(r1opp,(ni1,1,1),(ni1,r1opp[3],r1opp[4]))
    a31 = T.subzone(r2,(ni2,1,1),(ni1,r2[3],r2[4]))
    a41 = T.subzone(r2opp,(ni2,1,1),(ni1,r2opp[3],r2opp[4]))
    A = [a11,a21,a31,a41]; A = T.reorder(A,(3,1,2))
    r6 = G.TFI(A)
    m = G.TFI([r1,r1opp,r2,r2opp,r5,r6])
    return m


#=============================================================================
# retourne la surface projetee pour la premiere ligne + le contour associe
#=============================================================================
def projectExtensionOnSurface(extSurf,surfaces):
    # recuperation de l arete initiale
    nie = extSurf[2]; edge1 = T.subzone(extSurf,(1,1,1),(nie,1,1))
    # determination de ind1 : 0 ou nie
    d1 = DTW.distance2Walls([edge1], surfaces, type='ortho',loc='centers')[0]; d1 = C.center2Node(d1)
    ind1 = SW.getMinDistIndex__(d1)
    x1 = extSurf[1][0,ind1]; y1=extSurf[1][1,ind1]; z1=extSurf[1][2,ind1]
    # projection et preservation du premier pt
    contour1 = T.subzone(extSurf,(ind1+1,1,1),(ind1+1,extSurf[3],1))
    contour1 = T.projectOrtho(contour1,surfaces)
    extSurf = T.patch(extSurf, contour1, (ind1+1,1,1))
    extSurf[1][0,ind1] = x1; extSurf[1][1,ind1] = y1; extSurf[1][2,ind1] = z1
    contour1 = T.subzone(extSurf,(ind1+1,1,1),(ind1+1,extSurf[3],1))
    contour = T.reorder(contour1,(-2,1,3))
    return extSurf, contour

#-----------------------------------------------------------------------------
# order contour for walk on the difference surface
# IN : edge to reorder
# IN : s1, s2 : surfaces for which the boolean union is applied
# OUT : 0 : not reordered for s1, 1 : to be reordered for s1
#-----------------------------------------------------------------------------
def orderContourForUnionSurface__(edge, s1, s2):
    xt = C.extractVars(edge,['x'])[1]; xA = xt[0,0]; xB = xt[0,1]
    yt = C.extractVars(edge,['y'])[1]; yA = yt[0,0]; yB = yt[0,1]
    zt = C.extractVars(edge,['z'])[1]; zA = zt[0,0]; zB = zt[0,1]
    ksix = xB-xA; ksiy = yB-yA; ksiz = zB-zA
    ksin = 1./math.sqrt(ksix*ksix+ksiy*ksiy+ksiz*ksiz)
    ksix = ksix*ksin;  ksiy = ksiy*ksin; ksiz = ksiz*ksin
    #
    ptA = C.array('x,y,z',1,1,1); ptA[1][0,0]=xA; ptA[1][1,0]=yA; ptA[1][2,0]=zA
    n1 = G.getSmoothNormalMap(s1,niter=0); n1 = C.normalize(n1, ['sx','sy','sz'])
    pt1 = T.projectOrtho(ptA,[s1])
    x1 = pt1[1][0,0]; y1 = pt1[1][1,0]; z1 = pt1[1][2,0]
    n1 = C.addVars([s1,n1])
    [n1x,n1y,n1z] = P.extractPoint([n1],(x1,y1,z1))
    etax1 = n1y*ksiz-n1z*ksiy; etay1 = n1z*ksix-n1x*ksiz; etaz1 = n1x*ksiy-n1y*ksix
    normeta1 = max(1.e-12,math.sqrt(etax1*etax1+etay1*etay1+etaz1*etaz1)); normeta1 = 1./(normeta1)
    etax1 = etax1*normeta1; etay1 = etay1*normeta1; etaz1 = etaz1*normeta1
    #    
    n2 = G.getSmoothNormalMap(s2,niter=0); n2 = C.normalize(n2, ['sx','sy','sz'])
    pt2 = T.projectOrtho(ptA,[s2])
    x2 = pt2[1][0,0]; y2 = pt2[1][1,0]; z2 = pt2[1][2,0]
    n2 = C.addVars([s2,n2])
    [n2x,n2y,n2z] = P.extractPoint([n2],(x2,y2,z2))
    normn2 = max(1.e-12,math.sqrt(n2x*n2x+n2y*n2y+n2z*n2z)); normn2 = 1./normn2
    n2x = n2x*normn2; n2y = n2y*normn2; n2z = n2z*normn2
    ps = etax1*n2x+etay1*n2y+etaz1*n2z
    #print(ps)
    if ps < -0.25: return 1 # reordonner pour s1
    else: return 0

#-----------------------------------------------------------------------------
# order contour for walk on the difference surface
# IN : edge to reorder
# IN : s1, s2 : surfaces for which the difference is applied
# OUT : 0 : not reordered for s1, 1 : to be reordered for s1
#-----------------------------------------------------------------------------
def orderContourForDifferenceSurface__(edge, s1, s2):
    xt = C.extractVars(edge,['x'])[1]; xA = xt[0,0]; xB = xt[0,1]
    yt = C.extractVars(edge,['y'])[1]; yA = yt[0,0]; yB = yt[0,1]
    zt = C.extractVars(edge,['z'])[1]; zA = zt[0,0]; zB = zt[0,1]
    ksix = xB-xA; ksiy = yB-yA; ksiz = zB-zA
    ksin = 1./math.sqrt(ksix*ksix+ksiy*ksiy+ksiz*ksiz)
    ksix = ksix*ksin;  ksiy = ksiy*ksin; ksiz = ksiz*ksin
    #
    ptA = C.array('x,y,z',1,1,1); ptA[1][0,0]=xA; ptA[1][1,0]=yA; ptA[1][2,0]=zA
    n1 = G.getSmoothNormalMap(s1,niter=0); n1 = C.normalize(n1, ['sx','sy','sz'])
    pt1 = T.projectOrtho(ptA,[s1])
    x1 = pt1[1][0,0]; y1 = pt1[1][1,0]; z1 = pt1[1][2,0]
    n1 = C.addVars([s1,n1])
    [n1x,n1y,n1z] = P.extractPoint([n1],(x1,y1,z1))
    etax1 = n1y*ksiz-n1z*ksiy; etay1 = n1z*ksix-n1x*ksiz; etaz1 = n1x*ksiy-n1y*ksix
    #
    n2 = G.getSmoothNormalMap(s2,niter=0); n2 = C.normalize(n2, ['sx','sy','sz'])
    n2 = C.addVars([s2,n2])
    pt2 = T.projectOrtho(ptA,[s2])
    x2 = pt2[1][0,0]; y2 = pt2[1][1,0]; z2 = pt2[1][2,0]
    [n2x,n2y,n2z] = P.extractPoint([n2],(x2,y2,z2))
    ps = etax1*n2x + etay1*n2y + etaz1*n2z
    if ps < -0.8: return 0
    else: return 1
