# - OversetData -
# Module for internal functions used for overset info
from . import Connector
from . import connector
import KCore.Vector as Vector
import numpy
__version__ = Connector.__version__

import Converter.GhostCells as GhostCells

try:
    import Converter.Internal as Internal
    import Converter.PyTree as C
    import Converter
    import KCore
except:
    raise ImportError("Connector.OversetData requires Converter module.")

TypesOfIBC={}
TypesOfIBC["slip_cr"]=100 # slip + curvature radius
TypesOfIBC["slip"]=0
TypesOfIBC["noslip"]=1
TypesOfIBC["Log"]=2
TypesOfIBC["Musker"]=3 #also 31
TypesOfIBC["SA"]=32
TypesOfIBC["outpress"]=4
TypesOfIBC["inj"]=5
TypesOfIBC["TBLE"]=6
TypesOfIBC["MuskerMob"]=7 #Ivan
TypesOfIBC["Pohlhausen"]=8 #Lam fp
TypesOfIBC["Thwaites"]=9 #Lam fp inclined
TypesOfIBC["Mafzal"]=10 #Musker+gradP
TypesOfIBC["TBLE_FULL"]=11 #TBLE+gradP+conv+SA
TypesOfIBC["isothermal"]=12 #isothermal: set T_wall
TypesOfIBC["heatflux"]=13 #heatflux: set q_wall
TypesOfIBC["overlap"]=14 #TBLE+gradP+conv+SA
TypesOfIBC["MuskerLinear"]=331
TypesOfIBC["SALinear"]=332
TypesOfIBC["wiremodel"]=140 #wire mesh model

# Variables IBM pour le post traitement
__PRESSURE__= 'Pressure'
__UTAU__    = 'utau'
__YPLUS__   = 'yplus'
__YPLUSIP__ = 'yplusIP'
__DENSITY__ = 'Density'
__VELOCITYX__ = 'VelocityX'
__VELOCITYY__ = 'VelocityY'
__VELOCITYZ__ = 'VelocityZ'
__KCURV__ = 'KCurv' # curvature radius for IBCs
__TEMPERATURE__ = 'Temperature'
__GRADXPRESSURE__ = 'gradxPressure'
__GRADYPRESSURE__ = 'gradyPressure'
__GRADZPRESSURE__ = 'gradzPressure'
__CONV1__ = 'conv1' # convective terms for TBLE FULL
__CONV2__ = 'conv2'

# Revert dict
IBCTypes = {}
for ibctypename in TypesOfIBC:
    ibctypeint = TypesOfIBC[ibctypename]
    IBCTypes[ibctypeint] = ibctypename

#=============================================================================
# Ajout du cellN pour une zone qui ne le contient pas. celln=1
#=============================================================================
def addCellN__(t, loc='centers', cellNName='cellN'):
    tp = Internal.copyRef(t)
    _addCellN__(tp, loc=loc, cellNName=cellNName)
    return tp

def _addCellN__(t, loc='centers', cellNName='cellN'):
    if loc == 'centers': var = loc+':'+cellNName
    else: var = cellNName
    C._addVars(t, var)
    return None

#=============================================================================
# 1. CALCUL DES INTERSECTIONS ENTRE DOMAINES
#=============================================================================
def getIntersectingDomainsAABB(t, tol=1.e-10):
    """Return the intersection list of a list of bounding boxes."""
    m = C.getFields(Internal.__GridCoordinates__, t, api=3)
    ret = Connector.getIntersectingDomainsAABB(m, tol)
    dic = {}
    zoneNames = C.getZoneNames(t, False)
    c = 0
    for i in ret:
        dic[zoneNames[c]] = [zoneNames[j] for j in i]
        c += 1
    return dic

#------------------------------------------------------------------------------
def getBBIntersectingDomains(A, B, tol=1.e-6):
    """Returns the list zones of B that intersect zones of A.
    Usage: getBBIntersectingDomains(A, B, tol)"""
    try: import Generator.PyTree as G
    except: raise ImportError("getBBIntersectingDomains: requires Generator module.")

    BB_A = G.BB(A)
    BB_B = G.BB(B)

    doms = []
    zones2 = Internal.getZones(BB_B)
    for z2 in zones2:
        z2name = Internal.getName(z2)
        zones1 = Internal.getZones(BB_A)
        for z1 in zones1:
            z1name = Internal.getName(z1)
            if z1name != z2name and G.bboxIntersection(z1, z2, tol=tol, isBB=True) == 1:
                doms.append(z1)
    return doms

#------------------------------------------------------------------------------
def getIntersectingDomains(t, t2=None, method='AABB', taabb=None, tobb=None,
                           taabb2=None, tobb2=None, tol=1e-10):
    """Returns the dictionary of intersecting zones.
    Usage: getIntersectingDomains(t, t2, method, taabb, tobb, taabb2, tobb2)"""
    try: import Generator.PyTree as G
    except: raise ImportError("getIntersectingDomains: requires Generator module.")

    #TotInter = 0
    zones = Internal.getZones(t)

    # Gets all zone names and will iterate over them
    zonesNames = []
    for z in zones: zonesNames.append(z[0])
    N = len(zonesNames)
    if not t2:
        zonesNames2 = zonesNames
        zones2 = zones
        #N2 = N
    else:
        zones2 = Internal.getZones(t2)
        # Gets all zones names and will iterate over them
        zonesNames2 = []
        for z in zones2: zonesNames2.append(z[0])
        #N2 = len(zonesNames2)
    # Preallocates Dictionnary
    IntDict = dict.fromkeys(zonesNames)

    if method == 'AABB':
        if not taabb: taabb = G.BB(t)
        if not t2:
            # fast: t, taabb=t avec t bbox
            IntDict = getIntersectingDomainsAABB(taabb, tol=tol)
            #TotInter = 0
        else:
            if not taabb2: taabb2 = G.BB(t2)
            zones1 = Internal.getZones(taabb)
            zones2 = Internal.getZones(taabb2)
            for z1 in zones1:
                IntDict[z1[0]] = []
                for z2 in zones2:
                    if z1[0] != z2[0]:
                        if G.bboxIntersection(z1, z2, tol=tol, isBB=True, method='AABB') == 1:
                            IntDict[z1[0]].append(z2[0]) # saves the intersected zones names
                            #TotInter += 1

    elif method == 'OBB':
        if not tobb: tobb = G.BB(t, method='OBB')
        if not t2: tobb2 = tobb
        elif not tobb2: tobb2 = G.BB(t2, method='OBB')
        for z1 in zonesNames:
            obb1 = Internal.getZones(tobb)
            obb1 = Internal.getNodeFromName(obb1, z1)
            IntDict[z1] = []
            for z2 in zonesNames2:
                if z1 != z2:
                    obb2 = Internal.getZones(tobb2)
                    obb2 = Internal.getNodeFromName(obb2,z2)
                    if G.bboxIntersection(obb1, obb2, tol=tol, isBB=True, method='OBB') == 1:
                        IntDict[z1].append(z2) # saves the intersected zones names
                        #TotInter += 1

    elif method == 'hybrid':
        if not taabb: taabb = G.BB(t)
        if not tobb: tobb = G.BB(t, method='OBB')
        if not t2:
            taabb2 = taabb
            tobb2 = tobb
        if not not t2 and not taabb2: taabb2 = G.BB(t2)
        if not not t2 and not tobb2: tobb2 = G.BB(t2, method='OBB')
        for z1 in zonesNames:
            aabb1 = Internal.getZones(taabb)
            aabb1 = Internal.getNodeFromName(aabb1, z1)
            obb1 = Internal.getZones(tobb)
            obb1 = Internal.getNodeFromName(obb1, z1)
            IntDict[z1] = []
            for z2 in zonesNames2:
                if z1 != z2:
                    obb2 = Internal.getZones(tobb2)
                    obb2 = Internal.getNodeFromName(obb2,z2)
                    aabb2 = Internal.getZones(taabb2)
                    aabb2 = Internal.getNodeFromName(aabb2,z2)
                    if G.bboxIntersection(obb1, obb2, tol=tol, isBB=True,method='OBB') == 0:
                        continue
                    elif G.bboxIntersection(aabb1, aabb2, tol=tol, isBB=True,method='AABB') == 0:
                        continue
                    elif G.bboxIntersection(aabb1, obb2, tol=tol, isBB=True,method='AABBOBB') == 0:
                        continue
                    elif G.bboxIntersection(aabb2, obb1, tol=tol, isBB=True,method='AABBOBB') == 0:
                        continue
                    else:
                        IntDict[z1].append(z2) # saves the intersected zones names
                        #TotInter += 1
    else:
        print('Warning: getIntersectingDomains: method %s not implemented. Switched to AABB.'%method)
        return getIntersectingDomains(t, method='AABB', tol=tol, taabb=taabb, tobb=tobb)

    #print('Total zone/zone intersections: %d.'%TotInter)
    return IntDict

#------------------------------------------------------------------------------
def getCEBBIntersectingDomains(basis0, bases0, sameBase=0, tol=1e-6):
    """Return the list of interpolation domains defined in bases for any zone in basis.
    Usage: getCEBBIntersectingDomains(basis, bases, sameBase)"""
    try: import Generator.PyTree as G
    except: raise ImportError("getCEBBIntersectingDomains: requires Generator module.")

    basis = Internal.getBases(basis0)
    bases = Internal.getBases(bases0)
    if len(basis) != 1: raise TypeError("getCEBBIntersectingDomains: not a CGNS base.")
    if len(bases) < 1: raise TypeError("getCEBBIntersectingDomains: not a list of CGNS bases.")
    basis = basis[0]
    # precond : recherche des bases intersectantes
    [xmin1,ymin1,zmin1,xmax1,ymax1,zmax1] = G.bbox(basis)
    basename1 = basis[0]
    bases2 = []

    for b in bases:
        basename2 = b[0]
        if (sameBase == 1) or (sameBase == 0 and basename1 != basename2):
            [xmin2,ymin2,zmin2,xmax2,ymax2,zmax2] = G.bbox(b)
            if (xmax1 > xmin2-tol and xmin1 < xmax2+tol and
                ymax1 > ymin2-tol and ymin1 < ymax2+tol and
                    zmax1 > zmin2-tol and zmin1 < zmax2+tol):
                bases2.append(b)

    # recherche des zones intersectantes par base
    doms = []
    zones1 = Internal.getNodesFromType1(basis, 'Zone_t')
    for z1 in zones1:
        name1 = z1[0]
        doms1 = []
        for b2 in bases2:
            basename2 = b2[0]
            zones2 = Internal.getNodesFromType1(b2, 'Zone_t')
            if sameBase == 1 and basename1 == basename2:
                for z2 in zones2:
                    name2 = z2[0]
                    if name1 != name2:
                        if G.CEBBIntersection(z1, z2, tol=tol) == 1:
                            doms1.append(z2[0])
            else:
                for z2 in zones2:
                    if G.CEBBIntersection(z1, z2, tol=tol) == 1:
                        doms1.append(z2[0])
        doms.append(doms1)
    return doms

#==============================================================================
def getCEBBTimeIntersectingDomains(base0, func, bases0, funcs, \
                                   inititer=0, niter=1, dt=1, \
                                   sameBase=0, tol=1e-6):
    """Return the list domains defined in bases, that intersect in time
    with base. Func defines a python motion of base with respect to time.
    Usage: getCEBBTimeIntersectingDomains(base, func, bases, funcs,
                                          inititer,niter, dt, sameBase)"""
    try: import RigidMotion.PyTree as R
    except: raise ImportError("getCEBBTimeIntersectingDomains: requires RigidMotion module.")
    base = Internal.getBases(base0)
    bases = Internal.getBases(bases0)
    if len(base) != 1: raise TypeError("getCEBBIntersectingDomains: not a CGNS base.")
    if len(bases) < 1: raise TypeError("getCEBBIntersectingDomains: not a list of CGNS bases.")
    basem = base[0]
    zones = Internal.getNodesFromType1(basem,'Zone_t'); nzones = len(zones)
    zones = []
    t = 0.
    doms = []
    for i in range(inititer, niter):
        t = i*dt
        if func != []: basep = R.evalPosition(basem, t, func)
        else: basep = Internal.copyRef(basem)
        # mouvement relatif des bases
        nob = 0
        basesp = []
        for b in bases:
            if funcs[nob] != []: bp = R.evalPosition(b, t, funcs[nob])
            else: bp = Internal.copyRef(b)
            basesp.append(bp)
            nob += 1
        domsp = getCEBBIntersectingDomains(basep, basesp, sameBase, tol=tol)
        if domsp != [[]]*nzones:
            # premieres intersections
            if doms == []: doms = domsp
            else:
                nod = 0
                for dom1 in doms:
                    dom2 = domsp[nod]
                    doms1 = []
                    for name2 in dom2:
                        found = 0
                        for name1 in dom1:
                            if name1 == name2:
                                found = 1
                                break

                        if found == 0: doms1.append(name2)
                    dom1 += doms1
                    doms[nod] = dom1
                    nod += 1

    return doms


#=============================================================================
#=============================================================================
# 2. CALCUL DES DONNEES D INTERPOLATIONS
#=============================================================================
#=============================================================================
#=============================================================================
# Calcul des donnees pour les interpolations IBC
# IN: aR: arbre, zone, bases des receveurs
#     cellN=2 definit les points a corriger notes PC
#     vect(distance) definit le vecteur normal a la paroi * la distance du PC a la paroi
# IN: aD: arbre, zone, bases des donneurs
# IN: order: ordre des interpolations (2, 3, 5)
# IN: method: 'lagrangian', 'leastsquares'
# IN: penalty=1: penalise une cellule donneuse en terme de volume si elle est au bord
# IN: nature=0: aucun sommet de la cellule d'interpolation ne doit etre avec un cellN=0
#     nature=1: tous les sommets de la cellule d interpolation doivent etre de cellN=1
# IN: interpDataType: 0 pour recherche cartesienne, 1 pour recherche dans ADT
# IN: hook: hook sur l'arbre de recherche (pas reconstruit dans setInterpIBC), l'ordre doit suivre celui de zonesD
# IN: storage: type de stockage (direct: sur le bloc interpole, inverse: sur le bloc d'interpolation)
# IN: loc='nodes','cells','faces': interpolation appliquee pour les receveurs (localises en noeuds/centres/faces)
# IN: hi: decalage des pts IBC interieurs au corps pour avoir le pt interpole - hi peut etre constant ou variable
# IN: he: decalage des pts IBC exterieurs au corps pour avoir le pt interpole - he peut etre constant ou variable
# OUT: stockage direct: retourne tR avec les donnees IBC
# OUT: stockage indirect: retourne tD avec les donnees IBC
#=============================================================================
def setIBCData(tR, tD, order=2, penalty=0, nature=0,
               method='lagrangian', loc='nodes', storage='direct',
               interpDataType=1, hook=None, sameName=0, he=0., hi=0., tfront=None, dim=3, bcType=0,
               Reynolds=6.e6, yplus=100., Lref=1., frontType=0, tb=None):
    """Compute Immersed Boundary data."""
    aR = Internal.copyRef(tR)
    aD = Internal.copyRef(tD)
    _setIBCData(aR, aD, order=order, penalty=penalty, nature=nature, method=method, loc=loc, storage=storage, interpDataType=interpDataType,
                hook=hook, sameName=sameName, he=he, hi=hi, tfront=tfront, dim=dim, bcType=bcType, Reynolds=Reynolds, yplus=yplus, tb=tb, frontType=frontType)

    if storage == 'direct': return aR
    else: return aD

def _setIBCData(aR, aD, order=2, penalty=0, nature=0,
                method='lagrangian', loc='nodes', storage='direct',
                interpDataType=1, hook=None, sameName=0, he=0., hi=0., tfront=None, dim=3, bcType=0,
                Reynolds=6.e6, yplus=100., Lref=1., frontType=0, tb=None):
    """Compute Immersed Boundary data."""
    try: from . import ToolboxIBM as IBM
    except: raise ImportError("setIBCData: requires ToolboxIBM module.")

    locR = loc
    ReferenceState = Internal.getNodeFromType2(aR,'ReferenceState_t')
    model = Internal.getNodeFromName(aR, 'GoverningEquations')
    if model is not None: model = Internal.getValue(model)
    else:                 model = "Euler"

    # Si pas de cellN receveur, on retourne
    if loc == 'nodes': cellNPresent = C.isNamePresent(aR, 'cellN')
    else: cellNPresent = C.isNamePresent(aR, 'centers:cellN')
    if cellNPresent == -1:
        if storage == 'direct': return None
        else: return None

    locCellND = 'nodes'
    # pour l'enlever ensuite si addCellN le cree
    dictIsCellNPresent={}
    for zd in Internal.getZones(aD):
        dictIsCellNPresent[zd[0]]=C.isNamePresent(zd, 'cellN')
    _addCellN__(aD, loc=locCellND)

    zonesRcv = Internal.getZones(aR); nzonesRcv = len(zonesRcv)
    zonesDnr = Internal.getZones(aD); nzonesDnr = len(zonesDnr)

    #---------------------------------------------------------------------------
    # Extraction des pts IBC a corriger + pts parois + interpoles
    #---------------------------------------------------------------------------
    if hook is not None:
        allHooks = hook[:]
        allHooksL = allHooks[:]
    else: allHooks = None

    zonesDnrL = zonesDnr[:]
    nzonesDnrL = nzonesDnr
    #-------------------------------------------
    # 1. Get the list of IBC pts
    #-------------------------------------------
    res = IBM.getAllIBMPoints(zonesRcv, loc=locR, hi=hi, he=he, tfront=tfront, Reynolds=Reynolds, yplus=yplus, Lref=Lref, frontType=frontType, tb=tb)
    if tfront is None:
        correctedPts=res[0][3]; wallPts=res[1][3]; interpPts=res[2][3]
    else:
        correctedPts=list(res[0].values())[0]; wallPts=list(res[1].values())[0]; interpPts=list(res[2].values())[0]
    #-------------------------------------------
    # 2. Interpolation of IBC interp pts
    #-------------------------------------------
    nozr = -1
    for z in zonesRcv:
        nozr += 1
        if loc == 'nodes': cellNPresent = C.isNamePresent(z,'cellN')
        else: cellNPresent = C.isNamePresent(z, 'centers:cellN')
        if cellNPresent != -1:
            if sameName == 1:
                if hook is not None: allHooks = allHooksL[:]
                zonesDnr = zonesDnrL[:]
                nzonesDnr = nzonesDnrL
                cL = 0; found = 0
                for zo in zonesDnr:
                    if zo[0] == z[0]: found = 1; break
                    cL += 1
                if found == 1:
                    if hook is not None: allHooks.pop(cL)
                    zonesDnr.pop(cL); nzonesDnr = nzonesDnr-1
            if  correctedPts[nozr] != []:
                _setIBCDataForZone__(z, zonesDnr, correctedPts[nozr], wallPts[nozr], interpPts[nozr], loc=locR, order=order, penalty=penalty,nature=nature,method=method,\
                                     storage=storage, interpDataType=interpDataType, hook=hook, dim=dim, ReferenceState=ReferenceState, bcType=bcType,model=model)

    # fin parcours des zones receveuses
    for zd in Internal.getZones(aD):
        if dictIsCellNPresent[zd[0]]==-1: C._rmVars(zd,['cellN'])
    return None

def _setIBCDataForZone__(z, zonesDnr, correctedPts, wallPts, interpPts, loc='nodes', \
                         order=2, penalty=0, nature=0, extrap=1, method='lagrangian', storage='direct',\
                         interpDataType=1, hook=None, dim=3, bcType=-1, ReferenceState=None,model="Euler"):

    prefixIBCD ='IBCD_'
    bcName = None
    if isinstance(bcType,int):
        if bcType != -1:
            prefixIBCD += str(bcType)+'_'

    else:
        bcType0 = bcType
        ibcType = bcType.split("#")
        bcType = int(ibcType[0])
        prefixIBCD += ibcType[0]+'_'

        if len(ibcType)>1:
            bcName=ibcType[1]
            prefixIBCD += bcName+'_'

    arraysD = C.getFields(Internal.__GridCoordinates__, zonesDnr, api=1)
    cellND  = C.getField('cellN', zonesDnr, api=1); arraysD = Converter.addVars([arraysD, cellND])

    if model != "Euler" and bcType == -1: bcType = 3

    nzonesDnr = len(arraysD)
    #-------------------------------------------
    # 3. Interpolation of IBC interp pts
    #-------------------------------------------
    # resInterp = [rcvInd1D,donorInd1D,donorType,coefs,extrap,orphan, EXdirs]
    resInterp = Connector.setInterpData__(interpPts, arraysD, order=order, penalty=penalty, \
                                          nature=nature, extrap=extrap, method=method, interpDataType=interpDataType,\
                                          hook=hook, dim=dim)
    if resInterp is not None:
        # Bilan
        nborphan = 0; nbextrapolated = 0; nbinterpolated = 0
        nbinterpolated0 = interpPts[1].shape[1]
        if len(resInterp[4]) > 0:
            for r in resInterp[4]: nbextrapolated += r.shape[0]
        nborphan = resInterp[5].size
        nbinterpolated = nbinterpolated0-nbextrapolated-nborphan

        if bcType != -1:
            print('IBC zone %s: interpolated=%d; extrapolated=%d; orphan=%d for IBC type: %s'%(z[0],nbinterpolated,nbextrapolated,nborphan, IBCTypes[bcType]))
        else:
            print('IBC zone %s: interpolated=%d; extrapolated=%d; orphan=%d '%(z[0],nbinterpolated,nbextrapolated,nborphan))

        if  nborphan>0: print('Warning: zone %s has %d orphan points !'%(z[0], nborphan))

        #-------------------------------------------------------------------------
        # 4. Update the indices of interpolated pts -> indices of corrected points
        #-------------------------------------------------------------------------
        posindcell = KCore.isNamePresent(interpPts, 'indcell')
        indcells = interpPts[1][posindcell,:]
        #-------------------------------------------------------------------------
        # 5. Sort IBC coordinates wrt donor zones
        #-------------------------------------------------------------------------
        allCorrectedPts=[[]]*nzonesDnr; allWallPts=[[]]*nzonesDnr; allMirrorPts=[[]]*nzonesDnr
        xPC0 = correctedPts[1][0,:]; yPC0 = correctedPts[1][1,:]; zPC0 = correctedPts[1][2,:]
        xPW0 = wallPts[1][0,:]; yPW0 = wallPts[1][1,:]; zPW0 = wallPts[1][2,:]
        xPI0 = interpPts[1][0,:]; yPI0 = interpPts[1][1,:]; zPI0 = interpPts[1][2,:]

        #-------------------------------------------------------------------------
        # 6. Set the original indices of corrected pts in receptor zones
        #-------------------------------------------------------------------------
        # Orphan
        if nborphan>0:
            for noi in range(nborphan):
                noind = resInterp[5][noi]
                resInterp[5][noi] = indcells[noind]

        # Interpolated/Extrapolated
        for noz in range(nzonesDnr):
            ninterploc = resInterp[0][noz].size
            if ninterploc>0: # domaine d'interpolation
                correctedPtsL = Converter.array("CoordinateX,CoordinateY,CoordinateZ",ninterploc,1,1)
                xPC = correctedPtsL[1][0,:]; yPC = correctedPtsL[1][1,:]; zPC = correctedPtsL[1][2,:]
                wallPtsL = Converter.array("CoordinateX,CoordinateY,CoordinateZ",ninterploc,1,1)
                xPW = wallPtsL[1][0,:]; yPW = wallPtsL[1][1,:]; zPW = wallPtsL[1][2,:]
                mirrorPtsL = Converter.array("CoordinateX,CoordinateY,CoordinateZ",ninterploc,1,1)
                xPI = mirrorPtsL[1][0,:]; yPI = mirrorPtsL[1][1,:]; zPI = mirrorPtsL[1][2,:]

                for noi in range(ninterploc):
                    index = resInterp[0][noz][noi]
                    # Indices of receptor pts
                    resInterp[0][noz][noi] = indcells[index]
                    # coordinates of receptor/wall/mirror pts
                    xPC[noi] = xPC0[index]; yPC[noi] = yPC0[index]; zPC[noi] = zPC0[index]
                    xPW[noi] = xPW0[index]; yPW[noi] = yPW0[index]; zPW[noi] = zPW0[index]
                    xPI[noi] = xPI0[index]; yPI[noi] = yPI0[index]; zPI[noi] = zPI0[index]

                allCorrectedPts[noz] = correctedPtsL
                allWallPts[noz] = wallPtsL
                allMirrorPts[noz] = mirrorPtsL

        if len(resInterp[4]) > 0: # Sort extrapolated points wrt donor zones
            for noz in range(nzonesDnr):
                nextraploc = resInterp[4][noz].size
                if nextraploc > 0: # Extrapoles
                    for noi in range(nextraploc):
                        index = resInterp[4][noz][noi]
                        resInterp[4][noz][noi] = indcells[index]

        #----------------------------------
        # 7. Stockage dans l'arbre
        # direct: on stocke dans aR
        # inverse: on stocke dans aD
        #----------------------------------
        PL = numpy.array([],dtype=Internal.E_NpyInt)
        PLD = numpy.array([], dtype=Internal.E_NpyInt)
        EXTRAP = numpy.array([], dtype=Internal.E_NpyInt)
        INTERPTYPE = numpy.array([], dtype=Internal.E_NpyInt)
        CFS = numpy.array([],numpy.float64)
        VOL = numpy.array([],numpy.float64)
        ORPHAN = resInterp[5]
        for noz in range(nzonesDnr):
            # ajout des donnees d'interpolation
            ninterploc = resInterp[0][noz].size
            if ninterploc>0:# est un domaine donneur
                if resInterp[4][noz].size > 0:  EXTRAP = resInterp[4][noz]

                if storage == 'direct':
                    _createInterpRegion__(z, zonesDnr[noz][0],resInterp[0][noz],resInterp[1][noz],resInterp[3][noz],\
                                          resInterp[2][noz], VOL, EXTRAP, ORPHAN, \
                                          tag='Receiver', loc=loc,itype='ibc', prefix=prefixIBCD)
                    # add coordinates of corrected point, wall point, interpolated point
                    _addIBCCoords__(z, zonesDnr[noz][0], allCorrectedPts[noz], allWallPts[noz], allMirrorPts[noz], bcType, \
                                    bcName=bcName, ReferenceState=ReferenceState, prefix=prefixIBCD,model=model)
                else: # inverse
                    _createInterpRegion__(zonesDnr[noz],z[0],resInterp[1][noz],resInterp[0][noz],resInterp[3][noz],\
                                          resInterp[2][noz], VOL, EXTRAP, ORPHAN, \
                                          tag='Donor', loc=loc,itype='ibc', prefix=prefixIBCD)
                    # add coordinates of corrected point, wall point, interpolated point
                    _addIBCCoords__(zonesDnr[noz], z[0], allCorrectedPts[noz], allWallPts[noz], allMirrorPts[noz], bcType, \
                                    bcName=bcName, ReferenceState=ReferenceState, prefix=prefixIBCD,model=model)

            elif nborphan==nbinterpolated0: # Only orphan pts: orphan pt list is stored in any candidate donor zone
                if storage == 'direct':
                    _createInterpRegion__(z, zonesDnr[noz][0],PL, PLD, CFS, INTERPTYPE,VOL,EXTRAP,ORPHAN,
                                          tag='Receiver', loc=loc,itype='ibc', prefix=prefixIBCD)
                    # add coordinates of corrected point, wall point, interpolated point
                    _addIBCCoords__(z, zonesDnr[noz][0], correctedPts, wallPts, interpPts, bcType, \
                                    bcName=bcName, ReferenceState=ReferenceState, prefix=prefixIBCD,model=model)
                else: # inverse
                    _createInterpRegion__(zonesDnr[noz],z[0],PL, PLD, CFS, INTERPTYPE,VOL,EXTRAP,ORPHAN,
                                          tag='Donor', loc=loc,itype='ibc', prefix=prefixIBCD)
                    # add coordinates of corrected point, wall point, interpolated point
                    _addIBCCoords__(zonesDnr[noz], z[0], correctedPts, wallPts, interpPts, bcType, \
                                    bcName=bcName, ReferenceState=ReferenceState, prefix=prefixIBCD,model=model)

    return None

# 2nd PI
def _setIBCDataForZone2__(z, zonesDnr, correctedPts, wallPts, interpPts, interpPts2=None, loc='nodes', \
                          order=2, penalty=0, nature=0, extrap=1, method='lagrangian', storage='direct',\
                          interpDataType=1, hook=None, dim=3, bcType=-1, ReferenceState=None):

    prefixIBCD ='IBCD_'
    prefixIBC2D ='2_IBCD_'
    bcName = None
    if isinstance(bcType,int):
        if bcType != -1:
            prefixIBCD += str(bcType)+'_'
            prefixIBC2D += str(bcType)+'_'

    else:
        bcType0 = bcType
        ibcType = bcType.split("#")
        bcType = int(ibcType[0])
        prefixIBCD += ibcType[0]+'_'
        prefixIBC2D += ibcType[0]+'_'

        if len(ibcType)>1:
            bcName = ibcType[1]
            prefixIBCD += bcName+'_'
            prefixIBC2D += bcName+'_'

    arraysD = C.getFields(Internal.__GridCoordinates__, zonesDnr, api=1)
    cellND  = C.getField('cellN', zonesDnr, api=1); arraysD = Converter.addVars([arraysD,cellND])
    model = "Euler"
    a = Internal.getNodeFromName2(zonesDnr[0], 'model')
    if a is not None: model = Internal.getValue(a)
    if model != "Euler" and bcType == -1: bcType = 3

    nzonesDnr = len(arraysD)
    #-------------------------------------------
    # 3. Interpolation of IBC interp pts
    #-------------------------------------------
    # resInterp = [rcvInd1D,donorInd1D,donorType,coefs,extrap,orphan, EXdirs]
    if interpPts is not None:
        resInterp = Connector.setInterpData__(interpPts, arraysD, order=order, penalty=penalty, \
                                              nature=nature, extrap=extrap, method=method, interpDataType=interpDataType,\
                                              hook=hook, dim=dim)
    else:
        resInterp = None

    if interpPts2 is not None:
        resInterp2 = Connector.setInterpData__(interpPts2, arraysD, order=order, penalty=penalty, \
                                               nature=nature, extrap=extrap, method=method, interpDataType=interpDataType,\
                                               hook=hook, dim=dim)
    else:
        resInterp2 = None

    if resInterp is not None:
        # Bilan
        nborphan = 0; nbextrapolated = 0; nbinterpolated = 0
        nbinterpolated0 = interpPts[1].shape[1]
        if len(resInterp[4]) > 0:
            for r in resInterp[4]: nbextrapolated += r.shape[0]
        nborphan = resInterp[5].size
        nbinterpolated = nbinterpolated0-nbextrapolated-nborphan

        if bcType != -1:
            print('IBC zone %s: interpolated=%d; extrapolated=%d; orphan=%d for IBC type: %s'%(z[0],nbinterpolated,nbextrapolated,nborphan, IBCTypes[bcType]))
        else:
            print('IBC zone %s: interpolated=%d; extrapolated=%d; orphan=%d '%(z[0],nbinterpolated,nbextrapolated,nborphan))

        if  nborphan>0: print('Warning: zone %s has %d orphan points !'%(z[0], nborphan))

        #-------------------------------------------------------------------------
        # 4. Update the indices of interpolated pts -> indices of corrected points
        #-------------------------------------------------------------------------
        posindcell = KCore.isNamePresent(interpPts, 'indcell')
        indcells = interpPts[1][posindcell,:]
        #-------------------------------------------------------------------------
        # 5. Sort IBC coordinates wrt donor zones
        #-------------------------------------------------------------------------
        allCorrectedPts=[[]]*nzonesDnr; allWallPts=[[]]*nzonesDnr; allMirrorPts=[[]]*nzonesDnr
        xPC0 = correctedPts[1][0,:]; yPC0 = correctedPts[1][1,:]; zPC0 = correctedPts[1][2,:]
        xPW0 = wallPts[1][0,:]; yPW0 = wallPts[1][1,:]; zPW0 = wallPts[1][2,:]
        xPI0 = interpPts[1][0,:]; yPI0 = interpPts[1][1,:]; zPI0 = interpPts[1][2,:]

        #-------------------------------------------------------------------------
        # 6. Set the original indices of corrected pts in receptor zones
        #-------------------------------------------------------------------------
        # Orphan
        if nborphan>0:
            for noi in range(nborphan):
                noind = resInterp[5][noi]
                resInterp[5][noi] = indcells[noind]

        # Interpolated/Extrapolated
        for noz in range(nzonesDnr):
            ninterploc = resInterp[0][noz].size
            if ninterploc>0: # domaine d'interpolation
                correctedPtsL = Converter.array("CoordinateX,CoordinateY,CoordinateZ",ninterploc,1,1)
                xPC = correctedPtsL[1][0,:]; yPC = correctedPtsL[1][1,:]; zPC = correctedPtsL[1][2,:]
                wallPtsL = Converter.array("CoordinateX,CoordinateY,CoordinateZ",ninterploc,1,1)
                xPW = wallPtsL[1][0,:]; yPW = wallPtsL[1][1,:]; zPW = wallPtsL[1][2,:]
                mirrorPtsL = Converter.array("CoordinateX,CoordinateY,CoordinateZ",ninterploc,1,1)
                xPI = mirrorPtsL[1][0,:]; yPI = mirrorPtsL[1][1,:]; zPI = mirrorPtsL[1][2,:]

                for noi in range(ninterploc):
                    index = resInterp[0][noz][noi]
                    # Indices of receptor pts
                    resInterp[0][noz][noi] = indcells[index]
                    # coordinates of receptor/wall/mirror pts
                    xPC[noi] = xPC0[index]; yPC[noi] = yPC0[index]; zPC[noi] = zPC0[index]
                    xPW[noi] = xPW0[index]; yPW[noi] = yPW0[index]; zPW[noi] = zPW0[index]
                    xPI[noi] = xPI0[index]; yPI[noi] = yPI0[index]; zPI[noi] = zPI0[index]

                allCorrectedPts[noz] = correctedPtsL
                allWallPts[noz] = wallPtsL
                allMirrorPts[noz] = mirrorPtsL

        if len(resInterp[4]) > 0: # Sort extrapolated points wrt donor zones
            for noz in range(nzonesDnr):
                nextraploc = resInterp[4][noz].size
                if nextraploc > 0: # Extrapoles
                    for noi in range(nextraploc):
                        index = resInterp[4][noz][noi]
                        resInterp[4][noz][noi] = indcells[index]

        #----------------------------------
        # 7. Stockage dans l'arbre
        # direct: on stocke dans aR
        # inverse: on stocke dans aD
        #----------------------------------
        PL = numpy.array([], dtype=Internal.E_NpyInt)
        PLD = numpy.array([], dtype=Internal.E_NpyInt)
        EXTRAP = numpy.array([], dtype=Internal.E_NpyInt)
        INTERPTYPE = numpy.array([], dtype=Internal.E_NpyInt)
        CFS = numpy.array([],numpy.float64)
        VOL = numpy.array([],numpy.float64)
        ORPHAN = resInterp[5]
        for noz in range(nzonesDnr):
            # ajout des donnees d'interpolation
            ninterploc = resInterp[0][noz].size
            if ninterploc>0:# est un domaine donneur
                if resInterp[4][noz].size > 0:  EXTRAP = resInterp[4][noz]

                if storage == 'direct':
                    _createInterpRegion__(z, zonesDnr[noz][0],resInterp[0][noz],resInterp[1][noz],resInterp[3][noz],\
                                          resInterp[2][noz], VOL, EXTRAP, ORPHAN, \
                                          tag='Receiver', loc=loc,itype='ibc', prefix=prefixIBCD)
                    # add coordinates of corrected point, wall point, interpolated point
                    _addIBCCoords__(z, zonesDnr[noz][0], allCorrectedPts[noz], allWallPts[noz], allMirrorPts[noz], bcType, \
                                    bcName=bcName, ReferenceState=ReferenceState, prefix=prefixIBCD,model=model)
                else: # inverse
                    _createInterpRegion__(zonesDnr[noz],z[0],resInterp[1][noz],resInterp[0][noz],resInterp[3][noz],\
                                          resInterp[2][noz], VOL, EXTRAP, ORPHAN, \
                                          tag='Donor', loc=loc,itype='ibc', prefix=prefixIBCD)
                    # add coordinates of corrected point, wall point, interpolated point
                    _addIBCCoords__(zonesDnr[noz], z[0], allCorrectedPts[noz], allWallPts[noz], allMirrorPts[noz], bcType, \
                                    bcName=bcName, ReferenceState=ReferenceState, prefix=prefixIBCD,model=model)

            elif nborphan==nbinterpolated0: # Only orphan pts: orphan pt list is stored in any candidate donor zone
                if storage == 'direct':
                    _createInterpRegion__(z, zonesDnr[noz][0],PL, PLD, CFS, INTERPTYPE,VOL,EXTRAP,ORPHAN,
                                          tag='Receiver', loc=loc,itype='ibc', prefix=prefixIBCD)
                    # add coordinates of corrected point, wall point, interpolated point
                    _addIBCCoords__(z, zonesDnr[noz][0], correctedPts, wallPts, interpPts, bcType, \
                                    bcName=bcName, ReferenceState=ReferenceState, prefix=prefixIBCD,model=model)
                else: # inverse
                    _createInterpRegion__(zonesDnr[noz],z[0],PL, PLD, CFS, INTERPTYPE,VOL,EXTRAP,ORPHAN,
                                          tag='Donor', loc=loc,itype='ibc', prefix=prefixIBCD)
                    # add coordinates of corrected point, wall point, interpolated point
                    _addIBCCoords__(zonesDnr[noz], z[0], correctedPts, wallPts, interpPts, bcType, \
                                    bcName=bcName, ReferenceState=ReferenceState, prefix=prefixIBCD,model=model)

    if resInterp2 is not None:
        # Bilan
        nborphan = 0; nbextrapolated = 0; nbinterpolated = 0
        nbinterpolated0 = interpPts2[1].shape[1]
        if len(resInterp2[4]) > 0:
            for r in resInterp2[4]: nbextrapolated += r.shape[0]
        nborphan = resInterp2[5].size
        nbinterpolated = nbinterpolated0-nbextrapolated-nborphan

        if bcType != -1:
            print('2nd Front : IBC zone %s: interpolated=%d; extrapolated=%d; orphan=%d for IBC type: %s'%(z[0],nbinterpolated,nbextrapolated,nborphan, IBCTypes[bcType]))
        else:
            print('2nd Front : IBC zone %s: interpolated=%d; extrapolated=%d; orphan=%d '%(z[0],nbinterpolated,nbextrapolated,nborphan))

        if  nborphan>0: print('Warning: zone %s has %d orphan points !'%(z[0], nborphan))

        #-------------------------------------------------------------------------
        # 4. Update the indices of interpolated pts -> indices of corrected points
        #-------------------------------------------------------------------------
        posindcell = KCore.isNamePresent(interpPts2, 'indcell')
        indcells = interpPts2[1][posindcell,:]
        #-------------------------------------------------------------------------
        # 5. Sort IBC coordinates wrt donor zones
        #-------------------------------------------------------------------------
        allCorrectedPts=[[]]*nzonesDnr; allWallPts=[[]]*nzonesDnr; allMirrorPts=[[]]*nzonesDnr
        xPC0 = correctedPts[1][0,:]; yPC0 = correctedPts[1][1,:]; zPC0 = correctedPts[1][2,:]
        xPW0 = wallPts[1][0,:]; yPW0 = wallPts[1][1,:]; zPW0 = wallPts[1][2,:]
        xPI0 = interpPts2[1][0,:]; yPI0 = interpPts2[1][1,:]; zPI0 = interpPts2[1][2,:]

        #-------------------------------------------------------------------------
        # 6. Set the original indices of corrected pts in receptor zones
        #-------------------------------------------------------------------------
        # Orphan
        if nborphan>0:
            for noi in range(nborphan):
                noind = resInterp2[5][noi]
                resInterp2[5][noi] = indcells[noind]

        # Interpolated/Extrapolated
        for noz in range(nzonesDnr):
            ninterploc = resInterp2[0][noz].size
            if ninterploc>0: # domaine d'interpolation
                correctedPtsL = Converter.array("CoordinateX,CoordinateY,CoordinateZ",ninterploc,1,1)
                xPC = correctedPtsL[1][0,:]; yPC = correctedPtsL[1][1,:]; zPC = correctedPtsL[1][2,:]
                wallPtsL = Converter.array("CoordinateX,CoordinateY,CoordinateZ",ninterploc,1,1)
                xPW = wallPtsL[1][0,:]; yPW = wallPtsL[1][1,:]; zPW = wallPtsL[1][2,:]
                mirrorPtsL = Converter.array("CoordinateX,CoordinateY,CoordinateZ",ninterploc,1,1)
                xPI = mirrorPtsL[1][0,:]; yPI = mirrorPtsL[1][1,:]; zPI = mirrorPtsL[1][2,:]

                for noi in range(ninterploc):
                    index = resInterp2[0][noz][noi]
                    # Indices of receptor pts
                    resInterp2[0][noz][noi] = indcells[index]
                    # coordinates of receptor/wall/mirror pts
                    xPC[noi] = xPC0[index]; yPC[noi] = yPC0[index]; zPC[noi] = zPC0[index]
                    xPW[noi] = xPW0[index]; yPW[noi] = yPW0[index]; zPW[noi] = zPW0[index]
                    xPI[noi] = xPI0[index]; yPI[noi] = yPI0[index]; zPI[noi] = zPI0[index]

                allCorrectedPts[noz] = correctedPtsL
                allWallPts[noz] = wallPtsL
                allMirrorPts[noz] = mirrorPtsL

        if len(resInterp2[4]) > 0: # Sort extrapolated points wrt donor zones
            for noz in range(nzonesDnr):
                nextraploc = resInterp2[4][noz].size
                if nextraploc > 0: # Extrapoles
                    for noi in range(nextraploc):
                        index = resInterp2[4][noz][noi]
                        resInterp2[4][noz][noi] = indcells[index]

        #----------------------------------
        # 7. Stockage dans l'arbre
        # direct: on stocke dans aR
        # inverse: on stocke dans aD
        #----------------------------------
        PL = numpy.array([], dtype=Internal.E_NpyInt)
        PLD = numpy.array([], dtype=Internal.E_NpyInt)
        EXTRAP = numpy.array([], dtype=Internal.E_NpyInt)
        INTERPTYPE = numpy.array([], dtype=Internal.E_NpyInt)
        CFS = numpy.array([],numpy.float64)
        VOL = numpy.array([],numpy.float64)
        ORPHAN = resInterp2[5]
        for noz in range(nzonesDnr):
            # ajout des donnees d'interpolation
            ninterploc = resInterp2[0][noz].size
            if ninterploc>0:# est un domaine donneur
                if resInterp2[4][noz].size > 0:  EXTRAP = resInterp2[4][noz]

                if storage == 'direct':
                    _createInterpRegion__(z, zonesDnr[noz][0],resInterp2[0][noz],resInterp2[1][noz],resInterp2[3][noz],\
                                          resInterp2[2][noz], VOL, EXTRAP, ORPHAN, \
                                          tag='Receiver', loc=loc,itype='ibc2', prefix=prefixIBC2D)
                    # add coordinates of corrected point, wall point, interpolated point
                    _addIBCCoords2__(z, zonesDnr[noz][0], allCorrectedPts[noz], allWallPts[noz], allMirrorPts[noz], bcType, \
                                     bcName=bcName, ReferenceState=ReferenceState, prefix=prefixIBC2D)
                else: # inverse
                    _createInterpRegion__(zonesDnr[noz],z[0],resInterp2[1][noz],resInterp2[0][noz],resInterp2[3][noz],\
                                          resInterp2[2][noz], VOL, EXTRAP, ORPHAN, \
                                          tag='Donor', loc=loc,itype='ibc2', prefix=prefixIBC2D)
                    # add coordinates of corrected point, wall point, interpolated point
                    _addIBCCoords2__(zonesDnr[noz], z[0], allCorrectedPts[noz], allWallPts[noz], allMirrorPts[noz], bcType, \
                                     bcName=bcName, ReferenceState=ReferenceState, prefix=prefixIBC2D)
            elif nborphan==nbinterpolated0: # Only orphan pts: orphan pt list is stored in any candidate donor zone
                if storage == 'direct':
                    _createInterpRegion__(z, zonesDnr[noz][0],PL, PLD, CFS, INTERPTYPE,VOL,EXTRAP,ORPHAN,
                                          tag='Receiver', loc=loc,itype='ibc2', prefix=prefixIBC2D)
                    # add coordinates of corrected point, wall point, interpolated point
                    _addIBCCoords2__(z, zonesDnr[noz][0], correctedPts, wallPts, interpPts2, bcType, \
                                     bcName=bcName, ReferenceState=ReferenceState, prefix=prefixIBC2D)
                else: # inverse
                    _createInterpRegion__(zonesDnr[noz],z[0],PL, PLD, CFS, INTERPTYPE,VOL,EXTRAP,ORPHAN,
                                          tag='Donor', loc=loc,itype='ibc2', prefix=prefixIBC2D)
                    # add coordinates of corrected point, wall point, interpolated point
                    _addIBCCoords2__(zonesDnr[noz], z[0], correctedPts, wallPts, interpPts2, bcType, \
                                     bcName=bcName, ReferenceState=ReferenceState, prefix=prefixIBC2D)

    return None

def _addIBCCoords__(z, zname, correctedPts, wallPts, interpolatedPts, bcType, bcName=None, ReferenceState=None, prefix='IBCD_',model='NSLaminar'):
    nameSubRegion = prefix+zname
    zsr = Internal.getNodeFromName1(z, nameSubRegion)
    coordsPC = Converter.extractVars(correctedPts,['CoordinateX','CoordinateY','CoordinateZ'])
    coordsPW = Converter.extractVars(wallPts, ['CoordinateX','CoordinateY','CoordinateZ'])
    coordsPI = Converter.extractVars(interpolatedPts, ['CoordinateX','CoordinateY','CoordinateZ'])
    zsr[2].append(['CoordinateX_PC',coordsPC[1][0,:], [], 'DataArray_t'])
    zsr[2].append(['CoordinateY_PC',coordsPC[1][1,:], [], 'DataArray_t'])
    zsr[2].append(['CoordinateZ_PC',coordsPC[1][2,:], [], 'DataArray_t'])
    zsr[2].append(['CoordinateX_PI',coordsPI[1][0,:], [], 'DataArray_t'])
    zsr[2].append(['CoordinateY_PI',coordsPI[1][1,:], [], 'DataArray_t'])
    zsr[2].append(['CoordinateZ_PI',coordsPI[1][2,:], [], 'DataArray_t'])
    zsr[2].append(['CoordinateX_PW',coordsPW[1][0,:], [], 'DataArray_t'])
    zsr[2].append(['CoordinateY_PW',coordsPW[1][1,:], [], 'DataArray_t'])
    zsr[2].append(['CoordinateZ_PW',coordsPW[1][2,:], [], 'DataArray_t'])

    # Creation des numpy d'extraction
    nIBC = coordsPC[1].shape[1]
    if ReferenceState is None:
        densNP = numpy.zeros((nIBC),numpy.float64)
        pressNP = numpy.zeros((nIBC),numpy.float64)
    else:
        Pinf = Internal.getNodeFromName1(ReferenceState,'Pressure')
        if Pinf is not None:
            Pinf = Internal.getValue(Pinf)
            pressNP = Pinf*numpy.ones((nIBC),numpy.float64)
        else:
            pressNP = numpy.zeros((nIBC),numpy.float64)

        roinf = Internal.getNodeFromName1(ReferenceState,'Density')
        if roinf is not None:
            roinf = Internal.getValue(roinf)
            densNP = roinf*numpy.ones((nIBC),numpy.float64)
        else:
            densNP = numpy.zeros((nIBC),numpy.float64)

    zsr[2].append(['Pressure', pressNP, [], 'DataArray_t'])
    zsr[2].append(['Density' , densNP , [], 'DataArray_t'])

    vxNP = numpy.zeros((nIBC),numpy.float64)
    vyNP = numpy.zeros((nIBC),numpy.float64)
    vzNP = numpy.zeros((nIBC),numpy.float64)
    zsr[2].append(['VelocityX' , vxNP , [], 'DataArray_t'])
    zsr[2].append(['VelocityY' , vyNP , [], 'DataArray_t'])
    zsr[2].append(['VelocityZ' , vzNP , [], 'DataArray_t'])

    if bcType in [2, 3, 6, 7, 8, 9, 10, 11, 33]:
        utauNP  = numpy.zeros((nIBC),numpy.float64)
        yplusNP = numpy.zeros((nIBC),numpy.float64)
        zsr[2].append(['utau' , utauNP , [], 'DataArray_t'])
        zsr[2].append(['yplus', yplusNP, [], 'DataArray_t'])

    if bcType == 12:
        Temperature_local = numpy.zeros((nIBC),numpy.float64)
        Internal._createChild(zsr, 'Temperature'    , 'DataArray_t', value=Temperature_local)
        Internal._createChild(zsr, 'TemperatureWall', 'DataArray_t', value=Temperature_local)

    if bcType == 13:
        Temperature_local = numpy.zeros((nIBC),numpy.float64)
        Internal._createChild(zsr, 'Temperature' , 'DataArray_t', value=Temperature_local)
        Internal._createChild(zsr, 'WallHeatFlux', 'DataArray_t', value=Temperature_local)

    if bcType == 140:
        val_local = numpy.zeros((nIBC),numpy.float64)
        Internal._createChild(zsr, 'Density_WM'           , 'DataArray_t', value=val_local)
        Internal._createChild(zsr, 'VelocityX_WM'         , 'DataArray_t', value=val_local)
        Internal._createChild(zsr, 'VelocityY_WM'         , 'DataArray_t', value=val_local)
        Internal._createChild(zsr, 'VelocityZ_WM'         , 'DataArray_t', value=val_local)
        Internal._createChild(zsr, 'Temperature_WM'       , 'DataArray_t', value=val_local)
        Internal._createChild(zsr, 'TurbulentSANuTilde_WM', 'DataArray_t', value=val_local)

    if bcType == 33 and model=='NSLaminar':
        val_local = numpy.zeros((nIBC),numpy.float64)
        Internal._createChild(zsr, 't11_model', 'DataArray_t', value=val_local)
        Internal._createChild(zsr, 't12_model', 'DataArray_t', value=val_local)
        Internal._createChild(zsr, 't22_model', 'DataArray_t', value=val_local)
        Internal._createChild(zsr, 't13_model', 'DataArray_t', value=val_local)
        Internal._createChild(zsr, 't23_model', 'DataArray_t', value=val_local)
        Internal._createChild(zsr, 't33_model', 'DataArray_t', value=val_local)

    if bcType == 5:
        stagnationEnthalpy = numpy.zeros((nIBC),numpy.float64)
        Internal._createChild(zsr, 'StagnationEnthalpy', 'DataArray_t', value=stagnationEnthalpy)
        stagnationPressure = numpy.zeros((nIBC),numpy.float64)
        Internal._createChild(zsr, 'StagnationPressure', 'DataArray_t', value=stagnationPressure)
        dirx = numpy.zeros((nIBC),numpy.float64)
        Internal._createChild(zsr, 'dirx', 'DataArray_t', value=dirx)
        diry = numpy.zeros((nIBC),numpy.float64)
        Internal._createChild(zsr, 'diry', 'DataArray_t', value=diry)
        dirz = numpy.zeros((nIBC),numpy.float64)
        Internal._createChild(zsr, 'dirz', 'DataArray_t', value=dirz)

    if bcType == 100:
        KCurvNP = numpy.zeros((nIBC),numpy.float64)
        zsr[2].append([__KCURV__ , KCurvNP , [], 'DataArray_t'])

    if bcType == 10 or bcType == 11:
        gradxPressureNP  = numpy.zeros((nIBC),numpy.float64)
        zsr[2].append(['gradxPressure' , gradxPressureNP , [], 'DataArray_t'])
        gradyPressureNP  = numpy.zeros((nIBC),numpy.float64)
        zsr[2].append(['gradyPressure' , gradyPressureNP , [], 'DataArray_t'])
        gradzPressureNP  = numpy.zeros((nIBC),numpy.float64)
        zsr[2].append(['gradzPressure' , gradzPressureNP , [], 'DataArray_t'])

    if bcType == 11:
        gradxVelocityXNP  = numpy.zeros((nIBC),numpy.float64)
        zsr[2].append(['gradxVelocityX' , gradxVelocityXNP , [], 'DataArray_t'])
        gradyVelocityXNP  = numpy.zeros((nIBC),numpy.float64)
        zsr[2].append(['gradyVelocityX' , gradyVelocityXNP , [], 'DataArray_t'])
        gradzVelocityXNP  = numpy.zeros((nIBC),numpy.float64)
        zsr[2].append(['gradzVelocityX' , gradzVelocityXNP , [], 'DataArray_t'])

        gradxVelocityYNP  = numpy.zeros((nIBC),numpy.float64)
        zsr[2].append(['gradxVelocityY' , gradxVelocityYNP , [], 'DataArray_t'])
        gradyVelocityYNP  = numpy.zeros((nIBC),numpy.float64)
        zsr[2].append(['gradyVelocityY' , gradyVelocityYNP , [], 'DataArray_t'])
        gradzVelocityYNP  = numpy.zeros((nIBC),numpy.float64)
        zsr[2].append(['gradzVelocityY' , gradzVelocityYNP , [], 'DataArray_t'])

        gradxVelocityZNP  = numpy.zeros((nIBC),numpy.float64)
        zsr[2].append(['gradxVelocityZ' , gradxVelocityZNP , [], 'DataArray_t'])
        gradyVelocityZNP  = numpy.zeros((nIBC),numpy.float64)
        zsr[2].append(['gradyVelocityZ' , gradyVelocityZNP , [], 'DataArray_t'])
        gradzVelocityZNP  = numpy.zeros((nIBC),numpy.float64)
        zsr[2].append(['gradzVelocityZ' , gradzVelocityZNP , [], 'DataArray_t'])

    ##Moving IBM
    timeMotion = Internal.getNodeByName(z, 'TimeMotion')
    if timeMotion:
        zsr[2].append(['CoordinateX_PC#Init',coordsPC[1][0,:], [], 'DataArray_t'])
        zsr[2].append(['CoordinateY_PC#Init',coordsPC[1][1,:], [], 'DataArray_t'])
        zsr[2].append(['CoordinateZ_PC#Init',coordsPC[1][2,:], [], 'DataArray_t'])
        zsr[2].append(['CoordinateX_PI#Init',coordsPI[1][0,:], [], 'DataArray_t'])
        zsr[2].append(['CoordinateY_PI#Init',coordsPI[1][1,:], [], 'DataArray_t'])
        zsr[2].append(['CoordinateZ_PI#Init',coordsPI[1][2,:], [], 'DataArray_t'])
        zsr[2].append(['CoordinateX_PW#Init',coordsPW[1][0,:], [], 'DataArray_t'])
        zsr[2].append(['CoordinateY_PW#Init',coordsPW[1][1,:], [], 'DataArray_t'])
        zsr[2].append(['CoordinateZ_PW#Init',coordsPW[1][2,:], [], 'DataArray_t'])

        motionTypeLocal = numpy.ones((nIBC),numpy.float64)
        zsr[2].append(['MotionType'   ,  Internal.getNodeFromName(timeMotion,'MotionType')[1][0]*motionTypeLocal, [], 'DataArray_t'])

        zsr[2].append(['transl_speedX',  Internal.getNodeFromName(timeMotion,'transl_speed')[1][0]*motionTypeLocal, [], 'DataArray_t'])
        zsr[2].append(['transl_speedY',  Internal.getNodeFromName(timeMotion,'transl_speed')[1][1]*motionTypeLocal, [], 'DataArray_t'])
        zsr[2].append(['transl_speedZ',  Internal.getNodeFromName(timeMotion,'transl_speed')[1][2]*motionTypeLocal, [], 'DataArray_t'])

        zsr[2].append(['axis_pntX'    ,  Internal.getNodeFromName(timeMotion,'axis_pnt')[1][0]*motionTypeLocal, [], 'DataArray_t'])
        zsr[2].append(['axis_pntY'    ,  Internal.getNodeFromName(timeMotion,'axis_pnt')[1][1]*motionTypeLocal, [], 'DataArray_t'])
        zsr[2].append(['axis_pntZ'    ,  Internal.getNodeFromName(timeMotion,'axis_pnt')[1][2]*motionTypeLocal, [], 'DataArray_t'])

        zsr[2].append(['axis_vctX'    ,  Internal.getNodeFromName(timeMotion,'axis_vct')[1][0]*motionTypeLocal, [], 'DataArray_t'])
        zsr[2].append(['axis_vctY'    ,  Internal.getNodeFromName(timeMotion,'axis_vct')[1][1]*motionTypeLocal, [], 'DataArray_t'])
        zsr[2].append(['axis_vctZ'    ,  Internal.getNodeFromName(timeMotion,'axis_vct')[1][2]*motionTypeLocal, [], 'DataArray_t'])

        zsr[2].append(['omega'        ,  Internal.getNodeFromName(timeMotion,'omega')[1][0]*motionTypeLocal, [], 'DataArray_t'])


    if bcName is not None:
        Internal._createUniqueChild(zsr, 'FamilyName', 'FamilyName_t', value=bcName)
    return None

# 2nd PI
def _addIBCCoords2__(z, zname, correctedPts, wallPts, interpolatedPts, bcType, bcName=None, ReferenceState=None, prefix='2_IBCD_'):
    nameSubRegion = prefix+zname
    zsr = Internal.getNodeFromName1(z, nameSubRegion)
    coordsPC = Converter.extractVars(correctedPts,['CoordinateX','CoordinateY','CoordinateZ'])
    coordsPW = Converter.extractVars(wallPts, ['CoordinateX','CoordinateY','CoordinateZ'])
    coordsPI = Converter.extractVars(interpolatedPts, ['CoordinateX','CoordinateY','CoordinateZ'])
    zsr[2].append(['CoordinateX_PC',coordsPC[1][0,:], [], 'DataArray_t'])
    zsr[2].append(['CoordinateY_PC',coordsPC[1][1,:], [], 'DataArray_t'])
    zsr[2].append(['CoordinateZ_PC',coordsPC[1][2,:], [], 'DataArray_t'])
    zsr[2].append(['CoordinateX_PI',coordsPI[1][0,:], [], 'DataArray_t'])
    zsr[2].append(['CoordinateY_PI',coordsPI[1][1,:], [], 'DataArray_t'])
    zsr[2].append(['CoordinateZ_PI',coordsPI[1][2,:], [], 'DataArray_t'])
    zsr[2].append(['CoordinateX_PW',coordsPW[1][0,:], [], 'DataArray_t'])
    zsr[2].append(['CoordinateY_PW',coordsPW[1][1,:], [], 'DataArray_t'])
    zsr[2].append(['CoordinateZ_PW',coordsPW[1][2,:], [], 'DataArray_t'])

    # Creation des numpy d'extraction
    nIBC = coordsPC[1].shape[1]
    if ReferenceState is None:
        pressNP = numpy.zeros((nIBC),numpy.float64)
        densNP  = numpy.zeros((nIBC),numpy.float64)
    else:
        Pinf = Internal.getNodeFromName1(ReferenceState,'Pressure')
        if Pinf is not None:
            Pinf = Internal.getValue(Pinf)
            pressNP = Pinf*numpy.ones((nIBC),numpy.float64)
        else:
            pressNP = numpy.zeros((nIBC),numpy.float64)
        roinf = Internal.getNodeFromName1(ReferenceState,'Density')
        if roinf is not None:
            roinf = Internal.getValue(roinf)
            densNP = roinf*numpy.ones((nIBC),numpy.float64)
        else:
            densNP = numpy.zeros((nIBC),numpy.float64)

    zsr[2].append(['Pressure', pressNP, [], 'DataArray_t'])

    gradxPressureNP  = numpy.zeros((nIBC),numpy.float64)
    zsr[2].append(['gradxPressure' , gradxPressureNP , [], 'DataArray_t'])
    gradyPressureNP  = numpy.zeros((nIBC),numpy.float64)
    zsr[2].append(['gradyPressure' , gradyPressureNP , [], 'DataArray_t'])
    gradzPressureNP  = numpy.zeros((nIBC),numpy.float64)
    zsr[2].append(['gradzPressure' , gradzPressureNP , [], 'DataArray_t'])

    if bcName is not None:
        Internal._createUniqueChild(zsr, 'FamilyName', 'FamilyName_t', value=bcName)

    return None

def setInterpData2(aR, aD, order=2, loc='centers', cartesian=False, extrap=1, nature=1, penalty=1, verbose=2):
    """Compute interpolation data for 2 different trees."""
    tD = Internal.copyRef(aD)
    tR = Internal.copyRef(aR)
    _setInterpData2(tR, tD, order=order, loc=loc, cartesian=cartesian, extrap=extrap, nature=nature, penalty=penalty, verbose=verbose)
    return tD

def _setInterpData2(aR, aD, order=2, loc='centers', cartesian=False, extrap=1, nature=1, penalty=1, verbose=2):
    """Compute interpolation data for 2 different trees."""

    if loc == 'nodes': varcelln = 'cellN'
    else: varcelln = 'centers:cellN'

    # Clean previous IDs if necessary
    Internal._rmNodesFromType(aD, 'ZoneSubRegion_t')
    Internal._rmNodesFromName(aD, 'GridCoordinates#Init')

    if cartesian: interpDataType = 0
    else: interpDataType = 1

    cellNPresent = C.isNamePresent(aR, varcelln)
    if cellNPresent == -1: C._initVars(aR, varcelln, 2.) # interp all

    _setInterpData(aR, aD, order=order, penalty=penalty, nature=nature, extrap=extrap,
                   method='lagrangian', loc=loc, storage='inverse', verbose=verbose,
                   interpDataType=interpDataType, sameName=0, itype="chimera")

    if cellNPresent == -1: C._rmVars(aR, [varcelln])

    return None

#==============================================================================
# Calcul des donnees d'interpolation et stockage dans l'arbre CGNS/Python
# ----------------------------------------------------------------------------
# -------------------
# Si la variable cellN n'existe pas, toute la zone receveuse est interpolee
# Si la variable cellN existe, seuls les pts de cellN 2 sont interpoles (noeuds ou centres)
# Si loc='faces' et cellN existe pour des centres, interpolations aux pts EX effectues
# -------------------
# IN: aR: arbre, zone, bases des receveurs
# IN: aD: arbre, zone, bases des donneurs
# IN: order: ordre des interpolations (2, 3, 5)
# IN: method='lagrange', 'leastsquares','conservative'
# IN: penalty=1: penalise une cellule donneuse en terme de volume si elle est au bord
# IN: nature=0: aucun sommet de la cellule d'interpolation ne doit etre avec un cellN=0
#     nature=1: toutes les sommets de la cellule d'interpolation doivent etre de cellN=1
# IN: extrap=0: pas de calcul des point extrapoles
#     extrap=1: calcul et stockage des eventuels pts extrapoles
# IN: interpDataType: 1 for ADT, 0 if donor are cartesian (optimized)
# IN: hook: hook sur l'adt (pas reconstruit dans setInterpData), l'ordre doit suivre celui de zonesD
# IN: storage: type de stockage (direct: sur le bloc interpole, inverse: sur le bloc d'interpolation)
# IN: loc='nodes','cells','faces': interpolation appliquee pour les receveurs (localises en noeuds/centres/faces)
# IN: topTreeR: top tree des receveurs, sert pour extraire les FamilyBC du double wall
# IN: topTreeD: top tree des donneurs, sert pour extraire les FamilyBC du double wall
# IN: sameName: si 1, n'interpole pas a partir des zones donneuses portant le meme nom que les zones receveuses
# IN: itype='both','chimera','abutting': calcule toutes les donnees de transferts, seules les chimere, seules les match/nearmatch
# OUT: stockage direct: retourne tR avec les donnees d'interpolation
# OUT: stockage indirect: retourne tD avec les donnees d'interpolation
# verbose: 0 (rien), 1 (bilan interpolation), 2 (ecrit les indices de pts orphelins),
# RMQ: method='conservative' -> tout le domaine receveur est pour l'instant considere a interpoler (maquette)
#==============================================================================
def setInterpData(aR, aD, order=2, penalty=1, nature=0, extrap=1,
                  method='lagrangian', loc='nodes', storage='direct',
                  interpDataType=1, hook=None, verbose=2,
                  topTreeRcv=None, topTreeDnr=None, sameName=1, dim=3, itype='both', dictOfModels=None):
    """Compute interpolation data for abutting or chimera intergrid connectivity."""
    tR = Internal.copyRef(aR)
    tD = Internal.copyRef(aD)
    _setInterpData(tR, tD, order=order, penalty=penalty, nature=nature, extrap=extrap,
                   method=method, loc=loc, storage=storage, interpDataType=interpDataType,
                   hook=hook, verbose=verbose,
                   topTreeRcv=topTreeRcv, topTreeDnr=topTreeDnr, sameName=sameName, dim=dim, itype=itype, dictOfModels=dictOfModels)
    if storage == 'direct': return tR
    else: return tD

def _setInterpData(aR, aD, order=2, penalty=1, nature=0, extrap=1,
                   method='lagrangian', loc='nodes', storage='direct',
                   interpDataType=1, hook=None, verbose=2,
                   topTreeRcv=None, topTreeDnr=None, sameName=1, dim=3, itype='both', dictOfModels=None):
    """Compute interpolation data for abutting or chimera intergrid connectivity."""

    # Recherche pour les pts coincidents (base sur les GridConnectivity)
    if itype != 'chimera': # abutting
        if storage == 'direct':
            _setInterpDataForGhostCellsStruct__(aR, aD, storage, loc)
        else:
            _setInterpDataForGhostCellsStruct__(aR, aD, storage, loc)
            _setInterpDataForGhostCellsNGon__(aR, aD, storage, loc)

            # Determination du model pour RANS/LES
            #if itype == 'abutting': # SP : a mettre non ? sinon on le refait 2 fois
            _adaptForRANSLES__(aR, aD, dictOfModels=dictOfModels)

    if itype != 'abutting': # chimera
        _setInterpDataChimera__(aR, aD, order=order, penalty=penalty, nature=nature, extrap=extrap, verbose=verbose,
                                method=method, loc=loc, storage=storage, interpDataType=interpDataType, hook=hook,
                                topTreeRcv=topTreeRcv, topTreeDnr=topTreeDnr, sameName=sameName, dim=dim, itype=itype, dictOfModels=dictOfModels)

    # SP: pour l'instant adaptForRANSLES est appele 2 fois: pour les ghost cells et pour le chimere
    # peut on ne le mettre qu ici ?
    #if storage=='inverse': _adaptForRANSLES__(aR, aD)
    return None

def _setInterpDataChimera__(aR, aD, order=2, penalty=1, nature=0, extrap=1,
                            method='lagrangian', loc='nodes', storage='direct',
                            interpDataType=1, hook=None, verbose=2,
                            topTreeRcv=None, topTreeDnr=None, sameName=1, dim=3, itype='both', dictOfModels=None):
    locR = loc

    # Si pas de cellN receveur, on retourne
    if loc == 'nodes': cellNPresent = C.isNamePresent(aR, 'cellN')
    elif loc == 'centers': cellNPresent = C.isNamePresent(aR, 'centers:cellN')
    else: raise ValueError("setInterpData: invalid loc.")
    if cellNPresent == -1 or itype == 'abutting': return None

    locCellND = 'nodes'

    # pour l'enlever ensuite si addCellN le cree
    dictIsCellNPresent = {}
    for zd in Internal.getZones(aD):
        dictIsCellNPresent[zd[0]] = C.isNamePresent(zd, 'cellN')

    _addCellN__(aD, loc=locCellND)

    if method == 'conservative' and itype != 'abutting':
        if loc != 'centers':
            raise ValueError("setInterpData: conservative type is available only for loc='centers'.")
        else: return _setInterpDataConservative__(aR, aD, storage=storage)

    zonesRcv = Internal.getZones(aR); nzonesRcv = len(zonesRcv)
    zonesDnr = Internal.getZones(aD); nzonesDnr = len(zonesDnr)

    arraysD = C.getFields(Internal.__GridCoordinates__, zonesDnr, api=1)
    cellND = C.getField('cellN', zonesDnr, api=1)
    arraysD = Converter.addVars([arraysD,cellND])
    cellND = []

    #---------------------------------------------------------------------------
    # 1. Extraction des points a interpoler
    #    interpPts : un array avec les coordonnees des pts a interpoler modifiees
    # 2. Calcul des coefs et donneurs
    #---------------------------------------------------------------------------
    if hook is not None:
        allHooks = hook[:]
        allHooksL = allHooks[:]
    else: allHooks = None

    arraysDL = arraysD[:]
    zonesDnrL = zonesDnr[:]
    nzonesDnrL = nzonesDnr
    nozr = -1

    for z in zonesRcv:
        nozr += 1
        if loc == 'nodes': cellNPresent = C.isNamePresent(z, 'cellN')
        else: cellNPresent = C.isNamePresent(z, 'centers:cellN')
        if cellNPresent != -1:
            if sameName == 1:
                arraysD = arraysDL[:]
                if hook is not None: allHooks = allHooksL[:]
                zonesDnr = zonesDnrL[:]
                nzonesDnr = nzonesDnrL
                cL = 0; found = 0
                for zo in zonesDnr:
                    if zo[0] == z[0]: found = 1; break;
                    cL += 1
                if found == 1:
                    arraysD.pop(cL)
                    if hook is not None: allHooks.pop(cL)
                    zonesDnr.pop(cL); nzonesDnr = nzonesDnr-1

            #-------------------------------------------
            # Etape 1: on recupere les pts a interpoler
            #-------------------------------------------
            interpPts = []
            if locR == 'nodes':
                an = C.getFields(Internal.__GridCoordinates__, z, api=1)[0]
                cellN = C.getField('cellN', z, api=1)[0]
                an = Converter.addVars([an,cellN])
                interpPts = Connector.getInterpolatedPoints__(an)

            elif locR == 'centers':
                an = C.getFields(Internal.__GridCoordinates__, z, api=1)[0]
                ac = Converter.node2Center(an)
                cellN = C.getField('centers:cellN', z, api=1)[0]
                ac = Converter.addVars([ac,cellN])
                interpPts = Connector.getInterpolatedPoints__(ac)

            #---------------------------------------------
            # Etape 2 : calcul des donnees d'interpolation
            #---------------------------------------------
            # resInterp = [rcvInd1D,donorInd1D,donorType,coefs,extrap,orphan, EXdirs]
            resInterp = Connector.setInterpData__(interpPts, arraysD, order=order, penalty=penalty, nature=nature, extrap=extrap, method=method, interpDataType=interpDataType,
                                                  hook=allHooks, dim=dim)
            if resInterp is not None:
                # Bilan
                nbinterpolated=interpPts[1].shape[1]
                nbextrapolated=0
                if len(resInterp[4])>0:
                    for r in resInterp[4]: nbextrapolated += r.shape[0]
                nborphan = resInterp[5].size
                nbinterpolated = nbinterpolated-nbextrapolated-nborphan
                if verbose != 0:
                    print('Zone %s: interpolated=%d ; extrapolated=%d ; orphan=%d'%(z[0], nbinterpolated, nbextrapolated, nborphan))
                    if nborphan > 0: print('Warning: zone %s has %d orphan points !'%(z[0], nborphan))
                # on remet a une seule zone, attention si x,y,z sont necessaires ensuite
                # les coordonnees peuvent etre fausses a cause du double walls
                indcells=[]
                if loc == 'faces': vari = 'indcell1'
                else: vari = 'indcell'
                posindcell = KCore.isNamePresent(interpPts, vari)
                if posindcell != -1: indcells = interpPts[1][posindcell,:]

                # on recupere les bons indices de pts interpoles (indcell ds interpPts), de EXdir
                # pour ceux obtenus en sortie de setInterpData
                # Orphelins
                if nborphan > 0:
                    for noi in range(nborphan):
                        noind = resInterp[5][noi]
                        resInterp[5][noi] = indcells[noind]
                    if verbose == 3: # force cellN#Orphan=-1
                        listOrphan = resInterp[5]
                        if Internal.getNodeFromName2(z, 'cellN#Orphan') is None:
                            C._initVars(z, "{%s:cellN#Orphan} = {%s:cellN}"%(loc,loc))
                        cellNOrphan = Converter.array('cellN#Orphan', listOrphan.size, 1, 1)
                        cellNOrphan = Converter.initVars(cellNOrphan, 'cellN#Orphan', -1.)
                        C._setPartialFields(z, [cellNOrphan], [listOrphan], loc=locR)

                # Interpoles/Extrapoles
                for noz in range(nzonesDnr):
                    ninterploc = resInterp[0][noz].size
                    if ninterploc > 0: # domaine d'interpolation
                        # Indices des receveurs
                        for noi in range(ninterploc):
                            index = resInterp[0][noz][noi]
                            resInterp[0][noz][noi] = indcells[index]
                if len(resInterp[4])>0: # pts extrapoles par zone donneuse
                    for noz in range(nzonesDnr):
                        nextraploc = resInterp[4][noz].size
                        if nextraploc > 0: # Extrapoles
                            for noi in range(nextraploc):
                                index = resInterp[4][noz][noi]
                                resInterp[4][noz][noi] = indcells[index]
                            if verbose == 3: # force cellN#Orphan=-2
                                listExtrap = resInterp[4][noz]
                                if Internal.getNodeFromName2(z, 'cellN#Orphan') is None:
                                    C._initVars(z, "{%s:cellN#Orphan} = {%s:cellN}"%(loc,loc))
                                cellNOrphan = Converter.array('cellN#Orphan', nextraploc, 1, 1)
                                cellNOrphan = Converter.initVars(cellNOrphan, 'cellN#Orphan', -2.)
                                C._setPartialFields(z, [cellNOrphan], [listExtrap], loc=locR)

                #----------------------------------
                # Etape 3: Stockage dans l'arbre
                # direct: on stocke dans aR
                # inverse: on stocke dans aD
                #----------------------------------
                for noz in range(nzonesDnr):
                    # ajout des donnees d interpolation
                    ninterploc = resInterp[0][noz].size
                    if ninterploc>0: # domaine d'interpolation
                        if storage == 'direct':
                            if resInterp[4][noz].size == 0: extrapPts = numpy.array([], dtype=Internal.E_NpyInt)
                            else: extrapPts = resInterp[4][noz]
                            if resInterp[6] == []: EXdir = numpy.array([], dtype=Internal.E_NpyInt)
                            else: EXdir = resInterp[6][noz]
                            _createInterpRegion__(z, zonesDnr[noz][0], resInterp[0][noz], resInterp[1][noz], resInterp[3][noz], \
                                                  resInterp[2][noz], numpy.array([],numpy.float64), \
                                                  extrapPts, resInterp[5], tag='Receiver', loc=locR, EXDir=EXdir)

                        else: # inverse
                            if resInterp[4][noz].size == 0: extrapPts = numpy.array([], dtype=Internal.E_NpyInt)
                            else: extrapPts = resInterp[4][noz]
                            if resInterp[6] == []: EXdir = numpy.array([], dtype=Internal.E_NpyInt)
                            else: EXdir = resInterp[6][noz]
                            _createInterpRegion__(zonesDnr[noz], z[0], resInterp[1][noz], resInterp[0][noz], resInterp[3][noz], \
                                                  resInterp[2][noz], numpy.array([],numpy.float64), \
                                                  extrapPts, resInterp[5], tag='Donor', loc=locR, EXDir=EXdir)


    if storage != 'direct':
        _adaptForRANSLES__(aR, aD, dictOfModels=dictOfModels)

    # fin parcours des zones receveuses
    for zd in Internal.getZones(aD):
        if dictIsCellNPresent[zd[0]]==-1: C._rmVars(zd,['cellN'])
    return None

#-------------------------------------------------------------------------
# setInterpDataConservative__
#-------------------------------------------------------------------------
def setInterpDataConservative__(tR, tD, storage='direct'):
    aR = Internal.copyRef(tR)
    aD = Internal.copyRef(tD)
    _setInterpDataConservative__(aR,aD,storage=storage)
    if storage=='direct': return aR
    else: return aD

def _setInterpDataConservative__(aR, aD, storage='direct'):
    try: import Post.PyTree as P; import Generator.PyTree as G
    except: raise ImportError("setInterpDataConservative__: requires Post module.")
    locR = 'centers'
    tRBB = G.BB(aR,method='AABB')
    tDBB = G.BB(aD,method='AABB')
    intersectionDict = getIntersectingDomains(aR, aD, method='AABB',taabb=tRBB, taabb2=tDBB)
    # hook sur les zones donneuses
    allDnrCells={}; allDnrZones = {}; allIndicesDnrOrig={}
    for zd in Internal.getZones(aD):
        hookD = C.createHook(zd,'elementCenters')
        zdnr = P.selectCells(zd,'({centers:cellN}>0.)*({centers:cellN}<2.)>0.')
        zdnr = C.convertArray2NGon(zdnr)

        allDnrCells[zd[0]]=zdnr; allDnrZones[zd[0]]=zd
        indicesDnrOrig = C.identifyElements(hookD,zdnr)
        allIndicesDnrOrig[zd[0]]=indicesDnrOrig
        C.freeHook(hookD)

    for zr in Internal.getZones(aR):
        #1. recuperation des points a interpoler
        interpPts = P.selectCells(zr,'{centers:cellN}>1.')
        dimZR = Internal.getZoneDim(interpPts)
        if dimZR[1] >0:
            interpPts=C.convertArray2NGon(interpPts); interpPts = G.close(interpPts)
            #indices correspondant
            hookR = C.createHook(zr,'elementCenters')
            indicesRcvOrig = C.identifyElements(hookR,interpPts)
            C.freeHook(hookR)

            #2. calcul des donnees d interpolation
            # SP : pour optimiser on pourra prendre le dictionnaire d intersection des pts interpoles uniquement
            arraysD=[];  donorZoneNames=[]; listOfIndicesDnrOrig=[]
            for zdnrname in intersectionDict[zr[0]]:
                zdnr = allDnrCells[zdnrname]
                if zdnr != []:
                    ad = C.getFields(Internal.__GridCoordinates__, zdnr, api=1)[0]
                    arraysD.append(ad)
                    listOfIndicesDnrOrig.append(allIndicesDnrOrig[zdnrname])
                    donorZoneNames.append(zdnrname)
            nzonesDnr = len(arraysD)
            interpPtsA = C.getFields(Internal.__GridCoordinates__, interpPts, api=1)[0]

            resInterp = connector.setInterpDataCons(interpPtsA, arraysD, indicesRcvOrig, listOfIndicesDnrOrig)
            if resInterp is not None:
                # Bilan
                nborphan = 0; nbextrapolated = 0; nbinterpolated = 0
                for indicesR in resInterp[0]: nbinterpolated += indicesR.shape[0]
                nborphan = resInterp[4].size
                print('Zone %s: interpolated=%d ; orphan=%d'%(zr[0], nbinterpolated, nborphan))
                if  nborphan>0: print('Warning: zone %s has %d orphan points !'%(zr[0], nborphan))

                # Orphelins
                if nborphan > 0:
                    for noi in range(nborphan):
                        index = resInterp[4][noi]
                        resInterp[4][noi]=indicesRcvOrig[index]
                    orphanPts = resInterp[4]
                else: orphanPts = numpy.array([], dtype=Internal.E_NpyInt)

                extrapPts = numpy.array([], dtype=Internal.E_NpyInt)
                EXdir = numpy.array([], dtype=Internal.E_NpyInt)

                # 3. stockage dans l'arbre
                # on recupere les bons indices de pts interpoles (indcell ds interpPts)
                for noz in range(nzonesDnr):
                    dnrname=donorZoneNames[noz]
                    # ajout des donnees d interpolation
                    ninterploc = resInterp[0][noz].size
                    if ninterploc>0: # domaine d'interpolation
                        if storage == 'direct':
                            _createInterpRegion__(zr, dnrname, resInterp[0][noz], resInterp[1][noz], resInterp[3][noz], \
                                                  resInterp[2][noz], numpy.array([],numpy.float64), \
                                                  extrapPts, orphanPts, tag='Receiver', loc=locR, EXDir=EXdir)

                        else: # inverse
                            extrapPts = numpy.array([], dtype=Internal.E_NpyInt)
                            EXdir = numpy.array([], dtype=Internal.E_NpyInt)
                            zd = allDnrZones[dnrname]
                            _createInterpRegion__(zd, zr[0], resInterp[1][noz], resInterp[0][noz], resInterp[3][noz], \
                                                  resInterp[2][noz], numpy.array([],numpy.float64), \
                                                  extrapPts, orphanPts, tag='Donor', loc=locR, EXDir=EXdir)


    if storage != 'direct': _adaptForRANSLES__(aR, aD)
    return None

#===============================================================================
# Uses the topological information to set the interpolation data in the tree
# for 1-to-1 grid connectivity. Receptor points are ghost points (centers or nodes)
# IN: tR: pyTree with ghost cells already created
#          BCMatch ranges correspond to the original pytree without ghost cells
#          must contain a Rind_t node for each zone
# IN: tD: pyTree of donors (consistent with setInterpTransfers)
# IN: d : number of ghost cells that have been created
# IN: loc='nodes' or 'centers'
# IN: 'storage'='direct' or 'inverse'
# OUT: t: with interpolation data stored in ZoneSubRegion_t nodes of name 'ID*'
#===============================================================================
def _setInterpDataForGhostCellsNGon__(aR, aD, storage='inverse', loc='centers'):
    if loc != 'centers':
        print("WARNING: Connector.OversetData.__setInterpDataForGhostCellsNGon__: only valid for loc='centers'.")
    if storage!='inverse':
        print("WARNING: Connector.OversetData.__setInterpDataForGhostCellsNGon__: only valid for storage='inverse'.")

    indicesExtrap = numpy.array([], dtype=Internal.E_NpyInt)
    indicesOrphan = numpy.array([], dtype=Internal.E_NpyInt)
    vols =  numpy.array([], numpy.float64)
    EXdir = numpy.array([], dtype=Internal.E_NpyInt)
    prefix = 'ID_'
    for zp in Internal.getZones(aR):
        zname = zp[0]
        zoneDimR = Internal.getZoneDim(zp)
        if zoneDimR[0] == 'Structured':
            continue

        rind = -1
        nodes = Internal.getNodesFromType1(zp, 'Elements_t')
        for n in nodes:
            if n[1][0] == 23:
                p = Internal.getNodeFromName1(n, 'IntExt')[1]
                rind = p[0]
                break

        if rind != -1: # rind indices exist : ghost cell data to be computed
            # Array2
            a1 = C.getFields('GridCoordinates', zp, api=3)[0]
            # PE
            PE = Internal.getNodeFromName2(zp, 'ParentElements')[1]
            #
            for gcn in Internal.getNodesFromType2(zp,'GridConnectivity_t'):
                ctype = Internal.getNodeFromType1(gcn,'GridConnectivityType_t')
                if Internal.getValue(ctype)=='Abutting1to1':
                    FL = Internal.getNodeFromName1(gcn, 'PointList')[1]
                    FLd = Internal.getNodeFromName1(gcn, 'PointListDonor')[1]
                    zdonorname = Internal.getValue(gcn)
                    zdonor = Internal.getNodeFromNameAndType(aR, zdonorname, 'Zone_t')
                    if zdonor is None:
                        raise ValueError("setInterpDataForGhostCellsNGon: donor zone not found in aR.")

                    zdonorp = Internal.getNodeFromNameAndType(aD, zdonorname, 'Zone_t')
                    if zdonorp is None:
                        raise ValueError("setInterpDataForGhostCellsNGon: donor zone not found in aD.")

                    Periodic = Internal.getNodeFromType2(gcn,'Periodic_t')
                    RotationAngle=None; RotationCenter=None
                    if Periodic is not None:
                        RotationAngle = Internal.getNodeFromName1(Periodic,'RotationAngle')
                        RotationCenter = Internal.getNodeFromName1(Periodic,'RotationCenter')
                        if RotationAngle is not None:
                            RotationAngle[1][0]=-RotationAngle[1][0]
                            RotationAngle[1][1]=-RotationAngle[1][1]
                            RotationAngle[1][2]=-RotationAngle[1][2]

                    a2 = C.getFields('GridCoordinates', zdonor, api=3)[0]

                    # PE donor
                    PEd = Internal.getNodeFromName2(zdonor, 'ParentElements')[1]

                    # Rind donor
                    rindd = -1
                    nodes = Internal.getNodesFromType1(zdonor, 'Elements_t')
                    for n in nodes:
                        if n[1][0] == 23:
                            pd = Internal.getNodeFromName1(n, 'IntExt')[1]
                            rindd = pd[0]
                            break

                    ret = connector.setInterpDataForGCNGon(FL, FLd, p, pd, a1, a2, PE, PEd)
                    # Stockage
                    _createInterpRegion__(zdonorp, zp[0], ret[1], ret[0], ret[3], ret[2], vols, indicesExtrap,
                                          indicesOrphan, tag='Donor', loc='centers', EXDir=EXdir,
                                          itype='abutting', prefix=prefix,
                                          RotationAngle=RotationAngle, RotationCenter=RotationCenter)

    return None

def _setInterpDataForGhostCellsStruct__(aR, aD, storage='direct', loc='nodes'):

    # empty numpy arrays for zonesubregion nodes
    indicesExtrap = numpy.array([], dtype=Internal.E_NpyInt)
    indicesOrphan = numpy.array([], dtype=Internal.E_NpyInt)
    vols =  numpy.array([], numpy.float64)
    EXdir = numpy.array([], dtype=Internal.E_NpyInt)

    if loc == 'nodes': locR = 0; locS = 'Vertex'
    else: locR = 1; locS = 'CellCenter'

    for zp in Internal.getZones(aR):
        zname = zp[0]
        zoneDimR = Internal.getZoneDim(zp)
        if zoneDimR[0] == 'Unstructured': continue
        else: # Structured
            dimPb = zoneDimR[4]
            rindnode = Internal.getNodeFromType1(zp, 'Rind_t')
            if rindnode is not None: # rind indices exist : ghost cell data to be computed
                rindnode = rindnode[1]
                rindimin = rindnode[0]; rindimax = rindnode[1]
                rindjmin = 0; rindjmax = 0; rindkmin = 0; rindkmax = 0
                if dimPb > 1:
                    rindjmin = rindnode[2]; rindjmax = rindnode[3]
                    if dimPb == 3: rindkmin = rindnode[4]; rindkmax = rindnode[5]

                rindrcv = [rindimin,rindimax,rindjmin,rindjmax,rindkmin,rindkmax]
                imr = zoneDimR[1]; jmr = zoneDimR[2]; kmr = zoneDimR[3]
                if locR == 1: imr=imr-1; jmr=jmr-1; kmr=kmr-1
                listofjoins = Internal.getNodesFromType2(zp, 'GridConnectivity1to1_t')
                for join in listofjoins:
                    # receptor window ranges
                    prange = Internal.getNodeFromName1(join,'PointRange')[1]
                    # donor window ranges
                    prangedonor = Internal.getNodeFromName1(join,'PointRangeDonor')[1]

                    # trirac
                    if dimPb == 3: trirac = [1,2,3]
                    elif dimPb == 2: trirac = [1,2]
                    else: trirac = [1]
                    transfo = Internal.getNodeFromName1(join, 'Transform')
                    if transfo is not None:
                        trirac[0] = transfo[1][0]
                        if dimPb != 1: trirac[1] = transfo[1][1]
                        if dimPb == 3: trirac[2] = transfo[1][2]
                    Periodic = Internal.getNodeFromType2(join, 'Periodic_t')
                    RotationAngle=None; RotationCenter=None
                    if Periodic is not None:
                        RotationAngle = Internal.getNodeFromName1(Periodic,'RotationAngle')
                        RotationCenter = Internal.getNodeFromName1(Periodic,'RotationCenter')
                        if RotationAngle is not None:
                            RotationAngle[1][0]=-RotationAngle[1][0]
                            RotationAngle[1][1]=-RotationAngle[1][1]
                            RotationAngle[1][2]=-RotationAngle[1][2]
                    # donor zone name
                    zdonorname = Internal.getValue(join)
                    zdonor = Internal.getZones(aR)
                    zdonor = Internal.getNodeFromName(zdonor, zdonorname)
                    if zdonor is None: continue
                    # check if the donor zone is defined in aD
                    zdonorp = Internal.getNodesFromName2(aD, zdonorname)
                    if zdonorp == []:
                        raise ValueError("setInterpDataForGhostCells: donor zone not found in donor pyTree.")
                    zdonorp = Internal.getNodesFromType2(zdonorp, 'Zone_t')
                    if zdonorp == []:
                        raise ValueError("setInterpDataForGhostCells: donor zone not found in donor pyTree.")
                    if len(zdonorp)  > 1:
                        print('Warning: setInterpDataForGhostCells: zone name %s defined several times in donor pyTree.'%zdonorname)
                    zdonorp = zdonorp[0]
                    # donor zone dimensions
                    zoneDimD = Internal.getZoneDim(zdonor)
                    imd = zoneDimD[1]; jmd = zoneDimD[2]; kmd = zoneDimD[3]
                    if locR == 1: imd=imd-1; jmd=jmd-1; kmd=kmd-1

                    # rind for donor
                    rindnode = Internal.getNodeFromType1(zdonor, 'Rind_t')
                    if rindnode is not None: # rind indices exist: ghost cell data to be computed
                        rindnode = rindnode[1]
                        rindimin = rindnode[0]; rindimax = rindnode[1]
                        rindjmin = 0; rindjmax = 0; rindkmin = 0; rindkmax = 0
                        if dimPb > 1:
                            rindjmin = rindnode[2]; rindjmax = rindnode[3]
                            if dimPb == 3: rindkmin = rindnode[4]; rindkmax = rindnode[5]
                        #
                        rinddnr = [rindimin,rindimax,rindjmin,rindjmax,rindkmin,rindkmax]
                        # get directions of receptor and donor borders
                        dirR = GhostCells.getDirBorderStruct__(prange,dimPb)
                        dirD = GhostCells.getDirBorderStruct__(prangedonor,dimPb)
                        # get list of border nodes and direction of border
                        if dimPb != 1:
                            # arrayOfIndicesR : indices globaux des 1ers points pres de la frontiere
                            # definie par le prange dimensionnes dans le maillage ghost cells
                            shift = 0 # les indices sont ceux des pts en frontiere max GC
                            #arrayOfIndicesR = GhostCells.getBorderIndicesStruct__(prange,zoneDimR,dirR,0,locS,dimPb,shift)
                            [arrayOfIndicesR, dim1, dim2] = GhostCells.getBorderIndicesStruct__(prange,zoneDimR,dirR,0,locS,dimPb,shift)

                            # listOfIndicesD : indices globaux des 1ers pts donneurs associes a ceux definis par
                            # arrayOfIndicesR
                            shift = 0
                            if dirD  ==-1: shift += rinddnr[0]
                            elif dirD== 1: shift += rinddnr[1]
                            elif dirD==-2: shift += rinddnr[2]
                            elif dirD== 2: shift += rinddnr[3]
                            elif dirD==-3: shift += rinddnr[4]
                            else: shift += rinddnr[5]
                            if dirR  ==-1: shift += rindimin
                            elif dirR== 1: shift += rindimax
                            elif dirR==-2: shift += rindjmin
                            elif dirR== 2: shift += rindjmax
                            elif dirR==-3: shift += rindkmin
                            else: shift += rindkmax
                            if locR == 1: shift -= 1
                            listOfIndicesD = GhostCells.getJoinDonorIndicesStruct__(prange,prangedonor,zoneDimD,dirD, trirac,0,locS,dimPb, shift)

                    # Increments
                    absdirR=abs(dirR); absdirD=abs(dirD)
                    # list of directions of receptor zone
                    dir_rcv = [0]*dimPb; dir_rcv[absdirR-1] = 1
                    # list of direction of donor zone
                    dir_dnr = [0]*dimPb; dir_dnr[absdirD-1] = 1
                    if dimPb == 3:
                        incrR = GhostCells.increment__(dir_rcv[0],dir_rcv[1],dir_rcv[2], imr, jmr, kmr, 0)
                        incrD = GhostCells.increment__(dir_dnr[0],dir_dnr[1],dir_dnr[2], imd, jmd, kmd, 0)
                    elif dimPb == 2:
                        incrR = GhostCells.increment__(dir_rcv[0],dir_rcv[1], imr, jmr, 0)
                        incrD = GhostCells.increment__(dir_dnr[0],dir_dnr[1], imd, jmd, 0)
                    if dirR > 0: incrR = -incrR
                    if dirD < 0: incrD = -incrD

                    if dirR  ==-1: depth = rinddnr[0]
                    elif dirR== 1: depth = rinddnr[1]
                    elif dirR==-2: depth = rinddnr[2]
                    elif dirR== 2:  depth = rinddnr[3]
                    elif dirR==-3: depth = rinddnr[4]
                    else: depth = rinddnr[5]

                    res = connector.setInterpDataForGC(arrayOfIndicesR, listOfIndicesD, \
                                                       dimPb, locR, depth, incrR, incrD)
                    # Stockage
                    prefix = 'ID_'
                    if RotationAngle is not None:
                        val = Internal.getValue(RotationAngle)
                        if val[0]>0. or val[1]>0. or val[2]>0.: prefix='IDPERP_'
                        else: prefix = 'IDPERM_'
                    if storage == 'direct':
                        _createInterpRegion__(zp, zdonorname,res[0],res[1],res[3],res[2],vols,indicesExtrap,\
                                              indicesOrphan,tag='Receiver',loc=loc,EXDir=EXdir,itype='abutting',\
                                              prefix=prefix, RotationAngle=RotationAngle, RotationCenter=RotationCenter)

                    else:
                        _createInterpRegion__(zdonorp, zname, res[1], res[0], res[3], res[2], vols, indicesExtrap,\
                                              indicesOrphan, tag='Donor',loc=loc,EXDir=EXdir,itype='abutting',\
                                              prefix=prefix,RotationAngle=RotationAngle, RotationCenter=RotationCenter)
    return None

# Adapt aD pour RANS/LES
def _adaptForRANSLES__(tR, tD, dictOfModels=None):
    zrdict = {}
    if dictOfModels == None: #sequential behavior
        dictOfModels = {}
        for tp in [tR, tD]:
            bases = Internal.getBases(tp)
            if bases != []:
                for b in bases:
                    model_b = Internal.getNodeFromName2(b, 'GoverningEquations')
                    if model_b is not None: model_b = Internal.getValue(model_b)
                    else: model_b = 'None'
                    for z in Internal.getZones(b):
                        model = Internal.getNodeFromName2(z, 'GoverningEquations')
                        if model is None: model = model_b
                        else: model = Internal.getValue(model)
                        dictOfModels[z[0]] = model
            else:
                for z in Internal.getZones(tp):
                    model = Internal.getNodeFromName2(z, 'GoverningEquations')
                    if model is None: model = 'None'
                    else: model = Internal.getValue(model)
                    dictOfModels[z[0]] = model

    zonesR = Internal.getNodesFromType2(tR, 'Zone_t')
    for zr in zonesR: zrdict[zr[0]] = zr

    zonesD = Internal.getNodesFromType2(tD, 'Zone_t')
    for zd in zonesD:
        subRegions = Internal.getNodesFromType1(zd, 'ZoneSubRegion_t')
        model_zd = dictOfModels[zd[0]]
        for s in subRegions:
            zrcvname = Internal.getValue(s)
            if zrcvname in zrdict:
                zr = zrdict[zrcvname]
                model_zr = dictOfModels[zr[0]]

                if model_zr != 'None' and model_zd != 'None':

                    if (model_zr=='NSTurbulent' or model_zd=='NSTurbulent') and  model_zr != model_zd:
                        datap = numpy.ones(1, dtype=Internal.E_NpyInt)
                        Internal.createUniqueChild(s, 'RANSLES', 'DataArray_t', datap)
                        print("Info: _adaptForRANSLES__: {}/{} <-> {}/{}".format(zr[0], model_zr, zd[0], model_zd))

                    if (model_zr=='LBMLaminar' or model_zd=='LBMLaminar') and model_zr != model_zd:
                        datap = numpy.ones(1, dtype=Internal.E_NpyInt)
                        Internal.createUniqueChild(s, 'NSLBM', 'DataArray_t', datap)
                        print("Info: _adaptForNSLBM__: {}/{} <-> {}/{}".format(zr[0], model_zr, zd[0], model_zd))

    return None

def _createInterpRegion__(z, zname, pointlist, pointlistdonor, interpCoef, interpType, interpVol, indicesExtrap, indicesOrphan, \
                          EXDir=[], pointlistdonorExtC=None, tag='Receiver', loc='centers', itype='chimera',
                          prefix=None, RotationAngle=None, RotationCenter=None):
    if prefix is None:
        if itype == 'chimera': nameSubRegion = 'ID_'+zname
        elif itype == 'abutting': nameSubRegion = 'ID_'+zname
        elif itype == 'ibc': nameSubRegion='IBCD_'+zname
        elif itype == 'ibc2': nameSubRegion='2_IBCD_'+zname
    else: nameSubRegion=prefix+zname

    if RotationAngle is not None: RotationAngle[2]=[] # remove DimensionalUnits child node

    zsr = Internal.getNodesFromName1(z, nameSubRegion)
    # create new subregion for interpolations
    if zsr == []:
        Internal._createChild(z, nameSubRegion, 'ZoneSubRegion_t', value=zname)
        info = z[2][len(z[2])-1]
        Internal._createChild(info, 'ZoneRole', 'DataArray_t', value=tag)

        if loc == 'faces':
            info[2].append(['FaceList',      pointlist, [], 'IndexArray_t'])
            info[2].append(['FaceListDonor', pointlistdonor, [], 'IndexArray_t'])
            if pointlistdonorExtC is not None:
                if tag == 'Receiver':
                    info[2].append(['FaceListDonorExtC', pointlistdonorExtC , [], 'IndexArray_t'])
                else:
                    info[2].append(['FaceListExtC', pointlistdonorExtC , [], 'IndexArray_t'])
            info[2].append(['FaceDirection', EXDir, [], 'DataArray_t'])
            info[2].append(['FaceInterpolantsDonor', interpCoef , [], 'DataArray_t'])
            if interpVol.size > 0: info[2].append(['FaceInterpolantsVol',interpVol , [], 'DataArray_t'])
            info[2].append(['FaceInterpolantsType',interpType , [], 'DataArray_t'])
            if indicesOrphan.size>0: info[2].append(['OrphanFaceList',numpy.unique(indicesOrphan) , [], 'DataArray_t'])
            if indicesExtrap.size>0: info[2].append(['ExtrapFaceList',indicesExtrap , [], 'DataArray_t'])
        else:
            if loc == 'centers':
                Internal._createChild(info, 'GridLocation', 'GridLocation_t', value='CellCenter')
            info[2].append(['PointList',      pointlist , [], 'IndexArray_t'])
            info[2].append(['PointListDonor', pointlistdonor , [], 'IndexArray_t'])
            if pointlistdonorExtC is not None:
                if tag == 'Receiver':
                    info[2].append(['PointListDonorExtC', pointlistdonorExtC , [], 'IndexArray_t'])
                else:
                    info[2].append(['PointListExtC', pointlistdonorExtC , [], 'IndexArray_t'])
            info[2].append(['InterpolantsDonor',interpCoef, [], 'DataArray_t'])
            if interpVol.size>0: info[2].append(['InterpolantsVol',interpVol , [], 'DataArray_t'])
            info[2].append(['InterpolantsType',interpType , [], 'DataArray_t'])
            if indicesOrphan.size>0: info[2].append(['OrphanPointList',numpy.unique(indicesOrphan) , [], 'DataArray_t'])
            if indicesExtrap.size>0: info[2].append(['ExtrapPointList',indicesExtrap, [], 'DataArray_t'])
        if RotationAngle is not None: info[2].append(RotationAngle)
        if RotationCenter is not None: info[2].append(RotationCenter)
    else:
        if loc == 'faces':
            intd = Internal.getNodesFromName1(zsr[0],'FaceInterpolantsDonor')
            if intd == []:
                l = len(zsr[0][2]); info = zsr[0][2][l-1]
                zsr[0][2].append(['FaceList',      pointlist, [], 'IndexArray_t'])
                zsr[0][2].append(['FaceListDonor', pointlistdonor, [], 'IndexArray_t'])
                if pointlistdonorExtC is not None:
                    if tag == 'Receiver':
                        zsr[0][2].append(['FaceListDonorExtC', pointlistdonorExtC , [], 'IndexArray_t'])
                    else:
                        zsr[0][2].append(['FaceListExtC', pointlistdonorExtC , [], 'IndexArray_t'])
                zsr[0][2].append(['FaceInterpolantsDonor',interpCoef, [], 'DataArray_t'])
                zsr[0][2].append(['FaceDirection',EXDir, [], 'DataArray_t'])
                if interpVol.size>0: zsr[0][2].append(['FaceInterpolantsVol',interpVol, [], 'DataArray_t'])
                zsr[0][2].append(['FaceInterpolantsType',interpType, [], 'DataArray_t'])
                if indicesOrphan.size>0: zsr[0][2].append(['OrphanFaceList',numpy.unique(indicesOrphan),[],'DataArray_t'])
                if indicesExtrap.size>0: zsr[0][2].append(['ExtrapFaceList', indicesExtrap,[],'DataArray_t'])
                if RotationAngle is not None: zsr[0][2].append(RotationAngle)
                if RotationCenter is not None: zsr[0][2].append(RotationCenter)
            else:
                intd[0][1] = numpy.concatenate((intd[0][1],interpCoef))
                intv = Internal.getNodesFromName1(zsr[0],'FaceInterpolantsVol')
                if intv!= [] and interpVol.size>0: intv[0][1]=numpy.concatenate((intv[0][1],interpVol))
                inttype = Internal.getNodesFromName1(zsr[0],'FaceInterpolantsType')
                if inttype != []: inttype[0][1] = numpy.concatenate((inttype[0][1],interpType))
                node = Internal.getNodesFromName1(zsr[0],'FaceList')
                if node!=[]: node[0][1] = numpy.concatenate((node[0][1],pointlist))
                node = Internal.getNodesFromName1(zsr[0],'FaceListDonor')
                if node!=[]: node[0][1] = numpy.concatenate((node[0][1],pointlistdonor))
                if pointlistdonorExtC is not None:
                    if tag == 'Receiver':
                        node = Internal.getNodesFromName1(zsr[0],'FaceListDonorExtC')
                        if node!=[]: node[0][1]=numpy.concatenate((node[0][1],pointlistdonorExtC))
                    else:
                        node = Internal.getNodesFromName1(zsr[0],'FaceListExtC')
                        if node!=[]: node[0][1]=numpy.concatenate((node[0][1],pointlistdonorExtC))
                node = Internal.getNodesFromName1(zsr[0],'FaceDirection')
                if node!=[]: node[0][1] = numpy.concatenate((node[0][1],EXDir))
                node = Internal.getNodesFromName1(zsr[0],'OrphanFaceList')
                if node!=[]: node[0][1]=numpy.unique(indicesOrphan)
                node = Internal.getNodesFromName1(zsr[0],'ExtrapFaceList')
                if node!=[]: node[0][1]=numpy.concatenate((node[0][1],indicesExtrap))
        else:
            intd = Internal.getNodesFromName1(zsr[0],'InterpolantsDonor')
            if intd == []:
                l = len(zsr[0][2]); info = zsr[0][2][l-1]
                if loc == 'centers':
                    Internal._createChild(zsr[0], 'GridLocation', 'GridLocation_t', value='CellCenter')
                zsr[0][2].append(['PointList',      pointlist, [], 'IndexArray_t'])
                zsr[0][2].append(['PointListDonor', pointlistdonor, [], 'IndexArray_t'])
                if pointlistdonorExtC is not None:
                    if tag == 'Receiver':
                        zsr[0][2].append(['PointListDonorExtC', pointlistdonorExtC , [], 'IndexArray_t'])
                    else:
                        zsr[0][2].append(['PointListExtC', pointlistdonorExtC , [], 'IndexArray_t'])
                zsr[0][2].append(['InterpolantsDonor',interpCoef, [], 'DataArray_t'])
                if interpVol.size>0: zsr[0][2].append(['InterpolantsVol',interpVol, [], 'DataArray_t'])
                zsr[0][2].append(['InterpolantsType',interpType, [], 'DataArray_t'])
                if indicesOrphan.size>0: zsr[0][2].append(['OrphanPointList',numpy.unique(indicesOrphan),[], 'DataArray_t'])
                if indicesExtrap.size>0: zsr[0][2].append(['ExtrapPointList',indicesExtrap, [], 'DataArray_t'])
                if RotationAngle is not None: zsr[0][2].append(RotationAngle)
                if RotationCenter is not None: zsr[0][2].append(RotationCenter)
            else:
                intd[0][1] = numpy.concatenate((intd[0][1],interpCoef))
                intv = Internal.getNodesFromName1(zsr[0],'InterpolantsVol')
                if intv!= [] and interpVol.size>0: intv[0][1] = numpy.concatenate((intv[0][1],interpVol))
                inttype = Internal.getNodesFromName1(zsr[0],'InterpolantsType')
                if inttype != []: inttype[0][1] = numpy.concatenate((inttype[0][1],interpType))
                node = Internal.getNodesFromName1(zsr[0],'PointList')
                if node!=[]: node[0][1] = numpy.concatenate((node[0][1],pointlist))
                node = Internal.getNodesFromName1(zsr[0],'PointListDonor')
                if node!=[]: node[0][1] = numpy.concatenate((node[0][1],pointlistdonor))
                node = Internal.getNodesFromName1(zsr[0],'PointListExtC')
                if node!=[]: node[0][1] = numpy.concatenate((node[0][1],pointlistdonorExtC))
                node = Internal.getNodesFromName1(zsr[0],'PointListDonorExtC')
                if node!=[]: node[0][1] = numpy.concatenate((node[0][1],pointlistdonorExtC))
                node = Internal.getNodesFromName1(zsr[0],'OrphanPointList')
                if node!=[]: node[0][1] = numpy.unique(indicesOrphan)
                node = Internal.getNodesFromName1(zsr[0],'ExtrapPointList')
                if node!=[]: node[0][1] = numpy.concatenate((node[0][1],indicesExtrap))

    return None


#=============================================================================
#=============================================================================
# 3. TRANSFERTS
# interpDataType = 0 if donor is Cartesian
# interpDataType = otherwise (an ADT must define the hook)
#=============================================================================
#=============================================================================
def transferFields(z, interpXN, interpYN, interpZN, order=2, penalty=1, nature=0,
                   constraint=30., hook=None, variables=[], interpDataType=1):
    """Transfer given fields once interpolation data has been computed."""
    return connector.transferFields(z, interpXN, interpYN, interpZN,
                                    order, nature, penalty, constraint,
                                    hook, variables, interpDataType,
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)

#===============================================================================
# General transfers: Match + Chimera + IBC
# Interpolation is applied to aR
# Beware: variables must be defined in topTreeD at nodes in order to
# be consistent with the computation of connectivity by setInterpData and
# setIBCData
# loc='nodes','centers' defines the location in aR of transfered values
# IN: variablesID=['var1','var2',...]: variables to be used in Chimera transfers
#                =[] : the whole FlowSolutionNodes variables in topTreeD are transferred
# IN: variablesIBC=['var1','var2',...,'var5']: variables used in IBC transfers
# IN: bcType (IBC only)  see TypesOfIBC dictionary at top of this file
# IN: varType: defines the meaning of the variables IBC
#     varType = 2 : (ro,u,v,w,t)
#     varType = 21: (ro,u,v,w,t(,nutildeSA))
# CAUTION !!!!! for IBC transfers, compact=1 is mandatory : numpys in IBCD must be aligned (Density, Pressure, utau etc)
# IN: storage=-1/0/1: unknown/direct/inverse
# IN: loc = 'nodes' or 'centers': location of receiver zone field
# Pour les IBCs avec loi de paroi, il faut specifier Gamma, Cv, MuS, Cs, Ts
# compactD : 1 for IBC transfers. In that case, fields in IBCD zonesubregions must be aligned with density
#===============================================================================
def setInterpTransfers(aR, topTreeD, variables=[], cellNVariable='cellN',
                       variablesIBC=['Density','VelocityX','VelocityY','VelocityZ','Temperature'],
                       bcType=0, varType=2, storage=-1, compact=0, compactD=0,
                       Gamma=1.4, Cv=1.7857142857142865, MuS=1.e-08,
                       Cs=0.3831337844872463, Ts=1.0):
    """Transfer variables once interpolation data has been computed."""
    # if variablesIBC is not None:
    #     nvarIBC = len(variablesIBC)
    #     if nvarIBC != 5 and nvarIBC != 6:
    #         raise ValueError("setInterpTransfers: length of variablesIBC must be equal to 5.")

    tR = Internal.copyRef(aR)

    _setInterpTransfers(tR, topTreeD, variables=variables, variablesIBC=variablesIBC,
                        bcType=bcType, varType=varType, storage=storage, compact=compact, compactD=compactD,
                        cellNVariable=cellNVariable, Gamma=Gamma, Cv=Cv, MuS=MuS, Cs=Cs, Ts=Ts)
    return tR

#===============================================================================
# General transfers: Chimera + IBC - inplace version
# Interpolation is applied to aR
# Beware: variables must be defined in topTreeD at nodes in order to be
# consistent with the computation of connectivity by setInterpData and
# setIBCData
# loc='nodes','centers' defines the location in aR of transferred values
# IN: variablesI =['var1','var2',...]: variables to be used in Chimera transfers
#                =[]: the whole FlowSolutionNodes variables in topTreeD are transferred
# IN: variablesIBC=['var1','var2',...,'var5']: variables used in IBC transfers
# NB: if variables and/or variablesIBC are None, then no transfer (interp or IBC) is done
# IN: bcType (IBC only)  see TypesOfIBC dictionary at top of this file
# IN: varType: defines the meaning of the variables IBC
#     varType = 2 : (ro,u,v,w,t)
#     varType = 21: (ro,u,v,w,t(,nutildeSA))
# IN: storage=-1/0/1: unknown/direct/inverse
# Pour les IBCs avec loi de paroi, il faut specifier Gamma, Cv, MuS, Cs, Ts
# compact: 1 if fields in t are compacted
# compactD: 1 for IBC transfers. In that case, fields in IBCD zonesubregions must be aligned with density
#===============================================================================
def _setInterpTransfers(aR, topTreeD, variables=[], cellNVariable='',
                        variablesIBC=['Density','VelocityX','VelocityY','VelocityZ','Temperature'],
                        bcType=0, varType=2, storage=-1, compact=0, compactD=0,
                        Gamma=1.4, Cv=1.7857142857142865, MuS=1.e-08,Cs=0.3831337844872463, Ts=1.0):

    # Recup des donnees a partir des zones receveuses
    if storage < 1:
        # Dictionnaire pour optimisation
        znd = {}
        zones = Internal.getZones(topTreeD)
        for z in zones: znd[z[0]] = z

        zonesR = Internal.getZones(aR)
        for zr in zonesR:
            subRegions = Internal.getNodesFromType1(zr, 'ZoneSubRegion_t')
            for s in subRegions:
                sname = s[0][0:2]
                # test pour eviter parcours arbre inutile
                if ((sname == 'ID') and variables is not None) or (sname == 'IB' and variablesIBC is not None):
                    idn = Internal.getNodeFromName1(s, 'InterpolantsDonor')
                    if idn is not None: # la subRegion decrit des interpolations
                        zoneRole = Internal.getNodeFromName2(s, 'ZoneRole')
                        zoneRole = Internal.getValue(zoneRole)
                        if zoneRole == 'Receiver':
                            location = Internal.getNodeFromName1(s, 'GridLocation')
                            if location is not None: location = Internal.getValue(location)
                            Coefs = idn[1]
                            DonorType = Internal.getNodeFromName1(s, 'InterpolantsType')[1]
                            ListRcv   = Internal.getNodeFromName1(s, 'PointList')[1]
                            ListDonor = Internal.getNodeFromName1(s, 'PointListDonor')[1]
                            # Recup des champs du receveur
                            zdnrname = Internal.getValue(s)
                            zd = znd[zdnrname]

                            if location == 'CellCenter': loc = 1
                            else: loc = 0
                            # Transferts
                            if sname == 'ID':
                                RotAngleNode = Internal.getNodeFromName1(s, 'RotationAngle')
                                RotAngleX = 0.; RotAngleY = 0.; RotAngleZ = 0.
                                if RotAngleNode is not None:
                                    RotAngle = Internal.getValue(RotAngleNode)
                                    RotAngleX = RotAngle[0]
                                    RotAngleY = RotAngle[1]
                                    RotAngleZ = RotAngle[2]

                                connector._setInterpTransfers(zr, zd, variables, ListRcv, ListDonor, DonorType, Coefs, loc, varType, compact,
                                                              cellNVariable,
                                                              Internal.__GridCoordinates__,
                                                              Internal.__FlowSolutionNodes__,
                                                              Internal.__FlowSolutionCenters__,
                                                              RotAngleX, RotAngleY, RotAngleZ)

                            elif sname == 'IB' and compactD:
                                xPC = Internal.getNodeFromName1(s,'CoordinateX_PC')[1]
                                yPC = Internal.getNodeFromName1(s,'CoordinateY_PC')[1]
                                zPC = Internal.getNodeFromName1(s,'CoordinateZ_PC')[1]
                                xPW = Internal.getNodeFromName1(s,'CoordinateX_PW')[1]
                                yPW = Internal.getNodeFromName1(s,'CoordinateY_PW')[1]
                                zPW = Internal.getNodeFromName1(s,'CoordinateZ_PW')[1]
                                xPI = Internal.getNodeFromName1(s,'CoordinateX_PI')[1]
                                yPI = Internal.getNodeFromName1(s,'CoordinateY_PI')[1]
                                zPI = Internal.getNodeFromName1(s,'CoordinateZ_PI')[1]
                                density = Internal.getNodeFromName1(s,'Density')[1]
                                # Transferts
                                connector._setIBCTransfers(zr, zd, variablesIBC, ListRcv, ListDonor, DonorType, Coefs,
                                                           xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI,
                                                           density,
                                                           bcType, loc, varType, compactD, Gamma, Cv, MuS, Cs, Ts,
                                                           Internal.__GridCoordinates__,
                                                           Internal.__FlowSolutionNodes__,
                                                           Internal.__FlowSolutionCenters__)

    # Recup des donnees a partir des zones donneuses
    if storage != 0:
        # Dictionnaire pour optimisation
        znr = {}
        zones = Internal.getZones(aR)
        for z in zones: znr[z[0]] = z
        zonesD = Internal.getZones(topTreeD)
        for zd in zonesD:
            subRegions = Internal.getNodesFromType1(zd, 'ZoneSubRegion_t')
            for s in subRegions:
                sname = s[0][0:2]
                # test pour eviter parcours arbre inutile
                if (sname=='ID' and variables is not None) or (sname == 'IB' and variablesIBC is not None):
                    idn = Internal.getNodeFromName1(s, 'InterpolantsDonor')
                    if idn is not None:
                        zoneRole = Internal.getNodeFromName2(s, 'ZoneRole')
                        zoneRole = Internal.getValue(zoneRole)
                        if zoneRole == 'Donor':
                            location = Internal.getNodeFromName1(s, 'GridLocation') #localisation des donnees des rcvr
                            if location is not None: location = Internal.getValue(location)
                            Coefs = idn[1]
                            DonorType = Internal.getNodeFromName1(s,'InterpolantsType')[1]
                            ListDonor = Internal.getNodeFromName1(s,'PointList')[1]
                            ListRcv   = Internal.getNodeFromName1(s,'PointListDonor')[1]
                            # Recup des champs du receveur
                            zrcvname = Internal.getValue(s)
                            ##zr = Internal.getNodesFromName2(aR, zrcvname)
                            zr = znr.get(zrcvname, None)
                            if zr is not None:
                                if location == 'CellCenter': loc = 1
                                else: loc = 0
                                # Transferts
                                if sname == 'ID':
                                    RotAngleNode = Internal.getNodeFromName1(s,'RotationAngle')
                                    RotAngleX = 0.; RotAngleY = 0.; RotAngleZ = 0.
                                    if RotAngleNode is not None:
                                        RotAngleNode=RotAngleNode[1]
                                        RotAngleX = RotAngleNode[0]
                                        RotAngleY = RotAngleNode[1]
                                        RotAngleZ = RotAngleNode[2]
                                    connector._setInterpTransfers(zr, zd, variables, ListRcv, ListDonor, DonorType, Coefs,loc, varType, compact, cellNVariable,
                                                                  Internal.__GridCoordinates__,
                                                                  Internal.__FlowSolutionNodes__,
                                                                  Internal.__FlowSolutionCenters__,
                                                                  RotAngleX, RotAngleY, RotAngleZ)
                                elif sname == 'IB' and compactD:
                                    xPC = Internal.getNodeFromName1(s,'CoordinateX_PC')[1]
                                    yPC = Internal.getNodeFromName1(s,'CoordinateY_PC')[1]
                                    zPC = Internal.getNodeFromName1(s,'CoordinateZ_PC')[1]
                                    xPW = Internal.getNodeFromName1(s,'CoordinateX_PW')[1]
                                    yPW = Internal.getNodeFromName1(s,'CoordinateY_PW')[1]
                                    zPW = Internal.getNodeFromName1(s,'CoordinateZ_PW')[1]
                                    xPI = Internal.getNodeFromName1(s,'CoordinateX_PI')[1]
                                    yPI = Internal.getNodeFromName1(s,'CoordinateY_PI')[1]
                                    zPI = Internal.getNodeFromName1(s,'CoordinateZ_PI')[1]
                                    Density = Internal.getNodeFromName1(s,'Density')[1]

                                    connector._setIBCTransfers(zr, zd, variablesIBC, ListRcv, ListDonor, DonorType, Coefs,
                                                               xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI,
                                                               Density,
                                                               bcType, loc, varType, compact, Gamma, Cv, MuS, Cs, Ts,
                                                               Internal.__GridCoordinates__,
                                                               Internal.__FlowSolutionNodes__,
                                                               Internal.__FlowSolutionCenters__)
    return None

#===============================================================================
# General transfers: Chimera + IBC - inplace version optimiser par arbre tc compacte par zone donneuse
# Interpolation is applied to aR
# Beware: variables must be defined in topTreeD at nodes in order to be
# consistent with the computation of connectivity by setInterpData and
# setIBCData
# loc='nodes','centers' defines the location in aR of transferred values
# IN: variablesI =['var1','var2',...]: variables to be used in Chimera transfers
#                = None: the whole FlowSolutionNodes variables in topTreeD are transferred
# IN: variablesIBC=['var1','var2',...,'var5']: variables used in IBC transfers
# IN: bcType (IBC only) 0: glissement
#                       1: adherence
#                       2: loi de paroi log
#                       3: loi de paroi Musker
# IN: varType: defines the meaning of the variables IBC
#     varType = 1 : (ro,rou,rov,row,roE)
#     varType = 11: (ro,rou,rov,row,roE)+ronultideSA
#     varType = 2 : (ro,u,v,w,t)
#     varType = 21: (ro,u,v,w,t)+nultideSA
#     varType = 3 : (ro,u,v,w,p)
#     varType = 31: (ro,u,v,w,p)+nultideSA
# IN: storage=-1/0/1: unknown/direct/inverse
# Pour les IBCs avec loi de paroi, il faut specifier Gamma, Cv, MuS, Cs, Ts
#===============================================================================
def __setInterpTransfers(aR, topTreeD,
                         variables=[],
                         variablesIBC=['Density','MomentumX','MomentumY','MomentumZ','EnergyStagnationDensity'],
                         bcType=0, varType=1, storage=-1, compact=0,
                         Gamma=1.4, Cv=1.7857142857142865, MuS=1.e-08,
                         Cs=0.3831337844872463, Ts=1.0):

    # Recup des donnees a partir des zones receveuses
    if storage != 1 or compact == 0:
        raise ValueError("__setInterpTransfers: stockage chez le receveur non code. Mode compact obligatoire: compact=1 a imposer dans Fast.warmup.")

    # test pour savoir si appel transfert ou ibc
    flagibc = 1
    if variables is not None: flagibc = 0

    # Recup des donnees a partir des zones donneuses
    zones  = Internal.getZones(aR)
    zonesD = Internal.getZones(topTreeD)
    for zd in zonesD:
        param_int = Internal.getNodeFromName1(zd, 'Parameter_int')[1]
        param_real= Internal.getNodeFromName1(zd, 'Parameter_real')[1]
        if param_int[0] != 0 and flagibc == 0:
            connector.__setInterpTransfers(zones, zd, variables, param_int, param_real, varType, compact, flagibc, bcType, Gamma, Cv, MuS, Cs, Ts )
        elif param_int[1] != 0 and flagibc == 1:
            connector.__setInterpTransfers(zones, zd, variablesIBC, param_int, param_real, varType, compact, flagibc, bcType, Gamma, Cv, MuS, Cs, Ts )
    return None

#===============================================================================
# General transfers: Chimera + IBC parallele - getFromZone
# Beware: variables must be defined in topTreeD at nodes in order to be consistent with
# the computation of connectivity by setInterpData and setIBCData
# loc='nodes','centers' defines the location in aR of transferred values
# IN: variables=['var1','var2',...]: variables to be used in Chimera transfers
#              =[]: the whole FlowSolutionNodes variables in topTreeD are transferred
# IN: variablesIBC=['var1','var2',...,'var5']:
# IN: bcType (IBC only)  see TypesOfIBC dictionary at top of this file
# IN: varType: defines the meaning of the variables IBC
#     varType = 2 : (ro,u,v,w,t)
#     varType = 21: (ro,u,v,w,t(,nutildeSA))
# CAUTION !!!!! for IBC transfers, compact=1 is mandatory : numpys in IBCD must be aligned (Density, Pressure, utau etc)
# Adim: KCore.adim1 for Minf=0.1 (IBC only)
#===============================================================================
def setInterpTransfersD(topTreeD, variables=[], cellNVariable='',
                        variablesIBC=['Density','VelocityX','VelocityY','VelocityZ','Temperature'],
                        bcType=0, varType=2, compact=0, compactD=0,
                        Gamma=1.4, Cv=1.7857142857142865, MuS=1.e-08,
                        Cs=0.3831337844872463, Ts=1.0, extract=0, alpha=1.):
    tD = Internal.copyRef(topTreeD)
    return _setInterpTransfersD(tD, variables=variables, cellNVariable=cellNVariable,
                                variablesIBC=variablesIBC,
                                bcType=bcType, varType=varType, compact=compact, compactD=compactD,
                                Gamma=Gamma, Cv=Cv, MuS=MuS,
                                Cs=Cs, Ts=Ts, extract=extract)

def _setInterpTransfersD(topTreeD, variables=[], cellNVariable='',
                         variablesIBC=['Density','VelocityX','VelocityY','VelocityZ','Temperature'],
                         bcType=0, varType=2, compact=0, compactD=0,
                         Gamma=1.4, Cv=1.7857142857142865, MuS=1.e-08,
                         Cs=0.3831337844872463, Ts=1.0, extract=0, alpha=1.):

    # Recup des donnees a partir des zones donneuses
    zonesD = Internal.getZones(topTreeD)
    infos = []
    for zd in zonesD:
        subRegions = Internal.getNodesFromType1(zd, 'ZoneSubRegion_t')
        for s in subRegions:
            sname = s[0][0:2]
            #test pour eviter parcours arbre inutile
            if (sname=='ID' and variables is not None) or (sname == 'IB' and variablesIBC is not None):
                dname = Internal.getValue(s)
                idn = Internal.getNodeFromName1(s, 'InterpolantsDonor')
                if idn is not None: # la subRegion decrit des interpolations
                    zoneRole = Internal.getNodeFromName2(s, 'ZoneRole')
                    zoneRole = Internal.getValue(zoneRole)
                    if zoneRole == 'Donor':
                        location = Internal.getNodeFromName1(s, 'GridLocation') # localisation des donnees des receveurs
                        if location is not None: location = Internal.getValue(location)
                        if location == 'CellCenter': loc = 'centers'
                        else: loc = 'nodes'
                        Coefs     = idn[1]
                        DonorType = Internal.getNodeFromName1(s,'InterpolantsType')[1]
                        ListDonor = Internal.getNodeFromName1(s,'PointList')[1]
                        ListRcv   = Internal.getNodeFromName1(s,'PointListDonor')[1]
                        if sname == 'ID':
                            # print 'transfert ID: zd ', zd[0]
                            arrayT = connector._setInterpTransfersD(zd, variables, ListDonor, DonorType, Coefs, varType, compact,
                                                                    cellNVariable,
                                                                    Internal.__GridCoordinates__,
                                                                    Internal.__FlowSolutionNodes__,
                                                                    Internal.__FlowSolutionCenters__)
                            infos.append([dname,arrayT,ListRcv,loc])

                        elif sname == 'IB' and compactD and "gradxDensity" not in variablesIBC:
                            xPC = Internal.getNodeFromName1(s,'CoordinateX_PC')[1]
                            yPC = Internal.getNodeFromName1(s,'CoordinateY_PC')[1]
                            zPC = Internal.getNodeFromName1(s,'CoordinateZ_PC')[1]
                            xPW = Internal.getNodeFromName1(s,'CoordinateX_PW')[1]
                            yPW = Internal.getNodeFromName1(s,'CoordinateY_PW')[1]
                            zPW = Internal.getNodeFromName1(s,'CoordinateZ_PW')[1]
                            xPI = Internal.getNodeFromName1(s,'CoordinateX_PI')[1]
                            yPI = Internal.getNodeFromName1(s,'CoordinateY_PI')[1]
                            zPI = Internal.getNodeFromName1(s,'CoordinateZ_PI')[1]
                            Density = Internal.getNodeFromName1(s,'Density')[1]
                            arrayT = connector._setIBCTransfersD(zd, variablesIBC, ListDonor, DonorType, Coefs,
                                                                 xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI,
                                                                 Density,
                                                                 bcType, varType, compact, Gamma, Cv, MuS, Cs, Ts,
                                                                 Internal.__GridCoordinates__,
                                                                 Internal.__FlowSolutionNodes__,
                                                                 Internal.__FlowSolutionCenters__)
                            infos.append([dname,arrayT,ListRcv,loc])

                        elif sname == 'IB' and "gradxDensity" in variablesIBC:
                            pressure = Internal.getNodeFromName1(s,'Pressure')[1]
                            gradxP   = Internal.getNodeFromName1(s, 'gradxPressure')
                            gradyP   = Internal.getNodeFromName1(s, 'gradyPressure')
                            gradzP   = Internal.getNodeFromName1(s, 'gradzPressure')

                            if gradxP is not None:
                                gradxP = gradxP[1]
                                gradyP = gradyP[1]
                                gradzP = gradzP[1]
                            else:
                                nIBC = pressure.shape[0]
                                gradxPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                gradyPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                gradzPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                s[2].append(['gradxPressure' , gradxPressureNP , [], 'DataArray_t'])
                                s[2].append(['gradyPressure' , gradyPressureNP , [], 'DataArray_t'])
                                s[2].append(['gradzPressure' , gradzPressureNP , [], 'DataArray_t'])
                                gradxP = Internal.getNodeFromName1(s, 'gradxPressure')[1]
                                gradyP = Internal.getNodeFromName1(s, 'gradyPressure')[1]
                                gradzP = Internal.getNodeFromName1(s, 'gradzPressure')[1]

                            #print 'transfert IBC : zd ', zd[0]
                            arrayT = connector._setIBCTransfersD4GradP(zd, variablesIBC, ListDonor, DonorType, Coefs,
                                                                       pressure,
                                                                       gradxP, gradyP, gradzP,
                                                                       bcType, varType, compact, Gamma, Cv, MuS, Cs, Ts, alpha,
                                                                       Internal.__GridCoordinates__,
                                                                       Internal.__FlowSolutionNodes__,
                                                                       Internal.__FlowSolutionCenters__)
                            infos.append([dname,arrayT,ListRcv,loc])

    # Sortie
    return infos


#=============================================================================
#=============================================================================
# 1. Info chimere instationnaire
#=============================================================================
def getUnsteadyConnectInfos(t):
    """Return number of timelevel, min and max values present in the connectivity tree."""
    inst = {}
    numero_max =-100000000
    numero_min = 100000000

    zones = Internal.getZones(t)
    c        = 0
    for z in zones:
        subRegions  =  Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        for s in subRegions:
            #tri des pas de temps instationnaire
            if '#' in s[0]:
                numero_iter = int( s[0].split('#')[1].split('_')[0] )
                if numero_iter < numero_min : numero_min = numero_iter
                if numero_iter > numero_max : numero_max = numero_iter

                if numero_iter in inst.keys():
                    Noz = inst[ numero_iter ][0]
                    Noz = Noz + [c]
                    inst[ numero_iter ]=  [  Noz  ]
                else:
                    inst[ numero_iter ]= [ [c] ]

            TimeLevelNumber = len(inst)

    return [ TimeLevelNumber, numero_min, numero_max ]

#=============================================================================
#=============================================================================
# 4. INFO CHIMERE EXTRAITES
#=============================================================================
#=============================================================================


#=============================================================================
# Get overset info: if a cell is interpolated, orphan, extrapolated
#                   or the quality of the donor cell
# 'orphan' (0,1), 'extrapolated'(sum |cf|), 'interpolated'(0,1)
# 'cellRatio': max(volD/volR,volR/volD)
# 'donorAspect': edgeRatio of donor cell
# Compliant with setInterpData storage
#=============================================================================
def getOversetInfo(aR, topTreeD, loc='nodes', type='interpolated'):
    """Return information on overset connectivities."""
    tR = Internal.copyRef(aR)

    if type == 'interpolated':
        if loc == 'nodes': C._initVars(tR,'interpolated', 0.)
        else: C._initVars(tR,'centers:interpolated', 0.)
        tR = oversetNatureOfCells__(tR,topTreeD,'interpolated')
    elif type == 'extrapolated':
        if loc == 'nodes': C._initVars(tR,'extrapolated', 0.)
        else: C._initVars(tR,'centers:extrapolated', 0.)
        tR = oversetNatureOfCells__(tR,topTreeD,'extrapolated')
    elif type =='orphan':
        if loc == 'nodes': C._initVars(tR,'orphan', 0.)
        else: C._initVars(tR,'centers:orphan', 0.)
        tR = oversetNatureOfCells__(tR,topTreeD,'orphan')
    elif type =='cellRatio':
        if loc == 'nodes': C._initVars(tR,'cellRatio', 0.)
        else: C._initVars(tR,'centers:cellRatio', 0.)
        tR = oversetCellRatio__(tR,topTreeD)
    elif type == 'donorAspect':
        if loc == 'nodes': C._initVars(tR,'donorAspect', 0.)
        else: C._initVars(tR,'centers:donorAspect', 0.)
        tR = oversetDonorAspect__(tR,topTreeD)

    else: raise ValueError("getOversetInfo: type is invalid.")
    return tR


# =============================================================================
# extract info for interpolated/extrapolated/orphan cells
# Compliant with setInterpData
# =============================================================================
def oversetNatureOfCells__(aR,topTreeD,nature):
    tR = Internal.copyRef(aR)
    zonesR = Internal.getZones(tR)
    # First pass : direct storage
    for zr in zonesR:
        subRegions2 = Internal.getNodesFromType1(zr,'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0][0:2]
            if sname == 'ID':
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        for s in subRegions:
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Receiver': # direct storage ok
                location = Internal.getNodesFromName1(s, 'GridLocation')
                if location != []: location = Internal.getValue(location[0])
                locr = 'nodes'
                if location == 'CellCenter': locr = 'centers'

                if nature == 'interpolated':
                    ListRcv = Internal.getNodesFromName1(s, 'PointList')[0][1]
                    if ListRcv.size > 0:
                        field = Converter.array('interpolated', ListRcv.size,1,1)
                        field = Converter.initVars(field,'interpolated',1.)
                        zr = C.setPartialFields(zr, [field], [ListRcv], loc=locr)

                elif nature == 'extrapolated':
                    ListExtrap = Internal.getNodesFromName1(s, 'ExtrapPointList')
                    if ListExtrap != []:
                        ListRcv = Internal.getNodesFromName1(s, 'PointList')[0][1]
                        ListExtrap = ListExtrap[0][1]
                        DonorTypes = Internal.getNodesFromName1(s,'InterpolantsType')[0][1]
                        Coefs = Internal.getNodesFromName1(s,'InterpolantsDonor')[0][1]
                        # Somme des |coefs| : necessite ListRcv, ListExtrap, Coefs, DonorType
                        field = connector.getExtrapAbsCoefs(ListRcv, ListExtrap, DonorTypes, Coefs)
                        zr = C.setPartialFields(zr, [field], [ListExtrap], loc=locr)
                elif nature == 'orphan':
                    orphans = Internal.getNodesFromName(zr, 'OrphanPointList')
                    if orphans != []:
                        ListOrphan = orphans[0][1]
                        field = Converter.array('orphan', ListOrphan.size,1,1)
                        field = Converter.initVars(field, 'orphan', 1.)
                        zr = C.setPartialFields(zr, [field], [ListOrphan], loc=locr)

    # 2nd pass: storage in donor zones
    zones = Internal.getZones(topTreeD)
    for zd in zones:
        subRegions2 = Internal.getNodesFromType1(zd,'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0][0:2]
            if sname == 'ID':
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        for s in subRegions:
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Donor': # inverse storage ok
                location = Internal.getNodesFromName1(s, 'GridLocation')
                if location != []: location = Internal.getValue(location[0])
                locr = 'nodes'
                if location == 'CellCenter': locr = 'centers'

                zrcvname = Internal.getValue(s)
                zreceivers = Internal.getNodesFromName2(tR,zrcvname)
                zr = Internal.getNodesFromType1(zreceivers,'Zone_t')
                if zr != []:
                    zr = zr[0]
                    if nature == 'interpolated':
                        ListRcv = Internal.getNodesFromName1(s,'PointListDonor')[0][1]
                        if ListRcv.size > 0:
                            field = Converter.array('interpolated',ListRcv.size,1,1)
                            field = Converter.initVars(field,'interpolated',1.)
                            zr = C.setPartialFields(zr, [field], [ListRcv], loc=locr)
                    elif nature == 'extrapolated':
                        ListExtrap = Internal.getNodesFromName1(s,'ExtrapPointList')
                        if ListExtrap != []:
                            ListExtrap = ListExtrap[0][1]
                            DonorTypes = Internal.getNodesFromName1(s,'InterpolantsType')[0][1]
                            Coefs = Internal.getNodesFromName1(s,'InterpolantsDonor')[0][1]
                            ListRcv = Internal.getNodesFromName1(s,'PointListDonor')[0][1]
                            # Somme des |coefs| : necessite ListRcv, ListExtrap, Coefs, DonorType
                            field = connector.getExtrapAbsCoefs(ListRcv, ListExtrap, DonorTypes, Coefs)
                            zr = C.setPartialFields(zr, [field], [ListExtrap], loc=locr)
                    elif nature == 'orphan':
                        orphans = Internal.getNodesFromName(zd, 'OrphanPointList')
                        if orphans != []:
                            ListOrphan = orphans[0][1]
                            field = Converter.array('orphan',ListOrphan.size,1,1)
                            field = Converter.initVars(field,'orphan',1.)
                            zr = C.setPartialFields(zr, [field], [ListOrphan], loc=locr)
    return tR

#===============================================================================
# Information on overset quality : ratio of volumes between donors and receivers
# Compliant with setInterpolations
#===============================================================================
def oversetCellRatio__(aR, topTreeD):
    try: import Generator.PyTree as G
    except: raise ImportError("oversetInfo: requires Generator module.")

    tR = Internal.copyRef(aR)
    G._getVolumeMap(tR)
    tR = C.center2Node(tR,'centers:vol')

    tD = Internal.copyRef(topTreeD)
    G._getVolumeMap(tD)
    tD = C.center2Node(tD,'centers:vol')
    C._rmVars(tD,['centers:vol'])

    # First pass : direct storage
    zones = Internal.getZones(tR)
    for zr in zones:
        subRegions2 = Internal.getNodesFromType1(zr,'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0][0:2]
            if sname=='ID':
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        # parentr,dr = Internal.getParentOfNode(tR,zr)
        volNR = C.getField('vol', zr, api=1)[0][1]
        volCR = C.getField('centers:vol', zr, api=1)[0][1]
        for s in subRegions:
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Receiver': # direct storage ok
                # interpolated node or cell ?
                location = Internal.getNodesFromName1(s, 'GridLocation')
                if location != []: location = Internal.getValue(location[0])
                locr = 'nodes'; volRcv = []
                if location == 'CellCenter': locr = 'centers'; volRcv = volCR
                else: volRcv = volNR
                ListRcv = Internal.getNodesFromName1(s,'PointList')[0][1]
                # donor zone
                zdnrname = Internal.getValue(s)
                zdonors = Internal.getNodesFromName2(tD,zdnrname)
                zd = Internal.getNodesFromType1(zdonors,'Zone_t')
                if zd == []: raise ValueError("oversetInfo: donor zone %s not found."%zdnrname)
                else: zd = zd[0]
                volDnr = C.getField('vol', zd, api=1)[0][1]
                ListDnr = Internal.getNodesFromName1(s,'PointListDonor')[0][1]

                nindI = ListDnr.size
                if nindI > 0:
                    field = Converter.array('cellRatio',nindI,1,1)
                    for noind in range(nindI):
                        voldl = volDnr[0,ListDnr[noind]]
                        volrl = volRcv[0,ListRcv[noind]]
                        cr = max(voldl/volrl,volrl/voldl)
                        field[1][0,noind] = cr
                    zr = C.setPartialFields(zr, [field], [ListRcv], loc=locr)

    # 2nd pass : inverse storage
    zones = Internal.getZones(tD)
    for zd in zones:
        subRegions2 = Internal.getNodesFromType1(zd,'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0][0:2]
            if sname=='ID':
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        volDnr = C.getField('vol', zd, api=1)[0][1]
        for s in subRegions:
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Donor': # inverse storage ok
                zrcvname = Internal.getValue(s)
                zreceivers = Internal.getNodesFromName2(tR,zrcvname)
                zr = Internal.getNodesFromType1(zreceivers,'Zone_t')
                if zr != []:
                    zr = zr[0]
                    # parentr,dr = Internal.getParentOfNode(tR,zr)
                    # interpolated node or cell ?
                    location = Internal.getNodesFromName1(s, 'GridLocation')
                    if location != []: location = Internal.getValue(location[0])
                    locr = 'nodes'; volRcv = []
                    if location == 'CellCenter': locr = 'centers'; volRcv = C.getField('centers:vol', zr, api=1)[0][1]
                    else: volRcv = C.getField('vol', zr, api=1)[0][1]
                    ListRcv = Internal.getNodesFromName1(s, 'PointListDonor')[0][1]
                    ListDnr = Internal.getNodesFromName1(s, 'PointList')[0][1]
                    #
                    nindI = ListDnr.size
                    if nindI > 0:
                        field = Converter.array('cellRatio', nindI, 1, 1)
                        for noind in range(nindI):
                            voldl = volDnr[0,ListDnr[noind]]
                            volrl = volRcv[0,ListRcv[noind]]
                            cr = max(voldl/volrl,volrl/voldl)
                            field[1][0,noind] = cr
                        zr = C.setPartialFields(zr, [field], [ListRcv],loc=locr)
                        # parentr[2][dr] = zr

    C._rmVars(tR, 'centers:vol') # faut il la detruire ou non ? pas de test leger pour savoir si c'etait ds l arbre avant
    C._rmVars(tR, 'vol')
    C._rmVars(tD, 'vol')
    return tR

#===============================================================================
# Information on overset quality : ratio of volumes between donors and receivers
# Compliant with setInterpolations
#===============================================================================
def oversetDonorAspect__(aR, topTreeD):
    try: import Generator.PyTree as G
    except: raise ImportError("oversetInfo: requires Generator module.")
    tR = Internal.copyRef(aR)
    tD = Internal.copyRef(topTreeD)
    G._getEdgeRatio(tD)
    tD = C.center2Node(tD,'centers:EdgeRatio')
    C._rmVars(tD,['centers:EdgeRatio'])
    # First pass : direct storage
    zones = Internal.getZones(tR)
    for zr in zones:
        subRegions2 = Internal.getNodesFromType1(zr,'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0][0:2]
            if sname == 'ID':
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        # parentr,dr = Internal.getParentOfNode(tR,zr)
        for s in subRegions:
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Receiver': # direct storage ok
                # interpolated node or cell ?
                location = Internal.getNodesFromName1(s, 'GridLocation')
                if location != []: location = Internal.getValue(location[0])
                locr = 'nodes'
                if location == 'CellCenter': locr = 'centers';

                # donor zone
                zdnrname = Internal.getValue(s)
                zdonors = Internal.getNodesFromName2(tD,zdnrname)
                zd = Internal.getNodesFromType1(zdonors,'Zone_t')
                if zd == []: raise ValueError("oversetInfo: donor zone %s not found."%zdnrname)
                else: zd = zd[0]
                ER = C.getField('EdgeRatio', zd, api=1)[0][1]
                ListDnr = Internal.getNodesFromName1(s,'PointListDonor')[0][1]
                ListRcv = Internal.getNodesFromName1(s,'PointList')[0][1]
                nindI = ListDnr.size
                if nindI > 0:
                    field = Converter.array('donorAspect',nindI,1,1)
                    for noind in range(nindI):
                        field[1][0,noind] = ER[0,ListDnr[noind]]
                    zr = C.setPartialFields(zr, [field], [ListRcv],loc=locr)
        # parentr[2][dr] = zr

    # 2nd pass : inverse storage
    zones = Internal.getZones(tD)
    for zd in zones:
        subRegions2 = Internal.getNodesFromType1(zd,'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0][0:2]
            if sname == 'ID':
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        ER = C.getField('EdgeRatio',zd)[0][1]
        for s in subRegions:
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Donor': # inverse storage ok
                zrcvname = Internal.getValue(s)
                zreceivers = Internal.getNodesFromName2(tR,zrcvname)
                zr = Internal.getNodesFromType1(zreceivers,'Zone_t')
                if zr != []:
                    zr = zr[0]
                    # parentr,dr = Internal.getParentOfNode(tR,zr)
                    # interpolated node or cell ?
                    location = Internal.getNodesFromName1(s, 'GridLocation')
                    if location != []: location = Internal.getValue(location[0])
                    locr = 'nodes'
                    if location == 'CellCenter': locr = 'centers'
                    ListRcv = Internal.getNodesFromName1(s,'PointListDonor')[0][1]
                    ListDnr = Internal.getNodesFromName1(s,'PointList')[0][1]

                    nindI = ListDnr.size
                    if nindI > 0:
                        field = Converter.array('donorAspect',nindI,1,1)
                        for noind in range(nindI):
                            field[1][0,noind] = ER[0,ListDnr[noind]]
                        zr = C.setPartialFields(zr, [field], [ListRcv], loc=locr)
                        # parentr[2][dr] = zr

    C._rmVars(tD,'EdgeRatio')
    return tR


#=============================================================================
#=============================================================================
# 5. GESTION DU NOEUD OVERSETHOLES
#=============================================================================
#=============================================================================

#==============================================================================
# IN: ni, nj, nk
# IN: pointList: liste d'index i,j,k
# OUT: liste d'index ind = i + j*ni + k*ni*nj
#==============================================================================
def globalIndex__(ni, nj, nk, pointList):
    s = pointList.shape[1]
    a = numpy.empty((s), dtype=numpy.int32)
    nij = ni*nj
    a[:] = (pointList[0,:]-1) + (pointList[1,:]-1)*ni + (pointList[2,:]-1)*nij
    return a

#==============================================================================
# Cette fonction remplace le noeud oversetHoles
#==============================================================================
def cellN2OversetHoles(t, append=False):
    """Create OversetHole nodes from cellN field."""
    if append:
        try: import Compressor
        except: raise ImportError("cellN2OversetHoles: requires Compressor module.")
    a = Internal.copyRef(t)
    zones = Internal.getZones(a)
    for z in zones:
        # First find cellN
        cellN = C.getField('cellN', z, api=1)[0] # cellN en noeud
        loc = 'nodes'
        if cellN == []:
            loc = 'centers'
            cellN = C.getField('centers:cellN', z, api=1)[0] # cellN en centres
        if cellN != []:
            cellN = cellN[1]
            cellN = cellN.reshape(cellN.size, order='F' )
        # dim
        dim = Internal.getZoneDim(z); cellDim = dim[4]
        # Structured zone
        if dim[0] == 'Structured': # output im,jm,km
            im = dim[1]; jm = dim[2]; km = dim[3]
            if loc == 'centers': im = im-1; jm = jm-1; km = km-1
            nb = 0
            for i in range(cellN.size):
                if cellN[i] == 0: nb += 1
            pointList = numpy.empty((cellDim,nb), dtype=Internal.E_NpyInt, order='F' )
            connector.cellN2OversetHolesStruct(pointList, cellN, im, jm, cellDim, cellN.size)
        # Unstructured zone
        else:
            nb = 0
            for i in range(cellN.size):
                if cellN[i] == 0: nb += 1
            pointList = numpy.empty((1,nb), dtype=Internal.E_NpyInt)
            connector.cellN2OversetHolesUnStruct(pointList, cellN, cellN.size)
        if nb != 0:
            # Push OversetHoles Node
            GC = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
            if GC == []:
                z[2].append(['ZoneGridConnectivity', None, [],
                             'ZoneGridConnectivity_t'])
                info = z[2][len(z[2])-1]
            else: info = GC[0]

            OH = Internal.getNodesFromType1(info, 'OversetHoles_t')
            if OH == []:
                info[2].append(['OversetHoles', None, [], 'OversetHoles_t'])
                info = info[2][len(info[2])-1]
                if loc == 'centers':
                    Internal._createChild(info, 'GridLocation', 'GridLocation_t', value='CellCenter')
                info[2].append(['PointList', pointList, [], 'IndexArray_t'])
            else:
                PL = Internal.getNodesFromName(OH[0], 'PointList') # prev
                if PL == []:
                    OH[0][2].append(['PointList', pointList, [], 'IndexArray_t'])
                else:
                    if not append:
                        PL[0][1] = pointList
                        DL = Internal.getNodesFromName(OH[0], 'DeltaList')
                        if DL != []:
                            (p, c) = Internal.getParentOfNode(OH[0], DL[0])
                            del p[2][c]
                    else: # append
                        ref = PL[0][1]
                        PL[0][1] = pointList
                        if dim[0] == 'Structured': # conversion en indices globaux
                            im = dim[1]; jm = dim[2]; km = dim[3]
                            if loc == 'centers': im = im-1; jm = jm-1; km = km-1
                            ref = globalIndex__(im, jm, km, ref)
                            pointList = globalIndex__(im, jm, km, pointList)

                        delta = Compressor.deltaIndex(pointList, ref)
                        DL = Internal.getNodesFromName(OH[0], 'DeltaList')
                        if DL == []:
                            delta = numpy.concatenate([ref, delta])
                            OH[0][2].append(['DeltaList', delta, [],
                                             'IndexArray_t'])
                        else: # concatenate
                            prev = DL[0][1]
                            delta = numpy.concatenate([prev, delta])
                            DL[0][1] = delta
    return a

#==============================================================================
def setOversetHoles(z, indices):
    """Set an OversetHoles node in a zone for given mesh indices.
    Usage: setOversetHoles(z,ind)"""
    if z[3] != 'Zone_t': raise TypeError("setOversetHoles: only for a zone.")
    zp = Internal.copyRef(z)
    type = 0
    if isinstance(indices, numpy.ndarray): type = 1
    elif isinstance(indices, list):
        if indices == []: return zp
        if isinstance(indices[0], int): type = 2
        elif isinstance(indices[0], numpy.ndarray): type = 3
        elif isinstance(indices[0], list): type = 4
    else: raise TypeError("setOversetHoles: argument must be a list or a numpy of indices.")

    # cree la deltaList
    if type == 1: pointList = indices
    elif type == 2: pointList = numpy.array(indices, dtype=Internal.E_NpyInt)
    elif type == 3: pointList = numpy.concatenate(indices)
    elif type == 4:
        b = []
        for i in indices: b.append(numpy.array(indices, dtype=Internal.E_NpyInt))
        pointList = numpy.concatenate(b)

    # Push OversetHoles Node
    GC = Internal.getNodesFromType1(zp, 'ZoneGridConnectivity_t')
    if GC == []:
        zp[2].append(['ZoneGridConnectivity', None, [],
                      'ZoneGridConnectivity_t'])
        info = zp[2][len(z[2])-1]
    else: info = GC[0]
    OH = Internal.getNodesFromType(info, 'OversetHoles_t')
    if OH == []:
        info[2].append(['OversetHoles', None, [], 'OversetHoles_t'])
        info = info[2][len(info[2])-1]
        Internal._createChild(info, 'GridLocation', 'GridLocation_t', value='CellCenter')
        if type == 1 or type == 2:
            info[2].append(['PointList', pointList, [], 'IndexArray_t'])
        else: info[2].append(['DeltaList', pointList, [], 'IndexArray_t'])
    else:
        info = OH[0][2][len(info[2])-1]
        Internal._createChild(info, 'GridLocation', 'GridLocation_t', value='CellCenter')
        if type == 1 or type == 2:
            info[2].append(['PointList', pointList, [], 'IndexArray_t'])
        else: info[2].append(['DeltaList', pointList, [], 'IndexArray_t'])
    return zp

def extractChimeraInfo(a,type='interpolated',loc='centers'):
    """Extract interpolated/extrapolated/orphan points as zones.
    Usage: extractChimeraInfo(a,type,loc)"""
    FSCont = Internal.__FlowSolutionCenters__
    if loc == 'nodes': FSCont = Internal.__FlowSolutionNodes__
    try: import Post.PyTree as P
    except: raise ImportError("extractChimeraInfo requires Post.PyTree module.")

    typel = type.split('>')
    if len(typel) == 1:
        typel = type
        if typel == 'extrapolated':
            var = typel; prefix='ExtrapPts'
        elif typel == 'interpolated':
            var = typel; prefix='InterpPts'
        elif typel == 'orphan':
            var = typel; prefix='OrphanPts'
        else: raise AttributeError("extractChimeraInfo: not a valid type.")
        if loc == 'nodes': formula='{%s}>0.'%var
        else: formula='{centers:%s}>0.'%var

    elif len(typel)==2:
        val = float(typel[1])
        var='extrapolated'; prefix ='ExtrapCfPts'
        if loc == 'nodes': formula='{extrapolated}>%g'%val
        else: formula='{centers:extrapolated}>%g'%val
    else: raise AttributeError("extractChimeraInfo: not a valid type.")

    allChimPts=[]
    for z in Internal.getZones(a):
        FS = Internal.getNodeFromName1(z,FSCont)
        if FS is not None:
            chiminfo = Internal.getNodeFromName1(FS,var)
            if chiminfo is None:
                print('WARNING: extractChimeraInfo: chimera info %s cannot be extract for zone %s.'%(type,z[0]))
                print('You must apply oversetInfo or chimeraInfo to create the information in tree.')
            else:
                chimPts = P.selectCells(z,formula,strict=0)
                Internal._rmNodesByType(chimPts,"ZoneSubRegion_t")
                DIMS=Internal.getZoneDim(chimPts)
                if DIMS[1] > 0:
                    chimPts[0] = prefix+'_'+z[0]
                    allChimPts.append(chimPts)

    return allChimPts

#==============================================================================
# setIBCTransfersD for pressure gradients information (compact = 0)
# is secondOrder: update second order pressure gradients information as well
#==============================================================================
def _setIBCTransfersDForPressureGradients(topTreeD, ibctypes=[], secondOrder=False):
    variables = ["Pressure"]
    infos = []
    loc = 1

    for direction in ["x", "y", "z"]:
        var = "grad"+direction+"Pressure"
        variables.append(var)

    if secondOrder:
        for direction1 in ["x", "y", "z"]:
            for direction2 in ["x", "y", "z"]:
                var = "grad"+direction2+"grad"+direction1+"Pressure"
                variables.append(var)

    zonesD = Internal.getZones(topTreeD)
    for zd in zonesD:
        allIBCD = Internal.getNodesFromType(zd,"ZoneSubRegion_t")
        allIBCD = Internal.getNodesFromName(allIBCD,"IBCD_*")
        for IBCD in allIBCD:
            zrcvname = Internal.getValue(IBCD)
            ztype = int(IBCD[0].split("_")[1])
            if ibctypes != [] and ztype not in ibctypes: continue
            Coefs = Internal.getNodeFromName1(IBCD, 'InterpolantsDonor')[1]
            DonorType = Internal.getNodeFromName1(IBCD,'InterpolantsType')[1]
            ListDonor = Internal.getNodeFromName1(IBCD,'PointList')[1]
            ListRcv   = Internal.getNodeFromName1(IBCD,'PointListDonor')[1]

            Pressure = Internal.getNodeFromName1(IBCD,'Pressure')[1]
            nIBC = Pressure.shape[0]

            dictOfNumpys = {var: None for x in variables}

            for var in variables:
                field = Internal.getNodeFromName1(IBCD, var)
                if field is not None:
                    dictOfNumpys[var] = field[1]
                else:
                    fieldNP = numpy.zeros((nIBC),numpy.float64)
                    IBCD[2].append([var , fieldNP , [], 'DataArray_t'])
                    dictOfNumpys[var] = Internal.getNodeFromName1(IBCD, var)[1]

            if secondOrder:
                arrayT = connector._setIBCTransfersDForPressureGradientsOrder2(zd, variables, ListDonor, DonorType, Coefs,
                                                                               dictOfNumpys['Pressure'],
                                                                               dictOfNumpys['gradxPressure'] , dictOfNumpys['gradyPressure'] , dictOfNumpys['gradzPressure'],
                                                                               dictOfNumpys['gradxgradxPressure'], dictOfNumpys['gradxgradyPressure'], dictOfNumpys['gradxgradzPressure'],
                                                                               dictOfNumpys['gradygradxPressure'], dictOfNumpys['gradygradyPressure'], dictOfNumpys['gradygradzPressure'],
                                                                               dictOfNumpys['gradzgradxPressure'], dictOfNumpys['gradzgradyPressure'], dictOfNumpys['gradzgradzPressure'],
                                                                               Internal.__GridCoordinates__,
                                                                               Internal.__FlowSolutionNodes__,
                                                                               Internal.__FlowSolutionCenters__)
            else:
                arrayT = connector._setIBCTransfersDForPressureGradientsOrder1(zd, variables, ListDonor, DonorType, Coefs,
                                                                               dictOfNumpys['Pressure'],
                                                                               dictOfNumpys['gradxPressure'] , dictOfNumpys['gradyPressure'] , dictOfNumpys['gradzPressure'],
                                                                               Internal.__GridCoordinates__,
                                                                               Internal.__FlowSolutionNodes__,
                                                                               Internal.__FlowSolutionCenters__)

            infos.append([zrcvname,arrayT,ListRcv,loc])

    return infos

#==============================================================================
# setIBCTransfers for pressure gradients information (compact = 0)
# is secondOrder: update second order pressure gradients information as well
#==============================================================================
def _setIBCTransfersForPressureGradients(aR, topTreeD, ibctypes=[], secondOrder=False):

    variables = ["Pressure"]

    for direction in ["x", "y", "z"]:
        var = "grad"+direction+"Pressure"
        variables.append(var)

    if secondOrder:
        for direction1 in ["x", "y", "z"]:
            for direction2 in ["x", "y", "z"]:
                var = "grad"+direction2+"grad"+direction1+"Pressure"
                variables.append(var)

    loc = 1

    znr = {}
    zones = Internal.getZones(aR)
    for z in zones: znr[z[0]] = z

    zonesD = Internal.getZones(topTreeD)
    for zd in zonesD:
        allIBCD = Internal.getNodesFromType(zd,"ZoneSubRegion_t")
        allIBCD = Internal.getNodesFromName(allIBCD,"IBCD_*")
        for IBCD in allIBCD:
            zrcvname = Internal.getValue(IBCD)
            ztype = int(IBCD[0].split("_")[1])
            if ibctypes != [] and ztype not in ibctypes: continue
            Coefs = Internal.getNodeFromName1(IBCD, 'InterpolantsDonor')[1]
            DonorType = Internal.getNodeFromName1(IBCD,'InterpolantsType')[1]
            ListDonor = Internal.getNodeFromName1(IBCD,'PointList')[1]
            ListRcv   = Internal.getNodeFromName1(IBCD,'PointListDonor')[1]

            zr = znr.get(zrcvname, None)
            if zr is not None:
                Pressure = Internal.getNodeFromName1(IBCD,'Pressure')[1]
                nIBC = Pressure.shape[0]

                dictOfNumpys = {var: None for x in variables}

                for var in variables:
                    field = Internal.getNodeFromName1(IBCD, var)
                    if field is not None:
                        dictOfNumpys[var] = field[1]
                    else:
                        fieldNP = numpy.zeros((nIBC),numpy.float64)
                        IBCD[2].append([var , fieldNP , [], 'DataArray_t'])
                        dictOfNumpys[var] = Internal.getNodeFromName1(IBCD, var)[1]

                if secondOrder:
                    connector._setIBCTransfersForPressureGradientsOrder2(zr, zd, variables, ListRcv, ListDonor, DonorType, Coefs,
                                                                         dictOfNumpys['Pressure'],
                                                                         dictOfNumpys['gradxPressure'] , dictOfNumpys['gradyPressure'] , dictOfNumpys['gradzPressure'],
                                                                         dictOfNumpys['gradxgradxPressure'], dictOfNumpys['gradxgradyPressure'], dictOfNumpys['gradxgradzPressure'],
                                                                         dictOfNumpys['gradygradxPressure'], dictOfNumpys['gradygradyPressure'], dictOfNumpys['gradygradzPressure'],
                                                                         dictOfNumpys['gradzgradxPressure'], dictOfNumpys['gradzgradyPressure'], dictOfNumpys['gradzgradzPressure'],
                                                                         loc,
                                                                         Internal.__GridCoordinates__,
                                                                         Internal.__FlowSolutionNodes__,
                                                                         Internal.__FlowSolutionCenters__)
                else:
                    connector._setIBCTransfersForPressureGradientsOrder1(zr, zd, variables, ListRcv, ListDonor, DonorType, Coefs,
                                                                         dictOfNumpys['Pressure'],
                                                                         dictOfNumpys['gradxPressure'] , dictOfNumpys['gradyPressure'] , dictOfNumpys['gradzPressure'],
                                                                         loc,
                                                                         Internal.__GridCoordinates__,
                                                                         Internal.__FlowSolutionNodes__,
                                                                         Internal.__FlowSolutionCenters__)
    return None

# NOTE !!
## The functions below (_setIBCTransfers4GradP, _setIBCTransfers4GradP2, _setIBCTransfers4GradP3, _setIBCTransfers4FULLTBLE & _setIBCTransfers4FULLTBLE2)
## will be replaced by _setIBCTransfersDForPressureGradients & _setIBCTransfersForPressureGradients

#==============================================================================
# Copy of _setIBCTransfers for gradP info
# tc/tc2 -> RCV ZONES
#==============================================================================
def _setIBCTransfers4GradP(aR, topTreeD, variables=[], cellNVariable='',
                           variablesIBC=["Density", "Temperature", "gradxDensity", "gradyDensity", "gradzDensity", "gradxTemperature", "gradyTemperature", "gradzTemperature"],
                           bcType=3, varType=1, compact=0,
                           Gamma=1.4, Cv=1.7857142857142865, MuS=1.e-08,Cs=0.3831337844872463, Ts=1.0, alpha=1.0):

    # Recup des donnees a partir des zones donneuses
    # Dictionnaire pour optimisation
    znr = {}
    zones = Internal.getZones(aR)
    for z in zones: znr[z[0]] = z
    zonesD = Internal.getZones(topTreeD)
    for zd in zonesD:
        subRegions = Internal.getNodesFromType1(zd, 'ZoneSubRegion_t')
        for s in subRegions:
            sname = s[0][0:4]
            # test pour eviter parcours arbre inutile
            # if ((sname == 'IBCD' or sname == '2_IB') and variablesIBC is not None and not any(extension in s[0] for extension in ["_0_", "_1_", "_4_", "_5_"])):
            if ((sname == 'IBCD' or sname == '2_IB') and variablesIBC is not None):
                idn = Internal.getNodeFromName1(s, 'InterpolantsDonor')
                if idn is not None:
                    zoneRole = Internal.getNodeFromName2(s, 'ZoneRole')
                    zoneRole = Internal.getValue(zoneRole)
                    if zoneRole == 'Donor':
                        location = Internal.getNodeFromName1(s, 'GridLocation') #localisation des donnees des rcvr
                        if location is not None: location = Internal.getValue(location)
                        Coefs = idn[1]
                        DonorType = Internal.getNodeFromName1(s,'InterpolantsType')[1]
                        ListDonor = Internal.getNodeFromName1(s,'PointList')[1]
                        ListRcv   = Internal.getNodeFromName1(s,'PointListDonor')[1]
                        # Recup des champs du receveur
                        zrcvname = Internal.getValue(s)
                        ##zr = Internal.getNodesFromName2(aR, zrcvname)
                        zr = znr.get(zrcvname, None)
                        if zr is not None:
                            if location == 'CellCenter': loc = 1
                            else: loc = 0

                            Pressure = Internal.getNodeFromName1(s,'Pressure')[1]
                            gradxP   = Internal.getNodeFromName1(s, 'gradxPressure')
                            gradyP   = Internal.getNodeFromName1(s, 'gradyPressure')
                            gradzP   = Internal.getNodeFromName1(s, 'gradzPressure')

                            if gradxP is not None:
                                gradxP = gradxP[1]
                                gradyP = gradyP[1]
                                gradzP = gradzP[1]
                            else:
                                nIBC = Pressure.shape[0]
                                gradxPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                gradyPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                gradzPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                s[2].append(['gradxPressure' , gradxPressureNP , [], 'DataArray_t'])
                                s[2].append(['gradyPressure' , gradyPressureNP , [], 'DataArray_t'])
                                s[2].append(['gradzPressure' , gradzPressureNP , [], 'DataArray_t'])
                                gradxP = Internal.getNodeFromName1(s, 'gradxPressure')[1]
                                gradyP = Internal.getNodeFromName1(s, 'gradyPressure')[1]
                                gradzP = Internal.getNodeFromName1(s, 'gradzPressure')[1]

                            connector._setIBCTransfers4GradP(zr, zd, variablesIBC, ListRcv, ListDonor, DonorType, Coefs,
                                                             Pressure,
                                                             gradxP, gradyP, gradzP,
                                                             bcType, loc, varType, compact, Gamma, Cv, MuS, Cs, Ts, alpha,
                                                             Internal.__GridCoordinates__,
                                                             Internal.__FlowSolutionNodes__,
                                                             Internal.__FlowSolutionCenters__)
    return None

#==============================================================================
# Copy of _setIBCTransfers for gradP HO info
# tc/tc2 -> RCV ZONES
#==============================================================================
def _setIBCTransfers4GradP2(aR, topTreeD, variables=[], cellNVariable='',
                            variablesIBC=["Density", "Temperature", "gradxDensity", "gradyDensity", "gradzDensity", "gradxTemperature", "gradyTemperature", "gradzTemperature"],
                            bcType=3, varType=1, compact=0,
                            Gamma=1.4, Cv=1.7857142857142865, MuS=1.e-08,Cs=0.3831337844872463, Ts=1.0, alpha=1.0):

    # Recup des donnees a partir des zones donneuses
    # Dictionnaire pour optimisation
    znr = {}
    zones = Internal.getZones(aR)
    for z in zones: znr[z[0]] = z
    zonesD = Internal.getZones(topTreeD)
    for zd in zonesD:
        subRegions = Internal.getNodesFromType1(zd, 'ZoneSubRegion_t')
        for s in subRegions:
            sname = s[0][0:4]
            # test pour eviter parcours arbre inutile
            # if ((sname == 'IBCD' or sname == '2_IB') and variablesIBC is not None and not any(extension in s[0] for extension in ["_0_", "_1_", "_4_", "_5_"])):
            if ((sname == 'IBCD' or sname == '2_IB') and variablesIBC is not None):
                idn = Internal.getNodeFromName1(s, 'InterpolantsDonor')
                if idn is not None:
                    zoneRole = Internal.getNodeFromName2(s, 'ZoneRole')
                    zoneRole = Internal.getValue(zoneRole)
                    if zoneRole == 'Donor':
                        location = Internal.getNodeFromName1(s, 'GridLocation') #localisation des donnees des rcvr
                        if location is not None: location = Internal.getValue(location)
                        Coefs = idn[1]
                        DonorType = Internal.getNodeFromName1(s,'InterpolantsType')[1]
                        ListDonor = Internal.getNodeFromName1(s,'PointList')[1]
                        ListRcv   = Internal.getNodeFromName1(s,'PointListDonor')[1]
                        # Recup des champs du receveur
                        zrcvname = Internal.getValue(s)
                        ##zr = Internal.getNodesFromName2(aR, zrcvname)
                        zr = znr.get(zrcvname, None)
                        if zr is not None:
                            if location == 'CellCenter': loc = 1
                            else: loc = 0

                            Pressure = Internal.getNodeFromName1(s,'Pressure')[1]

                            gradxP   = Internal.getNodeFromName1(s, 'gradxPressure')
                            gradyP   = Internal.getNodeFromName1(s, 'gradyPressure')
                            gradzP   = Internal.getNodeFromName1(s, 'gradzPressure')

                            gradxxP   = Internal.getNodeFromName1(s, 'gradxgradxPressure')
                            gradxyP   = Internal.getNodeFromName1(s, 'gradxgradyPressure')
                            gradxzP   = Internal.getNodeFromName1(s, 'gradxgradzPressure')

                            gradyxP   = Internal.getNodeFromName1(s, 'gradygradxPressure')
                            gradyyP   = Internal.getNodeFromName1(s, 'gradygradyPressure')
                            gradyzP   = Internal.getNodeFromName1(s, 'gradygradzPressure')

                            gradzxP   = Internal.getNodeFromName1(s, 'gradzgradxPressure')
                            gradzyP   = Internal.getNodeFromName1(s, 'gradzgradyPressure')
                            gradzzP   = Internal.getNodeFromName1(s, 'gradzgradzPressure')

                            if gradxP is not None:
                                gradxP = gradxP[1]
                                gradyP = gradyP[1]
                                gradzP = gradzP[1]
                            else:
                                nIBC = Pressure.shape[0]
                                gradxPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                gradyPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                gradzPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                s[2].append(['gradxPressure' , gradxPressureNP , [], 'DataArray_t'])
                                s[2].append(['gradyPressure' , gradyPressureNP , [], 'DataArray_t'])
                                s[2].append(['gradzPressure' , gradzPressureNP , [], 'DataArray_t'])
                                gradxP = Internal.getNodeFromName1(s, 'gradxPressure')[1]
                                gradyP = Internal.getNodeFromName1(s, 'gradyPressure')[1]
                                gradzP = Internal.getNodeFromName1(s, 'gradzPressure')[1]

                            if gradxxP is not None:
                                gradxxP = gradxxP[1]
                                gradxyP = gradxyP[1]
                                gradxzP = gradxzP[1]
                            else:
                                nIBC = Pressure.shape[0]
                                gradxxPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                gradxyPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                gradxzPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                s[2].append(['gradxxPressure' , gradxxPressureNP , [], 'DataArray_t'])
                                s[2].append(['gradxyPressure' , gradxyPressureNP , [], 'DataArray_t'])
                                s[2].append(['gradxzPressure' , gradxzPressureNP , [], 'DataArray_t'])
                                gradxxP = Internal.getNodeFromName1(s, 'gradxxPressure')[1]
                                gradxyP = Internal.getNodeFromName1(s, 'gradxyPressure')[1]
                                gradxzP = Internal.getNodeFromName1(s, 'gradxzPressure')[1]

                            if gradyxP is not None:
                                gradyxP = gradyxP[1]
                                gradyyP = gradyyP[1]
                                gradyzP = gradyzP[1]
                            else:
                                nIBC = Pressure.shape[0]
                                gradyxPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                gradyyPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                gradyzPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                s[2].append(['gradyxPressure' , gradyxPressureNP , [], 'DataArray_t'])
                                s[2].append(['gradyyPressure' , gradyyPressureNP , [], 'DataArray_t'])
                                s[2].append(['gradyzPressure' , gradyzPressureNP , [], 'DataArray_t'])
                                gradyxP = Internal.getNodeFromName1(s, 'gradyxPressure')[1]
                                gradyyP = Internal.getNodeFromName1(s, 'gradyyPressure')[1]
                                gradyzP = Internal.getNodeFromName1(s, 'gradyzPressure')[1]

                            if gradzxP is not None:
                                gradzxP = gradzxP[1]
                                gradzyP = gradzyP[1]
                                gradzzP = gradzzP[1]
                            else:
                                nIBC = Pressure.shape[0]
                                gradzxPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                gradzyPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                gradzzPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                s[2].append(['gradzxPressure' , gradzxPressureNP , [], 'DataArray_t'])
                                s[2].append(['gradzyPressure' , gradzyPressureNP , [], 'DataArray_t'])
                                s[2].append(['gradzzPressure' , gradzzPressureNP , [], 'DataArray_t'])
                                gradzxP = Internal.getNodeFromName1(s, 'gradzxPressure')[1]
                                gradzyP = Internal.getNodeFromName1(s, 'gradzyPressure')[1]
                                gradzzP = Internal.getNodeFromName1(s, 'gradzzPressure')[1]

                            connector._setIBCTransfers4GradP2(zr, zd, variablesIBC, ListRcv, ListDonor, DonorType, Coefs,
                                                              Pressure,
                                                              gradxP, gradyP, gradzP,
                                                              gradxxP, gradxyP, gradxzP,
                                                              gradyxP, gradyyP, gradyzP,
                                                              gradzxP, gradzyP, gradzzP,
                                                              bcType, loc, varType, compact, Gamma, Cv, MuS, Cs, Ts, alpha,
                                                              Internal.__GridCoordinates__,
                                                              Internal.__FlowSolutionNodes__,
                                                              Internal.__FlowSolutionCenters__)
    return None

#==============================================================================
# Copy of _setIBCTransfers for gradP info
# RCV ZONES -> tc
# Parallel version + compact
#==============================================================================
def _setIBCTransfers4GradP3(zones, zonesD, graphInvIBCD, graphIBCD, procDict,
                            variablesIBC=['Density', 'Temperature', 'gradxDensity', 'gradyDensity', 'gradzDensity', 'gradxTemperature', 'gradyTemperature', 'gradzTemperature'],
                            bcType=0, varType=1, compact=1,
                            Gamma=1.4, Cv=1.7857142857142865, MuS=1.e-08,Cs=0.3831337844872463, Ts=1.0, alpha=1.0):

    import Converter.Mpi as Cmpi

    rank=Cmpi.rank
    datas = {}
    znr = {}
    for z in zones:
        znr[z[0]] = z

    for zd in zonesD:
        subRegions = Internal.getNodesFromType1(zd, 'ZoneSubRegion_t')
        for s in subRegions:
            sname = s[0][0:4]
            zname = s[0].split('_')[-1]
            dest = procDict[zname]
            if sname == 'IBCD' and dest != rank:
                Coefs     = numpy.copy(Internal.getNodeFromName1(s, 'InterpolantsDonor')[1])
                DonorType = numpy.copy(Internal.getNodeFromName1(s, 'InterpolantsType')[1])
                ListDonor = numpy.copy(Internal.getNodeFromName1(s, 'PointList')[1])
                ListRcv   = numpy.copy(Internal.getNodeFromName1(s, 'PointListDonor')[1])
                gradxP    = numpy.copy(Internal.getNodeFromName1(s, 'gradxPressure')[1])
                gradyP    = numpy.copy(Internal.getNodeFromName1(s, 'gradyPressure')[1])
                gradzP    = numpy.copy(Internal.getNodeFromName1(s, 'gradzPressure')[1])
                infos = [zd[0], zname, Coefs, DonorType, ListDonor, ListRcv, gradxP, gradyP, gradzP]
                rcvNode = dest
                if rcvNode not in datas: datas[rcvNode] = [infos]
                else: datas[rcvNode] += [infos]

    rcvDatas = Cmpi.sendRecv(datas, graphIBCD)

    datas = {}
    for dest in rcvDatas:
        for [name, zname, Coefs, DonorType, ListDonor, ListRcv, gradxP, gradyP, gradzP] in rcvDatas[dest]:
            zr = znr[zname]
            connector._setIBCTransfers4GradP3(zr, variablesIBC, ListRcv, ListDonor, DonorType, Coefs,
                                              gradxP, gradyP, gradzP,
                                              bcType, 1, varType, compact, Gamma, Cv, MuS, Cs, Ts, alpha,
                                              Internal.__GridCoordinates__,
                                              Internal.__FlowSolutionNodes__,
                                              Internal.__FlowSolutionCenters__)
            rcvNode = dest
            infos = [name, zname, gradxP, gradyP, gradzP]
            if rcvNode not in datas: datas[rcvNode] = [infos]
            else: datas[rcvNode] += [infos]

    rcvDatas = Cmpi.sendRecv(datas, graphInvIBCD)

    for dest in rcvDatas:
        for [name, zname, gradxP_new, gradyP_new, gradzP_new] in rcvDatas[dest]:
            for zd in zonesD:
                if zd[0] == name:
                    subRegions = Internal.getNodesFromType1(zd, 'ZoneSubRegion_t')
                    for s in subRegions:
                        sname = s[0][0:4]
                        znameD = s[0].split('_')[-1]
                        if sname == 'IBCD' and znameD == zname:
                            ListRcv   = Internal.getNodeFromName1(s, 'PointListDonor')[1]
                            gradxP    = Internal.getNodeFromName1(s, 'gradxPressure')[1]
                            gradyP    = Internal.getNodeFromName1(s, 'gradyPressure')[1]
                            gradzP    = Internal.getNodeFromName1(s, 'gradzPressure')[1]
                            connector._setIBCTransfers4GradP4(
                                gradxP_new, gradyP_new, gradzP_new,
                                gradxP, gradyP, gradzP)

    return None

#==============================================================================
# Copy of _setIBCTransfers for gradP + gradVelocity info
# tc/tc2 -> RCV ZONES
#==============================================================================

def _setIBCTransfers4FULLTBLE(aR, topTreeD, variables=[], cellNVariable='',
                              variablesIBC=['Density','MomentumX','MomentumY','MomentumZ','EnergyStagnationDensity'],
                              bcType=0, varType=1, storage=-1, compact=0,
                              Gamma=1.4, Cv=1.7857142857142865, MuS=1.e-08,Cs=0.3831337844872463, Ts=1.0, alpha=1.0):

    # Recup des donnees a partir des zones receveuses
    if storage < 1:
        # Dictionnaire pour optimisation
        znd = {}
        zones = Internal.getZones(topTreeD)
        for z in zones: znd[z[0]] = z

        zonesR = Internal.getZones(aR)
        for zr in zonesR:
            subRegions = Internal.getNodesFromType1(zr, 'ZoneSubRegion_t')
            for s in subRegions:
                sname = s[0][0:4]
                # test pour eviter parcours arbre inutile
                if ((sname == 'IBCD' or sname == '2_IB') and variablesIBC is not None and not any(extension in s[0] for extension in ["_0_", "_1_", "_4_", "_5_"])):
                    idn = Internal.getNodeFromName1(s, 'InterpolantsDonor')
                    if idn is not None: # la subRegion decrit des interpolations
                        zoneRole = Internal.getNodeFromName2(s, 'ZoneRole')
                        zoneRole = Internal.getValue(zoneRole)
                        if zoneRole == 'Receiver':
                            location = Internal.getNodeFromName1(s, 'GridLocation')
                            if location is not None: location = Internal.getValue(location)
                            Coefs = idn[1]
                            DonorType = Internal.getNodeFromName1(s,'InterpolantsType')[1]
                            ListRcv   = Internal.getNodeFromName1(s,'PointList')[1]
                            ListDonor = Internal.getNodeFromName1(s,'PointListDonor')[1]
                            # Recup des champs du receveur
                            zdnrname = Internal.getValue(s)
                            zd = znd[zdnrname]

                            if location == 'CellCenter': loc = 1
                            else: loc = 0

                            pressure  = Internal.getNodeFromName1(s,'Pressure')[1]
                            vx = Internal.getNodeFromName1(s,'VelocityX')
                            vy = Internal.getNodeFromName1(s,'VelocityY')
                            vz = Internal.getNodeFromName1(s,'VelocityZ')

                            gradxP   = Internal.getNodeFromName1(s, 'gradxPressure')
                            gradyP   = Internal.getNodeFromName1(s, 'gradyPressure')
                            gradzP   = Internal.getNodeFromName1(s, 'gradzPressure')

                            gradxU   = Internal.getNodeFromName1(s, 'gradxVelocityX')
                            gradyU   = Internal.getNodeFromName1(s, 'gradyVelocityX')
                            gradzU   = Internal.getNodeFromName1(s, 'gradzVelocityX')

                            gradxV   = Internal.getNodeFromName1(s, 'gradxVelocityY')
                            gradyV   = Internal.getNodeFromName1(s, 'gradyVelocityY')
                            gradzV   = Internal.getNodeFromName1(s, 'gradzVelocityY')

                            gradxW   = Internal.getNodeFromName1(s, 'gradxVelocityZ')
                            gradyW   = Internal.getNodeFromName1(s, 'gradyVelocityZ')
                            gradzW   = Internal.getNodeFromName1(s, 'gradzVelocityZ')


                            if vx is not None:
                                vx = vx[1]
                                vy = vy[1]
                                vz = vz[1]
                            else:
                                nIBC = pressure.shape[0]
                                vxNP  = numpy.zeros((nIBC),numpy.float64)
                                vyNP  = numpy.zeros((nIBC),numpy.float64)
                                vzNP  = numpy.zeros((nIBC),numpy.float64)
                                s[2].append(['VelocityX' , vxNP , [], 'DataArray_t'])
                                s[2].append(['VelocityY' , vyNP , [], 'DataArray_t'])
                                s[2].append(['VelocityZ' , vzNP , [], 'DataArray_t'])
                                vx = Internal.getNodeFromName1(s, 'VelocityX')[1]
                                vy = Internal.getNodeFromName1(s, 'VelocityY')[1]
                                vz = Internal.getNodeFromName1(s, 'VelocityZ')[1]

                            if gradxP is not None:
                                gradxP = gradxP[1]
                                gradyP = gradyP[1]
                                gradzP = gradzP[1]
                            else:
                                nIBC = pressure.shape[0]
                                gradxPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                gradyPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                gradzPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                s[2].append(['gradxPressure' , gradxPressureNP , [], 'DataArray_t'])
                                s[2].append(['gradyPressure' , gradyPressureNP , [], 'DataArray_t'])
                                s[2].append(['gradzPressure' , gradzPressureNP , [], 'DataArray_t'])
                                gradxP = Internal.getNodeFromName1(s, 'gradxPressure')[1]
                                gradyP = Internal.getNodeFromName1(s, 'gradyPressure')[1]
                                gradzP = Internal.getNodeFromName1(s, 'gradzPressure')[1]

                            if gradxU is not None:
                                gradxU = gradxU[1]
                                gradyU = gradyU[1]
                                gradzU = gradzU[1]
                            else:
                                nIBC = pressure.shape[0]
                                gradxVelocityXNP  = numpy.zeros((nIBC),numpy.float64)
                                gradyVelocityXNP  = numpy.zeros((nIBC),numpy.float64)
                                gradzVelocityXNP  = numpy.zeros((nIBC),numpy.float64)
                                s[2].append(['gradxVelocityX' , gradxVelocityXNP , [], 'DataArray_t'])
                                s[2].append(['gradyVelocityX' , gradyVelocityXNP , [], 'DataArray_t'])
                                s[2].append(['gradzVelocityX' , gradzVelocityXNP , [], 'DataArray_t'])
                                gradxU = Internal.getNodeFromName1(s, 'gradxVelocityX')[1]
                                gradyU = Internal.getNodeFromName1(s, 'gradyVelocityX')[1]
                                gradzU = Internal.getNodeFromName1(s, 'gradzVelocityX')[1]

                            if gradxV is not None:
                                gradxV = gradxV[1]
                                gradyV = gradyV[1]
                                gradzV = gradzV[1]
                            else:
                                nIBC = pressure.shape[0]
                                gradxVelocityYNP  = numpy.zeros((nIBC),numpy.float64)
                                gradyVelocityYNP  = numpy.zeros((nIBC),numpy.float64)
                                gradzVelocityYNP  = numpy.zeros((nIBC),numpy.float64)
                                s[2].append(['gradxVelocityY' , gradxVelocityYNP , [], 'DataArray_t'])
                                s[2].append(['gradyVelocityY' , gradyVelocityYNP , [], 'DataArray_t'])
                                s[2].append(['gradzVelocityY' , gradzVelocityYNP , [], 'DataArray_t'])
                                gradxV = Internal.getNodeFromName1(s, 'gradxVelocityY')[1]
                                gradyV = Internal.getNodeFromName1(s, 'gradyVelocityY')[1]
                                gradzV = Internal.getNodeFromName1(s, 'gradzVelocityY')[1]

                            if gradxW is not None:
                                gradxW = gradxW[1]
                                gradyW = gradyW[1]
                                gradzW = gradzW[1]
                            else:
                                nIBC = pressure.shape[0]
                                gradxVelocityZNP  = numpy.zeros((nIBC),numpy.float64)
                                gradyVelocityZNP  = numpy.zeros((nIBC),numpy.float64)
                                gradzVelocityZNP  = numpy.zeros((nIBC),numpy.float64)
                                s[2].append(['gradxVelocityZ' , gradxVelocityZNP , [], 'DataArray_t'])
                                s[2].append(['gradyVelocityZ' , gradyVelocityZNP , [], 'DataArray_t'])
                                s[2].append(['gradzVelocityZ' , gradzVelocityZNP , [], 'DataArray_t'])
                                gradxW = Internal.getNodeFromName1(s, 'gradxVelocityZ')[1]
                                gradyW = Internal.getNodeFromName1(s, 'gradyVelocityZ')[1]
                                gradzW = Internal.getNodeFromName1(s, 'gradzVelocityZ')[1]


                            #print 'transfert IBC : zr ', zr[0], ' et donor : ', zd[0]
                            connector._setIBCTransfers4FULLTBLE(zr, zd, variablesIBC, ListRcv, ListDonor, DonorType, Coefs,
                                                                pressure,
                                                                vx, vy, vz,
                                                                gradxP, gradyP, gradzP,
                                                                gradxU, gradyU, gradzU,
                                                                gradxV, gradyV, gradzV,
                                                                gradxW, gradyW, gradzW,
                                                                bcType, loc, varType, compact, Gamma, Cv, MuS, Cs, Ts, alpha,
                                                                Internal.__GridCoordinates__,
                                                                Internal.__FlowSolutionNodes__,
                                                                Internal.__FlowSolutionCenters__)

    # Recup des donnees a partir des zones donneuses
    if storage != 0:
        # Dictionnaire pour optimisation
        znr = {}
        zones = Internal.getZones(aR)
        for z in zones: znr[z[0]] = z
        zonesD = Internal.getZones(topTreeD)
        for zd in zonesD:
            subRegions = Internal.getNodesFromType1(zd, 'ZoneSubRegion_t')
            for s in subRegions:
                sname = s[0][0:4]
                # test pour eviter parcours arbre inutile
                if ((sname == 'IBCD' or sname == '2_IB') and variablesIBC is not None and not any(extension in s[0] for extension in ["_0_", "_1_", "_4_", "_5_"])):
                    idn = Internal.getNodeFromName1(s, 'InterpolantsDonor')
                    if idn is not None:
                        zoneRole = Internal.getNodeFromName2(s, 'ZoneRole')
                        zoneRole = Internal.getValue(zoneRole)
                        if zoneRole == 'Donor':
                            location = Internal.getNodeFromName1(s, 'GridLocation') #localisation des donnees des rcvr
                            if location is not None: location = Internal.getValue(location)
                            Coefs = idn[1]
                            DonorType = Internal.getNodeFromName1(s,'InterpolantsType')[1]
                            ListDonor = Internal.getNodeFromName1(s,'PointList')[1]
                            ListRcv   = Internal.getNodeFromName1(s,'PointListDonor')[1]
                            # Recup des champs du receveur
                            zrcvname = Internal.getValue(s)
                            ##zr = Internal.getNodesFromName2(aR, zrcvname)
                            zr = znr.get(zrcvname, None)
                            if zr is not None:
                                if location == 'CellCenter': loc = 1
                                else: loc = 0

                                pressure  = Internal.getNodeFromName1(s,'Pressure')[1]
                                vx = Internal.getNodeFromName1(s,'VelocityX')
                                vy = Internal.getNodeFromName1(s,'VelocityY')
                                vz = Internal.getNodeFromName1(s,'VelocityZ')

                                gradxP   = Internal.getNodeFromName1(s, 'gradxPressure')
                                gradyP   = Internal.getNodeFromName1(s, 'gradyPressure')
                                gradzP   = Internal.getNodeFromName1(s, 'gradzPressure')

                                gradxU   = Internal.getNodeFromName1(s, 'gradxVelocityX')
                                gradyU   = Internal.getNodeFromName1(s, 'gradyVelocityX')
                                gradzU   = Internal.getNodeFromName1(s, 'gradzVelocityX')

                                gradxV   = Internal.getNodeFromName1(s, 'gradxVelocityY')
                                gradyV   = Internal.getNodeFromName1(s, 'gradyVelocityY')
                                gradzV   = Internal.getNodeFromName1(s, 'gradzVelocityY')

                                gradxW   = Internal.getNodeFromName1(s, 'gradxVelocityZ')
                                gradyW   = Internal.getNodeFromName1(s, 'gradyVelocityZ')
                                gradzW   = Internal.getNodeFromName1(s, 'gradzVelocityZ')

                                if vx is not None:
                                    vx = vx[1]
                                    vy = vy[1]
                                    vz = vz[1]
                                else:
                                    nIBC = pressure.shape[0]
                                    vxNP  = numpy.zeros((nIBC),numpy.float64)
                                    vyNP  = numpy.zeros((nIBC),numpy.float64)
                                    vzNP  = numpy.zeros((nIBC),numpy.float64)
                                    s[2].append(['VelocityX' , vxNP , [], 'DataArray_t'])
                                    s[2].append(['VelocityY' , vyNP , [], 'DataArray_t'])
                                    s[2].append(['VelocityZ' , vzNP , [], 'DataArray_t'])
                                    vx = Internal.getNodeFromName1(s, 'VelocityX')[1]
                                    vy = Internal.getNodeFromName1(s, 'VelocityY')[1]
                                    vz = Internal.getNodeFromName1(s, 'VelocityZ')[1]

                                if gradxP is not None:
                                    gradxP = gradxP[1]
                                    gradyP = gradyP[1]
                                    gradzP = gradzP[1]
                                else:
                                    nIBC = pressure.shape[0]
                                    gradxPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                    gradyPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                    gradzPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                    s[2].append(['gradxPressure' , gradxPressureNP , [], 'DataArray_t'])
                                    s[2].append(['gradyPressure' , gradyPressureNP , [], 'DataArray_t'])
                                    s[2].append(['gradzPressure' , gradzPressureNP , [], 'DataArray_t'])
                                    gradxP = Internal.getNodeFromName1(s, 'gradxPressure')[1]
                                    gradyP = Internal.getNodeFromName1(s, 'gradyPressure')[1]
                                    gradzP = Internal.getNodeFromName1(s, 'gradzPressure')[1]

                                if gradxU is not None:
                                    gradxU = gradxU[1]
                                    gradyU = gradyU[1]
                                    gradzU = gradzU[1]
                                else:
                                    nIBC = pressure.shape[0]
                                    gradxVelocityXNP  = numpy.zeros((nIBC),numpy.float64)
                                    gradyVelocityXNP  = numpy.zeros((nIBC),numpy.float64)
                                    gradzVelocityXNP  = numpy.zeros((nIBC),numpy.float64)
                                    s[2].append(['gradxVelocityX' , gradxVelocityXNP , [], 'DataArray_t'])
                                    s[2].append(['gradyVelocityX' , gradyVelocityXNP , [], 'DataArray_t'])
                                    s[2].append(['gradzVelocityX' , gradzVelocityXNP , [], 'DataArray_t'])
                                    gradxU = Internal.getNodeFromName1(s, 'gradxVelocityX')[1]
                                    gradyU = Internal.getNodeFromName1(s, 'gradyVelocityX')[1]
                                    gradzU = Internal.getNodeFromName1(s, 'gradzVelocityX')[1]

                                if gradxV is not None:
                                    gradxV = gradxV[1]
                                    gradyV = gradyV[1]
                                    gradzV = gradzV[1]
                                else:
                                    nIBC = pressure.shape[0]
                                    gradxVelocityYNP  = numpy.zeros((nIBC),numpy.float64)
                                    gradyVelocityYNP  = numpy.zeros((nIBC),numpy.float64)
                                    gradzVelocityYNP  = numpy.zeros((nIBC),numpy.float64)
                                    s[2].append(['gradxVelocityY' , gradxVelocityYNP , [], 'DataArray_t'])
                                    s[2].append(['gradyVelocityY' , gradyVelocityYNP , [], 'DataArray_t'])
                                    s[2].append(['gradzVelocityY' , gradzVelocityYNP , [], 'DataArray_t'])
                                    gradxV = Internal.getNodeFromName1(s, 'gradxVelocityY')[1]
                                    gradyV = Internal.getNodeFromName1(s, 'gradyVelocityY')[1]
                                    gradzV = Internal.getNodeFromName1(s, 'gradzVelocityY')[1]

                                if gradxW is not None:
                                    gradxW = gradxW[1]
                                    gradyW = gradyW[1]
                                    gradzW = gradzW[1]
                                else:
                                    nIBC = pressure.shape[0]
                                    gradxVelocityZNP  = numpy.zeros((nIBC),numpy.float64)
                                    gradyVelocityZNP  = numpy.zeros((nIBC),numpy.float64)
                                    gradzVelocityZNP  = numpy.zeros((nIBC),numpy.float64)
                                    s[2].append(['gradxVelocityZ' , gradxVelocityZNP , [], 'DataArray_t'])
                                    s[2].append(['gradyVelocityZ' , gradyVelocityZNP , [], 'DataArray_t'])
                                    s[2].append(['gradzVelocityZ' , gradzVelocityZNP , [], 'DataArray_t'])
                                    gradxW = Internal.getNodeFromName1(s, 'gradxVelocityZ')[1]
                                    gradyW = Internal.getNodeFromName1(s, 'gradyVelocityZ')[1]
                                    gradzW = Internal.getNodeFromName1(s, 'gradzVelocityZ')[1]


                                # print 'transfert IBC : zr ', zr[0], ' et donor : ', zd[0]
                                connector._setIBCTransfers4FULLTBLE(zr, zd, variablesIBC, ListRcv, ListDonor, DonorType, Coefs,
                                                                    pressure,
                                                                    vx, vy, vz,
                                                                    gradxP, gradyP, gradzP,
                                                                    gradxU, gradyU, gradzU,
                                                                    gradxV, gradyV, gradzV,
                                                                    gradxW, gradyW, gradzW,
                                                                    bcType, loc, varType, compact, Gamma, Cv, MuS, Cs, Ts, alpha,
                                                                    Internal.__GridCoordinates__,
                                                                    Internal.__FlowSolutionNodes__,
                                                                    Internal.__FlowSolutionCenters__)
    return None

#==============================================================================
# Copy of _setIBCTransfers for gradP + gradVelocity info
# RCV ZONES -> tc
#==============================================================================

def _setIBCTransfers4FULLTBLE2(aR, topTreeD, variables=[], cellNVariable='',
                               variablesIBC=['Density','MomentumX','MomentumY','MomentumZ','EnergyStagnationDensity'],
                               bcType=0, varType=1, storage=-1, compact=0,
                               Gamma=1.4, Cv=1.7857142857142865, MuS=1.e-08,Cs=0.3831337844872463, Ts=1.0, alpha=1.0):

    # Recup des donnees a partir des zones receveuses
    if storage < 1:
        # Dictionnaire pour optimisation
        znd = {}
        zones = Internal.getZones(topTreeD)
        for z in zones: znd[z[0]] = z

        zonesR = Internal.getZones(aR)
        for zr in zonesR:
            subRegions = Internal.getNodesFromType1(zr, 'ZoneSubRegion_t')
            for s in subRegions:
                sname = s[0][0:4]
                # test pour eviter parcours arbre inutile
                if ((sname == 'IBCD' or sname == '2_IB') and variablesIBC is not None and not any(extension in s[0] for extension in ["_0_", "_1_", "_4_", "_5_"])):
                    idn = Internal.getNodeFromName1(s, 'InterpolantsDonor')
                    if idn is not None: # la subRegion decrit des interpolations
                        zoneRole = Internal.getNodeFromName2(s, 'ZoneRole')
                        zoneRole = Internal.getValue(zoneRole)
                        if zoneRole == 'Receiver':
                            location = Internal.getNodeFromName1(s, 'GridLocation')
                            if location is not None: location = Internal.getValue(location)
                            Coefs = idn[1]
                            DonorType = Internal.getNodeFromName1(s,'InterpolantsType')[1]
                            ListRcv   = Internal.getNodeFromName1(s,'PointList')[1]
                            ListDonor = Internal.getNodeFromName1(s,'PointListDonor')[1]
                            # Recup des champs du receveur
                            zdnrname = Internal.getValue(s)
                            zd = znd[zdnrname]

                            if location == 'CellCenter': loc = 1
                            else: loc = 0

                            pressure  = Internal.getNodeFromName1(s,'Pressure')[1]
                            vx = Internal.getNodeFromName1(s,'VelocityX')
                            vy = Internal.getNodeFromName1(s,'VelocityY')
                            vz = Internal.getNodeFromName1(s,'VelocityZ')

                            gradxP   = Internal.getNodeFromName1(s, 'gradxPressure')
                            gradyP   = Internal.getNodeFromName1(s, 'gradyPressure')
                            gradzP   = Internal.getNodeFromName1(s, 'gradzPressure')

                            gradxU   = Internal.getNodeFromName1(s, 'gradxVelocityX')
                            gradyU   = Internal.getNodeFromName1(s, 'gradyVelocityX')
                            gradzU   = Internal.getNodeFromName1(s, 'gradzVelocityX')

                            gradxV   = Internal.getNodeFromName1(s, 'gradxVelocityY')
                            gradyV   = Internal.getNodeFromName1(s, 'gradyVelocityY')
                            gradzV   = Internal.getNodeFromName1(s, 'gradzVelocityY')

                            gradxW   = Internal.getNodeFromName1(s, 'gradxVelocityZ')
                            gradyW   = Internal.getNodeFromName1(s, 'gradyVelocityZ')
                            gradzW   = Internal.getNodeFromName1(s, 'gradzVelocityZ')


                            if vx is not None:
                                vx = vx[1]
                                vy = vy[1]
                                vz = vz[1]
                            else:
                                nIBC = pressure.shape[0]
                                vxNP  = numpy.zeros((nIBC),numpy.float64)
                                vyNP  = numpy.zeros((nIBC),numpy.float64)
                                vzNP  = numpy.zeros((nIBC),numpy.float64)
                                s[2].append(['VelocityX' , vxNP , [], 'DataArray_t'])
                                s[2].append(['VelocityY' , vyNP , [], 'DataArray_t'])
                                s[2].append(['VelocityZ' , vzNP , [], 'DataArray_t'])
                                vx = Internal.getNodeFromName1(s, 'VelocityX')[1]
                                vy = Internal.getNodeFromName1(s, 'VelocityY')[1]
                                vz = Internal.getNodeFromName1(s, 'VelocityZ')[1]

                            if gradxP is not None:
                                gradxP = gradxP[1]
                                gradyP = gradyP[1]
                                gradzP = gradzP[1]
                            else:
                                nIBC = pressure.shape[0]
                                gradxPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                gradyPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                gradzPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                s[2].append(['gradxPressure' , gradxPressureNP , [], 'DataArray_t'])
                                s[2].append(['gradyPressure' , gradyPressureNP , [], 'DataArray_t'])
                                s[2].append(['gradzPressure' , gradzPressureNP , [], 'DataArray_t'])
                                gradxP = Internal.getNodeFromName1(s, 'gradxPressure')[1]
                                gradyP = Internal.getNodeFromName1(s, 'gradyPressure')[1]
                                gradzP = Internal.getNodeFromName1(s, 'gradzPressure')[1]

                            if gradxU is not None:
                                gradxU = gradxU[1]
                                gradyU = gradyU[1]
                                gradzU = gradzU[1]
                            else:
                                nIBC = pressure.shape[0]
                                gradxVelocityXNP  = numpy.zeros((nIBC),numpy.float64)
                                gradyVelocityXNP  = numpy.zeros((nIBC),numpy.float64)
                                gradzVelocityXNP  = numpy.zeros((nIBC),numpy.float64)
                                s[2].append(['gradxVelocityX' , gradxVelocityXNP , [], 'DataArray_t'])
                                s[2].append(['gradyVelocityX' , gradyVelocityXNP , [], 'DataArray_t'])
                                s[2].append(['gradzVelocityX' , gradzVelocityXNP , [], 'DataArray_t'])
                                gradxU = Internal.getNodeFromName1(s, 'gradxVelocityX')[1]
                                gradyU = Internal.getNodeFromName1(s, 'gradyVelocityX')[1]
                                gradzU = Internal.getNodeFromName1(s, 'gradzVelocityX')[1]

                            if gradxV is not None:
                                gradxV = gradxV[1]
                                gradyV = gradyV[1]
                                gradzV = gradzV[1]
                            else:
                                nIBC = pressure.shape[0]
                                gradxVelocityYNP  = numpy.zeros((nIBC),numpy.float64)
                                gradyVelocityYNP  = numpy.zeros((nIBC),numpy.float64)
                                gradzVelocityYNP  = numpy.zeros((nIBC),numpy.float64)
                                s[2].append(['gradxVelocityY' , gradxVelocityYNP , [], 'DataArray_t'])
                                s[2].append(['gradyVelocityY' , gradyVelocityYNP , [], 'DataArray_t'])
                                s[2].append(['gradzVelocityY' , gradzVelocityYNP , [], 'DataArray_t'])
                                gradxV = Internal.getNodeFromName1(s, 'gradxVelocityY')[1]
                                gradyV = Internal.getNodeFromName1(s, 'gradyVelocityY')[1]
                                gradzV = Internal.getNodeFromName1(s, 'gradzVelocityY')[1]

                            if gradxW is not None:
                                gradxW = gradxW[1]
                                gradyW = gradyW[1]
                                gradzW = gradzW[1]
                            else:
                                nIBC = pressure.shape[0]
                                gradxVelocityZNP  = numpy.zeros((nIBC),numpy.float64)
                                gradyVelocityZNP  = numpy.zeros((nIBC),numpy.float64)
                                gradzVelocityZNP  = numpy.zeros((nIBC),numpy.float64)
                                s[2].append(['gradxVelocityZ' , gradxVelocityZNP , [], 'DataArray_t'])
                                s[2].append(['gradyVelocityZ' , gradyVelocityZNP , [], 'DataArray_t'])
                                s[2].append(['gradzVelocityZ' , gradzVelocityZNP , [], 'DataArray_t'])
                                gradxW = Internal.getNodeFromName1(s, 'gradxVelocityZ')[1]
                                gradyW = Internal.getNodeFromName1(s, 'gradyVelocityZ')[1]
                                gradzW = Internal.getNodeFromName1(s, 'gradzVelocityZ')[1]


                            #print 'transfert IBC : zr ', zr[0], ' et donor : ', zd[0]
                            connector._setIBCTransfers4FULLTBLE2(zr, zd, variablesIBC, ListRcv, ListDonor, DonorType, Coefs,
                                                                 pressure,
                                                                 vx, vy, vz,
                                                                 gradxP, gradyP, gradzP,
                                                                 gradxU, gradyU, gradzU,
                                                                 gradxV, gradyV, gradzV,
                                                                 gradxW, gradyW, gradzW,
                                                                 bcType, loc, varType, compact, Gamma, Cv, MuS, Cs, Ts, alpha,
                                                                 Internal.__GridCoordinates__,
                                                                 Internal.__FlowSolutionNodes__,
                                                                 Internal.__FlowSolutionCenters__)

    # Recup des donnees a partir des zones donneuses
    if storage != 0:
        # Dictionnaire pour optimisation
        znr = {}
        zones = Internal.getZones(aR)
        for z in zones: znr[z[0]] = z
        zonesD = Internal.getZones(topTreeD)
        for zd in zonesD:
            subRegions = Internal.getNodesFromType1(zd, 'ZoneSubRegion_t')
            for s in subRegions:
                sname = s[0][0:4]
                # test pour eviter parcours arbre inutile
                if ((sname == 'IBCD' or sname == '2_IB') and variablesIBC is not None and not any(extension in s[0] for extension in ["_0_", "_1_", "_4_", "_5_"])):
                    idn = Internal.getNodeFromName1(s, 'InterpolantsDonor')
                    if idn is not None:
                        zoneRole = Internal.getNodeFromName2(s, 'ZoneRole')
                        zoneRole = Internal.getValue(zoneRole)
                        if zoneRole == 'Donor':
                            location = Internal.getNodeFromName1(s, 'GridLocation') #localisation des donnees des rcvr
                            if location is not None: location = Internal.getValue(location)
                            Coefs = idn[1]
                            DonorType = Internal.getNodeFromName1(s,'InterpolantsType')[1]
                            ListDonor = Internal.getNodeFromName1(s,'PointList')[1]
                            ListRcv   = Internal.getNodeFromName1(s,'PointListDonor')[1]
                            # Recup des champs du receveur
                            zrcvname = Internal.getValue(s)
                            ##zr = Internal.getNodesFromName2(aR, zrcvname)
                            zr = znr.get(zrcvname, None)
                            if zr is not None:
                                if location == 'CellCenter': loc = 1
                                else: loc = 0

                                pressure  = Internal.getNodeFromName1(s,'Pressure')[1]
                                vx = Internal.getNodeFromName1(s,'VelocityX')
                                vy = Internal.getNodeFromName1(s,'VelocityY')
                                vz = Internal.getNodeFromName1(s,'VelocityZ')

                                gradxP   = Internal.getNodeFromName1(s, 'gradxPressure')
                                gradyP   = Internal.getNodeFromName1(s, 'gradyPressure')
                                gradzP   = Internal.getNodeFromName1(s, 'gradzPressure')

                                gradxU   = Internal.getNodeFromName1(s, 'gradxVelocityX')
                                gradyU   = Internal.getNodeFromName1(s, 'gradyVelocityX')
                                gradzU   = Internal.getNodeFromName1(s, 'gradzVelocityX')

                                gradxV   = Internal.getNodeFromName1(s, 'gradxVelocityY')
                                gradyV   = Internal.getNodeFromName1(s, 'gradyVelocityY')
                                gradzV   = Internal.getNodeFromName1(s, 'gradzVelocityY')

                                gradxW   = Internal.getNodeFromName1(s, 'gradxVelocityZ')
                                gradyW   = Internal.getNodeFromName1(s, 'gradyVelocityZ')
                                gradzW   = Internal.getNodeFromName1(s, 'gradzVelocityZ')

                                if vx is not None:
                                    vx = vx[1]
                                    vy = vy[1]
                                    vz = vz[1]
                                else:
                                    nIBC = pressure.shape[0]
                                    vxNP  = numpy.zeros((nIBC),numpy.float64)
                                    vyNP  = numpy.zeros((nIBC),numpy.float64)
                                    vzNP  = numpy.zeros((nIBC),numpy.float64)
                                    s[2].append(['VelocityX' , vxNP , [], 'DataArray_t'])
                                    s[2].append(['VelocityY' , vyNP , [], 'DataArray_t'])
                                    s[2].append(['VelocityZ' , vzNP , [], 'DataArray_t'])
                                    vx = Internal.getNodeFromName1(s, 'VelocityX')[1]
                                    vy = Internal.getNodeFromName1(s, 'VelocityY')[1]
                                    vz = Internal.getNodeFromName1(s, 'VelocityZ')[1]

                                if gradxP is not None:
                                    gradxP = gradxP[1]
                                    gradyP = gradyP[1]
                                    gradzP = gradzP[1]
                                else:
                                    nIBC = pressure.shape[0]
                                    gradxPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                    gradyPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                    gradzPressureNP  = numpy.zeros((nIBC),numpy.float64)
                                    s[2].append(['gradxPressure' , gradxPressureNP , [], 'DataArray_t'])
                                    s[2].append(['gradyPressure' , gradyPressureNP , [], 'DataArray_t'])
                                    s[2].append(['gradzPressure' , gradzPressureNP , [], 'DataArray_t'])
                                    gradxP = Internal.getNodeFromName1(s, 'gradxPressure')[1]
                                    gradyP = Internal.getNodeFromName1(s, 'gradyPressure')[1]
                                    gradzP = Internal.getNodeFromName1(s, 'gradzPressure')[1]

                                if gradxU is not None:
                                    gradxU = gradxU[1]
                                    gradyU = gradyU[1]
                                    gradzU = gradzU[1]
                                else:
                                    nIBC = pressure.shape[0]
                                    gradxVelocityXNP  = numpy.zeros((nIBC),numpy.float64)
                                    gradyVelocityXNP  = numpy.zeros((nIBC),numpy.float64)
                                    gradzVelocityXNP  = numpy.zeros((nIBC),numpy.float64)
                                    s[2].append(['gradxVelocityX' , gradxVelocityXNP , [], 'DataArray_t'])
                                    s[2].append(['gradyVelocityX' , gradyVelocityXNP , [], 'DataArray_t'])
                                    s[2].append(['gradzVelocityX' , gradzVelocityXNP , [], 'DataArray_t'])
                                    gradxU = Internal.getNodeFromName1(s, 'gradxVelocityX')[1]
                                    gradyU = Internal.getNodeFromName1(s, 'gradyVelocityX')[1]
                                    gradzU = Internal.getNodeFromName1(s, 'gradzVelocityX')[1]

                                if gradxV is not None:
                                    gradxV = gradxV[1]
                                    gradyV = gradyV[1]
                                    gradzV = gradzV[1]
                                else:
                                    nIBC = pressure.shape[0]
                                    gradxVelocityYNP  = numpy.zeros((nIBC),numpy.float64)
                                    gradyVelocityYNP  = numpy.zeros((nIBC),numpy.float64)
                                    gradzVelocityYNP  = numpy.zeros((nIBC),numpy.float64)
                                    s[2].append(['gradxVelocityY' , gradxVelocityYNP , [], 'DataArray_t'])
                                    s[2].append(['gradyVelocityY' , gradyVelocityYNP , [], 'DataArray_t'])
                                    s[2].append(['gradzVelocityY' , gradzVelocityYNP , [], 'DataArray_t'])
                                    gradxV = Internal.getNodeFromName1(s, 'gradxVelocityY')[1]
                                    gradyV = Internal.getNodeFromName1(s, 'gradyVelocityY')[1]
                                    gradzV = Internal.getNodeFromName1(s, 'gradzVelocityY')[1]

                                if gradxW is not None:
                                    gradxW = gradxW[1]
                                    gradyW = gradyW[1]
                                    gradzW = gradzW[1]
                                else:
                                    nIBC = pressure.shape[0]
                                    gradxVelocityZNP  = numpy.zeros((nIBC),numpy.float64)
                                    gradyVelocityZNP  = numpy.zeros((nIBC),numpy.float64)
                                    gradzVelocityZNP  = numpy.zeros((nIBC),numpy.float64)
                                    s[2].append(['gradxVelocityZ' , gradxVelocityZNP , [], 'DataArray_t'])
                                    s[2].append(['gradyVelocityZ' , gradyVelocityZNP , [], 'DataArray_t'])
                                    s[2].append(['gradzVelocityZ' , gradzVelocityZNP , [], 'DataArray_t'])
                                    gradxW = Internal.getNodeFromName1(s, 'gradxVelocityZ')[1]
                                    gradyW = Internal.getNodeFromName1(s, 'gradyVelocityZ')[1]
                                    gradzW = Internal.getNodeFromName1(s, 'gradzVelocityZ')[1]


                                # print 'transfert IBC : zr ', zr[0], ' et donor : ', zd[0]
                                connector._setIBCTransfers4FULLTBLE2(zr, zd, variablesIBC, ListRcv, ListDonor, DonorType, Coefs,
                                                                     pressure,
                                                                     vx, vy, vz,
                                                                     gradxP, gradyP, gradzP,
                                                                     gradxU, gradyU, gradzU,
                                                                     gradxV, gradyV, gradzV,
                                                                     gradxW, gradyW, gradzW,
                                                                     bcType, loc, varType, compact, Gamma, Cv, MuS, Cs, Ts, alpha,
                                                                     Internal.__GridCoordinates__,
                                                                     Internal.__FlowSolutionNodes__,
                                                                     Internal.__FlowSolutionCenters__)
    return None

#==============================================================================
# Cette fonction retourne la transformation vectorielle d'une zone a l'autre
#==============================================================================
def getTransfo(zdonor, zrcv):

    transfo = numpy.zeros(3, dtype=numpy.int32)
    a = C.getFields(Internal.__GridCoordinates__, zdonor, api=1)[0]
    ni = a[2]; nj=a[3]; nk=a[4]

    if nk == 1: #2D
        i = ni//2; j = nj//2; k = nk//2
        ip1 = max(i+1,ni-1); jp1 = max(j+1,nj-1); kp1 = max(k+1,nk-1)
        ind = i + j*ni + k*ni*nj
        P0 = [ a[1][0][ind], a[1][1][ind],0.0 ]
        ind = ip1 + j*ni + k*ni*nj
        P1 = [ a[1][0][ind], a[1][1][ind],0.0 ]
        ind = i + jp1*ni + k*ni*nj
        P2 = [ a[1][0][ind], a[1][1][ind],0.0 ]
        l1 = Vector.sub(P1,P0)
        l2 = Vector.sub(P2,P0)
        l1=Vector.normalize(l1)
        l2=Vector.normalize(l2)
        x=[1.0,0.0,0.0]
        y=[0.0,1.0,0.0]
        a=Vector.dot(l1,x);b=Vector.dot(l1,y)
        c=Vector.dot(l2,x);d=Vector.dot(l2,y)
        mat_ = numpy.array([[a,b],
                            [c,d]])

        a = C.getFields(Internal.__GridCoordinates__, zrcv, api=1)[0]
        ni = a[2]; nj=a[3]; nk=a[4]
        i = ni//2; j = nj//2; k = nk//2
        ip1 = max(i+1,ni-1); jp1 = max(j+1,nj-1); kp1 = max(k+1,nk-1)
        #ip1=i+1;jp1=j+1
        ind = i + j*ni + k*ni*nj
        P0_  = [ a[1][0][ind], a[1][1][ind],0.0 ]
        ind = ip1 + j*ni + k*ni*nj
        P1_ = [ a[1][0][ind], a[1][1][ind],0.0 ]
        ind = i + jp1*ni + k*ni*nj
        P2_ = [ a[1][0][ind], a[1][1][ind],0.0 ]
        l1rcv = Vector.sub(P1_,P0_)
        l2rcv = Vector.sub(P2_,P0_)
        l1rcv=Vector.normalize(l1rcv)
        l2rcv=Vector.normalize(l2rcv)
        x=[1.0,0.0,0.0]
        y=[0.0,1.0,0.0]
        a=Vector.dot(l1rcv,x);b=Vector.dot(l1rcv,y)
        c=Vector.dot(l2rcv,x);d=Vector.dot(l2rcv,y)
        mat = numpy.array([[a,b],
                          [c,d]])
        mat=mat.T
        a=mat_[0,0]*mat[0,0] + mat_[0,1]*mat[1,0]
        b=mat_[0,0]*mat[0,1] + mat_[0,1]*mat[1,1]
        c=mat_[1,0]*mat[0,0] + mat_[1,1]*mat[1,0]
        d=mat_[1,0]*mat[0,1] + mat_[1,1]*mat[1,1]
        mat__ = numpy.array([[a,b],
                             [c,d]])
        transfo[0]=(numpy.nonzero(mat__[0,:])[0][0]+1)*numpy.sign(mat__[0,numpy.nonzero(mat__[0,:])[0][0]])
        transfo[1]=(numpy.nonzero(mat__[1,:])[0][0]+1)*numpy.sign(mat__[1,numpy.nonzero(mat__[1,:])[0][0]])
        transfo[2]=3

        return transfo

    else: #3D

        i = ni//2; j = nj//2; k = nk//2
        ip1 = max(i+1,ni-1); jp1 = max(j+1,nj-1); kp1 = max(k+1,nk-1)
        ind = i + j*ni + k*ni*nj
        P0 = [ a[1][0][ind], a[1][1][ind],a[1][2][ind] ]
        ind = ip1 + j*ni + k*ni*nj
        P1 = [ a[1][0][ind], a[1][1][ind],a[1][2][ind] ]
        ind = i + jp1*ni + k*ni*nj
        P2 = [ a[1][0][ind], a[1][1][ind],a[1][2][ind] ]
        ind = i + j*ni + kp1*ni*nj
        P3 = [ a[1][0][ind], a[1][1][ind],a[1][2][ind] ]
        l1 = Vector.sub(P1,P0)
        l2 = Vector.sub(P2,P0)
        l3 = Vector.sub(P3,P0)
        l1=Vector.normalize(l1)
        l2=Vector.normalize(l2)
        l3=Vector.normalize(l3)
        x=[1.0,0.0,0.0]
        y=[0.0,1.0,0.0]
        z=[0.0,0.0,1.0]
        a=Vector.dot(l1,x);b=Vector.dot(l1,y);c=Vector.dot(l1,z)
        d=Vector.dot(l2,x);e=Vector.dot(l2,y);f=Vector.dot(l2,z)
        g=Vector.dot(l3,x);h=Vector.dot(l3,y);i=Vector.dot(l3,z)
        mat_ = numpy.array([[a,b,c],
                            [d,e,f],
                            [g,h,i]])

        a = C.getFields(Internal.__GridCoordinates__, zrcv, api=1)[0]
        ni = a[2]; nj=a[3]; nk=a[4]
        i = ni//2; j = nj//2; k = nk//2
        ip1 = max(i+1,ni-1); jp1 = max(j+1,nj-1); kp1 = max(k+1,nk-1)
        ind = i + j*ni + k*ni*nj
        P0_  = [ a[1][0][ind], a[1][1][ind],a[1][2][ind] ]
        ind = ip1 + j*ni + k*ni*nj
        P1_ = [ a[1][0][ind], a[1][1][ind],a[1][2][ind] ]
        ind = i + jp1*ni + k*ni*nj
        P2_ = [ a[1][0][ind], a[1][1][ind],a[1][2][ind] ]
        ind = i + j*ni + kp1*ni*nj
        P3_ = [ a[1][0][ind], a[1][1][ind],a[1][2][ind] ]
        l1rcv = Vector.sub(P1_,P0_)
        l2rcv = Vector.sub(P2_,P0_)
        l3rcv = Vector.sub(P3_,P0_)
        l1rcv=Vector.normalize(l1rcv)
        l2rcv=Vector.normalize(l2rcv)
        l3rcv=Vector.normalize(l3rcv)
        x=[1.0,0.0,0.0]
        y=[0.0,1.0,0.0]
        z=[0.0,0.0,1.0]
        a=Vector.dot(l1rcv,x);b=Vector.dot(l1rcv,y);c=Vector.dot(l1rcv,z)
        d=Vector.dot(l2rcv,x);e=Vector.dot(l2rcv,y);f=Vector.dot(l2rcv,z)
        g=Vector.dot(l3rcv,x);h=Vector.dot(l3rcv,y);i=Vector.dot(l3rcv,z)
        mat = numpy.array([[a,b,c],
                           [d,e,f],
                           [g,h,i]])
        mat=mat.T
        a=mat_[0,0]*mat[0,0] + mat_[0,1]*mat[1,0] + mat_[0,2]*mat[2,0]
        b=mat_[0,0]*mat[0,1] + mat_[0,1]*mat[1,1] + mat_[0,2]*mat[2,1]
        c=mat_[0,0]*mat[0,2] + mat_[0,1]*mat[1,2] + mat_[0,2]*mat[2,2]
        d=mat_[1,0]*mat[0,0] + mat_[1,1]*mat[1,0] + mat_[1,2]*mat[2,0]
        e=mat_[1,0]*mat[0,1] + mat_[1,1]*mat[1,1] + mat_[1,2]*mat[2,1]
        f=mat_[1,0]*mat[0,2] + mat_[1,1]*mat[1,2] + mat_[1,2]*mat[2,2]
        g=mat_[2,0]*mat[0,0] + mat_[2,1]*mat[1,0] + mat_[2,2]*mat[2,0]
        h=mat_[2,0]*mat[0,1] + mat_[2,1]*mat[1,1] + mat_[2,2]*mat[2,1]
        i=mat_[2,0]*mat[0,2] + mat_[2,1]*mat[1,2] + mat_[2,2]*mat[2,2]

        mat__ = numpy.array([[a,b,c],
                             [d,e,f],
                             [g,h,i]])
        transfo[0]=(numpy.nonzero(mat__[0,:])[0][0]+1)*numpy.sign(mat__[0,numpy.nonzero(mat__[0,:])[0][0]])
        transfo[1]=(numpy.nonzero(mat__[1,:])[0][0]+1)*numpy.sign(mat__[1,numpy.nonzero(mat__[1,:])[0][0]])
        transfo[2]=(numpy.nonzero(mat__[2,:])[0][0]+1)*numpy.sign(mat__[2,numpy.nonzero(mat__[2,:])[0][0]])

        return transfo

#================================================================================
# Chimere: modif cellN BC pour interp chimere pour eviter interp inutile
# IN: t
# IN: depth= nombre de couche
#================================================================================
def _modifcellNBC(t, BCTarget=['BCWall'], depth=2,dim=3):

    for z in Internal.getZones(t):
        bc = Internal.getNodeFromName1(z,'ZoneBC')
        bcs=[]
        if bc is not None: bcs = Internal.getNodesFromType(bc,'BC_t')
        for bc in bcs:
            bctyp = Internal.getValue(bc)
            if bctyp in BCTarget:
                ptrg = Internal.getNodeFromName1(bc,'PointRange')
                dimZ = Internal.getZoneDim(z)
                #print(z[0], dimZ,'range', ptrg[1])

                idir = GhostCells.getDirection__(dim, [ptrg])

                #print('idir=',idir)
                cellN= Internal.getNodeFromName2(z,'cellN')[1]

                if idir==0:
                    jmin =  ptrg[1][1,0]-3; jmax =  ptrg[1][1,1]+1
                    kmin =  ptrg[1][2,0]-3; kmax =  ptrg[1][2,1]+1
                    for k in range(kmin,kmax):
                        for j in range(jmin,jmax):
                            for d in range(depth):
                                if cellN[depth-1-d,j,k]>=1.99 : cellN[depth-1-d,j,k]=1
                if idir==1:
                    i1 = ptrg[1][0,1]-depth-1
                    jmin =  ptrg[1][1,0]-3; jmax =  ptrg[1][1,1]+1
                    kmin =  ptrg[1][2,0]-3; kmax =  ptrg[1][2,1]+1
                    for k in range(kmin,kmax):
                        for j in range(jmin,jmax):
                            for d in range(depth):
                                if cellN[i1+d,j,k]>=1.99: cellN[i1+d,j,k]=1
                if idir==2:
                    imin =  ptrg[1][0,0]-3; imax =  ptrg[1][0,1]+1
                    kmin =  ptrg[1][2,0]-3; kmax =  ptrg[1][2,1]+1
                    #print(z[0],'range', imin, imax,kmin,kmax,depth-1)
                    for k in range(kmin,kmax):
                        for i in range(imin,imax):
                            for d in range(depth):
                                if  cellN[i,depth-1-d,k]>=1.99: cellN[i,depth-1-d,k]= 1
                if idir==3:
                    js = ptrg[1][1,1]-depth-1
                    imin =  ptrg[1][0,0]-3; imax =  ptrg[1][0,1]+1
                    kmin =  ptrg[1][2,0]-3; kmax =  ptrg[1][2,1]+1
                    for k in range(kmin,kmax):
                        for i in range(imin,imax):
                            for d in range(depth):
                                if cellN[i,js+d,k]>=1.99: cellN[i,js+d,k]= 1.
                if idir==4:
                    imin =  ptrg[1][0,0]-3; imax =  ptrg[1][0,1]+1
                    jmin =  ptrg[1][1,0]-3; jmax =  ptrg[1][1,1]+1
                    for j in range(jmin,jmax):
                        for i in range(imin,imax):
                            for d in range(depth):
                                if cellN[i,j,depth-1-d]>=1.99: cellN[i,j,depth-1-d]= 1.

                if idir==5:
                    ks = ptrg[1][2,1]-depth-1
                    imin =  ptrg[1][0,0]-3; imax =  ptrg[1][0,1]+1
                    jmin =  ptrg[1][1,0]-3; jmax =  ptrg[1][1,1]+1
                    for j in range(jmin,jmax):
                        for i in range(imin,imax):
                            for d in range(depth):
                                if cellN[i,j,ks+d]>=1.99: cellN[i,j,ks+d]= 1.
    return None
