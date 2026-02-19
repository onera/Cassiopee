"""Toolbox for IBM preprocessing
"""
from . import PyTree as X

try:
    import Converter.PyTree as C
    import Generator.PyTree as G
    import Transform.PyTree as T
    import Converter.Internal as Internal
    import Post.PyTree as P
    from .ToolboxIBM import *
except:
    raise ImportError("Connector.ToolboxIBM requires Converter, Generator, Transform, Dist2Walls and Post modules.")

varsn = ['gradxTurbulentDistance','gradyTurbulentDistance','gradzTurbulentDistance']
TOLDIST = 1.e-14
SHIFTF = 1.e-10
OPTFRONT = False
EPSCART = 1.e-6

#======================================================
# Returns NearBody meshes around each body component
# No BC is defined at external boundaries
#======================================================
def generateCompositeIBMMesh(tb, vmin, snears, dfar, dfarloc=0., DEPTH=2, NP=0, check=True, merged=1, sizeMax=4000000, symmetry=0, externalBCType='BCFarfield'):
    tbox = None; snearsf = None
    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    if dimPb is None: raise ValueError('EquationDimension is missing in input body tree.')
    dimPb = Internal.getValue(dimPb)
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    refstate = C.getState(tb)
    # list of near body meshes
    t = Internal.copyRef(tb)
    Internal._rmNodesFromType(t, "Zone_t")
    snearsExt=[]
    tblank = Internal.copyRef(tb)

    DEPTHEXT = 3*DEPTH+1
    tov = Internal.rmNodesFromType(tb,'Zone_t')
    for nob in range(len(tb[2])):
        base = tb[2][nob]
        if base[3] == 'CGNSBase_t':
            basename = base[0]
            res = generateIBMMesh_legacy(base, vmin, snears, dfarloc, DEPTH=DEPTH, NP=NP, tbox=None, snearsf=None,
                                         check=check, merged=merged, symmetry=symmetry, sizeMax=sizeMax, externalBCType='BCDummy', to=None,
                                         composite=1, mergeByParents=False)

            res = Internal.getZones(res)
            t[2][nob][2]=res
            # blanking box for off-body mesh
            extFacesL = C.extractBCOfType(res,"BCDummy")
            bbl = G.bbox(extFacesL)
            lMax = 0.
            for extf in extFacesL:
                if dimPb == 2: extf = T.subzone(extf,(1,1,1),(-1,1,1))
                extf = G.getVolumeMap(extf)
                dxloc=C.getValue(extf,'centers:vol',0)
                if dimPb == 3:
                    dxloc=dxloc**0.5
                lMax = max(lMax, dxloc)
                snearsExt.append(dxloc)
                tov[2][nob][2].append(extf)

            # ATTENTION : engendre une boite parallelepipedique - peut ne pas etre valide dans
            # les cas ou le maillage cartesien proche corps n est pas parallelepipedique
            # A ameliorer
            # IDEM pour les conditions aux limites dans octree2StructLoc et generateIBMMesh
            xminb = bbl[0]+DEPTHEXT*lMax; xmaxb = bbl[3]-DEPTHEXT*lMax
            yminb = bbl[1]+DEPTHEXT*lMax; ymaxb = bbl[4]-DEPTHEXT*lMax
            if dimPb == 3:
                zminb = bbl[2]+DEPTHEXT*lMax; zmaxb = bbl[5]-DEPTHEXT*lMax
                blankingBoxL = G.cart((xminb,yminb,zminb),(xmaxb-xminb,ymaxb-yminb,zmaxb-zminb),(2,2,2))
            else:
                blankingBoxL = G.cart((xminb,yminb,0.),(xmaxb-xminb,ymaxb-yminb,1.),(2,2,1))
            blankingBoxL = P.exteriorFaces(blankingBoxL)
            tblank[2][nob][2]=[blankingBoxL]

    tcart = generateIBMMesh_legacy(tov, vmin, snearsExt, dfar, DEPTH=DEPTH, NP=NP, tbox=tbox,
                                   snearsf=snearsf, check=check, merged=1, sizeMax=sizeMax,
                                   symmetry=symmetry, externalBCType='BCFarfield', to=None, mergeByParents=True)
    tcart[2][1][0] = 'OffBody'
    C._rmBCOfType(t,'BCDummy') # near body grids external borders must be BCOverlap
    t = C.fillEmptyBCWith(t,'ov_ext','BCOverlap',dim=dimPb)
    t = C.mergeTrees(t,tcart)
    model = Internal.getValue(model)
    C._addState(t,'GoverningEquations',model)
    C._addState(t,'EquationDimension',dimPb)
    C._addState(t,state=refstate)
    return t, tblank

#-----------------------------------------------------------
# Preprocessing pour les maillages proches corps en IBM
#-----------------------------------------------------------
def prepareCompositeIBMData(t,tb, DEPTH=2, loc='centers', frontType=1):
    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    if dimPb is None: raise ValueError('EquationDimension is missing in input body tree.')
    dimPb = Internal.getValue(dimPb)
    tc = C.newPyTree()
    for nob in range(len(t[2])):
        if Internal.getType(t[2][nob])=='CGNSBase_t':
            basename = t[2][nob][0]
            C._addBase2PyTree(tc,basename)
    for nob in range(len(tb[2])):
        if tb[2][nob][3]=='CGNSBase_t':
            basename = tb[2][nob][0]
            tloc = C.newPyTree([basename]); tloc[2][1]=t[2][nob]
            tbloc = C.newPyTree([basename]); tbloc[2][1]=tb[2][nob]
            # prepro IBM + interpolation entre blocs internes
            tloc,tcloc=prepareIBMData_legacy(tloc,tbloc,DEPTH=DEPTH, loc=loc, frontType=frontType, interp='composite')
            C._cpVars(tloc,'centers:cellN',tcloc,'cellN')
            tc[2][nob][2]+=Internal.getZones(tcloc)
            Internal._rmNodesFromType(t[2][nob],"Zone_t")
            t[2][nob][2]+=Internal.getZones(tloc)
    #
    offbodyZones = C.node2Center(t[2][-1][2])
    Internal._rmNodesByName(tc, "*TurbulentDistance*")
    Internal._rmNodesByName(tc, Internal.__FlowSolutionCenters__)
    tc[2][-1][2]+= offbodyZones
    return t, tc

def prepareCompositeChimeraData(t,tc,tblank,noBaseOff, DEPTH=2,loc='centers',
                                NIT=1, RotationCenter=None, RotationAngle=None, Translation=None):
    tBB = G.BB(tc)
    listOfOffBodyIntersectingNBZones=getListOfOffBodyIntersectingNBZones(tBB, noBaseOff, NIT=NIT, DEPTH=DEPTH,
                                                                         RotationCenter=RotationCenter,
                                                                         RotationAngle=RotationAngle,
                                                                         Translation=Translation)
    listOfSteadyOffBodyZones=getListOfSteadyOffBodyZones(tBB,noBaseOff,listOfOffBodyIntersectingNBZones)

    t,tc=prepareSteadyOffBodyChimeraData(t,tc,tblank, noBaseOff, tBB=tBB,DEPTH=2,loc='centers', NIT=NIT,
                                         RotationCenter=RotationCenter, RotationAngle=RotationAngle, Translation=Translation,
                                         listOfSteadyOffBodyZones=listOfSteadyOffBodyZones)

    t,tc=prepareMotionChimeraData(t,tc,tblank,noBaseOff, tBB=tBB,DEPTH=2,loc='centers', NIT=NIT,
                                  RotationCenter=RotationCenter, RotationAngle=RotationAngle, Translation=Translation,
                                  listOfSteadyOffBodyZones=listOfSteadyOffBodyZones,
                                  listOfOffBodyIntersectingNBZones=listOfOffBodyIntersectingNBZones)

    Internal._rmNodesByName(tc,Internal.__FlowSolutionNodes__)
    Internal._rmNodesByName(tc,Internal.__GridCoordinates__)
    return t,tc

#-------------------------------------------------------------------------------------------------------------------
# in : rotation angle en degres : a chaque pas de temps (dtheta)
# Soit rotation avec angle constant = RotationAngle=[dalphaX,dalphaY,dalphaZ]
# Soit rotation avec angle alpha(t) : RotationAngle=[[alphaX(t0),alphaY(t0),alphaZ(t0)],...,[alphaX(tNIT),alphaY(tNIT),alphaZ(tNIT)]]
# in : translation : translation a chq pas de temps - peut varier au cours des iterations comme la rotation
# in : noBaseOff : numero de la base correspondant au maillage de fond
# in : tblank : surface de masquage pour masquer le maillage de fond
#-------------------------------------------------------------------------------------------------------------------
def getListOfSteadyOffBodyZones(t, noBaseOff, listOfOffBodyIntersectingNBZones):

    listOfSteadyOffBodyZones=[]
    for z in Internal.getZones(t[2][noBaseOff]):
        if z[0] not in listOfOffBodyIntersectingNBZones: listOfSteadyOffBodyZones.append(z[0])
    return listOfSteadyOffBodyZones

#-------------------------------------------------------------------------------------------------------------------
# in : rotation angle en degres : a chaque pas de temps (dtheta)
# Soit rotation avec angle constant = RotationAngle=[dalphaX,dalphaY,dalphaZ]
# Soit rotation avec angle alpha(t) : RotationAngle=[[alphaX(t0),alphaY(t0),alphaZ(t0)],...,[alphaX(tNIT),alphaY(tNIT),alphaZ(tNIT)]]
# in : translation : translation a chq pas de temps - peut varier au cours des iterations comme la rotation
# in : noBaseOff : numero de la base correspondant au maillage de fond
# in : tBB : arbre de BB (Axis Aligned)
#-------------------------------------------------------------------------------------------------------------------
def getListOfOffBodyIntersectingNBZones(tBB, noBaseOff, NIT=1, DEPTH=2,
                                        RotationCenter=None, RotationAngle=None, Translation=None):

    rotation = True; translation = True
    if RotationAngle is None: rotation = False
    if Translation is None: translation = False
    if rotation is True and translation is True:
        raise ValueError("ToolboxIBM: translation and rotation not yet possible together.")
    constantMotion = False
    if rotation:
        if len(RotationAngle)==3 and isinstance(RotationAngle[0],list)==False:
            constantMotion = True
    if translation:
        if len(Translation)==3 and isinstance(Translation[0],list)==False:
            constantMotion=True

    xc0=RotationCenter[0];yc0=RotationCenter[1]; zc0=RotationCenter[2]

    listOfOffBodyIntersectingNBZones=[]
    C._initVars(tBB[2][noBaseOff],"{CoordinateXInit}={CoordinateX}")
    C._initVars(tBB[2][noBaseOff],"{CoordinateYInit}={CoordinateY}")
    C._initVars(tBB[2][noBaseOff],"{CoordinateZInit}={CoordinateZ}")

    for it in range(NIT):
        if constantMotion:
            if rotation:
                angleX = RotationAngle[0]*it
                angleY = RotationAngle[1]*it
                angleZ = RotationAngle[2]*it
                print('Rotation (degres) : (alphaX=%g,alphaY=%g,alphaZ=%g)'%(angleX,angleY,angleZ))
            elif translation:
                tx = Translation[0]
                ty = Translation[1]
                tz = Translation[2]
        else:
            if rotation:
                angleX = RotationAngle[it][0]
                angleY = RotationAngle[it][1]
                angleZ = RotationAngle[it][2]
                print('Rotation (degres) : (alphaX=%g,alphaY=%g,alphaZ=%g)'%(angleX,angleY,angleZ))
            elif translation:
                tx = Translation[it][0]
                ty = Translation[it][1]
                tz = Translation[it][2]

        # on fait bouger le maillage de fond dans le mvt oppose - pour ne pas bouger ts les maillages proches corps
        if rotation:
            tBB[2][noBaseOff]=T.rotate(tBB[2][noBaseOff],(xc0,yc0,zc0),(-angleX,-angleY,-angleZ))
        elif translation:
            tBB[2][noBaseOff]=T.translate(tBB[2][noBaseOff],(-tx,-ty,-tz))

        for nob in range(len(tBB[2])):
            if nob != noBaseOff and Internal.getType(tBB[2][nob])=='CGNSBase_t':
                dictOfIntersectingZones = X.getIntersectingDomains(tBB[2][nob],t2=tBB[2][noBaseOff],
                                                                   method='AABB',taabb=tBB[2][nob],taabb2=tBB[2][noBaseOff])
                for bodyZone in dictOfIntersectingZones:
                    listOfOffBodyIntersectingNBZones+=dictOfIntersectingZones[bodyZone]
                listOfOffBodyIntersectingNBZones=list(set(listOfOffBodyIntersectingNBZones))

        C._initVars(tBB[2][noBaseOff],"{CoordinateX}={CoordinateXInit}")
        C._initVars(tBB[2][noBaseOff],"{CoordinateY}={CoordinateYInit}")
        C._initVars(tBB[2][noBaseOff],"{CoordinateZ}={CoordinateZInit}")

    listOfOffBodyIntersectingNBZones=list(set(listOfOffBodyIntersectingNBZones))
    return listOfOffBodyIntersectingNBZones

#-------------------------------------------------------------------------------------------------------------------
# Calcul des donnees d interpolation pour les zones proches corps et de fond concernees par le mouvement
#-------------------------------------------------------------------------------------------------------------------
def prepareMotionChimeraData(t,tc,tblank,noBaseOff, tBB=None,DEPTH=2,loc='centers',
                             NIT=1, RotationCenter=None, RotationAngle=None, Translation=None,
                             listOfSteadyOffBodyZones=None,listOfOffBodyIntersectingNBZones=None):

    # arbre de BBox des zones donneuses
    if tBB is None: tBB = G.BB(tc)# BB des zones donneuses

    if listOfOffBodyIntersectingNBZones is None:
        listOfOffBodyIntersectingNBZones=getListOfOffBodyIntersectingNBZones(tBB, noBaseOff, NIT=NIT, DEPTH=DEPTH,
                                                                             RotationCenter=RotationCenter,
                                                                             RotationAngle=RotationAngle,
                                                                             Translation=Translation)
        listOfSteadyOffBodyZones=getListOfSteadyOffBodyZones(tBB,noBaseOff,listOfOffBodyIntersectingNBZones)

    if listOfSteadyOffBodyZones is None:
        listOfSteadyOffBodyZones=getListOfSteadyOffBodyZones(tBB,noBaseOff,listOfOffBodyIntersectingNBZones)

    dimPb = Internal.getNodeFromName(t, 'EquationDimension')
    if dimPb is None: raise ValueError('EquationDimension is missing in input body tree.')
    dimPb = Internal.getValue(dimPb)
    if dimPb == 2:
        z0 = Internal.getNodeFromType2(t, 'Zone_t')
        dims = Internal.getZoneDim(z0)
        npts = dims[1]*dims[2]*dims[3]
        zmin = C.getValue(z0,'CoordinateZ',0)
        zmax = C.getValue(z0,'CoordinateZ',npts-1)
        dz = zmax-zmin
        # Creation du corps 2D pour le preprocessing IBC
        tblank = T.addkplane(tblank)
        tblank = T.contract(tblank, (0,0,0), (1,0,0), (0,1,0), dz)

    rotation = True; translation = True
    if RotationAngle is None: rotation = False
    if Translation is None: translation = False
    if rotation is True and translation is True:
        raise ValueError("ToolboxIBM: translation and rotation not yet possible together.")
    constantMotion = False
    if rotation:
        if len(RotationAngle)==3 and isinstance(RotationAngle[0],list)==False:
            constantMotion = True
    if translation:
        if len(Translation)==3 and isinstance(Translation[0],list)==False:
            constantMotion=True

    xc0=RotationCenter[0];yc0=RotationCenter[1]; zc0=RotationCenter[2]

    C._initVars(t[2][noBaseOff],'{centers:cellNInit}={centers:cellN}')
    C._initVars(tc[2][noBaseOff],"{cellNInit}={cellN}")
    dictOfNearBodyBaseNb={}
    dictOfNearBodyZoneNb={}
    for nob in range(len(tc[2])):
        if nob != noBaseOff:
            base = tc[2][nob]
            if Internal.getType(base)=='CGNSBase_t':
                for noz in range(len(base[2])):
                    zone = base[2][noz]
                    if Internal.getType(zone)=='Zone_t':
                        zname=zone[0]
                        dictOfNearBodyZoneNb[zname]=noz
                        dictOfNearBodyBaseNb[zname]=nob

    dictOfOffBodyZoneNb={}
    dictOfOffBodyADT={} # preconditionnement
    for noz in range(len(tc[2][noBaseOff][2])):
        z = tc[2][noBaseOff][2][noz]
        zname = z[0]
        if Internal.getType(z)=='Zone_t':
            dictOfOffBodyZoneNb[zname] = noz
            if zname in listOfOffBodyIntersectingNBZones:
                zc = tc[2][noBaseOff][2][noz]
                hook0=C.createHook(zc,'adt')
                dictOfOffBodyADT[zname]=hook0

    C._initVars(tc,'{CoordinateXInit}={CoordinateX}')
    C._initVars(tc,'{CoordinateYInit}={CoordinateY}')
    C._initVars(tc,'{CoordinateZInit}={CoordinateZ}')
    C._initVars(t,'{CoordinateXInit}={CoordinateX}')
    C._initVars(t,'{CoordinateYInit}={CoordinateY}')
    C._initVars(t,'{CoordinateZInit}={CoordinateZ}')
    C._initVars(tBB,'{CoordinateXInit}={CoordinateX}')
    C._initVars(tBB,'{CoordinateYInit}={CoordinateY}')
    C._initVars(tBB,'{CoordinateZInit}={CoordinateZ}')

    if NIT > 1:
        for zc in Internal.getZones(tc):
            IDatas = Internal.getNodesFromName(zc,"ID_*")
            for ID in IDatas:
                name=ID[0].split('_')
                ID[0]='IDSteady_%s'%(name[1])

    tloc = C.newPyTree(['OffMotion'])
    for zr in t[2][noBaseOff][2]:
        if Internal.getType(zr)=='Zone_t':
            znamer=zr[0]
            if znamer not in listOfSteadyOffBodyZones:
                # print " Zone de fond en mvt %s"%znamer
                tloc[2][1][2].append(zr)
            else: tloc[2][1][2].append(zr)
        tBBloc=G.BB(tloc)


    # a remonter dans l interface
    intersectionsDictOffOff = X.getIntersectingDomains(tBB[2][noBaseOff], method='AABB', taabb=tBB[2][noBaseOff])
    listOfIntersectionDictsNBNB={}
    for nob in range(len(t[2])):
        if nob != noBaseOff and Internal.getType(t[2][nob])=='CGNSBase_t':
            intersectionsDictNB= X.getIntersectingDomains(tBB[2][nob], method='AABB', taabb=tBB[2][nob])
            listOfIntersectionDictsNBNB[nob]=intersectionsDictNB

    for it in range(NIT):
        print (' ------------------- Iteration %d ----------------------- '%it)
        if constantMotion:
            if rotation:
                angleX = RotationAngle[0]*it
                angleY = RotationAngle[1]*it
                angleZ = RotationAngle[2]*it
                print('Rotation (degres) : (alphaX=%g,alphaY=%g,alphaZ=%g)'%(angleX,angleY,angleZ))

            elif translation:
                tx = Translation[0]
                ty = Translation[1]
                tz = Translation[2]
        else:
            if rotation:
                angleX = RotationAngle[it][0]
                angleY = RotationAngle[it][1]
                angleZ = RotationAngle[it][2]
                print('Rotation (degres) : (alphaX=%g,alphaY=%g,alphaZ=%g)'%(angleX,angleY,angleZ))

            elif translation:
                tx = Translation[it][0]
                ty = Translation[it][1]
                tz = Translation[it][2]

        if rotation: tblankM = T.rotate(tblank,(xc0,yc0,zc0),(angleX,angleY,angleZ))
        elif translation: tblankM = T.translate(tblank,(tx,ty,tz))

        C._initVars(tloc,"{centers:cellN}={centers:cellNInit}")
        tloc = blankByIBCBodies(tloc,tblankM,'centers',dim=dimPb,gridType='composite')
        tloc = X.setHoleInterpolatedPoints(tloc,depth=DEPTH,loc='centers')
        for zloc in Internal.getZones(tloc):
            zname=zloc[0]
            nozloc = dictOfOffBodyZoneNb[zname]
            C._cpVars(zloc,'centers:cellN',tc[2][noBaseOff][2][nozloc],"cellN")

        dictOfMotionADT={} # ADT des blocs en mvt  a detruire a chq pas de temps

        # bases proches corps interpolees par le maillage de fond fixe + ses voisins
        for nob in range(len(t[2])):
            base = t[2][nob]
            if nob != noBaseOff and Internal.getType(base)=='CGNSBase_t':
                if rotation:
                    T._rotate(base,(xc0,yc0,zc0),(angleX,angleY,angleZ))
                    T._rotate(tBB[2][nob],(xc0,yc0,zc0),(angleX,angleY,angleZ))
                    T._rotate(tc[2][nob],(xc0,yc0,zc0),(angleX,angleY,angleZ))
                elif translation:
                    T._translate(base,(tx,ty,tz))
                    T._translate(tBB[2][nob],(tx,ty,tz))
                    T._translate(tc[2][nob],(tx,ty,tz))

        tBBNB = Internal.rmNodesByName(tBB, tBB[2][noBaseOff][0])
        intersectionsDictOffNB = X.getIntersectingDomains(tBB[2][noBaseOff], t2=tBBNB, method='AABB',
                                                          taabb=tBB[2][noBaseOff],taabb2=tBBNB)
        for nob in range(len(t[2])):
            base = t[2][nob]
            if nob != noBaseOff and Internal.getType(base)=='CGNSBase_t':
                # test intersection entre maillage proche corps et maillage de fond
                intersectionsDictNBO = X.getIntersectingDomains(tBB[2][nob], t2=tBB[2][noBaseOff], method='AABB',
                                                                taabb=tBB[2][nob],taabb2=tBB[2][noBaseOff])

                print('Near-body base %s in motion'%(base[0]))
                intersectionsDictNBNB=listOfIntersectionDictsNBNB[nob]
                for zr in Internal.getZones(base):
                    if Internal.getNodesFromName(zr,"ov_ext*")!=[]:
                        znamer = zr[0]
                        donorZones = []; hooks = []
                        for znamed in intersectionsDictNBO[znamer]:
                            if znamed in dictOfOffBodyZoneNb:# donneur=fond
                                nozd = dictOfOffBodyZoneNb[znamed]
                                zd = tc[2][noBaseOff][2][nozd]
                                hooks.append(dictOfOffBodyADT[znamed])
                                donorZones.append(zd)
                        for znamed in intersectionsDictNBNB[znamer]:
                            nozd = dictOfNearBodyZoneNb[znamed]
                            nobd = dictOfNearBodyBaseNb[znamed]
                            zd = tc[2][nobd][2][nozd]
                            if znamed not in dictOfMotionADT:
                                hook0 = C.createHook(zd,'adt')
                                dictOfMotionADT[znamed]=hook0
                            hooks.append(dictOfMotionADT[znamed])
                            donorZones.append(zd)

                        donorZones = X.setInterpData(zr,donorZones,nature=1,penalty=1,loc='centers',storage='inverse',sameName=1,\
                                                     hook=hooks, itype='chimera')
                        for zd in donorZones:
                            znamed = zd[0]
                            if znamed in dictOfOffBodyZoneNb:
                                nozd = dictOfOffBodyZoneNb[znamed]
                                tc[2][noBaseOff][2][nozd] = zd
                            elif znamed in dictOfNearBodyZoneNb:
                                nozd = dictOfNearBodyZoneNb[znamed]
                                nobd = dictOfNearBodyBaseNb[znamed]
                                tc[2][nobd][2][nozd]=zd

        # base de fond : interpolee depuis ses zones voisines + zones proches corps
        # les donneurs proches corps mobiles sont deja deplaces
        # intersectionsDict = X.getIntersectingDomains(tBBloc,t2=tBB,method='AABB',taabb=tBBloc,taabb2=tBB)
        for zr in Internal.getZones(tloc):
            znamer = zr[0]
            if znamer not in listOfSteadyOffBodyZones:
                donorZones = []; hooks = []
                for znamed in intersectionsDictOffOff[znamer]:# donneur=maillage de fond
                    nozd = dictOfOffBodyZoneNb[znamed]
                    zd = tc[2][noBaseOff][2][nozd]
                    if znamed not in dictOfOffBodyADT:
                        hook0 = C.createHook(zd,'adt')
                        dictOfOffBodyADT[znamed]=hook0
                    hooks.append(dictOfOffBodyADT[znamed])
                    donorZones.append(zd)
                for znamed in intersectionsDictOffNB[znamer]:
                    nozd = dictOfNearBodyZoneNb[znamed]
                    nobd = dictOfNearBodyBaseNb[znamed]
                    zd = tc[2][nobd][2][nozd]
                    if znamed not in dictOfMotionADT:
                        hook0 = C.createHook(zd,'adt')
                        dictOfMotionADT[znamed]=hook0
                    hooks.append(dictOfMotionADT[znamed])
                    donorZones.append(zd)

                # print 'Off-body motion zone %s'%znamer
                donorZones = X.setInterpData(zr,donorZones,nature=1,penalty=1,loc='centers',storage='inverse',sameName=1,\
                                             hook=hooks, itype='chimera')
                for zd in donorZones:
                    znamed = zd[0]
                    if znamed in dictOfOffBodyZoneNb:
                        nozd = dictOfOffBodyZoneNb[znamed]
                        tc[2][noBaseOff][2][nozd] = zd
                    elif znamed in dictOfNearBodyZoneNb:
                        nozd = dictOfNearBodyZoneNb[znamed]
                        nobd = dictOfNearBodyBaseNb[znamed]
                        tc[2][nobd][2][nozd]=zd

        for dnrname in dictOfMotionADT: C.freeHook(dictOfMotionADT[dnrname])

        # Reinit
        if NIT == 1:
            for dnrname in dictOfOffBodyADT: C.freeHook(dictOfOffBodyADT[dnrname])
            C._rmVars(tc,["CoordinateX","CoordinateY","CoordinateZ","cellNInit","CoordinateXInit","CoordinateYInit","CoordinateZInit"])
            C._initVars(t,'{CoordinateX}={CoordinateXInit}')
            C._initVars(t,'{CoordinateY}={CoordinateYInit}')
            C._initVars(t,'{CoordinateZ}={CoordinateZInit}')
            C._rmVars(t,["centers:cellNInit","CoordinateXInit","CoordinateYInit","CoordinateZInit"])
            C._cpVars(tc,"cellN",t,'centers:cellN')
            return t,tc

        else:
            for zc in Internal.getZones(tc):
                IDatas = Internal.getNodesFromName(zc,"ID_*")
                for ID in IDatas:
                    name=ID[0].split('_')
                    ID[0]='ID#%d_%s'%(it,name[1])

            if it == 5:
                C.convertPyTree2File(t,"t5.cgns")
                C.convertPyTree2File(tc,"tc5.cgns")
            C._initVars(tc[2][noBaseOff],"{cellN}={cellNInit}")# reinit
            C._initVars(t[2][noBaseOff],"{centers:cellN}={centers:cellNInit}")# reinit
            C._initVars(tc,'{CoordinateX}={CoordinateXInit}')
            C._initVars(tc,'{CoordinateY}={CoordinateYInit}')
            C._initVars(tc,'{CoordinateZ}={CoordinateZInit}')
            C._initVars(t,'{CoordinateX}={CoordinateXInit}')
            C._initVars(t,'{CoordinateY}={CoordinateYInit}')
            C._initVars(t,'{CoordinateZ}={CoordinateZInit}')
            C._initVars(tBB,'{CoordinateX}={CoordinateXInit}')
            C._initVars(tBB,'{CoordinateY}={CoordinateYInit}')
            C._initVars(tBB,'{CoordinateZ}={CoordinateZInit}')

    for dnrname in dictOfOffBodyADT: C.freeHook(dictOfOffBodyADT[dnrname])

    C._rmVars(t,["centers:cellNInit","CoordinateXInit","CoordinateYInit","CoordinateZInit"])
    C._initVars(tc[2][noBaseOff],"{cellN}={cellNInit}")
    C._cpVars(tc,"cellN",t,"centers:cellN")
    C._rmVars(tc,["cellNInit","CoordinateXInit","CoordinateYInit","CoordinateZInit"])
    return t, tc

#-------------------------------------------------------------------------------------------------------------------
# Calcul des donnees d interpolation pour les zones de fond non concernees par le mvt
#-------------------------------------------------------------------------------------------------------------------
def prepareSteadyOffBodyChimeraData(t,tc,tblank,noBaseOff, tBB=None,DEPTH=2,loc='centers',
                                    NIT=1, RotationCenter=[0,0,0], RotationAngle=[0,0,0.], Translation=[0,0,0],
                                    listOfSteadyOffBodyZones=None):

    if listOfSteadyOffBodyZones is None:
        raise ValueError("prepareSteadyOffBodyChimeraData: listOfSteadyOffBodyZones is None.")
    # arbre de BBox des zones donneuses
    if tBB is None: tBB = G.BB(tc)# BB des zones donneuses
    intersectionDict = X.getIntersectingDomains(tBB[2][noBaseOff])

    dictOfOffBodyZoneNb={}
    for noz in range(len(tc[2][noBaseOff][2])):
        z = tc[2][noBaseOff][2][noz]
        if Internal.getType(z)=='Zone_t': dictOfOffBodyZoneNb[z[0]] = noz

    dictOfOffBodyZoneNbRcv={}# t et tc peuvent ne pas avoir la meme structure
    for noz in range(len(t[2][noBaseOff][2])):
        z = t[2][noBaseOff][2][noz]
        if Internal.getType(z)=='Zone_t': dictOfOffBodyZoneNbRcv[z[0]] = noz

    C._initVars(t[2][noBaseOff],"centers:cellN", 1.)
    t[2][noBaseOff] = X.applyBCOverlaps(t[2][noBaseOff],depth=DEPTH,loc=loc)
    C._cpVars(t[2][noBaseOff],"centers:cellN",tc[2][noBaseOff],"cellN")
    dictOfADT={} # preconditionnement
    for zname in listOfSteadyOffBodyZones:
        noz = dictOfOffBodyZoneNbRcv[zname]
        z = t[2][noBaseOff][2][noz]
        intersectingZones = intersectionDict[zname]
        donorZones = []; hooks = []
        for znamed in intersectingZones:
            nozd = dictOfOffBodyZoneNb[znamed]
            zd = tc[2][noBaseOff][2][nozd]
            if znamed not in dictOfADT:
                hook0 = C.createHook(zd,'adt')
                dictOfADT[znamed] = hook0
            hooks.append(dictOfADT[znamed])
            donorZones.append(zd)

        donorZones = X.setInterpData(z,donorZones,nature=1,penalty=1,loc='centers',storage='inverse',sameName=1,\
                                     hook=hooks, itype='chimera')
        for zd in donorZones:
            znamed = zd[0]
            nozd = dictOfOffBodyZoneNb[znamed]
            tc[2][noBaseOff][2][nozd] = zd

    for dnrname in dictOfADT: C.freeHook(dictOfADT[dnrname])
    return t, tc
