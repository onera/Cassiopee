"""Find grid connectivity.
"""
from . import Connector
from . import connector
import numpy

from .OversetDataDtlocal import setInterpData3
from .OversetData import setInterpTransfers, _setInterpTransfers, __setInterpTransfers, _addCellN__, setInterpData, _setInterpData, _setIBCData, setInterpTransfersD, \
    _setInterpTransfersD, _setIBCTransfersForPressureGradients, _setIBCTransfersDForPressureGradients, _setIBCTransfers4GradP, _setIBCTransfers4GradP2, \
    _setIBCTransfers4GradP3,  _setIBCTransfers4FULLTBLE, _setIBCTransfers4FULLTBLE2, _modifcellNBC, getIntersectingDomains, getCEBBTimeIntersectingDomains, \
    getCEBBIntersectingDomains, cellN2OversetHoles, extractChimeraInfo, getOversetInfo, setIBCData, transferFields, setInterpData2, _setInterpData2
from .OversetDataElsA import _chimeraInfo, setInterpolations, chimeraInfo, chimeraTransfer
from .compactTransfers import ___setInterpTransfers, miseAPlatDonorTree__

__version__ = Connector.__version__
try:
    import KCore.kcore as KCore
    import Converter.Converter as Converter
    import Converter.PyTree as C
    import Converter.Internal as Internal
except: raise ImportError("Connector.PyTree: requires Converter module.")

__DEG2RAD__ = Internal.__DEG2RAD__
__RAD2DEG__ = Internal.__RAD2DEG__

#==============================================================================
# connectMatch between NGON zones
#==============================================================================
def _connectMatchNGON__(a, tol, dim, glob, allExtFaces=None, allExtIndices=None, periodic=0,
                        rotationCenter=None, rotationAngle=None, Translation=None, signT=0, signR=0,
                        unitAngle=None):
    # -------------------------
    # Exterior faces + indices
    # -------------------------
    # Ne conserve que les zones NGON (zonesp) avec leur no (indirZones)
    zones = Internal.getZones(a)
    zonesp = []; indirZones = []
    noz = 0
    for z in zones:
        dimZ = Internal.getZoneDim(z)
        if dimZ[3] == 'NGON':
            zonesp.append(z)
            indirZones.append(noz)
        noz += 1
    if len(zonesp) == 0: return glob

    # Get undefined faces
    if allExtFaces is None or allExtIndices is None:
        # - Par exteriorFaces
        #allExtIndices = []
        #allExtFaces = P.exteriorFaces(zonesp, indices=allExtIndices)
        #C._initVars(allExtFaces,'{centers:tag1}=-1.') # defines the opposite window
        #C._initVars(allExtFaces,'{centers:tag2}=-1.') # defines the opp index in opp window
        #indirBlkOfWins = []
        #for i in range(len(zonesp)): indirBlkOfWins.append(i)
        # - Par getEmptyWindows
        allExtFaces, indirBlkOfWins, allExtIndices = getEmptyWindowsInfoNGON__(zonesp, dim)
    else:
        # Cas periodique: allExtFaces vient d'un exteriorFaces
        indirBlkOfWins = []
        for i in range(len(allExtFaces)): indirBlkOfWins.append(i)

    # identify matching exterior faces
    if allExtFaces != []:
        nzones = len(zonesp)
        tagsF = C.node2Center(allExtFaces)
        tagsF = C.getAllFields(tagsF, 'nodes', api=1)
        tagsF = Connector.identifyMatching(tagsF, tol) # modifie tag1, tag2
        infos = Connector.gatherMatchingNGon__(tagsF, allExtIndices)
        rcvZones = infos[0]
        for i in range(rcvZones.size):
            rcvZones[i] = indirBlkOfWins[rcvZones[i]]
        dnrZones = infos[1]
        for i in range(dnrZones.size):
            dnrZones[i] = indirBlkOfWins[dnrZones[i]]
        allListRcvFaces = infos[2]
        allListDnrFaces = infos[3]
        nmatch = rcvZones.shape[0]
        for nm in range(nmatch):
            noz1p = rcvZones[nm]; noz2p = dnrZones[nm]
            isok = 1
            if periodic == 1:
                if noz1p >= nzones and noz2p < nzones: noz1p = noz1p-nzones
                elif noz2p >= nzones and noz1p < nzones: noz2p = noz2p-nzones
                else: isok = 0

            if isok == 1:
                noz1 = indirZones[noz1p]; noz2 = indirZones[noz2p]
                z1OppName = zones[noz1][0]
                z2OppName = zones[noz2][0]
                name1 = 'match%d_%d'%(noz1+1,glob); glob += 1
                name2 = 'match%d_%d'%(noz2+1,glob); glob += 1
                faceListR = allListRcvFaces[nm]
                faceListD = allListDnrFaces[nm]
                if periodic == 1:
                    tsign = numpy.array([signT,signT,signT], dtype=numpy.float64)
                    rsign = numpy.array([signR,signR,signR], dtype=numpy.float64)
                    C._addBC2Zone(zones[noz1], name1, 'BCMatch', faceList=faceListR,\
                                  zoneDonor=z2OppName, faceListDonor=faceListD,tol=tol,\
                                  rotationCenter=rotationCenter, rotationAngle=rsign*rotationAngle,\
                                  translation=tsign*Translation, unitAngle=unitAngle)

                    tsign = numpy.array([-signT,-signT,-signT], dtype=numpy.float64)
                    rsign = numpy.array([-signR,-signR,-signR], dtype=numpy.float64)
                    C._addBC2Zone(zones[noz2], name2, 'BCMatch', faceList=faceListD,\
                                  zoneDonor=z1OppName, faceListDonor=faceListR, tol=tol,\
                                  rotationCenter=rotationCenter, rotationAngle=rsign*rotationAngle,\
                                  translation=tsign*Translation, unitAngle=unitAngle)
                else:
                    C._addBC2Zone(zones[noz1], name1, 'BCMatch', faceList=faceListR,\
                                  zoneDonor=z2OppName, faceListDonor=faceListD, tol=tol)
                    C._addBC2Zone(zones[noz2], name2, 'BCMatch', faceList=faceListD,\
                                  zoneDonor=z1OppName, faceListDonor=faceListR,tol=tol)

        C.setFields(tagsF, allExtFaces, 'centers')
    # Sortie
    return glob

#==============================================================================
# connectMatch between NGON zones
#==============================================================================
def _connectMatchHybrid__(a, tol, dim, glob):

    # Tri des zones
    zones = []; indirZones = []
    noz = 0
    # identifyMatching runs over structured then unstructure. Same order must be kept
    # that is why 2 lists are first built and merged
    nzonesS = 0; nzonesU = 0
    for z in Internal.getZones(a):
        dimZ = Internal.getZoneDim(z)
        if dimZ[0] == 'Structured':
            zones.append(z)
            indirZones.append(noz)
            noz += 1
            nzonesS += 1
    for z in Internal.getZones(a):
        dimZ = Internal.getZoneDim(z)
        if dimZ[3] == 'NGON':
            zones.append(z)
            indirZones.append(noz)
            noz += 1
            nzonesU += 1
    if len(zones) == 0: return glob
    if nzonesS == 0 or nzonesU == 0: return glob

    # Get undefined faces
    allExtIndices=[]; allExtFaces=[]
    #for z in zones:
    #    extIndices=[]
    #    extFaces = P.exteriorFaces(z, indices=extIndices)
    #    Internal._rmNodesByType(extFaces,'ZoneBC_t')
    #    Internal._rmNodesByType(extFaces,'ZoneGridConnectivity_t')
    #    if Internal.getZoneType(z) == 1:
    #        extFaces = C.convertArray2NGon(extFaces)
    #        C._initVars(extFaces,'centers:tag1',-2.)
    #    else: C._initVars(extFaces,'centers:tag1',-1.)
    #    C._initVars(extFaces,'centers:tag2',-1.) # defines the opp index in opp window
    #    allExtFaces.append(extFaces); allExtIndices += extIndices
    allExtFaces, indirBlkOfWins, allExtIndices = getEmptyWindowsInfoHybrid__(zones, dim)

    # identify matching exterior faces
    if allExtFaces != []:
        tagsF = C.node2Center(allExtFaces)
        tagsF = C.getAllFields(tagsF, 'nodes', api=1)
        tagsF = Connector.identifyMatching(tagsF, tol)
        infos = Connector.gatherMatchingNGon__(tagsF, allExtIndices)
        rcvZones = infos[0]
        for i in range(rcvZones.size):
            rcvZones[i] = indirBlkOfWins[rcvZones[i]]
        dnrZones = infos[1]
        for i in range(dnrZones.size):
            dnrZones[i] = indirBlkOfWins[dnrZones[i]]
        allListRcvFaces = infos[2]
        allListDnrFaces = infos[3]
        nmatch = rcvZones.shape[0]
        for nm in range(nmatch):
            noz1p = rcvZones[nm]; noz2p = dnrZones[nm]
            noz1 = indirZones[noz1p]; noz2 = indirZones[noz2p]
            z1OppName = zones[noz1][0]
            z2OppName = zones[noz2][0]
            name1 = 'match%d_%d'%(noz1+1,glob); glob += 1
            name2 = 'match%d_%d'%(noz2+1,glob); glob += 1
            faceListR = allListRcvFaces[nm]
            faceListD = allListDnrFaces[nm]
            C._addBC2Zone(zones[noz1], name1, 'BCMatch', faceList=faceListR,\
                          zoneDonor=z2OppName, faceListDonor=faceListD, tol=tol)
            C._addBC2Zone(zones[noz2], name2, 'BCMatch', faceList=faceListD,\
                          zoneDonor=z1OppName, faceListDonor=faceListR, tol=tol)
        C.setFields(tagsF, allExtFaces, 'centers')
    # Sortie
    return glob

#==============================================================================
# Computes the 1 to 1 grid connectivity between structured zones
#==============================================================================
def _connectMatchStruct__(a, tol, dim, glob):
    zones = []
    dimPb = -1
    for z in Internal.getZones(a):
        if Internal.getZoneType(z) == 1:
            zones.append(z)
            if dimPb == -1: dimPb = Internal.getZoneDim(z)[4]
            else:
                if dimPb != Internal.getZoneDim(z)[4]:
                    print('Warning: some structured zones in connectMatch are not of same dimension. Function might fail...')
    # extract empty windows
    structTags,structWins,structIndirBlkOfWins,typeOfWins,dimsI,dimsJ,dimsK = \
        getEmptyWindowsInfoStruct__(zones, dim)

    model ='Euler'
    bases = Internal.getBases(a)
    if bases != []:
        c = Internal.getNodeFromName2(bases[0], 'GoverningEquations')
        if c is not None: model = Internal.getValue(c)

    # Identify matching cells for structured zones
    if structTags != []:
        structTags = Connector.identifyMatching(structTags, tol)
        structTags = Converter.extractVars(structTags, ['tag1','tag2'])
        # Gather into structured patches [[[noz1,noz2],[imin1,imax1,...],[imin2,imax2,...],trirac]]
        infos = Connector.gatherMatching(structWins, structTags, typeOfWins, structIndirBlkOfWins,
                                         dimsI, dimsJ, dimsK, dim, tol)
        for info in infos:
            noz1 = info[0][0]; noz2 = info[0][1]
            range1 = info[1]; range2 = info[2]
            topp0 = info[3]
            dimZ = Internal.getZoneDim(zones[noz1])
            dimzone = dimZ[4]
            if dimzone == 3: topp = [1,2,3]
            else:
                topp = [1,2]
                topp0 = [topp0[0], topp0[1]]

            if topp0[0] > 0: topp[topp0[0]-1] = 1
            else: topp[-topp0[0]-1] = -1
            if topp0[1] > 0: topp[topp0[1]-1] = 2
            else: topp[-topp0[1]-1] = -2
            if dimzone == 3:
                if topp0[2] > 0: topp[topp0[2]-1] = 3
                else: topp[-topp0[2]-1] = -3

            #print('match from ', zones[noz1][0], 'et ', zones[noz2][0])
            #print(topp, topp0, info[3])

            #------------------------------------------
            # addBC2Zone...
            name1 = 'match%d_%d'%(noz1+1,glob); glob += 1
            name2 = 'match%d_%d'%(noz2+1,glob); glob += 1
            C._addBC2Zone(zones[noz1], name1, 'BCMatch', range1, zones[noz2], range2, topp0)
            C._addBC2Zone(zones[noz2], name2, 'BCMatch', range2, zones[noz1], range1, topp)

            # couplage RANS/laminar ou euler
            model_z1 = model; model_z2 = model
            eq = Internal.getNodeFromName2(zones[noz1], 'GoverningEquations')
            if eq is not None: model_z1 = Internal.getValue(eq)
            eq = Internal.getNodeFromName2(zones[noz2], 'GoverningEquations')
            if eq is not None: model_z2 = Internal.getValue(eq)

            if model_z1 == 'NSTurbulent' and model_z1 != model_z2:
                # creation flag pour tranfert rans/LES
                datap1 = numpy.ones(1, dtype=Internal.E_NpyInt)
                datap2 = numpy.ones(1, dtype=Internal.E_NpyInt)
                Internal.createUniqueChild(Internal.getNodeFromName2(zones[noz1], name1), 'RANSLES', 'DataArray_t', datap1)
                Internal.createUniqueChild(Internal.getNodeFromName2(zones[noz2], name2), 'RANSLES', 'DataArray_t', datap2)
                name_extrap = 'RANS_LES%d_%d'%(noz1,noz2)
                C._addBC2Zone(zones[noz1],name_extrap,'BCExtrapolateRANS',range1)

            if model_z2 == 'NSTurbulent' and model_z1 != model_z2:
                datap1 = numpy.ones(1, dtype=Internal.E_NpyInt)
                datap2 = numpy.ones(1, dtype=Internal.E_NpyInt)
                Internal.createUniqueChild(Internal.getNodeFromName2(zones[noz2], name2), 'RANSLES', 'DataArray_t', datap2)
                Internal.createUniqueChild(Internal.getNodeFromName2(zones[noz1], name1), 'RANSLES', 'DataArray_t', datap1)
                name_extrap = 'RANS_LES%d_%d'%(noz2,noz1)
                C._addBC2Zone(zones[noz2], name_extrap, 'BCExtrapolateRANS', range2)

    return glob

#==============================================================================
# Returns the list of undefined windows as 2D zones
#         the tag defining matching blks and indices
#         the indirection tab to get the blk from which a window comes
#         the type of windows (i=1: 1, i=imax: 2, j=1: 3 etc)
#         the dimensions of the original blocks
#==============================================================================
def getEmptyWindowsInfoNGON__(t, dim=3):
    try: import Transform.PyTree as T
    except: raise ImportError("Connector.PyTree.getEmptyWindowsInfo__ requires Transform.PyTree module.")
    zones = Internal.getZones(t)
    nzones = len(zones)
    indirBlkOfWins=[]; allTags=[]; allExtIndices=[]
    for noz in range(nzones):
        z = zones[noz]
        zp = C.extractVars(z,['CoordinateX','CoordinateY','CoordinateZ']) # pour eviter des subzones de champs inutiles
        dimZ = Internal.getZoneDim(zp)
        if dimZ[3] == 'NGON':
            faceList = C.getEmptyBC(zp, dim=dim)
            for fl in faceList:
                winp = T.subzone(zp, fl, type='faces')
                C._initVars(winp,'centers:tag1',-1.) # defines the opposite window
                C._initVars(winp,'centers:tag2',-2.) # defines the opp index in opp window
                allTags += [winp]; indirBlkOfWins += [noz]
                allExtIndices += [fl]
    return allTags, indirBlkOfWins, allExtIndices

def getEmptyWindowsInfoHybrid__(t, dim=3):
    try: import Transform.PyTree as T
    except: raise ImportError("Connector.PyTree.getEmptyWindowsInfo__ requires Transform.PyTree module.")
    zones = Internal.getZones(t)
    nzones = len(zones)
    indirBlkOfWins=[]; allTags=[]; allExtIndices=[]
    for noz in range(nzones):
        z = zones[noz]
        zp = C.extractVars(z,['CoordinateX','CoordinateY','CoordinateZ']) # pour eviter des subzones de champs inutiles
        dimZ = Internal.getZoneDim(zp)
        if dimZ[3] == 'NGON':
            faceList = C.getEmptyBC(zp, dim=dim)
            for fl in faceList:
                winp = T.subzone(zp, fl, type='faces')
                C._initVars(winp,'centers:tag1',-1.) # defines the opposite window
                C._initVars(winp,'centers:tag2',-2.) # defines the opp index in opp window
                allTags += [winp]; indirBlkOfWins += [noz]
                allExtIndices += [fl]
        elif dimZ[0] == 'Structured':
            ranges = C.getEmptyBC(z, dim=dim)
            for r in ranges:
                winp = T.subzone(zp, (r[0],r[2],r[4]), (r[1],r[3],r[5]))
                winp = C.convertArray2NGon(winp)
                C._initVars(winp,'centers:tag1',-1.) # defines the opposite window
                C._initVars(winp,'centers:tag2',-2.) # defines the opp index in opp window
                allTags += [winp]; indirBlkOfWins += [noz]
                ind = Converter.converter.range2PointList(r[0],r[1],r[2],r[3],r[4],r[5],dimZ[1],dimZ[2],dimZ[3])
                allExtIndices += [ind]
    return allTags, indirBlkOfWins, allExtIndices

#==============================================================================
def getEmptyWindowsInfoStruct__(t, dim=3):
    try:import Transform.PyTree as T
    except:raise ImportError("getEmptyWindowsInfo__ requires Transform.PyTree modules.")
    zones = Internal.getZones(t)
    nzones = len(zones)
    allWins=[]; typeOfWins=[]; indirBlkOfWins=[]; allTags=[]
    dimsI=[]; dimsJ=[]; dimsK=[]
    for noz in range(nzones):
        z = zones[noz]
        dimsZ = Internal.getZoneDim(z)
        if dimsZ[0] == 'Structured':
            ni = dimsZ[1]; nj = dimsZ[2]; nk = dimsZ[3]
            nic = max(1,ni-1); njc = max(1,nj-1); nkc = max(1,nk-1)
            dimsI.append(ni); dimsJ.append(nj); dimsK.append(nk)
            ranges = C.getEmptyBC(z, dim=dim)
            if ranges != []:
                locWins=[]; locTypes=[]; locIndir=[]
                winp = T.subzone(z,(1,1,1),(1,nj,nk))
                win = C.getFields(Internal.__GridCoordinates__, winp, api=1)[0]
                locWins.append(win); locTypes.append(1); locIndir.append(noz)

                winp = T.subzone(z,(ni,1,1),(ni,nj,nk))
                win = C.getFields(Internal.__GridCoordinates__, winp, api=1)[0]
                locWins.append(win); locTypes.append(2); locIndir.append(noz)

                winp = T.subzone(z,(1,1,1),(ni,1,nk))
                win = C.getFields(Internal.__GridCoordinates__, winp, api=1)[0]
                locWins.append(win); locTypes.append(3); locIndir.append(noz)

                winp = T.subzone(z,(1,nj,1),(ni,nj,nk))
                win = C.getFields(Internal.__GridCoordinates__, winp, api=1)[0]
                locWins.append(win); locTypes.append(4); locIndir.append(noz)

                if dim == 3:
                    winp = T.subzone(z,(1,1,1),(ni,nj,1))
                    win = C.getFields(Internal.__GridCoordinates__, winp, api=1)[0]
                    locWins.append(win); locTypes.append(5); locIndir.append(noz)

                    winp = T.subzone(z,(1,1,nk),(ni,nj,nk))
                    win = C.getFields(Internal.__GridCoordinates__, winp, api=1)[0]
                    locWins.append(win); locTypes.append(6); locIndir.append(noz)

                locTags = Converter.node2Center(locWins)
                locTags = Converter.initVars(locTags,'tag1',-2.) # defines the opposite window
                locTags = Converter.initVars(locTags,'tag2',-2.) # defines the opposite index in opposite window
                for r in ranges:
                    imin = r[0]; jmin = r[2]; kmin = r[4]
                    imax = r[1]; jmax = r[3]; kmax = r[5]
                    now = 0
                    if imin == imax:
                        if imin == 1: now = 1
                        else: now = 2
                    elif jmin == jmax:
                        if jmin == 1: now = 3
                        else: now = 4
                    elif kmin == kmax:
                        if kmin == 1: now = 5
                        else: now = 6

                    tag = locTags[now-1]
                    postag = KCore.isNamePresent(tag, 'tag1')
                    taga = tag[1][postag,:]
                    imax = min(imax-1,nic)
                    jmax = min(jmax-1,njc)
                    kmax = min(kmax-1,nkc)
                    if dim == 2: kmax = max(1,kmax)

                    if now == 1 or now == 2:
                        for k in range(kmin-1,kmax):
                            taga[jmin-1+k*njc:jmax+k*njc] = -1.
                    elif now == 3 or now == 4:
                        for k in range(kmin-1,kmax):
                            taga[imin-1+k*nic:imax+k*nic] = -1.
                    elif now == 5 or now == 6:
                        for j in range(jmin-1,jmax):
                            taga[imin-1+j*nic:imax+j*nic] = -1.
                allWins += locWins
                allTags += locTags
                indirBlkOfWins += locIndir
                typeOfWins += locTypes
    return allTags, allWins, indirBlkOfWins, typeOfWins, dimsI, dimsJ, dimsK

#=============================================================================
# Computes the 1-to-1 connectivity between:
# 1. structured zones
# 2. NGON zones
# 3. structured/NGON
#=============================================================================
def connectMatch(t, tol=1.e-6, dim=3, type='all'):
    a,typen = Internal.node2PyTree(t)
    glob = 0
    if type == 'structured':
        glob = _connectMatchStruct__(a, tol, dim, glob)
    elif type == 'unstructured':
        glob = _connectMatchNGON__(a, tol, dim, glob)
    elif type == 'hybrid':
        glob = _connectMatchHybrid__(a, tol, dim, glob)
    else: # all
        glob = _connectMatchStruct__(a, tol, dim, glob)
        glob = _connectMatchNGON__(a, tol, dim, glob)
        glob = _connectMatchHybrid__(a, tol, dim, glob)
    return Internal.pyTree2Node(a, typen)

#===============================================================================
# Detection auto des raccords coincidents dans le cas periodique en
# rotation ou translation
# Standard: l'angle de rotation ou le vecteur de translation est oriente de
# l'interface courante a l'interface connectee.
# Dans le standard, si 3 angles de rotations non nuls, la composition se fait
# selon x, puis y, puis z
#===============================================================================
def connectMatchPeriodic(t, rotationCenter=[0.,0.,0.],
                         rotationAngle=[0.,0.,0.],
                         translation=[0.,0.,0.], tol=1.e-6, dim=3,
                         unitAngle=None):
    """Find periodic matching boundaries."""
    a = Internal.copyRef(t)
    removeFSN = []
    zones = Internal.getZones(a)
    for z in zones:
        if Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__) is None: removeFSN.append(z[0])

    a = C.tagDefinedBC(a)
    a = connectMatchPeriodicStruct__(a, rotationCenter, rotationAngle, translation, tol, dim, unitAngle)
    a = connectMatchPeriodicNGON__(a, rotationCenter, rotationAngle, translation, tol, dim, unitAngle)
    if len(zones) == len(removeFSN): C._rmNodes(a, Internal.__FlowSolutionNodes__)
    else: C._rmVars(a, 'definedBC')
    return a

def connectMatchPeriodicNGON__(a, rotationCenter, rotationAngle, translation, tol, dim, unitAngle):
    try: import Post.PyTree as P; import Transform.PyTree as T
    except: raise ImportError("connectMatchPeriodicNGON__: requires Transform and Post modules.")

    if unitAngle in ['Degree', None]: rotationAngleD=rotationAngle
    elif unitAngle == 'Radian': rotationAngleD=[v*__RAD2DEG__ for v in rotationAngle]
    else: raise ValueError('connectMatchPeriodicNGON__: value for unitAngle is not valid.')

    typeOfInputNode = Internal.typeOfNode(a)
    if typeOfInputNode == -1: raise ValueError("connectMatchPeriodicNGON__: invalid input node.")
    zones = Internal.getZones(a)
    indirZones = []; zonesU = []
    noz = 0
    for z in zones:
        dimZ = Internal.getZoneDim(z)
        if dimZ[3]=='NGON': zonesU.append(z); indirZones.append(noz)
        noz += 1
    if len(zonesU) == 0: return a

    # get exterior faces
    allExtIndices = []; allExtFaces0 = []
    for z in zonesU:
        indicesF = []
        f = P.exteriorFaces(z, indices=indicesF)
        indicesF = indicesF[0]
        C._initVars(f, 'centers:__tag__', 1.)

        # BC classique
        bnds = Internal.getNodesFromType2(z, 'BC_t')
        # BC Match
        bnds += Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
        # BC Overlap/NearMatch/NoMatch
        bnds += Internal.getNodesFromType2(z, 'GridConnectivity_t')

        zp = Internal.copyRef(z)
        C._deleteZoneBC__(zp); C._deleteFlowSolutions__(zp)

        defined = [] # BC deja definies
        for bc in bnds:
            flist = Internal.getNodeFromName1(bc, Internal.__FACELIST__)
            if flist is not None: defined.append(T.subzone(zp, flist[1], type='faces'))
            erange = Internal.getNodeFromName1(bc, Internal.__ELEMENTRANGE__)
            if erange is not None:
                r = erange[1]
                defined.append(C.selectOneConnectivity(zp,range=[r[0,0],r[0,1]]))

        hook = C.createHook(f, 'elementCenters')
        if defined != []:
            tag = Internal.getNodeFromName2(f, '__tag__')[1]
            defined = C.convertArray2NGon(defined)
            defined = T.join(defined)
            id0 = C.identifyElements(hook, defined)
            tag[id0[:]-1] = 0
        sel = P.selectCells2(f, 'centers:__tag__')
        id1 = C.identifyElements(hook,sel)
        id2 = numpy.empty(id1.size, dtype=Internal.E_NpyInt)
        id2[:] = indicesF[id1[:]-1]
        C.freeHook(hook)
        if id2.size != 0:
            allExtIndices.append(id2)
            C._initVars(sel,'centers:tag1',-1.) # defines the opposite window
            C._initVars(sel,'centers:tag2',-1.) # defines the opp index in opp window
            allExtFaces0.append(sel)
        else:
            allExtIndices.append(indicesF)
            C._rmVars(f, ['centers:__tag__'])
            C._initVars(f,'centers:tag1',-2.)# defines the opposite window
            C._initVars(f,'centers:tag2',-2.) # defines the opp index in opp window
            allExtFaces0.append(f)

    # duplicate exterior faces
    infoPer = duplicatePeriodicZones__(allExtFaces0,rotationCenter,rotationAngleD,translation,tol,dim)
    nzonesU = len(zonesU)
    typePeriodic = infoPer[0]
    if typePeriodic==1: signT = [-1,1]; signR = [0,0]
    elif typePeriodic==2: signT = [0,0]; signR = [-1,1]
    elif typePeriodic==3: signT = [-1,-1,1,1]; signR = [-1,1,-1,1]
    dupname = 'DUPPER_' # prefix for duplicated zones
    for i in range(1, len(infoPer)):
        # renommage des zones dupliquees
        for noz in range(nzonesU):
            zname = infoPer[i][noz][0]
            infoPer[i][noz][0] = dupname+zname
        glob = 0
        glob = _connectMatchNGON__(zonesU,tol,dim,glob,allExtFaces=allExtFaces0+infoPer[i],\
                                   allExtIndices=allExtIndices+allExtIndices, periodic=1,\
                                   rotationCenter=rotationCenter,rotationAngle=rotationAngle,\
                                   Translation=translation, signT=signT[i-1],signR=signR[i-1], \
                                   unitAngle=unitAngle)
        # update centers:tag1
        if glob > 0:
            tag1p = C.getField('centers:tag1', infoPer[i], api=1)
            tag1o = C.getField('centers:tag1', allExtFaces0, api=1)
            for noz in range(len(tag1p)): tag1p[noz][0]='tag1p'
            tag1o = Converter.addVars([tag1o,tag1p])
            tag1o = Converter.initVars(tag1o,'{tag1}=maximum({tag1},{tag1p})')
            C.setFields(tag1o, allExtFaces0, loc='centers')
            C._rmVars(allExtFaces0, 'centers:tag1p')

        infoPer[i] = []
    # update the original tree
    if typeOfInputNode == 1: # une zone
        a = zonesU[0]
    elif typeOfInputNode == 2: # une liste de zones
        noi = 0
        for indirZ in indirZones:
            zones[indirZ] = zonesU[noi]
            noi += 1
        return zones
    else: # base, liste de bases, arbre
        noi = 0
        for indirZ in indirZones:
            parent,d2 = Internal.getParentOfNode(a,zones[indirZ])
            parent[2][d2] = zonesU[noi]
            noi += 1
    return a

def connectMatchPeriodicStruct__(a,rotationCenter,rotationAngle,translation,tol,dim,unitAngle):
    typeOfInputNode = Internal.typeOfNode(a)
    if typeOfInputNode == -1: raise ValueError("connectMatchPeriodicStruct__: invalid input node.")
    zones = Internal.getZones(a)
    indirZones = []; zonesS = []
    noz = 0
    for z in zones:
        dimZ = Internal.getZoneDim(z)
        if dimZ[0] == 'Structured': zonesS.append(z); indirZones.append(noz)
        noz += 1
    if len(zonesS) == 0: return a

    if unitAngle in ['Degree',None]: rotationAngleD=rotationAngle
    elif unitAngle == 'Radian': rotationAngleD=[v*__RAD2DEG__ for v in rotationAngle]
    else: raise ValueError('connectMatchPeriodicStruct__: invalid value for unitAngle.')
    infoPer = duplicatePeriodicZones__(zonesS, rotationCenter, rotationAngleD, translation, tol, dim)
    nzonesS = len(zonesS)
    typePeriodic = infoPer[0]
    if typePeriodic == 1: signT = [-1,1]; signR = [0,0]
    elif typePeriodic == 2: signT = [0,0]; signR = [-1,1]
    elif typePeriodic == 3: signT = [-1,-1,1,1]; signR = [-1,1,-1,1]
    dupname = 'DUPPER_' # prefix for duplicated zones
    for i in range(1, len(infoPer)):
        # renommage des zones dupliquees
        for noz in range(nzonesS):
            zname = infoPer[i][noz][0]
            infoPer[i][noz][0] = dupname+zname
        glob = 0
        glob = _connectMatchStruct__(zonesS+infoPer[i], tol, dim, glob)
        for noz in range(nzonesS):
            gcnodes = Internal.getNodesFromType2(zonesS[noz], 'GridConnectivity1to1_t')
            _addPeriodicInfo__(gcnodes, rotationCenter, rotationAngle, translation,\
                               signT[i-1], signR[i-1], dupname, unitAngle=unitAngle)
            # CORRECTIF DANS LE CAS 180 DEG : LA PREMIERE PASSE FAIT LES 2 RACCORDS MIN ET MAX
            # MAIS L ANGLE N EST PAS BON POUR LE 2E RACCORD SI PAS CE CORRECTIF
            if i == 1:
                for angle in rotationAngle:
                    if angle == 180. or angle == -180.:
                        nogci = 0
                        for gc in gcnodes:
                            if Internal.getValue(gc) == Internal.getName(zonesS[noz]):
                                rotation_angle = Internal.getNodeFromName(gc, "RotationAngle")
                                if rotation_angle:
                                    nogci += 1
                                    if nogci == 2: # changement de signe
                                        rotation_angle = Internal.getValue(rotation_angle)
                                        for j in range(3): rotation_angle[j]=-rotation_angle[j]

        infoPer[i] = []
    # update the original tree
    if typeOfInputNode == 1: # une zone
        a = zonesS[0]
    elif typeOfInputNode == 2: # une liste de zones
        noi = 0
        for indirZ in indirZones:
            zones[indirZ] = zonesS[noi]
            noi += 1
        return zones
    else: # base, liste de bases, arbre
        noi = 0
        for indirZ in indirZones:
            parent,d2 = Internal.getParentOfNode(a,zones[indirZ])
            parent[2][d2] = zonesS[noi]
            noi += 1
    return a

#===============================================================================
def _addPeriodicInfo__(gcnodes,rotationCenter,rotationAngle,translation,signT,signR, dupname='DUPPER_',unitAngle=None):
    for info in gcnodes:
        if len(info[1]) > 8:
            donorNamePref = info[1][0:7]; donorName = info[1][7:]
            donorNamePref = Internal.getValue([0,donorNamePref])
            if donorNamePref == dupname: # cas periodique
                info[1] = donorName
                rotationAngleS = [v*signR for v in rotationAngle]
                translationS = [v*signT for v in translation]
                C._addPeriodicInfoInGC__(info, rotationCenter, rotationAngleS, translationS, unitAngle=unitAngle)
    return None

#============================================================================================
# duplique les zones par periodicite. En rotation, on suppose ici que l'angle est en degres
#============================================================================================
def duplicatePeriodicZones__(t, rotationCenter=[0.,0.,0.], rotationAngle=[0.,0.,0.],
                             translation=[0.,0.,0.], tol=1.e-6, dim=3):
    try: import Transform.PyTree as T
    except: raise ImportError("connectMatchPeriodic: requires Transform module.")
    a = Internal.copyRef(t)

    typePeriodic = 0 # 1: translation, 2: rotation, 3: les deux
    for i in range(3):
        if float(translation[i]) != 0.: typePeriodic = 1; break
    for i in range(3):
        if float(rotationAngle[i]) != 0.: typePeriodic += 2; break

    # Periodicite par rotation ou translation separement
    zones = Internal.getZones(a)
    dupZones = None
    if typePeriodic == 0: return None
    elif typePeriodic == 1: # translation pure
        zonesdupP = T.translate(zones,(translation[0],translation[1],translation[2]))# periodicite en +transVect
        zonesdupM = T.translate(zones,(-translation[0],-translation[1],-translation[2]))# periodicite en -transVect
        dupZones = [typePeriodic,zonesdupP,zonesdupM]
    elif typePeriodic == 2: # rotation pure
        zonesdupP = T.rotate(zones,(rotationCenter[0],rotationCenter[1],rotationCenter[2]),\
                             (rotationAngle[0],rotationAngle[1],rotationAngle[2]))
        zonesdupM = T.rotate(zones,(rotationCenter[0],rotationCenter[1],rotationCenter[2]),\
                             (-rotationAngle[0],-rotationAngle[1],-rotationAngle[2]))
        if 180. in rotationAngle or -180. in rotationAngle:
            dupZones = [typePeriodic,zonesdupP] # only one duplication
        else: dupZones = [typePeriodic,zonesdupP,zonesdupM]

    elif typePeriodic == 3: # rotation + translation
        zonesdup0 = T.translate(zones,(translation[0],translation[1],translation[2]))# periodicite en +transVect

        # 1. Premier sens de translation
        # 1.1 dans le premier sens de rotation
        zonesdupP = T.rotate(zonesdup0,(rotationCenter[0],rotationCenter[1],rotationCenter[2]),(rotationAngle[0],rotationAngle[1],rotationAngle[2]))
        # 1.2. dans le sens oppose
        zonesdupM = T.rotate(zonesdup0,(rotationCenter[0],rotationCenter[1],rotationCenter[2]),(-rotationAngle[0],-rotationAngle[1],-rotationAngle[2]))
        dupZones = [typePeriodic,zonesdupP,zonesdupM]

        zonesdup0 = T.translate(zones,(-translation[0],-translation[1],-translation[2]))# periodicite en -transVect
        zonesdupP = T.rotate(zonesdup0,(rotationCenter[0],rotationCenter[1],rotationCenter[2]),(rotationAngle[0],rotationAngle[1],rotationAngle[2]))
        # 2.2. dans le sens oppose
        zonesdupM = T.rotate(zonesdup0,(rotationCenter[0],rotationCenter[1],rotationCenter[2]),(-rotationAngle[0],-rotationAngle[1],-rotationAngle[2]))
        dupZones += [zonesdupP,zonesdupM]
    return dupZones

#==============================================================================
# Set degenerated BC
#==============================================================================
def setDegeneratedBC(t, dim=3, tol=1.e-10):
    """Find degenerated boundaries (lines)."""
    try: import Generator
    except: raise ImportError("setDegeneratedBC: requires Generator module.")

    a = Internal.copyRef(t)
    # Is t a topTree, a list of bases etc ?
    # listzones = 0: basis or toptree
    # listzones = 1: list of zones
    listzones = 1 # Is a  a list of zones or a basis or a toptree
    toptree = Internal.isTopTree(a)
    if toptree: listzones = 0
    else:
        # a base ou non
        base = Internal.getBases(a)
        if base != []: listzones = 0
        else: # liste zones ou zone ?
            stdNode = Internal.isStdNode(a)
            if stdNode != 0: return a

    allTags,allWins,indirBlkOfWins,typeOfWins,dimsI,dimsJ,dimsK=getEmptyWindowsInfoStruct__(a,dim)
    if allTags != []:
        volumes = Generator.getVolumeMap(allWins)
        allTags = Converter.addVars([allTags,volumes])
        allTags = Connector.identifyDegenerated(allTags, tol)
        allTags = Converter.extractVars(allTags,['tag1','tag2'])
        # Gather degenerated cells into structured patches [[[noz1],[imin1,imax1,...]],... ]
        infos = Connector.gatherDegenerated(allTags,typeOfWins,indirBlkOfWins,dimsI, dimsJ, dimsK, dim)
    else: infos = []
    zones = Internal.getZones(a)
    glob = 0
    for info in infos:
        noz1 = info[0]; range1 = info[1]
        # Get Parent Nodes
        if listzones == 0:
            r1 = Internal.getParentOfNode(a, zones[noz1])
            parent1 = r1[0]; d1 = r1[1]
        # addBC
        name1 = 'degen%d_%d'%(noz1+1,glob); glob+=1
        if dim == 2: zones[noz1] = C.addBC2Zone(zones[noz1],name1,'BCDegeneratePoint',range1)
        else: zones[noz1] = C.addBC2Zone(zones[noz1],name1,'BCDegenerateLine',range1)

        if listzones == 0: parent1[2][d1] = zones[noz1]
    if listzones == 1: a = zones
    return a

#=============================================================================
# Computes the n to p grid connectivity between zones defined in t
#=============================================================================
def connectNearMatch(t, ratio=2, tol=1.e-6, dim=3):
    """Find boundaries that matches with a given ratio."""
    try: import Generator.PyTree as G; import Transform.PyTree as T
    except: raise ImportError("connectNearMatch: requires Generator and Transform modules.")
    a,typen = Internal.node2PyTree(t)

    allRatios = []
    if not isinstance(ratio, list):
        if dim == 3:
            allRatios.append([ratio,1,1])#ij
            allRatios.append([1,ratio,1])#ik
            allRatios.append([1,1,ratio])#jk
            allRatios.append([ratio,ratio,1])#ij
            allRatios.append([ratio,1,ratio])#ik
            allRatios.append([1,ratio,ratio])#jk
            allRatios.append([ratio,ratio,ratio])#ijk
        else:
            allRatios.append([ratio,1,1])#ij
            allRatios.append([1,ratio,1])#ik
            allRatios.append([ratio,ratio,1])#ij
    else: allRatios = [ratio]

    model ='Euler'
    bases = Internal.getNodesFromType1(t, 'CGNSBase_t')
    if bases != []:
        eq = Internal.getNodeFromName2(bases[0], 'GoverningEquations')
        if eq is not None: model = Internal.getValue(eq)

    glob = 0
    zones = Internal.getZones(a)
    nzones = len(zones)
    for ratios in allRatios:
        zones2 = []
        for noz1 in range(nzones):
            z1 = zones[noz1]
            # modification des GC pour garder les frontieres definies par une BC lors du passage en oneovern
            bcmatch = Internal.getNodesFromType2(z1, 'GridConnectivity1to1_t')
            for i in bcmatch:
                r = Internal.getNodeFromName1(i, 'PointRange')
                if r is not None:
                    ranger = r[1]
                    w = Internal.range2Window(ranger)
                    z1 = C.addBC2Zone(z1, 'dummy', 'BCWall', w)

            bcnearmatch = Internal.getNodesFromType2(z1, 'GridConnectivity_t')
            for i in bcnearmatch:
                type = Internal.getNodeFromName1(i, 'GridConnectivityType')
                if type is not None:
                    val = Internal.getValue(type)
                    if val == 'Abutting':
                        r = Internal.getNodeFromName1(i, 'PointRange')
                        ranger = r[1]
                        w = Internal.range2Window(ranger)
                        z1 = C.addBC2Zone(z1, 'dummy', 'BCWall', w)
            # oneovern
            zones2.append(T.oneovern(z1,(ratios[0],ratios[1],ratios[2])))

        # Abutting info for oneovern zones
        allTags2, allWins2, indirBlkOfWins2, typeOfWins2,dimsI2,dimsJ2,dimsK2=getEmptyWindowsInfoStruct__(zones2,dim)
        #
        allTags1, allWins1, indirBlkOfWins1, typeOfWins1,dimsI1,dimsJ1,dimsK1=getEmptyWindowsInfoStruct__(zones,dim)
        #
        ntags1 = len(allTags1)
        for now1 in range(ntags1): indirBlkOfWins1[now1] +=nzones
        #
        if len(allTags2) < len(allTags1):
            for iblk in range(len(indirBlkOfWins2)):
                if indirBlkOfWins2[iblk] != indirBlkOfWins1[iblk]-len(zones):
                    allTags1.pop(iblk)
                    allWins1.pop(iblk)
                    indirBlkOfWins1.pop(iblk)
                    typeOfWins1.pop(iblk)

        # Identify matching windows of zones with windows of zones2
        if allTags1 != [] and len(allTags1) == len(allTags2):
            allTags = Connector.identifyMatchingNM(allTags2, allTags1,tol)
        else: allTags = []
        allWins = allWins2+allWins1
        typeOfWins = typeOfWins2+typeOfWins1
        indirBlkOfWins = indirBlkOfWins2+indirBlkOfWins1
        dimsI = dimsI2+dimsI1
        dimsJ = dimsJ2+dimsJ1
        dimsK = dimsK2+dimsK1
        #
        # Gather Matching cells into structured patches [ [[noz1,noz2],[imin1,imax1,...],[imin2,imax2,...],trirac] ]
        if allTags != []:
            infos = Connector.gatherMatchingNM(allWins, allTags, typeOfWins, indirBlkOfWins, dimsI, dimsJ, dimsK, dim, tol)
        else: infos = []
        for info in infos:
            noz1 = info[0][0]; noz2 = info[0][1]
            dimZ = Internal.getZoneDim(zones[noz1])[4]
            noz2 = noz2-nzones
            dims1 = Internal.getZoneDim(zones[noz1])
            ni1=dims1[1]; nj1=dims1[2]; nk1=dims1[3]
            dimsp1 = Internal.getZoneDim(zones2[noz1])
            nip1=dimsp1[1]; njp1=dimsp1[2]; nkp1=dimsp1[3]
            now1 = info[4][0]; now2 = info[4][1]
            if now1 != now2 and noz1 != noz2:
                range1 = info[1]; range2 = info[2]
                topp0 = info[3]
                now1 = info[4][0]; now2 = info[4][1]
                if dimZ == 3: topp = [1,2,3]
                else:
                    topp = [1,2]
                    topp0 = [topp0[0], topp0[1]]

                if topp0[0] > 0: topp[topp0[0]-1] = 1
                else: topp[-topp0[0]-1] = -1
                if topp0[1] > 0: topp[topp0[1]-1] = 2
                else: topp[-topp0[1]-1] = -2
                if dimZ == 3:
                    if topp0[2] > 0: topp[topp0[2]-1] = 3
                    else: topp[-topp0[2]-1] = -3

                # addBC2Zone...
                name1 = 'nmatch%d_%d'%(noz1+1,glob); glob+=1
                name2 = 'nmatch%d_%d'%(noz2+1,glob); glob+=1
                # Get Parent Nodes
                r1 = Internal.getParentOfNode(a, zones[noz1]); parent1 = r1[0]; d1 = r1[1]
                r2 = Internal.getParentOfNode(a, zones[noz2]); parent2 = r2[0]; d2 = r2[1]
                # addBCNearMatch
                rangenm1 = Internal.window2Range(range1)# copy
                if ratios[0] != 1: rangenm1 = G.refineBCRanges__(rangenm1, nip1, njp1, nkp1, ni1, nj1, nk1, 1, ratios[0])
                if ratios[1] != 1: rangenm1 = G.refineBCRanges__(rangenm1, nip1, njp1, nkp1, ni1, nj1, nk1, 2, ratios[1])
                if ratios[2] != 1: rangenm1 = G.refineBCRanges__(rangenm1, nip1, njp1, nkp1, ni1, nj1, nk1, 3, ratios[2])
                rangenm1 = Internal.range2Window(rangenm1)
                C._addBC2Zone(zones[noz1],name1,'BCNearMatch',rangenm1,zones[noz2],range2  , topp0)
                C._addBC2Zone(zones[noz2],name2,'BCNearMatch',range2  ,zones[noz1],rangenm1, topp )

                # couplage RANS/laminar ou euler
                model_z1 = model; model_z2 = model
                eq = Internal.getNodeFromName2(zones[noz1], 'GoverningEquations')
                if eq is not None: model_z1 = Internal.getValue( eq )
                eq = Internal.getNodeFromName2(zones[noz2], 'GoverningEquations')
                if eq is not None: model_z2 = Internal.getValue( eq )

                if model_z1 == 'NSTurbulent' and  model_z1 != model_z2:
                    #creation flag pour tranfert rans/LES
                    datap1 = numpy.ones(1, dtype=Internal.E_NpyInt)
                    datap2 = numpy.ones(1, dtype=Internal.E_NpyInt)
                    Internal.createUniqueChild( Internal.getNodeFromName2(zones[noz1], name1) , 'RANSLES', 'DataArray_t', datap1)
                    Internal.createUniqueChild( Internal.getNodeFromName2(zones[noz2], name2) , 'RANSLES', 'DataArray_t', datap2)
                    name1 = 'RANS_LES%d_%d'%(noz1,noz2)
                    C._addBC2Zone(zones[noz1],name1,'BCExtrapolateRANS',rangenm1)

                if model_z2 =='NSTurbulent' and  model_z1 != model_z2:
                    datap1 = numpy.ones(1, dtype=Internal.E_NpyInt)
                    datap2 = numpy.ones(1, dtype=Internal.E_NpyInt)
                    Internal.createUniqueChild( Internal.getNodeFromName2(zones[noz2], name2) , 'RANSLES', 'DataArray_t', datap2)
                    Internal.createUniqueChild( Internal.getNodeFromName2(zones[noz1], name1) , 'RANSLES', 'DataArray_t', datap1)
                    name2 = 'RANS_LES%d_%d'%(noz2,noz1)
                    C._addBC2Zone(zones[noz2],name2,'BCExtrapolateRANS',range2  )

                parent1[2][d1] = zones[noz1]; parent2[2][d2] = zones[noz2]

    return Internal.pyTree2Node(a,typen)

#=============================================================================
# blank intersecting and negative volume cells in a mesh. Useful for
# strand grids.
# Set the cellN to 0 for invalid cells, 2 for neighbouring cells
#=============================================================================
def blankIntersectingCells(t, tol=1.e-10, depth=2):
    """Set the cellN at 0 for intersecting cells in the normal direction
    to the wall.
    Usage: blankIntersectingCells(t, tol, depth)"""
    a = Internal.copyRef(t)
    _addCellN__(a, loc='centers')
    coords = C.getFields(Internal.__GridCoordinates__, a, api=1)
    cellN = C.getField('centers:cellN', a, api=1)
    res = Connector.blankIntersectingCells(coords, cellN, tol)
    C.setFields(res, a, 'centers')
    return a

#==============================================================================
# Masquage
#==============================================================================
def blankCells(t, bodies, blankingMatrix=[], depth=2,
               blankingType='cell_intersect', delta=1.e-10, dim=3,
               tol=1.e-8, XRaydim1=1000, XRaydim2=1000, cellNName='cellN'):
    try: import Transform as T
    except: raise ImportError("blankCells: requires Transform module.")
    if depth != 1 and depth != 2:
        print('Warning: blankCells: depth must be equal to 1 or 2. Set to default value (2).')
        depth = 2

    if blankingType != 'cell_intersect' and \
            blankingType != 'cell_intersect_opt' and \
            blankingType != 'center_in' and \
            blankingType != 'node_in':
        print('Warning: blankCells: blankingType must be cell_intersect, cell_intersect_opt, center_in or node_in.')
        print('Set to default (cell_intersect).')
        blankingType = 'cell_intersect'

    blankType = 1 # par defaut: cell_intersect
    if blankingType == 'node_in': blankType = 0
    elif blankingType == 'cell_intersect': blankType = 1
    elif blankingType == 'center_in': blankType = 2
    elif blankingType == 'cell_intersect_opt':
        if depth == 2: blankType = -2
        else: blankType = -1
    else:
        print('Warning: blankCells: blankingType must be cell_intersect, cell_intersect_opt, center_in or node_in.')
        print('Set to default (cell_intersect).')
        blankType = 1

    nb = 0
    a = Internal.copyRef(t)
    # ajout du celln aux centres si n'existe pas pour une zone
    loc = 'centers'
    if blankType == 0: loc = 'nodes'
    _addCellN__(a, loc=loc, cellNName=cellNName)
    bases = Internal.getBases(a)
    if bases == []: raise ValueError("blankCells: no CGNS base found in input tree.")

    if isinstance(blankingMatrix, list) and blankingMatrix == []: blankingMatrix = numpy.ones((len(bases), len(bodies)), Internal.E_NpyInt)
    for b in bases:
        coords = C.getFields(Internal.__GridCoordinates__, b, api=1)
        if coords != []:
            if loc == 'centers': cellN = C.getField('centers:'+cellNName, b, api=1)
            else: cellN = C.getField(cellNName, b, api=1)
            for nb2 in range(len(bodies)):
                blanking = blankingMatrix[nb, nb2]
                if bodies[nb2] != [] and (blanking == 1 or blanking == -1):
                    bc = []
                    for z in bodies[nb2]:
                        c = C.getFields(Internal.__GridCoordinates__, z, api=1)
                        if c != []:
                            c = c[0]
                            if len(c) == 5: # structure
                                # pour le 2D
                                if c[2] == 2: c = T.reorder(c, (-3,1,2))
                                elif c[3] == 2: c = T.reorder(c, (1,-3,2))
                            bc.append(c)
                    masknot = 0
                    if blanking == -1: masknot = 1
                    cellN = Connector.blankCells(
                        coords, cellN, bc, blankingType=blankType, \
                        delta=delta, dim=dim, masknot=masknot, tol=tol, \
                        XRaydim1=XRaydim1, XRaydim2=XRaydim2, cellNName=cellNName)
            C.setFields(cellN, b, loc, False)
        nb += 1
    return a

def _blankCells(a, bodies, blankingMatrix=[], depth=2,
                blankingType='cell_intersect', delta=1.e-10, dim=3,
                tol=1.e-8, XRaydim1=1000, XRaydim2=1000, cellNName='cellN'):
    try: import Transform as T
    except: raise ImportError("_blankCells: requires Transform module.")
    if depth != 1 and depth != 2:
        print('Warning: blankCells: depth must be equal to 1 or 2. Set to default value (2).')
        depth = 2
    if blankingType == 'center_in':
        print('Warning: blankCells: cannot be applied yet with center_in.')
        blankingType='cell_intersect'

    if blankingType != 'cell_intersect' and \
            blankingType != 'cell_intersect_opt' and \
            blankingType != 'center_in' and \
            blankingType != 'node_in':
        print('Warning: blankCells: blankingType must be cell_intersect, cell_intersect_opt, center_in or node_in.')
        print('Set to default (cell_intersect).')
        blankingType = 'cell_intersect'

    blankType = 1 # par defaut: cell_intersect
    if blankingType == 'node_in': blankType = 0
    elif blankingType == 'cell_intersect': blankType = 1
    elif blankingType == 'center_in': blankType = 2
    elif blankingType == 'cell_intersect_opt':
        if depth == 2: blankType = -2
        else: blankType = -1
    else:
        print('Warning: blankCells: blankingType must be cell_intersect, cell_intersect_opt, center_in or node_in.')
        print('Set to default (cell_intersect).')
        blankType = 1

    nb = 0

    # ajout du celln aux centres si n'existe pas pour une zone
    loc = 'centers'
    if blankType == 0: loc = 'nodes'
    _addCellN__(a, loc=loc, cellNName=cellNName)
    bases = Internal.getBases(a)
    if bases == []: raise ValueError("_blankCells: no CGNS base found in input tree.")

    if isinstance(blankingMatrix, list) and blankingMatrix == []: blankingMatrix = numpy.ones((len(bases), len(bodies)), dtype=Internal.E_NpyInt)
    for b in bases:
        coords = C.getFields(Internal.__GridCoordinates__, b, api=3)
        if coords != []:
            if loc == 'centers': cellN = C.getField('centers:'+cellNName, b, api=3)
            else: cellN = C.getField(cellNName, b, api=3)
            for nb2 in range(len(bodies)):
                blanking = blankingMatrix[nb, nb2]
                if bodies[nb2] != [] and (blanking == 1 or blanking == -1):
                    bc = []
                    for z in bodies[nb2]:
                        c = C.getFields(Internal.__GridCoordinates__, z, api=1)
                        if c != []:
                            c = c[0]
                            if len(c) == 5: # structure
                                # pour le 2D
                                if c[2] == 2: c = T.reorder(c, (-3,1,2))
                                elif c[3] == 2: c = T.reorder(c, (1,-3,2))
                            bc.append(c)
                    masknot = 0
                    if blanking == -1: masknot = 1
                    Connector._blankCells(coords, cellN, bc, blankingType=blankType, \
                                          delta=delta, dim=dim, masknot=masknot, tol=tol,\
                                          XRaydim1=XRaydim1, XRaydim2=XRaydim2, cellNName=cellNName)
        nb += 1
    return None

#==============================================================================
# Masquage par Tetra
#==============================================================================
def blankCellsTetra(t, mT4, blankingMatrix=[], blankingType='node_in',
                    tol=1.e-12, cellnval=0, overwrite=0, cellNName='cellN'):
    try: import Transform as T
    except: raise ImportError("blankCells: requires Transform module.")

    blankType = 1 # par defaut: cell_intersect
    if blankingType == 'node_in': blankType = 0
    elif blankingType == 'cell_intersect': blankType = 1
    elif blankingType == 'center_in': blankType = 2
    else:
        print('Warning: blankCellsTetra: blankingType must be cell_intersect, center_in or node_in.')
        print('Set to default (cell_intersect).')
        blankType = 1

    nb = -1
    a = Internal.copyRef(t)
    # ajout du celln aux centres si n'existe pas pour une zone
    loc = 'centers'
    if blankType == 0: loc = 'nodes'
    _addCellN__(a, loc=loc, cellNName=cellNName)
    bases = Internal.getBases(a)
    if bases == []: raise ValueError("blankCellsTetra: no basis found in input tree.")

    if isinstance(blankingMatrix, list) and blankingMatrix == []: blankingMatrix = numpy.ones((len(bases), len(mT4)), dtype=Internal.E_NpyInt)
    for b in bases:
        nb += 1
        coords = C.getFields(Internal.__GridCoordinates__, b, api=1)
        if coords == []: continue

        if len(coords[0]) == 5: coords = Converter.convertArray2Hexa(coords) # STRUCT -> HEXA

        if loc == 'centers': cellN = C.getField('centers:'+cellNName, b, api=1)
        else: cellN = C.getField(cellNName, b)
        bc = []
        for nb2 in range(len(mT4)):
            blanking = blankingMatrix[nb, nb2]
            if mT4[nb2] == [] or (blanking != 1 and blanking != -1): continue
            i = 0
            for z in mT4[nb2]:
                c = C.getFields(Internal.__GridCoordinates__, z, api=1)
                if c != []:
                    c = c[0]
                    bc.append(c)
                i += 1
        if bc == []:
            #print('Warning: nothing to mask for base %d'%(nb))
            continue
        bc = T.join(bc)
        cellN = Connector.blankCellsTetra(coords, cellN, bc, blankingType=blankType, tol=tol, \
                                          cellnval=cellnval, overwrite=overwrite, cellNName=cellNName)
        bc = None
        coords = None
        C.setFields(cellN, b, loc, False)
    return a

#==============================================================================
# Masquage par Tri (surface Tri)
#==============================================================================
def blankCellsTri(t, mT3, blankingMatrix=[], blankingType='node_in',
                  tol=1.e-12, cellnval=0, overwrite=0, cellNName='cellN'):
    try: import Transform as T
    except: raise ImportError("blankCellsTri: requires Transform module.")

    blankType = 1 # par defaut: cell_intersect
    if blankingType == 'node_in': blankType = 0
    elif blankingType == 'cell_intersect': blankType = 1
    elif blankingType == 'center_in': blankType = 2
    else:
        print('Warning: blankCellsTri: blankingType must be cell_intersect, center_in or node_in.')
        print('Set to cell_intersect.')
        blankType = 1

    nb = -1
    a = Internal.copyRef(t)
    # ajout du celln aux centres si n'existe pas pour une zone
    loc = 'centers'
    if blankType == 0: loc = 'nodes'
    _addCellN__(a, loc=loc,cellNName=cellNName)
    bases = Internal.getBases(a)
    if bases == []: raise ValueError("blankCellsTri: no basis found in input tree.")

    if isinstance(blankingMatrix, list) and blankingMatrix == []: blankingMatrix = numpy.ones((len(bases), len(mT3)), dtype=Internal.E_NpyInt)
    for b in bases:
        nb += 1
        coords = C.getFields(Internal.__GridCoordinates__, b, api=1)
        if coords == []: continue

        if len(coords[0]) == 5: coords = Converter.convertArray2Hexa(coords) # STRUCT -> HEXA

        if loc == 'centers': cellN = C.getField('centers:'+cellNName, b, api=1)
        else: cellN = C.getField(cellNName, b, api=1)
        bc = []
        for nb2 in range(len(mT3)):
            blanking = blankingMatrix[nb, nb2]
            if mT3[nb2] == [] or (blanking != 1 and blanking != -1): continue
            i = 0
            for z in mT3[nb2]:
                c = C.getFields(Internal.__GridCoordinates__, z, api=1)
                if c != []:
                    c = c[0]
                    bc.append(c)
                i += 1
        if bc == []:
            #print('Warning: nothing to mask for base %d'%(nb))
            continue
        bc = Converter.convertArray2Tetra(bc); bc = T.join(bc)
        cellN = Connector.blankCellsTri(coords, cellN, bc, blankingType=blankType, tol=tol, \
                                        cellnval=cellnval, overwrite=overwrite, cellNName=cellNName)
        bc = None
        coords = None
        C.setFields(cellN, b, loc, False)
    return a

def _blankCellsTri(a, mT3, blankingMatrix=[], blankingType='node_in',
                   tol=1.e-12, cellnval=0, overwrite=0, cellNName='cellN'):
    try: import Transform as T
    except: raise ImportError("blankCellsTri: requires Transform module.")

    blankType = 1 # par defaut: cell_intersect
    if blankingType == 'node_in': blankType = 0
    elif blankingType == 'cell_intersect': blankType = 1
    elif blankingType == 'center_in': blankType = 2
    else:
        print('Warning: blankCellsTri: blankingType must be cell_intersect, center_in or node_in.')
        print('Set to cell_intersect.')
        blankType = 1

    nb = -1
    # ajout du celln aux centres si n'existe pas pour une zone
    loc = 'centers'
    if blankType == 0: loc = 'nodes'
    _addCellN__(a, loc=loc, cellNName=cellNName)
    bases = Internal.getBases(a)
    if bases == []: raise ValueError("blankCellsTri: no basis found in input tree.")

    if isinstance(blankingMatrix, list) and blankingMatrix == []: blankingMatrix = numpy.ones((len(bases), len(mT3)), dtype=Internal.E_NpyInt)
    for b in bases:
        nb += 1
        coords = C.getFields(Internal.__GridCoordinates__, b, api=1)
        if coords == []: continue

        if len(coords[0]) == 5: coords = Converter.convertArray2Hexa(coords) # STRUCT -> HEXA

        if loc == 'centers': cellN = C.getField('centers:'+cellNName, b, api=1)
        else: cellN = C.getField(cellNName, b, api=1)
        bc = []
        for nb2 in range(len(mT3)):
            blanking = blankingMatrix[nb, nb2]
            if mT3[nb2] == [] or (blanking != 1 and blanking != -1): continue
            i = 0
            for z in mT3[nb2]:
                c = C.getFields(Internal.__GridCoordinates__, z, api=1)
                if c != []:
                    c = c[0]
                    bc.append(c)
                i += 1
        if bc == []:
            #print('Warning: nothing to mask for base %d'%(nb))
            continue
        bc = Converter.convertArray2Tetra(bc); bc = T.join(bc)
        cellN = Connector.blankCellsTri(coords, cellN, bc, blankingType=blankType, tol=tol, \
                                        cellnval=cellnval, overwrite=overwrite, cellNName=cellNName)
        bc = None; coords = None
        C.setFields(cellN, b, loc, False)
    return None

# cellN modifications
def _modCellN1(t, cellNName='cellN'):
    return C.__TZA3(t, 'centers', Connector._modCellN1, cellNName)
def _modCellN2(t, cellNName='cellN'):
    return C.__TZA3(t, 'centers', Connector._modCellN2, cellNName)

#=====================================================================================
# returns the numpys of indices of cellN=2 cell centers and corresponding coordinates
#=====================================================================================
def getInterpolatedPoints(z, loc='centers', cellNName='cellN'):
    if loc == 'centers':
        zc = C.node2Center(z)
        return connector.getInterpolatedPointsZ(zc, cellNName, Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__,Internal.__FlowSolutionCenters__)
    else:
        return connector.getInterpolatedPointsZ(z, cellNName, Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__,Internal.__FlowSolutionCenters__)

#==============================================================================
# optimisation du recouvrement
#==============================================================================
def optimizeOverlap(t, double_wall=0, priorities=[], planarTol=0., intersectionsDict=None):
    try: import Generator.PyTree as G
    except: raise ImportError('optimizeOverlap requires Generator module.')
    try: import Post.PyTree as P
    except: raise ImportError('optimizeOverlap requires Post module.')
    if double_wall == 1: from . import DoubleWall
    tol = 1.e-10
    a = Internal.copyRef(t)

    #=====================================================
    # 1-ajout du celln si n'existe pas pour une zone
    #=====================================================
    _addCellN__(a, loc='centers')
    G._getVolumeMap(a)
    bases = Internal.getBases(a)
    nbases = len(bases)

    #=====================================================
    # 2-traitement priorites
    #=====================================================
    nprios = len(priorities)//2
    prios = []
    size = 0
    if nprios == 0:
        for nob in range(nbases): prios.append(0)
    else:
        max_prio = 0
        for nob in range(nbases):
            baseName = bases[nob][0]
            prio = -1
            for nop in range(nprios):
                if priorities[2*nop] == baseName:
                    prio = priorities[2*nop+1]
                    max_prio = max(max_prio, prio)
                    break
            prios.append(prio)

        max_prio += 1
        nprios = len(prios)
        for nop in range(nprios):
            if prios[nop] == -1: prios[nop] = max_prio

    #=======================================================================
    # 2 - Recherche des periodicites:
    #     duplication des blocs periodiques dans les bases associees
    #     creation d'un noeud fils au niveau de la zone dupliquee de nom 'TemporaryPeriodicZone'
    #=======================================================================
    for nob in range(len(a[2])):
        if Internal.getType(a[2][nob]) == 'CGNSBase_t':
            C._addPeriodicZones__(a[2][nob])

    # Updates the intersection Dictionnary with new duplicated temporary zones:
    if not intersectionsDict:
        intersectionsDict = getIntersectingDomains(a, method='AABB')
    else:
        NewAndOldZones = Internal.getZones(a)
        #NAOZqty = len(NewAndOldZones)
        NewAndOldZonesNames = []
        #[NewAndOldZonesNames.append(NewAndOldZones[i][0]) for i in range(NAOZqty)]
        [NewAndOldZonesNames.append(z[0]) for z in NewAndOldZones]
        OldNames = intersectionsDict.keys()
        NewNames = []
        #for i in range(NAOZqty):
        #    if NewAndOldZonesNames[i] not in OldNames: NewNames.append(NewAndOldZonesNames[i])
        for zn in NewAndOldZonesNames:
            if zn not in OldNames: NewNames.append(zn)
        if not not NewNames:
            tDuplZones = C.newPyTree(['DuplicatedZones'])
            #for i in range(len(NewNames)): tDuplZones[2][1][2].append(Internal.getNodeFromName2(a,NewNames[i]))
            for newname in NewNames: tDuplZones[2][1][2].append(Internal.getNodeFromName2(a,newname))
            periodicZonesIntersectionDict = getIntersectingDomains(tDuplZones,a,method='hybrid')
            intersectionsDict.update(periodicZonesIntersectionDict)

    # Creation of extended centers meshes
    allExtCenters = []; allCenters = []; zones = []
    for nob1 in range(nbases):
        zones.append(Internal.getNodesFromType1(bases[nob1], 'Zone_t'))
        nodesPerBase = C.getFields(Internal.__GridCoordinates__, zones[nob1], api=1)
        if nodesPerBase == []:
            allExtCenters.append([[]]*len(zones[nob1]))
            allCenters.append([[]]*len(zones[nob1]))
        else:
            allExtCenters.append(Converter.node2ExtCenter(nodesPerBase))
            allCenters.append(Converter.node2Center(nodesPerBase))
        nodesPerBase = []

    #=======================================================
    # 4-Donor cell search: bbox intersection + adt creation
    #=======================================================
    # on cree par zone de chq base la liste des noms des domaines intersectants
    nobOfIntersectBasesAndZones=[]; allHooks=[]
    for nob1 in range(nbases):
        nobOfIntersectBasesAndZonesForBase1=[]
        allHooksForBase1 = []
        for noz1 in range(len(zones[nob1])):
            isIntersected = False
            nobOfIntersectBasesAndZonesForZone1=[]
            for nob2 in range(nbases):
                if nob2 != nob1:
                    nobOfIntersectZonesOfBase2 = [] # numero des zones de base2 intersectant z1
                    for noz2 in range(len(zones[nob2])):
                        if zones[nob2][noz2][0] in intersectionsDict[zones[nob1][noz1][0]]:
                            isIntersected = True
                            nobOfIntersectZonesOfBase2.append(noz2)
                    # par base2 opposee, dit si intersecte et les numeros dans base2 des zones intersectantes
                    nobOfIntersectBasesAndZonesForZone1 += [nob2, nobOfIntersectZonesOfBase2]
            nobOfIntersectBasesAndZonesForBase1.append(nobOfIntersectBasesAndZonesForZone1)
            if isIntersected:
                ae1 = allExtCenters[nob1][noz1]
                hook = Converter.createHook([ae1],'extractMesh')
            else: hook = None
            allHooksForBase1.append(hook)
        allHooks.append(allHooksForBase1)
        nobOfIntersectBasesAndZones.append(nobOfIntersectBasesAndZonesForBase1)

    #=====================================================
    # 5-optimisation du recouvrement
    #=====================================================
    isDW = 0
    if double_wall == 0:
        for nob1 in range(nbases-1):
            base1 = bases[nob1]
            zones1 = Internal.getNodesFromType1(base1, 'Zone_t')
            prio1 = prios[nob1]; noz1 = 0
            for noz1 in range(len(zones1)):
                z1 = zones1[noz1]
                isTempPeriodicZone1 = 0
                if Internal.getNodeFromName1(z1,'TempPeriodicZone') is None:
                    r1 = Internal.getParentOfNode(a, z1); parent1 = r1[0]; d1 = r1[1]

                else: isTempPeriodicZone1 = 1
                ae1 = allExtCenters[nob1][noz1]
                ac1 = allCenters[nob1][noz1]
                sol1 = C.getField('centers:cellN', z1, api=1)[0]
                vol1 = C.getField('centers:vol', z1, api=1)[0]
                ac1 = Converter.addVars([ac1,sol1,vol1])
                adt1 = allHooks[nob1][noz1]
                nobOfIntersectBasesAndZonesForZone1 = nobOfIntersectBasesAndZones[nob1][noz1]
                for nobi in range(len(nobOfIntersectBasesAndZonesForZone1)//2):
                    nob2 = nobOfIntersectBasesAndZonesForZone1[nobi*2]
                    if nob2 > nob1:
                        prio2 = prios[nob2]
                        base2 = bases[nob2]
                        zones2 = Internal.getNodesFromType1(base2,'Zone_t')
                        nobOfIntersectZones2 = nobOfIntersectBasesAndZonesForZone1[nobi*2+1]
                        for noz2 in nobOfIntersectZones2:
                            z2 = zones2[noz2]
                            adt2 = allHooks[nob2][noz2]
                            isTempPeriodicZone2 = 0
                            if Internal.getNodeFromName1(z2,'TempPeriodicZone') is None:
                                r2 = Internal.getParentOfNode(a, z2); parent2 = r2[0]; d2 = r2[1]
                            else: isTempPeriodicZone2 = 1
                            ae2 = allExtCenters[nob2][noz2]
                            ac2 = allCenters[nob2][noz2]
                            sol2 = C.getField('centers:cellN', z2, api=1)[0]
                            vol2 = C.getField('centers:vol', z2, api=1)[0]
                            ac2 = Converter.addVars([ac2, sol2, vol2])
                            res = Connector.optimizeOverlap__(ae1, ac1, ae2, ac2, prio1=prio1, prio2=prio2, \
                                                              isDW=isDW, hook1=adt1, hook2=adt2)
                            cellN1 = Converter.extractVars(res[0],['cellN'])
                            cellN2 = Converter.extractVars(res[1],['cellN'])
                            C.setFields([cellN1], z1, 'centers', False)
                            C.setFields([cellN2], z2, 'centers', False)
                            if isTempPeriodicZone1==0: parent1[2][d1] = z1
                            if isTempPeriodicZone2==0: parent2[2][d2] = z2

    else: #double wall
        dwInfo = DoubleWall.extractDoubleWallInfo__(a)
        firstWallCenters = dwInfo[0]; surfacesExtC = dwInfo[1]
        # liste des surfaces en centres etendus de toutes les zones de toutes les bases
        for nob1 in range(nbases-1):
            base1 = bases[nob1]
            zones1 = Internal.getNodesFromType1(base1, 'Zone_t')
            prio1 = prios[nob1]
            for noz1 in range(len(zones1)):
                z1 = zones1[noz1]
                isTempPeriodicZone1 = 0
                if Internal.getNodeFromName1(z1, 'TempPeriodicZone') is None:
                    r1 = Internal.getParentOfNode(a, z1); parent1 = r1[0]; d1 = r1[1]
                else: isTempPeriodicZone1 = 1
                ae1 = allExtCenters[nob1][noz1]
                ac1 = allCenters[nob1][noz1]
                sol1 = C.getField('centers:cellN', z1, api=1)[0]
                vol1 = C.getField('centers:vol', z1, api=1)[0]
                ac1 = Converter.addVars([ac1, sol1, vol1])
                adt1 = allHooks[nob1][noz1]

                firstWallCenters1 = firstWallCenters[nob1][noz1]
                surfacesExtC1 = surfacesExtC[nob1][noz1]

                # parcours des bases intersectantes
                nobOfIntersectBasesAndZonesForZone1 = nobOfIntersectBasesAndZones[nob1][noz1]
                for nobi in range(len(nobOfIntersectBasesAndZonesForZone1)//2):
                    nob2 = nobOfIntersectBasesAndZonesForZone1[nobi*2]
                    if nob2 > nob1:
                        prio2 = prios[nob2]
                        base2 = bases[nob2]
                        zones2 = Internal.getNodesFromType1(base2,'Zone_t')
                        nobOfIntersectZones2 = nobOfIntersectBasesAndZonesForZone1[nobi*2+1]
                        for noz2 in nobOfIntersectZones2:
                            z2 = zones2[noz2]
                            isTempPeriodicZone2 = 0
                            if Internal.getNodeFromName1(z2,'TempPeriodicZone') is None:
                                r2 = Internal.getParentOfNode(a, z2); parent2 = r2[0]; d2 = r2[1]
                            else: isTempPeriodicZone2 = 1
                            ae2 = allExtCenters[nob2][noz2]
                            ac2 = allCenters[nob2][noz2]
                            sol2 = C.getField('centers:cellN', z2, api=1)[0]
                            vol2 = C.getField('centers:vol', z2, api=1)[0]
                            ac2 = Converter.addVars([ac2,sol2,vol2])
                            adt2 = allHooks[nob2][noz2]
                            isDW = 0
                            firstWallCenters2 = firstWallCenters[nob2][noz2]
                            surfacesExtC2 = surfacesExtC[nob2][noz2]
                            if firstWallCenters1 != [] and firstWallCenters2 != []:
                                isDW = 1
                                acp1 = Converter.initVars(ac1,'{cellN}=minimum(2.,2*{cellN})')
                                acp2 = Converter.initVars(ac2,'{cellN}=minimum(2.,2*{cellN})')
                                acn1 = Connector.changeWall__(acp1, firstWallCenters1, surfacesExtC2, planarTol=planarTol)
                                acn2 = Connector.changeWall__(acp2, firstWallCenters2, surfacesExtC1, planarTol=planarTol)
                                cellN1 = Converter.extractVars(ac1,['cellN','vol'])
                                cellN2 = Converter.extractVars(ac2,['cellN','vol'])
                                acn1 = Converter.addVars([acn1,cellN1])
                                acn2 = Converter.addVars([acn2,cellN2])
                                res = Connector.optimizeOverlap__(ae1, acn1, ae2, acn2, prio1,prio2,isDW,\
                                                                  hook1=adt1, hook2=adt2)
                                cellN1 = Converter.extractVars(res[0],['cellN'])
                                cellN2 = Converter.extractVars(res[1],['cellN'])
                                ac1 = Converter.rmVars(ac1,['cellN']); ac1 = Converter.addVars([ac1,cellN1])
                                ac2 = Converter.rmVars(ac2,['cellN']); ac2 = Converter.addVars([ac2,cellN2])
                                C.setFields([cellN1], z1, 'centers', False)
                                C.setFields([cellN2], z2, 'centers', False)
                                if isTempPeriodicZone1==0: parent1[2][d1] = z1
                                if isTempPeriodicZone2==0: parent2[2][d2] = z2
                            else:
                                res = Connector.optimizeOverlap__(ae1, ac1, ae2, ac2, prio1,prio2,isDW,\
                                                                  hook1=adt1, hook2=adt2)
                                cellN1 = Converter.extractVars(res[0],['cellN'])
                                cellN2 = Converter.extractVars(res[1],['cellN'])
                                ac1 = Converter.rmVars(ac1,['cellN']); ac1 = Converter.addVars([ac1,cellN1])
                                ac2 = Converter.rmVars(ac2,['cellN']); ac2 = Converter.addVars([ac2,cellN2])
                                C.setFields([cellN1], z1, 'centers', False)
                                C.setFields([cellN2], z2, 'centers', False)

                                if isTempPeriodicZone1==0: parent1[2][d1] = z1
                                if isTempPeriodicZone2==0: parent2[2][d2] = z2
    # Delete duplicated periodic zones
    C._removeDuplicatedPeriodicZones__(a)

    C._rmVars(a, 'centers:vol')
    # remise a jour du celln : de 3 passe a 2
    C._initVars(a,'{centers:cellN}=minimum(2.,{centers:cellN})')
    #
    for nob1 in range(len(allHooks)):
        allHooksForBase1 = allHooks[nob1]
        for hookZ1 in allHooksForBase1:
            if hookZ1 is not None: Converter.freeHook(hookZ1)
    return a

#==============================================================================
def maximizeBlankedCells(t, depth=2, dir=1, loc='centers', cellNName='cellN', addGC=True):
    a = Internal.copyRef(t)
    _maximizeBlankedCells(a,depth=depth, dir=dir, loc=loc, cellNName=cellNName, addGC=addGC)
    return a

def _maximizeBlankedCells(t, depth=2, dir=1, loc='centers', cellNName='cellN', addGC=True):
    var = cellNName
    if loc == 'centers': var = 'centers:'+cellNName

    if addGC:
        ghost = Internal.getNodeFromName(t, 'ZoneRind')
        if ghost is None: Internal._addGhostCells(t, t, depth, adaptBCs=0, modified=[var])

    cellN = C.getField(var, t, api=1)
    cellN = Connector.maximizeBlankedCells(cellN, depth, dir, cellNName=cellNName)
    C.setFields(cellN, t, loc, False)

    if addGC:
        if ghost is None: Internal._rmGhostCells(t, t, depth, adaptBCs=0, modified=[var])
    return None

#==============================================================================
# Apply overlap BCs on cell nature field inside zone
# compatible avec une rangee de cellules d'interpolation
# Seulement pour les grilles structurees (no check)
#==============================================================================
# version in place (getFromArray2)
def _applyBCOverlapsStructured(z, depth, loc, val=2, cellNName='cellN', oversetFamNames=[]):
    varc = cellNName
    if loc == 'centers': varc = 'centers:'+varc; shift = 0
    else: shift = 1
    cellN = C.getField(varc, z, api=3)[0]
    ni = cellN[2]; nj = cellN[3]; nk = cellN[4]

    overlaps = Internal.getNodesFromType2(z, 'GridConnectivity_t')
    for o in overlaps:
        n = Internal.getNodeFromType(o, 'GridConnectivityType_t')
        if n is not None:
            v = Internal.getValue(n)
            if v == 'Overset':
                isDD = 0
                userDef = Internal.getNodesFromName(o, 'UserDefinedData')
                if userDef != []:
                    if len(userDef[0]) == 4:
                        info = userDef[0][2][0]
                        if info[0] == 'doubly_defined': isDD = 1
                if isDD == 0:
                    r = Internal.getNodesFromType1(o, 'IndexRange_t')
                    l = Internal.getNodesFromType1(o, 'IndexArray_t')
                    if r == [] and l == []:
                        print("Warning: applyBCOverlaps: BCOverlap is ill-defined.")
                    elif r != []:
                        rangew = r[0][1]
                        w = Internal.range2Window(rangew)
                        imin = w[0]; jmin = w[2]; kmin = w[4]
                        imax = w[1]; jmax = w[3]; kmax = w[5]
                        Connector._applyBCOverlapsStruct__(cellN,(imin,jmin,kmin),(imax,jmax,kmax),depth,loc,
                                                           val=val, cellNName=cellNName)
    # defined by a family with .Solver#Overlap
    # list of families of type overset
    for bc in Internal.getNodesFromType2(z, 'BC_t'):
        famName = Internal.getNodeFromType1(bc, 'FamilyName_t')
        if famName is not None:
            famName = Internal.getValue(famName)
            if famName in oversetFamNames:
                r = Internal.getNodesFromType1(bc, 'IndexRange_t')
                l = Internal.getNodesFromType1(bc, 'IndexArray_t')
                if r == [] and l == []:
                    print("Warning: applyBCOverlaps: BCOverlap is ill-defined.")
                elif r != []:
                    rangew = r[0][1]
                    w = Internal.range2Window(rangew)
                    imin = w[0]; jmin = w[2]; kmin = w[4]
                    imax = w[1]; jmax = w[3]; kmax = w[5]
                Connector._applyBCOverlapsStruct__(cellN,(imin,jmin,kmin),(imax,jmax,kmax),depth,loc, val=val, cellNName=cellNName)
    return None

# Version avec getField/setField - avec copie du tableau
def applyBCOverlapsStructured(z, depth, loc, val=2, cellNName='cellN',
                              oversetFamNames=[]):
    varc = cellNName
    if loc == 'centers': varc = 'centers:'+cellNName; shift = 0
    else: shift = 1
    cellN = C.getField(varc, z, api=1)[0]
    ni = cellN[2]; nj = cellN[3]; nk = cellN[4]

    overlaps = Internal.getNodesFromType2(z, 'GridConnectivity_t')
    for o in overlaps:
        n = Internal.getNodeFromType(o, 'GridConnectivityType_t')
        if n is not None:
            v = Internal.getValue(n)
            if v == 'Overset':
                isDD = 0
                userDef = Internal.getNodesFromName(o, 'UserDefinedData')
                if userDef != []:
                    if len(userDef[0]) == 4:
                        info = userDef[0][2][0]
                        if info[0] == 'doubly_defined': isDD = 1
                if isDD == 0:
                    r = Internal.getNodesFromType1(o, 'IndexRange_t')
                    l = Internal.getNodesFromType1(o, 'IndexArray_t')
                    if r == [] and l == []:
                        print("Warning: applyBCOverlaps: BCOverlap is ill-defined.")
                    elif r != []:
                        rangew = r[0][1]
                        w = Internal.range2Window(rangew)
                        imin = w[0]; jmin = w[2]; kmin = w[4]
                        imax = w[1]; jmax = w[3]; kmax = w[5]
                        cellN = Connector.applyBCOverlapsStruct__(cellN,(imin,jmin,kmin),(imax,jmax,kmax),depth,loc,
                                                                  val=val, cellNName=cellNName)
                        C.setFields([cellN], z, loc, False)
    # defined by a family with .Solver#Overlap
    # list of families of type overset
    for bc in Internal.getNodesFromType2(z,'BC_t'):
        famName = Internal.getNodeFromType1(bc,'FamilyName_t')
        if famName is not None:
            famName = Internal.getValue(famName)
            if famName in oversetFamNames:
                r = Internal.getNodesFromType1(bc, 'IndexRange_t')
                l = Internal.getNodesFromType1(bc, 'IndexArray_t')
                if r == [] and l == []:
                    print("Warning: applyBCOverlaps: BCOverlap is ill-defined.")
                elif r != []:
                    rangew = r[0][1]
                    w = Internal.range2Window(rangew)
                    imin = w[0]; jmin = w[2]; kmin = w[4]
                    imax = w[1]; jmax = w[3]; kmax = w[5]
                cellN = Connector.applyBCOverlapsStruct__(cellN,(imin,jmin,kmin),(imax,jmax,kmax),depth,loc, val=val, cellNName=cellNName)
                C.setFields([cellN], z, loc, False)
    return None

def applyBCOverlapsUnstructured(z, depth, loc, val=2, cellNName='cellN',oversetFamNames=[]):
    varc = cellNName
    if loc == 'centers': varc = 'centers:'+cellNName
    cellN = C.getField(varc, z, api=1)[0]
    zoneBC = Internal.getNodesFromType2(z, 'BC_t')
    for bc in zoneBC:
        v = Internal.getValue(bc)
        isOv = False
        if v == 'BCOverlap': isOv=True
        elif v=='FamilySpecified':
            famName = Internal.getNodeFromName1(bc,'FamilyName')
            if famName in oversetFamNames: isOv=True
        if isOv:
            faceListN = Internal.getNodesFromName(bc, Internal.__FACELIST__)
            if faceListN != []:
                faceList = Internal.getValue(faceListN[0])
                cellN = Connector.applyBCOverlapsNG__(cellN, faceList, depth, loc, val=val, cellNName=cellNName)
                C.setFields([cellN], z, loc, False)
    return None

def applyBCOverlaps(t, depth=2, loc='centers', val=2, cellNName='cellN'):
    a = Internal.copyRef(t)
    # ajout du celln si n'existe pas pour une zone
    _addCellN__(a, loc=loc, cellNName=cellNName)

    # only non doubly defined
    oversetFamNames = []
    for fam in Internal.getNodesFromType2(a, 'Family_t'):
        OV = Internal.getNodeFromName1(fam, '.Solver#Overlap')
        if OV is not None:
            dd = Internal.getNodeFromName1(OV, 'doubly_defined')
            if dd is None:
                oversetFamNames.append(Internal.getName(fam))

    zones = Internal.getZones(a)
    for z in zones:
        dimZ = Internal.getZoneDim(z)
        if dimZ[0] == 'Structured': applyBCOverlapsStructured(z, depth, loc, val=val, cellNName=cellNName, oversetFamNames=oversetFamNames)
        else:
            if dimZ[3] == 'NGON': applyBCOverlapsUnstructured(z, depth, loc, val=val, cellNName=cellNName, oversetFamNames=oversetFamNames)
            else:
                print('Warning: applyBCOverlaps: only for NGON unstructured zones.')
    return a

# VERSION getFromArray2 en structure
def _applyBCOverlaps(a, depth=2, loc='centers', val=2, cellNName='cellN', checkCellN=True):
    # ajout du celln si n'existe pas pour une zone
    if checkCellN: _addCellN__(a, loc=loc, cellNName=cellNName)
    oversetFamNames = []
    for fam in Internal.getNodesFromType2(a, 'Family_t'):
        OV = Internal.getNodeFromName1(fam, '.Solver#Overlap')
        if OV is not None:
            oversetFamNames.append(Internal.getName(fam))
    zones = Internal.getZones(a)
    for z in zones:
        dimZ = Internal.getZoneDim(z)
        if dimZ[0] == 'Structured': _applyBCOverlapsStructured(z, depth, loc, val, cellNName=cellNName, oversetFamNames=oversetFamNames)
        else:
            if dimZ[3] == 'NGON': applyBCOverlapsUnstructured(z, depth, loc, val, cellNName=cellNName, oversetFamNames=oversetFamNames)
            else:
                print('Warning: applyBCOverlaps: only for NGON unstructured zones.')
    return None

#==============================================================================
# IN: a: contains the cellN located at nodes or centers
# IN: depth can be positive or negative
# IN: dir=0 (directional), dir=1 (star), dir=2 (diamond), dir=3 (octaedre)
# Return depth layers of interpolated points at the fringe of blanked points
#==============================================================================
def setHoleInterpolatedPoints(a, depth=2, dir=0, loc='centers', cellNName='cellN'):
    """Set cellN=2 around cellN=0 points."""
    count = 0
    a = setHoleInterpolatedPoints__(a, depth, dir, count, loc, cellNName)
    count += 1

    ghost = Internal.getNodeFromName(a, 'ZoneRind')
    if loc == 'centers': varcelln = 'centers:'+cellNName
    else: varcelln = cellNName
    if ghost is None:
        a = Internal.addGhostCells(a, a, abs(depth), adaptBCs=0, modified=[varcelln])
        a = setHoleInterpolatedPoints__(a, depth, dir, count, loc, cellNName)
        a = Internal.rmGhostCells(a, a, abs(depth), adaptBCs=0, modified=[varcelln])
    return a

def setHoleInterpolatedPoints__(a, depth, dir, count, loc, cellNName='cellN'):
    if depth == 0: return a
    if loc == 'centers': varcelln = 'centers:'+cellNName
    else: varcelln = cellNName
    for z in Internal.getZones(a):
        dims = Internal.getZoneDim(z)
        if dims[0] == 'Unstructured' and count == 1: pass
        else: # passage ghost cells
            cellN = C.getField(varcelln, z, api=1)[0]
            if cellN != []:# cellN existe
                cellN = Connector.setHoleInterpolatedPoints(cellN, depth=depth, dir=dir, cellNName=cellNName)
                C.setFields([cellN], z, loc, False)
    return a

def _setHoleInterpolatedPoints__(a, depth, dir, count, loc, cellNName='cellN'):
    if depth == 0: return None
    if loc == 'centers': varcelln = 'centers:'+cellNName
    else: varcelln = cellNName
    count = 0
    for z in Internal.getZones(a):
        dims = Internal.getZoneDim(z)
        if dims[0] == 'Unstructured' and count == 1: pass
        else: # passage ghost cells
            cellN = C.getField(varcelln, z, api=3)[0]
            if cellN != []:
                Connector._setHoleInterpolatedPoints(cellN,depth=depth, dir=dir, cellNName=cellNName)
    return None

def _setHoleInterpolatedPoints(a, depth=2, dir=0, loc='centers',
                               cellNName='cellN', addGC=True):
    """Set cellN=2 around cellN=0 points."""
    count = 0
    _setHoleInterpolatedPoints__(a, depth, dir, count, loc, cellNName)

    if addGC:
        count += 1
        ghost = Internal.getNodeFromName(a, 'ZoneRind')
        if ghost is None:
            if loc == 'centers': varcelln = 'centers:'+cellNName
            else: varcelln = cellNName
            Internal._addGhostCells(a, a, abs(depth), adaptBCs=0, modified=[varcelln])
            _setHoleInterpolatedPoints__(a, depth, dir, count, loc, cellNName)
            Internal._rmGhostCells(a, a, abs(depth), adaptBCs=0, modified=[varcelln])
    return None

#=============================================================================
# Retourne la liste des zones donneuses definies dans la CL doubly defined
# ainsi que le cellN associe. Prend en compte les zones donneuses periodisees
#=============================================================================
def getDoublyDefinedDonorZones__(oversetgcnode, topTreeD):
    listOfDnrZones=[]; listOfDnrCellN=[]
    donorNames0 = []
    DNS = Internal.getValue(oversetgcnode)
    if DNS is not None: donorNames0 = DNS.split(",")# liste de noms de zones ou de familles
    else: # nom de famille a la Cannelle
        DFN = Internal.getNodeFromName2(oversetgcnode,'DonorFamilyName')
        if DFN is not None: donorNames0 = Internal.getValue(DFN).split(",")

    donorNames = []
    for donorName in donorNames0:
        famZones = C.getFamilyZones(topTreeD,donorName)
        if famZones != []:
            for zn in famZones: donorNames.append(Internal.getName(zn))
        else: donorNames.append(donorName)

    donorNames = list(set(donorNames))
    # Duplicated periodic zones
    donorNamesPer = []
    for zd in Internal.getZones(topTreeD):
        dupzoneinfo = Internal.getNodeFromName1(zd,'TempPeriodicZone')
        if dupzoneinfo is not None:
            donorName = Internal.getValue(dupzoneinfo)
            if donorName in donorNames: donorNamesPer.append(Internal.getName(zd))
    donorNames += donorNamesPer

    for donorName in donorNames:
        dnrZone = Internal.getNodeFromName2(topTreeD,donorName)
        if dnrZone is not None:
            coords = C.getFields(Internal.__GridCoordinates__, dnrZone, api=1)[0]
            listOfDnrZones.append(coords)
            cellN2 = C.getField('centers:cellN', dnrZone, api=1)[0]
            if cellN2 == []:
                print('Warning: setDoublyDefined: cellN init to 1 for zone %s.'%dnrZone[0])
                C._initVars(dnrZone,'centers:cellN',1.)
                cellN2 = C.getField('centers:cellN', dnrZone, api=1)[0]
            listOfDnrCellN.append(cellN2)

    return listOfDnrZones,listOfDnrCellN

#=============================================================================
# application de la CL doublement definie:
# remet a 1 le celln si cellule non interpolable
#=============================================================================
def setDoublyDefinedBC(t, depth=2):
    a = Internal.copyRef(t)
    _addCellN__(a, loc='centers')
    C._initVars(a, 'centers:cellN_dd', 1.)
    #=======================================================================
    # 2 - Recherche des periodicites :
    #     duplication des blocs periodiques dans les bases associees
    #     creation d un noeud fils au niveau de la zone dupliquee de nom 'TemporaryPeriodicZone'
    #=======================================================================
    for nob in range(len(a[2])):
        if Internal.getType(a[2][nob])=='CGNSBase_t':
            C._addPeriodicZones__(a[2][nob])

    zones = Internal.getZones(a)
    for z in zones:
        if Internal.getNodeFromName1(z,'TempPeriodicZone') is not None: pass
        else:
            cellNDD = C.getField('centers:cellN', z, api=1)[0]
            #(parent, d2) = Internal.getParentOfNode(a, z)
            overlaps = Internal.getNodesFromType2(z, 'GridConnectivity_t')
            coords = C.getFields(Internal.__GridCoordinates__, z, api=1)[0]
            for o in overlaps:
                n = Internal.getNodeFromType1(o, 'GridConnectivityType_t')
                if n is not None:
                    val = Internal.getValue(n)
                    if val == 'Overset':
                        userDef = Internal.getNodeFromName1(o, 'UserDefinedData')
                        if userDef is not None:
                            if len(userDef) == 4:
                                info = userDef[2][0]
                                if Internal.getName(info) == 'doubly_defined':
                                    r = Internal.getNodeFromType1(o,'IndexRange_t')
                                    win = Internal.range2Window(Internal.getValue(r))
                                    # recuperation des zones donneuses : a  partir de o[1] ou d une famille
                                    listOfInterpZones,cellns=getDoublyDefinedDonorZones__(o,a)
                                    print("Doubly defined: %s / %s : %i donor zones (periodic zones included)." % (z[0], o[0], len(listOfInterpZones)))

                                    # detection des pts non interpolables
                                    imin = win[0]; jmin = win[2]; kmin = win[4]
                                    imax = win[1]; jmax = win[3]; kmax = win[5]
                                    rangew = [int(imin),int(imax),int(jmin),int(jmax),int(kmin),int(kmax)]
                                    cellNDD = Connector.setDoublyDefinedBC(coords, cellNDD, listOfInterpZones,\
                                                                           cellns, rangew, depth)
            cellNDD[0] = 'cellN_dd'
            C.setFields([cellNDD], z, 'centers', False)
            #parent[2][d2] = z

    # Delete duplicated periodic zones
    C._removeDuplicatedPeriodicZones__(a)

    C._initVars(a,'{centers:cellN}=minimum({centers:cellN}*{centers:cellN_dd},2.)')
    C._rmVars(a,['centers:cellN_dd'])
    return a

#=============================================================================
# masque XRay: pierce pts
#=============================================================================
def maskXRay__(body, delta=0., dim=3, isNot=0, tol=1.e-8):
    """Create the pierce points of a X-Ray mask defined by body."""
    body = C.convertArray2Tetra(body)
    surf = C.getFields(Internal.__GridCoordinates__, body, api=1)
    pts = Connector.maskXRay__(surf, delta, dim, isNot, tol)
    return C.convertArrays2ZoneNode('XRayPts', [pts])

#==============================================================================
# Computes the connectivity & BC for NS/LBM interfaces
#==============================================================================
def connectNSLBM(t, tol=1.e-6, dim=3, type='all'):

    a,typen = Internal.node2PyTree(t)
    glob = 0
    #On recupere les zones structurees
    zones = []
    for z in Internal.getZones(a):
        if Internal.getZoneType(z) == 1: zones.append(z)
    # extract empty windows
    structTags,structWins,structIndirBlkOfWins,typeOfWins,dimsI,dimsJ,dimsK = \
        getEmptyWindowsInfoStruct__(zones, dim)

    model ='Euler'
    bases = Internal.getBases(a)
    if bases != []:
        c = Internal.getNodeFromName2(bases[0], 'GoverningEquations')
        if c is not None: model = Internal.getValue(c)

    # Identify matching cells for structured zones
    if structTags != []:
        structTags = Connector.identifyMatching(structTags,tol)
        structTags = Converter.extractVars(structTags,['tag1','tag2'])
        # Gather into structured patches [[[noz1,noz2],[imin1,imax1,...],[imin2,imax2,...],trirac]]
        infos = Connector.gatherMatching(structWins,structTags,typeOfWins,structIndirBlkOfWins,
                                         dimsI, dimsJ, dimsK, dim, tol)

        for info in infos:
            noz1 = info[0][0]; noz2 = info[0][1]
            range1 = info[1]; range2 = info[2]
            topp0 = info[3]
            dimZ = Internal.getZoneDim(zones[noz1])
            dimzone = dimZ[4]
            if dimzone == 3: topp = [1,2,3]
            else:
                topp = [1,2]
                topp0 = [topp0[0], topp0[1]]

            if topp0[0] > 0: topp[topp0[0]-1] = 1
            else: topp[-topp0[0]-1] = -1
            if topp0[1] > 0: topp[topp0[1]-1] = 2
            else: topp[-topp0[1]-1] = -2
            if dimzone == 3:
                if topp0[2] > 0: topp[topp0[2]-1] = 3
                else: topp[-topp0[2]-1] = -3
            #------------------------------------------
            # addBC2Zone...
            name1 = 'match%d_%d'%(noz1+1,glob); glob += 1
            name2 = 'match%d_%d'%(noz2+1,glob); glob += 1
            C._addBC2Zone(zones[noz1],name1,'BCMatch',range1,zones[noz2],range2,topp0)
            C._addBC2Zone(zones[noz2],name2,'BCMatch',range2,zones[noz1],range1,topp)

            # couplage  NS/LBM
            model_z1 = model; model_z2 = model
            eq= Internal.getNodeFromName2(zones[noz1], 'GoverningEquations')
            if eq is not None: model_z1 = Internal.getValue(eq)
            eq= Internal.getNodeFromName2(zones[noz2], 'GoverningEquations')
            if eq is not None: model_z2 = Internal.getValue(eq)

            if model_z1 == 'LBMLaminar' and model_z1 != model_z2:
                # creation flag pour transfert LBM/NS
                datap1 = numpy.ones(1, Internal.E_NpyInt)
                datap2 = numpy.ones(1, Internal.E_NpyInt)
                Internal.createUniqueChild(Internal.getNodeFromName2(zones[noz1],name1), 'NSLBM', 'DataArray_t', datap1)
                Internal.createUniqueChild(Internal.getNodeFromName2(zones[noz2],name2), 'NSLBM', 'DataArray_t', datap2)
                name_extrap = 'NS_LBM%d_%d'%(noz1,noz2)
                C._addBC2Zone(zones[noz1],name_extrap,'BCReconsLBM',range1) #Reconstruction des VDF cote LBM
                C._addBC2Zone(zones[noz2],name_extrap,'BCdimNS',range2)     #Permet de gerer l'adim cote NS

            if model_z2 == 'LBMLaminar' and model_z1 != model_z2:
                datap1 = numpy.ones(1, Internal.E_NpyInt)
                datap2 = numpy.ones(1, Internal.E_NpyInt)
                Internal.createUniqueChild(Internal.getNodeFromName2(zones[noz2], name2), 'NSLBM', 'DataArray_t', datap2)
                Internal.createUniqueChild(Internal.getNodeFromName2(zones[noz1], name1), 'NSLBM', 'DataArray_t', datap1)
                name_extrap = 'NS_LBM%d_%d'%(noz2,noz1)
                C._addBC2Zone(zones[noz2], name_extrap, 'BCReconsLBM', range2) #Idem : reconstruction des VDF cote LBM
                C._addBC2Zone(zones[noz1], name_extrap, 'BCdimNS', range1)     #Idem : gestion adim cote NS

    return Internal.pyTree2Node(a, typen)

#==============================================================================
# Double Wall treatment for chimera transfers (MPI friendly)
# Only modify tc
# familyBC* can either be a string or a list of strings
#==============================================================================
def _doubleWall(t, tc, familyBC1, familyBC2, ghostCells=False, check=False, surfaceCenters1=None):
    from . import DoubleWall

    if isinstance(familyBC1, str): familyBC1 = [familyBC1]
    if isinstance(familyBC2, str): familyBC2 = [familyBC2]

    listOfMismatch1 = []; listOfMismatch2 = []
    for b in Internal.getBases(t):
        for z in Internal.getZones(b):
            for f1 in familyBC1:
                wall1 = C.getFamilyBCs(z, f1)
                for w in wall1: listOfMismatch1.append(b[0]+'/'+z[0]+'/'+w[0])
            for f2 in familyBC2:
                wall2 = C.getFamilyBCs(z, f2)
                for w in wall2: listOfMismatch2.append(b[0]+'/'+z[0]+'/'+w[0])

    # project interpolated points (cellN=2) from listOfMismatch2 onto listOfMismatch1
    DoubleWall._changeWall2(t, tc, listOfMismatch1, listOfMismatch2, '_'.join(familyBC1), '_'.join(familyBC2), ghostCells, check, surfaceCenters1)

    return None

#==============================================================================
# Double Wall pre processing
# return surfaceCenters1 (array) to speed up _doubleWall calls
#==============================================================================
def initDoubleWall(t, familyBC1, check=False):
    from . import DoubleWall

    if isinstance(familyBC1, str): familyBC1 = [familyBC1]

    listOfMismatch1 = []
    for b in Internal.getBases(t):
        for z in Internal.getZones(b):
            for f1 in familyBC1:
                wall1 = C.getFamilyBCs(z, f1)
                for w in wall1: listOfMismatch1.append(b[0]+'/'+z[0]+'/'+w[0])

    return DoubleWall.getProjSurfaceForDoubleWall(t, listOfMismatch1, '_'.join(familyBC1), check)
