# ANCIENNES FONCTIONS DE POST.IBM A NE PLUS UTILISER
import Converter.PyTree as C
import Converter.Internal as Internal
from . import PyTree as P
import Generator.PyTree as G
import Converter.Mpi as Cmpi
import Transform.PyTree as T
import Converter.Distributed as Distributed
import Connector.PyTree as X
import math, numpy
import Connector.connector as connector
import Converter.converter
import Connector.OversetData as XOD
import Connector.IBM as X_IBM

#====================================================================================
#===========================================================
# IN: ts: geometry tree with solution
# IN: Sref: reference surface area
# IN: alpha: angle of attack (x-y plane)
# IN: beta: angle of attack (x-z plane)
# IN: dimPb: dimension of the problem
# IN: verbose: verbose for pressure & friction forces (& coeffcients) & total cl & total cd I/O to screen
# OUT: geometry tree with additional solution fields
# I/O screen:Cl & Cd
#===========================================================
def loads0(ts, Sref=None, alpha=0., beta=0., dimPb=3, verbose=False):
    if Sref is None:
        C._initVars(ts, '__ONE__',1.)
        Sref = P.integ(ts, '__ONE__')[0];
        C._rmVars(ts, ['__ONE__', 'centers:vol'])

    RefState = Internal.getNodeFromType(ts,'ReferenceState_t')
    PInf     = Internal.getValue(Internal.getNodeFromName(RefState,"Pressure"))
    RoInf    = Internal.getValue(Internal.getNodeFromName(RefState,"Density"))
    VxInf    = Internal.getValue(Internal.getNodeFromName(RefState,"VelocityX"))
    VyInf    = Internal.getValue(Internal.getNodeFromName(RefState,"VelocityY"))
    VzInf    = Internal.getValue(Internal.getNodeFromName(RefState,"VelocityZ"))
    VInf2    = VxInf*VxInf+VyInf*VyInf+VzInf*VzInf
    VInf     = math.sqrt(VInf2)

    q      = 0.5*RoInf*VInf2
    qinv   = 1./q
    alpha  = math.radians(alpha)
    beta   = math.radians(beta)
    calpha = math.cos(alpha); cbeta = math.cos(beta)
    salpha = math.sin(alpha); sbeta = math.sin(beta)

    #===========================
    # Compute pressure forces
    #===========================
    zw = Internal.getZones(ts)
    zw = T.join(zw)
    Internal._rmNodesFromType(ts,'Zone_t')
    ts[2][1][2].append(zw)
    FSN = Internal.getNodeFromName(ts,Internal.__FlowSolutionNodes__)
    isPresent=False
    if Internal.getNodeFromName(FSN,'Cp') is None:
        C._initVars(ts, '{Cp}=-({Pressure}-%g)*%g'%(PInf,qinv))
    res = P.integNorm(ts, 'Cp')[0]
    res = [i/Sref for i in res]
    calpha = math.cos(alpha); cbeta = math.cos(beta)
    salpha = math.sin(alpha); sbeta = math.sin(beta)
    if dimPb == 3:
        cd = res[0]*calpha*cbeta + res[1]*salpha*cbeta - res[2]*sbeta
        cl = res[1]*calpha       - res[0]*salpha
    else:
        cd = res[0]*calpha + res[1]*salpha
        cl = res[1]*calpha - res[0]*salpha
    if verbose:
        print("Normalized pressure drag = %g and lift = %g"%(cd, cl))
        print("Vector of pressure loads: (Fx_P,Fy_P,Fz_P)=(",res[0],res[1],res[2],")")
    cd_press = cd ; cl_press = cl

    #======================================
    # Compute viscous forces
    #======================================
    C._initVars(ts, '{tau_wall}=0.')
    G._getNormalMap(ts)

    # Calcul de yplus au point image
    C._initVars(ts, '{dist_cw}=sqrt(({CoordinateX_PW}-{CoordinateX_PC})**2 + ({CoordinateY_PW}-{CoordinateY_PC})**2 + ({CoordinateZ_PW}-{CoordinateZ_PC})**2)')
    C._initVars(ts, '{dist_iw}=sqrt(({CoordinateX_PW}-{CoordinateX_PI})**2 + ({CoordinateY_PW}-{CoordinateY_PI})**2 + ({CoordinateZ_PW}-{CoordinateZ_PI})**2)')
    C._initVars(ts, '{yplus_i}=({yplus}/{dist_cw})*{dist_iw}')

    variables = ['Density', 'Cp', 'tau_wall', 'Pressure','VelocityX','VelocityY','VelocityZ']
    FSN = Internal.getNodeFromType(ts,'FlowSolution_t')
    if Internal.getNodeFromName1(FSN,'utau') is not None:
        variables += ['utau','yplus','tau_wall']
        C._initVars(ts, '{tau_wall}={Density}*{utau}*{utau}')

    isGradP = False
    if Internal.getNodeFromName1(FSN,'gradxPressure') is not None:
        variables += ['gradxPressure','gradyPressure','gradzPressure']
        isGradP = True

    if Internal.getNodeFromName1(FSN,'yplus_i') is not None:
        variables += ['yplus_i']

    ts = C.node2Center(ts, variables)
    Internal._rmNodesFromName(ts,'FlowSolution')
    C._normalize(ts, ['centers:sx','centers:sy','centers:sz'])
    # Calc tangent vector
    C._initVars(ts, '{centers:tx}={centers:VelocityX}-({centers:VelocityX}*{centers:sx}+{centers:VelocityY}*{centers:sy}+{centers:VelocityZ}*{centers:sz})*{centers:sx}')
    C._initVars(ts, '{centers:ty}={centers:VelocityY}-({centers:VelocityX}*{centers:sx}+{centers:VelocityY}*{centers:sy}+{centers:VelocityZ}*{centers:sz})*{centers:sy}')
    C._initVars(ts, '{centers:tz}={centers:VelocityZ}-({centers:VelocityX}*{centers:sx}+{centers:VelocityY}*{centers:sy}+{centers:VelocityZ}*{centers:sz})*{centers:sz}')
    C._normalize(ts, ['centers:tx','centers:ty','centers:tz'])

    C._initVars(ts, '{centers:tauxx}=2*{centers:tau_wall}*{centers:tx}*{centers:sx}')
    C._initVars(ts, '{centers:tauyy}=2*{centers:tau_wall}*{centers:ty}*{centers:sy}')
    C._initVars(ts, '{centers:tauzz}=2*{centers:tau_wall}*{centers:tz}*{centers:sz}')
    C._initVars(ts, '{centers:tauxy}={centers:tau_wall}*({centers:tx}*{centers:sy}+{centers:ty}*{centers:sx})')
    C._initVars(ts, '{centers:tauxz}={centers:tau_wall}*({centers:tx}*{centers:sz}+{centers:tz}*{centers:sx})')
    C._initVars(ts, '{centers:tauyz}={centers:tau_wall}*({centers:ty}*{centers:sz}+{centers:tz}*{centers:sy})')

    # Cacl friction forces
    C._initVars(ts, '{centers:Fricx}={centers:tauxx}*{centers:sx}+{centers:tauxy}*{centers:sy}+{centers:tauxz}*{centers:sz}')
    C._initVars(ts, '{centers:Fricy}={centers:tauxy}*{centers:sx}+{centers:tauyy}*{centers:sy}+{centers:tauyz}*{centers:sz}')
    C._initVars(ts, '{centers:Fricz}={centers:tauxz}*{centers:sx}+{centers:tauyz}*{centers:sy}+{centers:tauzz}*{centers:sz}')

    # Calc coefficient of friction
    C._initVars(ts, '{centers:Cf}=(sqrt({centers:Fricx}**2+{centers:Fricy}**2+{centers:Fricz}**2))/%g'%q)

    # Update of pressure gradients
    if isGradP:
        C._initVars(ts, '{centers:gradnP}={centers:gradxPressure}*{centers:sx}+{centers:gradyPressure}*{centers:sy}+{centers:gradzPressure}*{centers:sz}')
        C._initVars(ts, '{centers:gradtP}={centers:gradxPressure}*{centers:tx}+{centers:gradyPressure}*{centers:ty}+{centers:gradzPressure}*{centers:tz}')

    G._getVolumeMap(ts)
    effortX = P.integ(ts, 'centers:Fricx')[0]
    effortY = P.integ(ts, 'centers:Fricy')[0]
    effortZ = P.integ(ts, 'centers:Fricz')[0]

    QADIMI = 1./(q*Sref)
    if dimPb == 3:
        cd = (effortX*calpha*cbeta + effortY*salpha*cbeta - effortZ*sbeta)*QADIMI
        cl = (effortY*calpha       - effortX*salpha)*QADIMI
    else:
        cd = (effortX*calpha + effortY*salpha)*QADIMI
        cl = (effortY*calpha - effortX*salpha)*QADIMI

    vars = ['centers:sx'   , 'centers:sy'   , 'centers:sz',
            'centers:tx'   , 'centers:ty'   , 'centers:tz',
            'centers:tauxx', 'centers:tauyy', 'centers:tauzz',
            'centers:tauxy', 'centers:tauxz', 'centers:tauyz',
            'centers:Fricx', 'centers:Fricy', 'centers:Fricz']
    C._rmVars(ts, vars)

    if verbose:
        print("Normalized skin friction drag = %g and lift = %g"%(cd, cl))
        print("Vector of skin friction loads: (Fx_f,Fy_f,Fz_f)=(",effortX*QADIMI, effortY*QADIMI, effortZ*QADIMI,")")

    ##################################################
    if verbose:
        print("****************************************")
        print("Total Drag :", cd_press+cd)
        print("Total Lift :", cl_press+cl)
        print("****************************************")

    return ts
#=============================================================================
# Post forces
# IN: tb: geometry file with solution projected onto it
# IN: Sref: Reference surface area
# OUT: wall_out ou None: fichier pour sortie des efforts sur la paroi aux centres
# IN: alpha: angle pour les efforts
# IN: beta: angle pour les efforts
#==============================================================================
def unsteadyLoads(tb, Sref=None, Pref=None, Qref=None, alpha=0., beta=0.):
    """Computes the viscous and pressure forces on the IB during the computation of the solution. 
    Usage: unsteadyLoads(tb, Sref, Pref, Qref, alpha, beta)"""
    tp = Internal.copyRef(tb)
    return _unsteadyLoads(tp, Sref=Sref, Pref=Pref, Qref=Qref, alpha=alpha, beta=beta)

def _unsteadyLoads(tb, Sref=None, Pref=None, Qref=None, alpha=0., beta=0.):
    """Computes the viscous and pressure forces on the IB during the computation of the solution. 
    Usage: _unsteadyLoads(tb, Sref, Pref, Qref, alpha, beta, dimPb)"""
    return loads0(tb, Sref=Sref, alpha=alpha, beta=beta, dimPb=3, verbose=False)


# reconstruit les champs parietaux a partir des infos stockees dans tw et les champs de tc
def _computeWallReconstruction(tw, tcw, tc, procDictR=None, procDictD=None, graph=None,
                               variables=['Pressure','Density','utau','yplus']):
    if procDictR is None:
        procDictR={}
        for z in Internal.getZones(tw): procDictR[z[0]]=0
    if procDictD is None:
        procDictD={}
        for z in Internal.getZones(tcw): procDictD[z[0]]=0

    if graph is None:
        graph = Cmpi.computeGraph(tcw, type='POST',procDict=procDictD, procDict2=procDictR, t2=tw)

    cellNVariable=''; varType=1; compact=0

    datas={}
    # interpolation from a cloud of points
    for zc_w in Internal.getZones(tcw):
        infos = []
        zcw_name = Internal.getName(zc_w)
        zcw_namel = zcw_name.split('#')
        zc_orig_name = zcw_namel[0]; zsrname = zcw_namel[1]

        FSN = Internal.getNodeFromName2(zc_w,Internal.__FlowSolutionNodes__)
        if FSN is None:
            FSN = Internal.newFlowSolution(name=Internal.__FlowSolutionNodes__,
                                           gridLocation='Vertex', parent=zc_w)

        # put the fields in corresponding zcw zone
        zc_orig = Internal.getNodeFromName2(tc,zc_orig_name)
        zsr_orig = Internal.getNodeFromName2(zc_orig,zsrname)
        for varname in variables:
            varnode = Internal.getNodeFromName(zsr_orig,varname)
            FSN[2].append(varnode)

        # search for IDW_zrname node
        for zsr in Internal.getNodesFromType1(zc_w, "ZoneSubRegion_t"):
            sname = zsr[0][0:3]
            if sname=='IDW' and variables is not None:
                dname = Internal.getValue(zsr)
                idn = Internal.getNodeFromName1(zsr, 'InterpolantsDonor')
                if idn is not None: # la subRegion decrit des interpolations
                    zoneRole = Internal.getNodeFromName2(zsr, 'ZoneRole')
                    zoneRole = Internal.getValue(zoneRole)
                    if zoneRole == 'Donor':
                        location = Internal.getNodeFromName1(zsr, 'GridLocation') # localisation des donnees des receveurs

                        if location is not None:
                            location = Internal.getValue(location)
                            if location == 'CellCenter': loc = 'centers'
                            else: loc = 'nodes'

                        else: loc = 'nodes'

                        Coefs     = idn[1]
                        DonorType = Internal.getNodeFromName1(zsr,'InterpolantsType')[1]
                        ListDonor = Internal.getNodeFromName1(zsr,'PointList')[1]
                        ListRcv   = Internal.getNodeFromName1(zsr,'PointListDonor')[1]

                        arrayT = connector._setInterpTransfersD(zc_w, variables, ListDonor, DonorType, Coefs, varType, compact,
                                                                cellNVariable,
                                                                Internal.__GridCoordinates__,
                                                                Internal.__FlowSolutionNodes__,
                                                                Internal.__FlowSolutionCenters__)
                        infos.append([dname,arrayT,ListRcv,loc])
        for n in infos:
            rcvName = n[0]
            proc = procDictR[rcvName]
            if proc == Cmpi.rank:
                fields = n[1]
                if fields != []:
                    listIndices = n[2]
                    zr = Internal.getNodeFromName2(tw,rcvName)
                    C._updatePartialFields(zr, [fields], [listIndices], loc=n[3])
            else:
                rcvNode = procDictR[rcvName]
                if rcvNode not in datas: datas[rcvNode] = [n]
                else: datas[rcvNode] += [n]

    #send numpys using graph
    rcvDatas = Cmpi.sendRecv(datas, graph)

    # put contribution of donor zone to the interpolated fields in the receptor zone
    for i in rcvDatas:
        for n in rcvDatas[i]:
            rcvName = n[0]
            field = n[1]
            if field != []:
                listIndices = n[2]
                z = Internal.getNodeFromName2(tw, rcvName)
                C._updatePartialFields(z,[field], [listIndices], loc=n[3])
    return None

def prepareWallReconstruction(tw, tc):
    # MLS interpolation order
    LSorder = 2

    dimPb = Internal.getNodeFromName(tw, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)
    rank = Cmpi.rank
    # GLOBAL
    # pour eviter les trous lies au split des zones IBM en nuages de points
    tolBB = 0.
    for snear in Internal.getNodesFromName(tw,"snear"):
        snearval = Internal.getValue(snear)
        tolBB = max(tolBB,2*snearval)

    tcw = X_IBM.createIBMWZones(tc,variables=[])

    nzones = len(Internal.getZones(tcw))

    dimPb = Internal.getNodeFromName(tw, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)

    if dimPb ==2: C._initVars(tcw,"CoordinateZ",0.0)
    #
    # Create the bbtree
    tcw_bb = C.newPyTree(["Base"])
    zonesl = []
    for base in Internal.getBases(tcw):
        for z in Internal.getZones(base):
            zbb = G.BB(z)

            # Clean up (zoneSubRegion)
            Internal._rmNodesFromType(zbb, 'ZoneSubRegion_t')
            valxm = C.getValue(zbb,'CoordinateX',0)-tolBB
            valxp = C.getValue(zbb,'CoordinateX',1)-tolBB
            C.setValue(zbb,'CoordinateX',0, valxm)
            C.setValue(zbb,'CoordinateX',1, valxp)
            valxm = C.getValue(zbb,'CoordinateY',0)-tolBB
            valxp = C.getValue(zbb,'CoordinateY',1)-tolBB
            C.setValue(zbb,'CoordinateY',0, valxm)
            C.setValue(zbb,'CoordinateY',1, valxp)
            if dimPb == 3:
                valxm = C.getValue(zbb,'CoordinateZ',0)-tolBB
                valxp = C.getValue(zbb,'CoordinateZ',1)-tolBB
                C.setValue(zbb,'CoordinateZ',0, valxm)
                C.setValue(zbb,'CoordinateZ',1, valxp)
            Cmpi._setProc(zbb,rank)
            zonesl.append(zbb)

    zonesBB = Cmpi.allgatherZones(zonesl)
    tcw_bb[2][1][2] += zonesBB
    graph={}
    for zw in Internal.getZones(tw):
        zwBB = G.BB(zw)
        for zdBB in Internal.getZones(tcw_bb):
            if G.bboxIntersection(zwBB,zdBB,isBB=True):
                popp = Cmpi.getProc(zdBB)
                zdname = zdBB[0]
                Distributed.updateGraph__(graph, popp, rank, zdname)

    allGraph = Cmpi.KCOMM.allgather(graph)
    graph = {}
    for i in allGraph:
        for k in i:
            if not k in graph: graph[k] = {}
            for j in i[k]:
                if not j in graph[k]: graph[k][j] = []
                graph[k][j] += i[k][j]
                graph[k][j] = list(set(graph[k][j]))

    Cmpi._addXZones(tcw, graph, subr=False)
    interDict = X.getIntersectingDomains(tw,tcw_bb)
    procDict = Cmpi.getProcDict(tcw)

    graphX = {}
    for zw in Internal.getZones(tw):
        zwname = Internal.getName(zw)
        for zdname in interDict[zwname]:
            zdbb = Internal.getNodeFromName2(tcw_bb,zdname)
            popp = Cmpi.getProc(zdbb)
            Distributed.updateGraph__(graphX, rank, popp, zdname)

    allGraph = Cmpi.KCOMM.allgather(graphX)
    graphX = {}
    for i in allGraph:
        for k in i:
            if not k in graphX: graphX[k] = {}
            for j in i[k]:
                if not j in graphX[k]: graphX[k][j] = []
                graphX[k][j] += i[k][j]
                graphX[k][j] = list(set(graphX[k][j]))
    del tcw_bb

    EXTRAP = numpy.array([],numpy.int32)
    VOL = numpy.array([],numpy.float64)
    ORPHAN = numpy.array([],numpy.float64)
    datas = {}

    for zw in Internal.getZones(tw):
        zwname = Internal.getName(zw)
        dnrZones=[]
        dnrZoneNames=[]
        for zdname in interDict[zwname]:
            zd = Internal.getNodeFromName2(tcw,zdname)
            dnrZones.append(zd)
            dnrZoneNames.append(zdname)

        coordsD = C.getFields(Internal.__GridCoordinates__, dnrZones, api=1)
        coordsR = C.getFields(Internal.__GridCoordinates__, zw, api=1)[0]
        if coordsR[1].shape[1]>0:
            coordsD = Converter.convertArray2Node(coordsD)
            ret = connector.setInterpData_IBMWall(coordsD, coordsR, dimPb, LSorder)

            allPLR = ret[0]; allPLD = ret[1]; allCOEFS=ret[3]; allITYPE=ret[2]

            nozd = 0
            for zd in dnrZones:
                zdname = zd[0]
                XOD._createInterpRegion__(zd, zw[0], allPLD[nozd], allPLR[nozd], allCOEFS[nozd], allITYPE[nozd], \
                                          VOL, EXTRAP, ORPHAN, \
                                          tag='Donor', loc='nodes',itype='chimera', prefix='IDW_')
                nozd+=1

                destProc = procDict[zdname]

                IDs = []
                for i in zd[2]:
                    if i[0][0:3] == 'IDW':
                        if Internal.getValue(i)==zwname: IDs.append(i)

                if IDs != []:
                    if destProc == rank:
                        zD = Internal.getNodeFromName2(tcw, zdname)
                        zD[2] += IDs
                    else:
                        if destProc not in datas: datas[destProc] = [[zdname,IDs]]
                        else: datas[destProc].append([zdname,IDs])
                else:
                    if destProc not in datas: datas[destProc] = []
    for i in range(Cmpi.size):
        if i not in datas: datas[i]=[]
    Cmpi._rmXZones(tcw)
    destDatas = Cmpi.sendRecv(datas, graphX)
    for i in destDatas:
        for n in destDatas[i]:
            zname = n[0]
            IDs = n[1]
            if IDs != []:
                zD = Internal.getNodeFromName2(tcw, zname)
                zD[2] += IDs
    datas = {}; destDatas = None; graph={}
    return tcw
