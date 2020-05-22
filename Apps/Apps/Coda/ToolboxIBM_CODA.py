import Converter.PyTree as C
import Converter.Internal as Internal
import Connector.ToolboxIBM as TIBM
import Transform.PyTree as T
import Converter.Internal as Internal
import Converter
import Generator.PyTree as G
import Post.PyTree as P
import Connector.OversetData as XOD
import numpy

#=============================================================================
# Extract info for skin post-processing
# INPUT : numpys of coordinates and fields to be projected onto the surface
# IN/OUT: surface defined by a CGNS/Python tree tb
#=============================================================================
def extractIBMWallFields(XCP, YCP, ZCP, arrayOfFields, tb, variables):
    VARLIST = Converter.getVarNames(arrayOfFields)
    dictOfVarNumber={}
    for var in variables:
        for nov in range(len(VARLIST)):
            if VARLIST[nov]==var:
                dictOfVarNumber[var]=nov
                break

    # 1. Creation of a CGNS zone O-D of cloud points
    zsize = numpy.empty((1,3), numpy.int32, order='F')
    zsize[0,0] = XCP.shape[0]; zsize[0,1] = 0; zsize[0,2] = 0
    z = Internal.newZone(name='IBW_Wall',zsize=zsize,ztype='Unstructured')
    gc = Internal.newGridCoordinates(parent=z)
    coordx = ['CoordinateX',XCP,[],'DataArray_t']
    coordy = ['CoordinateY',YCP,[],'DataArray_t']
    coordz = ['CoordinateZ',ZCP,[],'DataArray_t']
    gc[2] = [coordx,coordy,coordz]
    n = Internal.createChild(z, 'GridElements', 'Elements_t', [2,0])
    Internal.createChild(n, 'ElementRange', 'IndexRange_t', [1,0])
    Internal.createChild(n, 'ElementConnectivity', 'DataArray_t', None)
    FSN = Internal.newFlowSolution(name=Internal.__FlowSolutionNodes__,
                                   gridLocation='Vertex', parent=z)

    for varname in dictOfVarNumber:
        novar = dictOfVarNumber[varname]
        vararrayN = arrayOfFields[1][novar]
        FSN[2].append([varname,vararrayN, [],'DataArray_t'])


        
    dimPb = Internal.getNodeFromName(tb,'EquationDimension')
    if dimPb is None: 
        print('Warning: extractIBMWallFields: pb dimension is set to 3.')
        dimPb = 3
    else:
        dimPb = Internal.getValue(dimPb)
    # Force all the zones to be in a single CGNS base
    td = Internal.copyRef(tb)
    for nob in range(len(td[2])):
        b = td[2][nob]
        if b[3] == 'CGNSBase_t':                
            zones = Internal.getNodesFromType1(b, 'Zone_t')
            if zones != []:
                zones = C.convertArray2Tetra(zones)
                zones = T.join(zones); zones = G.close(zones)
                b[2] = [zones]
    for varname in dictOfVarNumber:
        C._initVars(td,varname,0.)
            
    td = P.projectCloudSolution(z, td, dim=dimPb)
    return td
   
# --------------------------------------------------------------------------------
# Creation of 0-D zones of name 'Zone#IBCD_*' such that their original zones can be
# retrieved in post processing
# Fonction is available in 
def createIBMWZones(tc,variables=[]):
    tw = C.newPyTree(['IBM_WALL'])
    for z in Internal.getZones(tc):
        ZSR = Internal.getNodesFromType2(z,'ZoneSubRegion_t')
        for IBCD in Internal.getNodesFromName(ZSR,"IBCD_*"):
            xPW = Internal.getNodesFromName(IBCD,"CoordinateX_PW")[0][1]
            yPW = Internal.getNodesFromName(IBCD,"CoordinateY_PW")[0][1]
            zPW = Internal.getNodesFromName(IBCD,"CoordinateZ_PW")[0][1]
            nptsW = xPW.shape[0]
            zw = G.cart((0,0,0),(1,1,1),(nptsW,1,1))
            COORDX = Internal.getNodeFromName2(zw,'CoordinateX'); COORDX[1]=xPW
            COORDY = Internal.getNodeFromName2(zw,'CoordinateY'); COORDY[1]=yPW
            COORDZ = Internal.getNodeFromName2(zw,'CoordinateZ'); COORDZ[1]=zPW
            if variables != [] and variables != None:
                FSN = Internal.newFlowSolution(parent=zw)
                for varo in variables:
                    fieldV = Internal.getNodeFromName2(IBCD,varo)
                    if fieldV is not None: 
                        C._initVars(zw,varo,0.)
                        fieldW = Internal.getNodeFromName2(FSN,varo)
                        fieldW[1] = fieldV[1]

            zw[0]=z[0]+"#"+IBCD[0]
            tw[2][1][2].append(zw)
    return tw

#------------------------------------------
# Creation of the subregions for IBM zones
#-----------------------------------------
varsRM = ["Density","utau","yplus","Pressure","Velocity*"]
def _createIB_ZSR(z, facelist, correctedPts, wallPts, imagePts, bctype, loc='faces'):
    zname = Internal.getName(z)
    nameSubRegion='IBCD_'+zname
    zsr = Internal.getNodesFromName1(z, nameSubRegion)
    # create new subregion for interpolations
    if zsr == []:
        dimZSR = numpy.zeros((1), numpy.int32)
        dimZSR[0] = correctedPts[1].shape[0]
        #v = numpy.fromstring(zname, 'c')
        z[2].append([nameSubRegion, dimZSR, [],'ZoneSubRegion_t'])
        info = z[2][len(z[2])-1]
        info[2].append(['PointList',facelist, [], 'IndexArray_t'])
        Internal.createChild(info,'GridLocation','GridLocation_t','FaceCenter')

    XOD._addIBCCoords__(z, z[0], correctedPts, wallPts, imagePts, bctype)
    # remove some nodes created by addIBCCoord but not useful for CODA
    zsr=Internal.getNodesFromName1(z,nameSubRegion)[0]
    for vrm in varsRM: Internal._rmNodesFromName(zsr,vrm)
    return None

# Convert Cartesian-octree structured mesh into an NS HEXA mesh
# hanging nodes information are computed too
# IN : t : structured mesh with nearmatch info
# OUT : tuple (z,hanging_elts_coarse, hanging_eltsf1, hanging_eltsf2)
def convertCart2NSMesh(t):
    import Converter
    # Conversion to NS HEXA with hanging nodes
    zones = Internal.getZones(t)
    nzones = len(zones)
    dictOfNobOfZones={}
    hashDict={}
    ncellsTot = 0
    for noz in range(nzones):
        zname = zones[noz][0]
        dictOfNobOfZones[zname]=noz
        if noz == 0:
            hashDict[zname] = 0
        else:
            znamep = zones[noz-1][0]
            dimZp = Internal.getZoneDim(zones[noz-1])
            if dimZp[4] == 2:
                sizezp = (dimZp[1]-1)*(dimZp[2]-1)
            else:
                sizezp = (dimZp[1]-1)*(dimZp[2]-1)*(dimZp[3]-1)
            hashDict[zname] = hashDict[znamep]+sizezp
            ncellsTot += sizezp

    dimZp = Internal.getZoneDim(zones[nzones-1])
    if dimZp[4] == 2:
        sizezp = (dimZp[1]-1)*(dimZp[2]-1)
    else:
        sizezp = (dimZp[1]-1)*(dimZp[2]-1)*(dimZp[3]-1)        
    ncellsTot += sizezp
    #print('ncellsTot = ', ncellsTot)

    HN_C = []
    HN_F1 = []
    HN_F2 = []

    for noz in range(nzones):
        z = zones[noz]
        dimZ = Internal.getZoneDim(z)
        niz = dimZ[1]; njz = dimZ[2]; nkz = dimZ[3]
        dimPb = dimZ[4]
        ni1 = niz-1; nj1 = njz-1; nk1 = nkz-1
        for gc in Internal.getNodesFromType(z,'GridConnectivity_t'):
            gctype = Internal.getNodeFromType(gc,'GridConnectivityType_t')
            if gctype is not None:
                gctype = Internal.getValue(gctype)
                if gctype=='Abutting':
                    PR = Internal.getNodeFromName(gc,'PointRange')
                    PRD = Internal.getNodeFromName(gc,'PointRangeDonor')
                    PR = Internal.getValue(PR)
                    PR = Internal.range2Window(PR)
                    PRD = Internal.getValue(PRD)
                    PRD = Internal.range2Window(PRD)
                    NMR = Internal.getNodeFromName(gc,'NMRatio')
                    NMR = Internal.getValue(NMR)
                    zdnrname = Internal.getValue(gc)
                    nozd = dictOfNobOfZones[zdnrname]
                    zdnr = zones[nozd]
                    dimZd = Internal.getZoneDim(zdnr)
                    nizd = dimZd[1]; njzd = dimZd[2]; nkzd = dimZd[3]
                    dimPb = dimZd[4]
                    nid1 = nizd-1; njd1 = njzd-1; nkd1 = nkzd-1
                    nid1njd1 = nid1*njd1
                    ni1nj1=ni1*nj1
                    i1 = PR[0]; i2 = PR[1]
                    j1 = PR[2]; j2 = PR[3]
                    k1 = PR[4]; k2 = PR[5]
                    id1 = PRD[0]; id2 = PRD[1]
                    jd1 = PRD[2]; jd2 = PRD[3]
                    kd1 = PRD[4]; kd2 = PRD[5]
                    #print(gc[0], PR, PRD, NMR)
                    coarseList=[]
                    fineList1 = []; fineList2 = []
                    shiftr = hashDict[z[0]]
                    shiftd = hashDict[zdnrname]

                    if dimPb == 2:
                        if NMR[0]==2:# z is coarse, zd is fine
                            idl = id1-1
                            if j1 > 1: j1 = j1-1
                            if jd1 > 1: jd1 = jd1-1
                            for i in range(i1-1,i2-1):
                                indr    = i    + (j1-1)*ni1 
                                indopp1 = idl  + (jd1-1)*nid1
                                indopp2 = idl+1+ (jd1-1)*nid1
                                idl = idl+2
                                coarseList.append(indr+shiftr)
                                fineList1.append(indopp1+shiftd)
                                fineList2.append(indopp2+shiftd)

                        elif NMR[1] == 2:
                            jdl = jd1-1
                            if i1 > 1: i1 = i1-1
                            if id1 > 1: id1 = id1-1
                            for j in range(j1-1,j2-1):
                                indr    = i1-1  + j*ni1 
                                indopp1 = id1-1 + jdl*nid1
                                indopp2 = id1-1 + (jdl+1)*nid1
                                jdl = jdl+2
                                coarseList.append(indr+shiftr)
                                fineList1.append(indopp1+shiftd)
                                fineList2.append(indopp2+shiftd)

                    else: #3D
                        if NMR[0]==2:# z is coarse, zd is fine
                            idl = id1-1
                            if j1 > 1: j1 = j1-1
                            if jd1 > 1: jd1 = jd1-1
                            if k1>1: k1=k1-1
                            if kd1>1: kd1 = kd1-1
                            for i in range(i1-1,i2-1):
                                indr    = i    + (j1-1)*ni1+(k1-1)*ni1nj1
                                indopp1 = idl  + (jd1-1)*nid1+(kd1-1)*nid1njd1
                                indopp2 = idl+1+ (jd1-1)*nid1+(kd1-1)*nid1njd1
                                idl = idl+2
                                coarseList.append(indr+shiftr)
                                fineList1.append(indopp1+shiftd)
                                fineList2.append(indopp2+shiftd)
                        elif NMR[1]==2:# z is coarse, zd is fine
                            jdl = jd1-1
                            if i1>1: i1 = i1-1
                            if id1>1: id1 = id1-1
                            if k1>1: k1=k1-1
                            if kd1>1: kd1 = kd1-1
                            for j in range(j1-1,j2-1):
                                indr    = i1-1  + j*ni1   + (k1-1)*ni1nj1
                                indopp1 = id1-1 + jdl*nid1+ (kd1-1)*nid1njd1
                                indopp2 = id1-1 + jdl*nid1+ (kd1-1)*nid1njd1
                                jdl = jdl+2
                                coarseList.append(indr+shiftr)
                                fineList1.append(indopp1+shiftd)
                                fineList2.append(indopp2+shiftd)

                        elif NMR[2]==2:# z is coarse, zd is fine
                            kdl = kd1-1
                            if i1>1: i1 = i1-1
                            if id1>1: id1 = id1-1
                            if j1>1: j1=j1-1
                            if jd1>1: jd1 = jd1-1
                            for k in range(k1-1,k2-1):
                                indr    = i1-1  + (j1-1)*ni1   + k*ni1nj1
                                indopp1 = id1-1 + (jd1-1)*nid1 + kdl*nid1njd1
                                indopp2 = id1-1 + (jd1-1)*nid1 + kdl*nid1njd1
                                kdl = kdl+2
                                coarseList.append(indr+shiftr)
                                fineList1.append(indopp1+shiftd)
                                fineList2.append(indopp2+shiftd)
                                
                    HN_C +=coarseList
                    HN_F1+=fineList1
                    HN_F2+=fineList2

    zones = C.convertArray2Hexa(zones)
    z = zones[0]
    for i in range(1, nzones):
        z = T.join([z,zones[i]])

    return (z, HN_C, HN_F1, HN_F2)
