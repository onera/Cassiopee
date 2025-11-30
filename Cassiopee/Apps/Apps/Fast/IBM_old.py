# Class for FastS "IBM" prepare and compute
import FastC.PyTree as FastC
import FastS.ToolboxChimera as TBX
import Converter.PyTree as C
import Geom.PyTree as D
import Generator.PyTree as G
import Transform.PyTree as T
import Post.PyTree as P
import Converter.Internal as Internal
import Connector.PyTree as X
import Dist2Walls.PyTree as DTW
import Distributor2.PyTree as D2
import Initiator.PyTree as I
import Compressor.PyTree as Compressor
import Converter.Mpi as Cmpi
import Connector.Mpi as Xmpi
import Post.Mpi as Pmpi
import Converter.Filter as Filter
import Converter.Distributed as Distributed
from Apps.Fast.Common import Common
import Connector.connector as connector
import Connector.OversetData as XOD
import Converter
import KCore.test as test
import Generator
import math
import numpy

import Geom.IBM as D_IBM
import Post.IBM as P_IBM
import Connector.IBM as X_IBM
import Generator.IBM as G_IBM
import Generator.IBMmodelHeight as G_IBM_Height

from mpi4py import MPI
COMM_WORLD = MPI.COMM_WORLD
KCOMM = COMM_WORLD

RENAMEIBCNODES=False # renommage des ibcd*

__IBCNameServer__={}

def getIBCDName(proposedName):
    global __IBCNameServer__
    (ibcname,__IBCNameServer__)=C.getUniqueName(proposedName, __IBCNameServer__)
    return ibcname

def _change_name_IBCD(tc2,NewIBCD):
    ZsubR=Internal.getNodesByType(tc2,'ZoneSubRegion_t')
    for z in ZsubR:
        zsplit=z[0].split('_')
        if zsplit[0]=='IBCD':
            zsplit[1]=str(NewIBCD)
            znew = '_'.join(zsplit)
            Internal.setName(z, znew)
    return None

# IN: maillage surfacique + reference State + snears
#================================================================================
# IBM prepare
# IN: t_case: fichier ou arbre body
# OUT: t_out, tc_out : fichier ou arbres de sorties
# snears: liste des snear, mais c'est mieux si ils sont dans t_case
# dfar, dfarList: liste des dfars, mais c'est mieux si ils sont dans t_case
# tbox: arbre de raffinement
# check: si true, fait des sorties
# NP: is the target number of processors for computation
# (maybe different from the number of processors the prep is run on)
# frontType=1,2,3: type de front
# expand=1,2,3: sorte d'expand des niveaux (1:classque,2:minimum,3:deux niveaux)
# tinit: arbre de champ d'avant pour la reprise
#================================================================================
def prepare(t_case, t_out, tc_out, snears=0.01, dfar=10., dfarList=[],
            tbox=None, snearsf=None, yplus=100.,
            vmin=21, check=False, NP=0, format='single',
            frontType=1, expand=3, tinit=None, initWithBBox=-1., wallAdapt=None,
            dfarDir=0, recomputeDist=False, redistribute=False):
    rank = Cmpi.rank; size = Cmpi.size
    ret = None
    ## sequential prep
    #if size == 1: ret = prepare0(t_case, t_out, tc_out, snears=snears, dfar=dfar, dfarList=dfarList,
    #                             tbox=tbox, snearsf=snearsf, yplus=yplus,
    #                             vmin=vmin, check=check, NP=NP, format=format, frontType=frontType, recomputeDist=recomputeDist,
    #                             expand=expand, tinit=tinit, initWithBBox=initWithBBox, wallAdapt=wallAdapt, dfarDir=dfarDir)
    # parallel prep
    #else:
    ret = prepare1(t_case, t_out, tc_out, snears=snears, dfar=dfar, dfarList=dfarList,
                   tbox=tbox, snearsf=snearsf, yplus=yplus,
                   vmin=vmin, check=check, NP=NP, format=format, frontType=frontType,
                   expand=expand, tinit=tinit, initWithBBox=initWithBBox,
                   wallAdapt=wallAdapt, dfarDir=dfarDir, redistribute=redistribute)

    return ret

#================================================================================
# IBM prepare - seq
#================================================================================
def prepare0(t_case, t_out, tc_out, snears=0.01, dfar=10., dfarList=[],
             tbox=None, snearsf=None, yplus=100.,
             vmin=21, check=False, NP=0, format='single', frontType=1, recomputeDist=False,
             expand=3, tinit=None, initWithBBox=-1., wallAdapt=None, dfarDir=0):
    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = t_case

    # list of dfars
    if dfarList == []:
        zones = Internal.getZones(tb)
        dfarList = [dfar*1.]*len(zones)
        for c, z in enumerate(zones):
            n = Internal.getNodeFromName2(z, 'dfar')
            if n is not None:
                dfarList[c] = Internal.getValue(n)*1.

    #-------------------------------------------------------
    # Refinement surfaces in the fluid
    #-------------------------------------------------------
    # snearsf: list of spacing required in the refinement surfaces
    if tbox is not None:
        if isinstance(tbox, str): tbox = C.convertFile2PyTree(tbox)
        else: tbox = tbox
        if snearsf is None:
            snearsf = []
            zones = Internal.getZones(tbox)
            for z in zones:
                sn = Internal.getNodeFromName2(z, 'snear')
                if sn is not None: snearsf.append(Internal.getValue(sn))
                else: snearsf.append(1.)

    #--------------------------------------------------------
    # Get Reference State and model from body pyTree
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('GoverningEquations is missing in input cgns.')
    # model: Euler, NSLaminar, NSTurbulent
    model = Internal.getValue(model)

    # check Euler non consistant avec Musker
    if model == 'Euler':
        for z in Internal.getZones(tb):
            ibctype = Internal.getNodeFromName2(z, 'ibctype')
            if ibctype is not None:
                ibctype = Internal.getValue(ibctype)
                if ibctype == 'Musker' or ibctype == 'Log':
                    raise ValueError("In tb: governing equations (Euler) not consistent with ibc type (%s)"%(ibctype))

    # reference state
    refstate = C.getState(tb)
    # dimension du pb
    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)
    if dimPb == 2: C._initVars(tb, 'CoordinateZ', 0.) # forced

    #--------------------------------------------------------
    # Generates the full Cartesian mesh
    t = G_IBM.generateIBMMesh_legacy(tb, vmin=vmin, snears=snears, dfar=dfar, dfarList=dfarList, DEPTH=2,
                                     tbox=tbox, snearsf=snearsf, check=check, sizeMax=1000000,
                                     expand=expand, dfarDir=dfarDir)
    test.printMem(">>> Build octree full [end]")

    #------------------------------------------------------
    # distribute the mesh over NP processors
    if NP > 0:
        print('distribution over %d processors'%NP)
        stats = D2._distribute(t, NP)
        if check: print(stats)

    #------------------------------------------------
    # Add reference state to the pyTree and init fields
    # Add viscosity if model is not Euler
    C._addState(t, state=refstate)
    C._addState(t, 'GoverningEquations', model)
    C._addState(t, 'EquationDimension', dimPb)
    # if check: C.convertPyTree2File(t, 'mesh1.cgns')

    #----------------------------------------
    # Computes distance field
    #----------------------------------------
    test.printMem(">>> wall distance [start]")
    if dimPb == 2:
        z0 = Internal.getZones(t)
        bb = G.bbox(z0); dz = bb[5]-bb[2]
        tb2 = C.initVars(tb, 'CoordinateZ', dz*0.5)
        DTW._distance2Walls(t, tb2, type='ortho', signed=0, dim=dimPb, loc='centers')
    else:
        DTW._distance2Walls(t, tb, type='ortho', signed=0, dim=dimPb, loc='centers')
    test.printMem(">>> wall distance [end]")

    #----------------------------------------
    # Create IBM info
    #----------------------------------------
    t,tc = X_IBM.prepareIBMData_legacy(t, tb, frontType=frontType, interpDataType=0, yplus=yplus, wallAdapt=wallAdapt)
    test.printMem(">>> ibm data [end]")

    # arbre donneur
    D2._copyDistribution(tc, t)
    if isinstance(tc_out, str): FastC.save(tc, tc_out, split=format, NP=-NP)

    #----------------------------------------
    # Extraction des coordonnees des pts IBM
    #----------------------------------------
    if check:
        tibm = P_IBM.extractIBMInfo(tc)
        C.convertPyTree2File(tibm, 'IBMInfo.cgns')
        del tibm

    #-----------------------------------------
    # Computes distance field for Musker only
    #-----------------------------------------
    if model != 'Euler' and recomputeDist:
        ibctypes = set()
        for node in Internal.getNodesFromName(tb,'ibctype'):
            ibctypes.add(Internal.getValue(node))
            if 'outpress' in ibctypes or 'inj' in ibctypes or 'slip' in ibctypes:
                test.printMem(">>> wall distance for viscous wall only [start]")
                for z in Internal.getZones(tb):
                    ibc = Internal.getNodeFromName(z,'ibctype')
                    if Internal.getValue(ibc)=='outpress' or Internal.getValue(ibc)=='inj' or Internal.getValue(ibc)=='slip':
                        Internal.rmNode(tb,z)

                if dimPb == 2:
                    z0 = Internal.getZones(t)
                    bb = G.bbox(z0); dz = bb[5]-bb[2]
                    tb2 = C.initVars(tb, 'CoordinateZ', dz*0.5)
                    DTW._distance2Walls(t,tb2,type='ortho', signed=0, dim=dimPb, loc='centers')
                else:
                    DTW._distance2Walls(t,tb,type='ortho', signed=0, dim=dimPb, loc='centers')
                test.printMem(">>> wall distance for viscous wall only [end]")

    # Initialisation
    if tinit is None:
        I._initConst(t, loc='centers')
        if model != "Euler": C._initVars(t, 'centers:ViscosityEddy', 0.)
    else:
        P._extractMesh(tinit, t, mode='accurate', constraint=40.)
        RefState = Internal.getNodeFromType(t, 'ReferenceState_t')
        ronutildeInf = Internal.getValue(Internal.getNodeFromName(RefState, 'TurbulentSANuTildeDensity'))
        vxInf = Internal.getValue(Internal.getNodeFromName(RefState, 'VelocityX'))
        vyInf = Internal.getValue(Internal.getNodeFromName(RefState, 'VelocityY'))
        vzInf = Internal.getValue(Internal.getNodeFromName(RefState, 'VelocityZ'))
        RhoInf = Internal.getValue(Internal.getNodeFromName(RefState, 'Density'))
        TInf = Internal.getValue(Internal.getNodeFromName(RefState, 'Temperature'))
        C._initVars(t,"{centers:VelocityX}=({centers:Density}<0.01)*%g+({centers:Density}>0.01)*{centers:VelocityX}"%vxInf)
        C._initVars(t,"{centers:VelocityY}=({centers:Density}<0.01)*%g+({centers:Density}>0.01)*{centers:VelocityY}"%vyInf)
        C._initVars(t,"{centers:VelocityZ}=({centers:Density}<0.01)*%g+({centers:Density}>0.01)*{centers:VelocityZ}"%vzInf)
        C._initVars(t,"{centers:Temperature}=({centers:Density}<0.01)*%g+({centers:Density}>0.01)*{centers:Temperature}"%TInf)
        C._initVars(t,"{centers:Density}=({centers:Density}<0.01)*%g+({centers:Density}>0.01)*{centers:Density}"%RhoInf)
        #C._initVars(t,"{centers:TurbulentSANuTildeDensity}=%g"%(ronutildeInf))

    # Init with BBox
    if initWithBBox>0.:
        print('initialisation par bounding box')
        bodybb = C.newPyTree(['Base'])
        for base in Internal.getBases(tb):
            bbox = G.bbox(base)
            bodybbz = D.box(tuple(bbox[:3]),tuple(bbox[3:]), N=2, ntype='STRUCT')
            Internal._append(bodybb,bodybbz,'Base')
        T._scale(bodybb, factor=(initWithBBox,initWithBBox,initWithBBox))
        tbb = G.BB(t)
        interDict = X.getIntersectingDomains(tbb,bodybb,taabb=tbb,taabb2=bodybb)
        for zone in Internal.getZones(t):
            zname = Internal.getName(zone)
            if interDict[zname] != []:
                C._initVars(zone, 'centers:MomentumX', 0.)
                C._initVars(zone, 'centers:MomentumY', 0.)
                C._initVars(zone, 'centers:MomentumZ', 0.)

    if isinstance(t_out, str): FastC.save(t, t_out, split=format, NP=-NP, cartesian=True)
    return t, tc


#================================================================================
# generateur grille Cart Ivan
#================================================================================
def generateCartesian(tb, dimPb=3, snears=0.01, dfar=10., dfarList=[], tbox=None, ext=3, snearsf=None, yplus=100.,
                      vmin=21, check=False, expand=3, dfarDir=0, to=None):

    rank = Cmpi.rank
    comm = Cmpi.COMM_WORLD
    refstate = C.getState(tb)
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('GoverningEquations is missing in input tree.')
    # model : Euler, NSLaminar, NSTurbulent
    model = Internal.getValue(model)


    # list of dfars
    if dfarList == []:
        zones = Internal.getZones(tb)
        dfarList = [dfar*1.]*len(zones)
        for c, z in enumerate(zones):
            n = Internal.getNodeFromName2(z, 'dfar')
            if n is not None: dfarList[c] = Internal.getValue(n)*1.

     # refinementSurfFile: surface meshes describing refinement zones
    if tbox is not None:
        if isinstance(tbox, str): tbox = C.convertFile2PyTree(tbox)
        else: tbox = tbox
        if snearsf is None:
            snearsf = []
            zones = Internal.getZones(tbox)
            for z in zones:
                sn = Internal.getNodeFromName2(z, 'snear')
                if sn is not None: snearsf.append(Internal.getValue(sn))
                else: snearsf.append(1.)
    symmetry = 0
    fileout = None
    if check: fileout = 'octree.cgns'
    # Octree identical on all procs
    test.printMem('>>> Octree unstruct [start]')
    if to is not None:
        if isinstance(to, str):
            o = C.convertFile2PyTree(to)
            o = Internal.getZones(o)[0]
        else:
            o = Internal.getZones(to)[0]
        parento = None
    else:
        o = G_IBM.buildOctree(tb, snears=snears, snearFactor=1., dfar=dfar, dfarList=dfarList,
                              to=to, tbox=tbox, snearsf=snearsf,
                              dimPb=dimPb, vmin=vmin, symmetry=symmetry, fileout=None, rank=rank,
                              expand=expand, dfarDir=dfarDir)

    if rank==0 and check: C.convertPyTree2File(o,fileout)
    # build parent octree 3 levels higher
    # returns a list of 4 octants of the parent octree in 2D and 8 in 3D
    parento = G_IBM.buildParentOctrees__(o, tb, snears=snears, snearFactor=4., dfar=dfar, dfarList=dfarList, to=to, tbox=tbox, snearsf=snearsf,
                                         dimPb=dimPb, vmin=vmin, symmetry=symmetry, fileout=None, rank=rank)
    test.printMem(">>> Octree unstruct [end]")

    # Split octree
    test.printMem(">>> Octree unstruct split [start]")
    bb = G.bbox(o)
    NPI = Cmpi.size
    if NPI == 1: p = Internal.copyRef(o) # keep reference
    else: p = T.splitNParts(o, N=NPI, recoverBC=False)[rank]
    del o
    test.printMem(">>> Octree unstruct split [end]")

    # fill vmin + merge in parallel
    test.printMem(">>> Octree struct [start]")
    res = G_IBM.octree2StructLoc__(p, vmin=vmin, ext=-1, optimized=0, parento=parento, sizeMax=1000000)
    del p
    if parento is not None:
        for po in parento: del po
    t = C.newPyTree(['CARTESIAN', res])
    zones = Internal.getZones(t)
    for z in zones: z[0] = z[0]+'X%d'%rank
    Cmpi._setProc(t, rank)

    C._addState(t, 'EquationDimension', dimPb)
    test.printMem(">>> Octree struct [end]")

    # Add xzones for ext
    test.printMem(">>> extended cart grids [start]")
    tbb = Cmpi.createBBoxTree(t)
    interDict = X.getIntersectingDomains(tbb)
    graph = Cmpi.computeGraph(tbb, type='bbox', intersectionsDict=interDict, reduction=False)
    del tbb
    Cmpi._addXZones(t, graph, variables=[], cartesian=True)
    test.printMem(">>> extended cart grids [after add XZones]")
    zones = Internal.getZones(t)
    coords = C.getFields(Internal.__GridCoordinates__, zones, api=3)
    coords, rinds = Generator.extendCartGrids(coords, ext=ext, optimized=1, extBnd=0)
    C.setFields(coords, zones, 'nodes')
    for noz in range(len(zones)):
        Internal.newRind(value=rinds[noz], parent=zones[noz])
    Cmpi._rmXZones(t)
    coords = None; zones = None
    test.printMem(">>> extended cart grids (after rmXZones) [end]")

    X_IBM._addBCOverlaps(t, bbox=bb)
    X_IBM._addExternalBCs(t, bbox=bb, dimPb=dimPb)

    if dimPb == 2:
        dz = 0.01
        T._addkplane(t)
        T._contract(t, (0,0,0), (1,0,0), (0,1,0), dz)

    # ReferenceState
    C._addState(t, state=refstate)
    C._addState(t, 'GoverningEquations', model)
    C._addState(t, 'EquationDimension', dimPb)
    return t

#================================================================================
# extrude mesh and case for z ou theta periodic cases
#
#================================================================================
def extrudeCartesian(t, tb, check=False, extrusion="cart", dz=0.01, NPas=10, span=1 , Ntranche=1, dict_Nz=None, ific=2,
                     isCartesianExtrude=False):

    rank = Cmpi.rank
    comm = Cmpi.COMM_WORLD

    if  extrusion == "cart": perio = span/float(Ntranche)
    else:                    perio = span/180.*math.pi/float(Ntranche)

    for z in Internal.getZones(t):

        cellN    = Internal.getNodeFromName(z, "cellN")[1]
        #cellNInit= Internal.getNodeFromName(z, "cellN#Init")[1]
        sh = numpy.shape(cellN)

        #modif CellN pour filtrer les cellule solide de l'interp chimere
        # et autoriser les ghost comme donneuse
        for k in range(sh[2]):
            for j in range(sh[1]):
                for i in range(sh[0]):
                    if  cellN[i,j,k] != 0:  cellN[i,j,k] =1
                    #if  cellNInit[i,j,k] != 0:  cellN[i,j,k] =1
                    #else: cellN[i,j,k] =0

    #Dim = 3D
    c = 0
    for tree in [t,tb]:
        for dim in Internal.getNodesFromName(tree,'EquationDimension'): dim[1]=3

        dz_loc={}; Nk={}
        if c ==0:
            for z in Internal.getZones(tree):
                name_zone = z[0]
                Nk_min = 1000000000
                if not isCartesianExtrude:
                    if dict_Nz is not None:
                        Nk[z[0]]     = int(dict_Nz[name_zone])
                        dz_loc[z[0]] = span/float(Ntranche*Nk[z[0]])
                    else:
                        Nk[z[0]]     = NPas-1
                        dz_loc[z[0]] = dz
                else:
                    if tree == t:
                        h = abs(C.getValue(z,'CoordinateX',0)-C.getValue(z,'CoordinateX',1))
                        NPas_local = int(round(span/h))
                        if NPas_local<2:
                            print("WARNING:: Zone %s has Nz=%d and is being clipped to Nz=2"%(z[0],NPas_local))
                            NPas_local=2
                        Nk[z[0]]     = NPas_local
                        dz_loc[z[0]] = span/float(Ntranche*Nk[z[0]])
                    else:
                        Nk[z[0]]     = NPas-1
                        dz_loc[z[0]] = dz

                if Nk[z[0]] < Nk_min: Nk_min =Nk[z[0]]
        else:
            for z in Internal.getZones(tree):
                Nk[z[0]]     =  Nk_min
                dz_loc[z[0]] = span/float(Ntranche*Nk_min)

        for z in Internal.getZones(tree):
            Nk[z[0]] +=2*ific-1  # -1, car un plan de maillage deja ajoute dans le cas 2D

        for z in Internal.getZones(tree):

            yy_2d   = Internal.getNodeFromName(z, "CoordinateY")[1]
            zz_2d   = Internal.getNodeFromName(z, "CoordinateZ")[1]

            sh_2d   = numpy.shape(yy_2d)

            T._addkplane(z   ,N=Nk[z[0]])

            zz   = Internal.getNodeFromName(z, "CoordinateZ")[1]
            yy   = Internal.getNodeFromName(z, "CoordinateY")[1]

            r       = numpy.sqrt( zz**2+yy**2)

            perio_loc = perio
            if c==1: period_loc= perio*1.5

            sh = numpy.shape(yy)
            if  len(sh_2d) ==1: #NGON
                if  extrusion == "cart":

                    for k in range( sh[1]):
                        zz[:,k] = (k-ific)* dz_loc[z[0]]
                else:
                    theta0  = numpy.arctan(zz_2d/yy_2d)
                    for k in range( sh[1]):
                        shift = perio_loc/float(sh[1]-5.)*(k-ific)
                        zz[:,k] = r[:,0]*numpy.sin(theta0[:] + shift)
                        yy[:,k] = r[:,0]*numpy.cos(theta0[:] + shift)
            else:
                if  extrusion == "cart":
                    for k in range(sh[2]):
                        zz[:,:,k] = (k-ific)* dz_loc[z[0]]
                else:
                    theta0  = numpy.arctan(zz_2d/yy_2d)
                    for k in range(sh[2]):
                        shift = perio_loc/float(sh[2]-5.)*(k-ific)
                        zz[:,:,k] = r[:,:,0]*numpy.sin(theta0[:,:,0] + shift)
                        yy[:,:,k] = r[:,:,0]*numpy.cos(theta0[:,:,0] + shift)

        if c==1: break
        c+=1

        for z in Internal.getZones(tree):

            zdim     = Internal.getValue(z)

            # Modif rind cellule GH en kmin et kmax
            rind = Internal.getNodeFromName1(z, "Rind")[1]
            rind[4]=ific; rind[5]=ific
            # Modif range BC
            BCs = Internal.getNodesFromType(z, "BC_t")
            for bc in BCs:
                ptrg = Internal.getNodeFromName(bc, "PointRange")[1]
                ptrg[2,0] = 3
                ptrg[2,1] = Nk[z[0]]

            # Creatioon connectivite perio dans t
            for idir in ['_kmax','_kmin']:
                if idir == '_kmax':
                    ktg  = zdim[2,0]
                    ktgD = 1
                    if extrusion =='cart':
                        angle = 0
                        trans = -perio
                    else:
                        angle = -perio
                        trans = 0.
                else:
                    ktgD  = zdim[2,0]
                    ktg   = 1
                    if extrusion =='cart':
                        angle = 0
                        trans = perio
                    else:
                        angle = perio
                        trans = 0.

                Conn = Internal.getNodeFromName(z, "ZoneGridConnectivity")
                name = 'match_'+z[0]+idir
                Internal.createUniqueChild(Conn, name, 'GridConnectivity1to1_t')
                tmp1     = Internal.getNodeFromName(Conn, name)
                Internal.setValue(tmp1, z[0])
                datap = numpy.empty( (3,2) , numpy.int32)
                datap[0,0]=1;datap[1,0]=1;datap[2,0]=ktg
                datap[0,1]=zdim[0,0];datap[1,1]=zdim[1,0];datap[2,1]= ktg
                Internal.createUniqueChild(tmp1 ,'PointRange', 'IndexRange_t',datap)
                datap = numpy.empty( (3,2) , numpy.int32)
                datap[0,0]=1;datap[1,0]=1;datap[2,0]=ktgD
                datap[0,1]=zdim[0,0];datap[1,1]=zdim[1,0];datap[2,1]= ktgD
                Internal.createUniqueChild(tmp1 ,'PointRangeDonor', 'IndexRange_t',datap)
                datap = numpy.empty( 3 , numpy.int32)
                datap[0]=1;datap[1]=2;datap[2]=3
                Internal.createUniqueChild(tmp1 ,'Transform', '"int[IndexDimension]"',datap)
                Internal.createUniqueChild(tmp1 ,'GridConnectivityProperty', 'GridConnectivityProperty_t')

                prop     = Internal.getNodeFromName(tmp1, 'GridConnectivityProperty')
                Internal.createUniqueChild(prop ,'Periodic', 'Periodic_t')
                period    = Internal.getNodeFromName(prop, 'Periodic')
                datap = numpy.zeros( 3 , numpy.float64)
                Internal.createUniqueChild(period ,'RotationCenter', 'DataArray_t',datap)
                datap = numpy.zeros( 3 , numpy.float64)
                datap[0]= angle
                Internal.createUniqueChild(period ,'RotationAngle' ,'DataArray_t',datap)
                datap = numpy.zeros( 3 , numpy.float64)
                datap[2]= trans
                Internal.createUniqueChild(period ,'Translation', 'DataArray_t',datap)

                rot     = Internal.getNodeFromName(period, 'RotationAngle')
                Units=['Kilogram','Meter','Second','Kelvin','Radian']
                Internal.createUniqueChild(rot ,'DimensionalUnits' ,'DimensionalUnits_t',Units)

    return t, tb

#================================================================================
# IBM prepare hybryde cart + curvi
# in : t_3d arbre cartesien   (issue de prepare1)
# in : tc_3d arbre cartesien connectivite   (issue de prepare1)
# in : t_curvi arbre curviligne   (avec bc et rac, mais sans ghost)
#==================================================
def setInterpData_Hybride(t_3d, tc_3d, t_curvi, extrusion=None, interpDataType=1):

    overlap     ='14'
    InterpOrder =2
    root_racHybride='joinIBC'

    #ajout 2 ghost maillag curvi
    Internal._addGhostCells(t_curvi, t_curvi, 2, adaptBCs=1, fillCorner=0)
    C._initVars(t_curvi,"centers:cellN",1.)
    t_curvi = TBX.modifyBCOverlapsForGhostMesh(t_curvi,2)

    dimPb = Internal.getValue(Internal.getNodeFromName(t_3d, 'EquationDimension'))

    for z in Internal.getZones(t_curvi):
        if root_racHybride in z[0]:
            C._fillEmptyBCWith(z,'inactive','BCExtrapolated',dim=dimPb)
            C._rmBCOfType(z,'BCOverlap')
            C._fillEmptyBCWith(z,'trou','BCOverlap',dim=dimPb)
            X._applyBCOverlaps(z,loc='centers',depth=2, val=0)
            X._applyBCOverlaps(z,loc='centers',depth=4, val=2)

    if extrusion =='cyl':
        T._cart2Cyl(t_3d , (0,0,0),(1,0,0))
        T._cart2Cyl(tc_3d, (0,0,0),(1,0,0))
        T._cart2Cyl(t_curvi   , (0,0,0),(1,0,0))

    tBB2 = G.BB(t_curvi)
    tBB  = G.BB(t_3d)

    tc2 = C.node2Center(t_curvi)

    X._setInterpData(t_curvi, tc2, nature=1, loc='centers', storage='inverse', sameName=1, dim=dimPb, itype='abutting')

    C._initVars(t_3d,"{centers:cellN#Init}={centers:cellN}")

    for var in ['wall','racChimer']:
        C._initVars(t_3d,"centers:cellN",1.)
        if var == 'wall': itype ="3"
        else: itype = overlap
        for zc in Internal.getZones(tc_3d):
            for zsr in Internal.getNodesFromType(zc, "ZoneSubRegion_t"):
                zsrname = Internal.getName(zsr)
                zsrname = zsrname.split('_')
                if zsrname[0]=='IBCD' and zsrname[1] == itype:
                    zrname = Internal.getValue(zsr)
                    ptlistD= Internal.getNodeFromName(zsr,'PointListDonor')[1]
                    zloc   = Internal.getNodeFromName(t_3d,zrname)
                    sol    = Internal.getNodeFromName(zloc,'FlowSolution#Centers')
                    cellN  = Internal.getNodeFromName(sol,'cellN')[1]
                    sh     = numpy.shape(cellN)
                    ni= sh[0]; ninj = sh[0]*sh[1]
                    for l0 in range(numpy.size(ptlistD)):
                        l = ptlistD[ l0]
                        k = l//ninj
                        j = (l-k*ninj)//ni
                        i =  l-k*ninj - j*ni
                        cellN[ i,j,k ]=2

        if var == 'wall': C._initVars(t_3d,"{centers:cellN#wall}={centers:cellN}")
        else:             C._initVars(t_3d,"{centers:cellN#racChim}={centers:cellN}")

    Internal._rmNodesFromName(tc_3d,'IBCD_'+overlap+'_*')

    for z in Internal.getZones(t_3d):
        sol            = Internal.getNodeFromName(z,'FlowSolution#Centers')
        cellNRac       = Internal.getNodeFromName(sol,'cellN#racChim')[1]
        cellN          = Internal.getNodeFromName(sol,'cellN')[1]
        C._initVars(z,"{centers:cellN}={centers:cellN#racChim}")
        sh             = numpy.shape(cellN)

        if sh[2]> 1:
            for k in [0,1, sh[2]-2, sh[2]-1]:
                for j in range(sh[1]):
                    for i in range(sh[0]):
                        if  cellN[i,j,k] != 0:  cellN[i,j,k] =1

    interDict = X.getIntersectingDomains(tBB, t2=tBB2, method='AABB', taabb=tBB, taabb2=tBB2)
    print(" Interp Cartesian from curvilinear")
    zonelist=[]
    for z in Internal.getZones(t_3d):
        if C.getMaxValue(z,'centers:cellN')==2:
            dnrZones = []
            for zdname in interDict[z[0]]:
                zd = Internal.getNodeFromName(tc2,zdname)
                dnrZones.append(zd)
            X._setInterpData(z,dnrZones, nature=0,penalty=1,loc='centers',storage='inverse',sameName=1,\
                             interpDataType=interpDataType, itype='chimera', order=InterpOrder)
            z = X.getOversetInfo(z, dnrZones, loc='center',type='extrapolated')
            zonelist.append(z)

    interDict_curvi = X.getIntersectingDomains(tBB2, t2=tBB, method='AABB', taabb=tBB2, taabb2=tBB)
    print(" Interp curvi from Cartesian")
    for z in Internal.getZones(t_curvi):
        if C.getMaxValue(z,'centers:cellN')==2:
            #if 'join' in z[0]:
            dnrZones = []
            for zdname in interDict_curvi[z[0]]:
                zd = Internal.getNodeFromName(tc_3d,zdname)
                dnrZones.append(zd)

            X._setInterpData(z,dnrZones,nature=0,penalty=1,loc='centers',storage='inverse',sameName=1,\
                             interpDataType=interpDataType, itype='chimera', order=InterpOrder)

            # to check orphans
            #z = X.getOversetInfo(z, dnrZones, loc='center',type='extrapolated')
            #z = X.getOversetInfo(z, dnrZones, loc='center',type='orphan')
            #zonelist.append(z)
            #C.convertPyTree2File(z,"rcv2.cgns")


    C._initVars(t_3d,'{centers:cellN}={centers:cellN#Init}')

    if extrusion =='cyl':
        T._cyl2Cart(t_3d,   (0,0,0),(1,0,0))
        T._cyl2Cart(t_curvi,  (0,0,0),(1,0,0))
        T._cyl2Cart(tc_3d,  (0,0,0),(1,0,0))
        T._cyl2Cart(tc2, (0,0,0),(1,0,0))

    for z in Internal.getZones(t_curvi):
        for bc in  Internal.getNodesFromType(z,'BC_t'):
            if 'inactive' in bc[0] and Internal.getValue(bc) == 'BCExtrapolated':
                Internal._rmNodesByName(z, bc[0])

    t  = C.mergeTrees(t_3d ,t_curvi )
    tc = C.mergeTrees(tc_3d,tc2)

    return t, tc

#================================================================================
# IBM prepare - para IM
#
# extrusion: make an extrusion from a 2D profile. ATTENTION, each zone of the profile must be joined in one single zone
# smoothing : smooth the front during the front 2 specific treatment in the cases of local refinements
# balancing ; balance the entire distribution after the octree generation, useful for symetries
# distrib : new distribution at the end of prepare1
#===================================================================================================================
def prepare1(t_case, t_out, tc_out, t_in=None, to=None, snears=0.01, dfar=10., dfarList=[],
             tbox=None, snearsf=None, yplus=100., Lref=1.,
             vmin=21, check=False, NP=0, format='single', interpDataType=0, order=2, ext=2,
             frontType=1, extrusion=None, smoothing=False, balancing=False, recomputeDist=False,
             distrib=True, expand=3, tinit=None, initWithBBox=-1., wallAdapt=None, yplusAdapt=100., dfarDir=0,
             correctionMultiCorpsF42=False, blankingF42=False, twoFronts=False, redistribute=False, IBCType=1,
             height_in=1.0,isoverideheight=False,isFilamentOnly=False,closedSolid=[],isWireModel=False, cleanCellN=True):

    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = Internal.copyTree(t_case)

    isOrtho_Project_First = isFilamentOnly
    if not closedSolid and not isFilamentOnly:
        closedSolid=[]
        for b in Internal.getBases(tb):
            closedSolid.append(b[0])

    tb_filament_list=[]
    if not isFilamentOnly:
        tb_filament     =Internal.copyTree(tb)
        for b in Internal.getBases(tb):
            if b[0] not in closedSolid:
                tb_filament_list.append(b[0])

        if tb_filament_list:
            list_delete=[]
            for b in Internal.getBases(tb):
                if b[0] in tb_filament_list:
                    list_delete.append(b)
            for b in list_delete:
                Internal._rmNode(tb,b)

            list_delete=[]
            for b in Internal.getBases(tb_filament):
                if b[0] not in tb_filament_list:
                    list_delete.append(b)
            for b in list_delete:
                Internal._rmNode(tb_filament,b)
            isOrtho_Project_First = True

    rank = Cmpi.rank
    comm = Cmpi.COMM_WORLD

    DEPTH=2

    cartesian = True
    if extrusion is not None: cartesian = False

    # reference state
    refstate = C.getState(tb)
    # dimension du pb
    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)

    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('GoverningEquations is missing in input tree.')
    # model : Euler, NSLaminar, NSTurbulent
    model = Internal.getValue(model)

    # check Euler non consistant avec Musker
    if model == 'Euler':
        for z in Internal.getZones(tb):
            ibctype = Internal.getNodeFromName2(z, 'ibctype')
            if ibctype is not None:
                ibctype = Internal.getValue(ibctype)
                if ibctype == 'Musker' or ibctype == 'Log':
                    raise ValueError("In tb: governing equations (Euler) not consistent with ibc type (%s)"%(ibctype))

    if dimPb == 2 and cleanCellN == False: C._initVars(tb, 'CoordinateZ', 0.) # forced
    if t_in is None:
        t = generateCartesian(tb, dimPb=dimPb, snears=snears, dfar=dfar, dfarList=dfarList, tbox=tbox, ext=ext+1,
                              snearsf=snearsf, yplus=yplus,vmin=vmin, check=check, expand=expand, dfarDir=dfarDir)
    else:
        t = t_in

    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)
    # Balancing
    if balancing:
        test.printMem(">>> balancing [start]")
        Cmpi.convertPyTree2File(t, t_out)
        # Need to wait for all the procs to write their parts before the new distribution
        comm.Barrier()
        ts = Cmpi.convertFile2SkeletonTree(t_out)
        D2._distribute(ts, Cmpi.size, algorithm='graph')
        t = Cmpi.readZones(ts, t_out, rank=rank)
        Cmpi._convert2PartialTree(t)
        zones = Internal.getZones(t)
        for z in zones: z[0] = z[0] + 'X%d'%rank
        del ts
        test.printMem(">>> balancing [end]")

    if extrusion == 'cyl':
        T._cart2Cyl(t, (0,0,0),(1,0,0))
        T._cart2Cyl(tb, (0,0,0),(1,0,0))

    # Distance a la paroi
    test.printMem(">>> Wall distance [start]")
    FSC = Internal.getNodeFromType(t,"FlowSolution_t")
    if (FSC is None or Internal.getNodeFromName(FSC,'TurbulentDistance') is None) and extrusion is None:
        C._initVars(t,'{centers:TurbulentDistance}=1e06')
        if dimPb == 2:
            z0 = Internal.getNodeFromType2(t, "Zone_t")
            bb0 = G.bbox(z0); dz = bb0[5]-bb0[2]

            tb2      = C.initVars(tb, 'CoordinateZ', dz*0.5)
            tbsave   = tb2
            t        = dist2wallNearBody(t, tb2, type='ortho', signed=0, dim=dimPb, loc='centers')

            if tb_filament_list:
                C._initVars(t,'{centers:TurbulentDistanceSolid}={centers:TurbulentDistance}')
                C._initVars(t,'{centers:TurbulentDistance}=1e06')

                tb2      = C.initVars(tb_filament, 'CoordinateZ', dz*0.5)
                tbsave   = tb2
                t        = dist2wallNearBody(t, tb2, type='ortho', signed=0, dim=dimPb, loc='centers')
                C._initVars(t,'{centers:TurbulentDistanceFilament}={centers:TurbulentDistance}')
                C._initVars(t,'{centers:TurbulentDistance}=minimum({centers:TurbulentDistanceSolid},{centers:TurbulentDistanceFilament})')

                C._initVars(t,"{centers:TurbulentDistanceSolid}=({centers:TurbulentDistanceSolid}>1e03)*0+({centers:TurbulentDistanceSolid}<1e03)*{centers:TurbulentDistanceSolid}")
                C._initVars(t,"{centers:TurbulentDistanceFilament}=({centers:TurbulentDistanceFilament}>1e03)*0+({centers:TurbulentDistanceFilament}<1e03)*{centers:TurbulentDistanceFilament}")

        else:
            t=dist2wallNearBody(t, tb, type='ortho', signed=0, dim=dimPb, loc='centers')

            if tb_filament_list:
                C._initVars(t,'{centers:TurbulentDistanceSolid}={centers:TurbulentDistance}')
                C._initVars(t,'{centers:TurbulentDistance}=1e06')

                t=dist2wallNearBody(t, tb_filament, type='ortho', signed=0, dim=dimPb, loc='centers')

                C._initVars(t,'{centers:TurbulentDistanceFilament}={centers:TurbulentDistance}')
                C._initVars(t,'{centers:TurbulentDistance}=minimum({centers:TurbulentDistanceSolid},{centers:TurbulentDistanceFilament})')

                C._initVars(t,"{centers:TurbulentDistanceSolid}=({centers:TurbulentDistanceSolid}>1e03)*0+({centers:TurbulentDistanceSolid}<1e03)*{centers:TurbulentDistanceSolid}")
                C._initVars(t,"{centers:TurbulentDistanceFilament}=({centers:TurbulentDistanceFilament}>1e03)*0+({centers:TurbulentDistanceFilament}<1e03)*{centers:TurbulentDistanceFilament}")

        C._initVars(t,"{centers:TurbulentDistance}=({centers:TurbulentDistance}>1e03)*0+({centers:TurbulentDistance}<1e03)*{centers:TurbulentDistance}")

        # Compute turbulentdistance wrt each body that is not a sym plan
        if correctionMultiCorpsF42 and frontType==42:
            test.printMem(">>> Individual wall distance [start]")
            # Keep track of the general turbulentDistance
            C._initVars(t,'{centers:TurbulentDistance_ori}={centers:TurbulentDistance}')

            Reynolds = Internal.getNodeFromName(tb, 'Reynolds')
            if Reynolds is not None:
                Reynolds = Internal.getValue(Reynolds)
            else:
                Reynolds = 6.e6

            if yplus > 0.:
                shiftDist = G_IBM_Height.computeModelisationHeight(Re=Reynolds, yplus=yplus, L=Lref)
            else:
                snears = Internal.getNodesFromName(tb, 'snear')
                h = max(snears, key=lambda x: x[1])[1]
                shiftDist = G_IBM_Height.computeBestModelisationHeight(Re=Reynolds, h=h) # meilleur compromis entre hauteur entre le snear et la hauteur de modelisation

            for z in Internal.getZones(t):
                cptBody = 1
                if dimPb == 3: tb2 = tb
                for body in Internal.getNodesFromType(tb2,'Zone_t'):
                    if body[0] != "sym" and ("closure" not in body[0]):
                        # Create extanded BBox around each body
                        bboxBody = G.BB(body)
                        coordX = Internal.getNodeFromName(bboxBody, 'CoordinateX')[1]
                        coordX[0] = coordX[0] - shiftDist
                        coordX[1] = coordX[1] + shiftDist
                        Internal.getNodeFromName(bboxBody, 'CoordinateX')[1] = coordX
                        coordY = Internal.getNodeFromName(bboxBody, 'CoordinateY')[1]
                        coordY[0][0] = coordY[0][0] - shiftDist
                        coordY[1][0] = coordY[1][0] - shiftDist
                        coordY[0][1] = coordY[0][1] + shiftDist
                        coordY[1][1] = coordY[1][1] + shiftDist
                        Internal.getNodeFromName(bboxBody, 'CoordinateY')[1] = coordY
                        coordZ = Internal.getNodeFromName(bboxBody, 'CoordinateZ')[1]
                        coordZ[0][0][0] = coordZ[0][0][0] - shiftDist
                        coordZ[0][1][0] = coordZ[0][1][0] - shiftDist
                        coordZ[1][0][0] = coordZ[1][0][0] - shiftDist
                        coordZ[1][1][0] = coordZ[1][1][0] - shiftDist
                        coordZ[0][0][1] = coordZ[0][0][1] + shiftDist
                        coordZ[0][1][1] = coordZ[0][1][1] + shiftDist
                        coordZ[1][0][1] = coordZ[1][0][1] + shiftDist
                        coordZ[1][1][1] = coordZ[1][1][1] + shiftDist
                        Internal.getNodeFromName(bboxBody, 'CoordinateZ')[1] = coordZ
                        bboxZone = G.BB(z)

                        # Compute new individual turbulentDistance when blocks are close enough
                        if G.bboxIntersection(bboxBody, bboxZone, isBB=True):
                            DTW._distance2Walls(z, body, type='ortho', signed=0, dim=dimPb, loc='centers')
                            C._initVars(z,'{centers:TurbulentDistance_body%i={centers:TurbulentDistance}'%cptBody)
                        else:
                            C._initVars(z,'{centers:TurbulentDistance_body%i=1000'%cptBody)
                        cptBody += 1
                if dimPb == 3: del tb2

            C._initVars(t,'{centers:TurbulentDistance}={centers:TurbulentDistance_ori}')
            C._rmVars(t,['centers:TurbulentDistance_ori'])

        #
        if dimPb == 2 and cleanCellN == False : C._initVars(t, '{centers:TurbulentDistanceAllBC}={centers:TurbulentDistance}')

    else:
        C._initVars(t, '{centers:TurbulentDistance}={centers:TurbulentDistanceAllBC}')

    test.printMem(">>> Wall distance [end]")

    X._applyBCOverlaps(t, depth=DEPTH, loc='centers', val=2, cellNName='cellN')

    # Blank des corps chimere
    # applyBCOverlap des maillages de corps
    # SetHoleInterpolated points

    C._initVars(t,'{centers:cellNChim}={centers:cellN}')
    C._initVars(t, 'centers:cellN', 1.)
    if dimPb == 2:
        z0 = Internal.getNodeFromType2(t, 'Zone_t')
        dims = Internal.getZoneDim(z0)
        npts = dims[1]*dims[2]*dims[3]
        zmin = C.getValue(z0,'CoordinateZ',0)
        zmax = C.getValue(z0,'CoordinateZ',npts-1)
        dz = zmax-zmin
        # Creation du corps 2D pour le preprocessing IBC
        T._addkplane(tb)
        T._contract(tb, (0,0,0), (1,0,0), (0,1,0), dz)
    test.printMem(">>> Blanking [start]")

    if extrusion is None:
        test.printMem(">>> Blanking [start]")
        if not isFilamentOnly:t = X_IBM.blankByIBCBodies(t, tb, 'centers', dimPb,closedSolid=closedSolid)
        if dimPb == 2 and cleanCellN == False: C._initVars(t, '{centers:cellNIBC_blank}={centers:cellN}')
    else:
        C._initVars(t, '{centers:cellN}={centers:cellNIBC_blank}')

    C._initVars(t, '{centers:cellNIBC}={centers:cellN}')


    if not isFilamentOnly:X_IBM._signDistance(t)

    if extrusion is not None:
        C._initVars(t,'{centers:cellN}={centers:cellNIBC_blank}')
    else:
        C._initVars(t,'{centers:cellN}={centers:cellNIBC}')

    if tb_filament_list or isFilamentOnly:
        if isFilamentOnly:
            C._initVars(t,'{centers:TurbulentDistanceFilament}={centers:TurbulentDistance}')
            maxy=C.getMaxValue(tb, ['CoordinateY']);
            miny=C.getMinValue(tb, ['CoordinateY']);
        if tb_filament_list:
            maxy=C.getMaxValue(tb_filament_list, ['CoordinateY']);
            miny=C.getMinValue(tb_filament_list, ['CoordinateY']);
        for z in Internal.getZones(t):
            if C.getMaxValue(z, 'centers:TurbulentDistanceFilament')> 1e-05:
                sol  = Internal.getNodeByName(z,"FlowSolution#Centers")
                cellN= Internal.getNodeByName(sol,'cellN')[1]
                dist = Internal.getNodeByName(sol,'TurbulentDistanceFilament')[1]
                ycord= Internal.getNodeByName(z,'CoordinateY')[1]
                h = abs(C.getValue(z,'CoordinateX',0)-C.getValue(z,'CoordinateX',1))
                sh=numpy.shape(dist)
                for k in range(sh[2]):
                    for j in range(sh[1]):
                        for i in range(sh[0]):
                            valy=0.5*(ycord[i,j,k]+ycord[i,j+1,k])
                            if dist[i,j]<numpy.sqrt(8)*h and valy<maxy and valy>miny :
                                cellN[i,j]=2
        C._rmVars(t,['centers:TurbulentDistanceFilament'])
        if tb_filament_list:C._rmVars(t,['centers:TurbulentDistanceSolid'])

    # determination des pts IBC
    Reynolds = Internal.getNodeFromName(tb, 'Reynolds')
    if Reynolds is not None: Reynolds = Internal.getValue(Reynolds)
    if Reynolds < 1.e5: frontType = 1
    if extrusion is None:
        if frontType != 42 :
            if IBCType == -1: X._setHoleInterpolatedPoints(t,depth=-DEPTH,dir=0,loc='centers',cellNName='cellN',addGC=False)
            elif IBCType == 1:
                X._setHoleInterpolatedPoints(t,depth=1,dir=1,loc='centers',cellNName='cellN',addGC=False) # pour les gradients
                if frontType < 2:
                    X._setHoleInterpolatedPoints(t,depth=DEPTH,dir=0,loc='centers',cellNName='cellN',addGC=False)
                else:
                    DEPTHL=DEPTH+1
                    X._setHoleInterpolatedPoints(t,depth=DEPTHL,dir=0,loc='centers',cellNName='cellN',addGC=False)
                    #cree des pts extrapoles supplementaires
                    # _blankClosestTargetCells(t,cellNName='cellN', depth=DEPTHL)
            else:
                raise ValueError('prepareIBMData: not valid IBCType. Check model.')
        else:
            # F42: tracking of IB points using distance information
            # the classical algorithm (front 1) is first used to ensure a minimum of two rows of target points around the geometry
            C._initVars(t,'{centers:cellNMin}={centers:cellNIBC}')
            if IBCType == -1: X._setHoleInterpolatedPoints(t,depth=-DEPTH,dir=0,loc='centers',cellNName='cellNMin',addGC=False)
            elif IBCType == 1: X._setHoleInterpolatedPoints(t,depth=1,dir=1,loc='centers',cellNName='cellNMin',addGC=False) # pour les gradients
            X._setHoleInterpolatedPoints(t,depth=DEPTH,dir=0,loc='centers',cellNName='cellNMin',addGC=False)

            for z in Internal.getZones(t):
                h = abs(C.getValue(z,'CoordinateX',0)-C.getValue(z,'CoordinateX',1))
                if yplus > 0.:
                    if isoverideheight:height = height_in
                    else:height = G_IBM_Height.computeModelisationHeight(Re=Reynolds, yplus=yplus, L=Lref)
                else:
                    if isoverideheight:height = height_in
                    else:height = G_IBM_Height.computeBestModelisationHeight(Re=Reynolds, h=h) # meilleur compromis entre hauteur entre le snear et la hauteur de modelisation
                    yplus  = G_IBM_Height.computeYplus(Re=Reynolds, height=height, L=Lref)
                C._initVars(z,'{centers:cellN}=({centers:TurbulentDistance}>%20.16g)+(2*({centers:TurbulentDistance}<=%20.16g)*({centers:TurbulentDistance}>0))'%(height,height))

                if correctionMultiCorpsF42:
                    # Prevent different body modeling from overlapping -> good projection of image points in the wall normal direction

                    epsilon_dist = 2*(abs(C.getValue(z,'CoordinateX',1)-C.getValue(z,'CoordinateX',0)))
                    max_dist = 2*0.1*Lref

                    # Try to find the best route between two adjacent bodies by finding optimal iso distances
                    def correctionMultiCorps(cellN, cellNF):
                        if cellN == 2 and cellNF == 2:
                            return 1
                        return cellN

                    def findIsoFront(cellNFront, Dist_1, Dist_2):
                        if Dist_1 < max_dist and Dist_2 < max_dist:
                            if abs(Dist_1-Dist_2) < epsilon_dist:
                                return 2
                        return max(cellNFront,1)

                    for i in range(1, cptBody):
                        for j in range(1, cptBody):
                            if j != i:
                                C._initVars(z,'centers:cellNFrontIso', findIsoFront, ['centers:cellNFrontIso', 'centers:TurbulentDistance_body%i'%i, 'centers:TurbulentDistance_body%i'%j])

                    C._initVars(z,'centers:cellN', correctionMultiCorps, ['centers:cellN', 'centers:cellNFrontIso'])

                    for i in range(1,cptBody):
                        C._rmVars(z,['centers:cellN_body%i'%i, 'centers:TurbulentDistance_body%i'%i])

            if wallAdapt is not None:
                # Use previous computation to adapt the positioning of IB points around the geometry (impose y+PC <= y+ref)
                # Warning: the wallAdapt file has to be obtained with TIBM.createWallAdapt(tc)
                C._initVars(t,'{centers:yplus}=100000.')
                w = C.convertFile2PyTree(wallAdapt)
                total = len(Internal.getZones(t))
                cpt = 1
                for z in Internal.getZones(t):
                    print("{} / {}".format(cpt, total))
                    cellN = Internal.getNodeFromName(z,'cellN')[1]
                    if 2 in cellN:
                        zname = z[0]
                        zd = Internal.getNodeFromName(w, zname)
                        if zd is not None:
                            yplus_w = Internal.getNodeFromName(zd, 'yplus')[1]
                            listIndices = Internal.getNodeFromName(zd, 'PointListDonor')[1]

                            n = numpy.shape(yplus_w)[0]
                            yplusA = Converter.array('yplus', n, 1, 1)
                            yplusA[1][:] = yplus_w

                            C._setPartialFields(z, [yplusA], [listIndices], loc='centers')

                    cpt += 1

                C._initVars(t,'{centers:cellN}=({centers:cellN}>0) * ( (({centers:cellN}) * ({centers:yplus}<=%20.16g)) + ({centers:yplus}>%20.16g) )'%(yplus,yplus))

            # final security gate, we ensure that we have at least to layers of target points
            C._initVars(t, '{centers:cellN} = maximum({centers:cellN}, {centers:cellNMin})')
            C._rmVars(t,['centers:yplus', 'centers:cellNMin'])

            # propagate max yplus between procs
            yplus = numpy.array([float(yplus)])
            yplus_max = numpy.zeros(1)
            comm.Allreduce(yplus, yplus_max, MPI.MAX)
            yplus = int(yplus_max[0])

            # Only keep the layer of target points useful for solver iterations, particularly useful in 3D
            if blankingF42: X._maximizeBlankedCells(t, depth=2, addGC=False)

        if dimPb == 2 and cleanCellN == False: C._initVars(t, '{centers:cellNIBC_hole}={centers:cellN}')

    else:  # extrusion
        C._initVars(t, '{centers:cellN}={centers:cellNIBC_hole}')

    if extrusion is None:
        if not isFilamentOnly:X_IBM._removeBlankedGrids(t, loc='centers')

    test.printMem(">>> Blanking [end]")

    print('Nb of Cartesian grids=%d.'%len(Internal.getZones(t)))
    npts = 0
    for i in Internal.getZones(t):
        dims = Internal.getZoneDim(i)
        npts += dims[1]*dims[2]*dims[3]
    print('Final number of points=%5.4f millions.'%(npts/1000000.))

    C._initVars(t, '{centers:cellNIBC}={centers:cellN}')

    if IBCType==-1:
        #print('Points IBC interieurs: on repousse le front un peu plus loin.')
        C._initVars(t,'{centers:cellNDummy}=({centers:cellNIBC}>0.5)*({centers:cellNIBC}<1.5)')
        X._setHoleInterpolatedPoints(t,depth=1,dir=1,loc='centers',cellNName='cellNDummy',addGC=False)
        C._initVars(t,'{centers:cellNFront}=logical_and({centers:cellNDummy}>0.5, {centers:cellNDummy}<1.5)')
        C._rmVars(t, ['centers:cellNDummy'])
        for z in Internal.getZones(t):
            connector._updateNatureForIBM(z, IBCType,
                                          Internal.__GridCoordinates__,
                                          Internal.__FlowSolutionNodes__,
                                          Internal.__FlowSolutionCenters__)
    else:
        C._initVars(t,'{centers:cellNFront}=logical_and({centers:cellNIBC}>0.5, {centers:cellNIBC}<1.5)')

        for z in Internal.getZones(t):
            if twoFronts:
                epsilon_dist = abs(C.getValue(z,'CoordinateX',1)-C.getValue(z,'CoordinateX',0))
                dmin = math.sqrt(3)*4*epsilon_dist
                if frontType == 42:
                    SHIFTB = G_IBM_Height.computeModelisationHeight(Re=Reynolds, yplus=yplus, L=Lref)
                    dmin = max(dmin, SHIFTB+math.sqrt(3)*2*epsilon_dist) # where shiftb = hmod
                C._initVars(z,'{centers:cellNIBC_2}=({centers:TurbulentDistance}>%20.16g)+(2*({centers:TurbulentDistance}<=%20.16g)*({centers:TurbulentDistance}>0))'%(dmin,dmin))
                C._initVars(z,'{centers:cellNFront_2}=logical_and({centers:cellNIBC_2}>0.5, {centers:cellNIBC_2}<1.5)')

            connector._updateNatureForIBM(z, IBCType,
                                          Internal.__GridCoordinates__,
                                          Internal.__FlowSolutionNodes__,
                                          Internal.__FlowSolutionCenters__)
    ##[AJ] - begin
    ##Ghost kmin et kmax donneuse potentiel
    if extrusion is not None:
        listvars_local =['cellN','cellNChim','cellNIBC','cellNFront']
        for z in Internal.getZones(t):
            sol            = Internal.getNodeFromName(z,'FlowSolution#Centers')
            for var in listvars_local:
                cellN          = Internal.getNodeFromName(sol,var)[1]
                sh             = numpy.shape(cellN)
                if var== 'cellNChim' or var== 'cellNIBC':
                    for k in [0,1, sh[2]-2, sh[2]-1]:
                        for j in range(sh[1]):
                            for i in range(sh[0]):
                                if  cellN[i,j,k] != 0:  cellN[i,j,k] =1
    ##[AJ] - end

    # setInterpData - Chimere
    C._initVars(t,'{centers:cellN}=maximum(0.,{centers:cellNChim})')# vaut -3, 0, 1, 2 initialement

    # maillage donneur: on MET les pts IBC comme donneurs
    # tp = Internal.copyRef(t)
    # FSN = Internal.getNodesFromName3(tp, Internal.__FlowSolutionNodes__)
    # Internal._rmNodesByName(FSN, 'cellNFront')
    # Internal._rmNodesByName(FSN, 'cellNIBC')
    # Internal._rmNodesByName(FSN, 'TurbulentDistance')
    # tc = C.node2Center(tp); del tp

    test.printMem(">>> Interpdata [start]")
    tc = C.node2Center(t)

    # abutting ?
    if Internal.getNodeFromType(t,"GridConnectivity1to1_t") is not None:
        test.printMem("setInterpData abutting")
        Xmpi._setInterpData(t, tc, nature=1, loc='centers', storage='inverse', sameName=1, dim=3, itype='abutting')
        test.printMem("setInterpData abutting done.")

    # setInterpData parallel pour le chimere
    tbbc = Cmpi.createBBoxTree(tc)
    interDict = X.getIntersectingDomains(tbbc)
    graph = Cmpi.computeGraph(tbbc, type='bbox', intersectionsDict=interDict, reduction=False)
    Cmpi._addXZones(tc, graph, variables=['cellN'], cartesian=cartesian)
    test.printMem(">>> Interpdata [after addXZones]")

    procDict = Cmpi.getProcDict(tc)
    datas = {}
    for zrcv in Internal.getZones(t):
        zrname = zrcv[0]
        dnrZones = []
        for zdname in interDict[zrname]:
            zd = Internal.getNodeFromName2(tc, zdname)
            dnrZones.append(zd)
        X._setInterpData(zrcv, dnrZones, nature=1, penalty=1, loc='centers', storage='inverse',
                         sameName=1, interpDataType=interpDataType, order=order, itype='chimera')
        for zd in dnrZones:
            zdname = zd[0]
            destProc = procDict[zdname]

            #allIDs = Internal.getNodesFromName(zd, 'ID*')
            #IDs = []
            #for zsr in allIDs:
            #    if Internal.getValue(zsr)==zrname: IDs.append(zsr)
            IDs = []
            for i in zd[2]:
                if i[0][0:2] == 'ID':
                    if Internal.getValue(i)==zrname: IDs.append(i)

            if IDs != []:
                if destProc == rank:
                    zD = Internal.getNodeFromName2(tc, zdname)
                    zD[2] += IDs
                else:
                    if destProc not in datas: datas[destProc] = [[zdname,IDs]]
                    else: datas[destProc].append([zdname,IDs])
            else:
                if destProc not in datas: datas[destProc] = []
    Cmpi._rmXZones(tc)
    test.printMem(">>> Interpdata [after rmXZones]")
    destDatas = Cmpi.sendRecv(datas, graph)
    for i in destDatas:
        for n in destDatas[i]:
            zname = n[0]
            IDs = n[1]
            if IDs != []:
                zD = Internal.getNodeFromName2(tc, zname)
                zD[2] += IDs
    datas = {}; destDatas = None; graph={}
    test.printMem(">>> Interpdata [after free]")
    test.printMem(">>> Interpdata [end]")

    # fin interpData
    C._initVars(t,'{centers:cellNIBCDnr}=minimum(2.,abs({centers:cellNIBC}))')
    C._initVars(t,'{centers:cellNIBC}=maximum(0.,{centers:cellNIBC})')# vaut -3, 0, 1, 2, 3 initialement
    C._initVars(t,'{centers:cellNIBC}={centers:cellNIBC}*({centers:cellNIBC}<2.5)')
    C._cpVars(t,'centers:cellNIBC',t,'centers:cellN')
    C._cpVars(t,'centers:cellN',tc,'cellN')

    # Transfert du cellNFront
    C._cpVars(t,'centers:cellNFront',tc,'cellNFront')


    # propager cellNVariable='cellNFront'
    Xmpi._setInterpTransfers(t,tc,variables=['cellNFront'], cellNVariable='cellNFront', compact=0)

    if twoFronts:
        C._cpVars(t,'centers:cellNFront_2',tc,'cellNFront_2')
        Xmpi._setInterpTransfers(t,tc,variables=['cellNFront_2'], cellNVariable='cellNFront_2', compact=0)

    ############################################################
    # Specific treatment for front 2
    ############################################################
    if frontType == 2:
        test.printMem(">>> pushBackImageFront2 [start]")

        # bboxDict needed for optimised AddXZones (i.e. "layers" not None)
        # Return a dict with the zones of t as keys and their specific bboxes as key values
        bboxDict = Cmpi.createBboxDict(t)
        tbbc = Cmpi.createBBoxTree(tc)
        interDict = X.getIntersectingDomains(tbbc)
        graph = Cmpi.computeGraph(tbbc, type='bbox', intersectionsDict=interDict, reduction=False)

        # if subr, the tree subregions are kept during the exchange
        # if layers not None, only communicate the desired number of layers
        Cmpi._addLXZones(tc, graph, variables=['cellNIBC','cellNChim','cellNFront'], cartesian=cartesian, interDict=interDict, bboxDict=bboxDict, layers=4, subr=False)
        Cmpi._addLXZones(t, graph, variables=['centers:cellNIBC', 'centers:cellNChim', 'centers:cellNFront'], cartesian=cartesian, interDict=interDict, bboxDict=bboxDict, layers=4, subr=False)

        # Zones of tc are modified after addXZones, new tbbc, interDict and intersectionDict
        tbbcx = G.BB(tc)
        interDict = X.getIntersectingDomains(tbbcx)
        intersectionsDict = X.getIntersectingDomains(tbbcx, method='AABB', taabb=tbbcx)

        # Reconstruction of cellNFront and cellN from cellNIBC (reduce the communications)
        # cellNFront_origin and cellNIBC_origin are initialised to store the Data of cellNFront and cellNIBC before the transfers
        C._initVars(t,'{centers:cellN}={centers:cellNIBC}')
        C._initVars(t,'{centers:cellNFront_origin}={centers:cellNFront}')
        C._initVars(t,'{centers:cellNIBC_origin}={centers:cellNIBC}')
        C._initVars(t,'{centers:cellN_interp}=maximum(0.,{centers:cellNChim})') # Second way of building the cellN field, see above

        C._cpVars(t,'centers:cellNFront',tc,'cellNFront')
        C._cpVars(t,'centers:cellNIBC',tc,'cellNIBC')
        C._cpVars(t,'centers:cellN',tc,'cellN')
        C._cpVars(t,'centers:cellN_interp',tc,'cellN_interp')
        C._cpVars(t,'centers:cellNFront_origin',tc,'cellNFront_origin')
        C._cpVars(t,'centers:cellNIBC_origin',tc,'cellNIBC_origin')

        # Find each zone that require the specific treatment
        C._initVars(t,'{centers:cellNFront2}=1.-({centers:cellNFront}<1.)*(abs({centers:cellNChim})>1.)')
        # i.e. if cellNFront_origin == 2 and cellNFront == 1 ou -3 => cellNFront2 = 1

        # Transfers the information at each grid connection
        for z in Internal.getZones(t):
            cellNFront = Internal.getNodeFromName2(z,'cellNFront2')
            if cellNFront != []:
                cellNFront = cellNFront[1]
                sizeTot = cellNFront.shape[0]*cellNFront.shape[1]*cellNFront.shape[2]
                sizeOne =  int(numpy.sum(cellNFront))
                if sizeOne < sizeTot:
                    X._setHoleInterpolatedPoints(z, depth=1, dir=0, loc='centers',cellNName='cellNFront2',addGC=False)
                    res = X.getInterpolatedPoints(z,loc='centers', cellNName='cellNFront2') # indices,X,Y,Z
                    if res is not None:
                        indicesI = res[0]
                        XI = res[1]; YI = res[2]; ZI = res[3]
                        allInterpFields=[]
                        for zc in Internal.getZones(tc):
                            if zc[0] in intersectionsDict[z[0]]:
                                C._cpVars(zc,'cellN_interp',zc,'cellN')
                                fields = X.transferFields(zc, XI, YI, ZI, hook=None, variables=['cellNFront_origin','cellNIBC_origin'], interpDataType=interpDataType, nature=1)
                                allInterpFields.append(fields)
                        if allInterpFields!=[]:
                            C._filterPartialFields(z, allInterpFields, indicesI, loc='centers', startFrom=0, filterName='donorVol',verbose=False)

        Cmpi._rmXZones(tc)
        Cmpi._rmXZones(t)

        # Update the cellNFront, cellNIBC and cellNIBCDnr fields
        for z in Internal.getZones(t):
            cellNFront = Internal.getNodeFromName2(z,'cellNFront2')
            if cellNFront != []:
                cellNFront = cellNFront[1]
                sizeTot = cellNFront.shape[0]*cellNFront.shape[1]*cellNFront.shape[2]
                sizeOne =  int(numpy.sum(cellNFront))
                if sizeOne < sizeTot:
                    C._initVars(z,'{centers:cellNFront}={centers:cellNFront}*({centers:cellNFront_origin}>0.5)') # Modification du Front uniquement lorsque celui-ci est repousse
                    # i.e. if cellNFront_origin == 0 and cellNFront == 1 => cellNfront = 0

                    C._initVars(z,'{centers:cellNIBC}={centers:cellNIBC}*(1.-({centers:cellNChim}==1.)*({centers:cellNIBC_origin}>1.5)*({centers:cellNIBC_origin}<2.5)) \
                        + 2.*({centers:cellNChim}==1.)*({centers:cellNIBC_origin}>1.5)*({centers:cellNIBC_origin}<2.5)')
                    # i.e. if cellNChim == 1 and cellNIBC_origin == 2 => cellNIBC = 2

                    C._initVars(z,'{centers:cellNIBCDnr}={centers:cellNIBCDnr}*(1.-({centers:cellNChim}==1.)*({centers:cellNIBC_origin}>1.5)*({centers:cellNIBC_origin}<2.5)) \
                        + 2.*({centers:cellNChim}==1.)*({centers:cellNIBC_origin}>1.5)*({centers:cellNIBC_origin}<2.5)')

        C._cpVars(t,'centers:cellNIBC',tc,'cellNIBC')
        C._cpVars(t,'centers:cellNIBC',t,'centers:cellN')
        C._cpVars(t,'centers:cellN',tc,'cellN')

        C._rmVars(t,['centers:cellNFront2'])
        C._rmVars(t,['centers:cellNFront_origin'])
        C._rmVars(t,['centers:cellNIBC_origin'])
        C._rmVars(t,['centers:cellN_interp'])

        # Smooth the front in case of a local refinement - only work in 2D
        if smoothing and dimPb == 2: X_IBM._smoothImageFront(t, tc)

        C._cpVars(t,'centers:cellNFront',tc,'cellNFront')

        Xmpi._setInterpTransfers(t,tc,variables=['cellNFront'], cellNVariable='cellNFront', compact=0)
        test.printMem(">>> pushBackImageFront2 [end]")
    ############################################################

    C._rmVars(t,['centers:cellNFront'])
    if twoFronts:
        C._rmVars(t,['centers:cellNFront_2', 'centers:cellNIBC_2'])

    C._cpVars(t,'centers:TurbulentDistance',tc,'TurbulentDistance')

    print('Minimum distance: %f.'%C.getMinValue(t,'centers:TurbulentDistance'))
    P._computeGrad2(t, 'centers:TurbulentDistance', ghostCells=True, withCellN=False)


    test.printMem(">>> Building IBM front [start]")
    front = X_IBM.getIBMFront(tc, 'cellNFront', dim=dimPb, frontType=frontType)
    front = X_IBM.gatherFront(front)

    if twoFronts:
        front2 = X_IBM.getIBMFront(tc, 'cellNFront_2', dim=dimPb, frontType=frontType)
        front2 = X_IBM.gatherFront(front2)

    if check and rank == 0:
        C.convertPyTree2File(front, 'front.cgns')
        if twoFronts: C.convertPyTree2File(front2, 'front2.cgns')

    zonesRIBC = []
    for zrcv in Internal.getZones(t):
        if C.getMaxValue(zrcv, 'centers:cellNIBC')==2.:
            zrcvname = zrcv[0]; zonesRIBC.append(zrcv)

    nbZonesIBC = len(zonesRIBC)
    if nbZonesIBC == 0:
        res = [{},{},{}]
        if twoFronts: res2 = [{},{},{}]
    else:
        print("isOrtho_Project_First=",isOrtho_Project_First)
        res    = X_IBM.getAllIBMPoints(zonesRIBC, loc='centers',tb=tb, tfront=front, frontType=frontType,
                                       cellNName='cellNIBC', depth=DEPTH, IBCType=IBCType, Reynolds=Reynolds, yplus=yplus, Lref=Lref, isOrthoFirst=isOrtho_Project_First)
        if twoFronts:
            res2 = X_IBM.getAllIBMPoints(zonesRIBC, loc='centers',tb=tb, tfront=front2, frontType=frontType,
                                         cellNName='cellNIBC', depth=DEPTH, IBCType=IBCType, Reynolds=Reynolds, yplus=yplus, Lref=Lref)
        if isWireModel:
            res2 = X_IBM.getAllIBMPoints(zonesRIBC, loc='centers',tb=tb, tfront=front, frontType=frontType,
                                         cellNName='cellNIBC', depth=DEPTH, IBCType=IBCType, Reynolds=Reynolds, yplus=yplus, Lref=Lref,isWireModel=isWireModel,
                                         isOrthoFirst=isOrtho_Project_First)

    # cleaning
    C._rmVars(tc,['cellNChim','cellNIBC','TurbulentDistance','cellNFront'])
    # dans t, il faut cellNChim et cellNIBCDnr pour recalculer le cellN a la fin
    varsRM = ['centers:gradxTurbulentDistance','centers:gradyTurbulentDistance','centers:gradzTurbulentDistance','centers:cellNFront','centers:cellNIBC']
    C._rmVars(t, varsRM)
    front = None
    if twoFronts: front2 = None
    test.printMem(">>> Building IBM front [end]")

    # Interpolation IBC (front, tbbc)

    # graph d'intersection des pts images de ce proc et des zones de tbbc
    zones = Internal.getZones(tbbc)
    allBBs = []
    dictOfCorrectedPtsByIBCType = res[0]
    dictOfWallPtsByIBCType = res[1]
    dictOfInterpPtsByIBCType = res[2]
    interDictIBM={}

    if twoFronts or isWireModel:
        dictOfCorrectedPtsByIBCType2 = res2[0]
        dictOfWallPtsByIBCType2 = res2[1]
        dictOfInterpPtsByIBCType2 = res2[2]
        interDictIBM2={}
    else:
        dictOfCorrectedPtsByIBCType2={}
        dictOfWallPtsByIBCType2={}
        dictOfInterpPtsByIBCType2={}
        interDictIBM2={}

    if dictOfCorrectedPtsByIBCType!={}:
        for ibcTypeL in dictOfCorrectedPtsByIBCType:
            allCorrectedPts = dictOfCorrectedPtsByIBCType[ibcTypeL]
            allWallPts = dictOfWallPtsByIBCType[ibcTypeL]
            allInterpPts = dictOfInterpPtsByIBCType[ibcTypeL]
            for nozr in range(nbZonesIBC):
                if allCorrectedPts[nozr] != []:
                    zrname = zonesRIBC[nozr][0]
                    interpPtsBB = Generator.BB(allInterpPts[nozr])
                    for z in zones:
                        bba = C.getFields('GridCoordinates', z, api=1)[0]
                        if Generator.bboxIntersection(interpPtsBB,bba,isBB=True):
                            zname = z[0]
                            popp = Cmpi.getProc(z)
                            Distributed.updateGraph__(graph, popp, rank, zname)
                            if zrname not in interDictIBM: interDictIBM[zrname]=[zname]
                            else:
                                if zname not in interDictIBM[zrname]: interDictIBM[zrname].append(zname)
        if twoFronts or isWireModel:
            for ibcTypeL in dictOfCorrectedPtsByIBCType2:
                allCorrectedPts2 = dictOfCorrectedPtsByIBCType2[ibcTypeL]
                allWallPts2 = dictOfWallPtsByIBCType2[ibcTypeL]
                allInterpPts2 = dictOfInterpPtsByIBCType2[ibcTypeL]
                for nozr in range(nbZonesIBC):
                    if allCorrectedPts2[nozr] != []:
                        zrname = zonesRIBC[nozr][0]
                        interpPtsBB2 = Generator.BB(allInterpPts2[nozr])
                        for z in zones:
                            bba = C.getFields('GridCoordinates', z, api=1)[0]
                            if Generator.bboxIntersection(interpPtsBB2,bba,isBB=True):
                                zname = z[0]
                                popp = Cmpi.getProc(z)
                                Distributed.updateGraph__(graph, popp, rank, zname)
                                if zrname not in interDictIBM2: interDictIBM2[zrname]=[zname]
                                else:
                                    if zname not in interDictIBM2[zrname]: interDictIBM2[zrname].append(zname)
    else: graph={}
    del tbbc
    allGraph = Cmpi.KCOMM.allgather(graph)
    #if rank == 0: print allGraph

    graph = {}
    for i in allGraph:
        for k in i:
            if not k in graph: graph[k] = {}
            for j in i[k]:
                if not j in graph[k]: graph[k][j] = []
                graph[k][j] += i[k][j]
                graph[k][j] = list(set(graph[k][j])) # pas utile?

    test.printMem(">>> Interpolating IBM [start]")
    # keyword subr=False to avoid memory overflow
    Cmpi._addXZones(tc, graph, variables=['cellN'], cartesian=cartesian, subr=False)
    test.printMem(">>> Interpolating IBM [after addXZones]")

    ReferenceState = Internal.getNodeFromType2(t, 'ReferenceState_t')
    nbZonesIBC = len(zonesRIBC)

    for i in range(Cmpi.size): datas[i] = [] # force

    if dictOfCorrectedPtsByIBCType!={}:
        for ibcTypeL in dictOfCorrectedPtsByIBCType:
            allCorrectedPts = dictOfCorrectedPtsByIBCType[ibcTypeL]
            allWallPts = dictOfWallPtsByIBCType[ibcTypeL]
            allInterpPts = dictOfInterpPtsByIBCType[ibcTypeL]
            for nozr in range(nbZonesIBC):
                if allCorrectedPts[nozr] != []:
                    zrcv = zonesRIBC[nozr]
                    zrname = zrcv[0]
                    dnrZones = []
                    for zdname in interDictIBM[zrname]:
                        zd = Internal.getNodeFromName2(tc, zdname)
                        #if zd is not None: dnrZones.append(zd)
                        if zd is None: print('!!!Zone None', zrname, zdname)
                        else: dnrZones.append(zd)
                    XOD._setIBCDataForZone__(zrcv, dnrZones, allCorrectedPts[nozr], allWallPts[nozr], allInterpPts[nozr],
                                             nature=1, penalty=1, loc='centers', storage='inverse', dim=dimPb,
                                             interpDataType=interpDataType, ReferenceState=ReferenceState, bcType=ibcTypeL)

                    nozr += 1
                    for zd in dnrZones:
                        zdname = zd[0]
                        destProc = procDict[zdname]

                        #allIDs = Internal.getNodesFromName(zd, 'IBCD*')
                        #IDs = []
                        #for zsr in allIDs:
                        #    if Internal.getValue(zsr)==zrname: IDs.append(zsr)

                        IDs = []
                        for i in zd[2]:
                            if i[0][0:4] == 'IBCD':
                                if Internal.getValue(i)==zrname: IDs.append(i)

                        if IDs != []:
                            if destProc == rank:
                                zD = Internal.getNodeFromName2(tc,zdname)
                                zD[2] += IDs
                            else:
                                if destProc not in datas: datas[destProc]=[[zdname,IDs]]
                                else: datas[destProc].append([zdname,IDs])
                        else:
                            if destProc not in datas: datas[destProc] = []

    if dictOfCorrectedPtsByIBCType2!={}:
        for ibcTypeL in dictOfCorrectedPtsByIBCType2:
            allCorrectedPts2 = dictOfCorrectedPtsByIBCType2[ibcTypeL]
            allWallPts2 = dictOfWallPtsByIBCType2[ibcTypeL]
            allInterpPts2 = dictOfInterpPtsByIBCType2[ibcTypeL]
            for nozr in range(nbZonesIBC):
                if allCorrectedPts2[nozr] != []:
                    zrcv = zonesRIBC[nozr]
                    zrname = zrcv[0]
                    dnrZones = []
                    for zdname in interDictIBM2[zrname]:
                        zd = Internal.getNodeFromName2(tc, zdname)
                        #if zd is not None: dnrZones.append(zd)
                        if zd is None: print('!!!Zone None', zrname, zdname)
                        else: dnrZones.append(zd)
                    XOD._setIBCDataForZone2__(zrcv, dnrZones, allCorrectedPts2[nozr], allWallPts2[nozr], None, allInterpPts2[nozr],
                                              nature=1, penalty=1, loc='centers', storage='inverse', dim=dimPb,
                                              interpDataType=interpDataType, ReferenceState=ReferenceState, bcType=ibcTypeL)

                    nozr += 1
                    for zd in dnrZones:
                        zdname = zd[0]
                        destProc = procDict[zdname]

                        IDs = []
                        for i in zd[2]:
                            if i[0][0:6] == '2_IBCD':
                                if Internal.getValue(i)==zrname: IDs.append(i)

                        if IDs != []:
                            if destProc == rank:
                                zD = Internal.getNodeFromName2(tc,zdname)
                                zD[2] += IDs
                            else:
                                if destProc not in datas: datas[destProc]=[[zdname,IDs]]
                                else: datas[destProc].append([zdname,IDs])
                        else:
                            if destProc not in datas: datas[destProc] = []

    test.printMem(">>> Interpolating IBM [end]")
    Cmpi._rmXZones(tc)
    dictOfCorrectedPtsByIBCType = None
    dictOfWallPtsByIBCType = None
    dictOfInterpPtsByIBCType = None
    interDictIBM = None
    if twoFronts or isWireModel:
        dictOfCorrectedPtsByIBCType2 = None
        dictOfWallPtsByIBCType2 = None
        dictOfInterpPtsByIBCType2 = None
        interDictIBM2 = None
    test.printMem(">>> Interpolating IBM [after rm XZones]")

    Internal._rmNodesByName(tc, Internal.__FlowSolutionNodes__)
    #Internal._rmNodesByName(tc, Internal.__GridCoordinates__)
    destDatas = Cmpi.sendRecv(datas, graph)
    for i in destDatas:
        for n in destDatas[i]:
            zname = n[0]
            IBCDs = n[1]
            if IBCDs != []:
                zD = Internal.getNodeFromName2(tc, zname)
                zD[2] += IBCDs

    datas = {}; graph = {}
    C._initVars(t,'{centers:cellN}=minimum({centers:cellNChim}*{centers:cellNIBCDnr},2.)')
    varsRM = ['centers:cellNChim','centers:cellNIBCDnr']
    if model == 'Euler': varsRM += ['centers:TurbulentDistance']
    C._rmVars(t, varsRM)

    #-----------------------------------------
    # Computes distance field for Musker only
    #-----------------------------------------
    # zones = Internal.getZones(t)
    # npts = 0
    # for z in zones:
    #     dims = Internal.getZoneDim(z)
    #     npts += dims[1]*dims[2]*dims[3]
    # Cmpi.barrier()
    # print('proc {} has {} blocks and {} Millions points'.format(rank, len(zones), npts/1.e6))

    if model != 'Euler' and recomputeDist and (extrusion!='cyl' and extrusion !='cart'):
        ibctypes = set()
        for node in Internal.getNodesFromName(tb,'ibctype'):
            ibctypes.add(Internal.getValue(node))
        ibctypes = list(ibctypes)
        if 'outpress' in ibctypes or 'inj' in ibctypes or 'slip' in ibctypes:
            test.printMem(">>> wall distance for viscous wall only [start]")
            for z in Internal.getZones(tb):
                ibc = Internal.getNodeFromName(z,'ibctype')
                if Internal.getValue(ibc)=='outpress' or Internal.getValue(ibc)=='inj' or Internal.getValue(ibc)=='slip':
                    Internal._rmNode(tb,z)
            if dimPb == 2:
                z0 = Internal.getZones(t)
                bb = G.bbox(z0); dz = bb[5]-bb[2]
                tb2 = C.initVars(tb, 'CoordinateZ', dz*0.5)
                DTW._distance2Walls(t,tb2,type='ortho', signed=0, dim=dimPb, loc='centers')
            else:
                DTW._distance2Walls(t,tb,type='ortho', signed=0, dim=dimPb, loc='centers')
            test.printMem(">>> wall distance for viscous wall only [end]")

            if dimPb == 2 and cleanCellN == False: C._initVars(t, '{centers:TurbulentDistanceWallBC}={centers:TurbulentDistance}')
    else:
        for z in Internal.getZones(t):
            dist = Internal.getNodeFromName2(z,'TurbulentDistanceWallBC')
            if dist is not None:  C._initVars(t, '{centers:TurbulentDistance}={centers:TurbulentDistanceWallBC}')

    # Sauvegarde des infos IBM
    if check:
        test.printMem(">>> Saving IBM infos [start]")
        tibm = X_IBM.extractIBMInfo(tc)

        # Avoid that two procs write the same information
        for z in Internal.getZones(tibm):
            if int(z[0][-1]) != rank:
                # Internal._rmNodesByName(tibm, z[0])
                z[0] = z[0]+"%{}".format(rank)

        Cmpi.convertPyTree2File(tibm, 'IBMInfo.cgns')


        if twoFronts or isWireModel:
            tibm2 = X_IBM.extractIBMInfo2(tc)

            # Avoid that two procs write the same information
            for z in Internal.getZones(tibm2):
                if int(z[0][-1]) != rank:
                    # Internal._rmNodesByName(tibm, z[0])
                    z[0] = z[0]+"%{}".format(rank)

            Cmpi.convertPyTree2File(tibm2, 'IBMInfo2.cgns')

        test.printMem(">>> Saving IBM infos [end]")
        del tibm
        if twoFronts or isWireModel: del tibm2

    # distribution par defaut (sur NP)
    tbbc = Cmpi.createBBoxTree(tc)

    # Perform the final distribution
    if distrib:
        if NP == 0: NP = Cmpi.size
        stats = D2._distribute(tbbc, NP, algorithm='graph', useCom='ID')
        D2._copyDistribution(tc, tbbc)
        D2._copyDistribution(t, tbbc)

    del tbbc

    if extrusion == 'cyl':
        T._cyl2Cart(t, (0,0,0),(1,0,0))
        T._cyl2Cart(tc,(0,0,0),(1,0,0))
        # modif info maillage des zonesubregion_t
        for z in Internal.getZones(tc):
            for zsr in Internal.getNodesFromType(z, "ZoneSubRegion_t"):
                zsrname = Internal.getName(zsr)
                zsrname = zsrname.split('_')
                if zsrname[0]=='IBCD':
                    for var in ['C','W','I']:
                        r     = Internal.getNodeFromName(zsr,'CoordinateY_P'+var)[1]
                        theta = Internal.getNodeFromName(zsr,'CoordinateZ_P'+var)[1]
                        #print('cyl2cart', var, numpy.size(r), zsr[0], z[0], toto,tutu )
                        for l in range(numpy.size(r)):
                            yy  = r[l]*numpy.cos( theta[l] )
                            zz  = r[l]*numpy.sin( theta[l] )
                            r[l]= yy; theta[l] = zz
    if redistribute:
        import Distributor2.Mpi as D2mpi
        tcs    = Cmpi.allgatherTree(Cmpi.convert2SkeletonTree(tc))
        stats  = D2._distribute(tcs, NP, algorithm='graph')
        if rank == 0: checkNcellsNptsPerProc(tcs, Cmpi.size, isAtCenter=True)
        D2._copyDistribution(tc, tcs)
        D2._copyDistribution(t , tcs)
        D2mpi._redispatch(tc)
        D2mpi._redispatch(t)

    #-----------------------------------------
    # Computes distance field for Musker only
    #-----------------------------------------
    # zones = Internal.getZones(t)
    # npts = 0
    # for z in zones:
    #     dims = Internal.getZoneDim(z)
    #     npts += dims[1]*dims[2]*dims[3]
    # Cmpi.barrier()
    # print('proc {} has {} blocks and {} Millions points'.format(rank, len(zones), npts/1.e6))

    if model == 'NSTurbulent':
        test.printMem(">>> wall distance for viscous wall only - RANS [start]")

        ibctypes = set()
        for node in Internal.getNodesFromName(tb,'ibctype'):
            ibctypes.add(Internal.getValue(node))
        ibctypes = list(ibctypes)
        if 'outpress' in ibctypes or 'inj' in ibctypes or 'slip' in ibctypes:

            for z in Internal.getZones(tb):
                ibc = Internal.getNodeFromName(z,'ibctype')
                if Internal.getValue(ibc)=='outpress' or Internal.getValue(ibc)=='inj' or Internal.getValue(ibc)=='slip':
                    Internal._rmNode(tb,z)

        if dimPb == 2:
            DTW._distance2Walls(t,tbsave,type='ortho', signed=0, dim=dimPb, loc='centers')
        else:
            DTW._distance2Walls(t,tb,type='ortho', signed=0, dim=dimPb, loc='centers')
        C._initVars(t, '{centers:TurbulentDistance}={centers:TurbulentDistance}*({centers:cellN}>0.)+(-1.)*{centers:TurbulentDistance}*({centers:cellN}<1.)')

        test.printMem(">>> wall distance for viscous wall only - RANS [end]")
    # Save tc
    if twoFronts or isWireModel:
        tc2 = Internal.copyTree(tc)
        tc2 = Internal.rmNodesByName(tc2, 'IBCD*')
        tc  = Internal.rmNodesByName(tc, '2_IBCD*')

        if isWireModel:tc2 = Internal.rmNodesByName(tc2, 'ID*')

    if RENAMEIBCNODES:
        for zc in Internal.getZones(tc):
            for ibcd in Internal.getNodesFromName1(zc,'IBCD_*'):
                proposedName = Internal.getName(ibcd)[0:6]+'_X%d'%(rank)
                ibcd[0]=getIBCDName(proposedName)

        if twoFronts or isWireModel:
            for zc in Internal.getZones(tc2):
                for ibcd in Internal.getNodesFromName1(zc,'2_IBCD_*'):
                    proposedName = Internal.getName(ibcd)[0:8]+'_X%d'%(rank)
                    ibcd[0]=getIBCDName(proposedName)

    if isinstance(tc_out, str):
        if isWireModel:
            tc2  = transformTc2(tc2)
            NewIBCD=140
            _change_name_IBCD(tc,NewIBCD)
            NewIBCD=141
            _change_name_IBCD(tc2,NewIBCD)
            tc=Internal.merge([tc,tc2])
            del tc2

        tcp = Compressor.compressCartesian(tc)
        Cmpi.convertPyTree2File(tcp, tc_out, ignoreProcNodes=True)

        if twoFronts:
            tc2  = transformTc2(tc2)
            tcp2 = Compressor.compressCartesian(tc2)
            Cmpi.convertPyTree2File(tcp2, 'tc2.cgns', ignoreProcNodes=True)
            del tc2


    # Initialisation
    if tinit is None: I._initConst(t, loc='centers')
    else:
        t = Pmpi.extractMesh(tinit, t, mode='accurate')
    if model != "Euler": C._initVars(t, 'centers:ViscosityEddy', 0.)

    if isWireModel:
        vars_wm = ['Density','VelocityX','VelocityY','VelocityZ','Temperature']
        if model == 'NSTurbulent':vars_wm.append('TurbulentSANuTilde')
        vars_wm.append('DistWall2IP')
        for z in Internal.getZones(t):
            #if C.getMaxValue(C.rmGhostCells(z, z, 2, adaptBCs=1), 'centers:cellN')>1:
            for v_local in vars_wm:
                C._initVars(z,'{centers:'+v_local+'_WM}=0.')

    # Init with BBox
    if initWithBBox>0.:
        print('initialisation par bounding box')
        bodybb = C.newPyTree(['Base'])
        for base in Internal.getBases(tb):
            bbox = G.bbox(base)
            bodybbz = D.box(tuple(bbox[:3]),tuple(bbox[3:]), N=2, ntype='STRUCT')
            Internal._append(bodybb,bodybbz,'Base')
        T._scale(bodybb, factor=(initWithBBox,initWithBBox,initWithBBox))
        tbb = G.BB(t)
        interDict = X.getIntersectingDomains(tbb,bodybb,taabb=tbb,taabb2=bodybb)
        for zone in Internal.getZones(t):
            zname = Internal.getName(zone)
            if interDict[zname] != []:
                C._initVars(zone, 'centers:MomentumX', 0.)
                C._initVars(zone, 'centers:MomentumY', 0.)
                C._initVars(zone, 'centers:MomentumZ', 0.)

    if not redistribute:
        tc_tmp = Cmpi.allgatherTree(Cmpi.convert2SkeletonTree(tc))
        if rank == 0: checkNcellsNptsPerProc(tc_tmp, Cmpi.size, isAtCenter=True)
        del tc_tmp

    # Save t
    if isinstance(t_out, str):
        tp = Compressor.compressCartesian(t)
        Cmpi.convertPyTree2File(tp, t_out, ignoreProcNodes=True)

    #Nettoyage arbre
    if extrusion is not None:
        vars = ['centers:TurbulentDistanceAllBC','centers:TurbulentDistanceWallBC', 'centers:cellNIBC_hole']
        C._rmVars(t, vars)


    if Cmpi.size > 1: Cmpi.barrier()
    return t, tc

#====================================================================================
# Redistrib on NP processors
#====================================================================================
def _distribute(t_in, tc_in, NP, algorithm='graph', tc2_in=None):
    if isinstance(tc_in, str):
        tcs = Cmpi.convertFile2SkeletonTree(tc_in, maxDepth=3)
    else: tcs = tc_in
    stats = D2._distribute(tcs, NP, algorithm=algorithm, useCom='ID')
    if isinstance(tc_in, str):
        paths = []; ns = []
        bases = Internal.getBases(tcs)
        for b in bases:
            zones = Internal.getZones(b)
            for z in zones:
                nodes = Internal.getNodesFromName2(z, 'proc')
                for n in nodes:
                    p = 'CGNSTree/%s/%s/.Solver#Param/proc'%(b[0],z[0])
                    paths.append(p); ns.append(n)
        Filter.writeNodesFromPaths(tc_in, paths, ns, maxDepth=0, mode=1)

    if isinstance(t_in, str):
        ts = Cmpi.convertFile2SkeletonTree(t_in, maxDepth=3)
    else: ts = t_in
    D2._copyDistribution(ts, tcs)

    if isinstance(t_in, str):
        paths = []; ns = []
        bases = Internal.getBases(ts)
        for b in bases:
            zones = Internal.getZones(b)
            for z in zones:
                nodes = Internal.getNodesFromName2(z, 'proc')
                for n in nodes:
                    p = 'CGNSTree/%s/%s/.Solver#Param/proc'%(b[0],z[0])
                    paths.append(p); ns.append(n)
        Filter.writeNodesFromPaths(t_in, paths, ns, maxDepth=0, mode=1)

    if tc2_in is not None:
        if isinstance(tc2_in, str):
            tc2s = Cmpi.convertFile2SkeletonTree(tc2_in, maxDepth=3)
        else: tc2s = tc2_in
        D2._copyDistribution(tc2s, tcs)

        if isinstance(tc2_in, str):
            paths = []; ns = []
            bases = Internal.getBases(tc2s)
            for b in bases:
                zones = Internal.getZones(b)
                for z in zones:
                    nodes = Internal.getNodesFromName2(z, 'proc')
                    for n in nodes:
                        p = 'CGNSTree/%s/%s/.Solver#Param/proc'%(b[0],z[0])
                        paths.append(p); ns.append(n)
            Filter.writeNodesFromPaths(tc2_in, paths, ns, maxDepth=0, mode=1)

    checkNcellsNptsPerProc(ts,NP)
    return None


def checkNcellsNptsPerProc(ts,NP,isAtCenter=False):
    # Affichage du nombre de points par proc - equilibrage ou pas
    NptsTot = 0
    NcellsTot = 0
    ncellslocal=[]
    nptslocal  =[]
    for i in range(NP):
        NPTS = 0
        NCELLS = 0
        for z in Internal.getZones(ts):
            if Cmpi.getProc(z) == i:
                NPTS += C.getNPts(z)
                if not isAtCenter:
                    NCELLS += C.getNCells(z)
                else:
                    NCELLS=NPTS
        ncellslocal.append(NCELLS)
        nptslocal.append(NPTS)

        NptsTot   += NPTS
        NcellsTot += NCELLS
        if isAtCenter:
            print('Rank {} has {} cells'.format(i,NCELLS))
        else:
            print('Rank {} has {} points & {} cells'.format(i,NPTS,NCELLS))
    print('All points: {} million points & {} million cells'.format(NptsTot/1.e6,NcellsTot/1.e6))

    for i in range(NP):
        ncellslocal[i] = ncellslocal[i]/NcellsTot*100
        print('Rank {} :: {} % of cells'.format(i,ncellslocal[i]))
    print('Range of % of cells: {} - {}'.format(min(ncellslocal),max(ncellslocal)))

    return None


class IBM(Common):
    """Preparation et calculs IBM avec le module FastS."""
    def __init__(self, format=None, numb=None, numz=None):
        Common.__init__(self, format, numb, numz)
        self.__version__ = "0.0"
        self.authors = ["ash@onera.fr"]
        self.cartesian = True

    # Prepare
    def prepare(self, t_case, t_out, tc_out, snears=0.01, dfar=10., dfarList=[],
                tbox=None, snearsf=None, yplus=100.,
                vmin=21, check=False, frontType=1, NP=None, expand=3, tinit=None,
                initWithBBox=-1., wallAdapt=None, dfarDir=0, redistribute=False):
        if NP is None: NP = Cmpi.size
        if NP == 0: print('Preparing for a sequential computation.')
        else: print('Preparing for an IBM computation on %d processors.'%NP)
        ret = prepare(t_case, t_out, tc_out, snears=snears, dfar=dfar, dfarList=dfarList,
                      tbox=tbox, snearsf=snearsf, yplus=yplus,
                      vmin=vmin, check=check, NP=NP, format=self.data['format'],
                      frontType=frontType, expand=expand, tinit=tinit, dfarDir=dfarDir,
                      redistribute=redistribute)
        return ret

    # post-processing: extrait la solution aux noeuds + le champs sur les surfaces
    def post(self, t_case, t_in, tc_in, t_out, wall_out):
        return post(t_case, t_in, tc_in, t_out, wall_out)

    # post-processing: extrait les efforts sur les surfaces
    def loads(self, t_case, tc_in=None, wall_out=None, alpha=0., beta=0., Sref=None, famZones=[]):
        return loads(t_case, tc_in=tc_in, wall_out=wall_out, alpha=alpha, beta=beta, Sref=Sref, famZones=famZones)



## IMPORTANT NOTE !!
## FUNCTIONS MIGRATED TO $CASSIOPEE/Cassiopee/Geom/Geom/IBM.py
## The functions below will become deprecated after Jan. 1 2023
#====================================================================================
def setSnear(t, value):
    tp=D_IBM.setSnear(t, value)
    return tp

def _setSnear(t, value):
    D_IBM._setSnear(t, value)
    return None

def setDfar(t, value):
    tp=D_IBM.setDfar(t, value)
    return tp

def _setDfar(t, value):
    D_IBM._setDfar(t, value)
    return None

def snearFactor(t, factor):
    tp=D_IBM.snearFactor(t, factor)
    return tp

def _snearFactor(t, factor):
    D_IBM._snearFactor(t, factor)
    return None

def setFluidInside(t):
    return D_IBM.setFluidInside(t)

def _setFluidInside(t):
    return D_IBM._setFluidInside(t)

def setIBCType(t, value):
    tp=D_IBM.setIBCType(t, value)
    return tp

def _setIBCType(t, value):
    D_IBM._setIBCType(t, value)
    return None

def changeBCType(tc, oldBCType, newBCType):
    tc=D_IBM.changeIBCType(tc, oldBCType, newBCType)
    return tc

def initOutflow(tc, familyNameOutflow, P_tot):
    tc=D_IBM.initOutflow(tc, familyNameOutflow, P_tot)
    return tc

def initInj(tc, familyNameInj, P_tot, H_tot, injDir=[1.,0.,0.]):
    tc=D_IBM.initInj(tc, familyNameInj, P_tot, H_tot, injDir)
    return tc

def _initOutflow(tc, familyNameOutflow, P_tot):
    return D_IBM._initOutflow(tc, familyNameOutflow, P_tot)

def _initInj(tc, familyNameInj, P_tot, H_tot, injDir=[1.,0.,0.]):
    return D_IBM._initInj(tc, familyNameInj, P_tot, H_tot, injDir)

def transformTc2(tc2):
    tc2=D_IBM.transformTc2(tc2)
    return tc2


#=============================================================================
# Post - General
# IN: t_case: geometry file name or tree
# IN: t_in: result file name or tree
# IN: tc_in: connectivity file name or tree
# OUT: t_out ou None: output file name - values at nodes
# OUT: wall_out ou None: output file name - wall values
#==============================================================================
def post(t_case, t_in, tc_in, t_out, wall_out):
    if isinstance(t_in, str): t = C.convertFile2PyTree(t_in)
    else: t = t_in
    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = t_case

    #=============================
    # Deleting unnecessary fields
    #=============================
    vars =['centers:TurbulentDistance',
           'centers:Density_M1'       , 'centers:Temperature_M1',
           'centers:VelocityX_M1'     , 'centers:VelocityY_M1'  , 'centers:VelocityZ_M1',
           'centers:Density_P1'       , 'centers:Temperature_P1',
           'centers:VelocityX_P1'     , 'centers:VelocityY_P1'  , 'centers:VelocityZ_P1']
    C._rmVars(t, vars)

    #=============================
    # Connectivity tree
    #=============================
    if isinstance(tc_in, str): tc = C.convertFile2PyTree(tc_in)
    else: tc = tc_in
    Internal._rmNodesByName(tc, 'GridCoordinates')

    #==========================================================
    # Compute Cp, Cf, ... on the geometry surface (interpolation)
    #==========================================================
    tb = C.convertArray2Tetra(tb)

    #--------------------------------------------------------
    # Get Reference State and model from body pyTree
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('GoverningEquations is missing in input cgns.')
    model = Internal.getValue(model)

    if model == 'Euler': bcType = 0 # slip
    elif model =='NSLaminar': bcType = 1 # noslip
    else: bcType = 3 # Musker

    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)

    [RoInf, RouInf, RovInf, RowInf, RoeInf, PInf, TInf, cvInf, MInf,
     ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,Mus, Cs, Ts, Pr] = C.getState(tb)

    varType = 2 # IBM updated variables (rho,u,t)
    varsIBC = ['Density', 'VelocityX', 'VelocityY', 'VelocityZ', 'Temperature']
    vars    = ['Density', 'VelocityX', 'VelocityY', 'VelocityZ', 'Temperature']


    if model != 'Euler':
        vars += ['ViscosityEddy']
        if model == 'NSTurbulent':
            vars += ['TurbulentSANuTilde']
            varsIBC += ['TurbulentSANuTilde']
            varType = 21

    for z in Internal.getNodesFromType2(t, "Zone_t"):
        zc = Internal.getNodeFromName(tc, z[0])
        for v in varsIBC: C._cpVars(z, 'centers:'+v, zc, v)

    X._setInterpTransfers(t, tc, variables=vars,
                          variablesIBC=varsIBC, bcType=bcType,
                          varType=varType, storage=1,
                          Gamma=Gamma, Cv=cvInf, MuS=Mus,
                          Cs=Cs, Ts=Ts)

    zw = P_IBM.extractIBMWallFields(tc, tb=tb)
    RoUInf2I = 1./(RouInf*RouInf+RovInf*RovInf+RowInf*RowInf)
    C._initVars(zw,'{Cp}=2*%f*({Pressure}-%f)*%f'%(RoInf,PInf,RoUInf2I))
    if model != 'Euler':
        C._initVars(zw,'{Cf}=2*%f*{Density}*{utau}**2*%f'%(RoInf,RoUInf2I))

    Internal._rmNodesByName(zw, '.Solver#Param')
    Internal._rmNodesByName(zw, '.Solver#ownData')

    if isinstance(wall_out, str): C.convertPyTree2File(zw, wall_out)

    #================================
    # For 2D, extract a single k plane
    #================================
    if dimPb == 2:
        t = T.subzone(t, (1,1,1), (-1,-1,1))
        C._initVars(t, 'CoordinateZ', 0.) # forced

    #=================================
    # Calc. mu_t/mu in the flow field
    #=================================
    if model != 'Euler':
        betas = Mus*(Ts+Cs)/(Ts**(3./2.))
        C._initVars(t, '{centers:ViscosityMolecular} = %20.16g*sqrt({centers:Temperature})/(1.+%20.16g/{centers:Temperature})'%(betas,Cs))
        C._initVars(t, '{centers:mutsmu}=({centers:ViscosityEddy})/({centers:ViscosityMolecular})-1.')

    #======================================
    # Output of flow solution at cell nodes
    #======================================
    vars = ['centers:Density','centers:VelocityX', 'centers:VelocityY', 'centers:VelocityZ', 'centers:Temperature','centers:ViscosityEddy',
            'centers:TurbulentSANuTilde','centers:ViscosityMolecular', 'centers:mutsmu', 'centers:cellN']
    #vars = ['centers:Density','centers:VelocityX', 'centers:Temperature','centers:ViscosityEddy',
    #        'centers:TurbulentSANuTilde','centers:ViscosityMolecular', 'centers:mutsmu', 'centers:cellN']
    t = C.center2Node(t, vars)
    Internal._rmNodesByName(t, 'FlowSolution#Centers')
    if isinstance(t_out, str): C.convertPyTree2File(t, t_out)

    return t, zw


#=============================================================================
# Post efforts
# IN: t_case: geometry tree
# IN: tc_in: connectivity tree
# IN: tc2_in: second connectivity tree (when using 2 image points)
# OUT: wall_out or None: file for the output of the forces on the wall at the centers
# IN: alpha: angle for the computation of the forces
# IN: beta: angle for the computation of the forces
# IN: gradP: calculate the pressure gradient
# IN: order: order of the extrapolation of pressure
# IN: Sref: reference area
# IN : famZones : list of family names of surface zones on which the solution is projected
# NOTE: if tc_in = None, t_case is the geometry tree with the projected solution
#==============================================================================
def loads(t_case, tc_in=None, tc2_in=None, wall_out=None, alpha=0., beta=0., gradP=False, order=1, Sref=None, famZones=[]):
    """Computes the viscous and pressure forces on the IB. If tc_in=None, t_case must also contain the projection of the flow field solution onto the IB.
    Usage: loads(t_case, tc_in, tc2_in, wall_out, alpha, beta, gradP, order, Sref, famZones)"""
    if tc_in is not None:
        if isinstance(tc_in, str):
            tc = C.convertFile2PyTree(tc_in)
        else: tc = tc_in
    else: tc = None

    if tc2_in is not None:
        if isinstance(tc2_in, str):
            tc2 = C.convertFile2PyTree(tc2_in)
        else: tc2 = tc2_in
    else: tc2 = None

    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = t_case

    if Sref is None:
        C._initVars(tb, '__ONE__',1.)
        Sref = P.integ(tb, '__ONE__')[0]; print(Sref)
        C._rmVars(tb, ['__ONE__', 'centers:vol'])

    #====================================
    # Wall pressure correction
    #====================================
    if gradP:
        # add gradP fields in tc if necessary
        if tc is not None:
            for z in Internal.getZones(tc):
                P_IBM._addGradxiP__(z)

            if tc2 is None:
                if order < 2:
                    tc = P_IBM.extractPressureHO(tc)
                else:
                    tc = P_IBM.extractPressureHO2(tc)

        # add gradP fields in tc2 if necessary
        if tc2 is not None:

            for z in Internal.getZones(tc2):
                P_IBM._addGradxiP__(z)

            if order < 2:
                tc2 = P_IBM.extractPressureHO(tc2)
            else:
                tc2 = P_IBM.extractPressureHO2(tc2)


    #====================================
    # Extraction des grandeurs a la paroi
    #====================================
    if tc is None:
        zw = Internal.getZones(tb)
        zw = T.join(zw)
    else:
        zw = P_IBM.extractIBMWallFields(tc, tb=tb, famZones=famZones)

    #====================================
    # Extract pressure info from tc2 to tc
    #====================================
    if tc2 is not None:
        zw2 = P_IBM.extractIBMWallFields(tc2, tb=tb, famZones=famZones, front=1)
        zones_zw  = []
        zones_zw2 = []
        for zone in Internal.getZones(zw): zones_zw.append(zone[0])
        for zone in Internal.getZones(zw2): zones_zw2.append(zone[0])
        nbZones = len(zones_zw)

        for i in range(nbZones): # for multi corps
            szw  = Internal.getNodeFromName(zw, zones_zw[i])
            szw2 = Internal.getNodeFromName(zw2, zones_zw2[i])

            Internal.getNodeFromName(szw, 'Pressure')[1] = Internal.getNodeFromName(szw2, 'Pressure')[1]
            Internal.getNodeFromName(szw, 'Density')[1]  = Internal.getNodeFromName(szw2, 'Density')[1]

            Internal.getNodeFromName(szw, 'gradxPressure')[1] = Internal.getNodeFromName(szw2, 'gradxPressure')[1]
            Internal.getNodeFromName(szw, 'gradyPressure')[1] = Internal.getNodeFromName(szw2, 'gradyPressure')[1]
            Internal.getNodeFromName(szw, 'gradzPressure')[1] = Internal.getNodeFromName(szw2, 'gradzPressure')[1]

    dimPb = Internal.getValue(Internal.getNodeFromName(tb, 'EquationDimension'))

    if dimPb == 2: T._addkplane(zw)

    zw = C.convertArray2Tetra(zw)
    zw = T.reorderAll(zw, 1)

    ts = C.newPyTree(['SKIN']);
    if famZones:ts[2][1][2]=zw
    else:ts[2][1][2]=zw[2][1][2]
    #==============================
    # Reference state
    #==============================
    RefState = Internal.getNodeFromType(tb,'ReferenceState_t')
    ts[2][1][2].append(RefState)

    ts=C.node2Center(ts,'FlowSolution')
    C._rmVars(ts, 'FlowSolution')

    P_IBM._loads0(ts, Sref=Sref, Pref=None, Qref=None, alpha=alpha, beta=beta, dimPb=dimPb, verbose=True)

    if dimPb == 2: # reextrait en 2D
        ts = P.isoSurfMC(ts, "CoordinateZ", 0.)
        nodes = Internal.getNodesFromName(ts, 'CoordinateX')
        xmin = numpy.min(nodes[0][1])
        xmax = numpy.max(nodes[0][1])
        dxi = 1./(xmax-xmin)
        C._initVars(ts, 'xsc=({CoordinateX}-%g)*%g'%(xmin, dxi))

    if isinstance(wall_out, str): C.convertPyTree2File(ts, wall_out)
    return ts

#====================================================================================

## IMPORTANT NOTE !!
## FUNCTIONS MIGRATED TO $CASSIOPEE/Cassiopee/Post/Post/IBM.py
## The functions below will become deprecated after Jan. 1 2023
#====================================================================================
def extractIBMInfo(tc_in, t_out='IBMInfo.cgns'):
    tibm=P_IBM.extractIBMInfo(tc_in, t_out=t_out)
    return tibm

def _loads0(ts, Sref=None, Pref=None, Qref=None, alpha=0., beta=0., dimPb=3, verbose=False):
    return P_IBM.loads0(ts, Sref=Sref, Pref=Pref, Qref=Qref, alpha=alpha, beta=beta, dimPb=dimPb, verbose=verbose)

def loads0(ts, Sref=None, alpha=0., beta=0., dimPb=3, verbose=False):
    return P_IBM.loads0(ts, Sref=Sref, alpha=alpha, beta=beta, dimPb=dimPb, verbose=verbose)

def extractPressureHO(tc):
    tp=P_IBM.extractPressureHO(tc)
    return tp

def extractPressureHO2(tc):
    tp=P_IBM.extractPressureHO2(tc)
    return tp

def extractConvectiveTerms(tc):
    return P_IBM.extractConvectiveTerms(tc)

def unsteadyLoads(tb, Sref=None, Pref=None, Qref=None, alpha=0., beta=0.):
    return P_IBM.unsteadyLoads(tb, Sref=Sref, Pref=Pref, Qref=Qref, alpha=alpha, beta=beta)

def _unsteadyLoads(tb, Sref=None, Pref=None, Qref=None, alpha=0., beta=0.):
    return P_IBM._unsteadyLoads(tb, Sref=Sref, Pref=Pref, Qref=Qref, alpha=alpha, beta=beta)

def _prepareSkinReconstruction(ts, tc):
    return P_IBM._prepareSkinReconstruction(ts,tc)

def _computeSkinVariables(ts, tc, tl, graphWPOST):
    return P_IBM._computeSkinVariables(ts, tc, tl, graphWPOST)

def _modifIBCD(tc):
    raise NotImplementedError("_modifyIBCD is obsolete. Use _initOutflow and _initInj functions.")

#====================================================================================

## IMPORTANT NOTE !!
## FUNCTIONS MIGRATED TO $CASSIOPEE/Cassiopee/Generator/Generator/IBMmodelHeight.py
## The functions below will become deprecated after Jan. 1 2023
#====================================================================================
def compute_Cf(Re, Cf_law='ANSYS'):
    val=G_IBM_Height.compute_Cf(Re, Cf_law=Cf_law)
    return val

def computeYplusOpt(Re=None,tb=None,Lref=1.,q=1.2,snear=None,Cf_law='ANSYS'):
    val=G_IBM_Height.computeYplusOpt(Re=Re,tb=tb,Lref=Lref,q=q,snear=snear,Cf_law=Cf_law)
    return val

def computeSnearOpt(Re=None,tb=None,Lref=1.,q=1.2,yplus=300.,Cf_law='ANSYS'):
    val=G_IBM_Height.computeSnearOpt(Re=Re,tb=tb,Lref=Lref,q=q,yplus=yplus,Cf_law=Cf_law)
    return val


def dist2wallNearBody(t, tb, type='ortho', signed=0, dim=3, loc='centers'):
    list_final_zones=[]
    for z in Internal.getZones(t):
        list_final_zones.append(z[0])

    tBB =G.BB(t)
    tbBB=G.BB(tb)

    interDict = X.getIntersectingDomains(tBB, tbBB)

    #FULL TB
    zt       = []
    zt_names = []
    for i in interDict:
        if interDict[i]:
            zt.append(Internal.getNodeByName(t,i))
            zt_names.append(i)

    if zt_names:
        DTW._distance2Walls(zt, tb, type=type, signed=signed, dim=dim, loc=loc)

    ##PRT1
    list_additional_zones = get_zones_scale_up_down(tbBB,tBB,zt_names,dim=dim)

    ###PRT2
    if list_additional_zones:
        zt=[]
        for i in list_additional_zones:
            zt.append(Internal.getNodeByName(t,i))

        DTW._distance2Walls(zt, tb, type=type, signed=signed, dim=dim, loc=loc)
    return t


def get_zones_scale_up_down(tbBB,tBB,zt_names,diff_percent=0.15,sweep_num=4,scaleDirection=0,dim=2):
    minval_tb = C.getMinValue(tbBB, ['CoordinateX', 'CoordinateY','CoordinateZ']);
    maxval_tb = C.getMaxValue(tbBB, ['CoordinateX', 'CoordinateY','CoordinateZ']);
    mean_tb   = get_mean(maxval_tb,minval_tb)
    diff_percentz=diff_percent
    if dim==2: diff_percentz=0

    list_additional_zones=[]
    for i in range(1,sweep_num+1):
        if scaleDirection>=0:
            tbBB_scale    = T.scale(tbBB, factor=(1.0+i*diff_percent,1.0+i*diff_percent,1.0+i*diff_percentz))
            add2listAdditionalZones(list_additional_zones,tbBB_scale,tBB,mean_tb,zt_names)

        if scaleDirection<=0:
            tbBB_scale    = T.scale(tbBB, factor=(1.0-i*diff_percent,1.0-i*diff_percent,1.0-i*diff_percentz))
            add2listAdditionalZones(list_additional_zones,tbBB_scale,tBB,mean_tb,zt_names)


    return list_additional_zones


def get_mean(max_local,min_local):
    mean_local=[]
    for i in range(len(max_local)):
        mean_local.append((max_local[i]+min_local[i])/2)
    return mean_local


def add2listAdditionalZones(list_additional_zones,tbBB_scale,tBB,mean_tb,zt_names):
    minval_tbscale = C.getMinValue(tbBB_scale, ['CoordinateX', 'CoordinateY','CoordinateZ']);
    maxval_tbscale = C.getMaxValue(tbBB_scale, ['CoordinateX', 'CoordinateY','CoordinateZ']);
    mean_tbscale   = get_mean(maxval_tbscale,minval_tbscale)
    T._translate(tbBB_scale, (mean_tb[0]-mean_tbscale[0],mean_tb[1]-mean_tbscale[1],mean_tb[2]-mean_tbscale[2]))
    interDict_scale = X.getIntersectingDomains(tBB, tbBB_scale)
    for i in interDict_scale:
        if interDict_scale[i] and i not in list_additional_zones and i not in zt_names:
            list_additional_zones.append(i)
    return None
