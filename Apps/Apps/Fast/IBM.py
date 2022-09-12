# Class for FastS "IBM" prepare and compute
import FastC.PyTree as FastC
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

import timeit
import getpass
import socket
import os


try: range = xrange
except: pass
from mpi4py import MPI
COMM_WORLD = MPI.COMM_WORLD
KCOMM = COMM_WORLD

RENAMEIBCNODES=False # renommage des ibcd*

__IBCNameServer__={}

def getIBCDName(proposedName):
    global __IBCNameServer__
    (ibcname,__IBCNameServer__)=C.getUniqueName(proposedName, __IBCNameServer__)
    return ibcname


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
            frontType=1, expand=3, tinit=None, initWithBBox=-1., wallAdapt=None, dfarDir=0, recomputeDist=True):
    rank = Cmpi.rank; size = Cmpi.size
    ret = None
    # sequential prep
    if size == 1: ret = prepare0(t_case, t_out, tc_out, snears=snears, dfar=dfar, dfarList=dfarList,
                                 tbox=tbox, snearsf=snearsf, yplus=yplus,
                                 vmin=vmin, check=check, NP=NP, format=format, frontType=frontType, recomputeDist=recomputeDist,
                                 expand=expand, tinit=tinit, initWithBBox=initWithBBox, wallAdapt=wallAdapt, dfarDir=dfarDir)
    # parallel prep
    else: ret = prepare1(t_case, t_out, tc_out, snears=snears, dfar=dfar, dfarList=dfarList,
                         tbox=tbox, snearsf=snearsf, yplus=yplus, 
                         vmin=vmin, check=check, NP=NP, format=format, frontType=frontType,
                         expand=expand, tinit=tinit, initWithBBox=initWithBBox, wallAdapt=wallAdapt, dfarDir=dfarDir)

    return ret

#================================================================================
# IBM prepare - seq
#================================================================================
def prepare0(t_case, t_out, tc_out, snears=0.01, dfar=10., dfarList=[],
             tbox=None, snearsf=None, yplus=100.,
             vmin=21, check=False, NP=0, format='single', frontType=1, recomputeDist=True,
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
    t = G_IBM.generateIBMMesh(tb, vmin=vmin, snears=snears, dfar=dfar, dfarList=dfarList, DEPTH=2,
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
    t,tc = X_IBM.prepareIBMData(t, tb, frontType=frontType, interpDataType=0, yplus=yplus, wallAdapt=wallAdapt)
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


def generateCartesian(tb, to=None, dimPb=3, snears=0.01, dfar=10., dfarList=[], tbox=None, ext=3, snearsf=None, yplus=100.,
                      vmin=21, check=False, expand=3, dfarDir=0, extrusion=False,dz_in=0.01,NPas_in=200,span_in=1,check2Donly=False,dict_Nz={},isCartesianExtrude=False,isExtrudeByZone=False,directory_tmp_files='./'):
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

    if rank==0 and check: C.convertPyTree2File(o, directory_tmp_files+fileout)
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
    coords = C.getFields(Internal.__GridCoordinates__, zones, api=2)
    coords, rinds = Generator.extendCartGrids(coords, ext=ext, optimized=1, extBnd=0)
    C.setFields(coords, zones, 'nodes')
    for noz in range(len(zones)):
        Internal.newRind(value=rinds[noz], parent=zones[noz])
    Cmpi._rmXZones(t)
    coords = None; zones = None
    test.printMem(">>> extended cart grids (after rmXZones) [end]")

    if not extrusion:
        X_IBM._addBCOverlaps(t, bbox=bb)
        X_IBM._addExternalBCs(t, bbox=bb, dimPb=dimPb)

    if check and check2Donly:Cmpi.convertPyTree2File(t, directory_tmp_files+'2Dmesh.cgns')
    if check2Donly: exit()
        
    dz = dz_in
    if dimPb == 2:
        if not extrusion:
            T._addkplane(t)
            T._contract(t, (0,0,0), (1,0,0), (0,1,0), dz)
        if extrusion:
            NSplit = 1
            NPas   = NPas_in
            span   = span_in
            dimPb = 3
            # Extrude 2D case
            T._addkplane(tb,N=NPas+4)
            for node in Internal.getNodesFromName(tb,'EquationDimension'): Internal.setValue(node,dimPb)
            T._contract(tb, (0.,0.,0.), (1,0,0), (0,1,0), span/NPas)
            zmax = C.getMaxValue(tb,'CoordinateZ')
            T._translate(tb,(0.,0.,-0.5*zmax))

                        
            # Close new 3D case
            for b in Internal.getBases(tb):
                name = Internal.getName(b)
                b = C.convertArray2Tetra(b)
                b = G.close(b)
                b = P.exteriorFaces(b)
                b = T.splitConnexity(b)
                for line in Internal.getZones(b):
                    closure = G.tetraMesher(line, algo=1)
                    tb = Internal.append(tb, closure, name)

            ##[AJ] A random unstructured mesh was created (NOT SURE WHY)
            for z in Internal.getZones(tb):
                ztype = Internal.getNodeFromType2(z, 'ZoneType_t');
                if Internal.getValue(ztype)!='Structured':
                    Internal._rmNode(tb,z)

                    
            if rank == 0: C.convertPyTree2File(tb, directory_tmp_files+'3Dcase.cgns')
            
            # create new 3D tree
            t = T.subzone(t, (1,1,1), (-1,-1,1))
            bbox = G.bbox(t); bbox = [round(i,1) for i in bbox]
            bbox = numpy.array(bbox)
            
            # Share the boundaries of the entire mesh for BCFarfield
            comm.Barrier()
            minbox = numpy.zeros(3)
            maxbox = numpy.zeros(3)
            comm.Allreduce([bbox[0:3], MPI.DOUBLE], [minbox, MPI.DOUBLE], MPI.MIN)
            comm.Allreduce([bbox[3:], MPI.DOUBLE], [maxbox, MPI.DOUBLE], MPI.MAX)
            comm.Barrier()
            bbox[0:3] = minbox
            bbox[3:]  = maxbox
            C._rmBCOfType(t, 'BCFarfield')
            C._rmBCOfType(t, 'BCOverlap')
            Internal._rmNodesByType(t,'FlowSolution_t')
            for z in Internal.getZones(t):
                xmin = C.getValue( z, 'CoordinateX', (1,1,1))
                xmax = C.getValue( z, 'CoordinateX', (0,1,1))
                ymin = C.getValue( z, 'CoordinateY', (1,1,1))
                ymax = C.getValue( z, 'CoordinateY', (1,0,1))
                if abs(round(xmin-bbox[0]))==0.: C._addBC2Zone(z, 'external', 'BCFarfield', 'imin')
                if abs(round(xmax-bbox[3]))==0.: C._addBC2Zone(z, 'external', 'BCFarfield', 'imax')
                if abs(round(ymin-bbox[1]))==0.: C._addBC2Zone(z, 'external', 'BCFarfield', 'jmin')
                if abs(round(ymax-bbox[4]))==0.: C._addBC2Zone(z, 'external', 'BCFarfield', 'jmax')
            C._fillEmptyBCWith(t,'overlap','BCOverlap')

            if isCartesianExtrude:
                for z in Internal.getZones(t):
                    h = abs(C.getValue(z,'CoordinateX',0)-C.getValue(z,'CoordinateX',1))
                    NPas_local = int(round(span/h))
                    if NPas_local<2:
                        print("WARNING:: Zone %s has Nz=%d and is being clipped to Nz=2"%(z[0],NPas_local))
                        NPas_local=2
                    T._addkplane(z,N=NPas_local+4)
                    T._contract(z, (0.,0.,0.), (1,0,0), (0,1,0), span/NPas_local)
                    zmax_local = C.getMaxValue(z, 'CoordinateZ')
                    T._translate(z,(0.,0.,-0.5*zmax_local))
                    
            if isExtrudeByZone and len(dict_Nz)>0:
                for z in Internal.getZones(t):
                    name_zone = z[0][0:len(z[0])-2]
                    NPas_local = int(dict_Nz[name_zone])
                    T._addkplane(z,N=NPas_local+4)                    
                    T._contract(z, (0.,0.,0.), (1,0,0), (0,1,0), span/NPas_local)
                    zmax_local = C.getMaxValue(z, 'CoordinateZ')
                    T._translate(z,(0.,0.,-0.5*zmax_local))
                    
            if not isCartesianExtrude and not isExtrudeByZone:
                T._addkplane(t,N=NPas+4)
                T._contract(t, (0.,0.,0.), (1,0,0), (0,1,0), span/NPas)
                T._translate(t,(0.,0.,-0.5*zmax))
                
            for node in Internal.getNodesFromName(t,'EquationDimension'): Internal.setValue(node,3)
            C._addBC2Zone(t, 'period', 'BCautoperiod', 'kmin')
            C._addBC2Zone(t, 'period', 'BCautoperiod', 'kmax')
            
            if check: Cmpi.convertPyTree2File(t, directory_tmp_files+'3Dmesh.cgns')

    
    # ReferenceState
    C._addState(t, state=refstate)
    C._addState(t, 'GoverningEquations', model)
    C._addState(t, 'EquationDimension', dimPb)            
    return t,tb

#================================================================================
# IBM prepare - para
#
# extrusion: make an extrusion from a 2D profile. ATTENTION, each zone of the profile must be joined in one single zone
# smoothing : smooth the front during the front 2 specific treatment in the cases of local refinements
# balancing ; balance the entire distribution after the octree generation, useful for symetries
# distrib : new distribution at the end of prepare1
#===================================================================================================================
def prepare1(t_case, t_out, tc_out, t_in=None, to=None, snears=0.01, dfar=10., dfarList=[],
             tbox=None, snearsf=None, yplus=100., Lref=1.,
             vmin=21, check=False, NP=0, format='single',
             frontType=1, extrusion=False, smoothing=False, balancing=False,
             distrib=True, expand=3, tinit=None, initWithBBox=-1., wallAdapt=None, yplusAdapt=100.,
             dfarDir=0,dz_in=0.01,span_in=0.25,NPas_in=10,height_in=0.1,
             correctionMultiCorpsF42=False, blankingF42=False, twoFronts=False,
             redistribute=False,isDist2WallNearBodyOnly=False,isoverideheight=False,check2Donly=False,
             dict_Nz={},isCartesianExtrude=False,isExtrudeByZone=False,directory_tmp_files='./'):
    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = t_case


    
    if isCartesianExtrude and isExtrudeByZone:
        FastC._print2screen("isCartesianExtrude and isExtrudeByZone are both True. ONLY one can be true.",0)
        FastC._print2screen("Exiting...",1)
    if isExtrudeByZone and len(dict_Nz)==0:
        FastC._print2screen("isExtrudeByZone is True but the Nz per zone dictionary is empty.",0)
        FastC._print2screen("Exiting...",1)
    
    rank = Cmpi.rank
    comm = Cmpi.COMM_WORLD

    DEPTH=2
    IBCType=1

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

    if dimPb == 2: C._initVars(tb, 'CoordinateZ', 0.) # forced
    if t_in is None:
        t,tb = generateCartesian(tb, to=to, dimPb=dimPb, snears=snears, dfar=dfar, dfarList=dfarList, tbox=tbox, ext=DEPTH+1,
                              snearsf=snearsf, yplus=yplus,vmin=vmin, check=check, expand=expand, dfarDir=dfarDir, extrusion=extrusion,
                                 dz_in=dz_in,NPas_in=NPas_in,span_in=span_in,check2Donly=check2Donly,dict_Nz=dict_Nz,isCartesianExtrude=isCartesianExtrude,isExtrudeByZone=isExtrudeByZone,directory_tmp_files=directory_tmp_files)                    
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

    # Distance a la paroi
    test.printMem(">>> Wall distance [start]")
    FSC = Internal.getNodeFromType(t,"FlowSolution_t")
    
    if FSC is None or Internal.getNodeFromName(FSC,'TurbulentDistance') is None:
        if dimPb == 2:
            z0 = Internal.getNodeFromType2(t, "Zone_t")
            bb0 = G.bbox(z0); dz = bb0[5]-bb0[2]
            tb2 = C.initVars(tb, 'CoordinateZ', dz*0.5)
            tbsave =tb2
            t=dist2wallNearBody(t, tb2, type='ortho', signed=0, dim=dimPb, loc='centers')
        else:
            t=dist2wallNearBody(t, tb, type='ortho', signed=0, dim=dimPb, loc='centers')

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
    t = X_IBM.blankByIBCBodies(t, tb, 'centers', dimPb)
    C._initVars(t, '{centers:cellNIBC}={centers:cellN}')

    X_IBM._signDistance(t)

    C._initVars(t,'{centers:cellN}={centers:cellNIBC}')
    # determination des pts IBC
    Reynolds = Internal.getNodeFromName(tb, 'Reynolds')
    if Reynolds is not None: Reynolds = Internal.getValue(Reynolds)
    if Reynolds < 1.e5: frontType = 1
    if frontType != 42:
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
                
            print("height %f & heigh/h %f ="%(height,height/h))
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

    X_IBM._removeBlankedGrids(t, loc='centers')
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

    # setInterpData - Chimere
    C._initVars(t,'{centers:cellN}=maximum(0.,{centers:cellNChim})')# vaut -3, 0, 1, 2 initialement

    # maillage donneur: on MET les pts IBC comme donneurs
    # tp = Internal.copyRef(t)
    # FSN = Internal.getNodesFromName3(tp, Internal.__FlowSolutionNodes__)
    # Internal._rmNodesByName(FSN, 'cellNFront')
    # Internal._rmNodesByName(FSN, 'cellNIBC')
    # Internal._rmNodesByName(FSN, 'TurbulentDistance')
    # tc = C.node2Center(tp); del tp

    if extrusion:
        listvars_local =['cellN','cellNChim','cellNIBC','cellNFront']
        for z in Internal.getZones(t):
            sol            = Internal.getNodeFromName(z,'FlowSolution#Centers')
            for i in listvars_local:
                cellN          = Internal.getNodeFromName(sol,i)[1]
                sh             = numpy.shape(cellN)
                if i== 'cellNChim':
                    cellN[:,:,0:2]            = -3
                    cellN[:,:,sh[2]-2:sh[2]]  = -3
                else:
                    cellN[:,:,0:2]            = 0
                    cellN[:,:,sh[2]-2:sh[2]]  = 0

    test.printMem(">>> Interpdata [start]")
    tc = C.node2Center(t)

    # abutting ? 
    if Internal.getNodeFromType(t,"GridConnectivity1to1_t") is not None:
        test.printMem("setInterpData abutting")
        Xmpi._setInterpData(t, tc, 
                            nature=1, loc='centers', storage='inverse', 
                            sameName=1, dim=3, itype='abutting')
        test.printMem("setInterpData abutting done.")

    # setInterpData parallel pour le chimere
    tbbc = Cmpi.createBBoxTree(tc)
    interDict = X.getIntersectingDomains(tbbc)
    graph = Cmpi.computeGraph(tbbc, type='bbox', intersectionsDict=interDict, reduction=False)
    Cmpi._addXZones(tc, graph, variables=['cellN'], cartesian=True)
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
                         sameName=1, interpDataType=0, itype='chimera')
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
        Cmpi._addLXZones(tc, graph, variables=['cellNIBC','cellNChim','cellNFront'], cartesian=True, interDict=interDict, bboxDict=bboxDict, layers=4, subr=False)
        Cmpi._addLXZones(t, graph, variables=['centers:cellNIBC', 'centers:cellNChim', 'centers:cellNFront'], cartesian=True, interDict=interDict, bboxDict=bboxDict, layers=4, subr=False)

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
                                fields = X.transferFields(zc, XI, YI, ZI, hook=None, variables=['cellNFront_origin','cellNIBC_origin'], interpDataType=0, nature=1)
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
        C.convertPyTree2File(front, directory_tmp_files+'front.cgns')
        if twoFronts: C.convertPyTree2File(front2, directory_tmp_files+'front2.cgns')

    zonesRIBC = []
    for zrcv in Internal.getZones(t):
        if C.getMaxValue(zrcv, 'centers:cellNIBC')==2.:
            zrcvname = zrcv[0]; zonesRIBC.append(zrcv)

    nbZonesIBC = len(zonesRIBC)
    if nbZonesIBC == 0:
        res = [{},{},{}]
        if twoFronts: res2 = [{},{},{}]
    else:
        res = X_IBM.getAllIBMPoints(zonesRIBC, loc='centers',tb=tb, tfront=front, frontType=frontType,
                                   cellNName='cellNIBC', depth=DEPTH, IBCType=IBCType, Reynolds=Reynolds, yplus=yplus, Lref=Lref)
        if twoFronts:
            res2 = X_IBM.getAllIBMPoints(zonesRIBC, loc='centers',tb=tb, tfront=front2, frontType=frontType,
                                        cellNName='cellNIBC', depth=DEPTH, IBCType=IBCType, Reynolds=Reynolds, yplus=yplus, Lref=Lref)

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
    if twoFronts:
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
                        bba = C.getFields('GridCoordinates', z)[0]
                        if Generator.bboxIntersection(interpPtsBB,bba,isBB=True):
                            zname = z[0]
                            popp = Cmpi.getProc(z)
                            Distributed.updateGraph__(graph, popp, rank, zname)
                            if zrname not in interDictIBM: interDictIBM[zrname]=[zname]
                            else:
                                if zname not in interDictIBM[zrname]: interDictIBM[zrname].append(zname)
        if twoFronts:
            for ibcTypeL in dictOfCorrectedPtsByIBCType2:
                    allCorrectedPts2 = dictOfCorrectedPtsByIBCType2[ibcTypeL]
                    allWallPts2 = dictOfWallPtsByIBCType2[ibcTypeL]
                    allInterpPts2 = dictOfInterpPtsByIBCType2[ibcTypeL]
                    for nozr in range(nbZonesIBC):
                        if allCorrectedPts2[nozr] != []:
                            zrname = zonesRIBC[nozr][0]
                            interpPtsBB2 = Generator.BB(allInterpPts2[nozr])
                            for z in zones:
                                bba = C.getFields('GridCoordinates', z)[0]
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
    Cmpi._addXZones(tc, graph, variables=['cellN'], cartesian=True, subr=False)
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
                                             interpDataType=0, ReferenceState=ReferenceState, bcType=ibcTypeL)

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
                                                     interpDataType=0, ReferenceState=ReferenceState, bcType=ibcTypeL)

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
    if twoFronts:
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


        if twoFronts:
            tibm2 = X_IBM.extractIBMInfo2(tc)

            # Avoid that two procs write the same information
            for z in Internal.getZones(tibm2):
               if int(z[0][-1]) != rank:
                  # Internal._rmNodesByName(tibm, z[0])
                  z[0] = z[0]+"%{}".format(rank)

            Cmpi.convertPyTree2File(tibm2, 'IBMInfo2.cgns')

        test.printMem(">>> Saving IBM infos [end]")
        del tibm
        if twoFronts: del tibm2

    # distribution par defaut (sur NP)
    tbbc = Cmpi.createBBoxTree(tc)

    # Perform the final distribution
    if distrib:
        if NP == 0: NP = Cmpi.size
        stats = D2._distribute(tbbc, NP, algorithm='graph', useCom='ID')
        D2._copyDistribution(tc, tbbc)
        D2._copyDistribution(t, tbbc)

    del tbbc

    
    tmpFilename = 'tmp_tc_IBM_prep'    
    localUser= getpass.getuser()
    localPc  = socket.gethostname().split('-')
    if 'sator' not in localPc:
        path2writetmp='/stck/'+localUser+'/'
    else:
        path2writetmp='/tmp_user/sator/'+localUser+'/'
    tmpFilename = path2writetmp+tmpFilename+'.cgns'
    
    if redistribute:
        # Distribute over NP procs
        if rank == 0: print("REDISTRIBUTE - Final")
        tcp = Compressor.compressCartesian(tc)
        Cmpi.convertPyTree2File(tcp, tmpFilename, ignoreProcNodes=True)
        del tcp

        tc    = Cmpi.convertFile2PyTree(tmpFilename)
        stats = D2._distribute(tc, NP, algorithm='graph', useCom='ID')
        D2._copyDistribution(t, tc)    
        if rank==0 and os.path.exists(tmpFilename): os.remove(tmpFilename)

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
    if twoFronts:
        tc2 = Internal.copyTree(tc)
        tc2 = Internal.rmNodesByName(tc2, 'IBCD*')
        tc  = Internal.rmNodesByName(tc, '2_IBCD*')

    if RENAMEIBCNODES:
        for zc in Internal.getZones(tc):
            for ibcd in Internal.getNodesFromName1(zc,'IBCD_*'):            
                proposedName = Internal.getName(ibcd)[0:6]+'_X%d'%(rank)
                ibcd[0]=getIBCDName(proposedName)

        if twoFronts:
            for zc in Internal.getZones(tc2):
                for ibcd in Internal.getNodesFromName1(zc,'2_IBCD_*'):            
                    proposedName = Internal.getName(ibcd)[0:8]+'_X%d'%(rank)
                    ibcd[0]=getIBCDName(proposedName)

    if isinstance(tc_out, str):
        tcp = Compressor.compressCartesian(tc)
        Cmpi.convertPyTree2File(tcp, tc_out, ignoreProcNodes=True)

        if twoFronts:
            tc2 = transformTc2(tc2)
            tcp2 = Compressor.compressCartesian(tc2)
            Cmpi.convertPyTree2File(tcp2, 'tc2.cgns', ignoreProcNodes=True)
            del tc2

    # Initialisation
    if tinit is None: I._initConst(t, loc='centers')
    else:
        t = Pmpi.extractMesh(tinit, t, mode='accurate')
    if model != "Euler": C._initVars(t, 'centers:ViscosityEddy', 0.)

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

    # Save t    
    if not redistribute:
        tc_local = Internal.rmNodesByName(tc,'GridCoordinates')
        Internal._rmNodesByName(tc_local,'ID*')
        Internal._rmNodesByName(tc_local,'IBCD*')
        Cmpi.convertPyTree2File(tc_local, tmpFilename, ignoreProcNodes=True)
        if rank ==0:
            tc_local= C.convertFile2PyTree(tmpFilename)
            checkNcellsNptsPerProc(tc_local,Cmpi.size,isAtCenter=True)
            del tc_local
            if rank==0 and os.path.exists(tmpFilename): os.remove(tmpFilename)
    else:
        if rank ==0:checkNcellsNptsPerProc(tc,Cmpi.size,isAtCenter=True)
     
    if isinstance(t_out, str):
        tp = Compressor.compressCartesian(t)
        Cmpi.convertPyTree2File(tp, t_out, ignoreProcNodes=True)

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
    print(stats)
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
                if isAtCenter:
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
                initWithBBox=-1., wallAdapt=None,dfarDir=0):
        if NP is None: NP = Cmpi.size
        if NP == 0: print('Preparing for a sequential computation.')
        else: print('Preparing for an IBM computation on %d processors.'%NP)
        ret = prepare(t_case, t_out, tc_out, snears=snears, dfar=dfar, dfarList=dfarList,
                      tbox=tbox, snearsf=snearsf, yplus=yplus,
                      vmin=vmin, check=check, NP=NP, format=self.data['format'],
                      frontType=frontType, expand=expand, tinit=tinit, dfarDir=dfarDir)
        return ret

    # post-processing: extrait la solution aux noeuds + le champs sur les surfaces
    def post(self, t_case, t_in, tc_in, t_out, wall_out):
        return post(t_case, t_in, tc_in, t_out, wall_out)

    # post-processing: extrait les efforts sur les surfaces
    def loads(self, t_case, tc_in=None, wall_out=None, alpha=0., beta=0., Sref=None, famZones=[]):
        return loads(t_case, tc_in=tc_in, wall_out=wall_out, alpha=alpha, beta=beta, Sref=Sref, famZones=famZones)



## IMPORTANT NOTE !!
## FUNCTIONS MIGRATED TO $CASSIOPEE/Apps/Modules/Geom/Geom/IBM.py
## The functions below will become decrepit after Jan. 1 2023
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
## FUNCTIONS MIGRATED TO $CASSIOPEE/Apps/Modules/Post/Post/IBM.py
## The functions below will become decrepit after Jan. 1 2023
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
## FUNCTIONS MIGRATED TO $CASSIOPEE/Apps/Modules/Generator/Generator/IBMmodelHeight.py
## The functions below will become decrepit after Jan. 1 2023
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





