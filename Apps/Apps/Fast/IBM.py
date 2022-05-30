# Class for FastS "IBM" prepare and compute
import FastC.PyTree as FastC
import Converter.PyTree as C
import Geom.PyTree as D
import Generator.PyTree as G
import Transform.PyTree as T
import Post.PyTree as P
import Converter.Internal as Internal
import Connector.PyTree as X
import Connector.ToolboxIBM as TIBM
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
try: range = xrange
except: pass
from mpi4py import MPI
COMM_WORLD = MPI.COMM_WORLD
KCOMM = COMM_WORLD

def compute_Cf(Re, Cf_law='ANSYS'):
    if Cf_law == 'ANSYS':
        return 0.058*Re**(-0.2)
    elif Cf_law == 'PW':
        return 0.026*Re**(-1/7.)
    elif Cf_law == 'PipeDiameter':
        return 0.079*Re**(-0.25)
    elif Cf_law == 'Laminar':
        return 1.328*Re**(-0.5)

def computeYplusOpt(Re=None,tb=None,Lref=1.,q=1.2,snear=None,Cf_law='ANSYS'):
    fail=0
    if Re is None:
        if tb is not None:
            Re = Internal.getNodeFromName(tb,"Reynolds")
            if Re is None: fail=1
            else:
                Re = Internal.getValue(Re)
        else: fail = 1
    if fail: 
        raise ValueError("computeYplusOpt: requires Reynolds number as a float or in tb.")
    fail = 0
    if snear is None:
        snear = Internal.getNodeFromName(tb,"snear")
        if snear is None: fail=1
        else: snear = Internal.getValue(snear)
    if fail:
        raise ValueError("computeYlusOpt: requires snear as a float or in tb.")

    print("Warning: estimation of the optimum y+ at Reynolds number ", Re, " and snear target at image point ", snear)
    h0 = (1.*Lref*math.sqrt(2.))/(Re*math.sqrt(compute_Cf(Re,Cf_law))) #Taille de maille pour y+1
    h_opti = (h0-q*snear)/(1.-q) #Hauteur de modelisation opti
    yplus_opti = h_opti/h0 #yplus opti

    # print('\nInformation for the body-fitted mesh :')
    # print('h_opti     = {:.2e}'.format(h_opti))
    # print('h0         = {:.2e}\n'.format(h0))
    # print('Information for the Cartesian mesh :')
    # print('yplus_opti = {}\n'.format(math.ceil(yplus_opti)))
    return yplus_opti

# compute the near wall spacing in agreement with the yplus target at image points - front42
def computeSnearOpt(Re=None,tb=None,Lref=1.,q=1.2,yplus=300.,Cf_law='ANSYS'):
    fail=0
    if Re is None:
        if tb is not None:
            Re = Internal.getNodeFromName(tb,"Reynolds")
            if Re is None: fail=1
            else: Re = Internal.getValue(Re)
        else: fail = 1
    if fail: 
        raise ValueError("computeSnearOpt: requires Reynolds number as a float or in tb.")


    print("Estimation of the optimum near-wall spacing at Reynolds number ", Re, " and yplus target at image point ", yplus)
    h_mod = (yplus*Lref*math.sqrt(2.))/(Re*math.sqrt(compute_Cf(Re,Cf_law)))
    h0    = (Lref*math.sqrt(2.))/(Re*math.sqrt(compute_Cf(Re,Cf_law))) #Taille de maille pour y+=1
    n     = int(math.ceil(math.log(1-yplus*(1-q))/math.log(q))) # number of cells in the BF mesh for the height h
    snear_opti = q**(n-1)*h0 # best snear for the target yplus
    print('\nInformation for the body-fitted mesh :')
    print('h           = {:.2e}'.format(h_mod))
    print('h0          = {:.2e}\n'.format(h0))
    print('Information for the Cartesian mesh :')
    print('snear_opti  = {:.3e}\n'.format(snear_opti))
    return snear_opti

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
                         vmin=vmin, check=check, NP=NP, format=format, frontType=frontType, recomputeDist=recomputeDist,
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
    t = TIBM.generateIBMMesh(tb, vmin=vmin, snears=snears, dfar=dfar, dfarList=dfarList, DEPTH=2,
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
    t,tc = TIBM.prepareIBMData(t, tb, frontType=frontType, interpDataType=0, yplus=yplus, wallAdapt=wallAdapt)
    test.printMem(">>> ibm data [end]")

    # arbre donneur
    D2._copyDistribution(tc, t)
    if isinstance(tc_out, str): FastC.save(tc, tc_out, split=format, NP=-NP)

    #----------------------------------------
    # Extraction des coordonnees des pts IBM
    #----------------------------------------
    if check:
        tibm = TIBM.extractIBMInfo(tc)
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

def generateCartesian(tb, dimPb=3, snears=0.01, dfar=10., dfarList=[], tbox=None, ext=3, snearsf=None, yplus=100.,
                      vmin=21, check=False, expand=3, dfarDir=0, extrusion=False):
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

    # a mettre dans la classe ou en parametre de prepare1 ???
    to = None

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

    o = TIBM.buildOctree(tb, snears=snears, snearFactor=1., dfar=dfar, dfarList=dfarList,
                         to=to, tbox=tbox, snearsf=snearsf,
                         dimPb=dimPb, vmin=vmin, symmetry=symmetry, fileout=None, rank=rank,
                         expand=expand, dfarDir=dfarDir)

    if rank==0 and check: C.convertPyTree2File(o, fileout)
    # build parent octree 3 levels higher
    # returns a list of 4 octants of the parent octree in 2D and 8 in 3D
    parento = TIBM.buildParentOctrees__(o, tb, snears=snears, snearFactor=4., dfar=dfar, dfarList=dfarList, to=to, tbox=tbox, snearsf=snearsf,
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
    res = TIBM.octree2StructLoc__(p, vmin=vmin, ext=-1, optimized=0, parento=parento, sizeMax=1000000)
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
        TIBM._addBCOverlaps(t, bbox=bb)
        TIBM._addExternalBCs(t, bbox=bb, dimPb=dimPb)

    dz = 0.01
    if dimPb == 2:
        if not extrusion:
            T._addkplane(t)
            T._contract(t, (0,0,0), (1,0,0), (0,1,0), dz)
        if extrusion:
            chord = 1.
            NSplit = 1
            NPas = 200
            span = 0.25*chord
            dimPb = 3
            # Extrude 2D case
            T._addkplane(tb,N=NPas+4)
            for node in Internal.getNodesFromName(tb,'EquationDimension'): Internal.setValue(node,3)
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
            if rank == 0: C.convertPyTree2File(tb, '3Dcase.cgns')
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
            T._addkplane(t,N=NPas+4)
            for node in Internal.getNodesFromName(t,'EquationDimension'): Internal.setValue(node,3)
            T._contract(t, (0.,0.,0.), (1,0,0), (0,1,0), span/NPas)
            T._translate(t,(0.,0.,-0.5*zmax))
            C._addBC2Zone(t, 'period', 'BCautoperiod', 'kmin')
            C._addBC2Zone(t, 'period', 'BCautoperiod', 'kmax')
            if check: Cmpi.convertPyTree2File(t, '3Dmesh.cgns')


    # ReferenceState
    C._addState(t, state=refstate)
    C._addState(t, 'GoverningEquations', model)
    C._addState(t, 'EquationDimension', dimPb)            
    return t

#================================================================================
# IBM prepare - para
#
# extrusion: make an extrusion from a 2D profile. ATTENTION, each zone of the profile must be joined in one single zone
# smoothing : smooth the front during the front 2 specific treatment in the cases of local refinements
# balancing ; balance the entire distribution after the octree generation, useful for symetries
# distrib : new distribution at the end of prepare1
#===================================================================================================================
def prepare1(t_case, t_out, tc_out, t_in=None, snears=0.01, dfar=10., dfarList=[],
             tbox=None, snearsf=None, yplus=100., Lref=1.,
             vmin=21, check=False, NP=0, format='single',
             frontType=1, extrusion=False, smoothing=False, balancing=False, recomputeDist=True,
             distrib=True, expand=3, tinit=None, initWithBBox=-1., wallAdapt=None, yplusAdapt=100., dfarDir=0, 
             correctionMultiCorpsF42=False, blankingF42=False, twoFronts=False):
    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = t_case

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
        t = generateCartesian(tb, dimPb=dimPb, snears=snears, dfar=dfar, dfarList=dfarList, tbox=tbox, ext=DEPTH+1,
                              snearsf=snearsf, yplus=yplus,vmin=vmin, check=check, expand=expand, dfarDir=dfarDir, extrusion=extrusion)                    
    else: 
        t = t_in

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
            DTW._distance2Walls(t, tb2, type='ortho', signed=0, dim=dimPb, loc='centers')
        else:
            DTW._distance2Walls(t, tb, type='ortho', signed=0, dim=dimPb, loc='centers')


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
            shiftDist = TIBM.computeModelisationHeight(Re=Reynolds, yplus=yplus, L=Lref)
        else:
            snears = Internal.getNodesFromName(tb, 'snear')
            h = max(snears, key=lambda x: x[1])[1]
            shiftDist = TIBM.computeBestModelisationHeight(Re=Reynolds, h=h) # meilleur compromis entre hauteur entre le snear et la hauteur de modelisation

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
    t = TIBM.blankByIBCBodies(t, tb, 'centers', dimPb)
    C._initVars(t, '{centers:cellNIBC}={centers:cellN}')

    TIBM._signDistance(t)

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
                height = TIBM.computeModelisationHeight(Re=Reynolds, yplus=yplus, L=Lref)
            else:
                height = TIBM.computeBestModelisationHeight(Re=Reynolds, h=h) # meilleur compromis entre hauteur entre le snear et la hauteur de modelisation
                yplus  = TIBM.computeYplus(Re=Reynolds, height=height, L=Lref)
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

    TIBM._removeBlankedGrids(t, loc='centers')
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
                    SHIFTB = TIBM.computeModelisationHeight(Re=Reynolds, yplus=yplus, L=Lref)
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
        if smoothing and dimPb == 2: TIBM._smoothImageFront(t, tc)

        C._cpVars(t,'centers:cellNFront',tc,'cellNFront')

        Xmpi._setInterpTransfers(t,tc,variables=['cellNFront'], cellNVariable='cellNFront', compact=0)
        test.printMem(">>> pushBackImageFront2 [end]")
    ############################################################

    C._rmVars(t,['centers:cellNFront'])
    if twoFronts:
        C._rmVars(t,['centers:cellNFront_2', 'centers:cellNIBC_2'])
    C._cpVars(t,'centers:TurbulentDistance',tc,'TurbulentDistance')

    print('Minimum distance: %f.'%C.getMinValue(t,'centers:TurbulentDistance'))
    P._computeGrad2(t, 'centers:TurbulentDistance',ghostCells=True)

    test.printMem(">>> Building IBM front [start]")
    front = TIBM.getIBMFront(tc, 'cellNFront', dim=dimPb, frontType=frontType)
    front = TIBM.gatherFront(front)

    if twoFronts:
        front2 = TIBM.getIBMFront(tc, 'cellNFront_2', dim=dimPb, frontType=frontType)
        front2 = TIBM.gatherFront(front2)

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
        res = TIBM.getAllIBMPoints(zonesRIBC, loc='centers',tb=tb, tfront=front, frontType=frontType,
                                   cellNName='cellNIBC', depth=DEPTH, IBCType=IBCType, Reynolds=Reynolds, yplus=yplus, Lref=Lref)
        if twoFronts:
            res2 = TIBM.getAllIBMPoints(zonesRIBC, loc='centers',tb=tb, tfront=front2, frontType=frontType,
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
    
    if model != 'Euler' and recomputeDist:
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


    # Sauvegarde des infos IBM
    if check:
        test.printMem(">>> Saving IBM infos [start]")
        tibm = TIBM.extractIBMInfo(tc)

        # Avoid that two procs write the same information
        for z in Internal.getZones(tibm):
           if int(z[0][-1]) != rank:
              # Internal._rmNodesByName(tibm, z[0])
              z[0] = z[0]+"%{}".format(rank)

        Cmpi.convertPyTree2File(tibm, 'IBMInfo.cgns')


        if twoFronts:
            tibm2 = TIBM.extractIBMInfo2(tc)

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

    # Save tc
    if twoFronts:
        tc2 = Internal.copyTree(tc)
        tc2 = Internal.rmNodesByName(tc2, 'IBCD*')
        tc  = Internal.rmNodesByName(tc, '2_IBCD*')

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
    if isinstance(t_out, str):
        tp = Compressor.compressCartesian(t)
        Cmpi.convertPyTree2File(tp, t_out, ignoreProcNodes=True)

    if Cmpi.size > 1: Cmpi.barrier()
    return t, tc

def extractIBMInfo(tc_in, t_out='IBMInfo.cgns'):
    if isinstance(tc_in, str): tc = Cmpi.convertFile2PyTree(tc_in)
    else: tc = tc_in

    tibm = TIBM.extractIBMInfo(tc)
    rank = Cmpi.rank
    Distributed._setProc(tibm,rank)
    if isinstance(t_out, str): Cmpi.convertPyTree2File(tibm, t_out)
    return tibm

#=============================================================================
# Post
# IN: t_case: fichier ou arbre du cas
# IN: t_in: fichier ou arbre de resultat
# IN: tc_in: fichier ou arbre de connectivite
# OUT: t_out ou None: fichier pour sortie du champ aux noeuds
# OUT: wall_out ou None: fichier pour sortie du champ sur la paroi
#==============================================================================
def post(t_case, t_in, tc_in, t_out, wall_out):
    if isinstance(t_in, str): t = C.convertFile2PyTree(t_in)
    else: t = t_in
    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = t_case

    #=============================
    # Supprime les champs inutiles
    #=============================
    vars = ['centers:Density_M1', 'centers:VelocityX_M1', 'centers:VelocityY_M1', 'centers:VelocityZ_M1', 'centers:Temperature_M1', 'centers:Density_P1', 'centers:VelocityX_P1', 'centers:VelocityY_P1', 'centers:VelocityZ_P1', 'centers:Temperature_P1','centers:TurbulentDistance']
    C._rmVars(t, vars)

    #=============================
    # Arbre de connectivite
    #=============================
    if isinstance(tc_in, str): tc = C.convertFile2PyTree(tc_in)
    else: tc = tc_in
    Internal._rmNodesByName(tc, 'GridCoordinates')

    #==========================================================
    # Extraction Cp, Cf, ... sur les surfaces par interpolation
    #==========================================================
    tb = C.convertArray2Tetra(tb)

    #--------------------------------------------------------
    # Get Reference State and model from body pyTree
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('GoverningEquations is missing in input cgns.')
    # model : Euler, NSLaminar, NSTurbulent
    model = Internal.getValue(model)

    if model == 'Euler': bcType = 0
    elif model =='NSLaminar': bcType = 1
    else: bcType = 3 # Musker

    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)

    [RoInf, RouInf, RovInf, RowInf, RoeInf, PInf, TInf, cvInf, MInf,
    ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
    Mus, Cs, Ts, Pr] = C.getState(tb)

    varType = 2 # IBM updated variables (rho,u,t)
    varsIBC = ['Density', 'VelocityX', 'VelocityY', 'VelocityZ', 'Temperature']
    vars = ['Density', 'VelocityX', 'VelocityY', 'VelocityZ', 'Temperature']
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
    zw = TIBM.extractIBMWallFields(tc, tb=tb)
    RoUInf2I = 1./(RouInf*RouInf+RovInf*RovInf+RowInf*RowInf)
    C._initVars(zw,'{Cp}=2*%f*({Pressure}-%f)*%f'%(RoInf,PInf,RoUInf2I))
    if model != 'Euler':
        C._initVars(zw,'{Cf}=2*%f*{Density}*{utau}**2*%f'%(RoInf,RoUInf2I))

    Internal._rmNodesByName(zw, '.Solver#Param')
    Internal._rmNodesByName(zw, '.Solver#ownData')

    if isinstance(wall_out, str): C.convertPyTree2File(zw, wall_out)

    #===============================
    # En 2D, extrait un seul plan k
    #================================
    if dimPb == 2:
        t = T.subzone(t, (1,1,1), (-1,-1,1))
        C._initVars(t, 'CoordinateZ', 0.) # forced

    #=================================
    # Calcul de mut/mu dans le volume
    #=================================
    if model != 'Euler':
        betas = Mus*(Ts+Cs)/(Ts**(3./2.))
        C._initVars(t,'{centers:ViscosityMolecular} = %20.16g*sqrt({centers:Temperature})/(1.+%20.16g/{centers:Temperature})'%(betas,Cs))
        C._initVars(t,'{centers:mutsmu}=({centers:ViscosityEddy})/({centers:ViscosityMolecular})-1.')

    #==============================
    # Sortie champs aux noeuds
    #==============================
    # vars = ['centers:Density','centers:VelocityX', 'centers:VelocityY', 'centers:VelocityZ', 'centers:Temperature','centers:ViscosityEddy',
    # 'centers:TurbulentSANuTilde','centers:ViscosityMolecular', 'centers:mutsmu', 'centers:cellN']
    vars = ['centers:Density','centers:VelocityX', 'centers:Temperature','centers:ViscosityEddy',
    'centers:TurbulentSANuTilde','centers:ViscosityMolecular', 'centers:mutsmu', 'centers:cellN']
    t = C.center2Node(t, vars)
    Internal._rmNodesByName(t, 'FlowSolution#Centers')
    if isinstance(t_out, str): C.convertPyTree2File(t, t_out)

    return t, zw

#===========================================================
# compute [Cl, Cd]
# return wall pytree with Cp/Cf and gradtP/gradnP
# alpha, beta are angles in degrees
#===========================================================
def loads0(ts, Sref=None, alpha=0., beta=0., dimPb=3, verbose=False):
    if Sref is None:
        C._initVars(ts, '__ONE__',1.)
        Sref = P.integ(ts, '__ONE__')[0];
        C._rmVars(ts, ['__ONE__', 'centers:vol'])

    RefState = Internal.getNodeFromType(ts,'ReferenceState_t')
    PInf  = Internal.getValue(Internal.getNodeFromName(RefState,"Pressure"))
    RoInf = Internal.getValue(Internal.getNodeFromName(RefState,"Density"))
    VxInf = Internal.getValue(Internal.getNodeFromName(RefState,"VelocityX"))
    VyInf = Internal.getValue(Internal.getNodeFromName(RefState,"VelocityY"))
    VzInf = Internal.getValue(Internal.getNodeFromName(RefState,"VelocityZ"))
    VInf2 = VxInf*VxInf+VyInf*VyInf+VzInf*VzInf
    VInf  = math.sqrt(VInf2)

    q = 0.5*RoInf*VInf2
    qinv = 1./q
    alpha  = math.radians(alpha)
    beta   = math.radians(beta)
    calpha = math.cos(alpha); cbeta = math.cos(beta)
    salpha = math.sin(alpha); sbeta = math.sin(beta)
    #===========================
    # Calcul efforts de pression
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
    if dimPb==3:
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
    # Calcul frottement et efforts visqueux
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
    # calcul du vecteur tangent
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

    # calcul forces de frottement
    C._initVars(ts, '{centers:Fricx}={centers:tauxx}*{centers:sx}+{centers:tauxy}*{centers:sy}+{centers:tauxz}*{centers:sz}')
    C._initVars(ts, '{centers:Fricy}={centers:tauxy}*{centers:sx}+{centers:tauyy}*{centers:sy}+{centers:tauyz}*{centers:sz}')
    C._initVars(ts, '{centers:Fricz}={centers:tauxz}*{centers:sx}+{centers:tauyz}*{centers:sy}+{centers:tauzz}*{centers:sz}')

    # calcul coefficient de frottement
    C._initVars(ts, '{centers:Cf}=(sqrt({centers:Fricx}**2+{centers:Fricy}**2+{centers:Fricz}**2))/%g'%q)

    # maj des gradients de pression (norm/tang)
    if isGradP:
        C._initVars(ts, '{centers:gradnP}={centers:gradxPressure}*{centers:sx}+{centers:gradyPressure}*{centers:sy}+{centers:gradzPressure}*{centers:sz}')
        C._initVars(ts, '{centers:gradtP}={centers:gradxPressure}*{centers:tx}+{centers:gradyPressure}*{centers:ty}+{centers:gradzPressure}*{centers:tz}')

    G._getVolumeMap(ts)
    effortX = P.integ(ts, 'centers:Fricx')[0]
    effortY = P.integ(ts, 'centers:Fricy')[0]
    effortZ = P.integ(ts, 'centers:Fricz')[0]

    QADIMI = 1./(q*Sref)
    if dimPb==3:
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
# Post efforts
# IN: t_case: fichier ou arbre du cas
# IN: tc_in: fichier ou arbre de connectivite contenant les IBCD
# si tc_in =None, t_case est la surface avec la solution deja projetee
# OUT: wall_out ou None: fichier pour sortie des efforts sur la paroi aux centres
# IN: alpha: angle pour les efforts
# IN: beta: angle pour les efforts
#==============================================================================
def loads(t_case, tc_in=None, tc2_in=None, wall_out=None, alpha=0., beta=0., gradP=False, order=1, Sref=None, famZones=[]):
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
        for z in Internal.getZones(tc):
            subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
            for zsr in subRegions:
                nameSubRegion = zsr[0]
                if nameSubRegion[:4] == "IBCD":
                    pressure = Internal.getNodeFromName(zsr, 'Pressure')[1]
                    gradxP   = Internal.getNodeFromName(zsr, 'gradxPressure')
                    gradyP   = Internal.getNodeFromName(zsr, 'gradyPressure')
                    gradzP   = Internal.getNodeFromName(zsr, 'gradzPressure')
                    nIBC = pressure.shape[0]

                    if gradxP is  None:
                        gradxPressureNP = numpy.zeros((nIBC),numpy.float64)
                        gradyPressureNP = numpy.zeros((nIBC),numpy.float64)
                        gradzPressureNP = numpy.zeros((nIBC),numpy.float64)
                        zsr[2].append(['gradxPressure' , gradxPressureNP , [], 'DataArray_t'])
                        zsr[2].append(['gradyPressure' , gradyPressureNP , [], 'DataArray_t'])
                        zsr[2].append(['gradzPressure' , gradzPressureNP , [], 'DataArray_t'])

        if tc2 is not None: 
            # add gradP fields in tc2 if necessary
            for z in Internal.getZones(tc2):
                subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
                for zsr in subRegions:
                    nameSubRegion = zsr[0]
                    if nameSubRegion[:4] == "IBCD":
                        pressure = Internal.getNodeFromName(zsr, 'Pressure')[1]
                        gradxP   = Internal.getNodeFromName(zsr, 'gradxPressure')
                        gradyP   = Internal.getNodeFromName(zsr, 'gradyPressure')
                        gradzP   = Internal.getNodeFromName(zsr, 'gradzPressure')
                        nIBC = pressure.shape[0]

                        if gradxP is  None:
                            gradxPressureNP = numpy.zeros((nIBC),numpy.float64)
                            gradyPressureNP = numpy.zeros((nIBC),numpy.float64)
                            gradzPressureNP = numpy.zeros((nIBC),numpy.float64)
                            zsr[2].append(['gradxPressure' , gradxPressureNP , [], 'DataArray_t'])
                            zsr[2].append(['gradyPressure' , gradyPressureNP , [], 'DataArray_t'])
                            zsr[2].append(['gradzPressure' , gradzPressureNP , [], 'DataArray_t'])

        if tc2 is not None: 
            if order < 2:
                tc2 = extractPressureHO(tc2)
            else:
                tc2 = extractPressureHO2(tc2)
        else: 
            if order < 2:
                tc = extractPressureHO(tc)
            else:
                tc = extractPressureHO2(tc)

    #====================================
    # Extraction des grandeurs a la paroi
    #====================================
    if tc is None: 
        zw = Internal.getZones(tb)
        zw = T.join(zw)
    else:
        zw = TIBM.extractIBMWallFields(tc, tb=tb, famZones=famZones)

    #====================================
    # Extract pressure info from tc2 to tc
    #====================================
    if tc2 is not None:
        zw2 = TIBM.extractIBMWallFields(tc2, tb=tb, famZones=famZones, front=1)

        zones_zw = []
        zones_zw2 = []
        for zone in Internal.getZones(zw): zones_zw.append(zone[0])
        for zone in Internal.getZones(zw2): zones_zw2.append(zone[0])
        nbZones = len(zones_zw)

        for i in range(nbZones): # for multi corps
            szw = Internal.getNodeFromName(zw, zones_zw[i])
            szw2 = Internal.getNodeFromName(zw2, zones_zw2[i])

            Pressure2 = Internal.getNodeFromName(szw2, 'Pressure')[1]
            Internal.getNodeFromName(szw, 'Pressure')[1] = Pressure2

            Density2 = Internal.getNodeFromName(szw2, 'Density')[1]
            Internal.getNodeFromName(szw, 'Density')[1] = Density2

            gradxPressure2 = Internal.getNodeFromName(szw2, 'gradxPressure')[1]
            Internal.getNodeFromName(szw, 'gradxPressure')[1] = gradxPressure2
            gradyPressure2 = Internal.getNodeFromName(szw2, 'gradyPressure')[1]
            Internal.getNodeFromName(szw, 'gradyPressure')[1] = gradyPressure2
            gradzPressure2 = Internal.getNodeFromName(szw2, 'gradzPressure')[1]
            Internal.getNodeFromName(szw, 'gradzPressure')[1] = gradzPressure2

    dimPb = Internal.getValue(Internal.getNodeFromName(tb, 'EquationDimension'))

    if dimPb == 2: T._addkplane(zw)

    zw = C.convertArray2Tetra(zw)
    zw = T.reorderAll(zw, 1)

    ts = C.newPyTree(['SKIN']); ts[2][1][2]=zw[2][1][2]

    #==============================
    # Reference state
    #==============================
    RefState = Internal.getNodeFromType(tb,'ReferenceState_t')
    ts[2][1][2].append(RefState)
    ts = loads0(ts, Sref=Sref, alpha=alpha, beta=beta, dimPb=dimPb, verbose=True)

    if dimPb == 2: # reextrait en 2D
        ts = P.isoSurfMC(ts, "CoordinateZ", 0.)
        nodes = Internal.getNodesFromName(ts, 'CoordinateX')
        xmin = numpy.min(nodes[0][1])
        xmax = numpy.max(nodes[0][1])
        dxi = 1./(xmax-xmin)
        C._initVars(ts, 'xsc=({CoordinateX}-%g)*%g'%(xmin, dxi))

    if isinstance(wall_out, str): C.convertPyTree2File(ts, wall_out)
    return ts

#==========================================================================================
# In: ts: skin (TRI zones) distributed already (partial tree here)
# tc : transfer tree
#out: tl, graphWPOST : NODE-type zones of IBM points to be projected locally on ts
#out: graphWPOST: graph of coms between tc and tl 
#==========================================================================================
def _prepareSkinReconstruction(ts, tc):
    tBBs=Cmpi.createBBoxTree(ts)
    procDictBBs = Cmpi.getProcDict(tBBs)

    basename=Internal.getName(Internal.getBases(ts)[0])
    tl = C.newPyTree([basename])
    utauPresent = 0; vxPresent=0; yplusPresent = 0
    hmin = 0.
    for zc in Internal.getZones(tc):
        allIBCD = Internal.getNodesFromType(zc,"ZoneSubRegion_t")
        allIBCD = Internal.getNodesFromName(allIBCD,"IBCD_*")                  
        GCnode = Internal.getNodeFromType(zc,"GridCoordinates_t")
        XN = Internal.getNodeFromName(GCnode,'CoordinateX')

        for IBCD in allIBCD:
            if XN is not None:
                if XN[1].shape[0]>1:
                    hx = C.getValue(zc,'CoordinateX',1)-C.getValue(zc,'CoordinateX',0)
                    hy = C.getValue(zc,'CoordinateY',1)-C.getValue(zc,'CoordinateY',0)
                    hloc = max(abs(hx),abs(hy))
                    hmin = max(hloc,hmin)

            zname = Internal.getValue(IBCD)
            XW = Internal.getNodeFromName(IBCD,'CoordinateX_PW')[1]
            YW = Internal.getNodeFromName(IBCD,'CoordinateY_PW')[1]
            ZW = Internal.getNodeFromName(IBCD,'CoordinateZ_PW')[1]

            zsize = numpy.empty((1,3), numpy.int32, order='F')
            zsize[0,0] = XW.shape[0]; zsize[0,1] = 0; zsize[0,2] = 0
            z = Internal.newZone(name='IBW_Wall_%s_%s'%(zc[0],zname),zsize=zsize,
                                 ztype='Unstructured')
            gc = Internal.newGridCoordinates(parent=z)
            coordx = ['CoordinateX',XW,[],'DataArray_t']
            coordy = ['CoordinateY',YW,[],'DataArray_t']
            coordz = ['CoordinateZ',ZW,[],'DataArray_t']
            gc[2] = [coordx,coordy,coordz]
            n = Internal.createChild(z, 'GridElements', 'Elements_t', [2,0])
            Internal.createChild(n, 'ElementRange', 'IndexRange_t', [1,0])
            Internal.createChild(n, 'ElementConnectivity', 'DataArray_t', None)
            FSN = Internal.newFlowSolution(name=Internal.__FlowSolutionNodes__,
                                           gridLocation='Vertex', parent=z)
            pressNP = []; utauNP = []; yplusNP = []; densNP = []
            vxNP = []; vyNP = []; vzNP = []

            PW = Internal.getNodeFromName1(IBCD,XOD.__PRESSURE__)
            if PW is not None: pressNP.append(PW[1])
            RHOW = Internal.getNodeFromName1(IBCD,XOD.__DENSITY__)
            if RHOW is not None: densNP.append(RHOW[1])
            UTAUW = Internal.getNodeFromName1(IBCD,XOD.__UTAU__)
            if UTAUW is not None: utauNP.append(UTAUW[1])
            YPLUSW = Internal.getNodeFromName1(IBCD, XOD.__YPLUS__)
            if YPLUSW is not None: yplusNP.append(YPLUSW[1])

            VXW = Internal.getNodeFromName1(IBCD, XOD.__VELOCITYX__)
            if VXW is not None: vxNP.append(VXW[1])
            VYW = Internal.getNodeFromName1(IBCD, XOD.__VELOCITYY__)
            if VYW is not None: vyNP.append(VYW[1])
            VZW = Internal.getNodeFromName1(IBCD, XOD.__VELOCITYZ__)
            if VZW is not None: vzNP.append(VZW[1])

            FSN[2].append([XOD.__PRESSURE__,pressNP[0], [],'DataArray_t'])
            FSN[2].append([XOD.__DENSITY__,densNP[0], [],'DataArray_t'])
            if utauNP != []:
                utauPresent = 1
                FSN[2].append([XOD.__UTAU__,utauNP[0], [],'DataArray_t'])
            if yplusNP != []:
                yplusPresent = 1
                FSN[2].append([XOD.__YPLUS__,yplusNP[0], [],'DataArray_t'])

            if vxNP != []:
                vxPresent = 1
                FSN[2].append([XOD.__VELOCITYX__,vxNP[0], [],'DataArray_t'])
                FSN[2].append([XOD.__VELOCITYY__,vyNP[0], [],'DataArray_t'])
                FSN[2].append([XOD.__VELOCITYZ__,vzNP[0], [],'DataArray_t'])

            Cmpi._setProc(z,Cmpi.rank)          
            tl[2][1][2].append(z)

    tlBB=Cmpi.createBBoxTree(tl, tol=hmin)
    procDictWPOST = Cmpi.getProcDict(tlBB)
    interDictWPOST = X.getIntersectingDomains(tlBB, tBBs)
    graphWPOST = Cmpi.computeGraph(tlBB, type='bbox3',intersectionsDict=interDictWPOST,
                                   procDict=procDictWPOST, procDict2=procDictBBs, t2=tBBs)
 
    RefStateNode = Internal.getNodeFromName(ts,'ReferenceState')
    tl[2][1][2].append(RefStateNode)
    FES =  Internal.getNodeFromName(ts,'FlowEquationSet')
    tl[2][1][2].append(FES)

    C._initVars(ts,XOD.__PRESSURE__,0.)
    C._initVars(ts,XOD.__DENSITY__,0.)
    C._initVars(ts,XOD.__VELOCITYX__,0.)
    C._initVars(ts,XOD.__VELOCITYY__,0.)
    C._initVars(ts,XOD.__VELOCITYZ__,0.)    
    if Internal.getValue(Internal.getNodeFromType1(FES,'GoverningEquations_t'))!= 'Euler':
        C._initVars(ts,XOD.__UTAU__,0.)
        C._initVars(ts,XOD.__YPLUS__,0.)
    
    return tl, graphWPOST, interDictWPOST

# Distributed skin reconstruction (unsteady)
def _computeSkinVariables(ts, tc, tl, graphWPOST, interDictWPOST):
    for zc in Internal.getZones(tc):
        allIBCD = Internal.getNodesFromType(zc,"ZoneSubRegion_t")
        allIBCD = Internal.getNodesFromName(allIBCD,"IBCD_*")
        for IBCD in allIBCD:
            PW = Internal.getNodeFromName1(IBCD,XOD.__PRESSURE__)
            RHOW = Internal.getNodeFromName1(IBCD,XOD.__DENSITY__)
            UTAUW = Internal.getNodeFromName1(IBCD,XOD.__UTAU__)
            YPLUSW = Internal.getNodeFromName1(IBCD, XOD.__YPLUS__)
            VXW = Internal.getNodeFromName1(IBCD, XOD.__VELOCITYX__)
            VYW = Internal.getNodeFromName1(IBCD, XOD.__VELOCITYY__)
            VZW = Internal.getNodeFromName1(IBCD, XOD.__VELOCITYZ__)
            
            zname = Internal.getValue(IBCD)
            znamepostw = 'IBW_Wall_%s_%s'%(zc[0],zname)
            zpostw = Internal.getNodeFromName(tl,znamepostw)
            FSP = Internal.getNodeFromType(zpostw,'FlowSolution_t')
            PW2 = Internal.getNodeFromName1(FSP,XOD.__PRESSURE__)
            RHOW2 = Internal.getNodeFromName1(FSP,XOD.__DENSITY__)
            PW2[1]=PW[1]; RHOW2[1]=RHOW[1]

            UTAUW2 = Internal.getNodeFromName1(FSP,XOD.__UTAU__)
            if UTAUW2 is not None:
                YPLUSW2 = Internal.getNodeFromName1(FSP, XOD.__YPLUS__)
                UTAUW2[1]=UTAUW[1]; YPLUSW2[1]=YPLUSW[1]
            VXW2 = Internal.getNodeFromName1(FSP, XOD.__VELOCITYX__)     
            if VXW2 is not None:
                VYW2 = Internal.getNodeFromName1(FSP, XOD.__VELOCITYY__)
                VZW2 = Internal.getNodeFromName1(FSP, XOD.__VELOCITYZ__)
                VXW2[1]=VXW[1]
                VYW2[1]=VYW[1]
                VZW2[1]=VZW[1]


    tdl = Cmpi.addXZones(tl, graphWPOST)
    tdl = Cmpi.convert2PartialTree(tdl)
    for nobs in range(len(ts[2])):
        if Internal.getType(ts[2][nobs])=='CGNSBase_t':
            for nozs in range(len(ts[2][nobs][2])):
                zs = ts[2][nobs][2][nozs]
                if Internal.getType(zs)=='Zone_t':
                    cloud = []
                    for zl in Internal.getZones(tdl):
                        if zl != [] and zl is not None and zs[0] in interDictWPOST[zl[0]]:
                            zl = C.convertArray2Node(zl)
                            cloud.append(zl)
    
                    if cloud != []:
                        cloud = T.join(cloud)
                        
                        ts[2][nobs][2][nozs] = P.projectCloudSolution(cloud, zs, dim=3)
                        
    return None

#=============================================================================
# Post efforts
# IN: t_case: fichier ou arbre du cas
# IN: tc_in: fichier ou arbre de connectivite contenant les IBCD
# si tc_in =None, t_case est la surface avec la solution deja projetee
# OUT: wall_out ou None: fichier pour sortie des efforts sur la paroi aux centres
# IN: alpha: angle pour les efforts
# IN: beta: angle pour les efforts
#==============================================================================
def _unsteadyLoads(tb, Sref=None, alpha=0., beta=0.):
    zones = KCOMM.allgather(Internal.getZones(tb))
    ts = Distributed.setZonesInTree(tb, zones)
    dimPb = Internal.getValue(Internal.getNodeFromName(ts, 'EquationDimension'))
    return _loads0(ts, Sref=Sref, alpha=alpha, beta=beta, dimPb=dimPb, verbose =False)

#=============================================================================
# Correction de la pression en post-traitement
#=============================================================================
def extractPressureHO(tc):

    for z in Internal.getZones(tc):
        subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        for zsr in subRegions:
            nameSubRegion = zsr[0]
            if (nameSubRegion[:4] == 'IBCD' or nameSubRegion[:4] == '2_IB'):
                gradxPressure = Internal.getNodeFromName(zsr, 'gradxPressure')[1]
                gradyPressure = Internal.getNodeFromName(zsr, 'gradyPressure')[1]
                gradzPressure = Internal.getNodeFromName(zsr, 'gradzPressure')[1]

                CoordinateX = Internal.getNodeFromName(zsr, 'CoordinateX_PW')[1]
                CoordinateY = Internal.getNodeFromName(zsr, 'CoordinateY_PW')[1]
                CoordinateZ = Internal.getNodeFromName(zsr, 'CoordinateZ_PW')[1]

                CoordinateX_PC = Internal.getNodeFromName(zsr, 'CoordinateX_PC')[1]
                CoordinateY_PC = Internal.getNodeFromName(zsr, 'CoordinateY_PC')[1]
                CoordinateZ_PC = Internal.getNodeFromName(zsr, 'CoordinateZ_PC')[1]

                CoordinateX_PI = Internal.getNodeFromName(zsr, 'CoordinateX_PI')[1]
                CoordinateY_PI = Internal.getNodeFromName(zsr, 'CoordinateY_PI')[1]
                CoordinateZ_PI = Internal.getNodeFromName(zsr, 'CoordinateZ_PI')[1]

                Pressure = Internal.getNodeFromName(zsr, 'Pressure')[1]
                Density = Internal.getNodeFromName(zsr, 'Density')[1]

                nIBC = numpy.shape(CoordinateX)[0]
                for i in range(nIBC):
                    nx = CoordinateX_PC[i] - CoordinateX[i]
                    ny = CoordinateY_PC[i] - CoordinateY[i]
                    nz = CoordinateZ_PC[i] - CoordinateZ[i]
                    norm = math.sqrt(nx*nx + ny*ny + nz*nz)
                    nx = nx/norm
                    ny = ny/norm
                    nz = nz/norm

                    nGradP = nx*gradxPressure[i] + ny*gradyPressure[i] + nz*gradzPressure[i]

                    bx = CoordinateX_PI[i] - CoordinateX[i]
                    by = CoordinateY_PI[i] - CoordinateY[i]
                    bz = CoordinateZ_PI[i] - CoordinateZ[i]
                    beta = math.sqrt(bx*bx + by*by + bz*bz)

                    Density[i] = Density[i]/Pressure[i]*(Pressure[i] - nGradP*beta)
                    Pressure[i] = Pressure[i] - nGradP*beta

                Internal.getNodeFromName(zsr, 'Pressure')[1] = Pressure
                Internal.getNodeFromName(zsr, 'Density')[1]  = Density

    return tc

#=============================================================================
# Correction de la pression en post-traitement
#=============================================================================    
def extractPressureHO2(tc):

    for z in Internal.getZones(tc):
        subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        for zsr in subRegions:
            nameSubRegion = zsr[0]
            if (nameSubRegion[:4] == 'IBCD' or nameSubRegion[:4] == '2_IB'):

                CoordinateX = Internal.getNodeFromName(zsr, 'CoordinateX_PW')[1]
                CoordinateY = Internal.getNodeFromName(zsr, 'CoordinateY_PW')[1]
                CoordinateZ = Internal.getNodeFromName(zsr, 'CoordinateZ_PW')[1]

                CoordinateX_PC = Internal.getNodeFromName(zsr, 'CoordinateX_PC')[1]
                CoordinateY_PC = Internal.getNodeFromName(zsr, 'CoordinateY_PC')[1]
                CoordinateZ_PC = Internal.getNodeFromName(zsr, 'CoordinateZ_PC')[1]

                CoordinateX_PI = Internal.getNodeFromName(zsr, 'CoordinateX_PI')[1]
                CoordinateY_PI = Internal.getNodeFromName(zsr, 'CoordinateY_PI')[1]
                CoordinateZ_PI = Internal.getNodeFromName(zsr, 'CoordinateZ_PI')[1]

                gradxPressure  = Internal.getNodeFromName(zsr, 'gradxPressure')[1]
                gradyPressure  = Internal.getNodeFromName(zsr, 'gradyPressure')[1]
                gradzPressure  = Internal.getNodeFromName(zsr, 'gradzPressure')[1]

                gradxxPressure = Internal.getNodeFromName(zsr, 'gradxxPressure')[1]
                gradxyPressure = Internal.getNodeFromName(zsr, 'gradxyPressure')[1]
                gradxzPressure = Internal.getNodeFromName(zsr, 'gradxzPressure')[1]

                gradyxPressure = Internal.getNodeFromName(zsr, 'gradyxPressure')[1]
                gradyyPressure = Internal.getNodeFromName(zsr, 'gradyyPressure')[1]
                gradyzPressure = Internal.getNodeFromName(zsr, 'gradyzPressure')[1]

                gradzxPressure = Internal.getNodeFromName(zsr, 'gradzxPressure')[1]
                gradzyPressure = Internal.getNodeFromName(zsr, 'gradzyPressure')[1]
                gradzzPressure = Internal.getNodeFromName(zsr, 'gradzzPressure')[1]

                Pressure = Internal.getNodeFromName(zsr, 'Pressure')[1]
                Density = Internal.getNodeFromName(zsr, 'Density')[1]

                nIBC = numpy.shape(CoordinateX)[0]
                for i in range(nIBC):
                    nx = CoordinateX_PC[i] - CoordinateX[i]
                    ny = CoordinateY_PC[i] - CoordinateY[i]
                    nz = CoordinateZ_PC[i] - CoordinateZ[i]
                    norm = math.sqrt(nx*nx + ny*ny + nz*nz)
                    nx = nx/norm
                    ny = ny/norm
                    nz = nz/norm

                    nGradP   = nx*gradxPressure[i] + ny*gradyPressure[i] + nz*gradzPressure[i]

                    nnxGradP = nx*gradxxPressure[i] + ny*gradxyPressure[i] + nz*gradxzPressure[i]
                    nnyGradP = nx*gradyxPressure[i] + ny*gradyyPressure[i] + nz*gradyzPressure[i]
                    nnzGradP = nx*gradzxPressure[i] + ny*gradzyPressure[i] + nz*gradzzPressure[i]
                    nnGradP  = nx*nnxGradP + ny*nnyGradP + nz*nnzGradP

                    bx = CoordinateX_PI[i] - CoordinateX[i]
                    by = CoordinateY_PI[i] - CoordinateY[i]
                    bz = CoordinateZ_PI[i] - CoordinateZ[i]
                    beta = math.sqrt(bx*bx + by*by + bz*bz)

                    Density[i] = Density[i]/Pressure[i]*(Pressure[i] - nGradP*beta + 0.5*nnGradP*beta**2)
                    Pressure[i] = Pressure[i] - nGradP*beta + 0.5*nnGradP*beta**2

                Internal.getNodeFromName(zsr, 'Pressure')[1] = Pressure
                Internal.getNodeFromName(zsr, 'Density')[1]  = Density

    return tc

#=============================================================================
# Calcul des termes convectifs (TBLE FULL)
#=============================================================================
def extractConvectiveTerms(tc):

    for z in Internal.getZones(tc):
        subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        for zsr in subRegions:
            nameSubRegion = zsr[0]
            if (nameSubRegion[:4] == 'IBCD' or nameSubRegion[:4] == '2_IB'):

                CoordinateX    = Internal.getNodeFromName(zsr, 'CoordinateX_PW')[1]
                CoordinateY    = Internal.getNodeFromName(zsr, 'CoordinateY_PW')[1]
                CoordinateZ    = Internal.getNodeFromName(zsr, 'CoordinateZ_PW')[1]

                CoordinateX_PC = Internal.getNodeFromName(zsr, 'CoordinateX_PC')[1]
                CoordinateY_PC = Internal.getNodeFromName(zsr, 'CoordinateY_PC')[1]
                CoordinateZ_PC = Internal.getNodeFromName(zsr, 'CoordinateZ_PC')[1]

                CoordinateX_PI = Internal.getNodeFromName(zsr, 'CoordinateX_PI')[1]
                CoordinateY_PI = Internal.getNodeFromName(zsr, 'CoordinateY_PI')[1]
                CoordinateZ_PI = Internal.getNodeFromName(zsr, 'CoordinateZ_PI')[1]

                Pressure       = Internal.getNodeFromName(zsr, 'Pressure')[1]
                Density        = Internal.getNodeFromName(zsr, 'Density')[1]

                gradxPressure  = Internal.getNodeFromName(zsr, 'gradxPressure')[1]
                gradyPressure  = Internal.getNodeFromName(zsr, 'gradyPressure')[1]
                gradzPressure  = Internal.getNodeFromName(zsr, 'gradzPressure')[1]

                gradxVelocityX = Internal.getNodeFromName(zsr, 'gradxVelocityX')[1]
                gradyVelocityX = Internal.getNodeFromName(zsr, 'gradyVelocityX')[1]
                gradzVelocityX = Internal.getNodeFromName(zsr, 'gradzVelocityX')[1]

                gradxVelocityY = Internal.getNodeFromName(zsr, 'gradxVelocityY')[1]
                gradyVelocityY = Internal.getNodeFromName(zsr, 'gradyVelocityY')[1]
                gradzVelocityY = Internal.getNodeFromName(zsr, 'gradzVelocityY')[1]

                gradxVelocityZ = Internal.getNodeFromName(zsr, 'gradxVelocityZ')[1]
                gradyVelocityZ = Internal.getNodeFromName(zsr, 'gradyVelocityZ')[1]
                gradzVelocityZ = Internal.getNodeFromName(zsr, 'gradzVelocityZ')[1]

                VelocityX      = Internal.getNodeFromName(zsr, 'VelocityX')[1]
                VelocityY      = Internal.getNodeFromName(zsr, 'VelocityY')[1]
                VelocityZ      = Internal.getNodeFromName(zsr, 'VelocityZ')[1]

                nIBC = numpy.shape(CoordinateX)[0]
                conv1  = numpy.zeros((nIBC),numpy.float64)
                conv2  = numpy.zeros((nIBC),numpy.float64)

                for i in range(nIBC):
                    nx = CoordinateX_PC[i] - CoordinateX[i]
                    ny = CoordinateY_PC[i] - CoordinateY[i]
                    nz = CoordinateZ_PC[i] - CoordinateZ[i]
                    norm_n = math.sqrt(nx*nx + ny*ny + nz*nz)
                    nx = nx/norm_n
                    ny = ny/norm_n
                    nz = nz/norm_n

                    uscaln = (VelocityX[i]*nx + VelocityY[i]*ny + VelocityZ[i]*nz)

                    tx = VelocityX[i] - uscaln*nx
                    ty = VelocityY[i] - uscaln*ny
                    tz = VelocityZ[i] - uscaln*nz
                    norm_t = math.sqrt(tx*tx + ty*ty + tz*tz)
                    tx = tx/norm_t
                    ty = ty/norm_t
                    tz = tz/norm_t

                    tgradU =  (tx*gradxVelocityX[i] + ty*gradyVelocityX[i] + tz*gradzVelocityX[i])*tx
                    tgradU += (tx*gradxVelocityY[i] + ty*gradyVelocityY[i] + tz*gradzVelocityY[i])*ty
                    tgradU += (tx*gradxVelocityZ[i] + ty*gradyVelocityZ[i] + tz*gradzVelocityZ[i])*tz

                    ngradU =  (nx*gradxVelocityX[i] + ny*gradyVelocityX[i] + nz*gradzVelocityX[i])*tx
                    ngradU += (nx*gradxVelocityY[i] + ny*gradyVelocityY[i] + nz*gradzVelocityY[i])*ty
                    ngradU += (nx*gradxVelocityZ[i] + ny*gradyVelocityZ[i] + nz*gradzVelocityZ[i])*tz

                    norm_n = math.sqrt((uscaln*nx)**2 + (uscaln*ny)**2 + (uscaln*nz)**2)

                    conv1[i] = norm_t*tgradU
                    conv2[i] = norm_n*ngradU

                zsr[2].append(['conv1'  , conv1  , [], 'DataArray_t'])
                zsr[2].append(['conv2'  , conv2  , [], 'DataArray_t'])

    return tc

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

    # Affichage du nombre de points par proc - equilibrage ou pas
    NptsTot = 0
    for i in range(NP):
        NPTS = 0
        for z in Internal.getZones(ts):
            if Cmpi.getProc(z) == i: NPTS += C.getNPts(z)
        NptsTot += NPTS
        print('Rank {} has {} points'.format(i,NPTS))
    print('All points: {} million points'.format(NptsTot/1.e6))
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

    tcw = TIBM.createIBMWZones(tc,variables=[])

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


#====================================================================================
# Prend les snears dans t, les multiplie par factor
def snearFactor(t, factor=1.):
    tp = Internal.copyRef(t)
    _snearFactor(t, factor)
    return tp

def _snearFactor(t, factor=1.):
    zones = Internal.getZones(t)
    for z in zones:
        nodes = Internal.getNodesFromName2(z, 'snear')
        for n in nodes:
            Internal._setValue(n, factor*Internal.getValue(n))
    return None

# Set IBC type in zones
def setIBCType(t, value):
    tp = Internal.copyRef(t)
    _setIBCType(t, value)
    return tp

def _setIBCType(z, value):
    zones = Internal.getZones(z)
    for z in zones:
        Internal._createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
        n = Internal.getNodeFromName1(z, '.Solver#define')
        Internal._createUniqueChild(n, 'ibctype', 'DataArray_t', value)
    return None

# Set snear in zones
def setSnear(t, value):
    tp = Internal.copyRef(t)
    _setSnear(t, value)
    return tp

def _setSnear(z, value):
    zones = Internal.getZones(z)
    for z in zones:
        Internal._createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
        n = Internal.getNodeFromName1(z, '.Solver#define')
        Internal._createUniqueChild(n, 'snear', 'DataArray_t', value)
    return None

# Set dfar in zones
def setDfar(t, value):
    tp = Internal.copyRef(t)
    _setDfar(t, value)
    return tp

def _setDfar(z, value):
    zones = Internal.getZones(z)
    for z in zones:
        Internal._createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
        n = Internal.getNodeFromName1(z, '.Solver#define')
        Internal._createUniqueChild(n, 'dfar', 'DataArray_t', value)
    return None

def _modifIBCD(tc):
    raise NotImplementedError("_modifyIBCD is obsolete. Use _initOutflow and _initInj functions.")

# set Pressure to P_tot for a IBC of type outpress of family name FamilyName
def _initOutflow(tc, familyNameOutflow, P_tot):
    for zc in Internal.getZones(tc):
        for zsr in Internal.getNodesFromName(zc,'IBCD_4_*'):
            FamNode = Internal.getNodeFromType1(zsr,'FamilyName_t')
            if FamNode is not None:
                FamName = Internal.getValue(FamNode)
                if FamName==familyNameOutflow:
                    stagPNode =  Internal.getNodeFromName(zsr,'Pressure')    
                    sizeIBC = numpy.shape(stagPNode[1])
                    Internal.setValue(stagPNode,P_tot*numpy.ones(sizeIBC))
    return None

def _initInj(tc, familyNameInj, P_tot, H_tot, injDir=[1.,0.,0.]):
    for zc in Internal.getZones(tc):
        for zsr in Internal.getNodesFromName(zc,'IBCD_5_*'):
            FamNode = Internal.getNodeFromType1(zsr,'FamilyName_t')
            if FamNode is not None:
                FamName = Internal.getValue(FamNode)
                if FamName==familyNameInj:
                    stagPNode =  Internal.getNodeFromName(zsr,'StagnationPressure')
                    stagHNode =  Internal.getNodeFromName(zsr,'StagnationEnthalpy')
                    dirxNode = Internal.getNodeFromName(zsr,'dirx')
                    diryNode = Internal.getNodeFromName(zsr,'diry')
                    dirzNode = Internal.getNodeFromName(zsr,'dirz')
                    sizeIBC = numpy.shape(stagHNode[1])
                    Internal.setValue(stagHNode,H_tot*numpy.ones(sizeIBC))
                    Internal.setValue(stagPNode,P_tot*numpy.ones(sizeIBC))
                    if injDir[0] != 0.: 
                        Internal.setValue(dirxNode, injDir[0]*numpy.ones(sizeIBC))
                    else:
                        Internal.setValue(dirxNode, numpy.zeros(sizeIBC)) 

                    if injDir[1] != 0.: 
                        Internal.setValue(diryNode, injDir[1]*numpy.ones(sizeIBC))
                    else:
                        Internal.setValue(diryNode, numpy.zeros(sizeIBC))

                    if injDir[2] != 0.: 
                        Internal.setValue(dirzNode, injDir[2]*numpy.ones(sizeIBC))
                    else:
                        Internal.setValue(dirzNode, numpy.zeros(sizeIBC)) 
                    
    return None

def changeBCType(tc, oldBCType, newBCType):
    for z in Internal.getZones(tc):
        subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        for zsr in subRegions:
            nameSubRegion = zsr[0]
            if nameSubRegion[:4] == "IBCD":
                bcType = int(nameSubRegion.split("_")[1])
                if bcType == oldBCType:
                    zsr[0] = "IBCD_{}_".format(newBCType)+"_".join(nameSubRegion.split("_")[2:])

                    pressure = Internal.getNodeFromName(zsr, 'Pressure')[1]
                    nIBC = pressure.shape[0]

                    Internal._rmNodesByName(zsr, 'utau')
                    Internal._rmNodesByName(zsr, 'yplus')

                    Internal._rmNodesByName(zsr, 'StagnationEnthalpy')
                    Internal._rmNodesByName(zsr, 'StagnationPressure')
                    Internal._rmNodesByName(zsr, 'dirx')
                    Internal._rmNodesByName(zsr, 'diry')
                    Internal._rmNodesByName(zsr, 'dirz')

                    Internal._rmNodesByName(zsr, 'gradxPressure')
                    Internal._rmNodesByName(zsr, 'gradyPressure')
                    Internal._rmNodesByName(zsr, 'gradzPressure')

                    Internal._rmNodesByName(zsr, 'gradxVelocityX')
                    Internal._rmNodesByName(zsr, 'gradyVelocityX')
                    Internal._rmNodesByName(zsr, 'gradzVelocityX')

                    Internal._rmNodesByName(zsr, 'gradxVelocityY')
                    Internal._rmNodesByName(zsr, 'gradyVelocityY')
                    Internal._rmNodesByName(zsr, 'gradzVelocityY')

                    Internal._rmNodesByName(zsr, 'gradxVelocityZ')
                    Internal._rmNodesByName(zsr, 'gradyVelocityZ')
                    Internal._rmNodesByName(zsr, 'gradzVelocityZ')

                    Internal._rmNodesByName(zsr, 'KCurv')

                    if newBCType in [2, 3, 6, 10, 11]:
                        utauNP  = numpy.zeros((nIBC),numpy.float64)
                        yplusNP = numpy.zeros((nIBC),numpy.float64)
                        zsr[2].append(['utau' , utauNP , [], 'DataArray_t'])
                        zsr[2].append(['yplus', yplusNP, [], 'DataArray_t'])

                    if newBCType == 5:
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

                    if newBCType == 100:
                        KCurvNP = numpy.zeros((nIBC),numpy.float64)
                        zsr[2].append(["KCurv" , KCurvNP , [], 'DataArray_t'])

                    if newBCType == 10 or newBCType == 11:
                        gradxPressureNP  = numpy.zeros((nIBC),numpy.float64)
                        zsr[2].append(['gradxPressure' , gradxPressureNP , [], 'DataArray_t'])
                        gradyPressureNP  = numpy.zeros((nIBC),numpy.float64)
                        zsr[2].append(['gradyPressure' , gradyPressureNP , [], 'DataArray_t'])
                        gradzPressureNP  = numpy.zeros((nIBC),numpy.float64)
                        zsr[2].append(['gradzPressure' , gradzPressureNP , [], 'DataArray_t'])

                    if newBCType == 11:
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

    return tc

def transformTc2(tc2):
    for z in Internal.getZones(tc2):
        subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        for zsr in subRegions:
            nameSubRegion = zsr[0]
            if nameSubRegion[:6] == "2_IBCD":
                ibctype = int(nameSubRegion.split("_")[2])
                zsr[0] = "IBCD_{}_".format(ibctype)+"_".join(nameSubRegion.split("_")[3:])

                pressure = Internal.getNodeFromName(zsr, 'Pressure')[1]
                nIBC = pressure.shape[0]

                Internal._rmNodesByName(zsr, 'Density')

                Internal._rmNodesByName(zsr, 'VelocityX')
                Internal._rmNodesByName(zsr, 'VelocityY')
                Internal._rmNodesByName(zsr, 'VelocityZ')

                Internal._rmNodesByName(zsr, 'utau')
                Internal._rmNodesByName(zsr, 'yplus')

                Internal._rmNodesByName(zsr, 'StagnationEnthalpy')
                Internal._rmNodesByName(zsr, 'StagnationPressure')
                Internal._rmNodesByName(zsr, 'dirx')
                Internal._rmNodesByName(zsr, 'diry')
                Internal._rmNodesByName(zsr, 'dirz')

                Internal._rmNodesByName(zsr, 'gradxPressure')
                Internal._rmNodesByName(zsr, 'gradyPressure')
                Internal._rmNodesByName(zsr, 'gradzPressure')

                Internal._rmNodesByName(zsr, 'gradxVelocityX')
                Internal._rmNodesByName(zsr, 'gradyVelocityX')
                Internal._rmNodesByName(zsr, 'gradzVelocityX')

                Internal._rmNodesByName(zsr, 'gradxVelocityY')
                Internal._rmNodesByName(zsr, 'gradyVelocityY')
                Internal._rmNodesByName(zsr, 'gradzVelocityY')

                Internal._rmNodesByName(zsr, 'gradxVelocityZ')
                Internal._rmNodesByName(zsr, 'gradyVelocityZ')
                Internal._rmNodesByName(zsr, 'gradzVelocityZ')

                Internal._rmNodesByName(zsr, 'KCurv')

                DensityNP = numpy.zeros((nIBC),numpy.float64)
                zsr[2].append(['Density' , DensityNP , [], 'DataArray_t'])

                VelocityXNP = numpy.zeros((nIBC),numpy.float64)
                zsr[2].append(['VeloicityX' , VelocityXNP , [], 'DataArray_t'])
                VelocityYNP= numpy.zeros((nIBC),numpy.float64)
                zsr[2].append(['VeloicityY' , VelocityYNP , [], 'DataArray_t'])
                VelocityZNP = numpy.zeros((nIBC),numpy.float64)
                zsr[2].append(['VeloicityZ' , VelocityZNP , [], 'DataArray_t'])

                if ibctype in [2, 3, 6, 10, 11]:
                    utauNP  = numpy.zeros((nIBC),numpy.float64)
                    yplusNP = numpy.zeros((nIBC),numpy.float64)
                    zsr[2].append(['utau' , utauNP , [], 'DataArray_t'])
                    zsr[2].append(['yplus', yplusNP, [], 'DataArray_t'])

                if ibctype == 5:
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

                if ibctype == 100:
                    KCurvNP = numpy.zeros((nIBC),numpy.float64)
                    zsr[2].append(["KCurv" , KCurvNP , [], 'DataArray_t'])

                if ibctype == 10 or ibctype == 11:
                    gradxPressureNP  = numpy.zeros((nIBC),numpy.float64)
                    zsr[2].append(['gradxPressure' , gradxPressureNP , [], 'DataArray_t'])
                    gradyPressureNP  = numpy.zeros((nIBC),numpy.float64)
                    zsr[2].append(['gradyPressure' , gradyPressureNP , [], 'DataArray_t'])
                    gradzPressureNP  = numpy.zeros((nIBC),numpy.float64)
                    zsr[2].append(['gradzPressure' , gradzPressureNP , [], 'DataArray_t'])

                if ibctype == 11:
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

    return tc2

#====================================================================================
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
