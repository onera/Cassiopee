"""Connector for IBM preprocessing"""
from . import PyTree as X
from . import OversetData as XOD
from . import Connector
from . import connector

import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Connector.Mpi as Xmpi
import Post.Mpi as Pmpi
import Generator.PyTree as G
import Initiator.PyTree as I
import Converter.Distributed as Distributed
import Generator.IBMmodelHeight as G_IBM_Height
import Geom.IBM as D_IBM
import Transform.PyTree as T
import Converter.Internal as Internal
import Compressor.PyTree as Compressor
import Dist2Walls.PyTree as DTW
import Distributor2.PyTree as D2
import Post.PyTree as P
import KCore.test as test
import Converter
import Generator
import Transform
import KCore
import numpy
import math

varsn       = ['gradxTurbulentDistance','gradyTurbulentDistance','gradzTurbulentDistance']
varsnDouble = ['gradxTurbulentDistanceDouble','gradyTurbulentDistanceDouble','gradzTurbulentDistanceDouble']
TOLDIST     = 1.e-14
SHIFTF      = 1.e-10
EPSCART     = 1.e-6
TOLCELLN    = 0.01

LBM_IBC_NUM = 113 #from param_solver.h

TypesOfIBC  = XOD.TypesOfIBC

# computes the friction velocity
def _computeFrictionVelocity(a):
    for z in Internal.getZones(a):
        connector._computeFrictionVelocity(z,
                                           Internal.__GridCoordinates__,
                                           Internal.__FlowSolutionNodes__,
                                           Internal.__FlowSolutionCenters__)
    return None

RENAMEIBCNODES=False

__IBCNameServer__={}

def _changeNameIBCD__(tc,NewIBCD):
    ZsubR = Internal.getNodesByType(tc, 'ZoneSubRegion_t')
    for z in ZsubR:
        zsplit = z[0].split('_')
        if zsplit[0] == 'IBCD':
            zsplit[1] = str(NewIBCD)
            znew = '_'.join(zsplit)
            Internal.setName(z, znew)
    return None

def getIBCDName(proposedName):
    global __IBCNameServer__
    (ibcname,__IBCNameServer__)=C.getUniqueName(proposedName, __IBCNameServer__)
    return ibcname


#====================================================================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## MACRO FUNCTIONS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def printTimeAndMemory__(message, time=-1, functionName='prepareIBMData'):
    Cmpi.barrier()
    test.printMem("Info: %s: %s [%s]"%(functionName,message, 'end' if time > 0 else 'start'))
    Cmpi.barrier()
    if time > 0 and Cmpi.rank == 0: print("Info: %s: %s running time = %4.4fs"%(functionName,message, time))
    Cmpi.barrier()

    return None

def computeMeshInfo__(z, dim):
    if dim[0] == 'Structured':
        # Compute nodes number
        np = dim[1]*dim[2]*dim[3]
        # Compute cell number
        if dim[4] == 1:
            ncells = dim[1]-1
        elif dim[4] == 2:
            ncells = (dim[1]-1)*(dim[2]-1)
        else:
            ncells = (dim[1]-1)*(dim[2]-1)*(dim[3]-1)
        # Compute face number
        if dim[4] == 1:
            nfaces = dim[1]
        elif dim[4] == 2:
            nfaces = dim[1]*(dim[2]-1)+dim[2]*(dim[1]-1)
        elif dim[4] == 3:
            nfaces = dim[1]*(dim[2]-1)*(dim[3]-1)+ \
                dim[2]*(dim[3]-1)*(dim[1]-1)+ \
                dim[3]*(dim[1]-1)*(dim[2]-1)
        return np, ncells, nfaces
    return 0, 0, 0

def _computeMeshInfo(t):
    np_total, nc_total, nf_total, nzones = 0,0,0,0
    nc0, nc1, nc2 = 0,0,0

    NP     = Cmpi.size
    rank   = Cmpi.rank
    NPTS   = numpy.zeros(NP, dtype=Internal.E_NpyInt)
    NCELLS = numpy.zeros(NP, dtype=Internal.E_NpyInt)
    for z in Internal.getZones(t):
        cellN = Internal.getNodeFromName(z, 'cellN')[1]
        c0 = numpy.count_nonzero(cellN == 0); nc0 += c0
        c1 = numpy.count_nonzero(cellN == 1); nc1 += c1
        c2 = numpy.count_nonzero(cellN == 2); nc2 += c2

        dim = Internal.getZoneDim(z)
        np, nc, nf = computeMeshInfo__(z, dim)
        np_total += np
        nc_total += nc
        nf_total += nf
        nzones   += 1

        zNoGhost      = C.rmGhostCells(z, z, 2)
        NPTS[rank]   += C.getNPts(zNoGhost)
        NCELLS[rank] += C.getNCells(zNoGhost)

    NPTS     = Cmpi.allreduce(NPTS  ,op=Cmpi.SUM)
    NCELLS   = Cmpi.allreduce(NCELLS,op=Cmpi.SUM)
    print("Info: mesh info for rank {}: number of points: {:.2f}M / number of cells: {:.2f}M".format(rank, np_total/1.e6, nc_total/1.e6))
    if rank == 0:
        NcellsTot     = numpy.sum(NCELLS)
        ncells_percent= []
        for i in range(NP):
            ncells_percent.append(NCELLS[i]/NcellsTot*100)
            print('Info: Rank {} has {:.3f}e06 points & {:.3f}e06 cells & {} % of cells - no ghost cells'.format(i,int(NPTS[i])/1e06,int(NCELLS[i])/1e06,round(ncells_percent[i],2)))
        print('Info: Range of % of cells: {} - {}'.format(round(min(ncells_percent),2),round(max(ncells_percent),2)))
    Cmpi.barrier()

    np_total = Cmpi.allreduce(np_total, op=Cmpi.SUM)
    nc_total = Cmpi.allreduce(nc_total, op=Cmpi.SUM)
    nc0 = Cmpi.allreduce(nc0, op=Cmpi.SUM)
    nc1 = Cmpi.allreduce(nc1, op=Cmpi.SUM)
    nc2 = Cmpi.allreduce(nc2, op=Cmpi.SUM)

    natures = [ncx/float(nc_total)*100. for ncx in [nc2, nc1, nc0]]

    if Cmpi.rank == 0:
        print("Info: global mesh info: Interpolated cells (cellN 2): {:.2f}%, Computed cells (cellN 1): {:.2f}%, Blanked cells (cellN 0): {:.2f}%".format(*natures))
        print("Info: global mesh info: total number of points: {:.2f}M / total number of cells: {:.2f}M".format(np_total/1.e6, nc_total/1.e6))
    Cmpi.barrier()

    return None

def prepareIBMData(t_case, t_out, tc_out, t_in=None, to=None, tbox=None, tinit=None, tbCurvi=None,
                   snears=0.01, snearsf=None, dfars=10., dfarDir=0, vmin=21, depth=2, frontType=1, octreeMode=0,
                   IBCType=1, verbose=True, expand=3, ext=-1, optimized=1,
                   check=False, twoFronts=False, cartesian=True, cleanCellN=True,
                   yplus=100., Lref=1., correctionMultiCorpsF42=False, blankingF42=False, wallAdaptF42=None, heightMaxF42=-1.):

    import Generator.IBM as G_IBM
    import time as python_time

    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = Internal.copyTree(t_case)

    ## Note: cartesian = True is left as an input argument to avoid regressing  during the non-regression test.
    ##       In the near future the ref. values for the non-regression tests will be updated with cartesian=True.
    ##       At this point, cartesian=True input argument can be deleted.
    ## Note: when cartesian = True is deleted as an input argument the line below must be uncommented and the cartesian in the if statement but be deleted.
    #cartesian = True
    if t_in and cartesian:
        cartesian = G_IBM.checkCartesian(t_in, nghost=2)
        if cartesian:
            RED  = "\033[1;31;40m"
            END  = "\033[0m"
            print("===========================================")
            print("Note: t_in is a " + RED + "CARTESIAN " + END + "grid")
            print("===========================================")
        else:
            RED  = "\033[1;31;40m"
            END  = "\033[0m"
            print("===========================================")
            print("Note: t_in is " + RED + "NOT" + END + " a " + RED + "CARTESIAN " + END + "grid")
            print("===========================================")

    ##[AJ] keep for now...will delete in the near future
    ##if isinstance(tc_out, str):
    ##    if '/' in tc_out: fileoutpre = tc_out.split('/')
    ##    else: fileoutpre = ['.','template.cgns']
    ##else: fileoutpre = ['.','template.cgns']

    refstate = Internal.getNodeFromName(tb, 'ReferenceState')
    flowEqn  = Internal.getNodeFromName(tb, 'FlowEquationSet')

    Reynolds = Internal.getNodeFromName(tb, 'Reynolds')
    if Reynolds is not None:
        Reynolds = Internal.getValue(Reynolds)
        if Reynolds < 1.e5: frontType = 1
    else: Reynolds = 1.e6

    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    if dimPb is None: raise ValueError('prepareIBMData: EquationDimension is missing in input geometry tree.')
    dimPb = Internal.getValue(dimPb)
    if dimPb == 2: C._initVars(tb, 'CoordinateZ', 0.)
    else: cleanCellN = True

    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('prepareIBMData: GoverningEquations is missing in input geometry tree.')
    model = Internal.getValue(model)

    ibctypes = Internal.getNodesFromName(tb, 'ibctype')
    if ibctypes is None: raise ValueError('prepareIBMData: ibc type is missing in input geometry tree.')
    ibctypes = list(set(Internal.getValue(ibc) for ibc in ibctypes))

    if model == 'Euler':
        if any(ibc in ['Musker', 'MuskerMob', 'Mafzal', 'Log', 'TBLE', 'TBLE_FULL'] for ibc in ibctypes):
            raise ValueError("prepareIBMData: governing equations (Euler) not consistent with ibc types %s"%(ibctypes))

    if frontType == 42 and tbox is None:
        print("Info: prepareIBMData: frontType 42 is used, but no tbox has been provided to ensure that the near-wall resolution is sufficiently propagated. Forcing expand=4.")
        expand = 4

    if ext == -1: ext = depth+1

    #===================
    # STEP 0 : GET FILAMENT BODIES
    #===================
    tb, tbFilament              = D_IBM.determineClosedSolidFilament__(tb)
    isFilamentOnly, isWireModel = D_IBM.localWMMFlags__(tb, tbFilament)

    #===================
    # STEP 1 : GENERATE MESH
    #===================
    if isFilamentOnly: tbLocal=tbFilament
    elif tbFilament: tbLocal=Internal.merge([tb,tbFilament])
    else: tbLocal=tb


    if t_in is None:
        if verbose: pt0 = python_time.time(); printTimeAndMemory__('generate Cartesian mesh', time=-1)
        test.printMem("Info: prepareIBMData: generate Cartesian mesh [start]")
        t = G_IBM.generateIBMMesh(tbLocal, vmin=vmin, snears=snears, dimPb=dimPb, dfars=dfars, tbox=tbox,
                                  snearsf=snearsf, check=check, to=to, ext=ext, optimized=optimized,
                                  expand=expand, dfarDir=dfarDir, octreeMode=octreeMode)
        Internal._rmNodesFromName(tb,"SYM")

        _redispatch__(t=t)
        if verbose: printTimeAndMemory__('generate Cartesian mesh', time=python_time.time()-pt0)

    else:
        t = t_in

    for b in Internal.getBases(t):
        Internal.addChild(b, refstate, pos=0)
        Internal.addChild(b, flowEqn , pos=0)

    #===================
    # STEP 2 : DIST2WALL
    #===================
    if verbose: pt0 = python_time.time(); printTimeAndMemory__('compute wall distance', time=-1)
    _dist2wallIBM(t, tb, dimPb=dimPb, frontType=frontType, Reynolds=Reynolds, yplus=yplus, Lref=Lref,
                  correctionMultiCorpsF42=correctionMultiCorpsF42, heightMaxF42=heightMaxF42,
                  tbFilament=tbFilament, cleanCellN=cleanCellN)
    if verbose: printTimeAndMemory__('compute wall distance', time=python_time.time()-pt0)

    #===================
    # STEP 3 : BLANKING IBM
    #===================
    if verbose: pt0 = python_time.time(); printTimeAndMemory__('blank by IBC bodies', time=-1)
    _blankingIBM(t, tb, dimPb=dimPb, frontType=frontType, IBCType=IBCType, depth=depth,
                 Reynolds=Reynolds, yplus=yplus, Lref=Lref, twoFronts=twoFronts,
                 heightMaxF42=heightMaxF42, correctionMultiCorpsF42=correctionMultiCorpsF42,
                 wallAdaptF42=wallAdaptF42, blankingF42=blankingF42,
                 tbFilament=tbFilament, cleanCellN=cleanCellN)
    ##Keep for now to check
    #filamentBases=filamentBases, isFilamentOnly=isFilamentOnly, tbFilament=tbFilament,
    #isWireModel=isWireModel)
    Cmpi.barrier()
    _redispatch__(t=t)
    if verbose: printTimeAndMemory__('blank by IBC bodies', time=python_time.time()-pt0)
    #===================
    # STEP 4 : INTERP DATA CHIM
    #===================
    if verbose: pt0 = python_time.time(); printTimeAndMemory__('compute interpolation data (Abutting & Chimera)', time=-1)
    tc = C.node2Center(t)

    if Internal.getNodeFromType(t, "GridConnectivity1to1_t") is not None:
        Xmpi._setInterpData(t, tc, nature=1, loc='centers', storage='inverse', sameName=1, dim=dimPb, itype='abutting', order=2, cartesian=cartesian)
    Xmpi._setInterpData(t, tc, nature=1, loc='centers', storage='inverse', sameName=1, sameBase=1, dim=dimPb, itype='chimera', order=2, cartesian=cartesian)
    if verbose: printTimeAndMemory__('compute interpolation data (Abutting & Chimera)', time=python_time.time()-pt0)

    #===================
    # STEP 4 : BUILD FRONT
    #===================
    if verbose: pt0 = python_time.time(); printTimeAndMemory__('build IBM front', time=-1)
    t, tc, front, front2, frontWMM = buildFrontIBM(t, tc, tb=tb, dimPb=dimPb, frontType=frontType,
                                                   cartesian=cartesian, twoFronts=twoFronts, check=check,
                                                   tbFilament=tbFilament)
    if verbose: printTimeAndMemory__('build IBM front', time=python_time.time()-pt0)

    #===================
    # STEP 5 : INTERP DATA IBM
    #===================
    if verbose: pt0 = python_time.time(); printTimeAndMemory__('compute interpolation data (IBM)', time=-1)
    _setInterpDataIBM(t, tc, tb, front, front2=front2, dimPb=dimPb, frontType=frontType, IBCType=IBCType, depth=depth,
                      Reynolds=Reynolds, yplus=yplus, Lref=Lref,
                      cartesian=cartesian, twoFronts=twoFronts, check=check,
                      tbFilament=tbFilament, frontWMM=frontWMM)
    if verbose: printTimeAndMemory__('compute interpolation data (IBM)', time=python_time.time()-pt0)
    #===================
    # STEP 6 : INIT IBM
    #===================
    if verbose: pt0 = python_time.time(); printTimeAndMemory__('initialize and clean', time=-1)

    t, tc, tc2 = initializeIBM(t, tc, tb, tinit=tinit, tbCurvi=tbCurvi, dimPb=dimPb, twoFronts=twoFronts,
                               tbFilament=tbFilament, cleanCellN=cleanCellN)

    _redispatch__(t=t, tc=tc, tc2=tc2)

    _setInjOutlet__(tc, tb)

    if isinstance(tc_out, str):
        if cartesian: tcp = Compressor.compressCartesian(tc)
        else: tcp = tc
        Cmpi.convertPyTree2File(tcp, tc_out, ignoreProcNodes=True)

        if tc2:
            if cartesian: tcp2 = Compressor.compressCartesian(tc2)
            else: tcp2 = tc2
            tc2_out = tc_out.replace('tc', 'tc2') if 'tc' in tc_out else 'tc2.cgns'
            Cmpi.convertPyTree2File(tcp2, tc2_out, ignoreProcNodes=True)

    if isinstance(t_out, str):
        if cartesian: tp = Compressor.compressCartesian(t)
        else: tp = t
        Cmpi.convertPyTree2File(tp, t_out, ignoreProcNodes=True)

    _computeMeshInfo(t)

    Cmpi.barrier()
    if verbose: printTimeAndMemory__('initialize and clean', time=python_time.time()-pt0)

    if tc2 is not None: return t, tc, tc2
    else: return t, tc

def prepareIBMDataExtrude(t_case, t_out, tc_out, t, to=None,
                          depth=2, frontType=1, octreeMode=0, IBCType=1,
                          verbose=True, check=False, twoFronts=False, cartesian=True,
                          yplus=100., Lref=1., correctionMultiCorpsF42=False, blankingF42=False, wallAdaptF42=None, heightMaxF42=-1.,
                          tbox=None, extrusion='cart'):
    import Generator.IBM as G_IBM
    import time as python_time

    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = Internal.copyTree(t_case)

    ##Will be changed in the near future - pending commit #165
    #cartesian = True
    #if t:
    #    cartesian = G_IBM.checkCartesian(t, nghost=2)
    #    if cartesian:
    #        RED  = "\033[1;31;40m"
    #        END  = "\033[0m"
    #        print("===========================================")
    #        print("Note: t_in is a " + RED + "CARTESIAN " + END + "grid")
    #        print("===========================================")
    #    else:
    #        RED  = "\033[1;31;40m"
    #        END  = "\033[0m"
    #        print("===========================================")
    #        print("Note: t_in is " + RED + "NOT" + END + " a " + RED + "CARTESIAN " + END + "grid")
    #        print("===========================================")

    refstate = Internal.getNodeFromName(tb, 'ReferenceState')
    flowEqn  = Internal.getNodeFromName(tb, 'FlowEquationSet')

    Reynolds = Internal.getNodeFromName(tb, 'Reynolds')
    if Reynolds is not None:
        Reynolds = Internal.getValue(Reynolds)
        if Reynolds < 1.e5: frontType = 1
    else: Reynolds = 1.e6

    expand = 3 if frontType != 42 else 4

    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    if dimPb is None: raise ValueError('prepareIBMDataPara: EquationDimension is missing in input geometry tree.')
    dimPb = Internal.getValue(dimPb)
    if dimPb == 2: C._initVars(tb, 'CoordinateZ', 0.)

    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('prepareIBMDataPara: GoverningEquations is missing in input geometry tree.')
    model = Internal.getValue(model)

    ibctypes = Internal.getNodesFromName(tb, 'ibctype')
    if ibctypes is None: raise ValueError('prepareIBMDataPara: ibc type is missing in input geometry tree.')
    ibctypes = list(set(Internal.getValue(ibc) for ibc in ibctypes))

    if model == 'Euler':
        if any(ibc in ['Musker', 'MuskerMob', 'Mafzal', 'Log', 'TBLE', 'TBLE_FULL'] for ibc in ibctypes):
            raise ValueError("prepareIBMDataPara: governing equations (Euler) not consistent with ibc types %s"%(ibctypes))

    #===================
    # STEP 0 : GET FILAMENT BODIES
    #===================
    tb, tbFilament              = D_IBM.determineClosedSolidFilament__(tb)
    isFilamentOnly, isWireModel = D_IBM.localWMMFlags__(tb, tbFilament)

    if extrusion=='cyl': cartesian = False
    cellN = Internal.getNodeFromName(t, 'cellN')
    if cellN is not None: C._initVars(t, '{centers:cellN}=minimum(1, {centers:cellN})') # modification needed for extrude

    #===================
    # STEP 1 : GENERATE MESH
    #===================
    ##SKIPPED - mesh is provided as an input

    #===================
    # STEP 2 : DIST2WALL
    #===================
    ##SKIPPED - mesh is provided as an input & has TurbulentDistance already

    #===================
    # STEP 3 : BLANKING IBM
    #===================
    if verbose: pt0 = python_time.time(); printTimeAndMemory__('blank by IBC bodies', time=-1, functionName='prepareIBMDataExtrude')
    _blankingIBM(t, tb, dimPb=dimPb, frontType=frontType, IBCType=IBCType, depth=depth,
                 Reynolds=Reynolds, yplus=yplus, Lref=Lref, twoFronts=twoFronts,
                 heightMaxF42=heightMaxF42, correctionMultiCorpsF42=correctionMultiCorpsF42,
                 wallAdaptF42=wallAdaptF42, blankingF42=blankingF42,
                 tbFilament=tbFilament)

    ##set the kmin et kmax Ghost cells are potential donors                                          #__
    listvars_local =['cellNChim','cellNIBC']                                                         #  |
    for z in Internal.getZones(t):                                                                   #  |
        sol            = Internal.getNodeFromName(z,'FlowSolution#Centers')                          #  |
        for var in listvars_local:                                                                   #  |
            cellN          = Internal.getNodeFromName(sol,var)[1]                                    #  | Modification needed for extrude &
            sh             = numpy.shape(cellN)                                                      #  | to replicate previous behavior
            for k in [0,1, sh[2]-2, sh[2]-1]:                                                        #  |
                for j in range(sh[1]):                                                               #  |
                    for i in range(sh[0]):                                                           #  |
                        if  cellN[i,j,k] != 0:  cellN[i,j,k] =1                                      #  |
    C._initVars(t,'{centers:cellN}=maximum(0.,{centers:cellNChim})')# vaut -3, 0, 1, 2 initialement  #__
    Cmpi.barrier()
    _redispatch__(t=t)
    if verbose: printTimeAndMemory__('blank by IBC bodies', time=python_time.time()-pt0, functionName='prepareIBMDataExtrude')
    #===================
    # STEP 4 : INTERP DATA CHIM
    #===================
    ## REQUIREMENT:: cellN mush be correct here      --> _setInterpData uses cellN
    ##               if cellN* is correct henceforth --> correct values at the end of prepareIBMDataExtrude
    if verbose: pt0 = python_time.time(); printTimeAndMemory__('compute interpolation data (Abutting & Chimera)', time=-1, functionName='prepareIBMDataExtrude')
    tc = C.node2Center(t)
    if Internal.getNodeFromType(t, "GridConnectivity1to1_t") is not None:
        Xmpi._setInterpData(t, tc, nature=1, loc='centers', storage='inverse', sameName=1, dim=dimPb, itype='abutting', order=2, cartesian=cartesian)
    Xmpi._setInterpData(t, tc, nature=1, loc='centers', storage='inverse', sameName=1, sameBase=1, dim=dimPb, itype='chimera', order=2, cartesian=cartesian)
    if verbose: printTimeAndMemory__('compute interpolation data (Abutting & Chimera)', time=python_time.time()-pt0, functionName='prepareIBMDataExtrude')
    #===================
    # STEP 5 : BUILD FRONT
    #===================
    if verbose: pt0 = python_time.time(); printTimeAndMemory__('build IBM front', time=-1, functionName='prepareIBMDataExtrude')
    t, tc, front, front2, frontWMM = buildFrontIBM(t, tc, tb=tb, dimPb=dimPb, frontType=frontType,
                                                   cartesian=cartesian, twoFronts=twoFronts, check=check,
                                                   tbFilament=tbFilament)
    if verbose: printTimeAndMemory__('build IBM front', time=python_time.time()-pt0, functionName='prepareIBMDataExtrude')
    #===================
    # STEP 6 : INTERP DATA IBM
    #===================
    if verbose: pt0 = python_time.time(); printTimeAndMemory__('compute interpolation data (IBM)', time=-1, functionName='prepareIBMDataExtrude')
    _setInterpDataIBM(t, tc, tb, front, front2=front2, dimPb=dimPb, frontType=frontType, IBCType=IBCType, depth=depth,
                      Reynolds=Reynolds, yplus=yplus, Lref=Lref,
                      cartesian=cartesian, twoFronts=twoFronts, check=check,
                      tbFilament=tbFilament, frontWMM=frontWMM)
    if verbose: printTimeAndMemory__('compute interpolation data (IBM)', time=python_time.time()-pt0, functionName='prepareIBMDataExtrude')
    #===================
    # STEP 7 : INIT IBM
    #===================
    if verbose: pt0 = python_time.time(); printTimeAndMemory__('initialize and clean', time=-1, functionName='prepareIBMDataExtrude')
    tsave = Internal.copyTree(t)    # Modification needed to by pass the initialization of t in the macro function initializeIBM
    t     = None                    #
    t, tc, tc2 = initializeIBM(t, tc, tb, dimPb=dimPb, twoFronts=twoFronts, tbFilament=tbFilament)
    t = Internal.copyTree(tsave)    # Modification needed to by pass the initialization of t in the macro function initializeIBM
    _redispatch__(t=t, tc=tc, tc2=tc2)

    if extrusion == 'cyl':                                                                              #__
        T._cyl2Cart(t, (0,0,0),(1,0,0))                                                                 #  |
        T._cyl2Cart(tc,(0,0,0),(1,0,0))                                                                 #  |
        # modif info mesh in zonesubregion_t                                                            #  |
        for z in Internal.getZones(tc):                                                                 #  |
            for zsr in Internal.getNodesFromType(z, "ZoneSubRegion_t"):                                 #  |
                zsrname = Internal.getName(zsr)                                                         #  |
                zsrname = zsrname.split('_')                                                            #  |Modification needed for extrude in cylindrical coordinates
                if zsrname[0]=='IBCD':                                                                  #  |
                    for var in ['C','W','I']:                                                           #  |
                        r     = Internal.getNodeFromName(zsr,'CoordinateY_P'+var)[1]                    #  |
                        theta = Internal.getNodeFromName(zsr,'CoordinateZ_P'+var)[1]                    #  |
                        for l in range(numpy.size(r)):                                                  #  |
                            yy  = r[l]*numpy.cos( theta[l] )                                            #  |
                            zz  = r[l]*numpy.sin( theta[l] )                                            #  |
                            r[l]= yy; theta[l] = zz                                                     #  |
                            #__

    if isinstance(tc_out, str):
        tcp = Compressor.compressCartesian(tc)
        Cmpi.convertPyTree2File(tcp, tc_out, ignoreProcNodes=True)

        if twoFronts:
            tcp2 = Compressor.compressCartesian(tc2)
            tc2_out = tc_out.replace('tc', 'tc2') if 'tc' in tc_out else 'tc2.cgns'
            Cmpi.convertPyTree2File(tcp2, tc2_out, ignoreProcNodes=True)

    if isinstance(t_out, str):
        tp = Compressor.compressCartesian(t)
        Cmpi.convertPyTree2File(tp, t_out, ignoreProcNodes=True)

    _computeMeshInfo(t)

    if Cmpi.size > 1: Cmpi.barrier()
    if verbose: printTimeAndMemory__('initialize and clean', time=python_time.time()-pt0, functionName='prepareIBMDataExtrude')

    if tc2 is not None: return t, tc, tc2
    else: return t, tc

def prepareIBMDataAdapt(t_case, t_out, tc_out, t_in,
                        depth=2, IBCType=1, verbose=True,
                        check=False, twoFronts=False, cartesian=True,
                        yplus=100., Lref=1., correctionMultiCorpsF42=False, blankingF42=False, wallAdaptF42=None, heightMaxF42=-1.):

    import time as python_time

    if t_in is None: raise ValueError('prepareIBMDataAdapt: input tree is missing.')

    if isinstance(t_in, str): t = Cmpi.convertFile2PyTree(t_in, proc=Cmpi.rank)
    else: t = Internal.copyTree(t_in)

    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = Internal.copyTree(t_case)

    frontType = 42

    refstate = Internal.getNodeFromName(tb, 'ReferenceState')
    flowEqn  = Internal.getNodeFromName(tb, 'FlowEquationSet')

    Reynolds = Internal.getNodeFromName(tb, 'Reynolds')
    if Reynolds is not None:
        Reynolds = Internal.getValue(Reynolds)
        if Reynolds < 1.e5: frontType = 1
    else: Reynolds = 1.e6

    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    if dimPb is None: raise ValueError('prepareIBMDataAdapt: EquationDimension is missing in input geometry tree.')
    dimPb = Internal.getValue(dimPb)
    if dimPb == 2: C._initVars(tb, 'CoordinateZ', 0.)

    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('prepareIBMDataAdapt: GoverningEquations is missing in input geometry tree.')
    model = Internal.getValue(model)

    ibctypes = Internal.getNodesFromName(tb, 'ibctype')
    if ibctypes is None: raise ValueError('prepareIBMDataAdapt: ibc type is missing in input geometry tree.')
    ibctypes = list(set(Internal.getValue(ibc) for ibc in ibctypes))

    if model == 'Euler':
        if any(ibc in ['Musker', 'MuskerMob', 'Mafzal', 'Log', 'TBLE', 'TBLE_FULL'] for ibc in ibctypes):
            raise ValueError("prepareIBMDataAdapt: governing equations (Euler) not consistent with ibc types %s"%(ibctypes))

    C._initVars(t, '{centers:TurbulentDistance}=(({centers:cellN} > 0)*1 + ({centers:cellN} == 0)*-1)*{centers:TurbulentDistance}')
    C._rmVars(t, 'centers:cellN')

    Internal._rmNodesFromName(tb,"SYM")

    #===================
    # STEP 3 : BLANKING IBM
    #===================
    if verbose: pt0 = python_time.time(); printTimeAndMemory__('blank by IBC bodies', time=-1)
    _blankingIBM(t, tb, dimPb=dimPb, frontType=frontType, IBCType=IBCType, depth=depth,
                 Reynolds=Reynolds, yplus=yplus, Lref=Lref, twoFronts=twoFronts,
                 heightMaxF42=heightMaxF42, correctionMultiCorpsF42=correctionMultiCorpsF42,
                 wallAdaptF42=wallAdaptF42, blankingF42=blankingF42)
    Cmpi.barrier()
    _redispatch__(t=t)
    if verbose: printTimeAndMemory__('blank by IBC bodies', time=python_time.time()-pt0)

    #===================
    # STEP 4 : INTERP DATA CHIM
    #===================
    if verbose: pt0 = python_time.time(); printTimeAndMemory__('compute interpolation data (Abutting & Chimera)', time=-1)
    tc = C.node2Center(t)

    if Internal.getNodeFromType(t, "GridConnectivity1to1_t") is not None:
        Xmpi._setInterpData(t, tc, nature=1, loc='centers', storage='inverse', sameName=1, dim=dimPb, itype='abutting', order=2, cartesian=cartesian)
    Xmpi._setInterpData(t, tc, nature=1, loc='centers', storage='inverse', sameName=1, sameBase=1, dim=dimPb, itype='chimera', order=2, cartesian=cartesian)
    if verbose: printTimeAndMemory__('compute interpolation data (Abutting & Chimera)', time=python_time.time()-pt0)

    #===================
    # STEP 4 : BUILD FRONT
    #===================
    if verbose: pt0 = python_time.time(); printTimeAndMemory__('build IBM front', time=-1)
    t, tc, front, front2, frontWMM = buildFrontIBM(t, tc, tb=tb, dimPb=dimPb, frontType=frontType,
                                                   cartesian=cartesian, twoFronts=twoFronts, check=check)
    if verbose: printTimeAndMemory__('build IBM front', time=python_time.time()-pt0)

    #===================
    # STEP 5 : INTERP DATA IBM
    #===================
    if verbose: pt0 = python_time.time(); printTimeAndMemory__('compute interpolation data (IBM)', time=-1)
    _setInterpDataIBM(t, tc, tb, front, front2=front2, dimPb=dimPb, frontType=frontType, IBCType=IBCType, depth=depth,
                      Reynolds=Reynolds, yplus=yplus, Lref=Lref,
                      cartesian=cartesian, twoFronts=twoFronts, check=check)
    if verbose: printTimeAndMemory__('compute interpolation data (IBM)', time=python_time.time()-pt0)

    #===================
    # STEP 6 : INIT IBM
    #===================
    if verbose: pt0 = python_time.time(); printTimeAndMemory__('initialize and clean', time=-1)

    _, tc, tc2 = initializeIBM(None, tc, tb, dimPb=dimPb, twoFronts=twoFronts)

    _redispatch__(t=t, tc=tc, tc2=tc2)

    _setInjOutlet__(tc, tb)

    if isinstance(tc_out, str):
        if cartesian: tcp = Compressor.compressCartesian(tc)
        else: tcp = tc
        Cmpi.convertPyTree2File(tcp, tc_out, ignoreProcNodes=True)

        if tc2:
            if cartesian: tcp2 = Compressor.compressCartesian(tc2)
            else: tcp2 = tc2
            tc2_out = tc_out.replace('tc', 'tc2') if 'tc' in tc_out else 'tc2.cgns'
            Cmpi.convertPyTree2File(tcp2, tc2_out, ignoreProcNodes=True)

    if isinstance(t_out, str):
        if cartesian: tp = Compressor.compressCartesian(t)
        else: tp = t
        Cmpi.convertPyTree2File(tp, t_out, ignoreProcNodes=True)

    _computeMeshInfo(t)

    Cmpi.barrier()
    if verbose: printTimeAndMemory__('initialize and clean', time=python_time.time()-pt0)

    if Cmpi.size > 1: Cmpi.barrier()

    if tc2 is not None: return t, tc, tc2
    else: return t, tc

def _redispatch__(t=None, tc=None, tc2=None):
    import Distributor2.Mpi as D2mpi

    if tc:
        algo = 'graph'
        tskel = Cmpi.convert2SkeletonTree(tc)
        Internal._rmNodesByType(tskel, 'ZoneSubRegion_t')
        tcs    = Cmpi.allgatherTree(tskel)
        stats  = D2._distribute(tcs, Cmpi.size, algorithm=algo)

        D2._copyDistribution(tc, tcs)
        D2mpi._redispatch(tc, verbose=1)

        if t is not None:
            D2._copyDistribution(t, tcs)
            D2mpi._redispatch(t)
        if tc2:
            D2._copyDistribution(tc2, tcs)
            D2mpi._redispatch(tc2)
        del tcs
    else:
        algo = 'fast'
        tskel = Cmpi.convert2SkeletonTree(t)
        ts    = Cmpi.allgatherTree(tskel)
        stats = D2._distribute(ts, Cmpi.size, algorithm=algo)
        D2._copyDistribution(t, ts)
        del ts
        D2mpi._redispatch(t, verbose=1)

    return None

# alias prepareIBMData new version
prepareIBMDataPara = prepareIBMData

#=========================================================================
# Compute the wall distance for IBM pre-processing.
# IN: t (tree): computational tree
# IN: tb (tree): geometry tree (IBM bodies)
# IN: dimPb (2 or 3): problem dimension
# IN: frontType (0,1,2 or 42): type of IBM front
# IN: Reynolds (float): Reynolds number (F42)
# IN: yplus (float): estimated yplus at the first computed cells (F42)
# IN: Lref (float): characteristic length of the geometry (F42)
# IN: correctionMultiCorpsF42 (boolean): if True, computes the wall distance
# w.r.t each body that is not a symmetry plane (F42)
# IN: heightMaxF42 (float): if heightMaxF42 > 0: uses a maximum modeling height
# to speed up individual wall distance calculations when correctionMultiCorpsF42 is active (F42)
# OUT: centers:TurbulentDistance field
# OUT: (optional) centers:TurbulentDistance_body%i fields
#=========================================================================
def dist2wallIBM(t, tb, dimPb=3, frontType=1, Reynolds=1.e6, yplus=100, Lref=1.,
                 correctionMultiCorpsF42=False, heightMaxF42=-1., tbFilament=None, cleanCellN=True):
    """Compute the wall distance for IBM pre-processing."""
    tp = Internal.copyRef(t)
    _dist2wallIBM(t, tb, dimPb=dimPb, frontType=frontType, Reynolds=Reynolds, yplus=yplus, Lref=Lref,
                  correctionMultiCorpsF42=correctionMultiCorpsF42, heightMaxF42=heightMaxF42,
                  tbFilament=tbFilament, cleanCellN=cleanCellN)
    return tp

def _dist2wallIBMFilamentWMM__(t, tb2, tbsave, dimPb):
    tbFilamentnoWMM = []
    tbFilamentWMM   = []
    for z in Internal.getZones(tb2):
        soldef  = Internal.getNodeFromName(z,'.Solver#define')
        ibctype = Internal.getNodeFromName(soldef,'ibctype')
        if Internal.getValue(ibctype) != "wiremodel":
            tbFilamentnoWMM.append(z)
        else:
            tbFilamentWMM.append(z)

    if tbFilamentnoWMM:
        tbFilamentnoWMM = C.newPyTree(['BasenoWMM', tbFilamentnoWMM])
        tbsave          = Internal.merge([tbsave,tbFilamentnoWMM])
        C._initVars(t,'{centers:TurbulentDistance}=1e06')
        DTW._distance2Walls(t, tbFilamentnoWMM, type='ortho', signed=0, dim=dimPb, loc='centers')
        C._initVars(t,'{centers:TurbulentDistanceFilament}={centers:TurbulentDistance}')
        C._initVars(t,'{centers:TurbulentDistance}=minimum({centers:TurbulentDistanceSolid},{centers:TurbulentDistanceFilament})')

    if tbFilamentWMM:
        C._initVars(t,'{centers:TurbulentDistance}=1e06')
        DTW._distance2Walls(t, tbFilamentWMM, type='ortho', signed=0, dim=dimPb, loc='centers')
        C._initVars(t,'{centers:TurbulentDistanceFilamentWMM}={centers:TurbulentDistance}')
        C._initVars(t,'{centers:TurbulentDistance}=minimum({centers:TurbulentDistanceSolid},{centers:TurbulentDistanceFilamentWMM})')
        if tbFilamentnoWMM:
            C._initVars(t,'{centers:TurbulentDistance}=minimum({centers:TurbulentDistance},{centers:TurbulentDistanceFilament})')

    C._initVars(t,"{centers:TurbulentDistanceSolid}=({centers:TurbulentDistanceSolid}>1e03)*0+({centers:TurbulentDistanceSolid}<1e03)*{centers:TurbulentDistanceSolid}")
    C._initVars(t,"{centers:TurbulentDistanceFilament}=({centers:TurbulentDistanceFilament}>1e03)*0+({centers:TurbulentDistanceFilament}<1e03)*{centers:TurbulentDistanceFilament}")
    C._initVars(t,"{centers:TurbulentDistanceFilamentWMM}=({centers:TurbulentDistanceFilamentWMM}>1e03)*0+({centers:TurbulentDistanceFilamentWMM}<1e03)*{centers:TurbulentDistanceFilamentWMM}")
    return None

def _dist2wallIBM(t, tb, dimPb=3, frontType=1, Reynolds=1.e6, yplus=100, Lref=1.,
                  correctionMultiCorpsF42=False, heightMaxF42=-1., tbFilament=None, cleanCellN=True):
    """Compute the wall distance for IBM pre-processing."""

    isFilamentOnly, isWireModel = D_IBM.localWMMFlags__(tb, tbFilament)

    if dimPb == 2:
        # Set CoordinateZ to a fixed value
        z0 = Internal.getNodeFromType2(t, "Zone_t")
        bb0 = G.bbox(z0); dz = (bb0[5]-bb0[2])*0.5
        tb2 = Internal.copyRef(tb)
        tb2 = Internal.copyValue(tb2, byName='CoordinateZ')
        C._initVars(tb2, 'CoordinateZ', dz)
    else:
        tb2 = tb

    DTW._distance2Walls(t, tb2, type='ortho', signed=0, dim=dimPb, loc='centers')

    tbsave = tb2

    if tbFilament and not isFilamentOnly:
        C._initVars(t,'{centers:TurbulentDistanceSolid}={centers:TurbulentDistance}')
        if dimPb ==2: tb2 = C.initVars(tbFilament, 'CoordinateZ', dz)
        _dist2wallIBMFilamentWMM__(t, tb2, tbsave, dimPb)
    elif isFilamentOnly:
        C._initVars(t,"{centers:TurbulentDistanceSolid}=0")
        C._initVars(t,"{centers:TurbulentDistanceFilament}={centers:TurbulentDistance}")
        C._initVars(t,"{centers:TurbulentDistanceFilamentWMM}={centers:TurbulentDistance}")


    # Compute turbulentdistance wrt each body that is not a sym plan (centers:TurbulentDistance_bodyX)
    if correctionMultiCorpsF42 and frontType == 42:
        print("correctionMultiCorpsF42")
        # Keep track of the general turbulentDistance
        C._initVars(t,'{centers:TurbulentDistance_ori}={centers:TurbulentDistance}')

        if yplus > 0.:
            shiftDist = G_IBM_Height.computeModelisationHeight(Re=Reynolds, yplus=yplus, L=Lref)
        else:
            snears    = Internal.getNodesFromName(tb, 'snear')
            h         = max(snears, key=lambda x: x[1])[1]
            shiftDist = G_IBM_Height.computeBestModelisationHeight(Re=Reynolds, h=h) # meilleur compromis entre hauteur entre le snear et la hauteur de modelisation

        if heightMaxF42 > 0.:
            if shiftDist > heightMaxF42: shiftDist = heightMaxF42

        for z in Internal.getZones(t):
            cptBody = 1
            if dimPb == 3: tb2 = tb
            for body in Internal.getNodesFromType(tb2, 'Zone_t'):
                if body[0] != "sym" and ("closure" not in body[0]):
                    # Create extended BBox around each body
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

        C._initVars(t, '{centers:TurbulentDistance}={centers:TurbulentDistance_ori}')
        C._rmVars(t, ['centers:TurbulentDistance_ori'])

    if not cleanCellN: C._initVars(t, '{centers:TurbulentDistanceAllBC}={centers:TurbulentDistance}')

    return None

#=========================================================================
# Blank t by IBC bodies for IBM pre-processing.
# IN: t (tree): computational tree
# IN: tb (tree): geometry tree (IBM bodies)
# IN: dimPb (2 or 3): problem dimension
# IN: frontType (0,1,2 or 42): type of IBM front
# IN: IBCType (-1 or 1): type of IBM, -1: IB target points are located inside the solid, 1: IB target points are located in the fluid
# IN: depth (int): depth of overlapping regions
# IN: Reynolds (float): Reynolds number (F42)
# IN: yplus (float): estimated yplus at the first computed cells (F42)
# IN: Lref (float): characteristic length of the geometry (F42)
# IN: twoFronts (boolean): if True, performs the IBM pre-processing for an additional image point positioned farther away
# IN: correctionMultiCorpsF42 (boolean): if True, ensures that there are calculated points between
# the immersed bodies by using individual wall distances (F42)
# IN: blankingF42 (boolean): if True, reduces as much as possible the number of IB target points inside the boundary layer (F42)
# IN: wallAdaptF42 (cloud of IB target points with yplus information):
# use a previous computation to adapt the positioning of IB target points around the geometry according to a target yplus (F42)
# IN: heightMaxF42 (float): maximum modeling height for the location of IB target points around the geometry (F42)
# IN: filamentBases: list of bases that are IBC filaments
# IN: isFilamentOnly: boolean on whether there is only a IBC filament
# IN: tbFilament: PyTree of the IBC filaments
# OUT: centers:cellN, centers:cellNIBC, centers:cellNChim, centers:cellNFront fields
# OUT: (optional) centers:cellNIBC_2, centers:cellNFront_2 fields for second image points
#=========================================================================
def _blankingIBM__(t, tb, dimPb=3, frontType=1, IBCType=1, depth=2, Reynolds=1.e6, yplus=100, Lref=1.,
                   correctionMultiCorpsF42=False, blankingF42=False, wallAdaptF42=None, heightMaxF42=-1.,
                   tbFilament=None, cleanCellN=True):

    isFilamentOnly, isWireModel = D_IBM.localWMMFlags__(tb, tbFilament)

    isSkipDist = False
    minval     = C.getMinValue(t, 'centers:TurbulentDistance');
    minval     = Cmpi.allreduce(minval, op=Cmpi.MIN)
    if minval<0: isSkipDist=True

    snear_min = 10.e10
    for z in Internal.getZones(tb):
        sdd = Internal.getNodeFromName1(z, ".Solver#define")
        if sdd is not None:
            snearl = Internal.getNodeFromName1(sdd, "snear")
            if snearl is not None:
                snearl = Internal.getValue(snearl)
                snear_min = min(snear_min, snearl)
    snear_min = Cmpi.allreduce(snear_min, op=Cmpi.MIN)

    cellNIBC_blank = Internal.getNodeFromName(t, 'cellNIBC_blank')
    if cellNIBC_blank is None:
        if not isFilamentOnly: _blankByIBCBodies(t, tb, 'centers', dimPb) #cellN -> 0 outside domain ; 1 inside domain
        if not cleanCellN: C._initVars(t, '{centers:cellNIBC_blank}={centers:cellN}')
    else:
        print("Info: blankingIBM: cellNIBC_blank field found in t. blankByIBCBodies is skipped.")
        C._initVars(t, '{centers:cellN}={centers:cellNIBC_blank}')

    C._initVars(t, '{centers:cellNIBC}={centers:cellN}')
    if not isFilamentOnly and not isSkipDist: _signDistance(t)
    C._initVars(t,'{centers:cellN}={centers:cellNIBC}')

    if tbFilament:
        C._initVars(t,'{centers:cellNFil}={centers:cellN}')
        C._initVars(t,'{centers:cellNFilWMM}={centers:cellN}')

    # determination des pts IBC
    cellNIBC_hole = Internal.getNodeFromName(t, 'cellNIBC_hole')
    if cellNIBC_hole is None:
        if frontType != 42:
            if IBCType == -1: X._setHoleInterpolatedPoints(t,depth=-depth,dir=0,loc='centers',cellNName='cellN',addGC=False)
            elif IBCType == 1:
                X._setHoleInterpolatedPoints(t,depth=1,dir=1,loc='centers',cellNName='cellN',addGC=False) # pour les gradients
                if frontType < 2:
                    X._setHoleInterpolatedPoints(t,depth=depth,dir=0,loc='centers',cellNName='cellN',addGC=False)
                else:
                    depthL = depth+1
                    X._setHoleInterpolatedPoints(t,depth=depthL,dir=0,loc='centers',cellNName='cellN',addGC=False)
                    #cree des pts extrapoles supplementaires
                    # _blankClosestTargetCells(t,cellNName='cellN', depth=depthL)
            else:
                raise ValueError('prepareIBMData: not valid IBCType. Check model.')
        else:
            # F42: tracking of IB points using distance information
            # the classical algorithm (front 1) is first used to ensure a minimum of two rows of target points around the geometry
            C._initVars(t,'{centers:cellNMin}={centers:cellNIBC}')
            if IBCType == -1: X._setHoleInterpolatedPoints(t,depth=-depth,dir=0,loc='centers',cellNName='cellNMin',addGC=False)
            elif IBCType == 1: X._setHoleInterpolatedPoints(t,depth=1,dir=1,loc='centers',cellNName='cellNMin',addGC=False) # pour les gradients
            X._setHoleInterpolatedPoints(t,depth=depth,dir=0,loc='centers',cellNName='cellNMin',addGC=False)

            for z in Internal.getZones(t):
                if Internal.getNodeFromName(z, 'SaveF1'):
                    Internal._rmNode(z, Internal.getNodeFromName(z, '.Solver#defineTMP'))
                    continue

                h = abs(C.getValue(z,'CoordinateX',0)-C.getValue(z,'CoordinateX',1))
                if yplus > 0.:
                    height = G_IBM_Height.computeModelisationHeight(Re=Reynolds, yplus=yplus, L=Lref)
                else:
                    height = G_IBM_Height.computeBestModelisationHeight(Re=Reynolds, h=h) # meilleur compromis entre hauteur entre le snear et la hauteur de modelisation
                    yplus  = G_IBM_Height.computeYplus(Re=Reynolds, height=height, L=Lref)
                if heightMaxF42 > 0.:
                    if height > heightMaxF42:
                        height = heightMaxF42
                        #print("Snear min (SM) = %g || Wall Modeling Height (WMH) = %g || WMH/SM = %g"%(snear_min,height,height/snear_min))
                C._initVars(z,'{centers:cellN}=({centers:TurbulentDistance}>%20.16g)+(2*({centers:TurbulentDistance}<=%20.16g)*({centers:TurbulentDistance}>0))'%(height,height))

                if correctionMultiCorpsF42:
                    # Prevent different body modeling from overlapping -> good projection of image points in the wall normal direction

                    cptBody = 1
                    for body in Internal.getNodesFromType(tb,'Zone_t'):
                        if body[0] != "sym" and ("closure" not in body[0]):
                            cptBody += 1

                    epsilon_dist = 2*(abs(C.getValue(z,'CoordinateX',1)-C.getValue(z,'CoordinateX',0)))
                    max_dist     = 2*0.1*Lref

                    # Try to find the best route between two adjacent bodies by finding optimal iso distances
                    def correctionMultiCorps(cellN, cellNF):
                        if cellN == 2 and cellNF == 2: return 1
                        return cellN

                    def findIsoFront(cellNFront, Dist_1, Dist_2):
                        if Dist_1 < max_dist and Dist_2 < max_dist:
                            if abs(Dist_1-Dist_2) < epsilon_dist: return 2
                        return max(cellNFront,1)

                    for i in range(1, cptBody):
                        for j in range(1, cptBody):
                            if j != i:
                                C._initVars(z, 'centers:cellNFrontIso', findIsoFront, ['centers:cellNFrontIso', 'centers:TurbulentDistance_body%i'%i, 'centers:TurbulentDistance_body%i'%j])

                    C._initVars(z, 'centers:cellN', correctionMultiCorps, ['centers:cellN', 'centers:cellNFrontIso'])

                    for i in range(1, cptBody):
                        C._rmVars(z,['centers:cellN_body%i'%i, 'centers:TurbulentDistance_body%i'%i])

            if wallAdaptF42 is not None:
                # Use previous computation to adapt the positioning of IB points around the geometry (impose y+PC <= y+ref)
                # Warning: the wallAdapt file has to be obtained with TIBM.createWallAdapt(tc)
                C._initVars(t,'{centers:yplus}=100000.')
                if isinstance(wallAdaptF42, str): w = C.convertFile2PyTree(wallAdaptF42)
                else: w = wallAdaptF42
                zones = Internal.getZones(t)
                total = len(zones)
                cpt = 1
                for z in zones:
                    print("Info: blankingIBM: modeling height adaptation: zone %d / %d"%(cpt, total))
                    cellN = Internal.getNodeFromName(z, 'cellN')[1]
                    if 2 in cellN:
                        hloc = abs(C.getValue(z,'CoordinateX',1)-C.getValue(z,'CoordinateX',0))
                        zname = z[0]
                        zd = Internal.getNodeFromName(w, zname)
                        if zd is not None:
                            yplus_w = Internal.getNodeFromName(zd, 'yplus')[1]
                            listIndices = Internal.getNodeFromName(zd, 'PointListDonor')[1]

                            n = numpy.shape(yplus_w)[0]
                            yplusA = Converter.array('yplus', n, 1, 1)
                            yplusA[1][:] = yplus_w

                            C._setPartialFields(z, [yplusA], [listIndices], loc='centers')
                            C._initVars(z,'{centers:yplus}={centers:yplus}*(1-%20.16g/{centers:TurbulentDistance})'%(hloc)) #safety measure

                    cpt += 1

                C._initVars(t,'{centers:cellN}=({centers:cellN}>0) * ( (({centers:cellN}) * ({centers:yplus}<=%20.16g)) + ({centers:yplus}>%20.16g) )'%(yplus,yplus))

            # final security gate, we ensure that we have at least to layers of target points
            C._initVars(t, '{centers:cellN} = maximum({centers:cellN}, {centers:cellNMin})')
            C._rmVars(t,['centers:yplus', 'centers:cellNMin'])

            # propagate max yplus between procs
            yplus = numpy.array([float(yplus)])
            yplus = Cmpi.allreduce(yplus, op=Cmpi.MAX)[0]

            # Only keep the layer of target points useful for solver iterations, particularly useful in 3D
            if blankingF42:
                if wallAdaptF42 is None: X._maximizeBlankedCells(t, depth=2, addGC=False)
                else: print("Info: blankingIBM: blankingF42 cannot operate with a local modeling height")
        if not cleanCellN: C._initVars(t, '{centers:cellNIBC_hole}={centers:cellN}')
    else:
        print("Info: blankingIBM: cellNIBC_hole field found in t. Determination of IBC points is skipped.")
        C._initVars(t, '{centers:cellN}={centers:cellNIBC_hole}')

    if tbFilament:
        tbFilamentWMM   = []
        for z in Internal.getZones(tbFilament):
            soldef  = Internal.getNodeFromName(z,'.Solver#define')
            ibctype = Internal.getNodeFromName(soldef,'ibctype')
            if Internal.getValue(ibctype) == "wiremodel":
                tbFilamentWMM.append(z)
        maxy=C.getMaxValue(tbFilamentWMM, ['CoordinateY']);
        miny=C.getMinValue(tbFilamentWMM, ['CoordinateY']);

        if dimPb == 3:
            maxz=C.getMaxValue(tbFilamentWMM, ['CoordinateZ']);
            minz=C.getMinValue(tbFilamentWMM, ['CoordinateZ']);

        for z in Internal.getZones(t):
            if C.getMaxValue(z, 'centers:TurbulentDistanceFilament')> 1e-05:
                sol     = Internal.getNodeByName(z,"FlowSolution#Centers")
                cellNFil= Internal.getNodeByName(sol,'cellNFil')[1]
                cellN   = Internal.getNodeByName(sol,'cellN')[1]
                dist    = Internal.getNodeByName(sol,'TurbulentDistanceFilament')[1]
                ycord   = Internal.getNodeByName(z,'CoordinateY')[1]
                h       = abs(C.getValue(z,'CoordinateX',4)-C.getValue(z,'CoordinateX',5))
                sh      = numpy.shape(dist)
                if dimPb == 3:
                    for k in range(sh[2]):
                        for j in range(sh[1]):
                            for i in range(sh[0]):
                                valy=0.5*(ycord[i,j,k]+ycord[i,j+1,k])
                                if dist[i,j,k]<numpy.sqrt(8)*h and cellN[i,j,k]==1:
                                    cellNFil[i,j,k]=2
                                    cellN[i,j,k]=2
                else:
                    for j in range(sh[1]):
                        for i in range(sh[0]):
                            valy=0.5*(ycord[i,j]+ycord[i,j+1])
                            if dist[i,j]<numpy.sqrt(8)*h and cellN[i,j]==1:
                                cellNFil[i,j]=2
                                cellN[i,j]=2
            if C.getMaxValue(z, 'centers:TurbulentDistanceFilamentWMM')> 1e-05:
                sol     = Internal.getNodeByName(z,"FlowSolution#Centers")
                cellNFil= Internal.getNodeByName(sol,'cellNFilWMM')[1]
                cellN   = Internal.getNodeByName(sol,'cellN')[1]
                dist    = Internal.getNodeByName(sol,'TurbulentDistanceFilamentWMM')[1]
                ycord   = Internal.getNodeByName(z,'CoordinateY')[1]
                h       = abs(C.getValue(z,'CoordinateX',4)-C.getValue(z,'CoordinateX',5))
                sh      = numpy.shape(dist)
                if dimPb == 3:
                    zcord   = Internal.getNodeByName(z,'CoordinateZ')[1]
                    for k in range(sh[2]):
                        for j in range(sh[1]):
                            for i in range(sh[0]):
                                valy=0.5*(ycord[i,j,k]+ycord[i,j+1,k  ])
                                valz=0.5*(zcord[i,j,k]+zcord[i,j  ,k+1])
                                if dist[i,j,k]<numpy.sqrt(8)*h and \
                                   cellN[i,j,k]==1 and \
                                   valy<maxy and valy>miny and \
                                   valz<maxz and valz>minz :
                                    cellNFil[i,j,k]=2
                                    cellN[i,j,k]=2
                else:
                    for j in range(sh[1]):
                        for i in range(sh[0]):
                            valy=0.5*(ycord[i,j,0]+ycord[i,j+1,0])
                            if dist[i,j,0]<numpy.sqrt(8)*h and valy<maxy and valy>miny and cellN[i,j]==1:
                                cellNFil[i,j]=2
                                cellN[i,j]=2
        C._rmVars(t,['centers:TurbulentDistanceFilament','centers:TurbulentDistanceFilamentWMM'])
        if tbFilament: C._rmVars(t,['centers:TurbulentDistanceSolid'])
        if isFilamentOnly: C._initVars(t,'{centers:cellNFilWMM}={centers:cellNFil}')

    if not isFilamentOnly: _removeBlankedGrids(t, loc='centers')
    return None

def blankingIBM(t, tb, dimPb=3, frontType=1, IBCType=1, depth=2, Reynolds=1.e6, yplus=100, Lref=1., twoFronts=False,
                correctionMultiCorpsF42=False, blankingF42=False, wallAdaptF42=None, heightMaxF42=-1.,
                tbFilament=None, cleanCellN=True):
    """Blank the computational tree by IBC bodies for IBM pre-processing."""
    tp = Internal.copyRef(t)
    _blankingIBM(t, tb, dimPb=dimPb, frontType=frontType, IBCType=IBCType, depth=depth,
                 Reynolds=Reynolds, yplus=yplus, Lref=Lref, twoFronts=twoFronts,
                 heightMaxF42=heightMaxF42, correctionMultiCorpsF42=correctionMultiCorpsF42,
                 wallAdaptF42=wallAdaptF42, blankingF42=blankingF42,
                 tbFilament=tbFilament, cleanCellN=cleanCellN)
    return tp

def _blankingIBM(t, tb, dimPb=3, frontType=1, IBCType=1, depth=2, Reynolds=1.e6, yplus=100, Lref=1., twoFronts=False,
                 correctionMultiCorpsF42=False, blankingF42=False, wallAdaptF42=None, heightMaxF42=-1.,
                 tbFilament=None, cleanCellN=True):
    """Blank the computational tree by IBC bodies for IBM pre-processing."""

    isFilamentOnly, isWireModel = D_IBM.localWMMFlags__(tb, tbFilament)

    if dimPb == 2:
        T._addkplane(tb)
        T._contract(tb, (0,0,0), (1,0,0), (0,1,0), 0.01)
        if tbFilament:
            T._addkplane(tbFilament)
            T._contract(tbFilament, (0,0,0), (1,0,0), (0,1,0), 0.01)

    X._applyBCOverlaps(t, depth=depth, loc='centers', val=2, cellNName='cellN')
    C._initVars(t,'{centers:cellNChim}={centers:cellN}')
    C._initVars(t, 'centers:cellN', 1.)

    _blankingIBM__(t, tb, dimPb=dimPb, frontType=frontType, IBCType=IBCType, depth=depth,
                   Reynolds=Reynolds, yplus=yplus, Lref=Lref,
                   heightMaxF42=heightMaxF42, correctionMultiCorpsF42=correctionMultiCorpsF42,
                   wallAdaptF42=wallAdaptF42, blankingF42=blankingF42,
                   tbFilament=tbFilament, cleanCellN=cleanCellN)

    C._initVars(t, '{centers:cellNIBC}={centers:cellN}')
    if IBCType == -1:
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
        if isWireModel:
            C._initVars(t,'{centers:cellNFrontFilWMM}={centers:cellNFilWMM}*({centers:cellNFilWMM}>0.5)+1*({centers:cellNFilWMM}<0.5)')
            C._initVars(t,'{centers:cellNFrontFilWMM}=logical_and({centers:cellNFrontFilWMM}>0.5, {centers:cellNFrontFilWMM}<1.5)')

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

    C._initVars(t,'{centers:cellN}=maximum(0.,{centers:cellNChim})') # vaut -3, 0, 1, 2 initialement

    if cleanCellN:
        vars = ['centers:TurbulentDistanceAllBC','centers:TurbulentDistanceWallBC', 'centers:cellNIBC_hole']
        C._rmVars(t, vars)

    return None

#=========================================================================
# Build the IBM front for IBM pre-processing.
# IN: t (tree): computational tree
# IN: tc (tree): connectivity tree
# IN: dimPb (2 or 3): problem dimension
# IN: frontType (0,1,2 or 42): type of IBM front
# IN: cartesian (boolean): if True, activates optimized algorithms for Cartesian meshes
# IN: twoFronts (boolean): if True, performs the IBM pre-processing for an additional image point positioned farther away
# IN: check (boolean): if True, saves front.cgns and front2.cgns if twoFronts is active
# OUT: updated centers:cellNIBC and centers:cellNFront fields
# OUT: (optional) updated centers:cellNIBC_2, centers:cellNFront_2 fields for second image points
# OUT: front: front of image points
# OUT: front2: (optional): front of second image points
#=========================================================================
def _pushBackImageFront2__(t, tc, tbbc, cartesian=True):
    # bboxDict needed for optimised AddXZones (i.e. "layers" not None)
    # Return a dict with the zones of t as keys and their specific bboxes as key values
    bboxDict  = Cmpi.createBboxDict(t)
    interDict = X.getIntersectingDomains(tbbc, taabb=tbbc)
    graph     = Cmpi.computeGraph(tbbc, type='bbox', intersectionsDict=interDict, reduction=False)

    interpDataType = 0 if cartesian else 1

    # if subr, the tree subregions are kept during the exchange
    # if layers not None, only communicate the desired number of layers
    Cmpi._addLXZones(tc, graph, variables=['cellNIBC','cellNChim','cellNFront'], cartesian=cartesian,
                     interDict=interDict,bboxDict=bboxDict, layers=4, subr=False)
    Cmpi._addLXZones(t, graph, variables=['centers:cellNIBC', 'centers:cellNChim', 'centers:cellNFront'], cartesian=cartesian,
                     interDict=interDict, bboxDict=bboxDict, layers=4, subr=False)

    # Zones of tc are modified after addXZones, new tbbc, interDict and intersectionDict
    tbbcx             = G.BB(tc)
    #interDict         = X.getIntersectingDomains(tbbcx)
    intersectionsDict = X.getIntersectingDomains(tbbcx, method='AABB', taabb=tbbcx)

    # Reconstruction of cellNFront and cellN from cellNIBC (reduce the communications)
    # cellNFront_origin and cellNIBC_origin are initialised to store the Data of cellNFront and cellNIBC before the transfers
    C._initVars(t, '{centers:cellN}={centers:cellNIBC}')
    C._initVars(t, '{centers:cellNFront_origin}={centers:cellNFront}')
    C._initVars(t, '{centers:cellNIBC_origin}={centers:cellNIBC}')
    C._initVars(t, '{centers:cellN_interp}=maximum(0.,{centers:cellNChim})') # Second way of building the cellN field, see above

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
        cellNFront = Internal.getNodeFromName2(z, 'cellNFront2')
        if cellNFront != []:
            cellNFront = cellNFront[1]
            sizeTot = cellNFront.shape[0]*cellNFront.shape[1]*cellNFront.shape[2]
            sizeOne =  int(numpy.sum(cellNFront))
            if sizeOne < sizeTot:
                X._setHoleInterpolatedPoints(z, depth=1, dir=0, loc='centers', cellNName='cellNFront2', addGC=False)
                res = X.getInterpolatedPoints(z, loc='centers', cellNName='cellNFront2') # indices,X,Y,Z
                if res is not None:
                    indicesI = res[0]
                    XI = res[1]; YI = res[2]; ZI = res[3]
                    allInterpFields=[]
                    for zc in Internal.getZones(tc):
                        if zc[0] in intersectionsDict[z[0]]:
                            C._cpVars(zc,'cellN_interp',zc,'cellN')
                            fields = X.transferFields(zc, XI, YI, ZI, hook=None, variables=['cellNFront_origin','cellNIBC_origin'],
                                                      interpDataType=interpDataType, nature=1)
                            allInterpFields.append(fields)
                    if allInterpFields!=[]:
                        C._filterPartialFields(z, allInterpFields, indicesI, loc='centers', startFrom=0, filterName='donorVol',
                                               verbose=False)

    Cmpi._rmXZones(tc)
    Cmpi._rmXZones(t)

    # Update the cellNFront, cellNIBC and cellNIBCDnr fields
    for z in Internal.getZones(t):
        cellNFront = Internal.getNodeFromName2(z, 'cellNFront2')
        if cellNFront != []:
            cellNFront = cellNFront[1]
            sizeTot = cellNFront.shape[0]*cellNFront.shape[1]*cellNFront.shape[2]
            sizeOne =  int(numpy.sum(cellNFront))
            if sizeOne < sizeTot:
                # Modification du Front uniquement lorsque celui-ci est repousse
                C._initVars(z,'{centers:cellNFront}={centers:cellNFront}*({centers:cellNFront_origin}>0.5)')
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

    C._cpVars(t,'centers:cellNFront',tc,'cellNFront')

    Xmpi._setInterpTransfers(t,tc,variables=['cellNFront'], cellNVariable='cellNFront', compact=0)

    return None

def buildFrontIBM(t, tc, tb=None, dimPb=3, frontType=1, cartesian=True, twoFronts=False, check=False,
                  tbFilament=None):
    """Build the IBM front for IBM pre-processing."""

    isFilamentOnly, isWireModel = D_IBM.localWMMFlags__(tb, tbFilament)

    tbbc = Cmpi.createBBoxTree(tc)

    interpDataType = 0 if cartesian else 1

    C._initVars(t,'{centers:cellNIBCDnr}=minimum(2.,abs({centers:cellNIBC}))')
    C._initVars(t,'{centers:cellNIBC}=maximum(0.,{centers:cellNIBC})')# vaut -3, 0, 1, 2, 3 initialement
    C._initVars(t,'{centers:cellNIBC}={centers:cellNIBC}*({centers:cellNIBC}<2.5)')
    C._cpVars(t,'centers:cellNIBC',t,'centers:cellN')
    C._cpVars(t,'centers:cellN',tc,'cellN')

    # Transfert du cellNFront
    C._cpVars(t,'centers:cellNFront',tc,'cellNFront')

    # propager cellNVariable='cellNFront'
    Xmpi._setInterpTransfers(t, tc, variables=['cellNFront'], cellNVariable='cellNFront', compact=0)

    if twoFronts:
        C._cpVars(t,'centers:cellNFront_2',tc,'cellNFront_2')
        Xmpi._setInterpTransfers(t, tc, variables=['cellNFront_2'], cellNVariable='cellNFront_2', compact=0)

    if frontType == 2: _pushBackImageFront2__(t, tc, tbbc, cartesian=cartesian)

    C._rmVars(t,['centers:cellNFront'])
    if twoFronts: C._rmVars(t,['centers:cellNFront_2', 'centers:cellNIBC_2'])

    C._cpVars(t,'centers:TurbulentDistance',tc,'TurbulentDistance')

    print('Minimum distance: %f.'%C.getMinValue(t,'centers:TurbulentDistance'))
    P._computeGrad2(t, 'centers:TurbulentDistance', ghostCells=True, withCellN=False)

    front = getIBMFront(tc, 'cellNFront', dim=dimPb, frontType=frontType)
    front = gatherFront(front)

    if twoFronts:
        front2 = getIBMFront(tc, 'cellNFront_2', dim=dimPb, frontType=frontType)
        front2 = gatherFront(front2)
    else:
        front2 = None

    if isWireModel:
        frontWMM = getIBMFront(tc, 'cellNFrontFilWMM', dim=dimPb, frontType=frontType)
        frontWMM = gatherFront(frontWMM)
    else:
        frontWMM = None

    if check and Cmpi.rank == 0:
        C.convertPyTree2File(front, 'front.cgns')
        if twoFronts:
            C.convertPyTree2File(front2, 'front2.cgns')
        if isWireModel:
            C.convertPyTree2File(frontWMM, 'frontWMM.cgns')

    return t, tc, front, front2, frontWMM

#=========================================================================
# Compute the transfer coefficients and data for IBM pre-processing.
# IN: t (tree): computational tree
# IN: tc (tree): connectivity tree
# IN: tb (tree): geometry tree (IBM bodies)
# IN: front (tree): front of image points
# IN: front2 (tree, optional): front of second image points
# IN: dimPb (2 or 3): problem dimension
# IN: frontType (0,1,2 or 42): type of IBM front
# IN: IBCType (-1 or 1): type of IBM, -1: IB target points are located inside the solid, 1: IB target points are located in the fluid
# IN: depth (int): depth of overlapping regions
# IN: Reynolds (float): Reynolds number (F42)
# IN: yplus (float): estimated yplus at the first computed cells (F42)
# IN: Lref (float): characteristic length of the geometry (F42)
# IN: cartesian (boolean): if True, activates optimized algorithms for Cartesian meshes
# IN: twoFronts (boolean): if True, performs the IBM pre-processing for an additional image point positioned farther away
# OUT: IBCD* zones inside tc
# OUT: (optional) 2_IBCD* zones inside tc
#=========================================================================
def setInterpDataIBM(t, tc, tb, front, front2=None, dimPb=3, frontType=1, IBCType=1, depth=2, Reynolds=1.e6,
                     yplus=100, Lref=1., cartesian=True, twoFronts=False, check=False,
                     tbFilament=None, frontWMM=None):
    """Compute the transfer coefficients and data for IBM pre-processing."""
    tp = Internal.copyRef(t)
    _setInterpDataIBM(t, tc, tb, front, front2=front2, dimPb=dimPb, frontType=frontType, IBCType=IBCType, depth=depth,
                      Reynolds=Reynolds, yplus=yplus, Lref=Lref,
                      cartesian=cartesian, twoFronts=twoFronts, check=check,
                      tbFilament=tbFilament, frontWMM=frontWMM)
    return tp

def _setInterpDataIBM(t, tc, tb, front, front2=None, dimPb=3, frontType=1, IBCType=1, depth=2, Reynolds=1.e6,
                      yplus=100, Lref=1., cartesian=True, twoFronts=False, check=False,
                      tbFilament=None, frontWMM=None):
    """Compute the transfer coefficients and data for IBM pre-processing."""

    isFilamentOnly, isWireModel = D_IBM.localWMMFlags__(tb, tbFilament)

    for zc in Internal.getZones(tc):
        proc = Cmpi.getProc(zc)
        if proc == -1: Cmpi._setProc(zc, 0)

    tbbc = Cmpi.createBBoxTree(tc)

    interpDataType = 0 if cartesian else 1

    zonesRIBC = []
    for zrcv in Internal.getZones(t):
        if C.getMaxValue(zrcv, 'centers:cellNIBC')==2.:
            zonesRIBC.append(zrcv)

    nbZonesIBC = len(zonesRIBC)
    if nbZonesIBC == 0:
        res = [{},{},{}]
        if twoFronts or isWireModel: res2 = [{},{},{}]
    else:
        tb_local = tb
        if tbFilament:
            if isFilamentOnly:tb_local= tbFilament
            else:             tb_local= Internal.merge([tb,tbFilament])
        res = getAllIBMPoints(zonesRIBC, loc='centers',tb=tb_local, tfront=front, frontType=frontType,
                              cellNName='cellNIBC', depth=depth, IBCType=IBCType, Reynolds=Reynolds, yplus=yplus, Lref=Lref,
                              isOrthoFirst=isFilamentOnly, check=check)
        if twoFronts:
            res2 = getAllIBMPoints(zonesRIBC, loc='centers',tb=tb, tfront=front2, frontType=frontType,
                                   cellNName='cellNIBC', depth=depth, IBCType=IBCType, Reynolds=Reynolds, yplus=yplus, Lref=Lref, check=check,
                                   twoFronts=twoFronts)
        if isWireModel:
            tb_filament_localWMM=Internal.copyTree(tbFilament)
            for z in Internal.getZones(tbFilament):
                ibctype = Internal.getNodeFromName(Internal.getNodeFromName(z,'.Solver#define'),'ibctype')
                if Internal.getValue(ibctype) != "wiremodel":
                    zlocal=Internal.getNodeFromName(tb_filament_localWMM,z[0])
                    Internal._rmNode(tb_filament_localWMM,zlocal)
            res2 = getAllIBMPoints(zonesRIBC, loc='centers',tb=tb_filament_localWMM, tfront=frontWMM, frontType=frontType,
                                   cellNName='cellNFilWMM', depth=depth, IBCType=IBCType, Reynolds=Reynolds, yplus=yplus, Lref=Lref,
                                   isWireModel=isWireModel, isOrthoFirst=isFilamentOnly, check=check)

            restmp = getAllIBMPoints(zonesRIBC, loc='centers',tb=tb_filament_localWMM, tfront=frontWMM, frontType=frontType,
                                     cellNName='cellNFilWMM', depth=depth, IBCType=IBCType, Reynolds=Reynolds,
                                     yplus=yplus, Lref=Lref, isOrthoFirst=isFilamentOnly, check=check)

            for j in range(3):
                ##delete in res
                item_del=[]
                for ii in res[j]:
                    if "140" in ii: item_del.append(ii)
                for ii in item_del: del res[j][ii]
                ##detele in res2
                item_del=[]
                for ii in res2[j]:
                    if "140" not in ii: item_del.append(ii)
                for ii in item_del: del res2[j][ii]
                ##add to res
                for ii in restmp[j]: res[j][ii] = restmp[j][ii]

    # cleaning
    C._rmVars(tc,['cellNChim','cellNIBC','TurbulentDistance','cellNFront'])
    # dans t, il faut cellNChim et cellNIBCDnr pour recalculer le cellN a la fin
    varsRM = ['centers:gradxTurbulentDistance','centers:gradyTurbulentDistance','centers:gradzTurbulentDistance','centers:cellNFront','centers:cellNIBC']
    C._rmVars(t, varsRM)
    front = None
    if twoFronts: front2 = None

    if isWireModel:
        C._rmVars(tc,['cellNFrontFilWMM'])
        C._rmVars(t,['cellNFrontFilWMM'])
        frontWMM = None

    graph = {}; datas = {}
    procDict = Cmpi.getProcDict(tc)

    # graph d'intersection des pts images de ce proc et des zones de tbbc
    zones  = Internal.getZones(tbbc)
    allBBs = []
    dictOfCorrectedPtsByIBCType = res[0]
    dictOfWallPtsByIBCType      = res[1]
    dictOfInterpPtsByIBCType    = res[2]
    interDictIBM={}

    if twoFronts or isWireModel:
        dictOfCorrectedPtsByIBCType2 = res2[0]
        dictOfWallPtsByIBCType2      = res2[1]
        dictOfInterpPtsByIBCType2    = res2[2]
    else:
        dictOfCorrectedPtsByIBCType2={}
        dictOfWallPtsByIBCType2     ={}
        dictOfInterpPtsByIBCType2   ={}
    interDictIBM2={}

    if dictOfCorrectedPtsByIBCType!={}:
        for ibcTypeL in dictOfCorrectedPtsByIBCType:
            allCorrectedPts = dictOfCorrectedPtsByIBCType[ibcTypeL]
            allWallPts      = dictOfWallPtsByIBCType[ibcTypeL]
            allInterpPts    = dictOfInterpPtsByIBCType[ibcTypeL]
            for nozr in range(nbZonesIBC):
                if allCorrectedPts[nozr] != []:
                    ##[AJ] Keep temporarily
                    ##if ibcTypeL=="140#filament":
                    ##    nlen  = numpy.shape(allInterpPts[nozr][1])[1]
                    ##    save2file = numpy.zeros((nlen,2),dtype=float)
                    ##    save2file[:,0]=allInterpPts[nozr][1][0][:]
                    ##    save2file[:,1]=allInterpPts[nozr][1][1][:]
                    ##    numpy.savetxt('allInterpPts.txt', save2file, delimiter=',')   # X is an array
                    zrname = zonesRIBC[nozr][0]
                    interpPtsBB = Generator.BB(allInterpPts[nozr])
                    for z in zones:
                        bba = C.getFields('GridCoordinates', z)[0]
                        if Generator.bboxIntersection(interpPtsBB, bba, isBB=True):
                            zname = z[0]
                            popp  = Cmpi.getProc(z)
                            if Cmpi.size > 1: Distributed.updateGraph__(graph, popp, Cmpi.rank, zname)

                            if zrname not in interDictIBM: interDictIBM[zrname]=[zname]
                            else:
                                if zname not in interDictIBM[zrname]: interDictIBM[zrname].append(zname)
        if twoFronts or isWireModel:
            for ibcTypeL in dictOfCorrectedPtsByIBCType2:
                allCorrectedPts2 = dictOfCorrectedPtsByIBCType2[ibcTypeL]
                allWallPts2      = dictOfWallPtsByIBCType2[ibcTypeL]
                allInterpPts2    = dictOfInterpPtsByIBCType2[ibcTypeL]
                for nozr in range(nbZonesIBC):
                    if allCorrectedPts2[nozr] != []:
                        ##[AJ] Keep temporarily
                        ##if ibcTypeL=="140#filament":
                        ##    nlen  = numpy.shape(allInterpPts2[nozr][1])[1]
                        ##    save2file = numpy.zeros((nlen,2),dtype=float)
                        ##    save2file[:,0]=allInterpPts2[nozr][1][0][:]
                        ##    save2file[:,1]=allInterpPts2[nozr][1][1][:]
                        ##    numpy.savetxt('allInterpPts2.txt', save2file, delimiter=',')   # X is an array
                        zrname = zonesRIBC[nozr][0]
                        interpPtsBB2 = Generator.BB(allInterpPts2[nozr])
                        for z in zones:
                            bba = C.getFields('GridCoordinates', z)[0]
                            if Generator.bboxIntersection(interpPtsBB2,bba,isBB=True):
                                zname = z[0]
                                popp  = Cmpi.getProc(z)
                                if Cmpi.size > 1: Distributed.updateGraph__(graph, popp, Cmpi.rank, zname)
                                if zrname not in interDictIBM2: interDictIBM2[zrname]=[zname]
                                else:
                                    if zname not in interDictIBM2[zrname]: interDictIBM2[zrname].append(zname)
    else: graph={}
    if Cmpi.KCOMM is not None: allGraph = Cmpi.KCOMM.allgather(graph)
    else: allGraph = [graph]

    graph = {}
    for i in allGraph:
        for k in i:
            if not k in graph: graph[k] = {}
            for j in i[k]:
                if not j in graph[k]: graph[k][j] = []
                graph[k][j] += i[k][j]
                graph[k][j] = list(set(graph[k][j])) # pas utile?

    # keyword subr=False to avoid memory overflow
    Cmpi._addXZones(tc, graph, variables=['cellN'], cartesian=cartesian, subr=False)

    ReferenceState = Internal.getNodeFromType2(t, 'ReferenceState_t')
    model = Internal.getNodeFromName(t, 'GoverningEquations')
    if model is not None: model = Internal.getValue(model)
    else:                 model = "Euler"

    for i in range(Cmpi.size): datas[i] = [] # force

    if dictOfCorrectedPtsByIBCType!={}:
        for ibcTypeL in dictOfCorrectedPtsByIBCType:
            if '#' in ibcTypeL: ibcNameL = 'IBCD_'+'_'.join(ibcTypeL.split('#')) #ibctype with familyname
            else: ibcNameL = 'IBCD_'+ibcTypeL #regular ibctype
            allCorrectedPts = dictOfCorrectedPtsByIBCType[ibcTypeL]
            allWallPts = dictOfWallPtsByIBCType[ibcTypeL]
            allInterpPts = dictOfInterpPtsByIBCType[ibcTypeL]
            for nozr in range(nbZonesIBC):
                if allCorrectedPts[nozr] != []:
                    zrcv = zonesRIBC[nozr]
                    zrname = zrcv[0]
                    dnrZones = []
                    for zdname in interDictIBM[zrname]:
                        zd = Internal.copyRef(Internal.getNodeFromName2(tc, zdname))
                        if zd is None: print('!!!Zone None', zrname, zdname)
                        else: dnrZones.append(zd)
                    XOD._setIBCDataForZone__(zrcv, dnrZones, allCorrectedPts[nozr], allWallPts[nozr], allInterpPts[nozr],
                                             nature=1, penalty=1, loc='centers', storage='inverse', dim=dimPb,
                                             interpDataType=interpDataType, ReferenceState=ReferenceState, bcType=ibcTypeL,model=model)

                    nozr += 1
                    for zd in dnrZones:
                        zdname = zd[0]
                        destProc = procDict[zdname]

                        IDs = []
                        for i in zd[2]:
                            if i[0][0:4] == 'IBCD':
                                #where i[0] == 'IBCD_#ibctype_#familyname_#zname' or 'IBCD_#ibctype_#zname'
                                if ibcNameL == '_'.join(i[0].split('_')[:-1]) and Internal.getValue(i)==zrname:
                                    IDs.append(i) #add only ibcd zones related to local ibcType to avoid doublons

                        if IDs != []:
                            if destProc == Cmpi.rank:
                                zD = Internal.getNodeFromName2(tc, zdname)
                                zD[2] += IDs
                            else:
                                if destProc not in datas: datas[destProc]=[[zdname,IDs]]
                                else: datas[destProc].append([zdname,IDs])
                        else:
                            if destProc not in datas: datas[destProc] = []

    if dictOfCorrectedPtsByIBCType2!={}:
        for ibcTypeL in dictOfCorrectedPtsByIBCType2:
            if '#' in ibcTypeL: ibcNameL = '2_IBCD_'+'_'.join(ibcTypeL.split('#'))
            else: ibcNameL = '2_IBCD_'+ibcTypeL
            allCorrectedPts2 = dictOfCorrectedPtsByIBCType2[ibcTypeL]
            allWallPts2      = dictOfWallPtsByIBCType2[ibcTypeL]
            allInterpPts2    = dictOfInterpPtsByIBCType2[ibcTypeL]
            for nozr in range(nbZonesIBC):
                if allCorrectedPts2[nozr] != []:
                    zrcv     = zonesRIBC[nozr]
                    zrname   = zrcv[0]
                    dnrZones = []
                    for zdname in interDictIBM2[zrname]:
                        zd = Internal.copyRef(Internal.getNodeFromName2(tc, zdname))
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
                                if ibcNameL == '_'.join(i[0].split('_')[:-1]) and Internal.getValue(i)==zrname:
                                    IDs.append(i) #add only ibcd zones related to local ibcType to avoid doublons

                        if IDs != []:
                            if destProc == Cmpi.rank:
                                zD = Internal.getNodeFromName2(tc, zdname)
                                zD[2] += IDs
                            else:
                                if destProc not in datas: datas[destProc]=[[zdname,IDs]]
                                else: datas[destProc].append([zdname,IDs])
                        else:
                            if destProc not in datas: datas[destProc] = []

    Cmpi._rmXZones(tc)
    dictOfCorrectedPtsByIBCType = None
    dictOfWallPtsByIBCType      = None
    dictOfInterpPtsByIBCType    = None
    interDictIBM = None
    if twoFronts or isWireModel:
        dictOfCorrectedPtsByIBCType2 = None
        dictOfWallPtsByIBCType2      = None
        dictOfInterpPtsByIBCType2    = None
        interDictIBM2 = None

    Internal._rmNodesByName(tc, Internal.__FlowSolutionNodes__)
    destDatas = Cmpi.sendRecv(datas, graph)
    for i in destDatas:
        for n in destDatas[i]:
            zname = n[0]
            IBCDs = n[1]
            if IBCDs != []:
                zD = Internal.getNodeFromName2(tc, zname)
                zD[2] += IBCDs

    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is not None: model = Internal.getValue(model)
    else: model = "NSTurbulent"

    C._initVars(t,'{centers:cellN}=minimum({centers:cellNChim}*{centers:cellNIBCDnr},2.)')
    varsRM = ['centers:cellNChim', 'centers:cellNIBCDnr']
    if model == 'Euler': varsRM += ['centers:TurbulentDistance']
    C._rmVars(t, varsRM)

    if check:
        extractIBMInfo(tc, IBCNames="IBCD_*", fileout='IBMInfo.cgns')

        if twoFronts:
            extractIBMInfo(tc, IBCNames="2_IBCD_*", fileout='IBMInfo2.cgns')

    return None

#=========================================================================
# Init the final trees for IBM pre-processing.
# IN: t (tree): computational tree
# IN: tc (tree): connectivity tree
# IN: tb (tree): geometry tree (IBM bodies)
# IN: tinit (tree): computational tree from previous computation used to initialize the flow solution in t
# IN: tbCurvi (tree): geometry tree (curvilign bodies) for hybrid computation cart + curvi
# IN: dimPb (2 or 3): problem dimension
# IN: twoFronts (boolean): if True, creates a new connectivity tree that contains second image points information
# OUT: updated t, tc
# OUT: (optional) new connectivity tree tc2
#=========================================================================
def _recomputeDistForViscousWall__(t, tb, tbCurvi=None, dimPb=3, tbFilament=None):

    isFilamentOnly, isWireModel = D_IBM.localWMMFlags__(tb, tbFilament)

    filteredBC=['outpress', 'inj', 'slip', 'overlap']
    recompute=False
    tb2 = Internal.copyRef(tb)
    for z in Internal.getZones(tb2):
        ibc = Internal.getValue( Internal.getNodeFromName(z, 'ibctype') )
        if ibc in filteredBC:
            Internal._rmNode(tb2, z)
            recompute=True
    if tbCurvi is not None:
        recompute=True # on recalcule distance si cas hybride cart + curvi
        bases = Internal.getBases(tbCurvi)
        for base in bases: base[0]='curvi'+base[0]
        tb2[2] += bases

    if recompute:
        if dimPb == 2:
            z0 = Internal.getNodeFromType2(t, "Zone_t")
            bb0 = G.bbox(z0); dz = (bb0[5]-bb0[2])*0.5
            tb2 = Internal.copyValue(tb2, byName='CoordinateZ')
            C._initVars(tb2, 'CoordinateZ', dz)
    else: return None

    tbsave = tb2

    if tbFilament and not isFilamentOnly:
        if dimPb == 2: tb2 = C.initVars(tbFilament, 'CoordinateZ', dz)
        tbFilamentnoWMM = []
        tbFilamentWMM   = []
        for z in Internal.getZones(tb2):
            soldef  = Internal.getNodeFromName(z,'.Solver#define')
            ibctype = Internal.getNodeFromName(soldef,'ibctype')
            if Internal.getValue(ibctype) != "wiremodel":
                tbFilamentnoWMM.append(z)
            else:
                tbFilamentWMM.append(z)

        if tbFilamentnoWMM:
            tbFilamentnoWMM = C.newPyTree(['BasenoWMM', tbFilamentnoWMM])
            tbsave          = Internal.merge([tbsave,tbFilamentnoWMM])

    DTW._distance2Walls(t, tbsave, type='ortho', signed=0, dim=dimPb, loc='centers')

    return None

def _tcInitialize__(tc, tc2=None, ibctypes=[], isWireModel=False):

    if tc2:
        Internal._rmNodesByName(tc2, 'IBCD*')
        Internal._rmNodesByName(tc, '2_IBCD*')

        D_IBM._transformTc2(tc2)

    if RENAMEIBCNODES:
        for zc in Internal.getZones(tc):
            for ibcd in Internal.getNodesFromName1(zc,'IBCD_*'):
                proposedName = Internal.getName(ibcd)[0:6]+'_X%d'%(Cmpi.rank)
                ibcd[0] = getIBCDName(proposedName)

        if tc2:
            for zc in Internal.getZones(tc2):
                for ibcd in Internal.getNodesFromName1(zc,'IBCD_*'):
                    proposedName = Internal.getName(ibcd)[0:8]+'_X%d'%(Cmpi.rank)
                    ibcd[0] = getIBCDName(proposedName)

    # Adding a userdefined node to the tc tree for the IBC conditions that are provided
    # to FastS solver to reduce the number of input arguments and to make a clear distinction
    # of the solver parameters and the IBC parameters
    base = Internal.getBases(tc)[0]
    Internal._createUniqueChild(base, '.Solver#IBCdefine', 'UserDefinedData_t')
    solverIBC = Internal.getNodeFromName1(base, '.Solver#IBCdefine')

    Internal._createUniqueChild(solverIBC, 'Reref', 'DataArray_t', -1)
    Internal._createUniqueChild(solverIBC, 'Lref' , 'DataArray_t',  1)

    Internal._createUniqueChild(solverIBC, 'isgradP'    , 'DataArray_t', 'False')
    Internal._createUniqueChild(solverIBC, 'isWireModel', 'DataArray_t', 'False')
    Internal._createUniqueChild(solverIBC, 'isTBLE'     , 'DataArray_t', 'False')

    # note: here alphagrad is the corrected nomenclature for alghagradp found in param_solver.h (B.Constant confirmed)
    #      changed some other variables names to be consistent with other options/coding "guidelines"
    if 'Mafzal' in ibctypes:
        Internal._createUniqueChild(solverIBC, 'isgradP'   , 'DataArray_t', 'True')
        Internal._createUniqueChild(solverIBC, 'mafzalMode', 'DataArray_t', 0)
        Internal._createUniqueChild(solverIBC, 'alphaGrad' , 'DataArray_t', 0)
    # if twoFronts:
    #     Internal._createUniqueChild(solverIBC, 'isgradP'    , 'DataArray_t', 'True')
    #     Internal._createUniqueChild(solverIBC, 'alphaGrad'  , 'DataArray_t', 0)
    if isWireModel:
        Internal._createUniqueChild(solverIBC, 'isWireModel' , 'DataArray_t', 'True')
        Internal._createUniqueChild(solverIBC, 'DeltaVWire'  , 'DataArray_t', 0)
        Internal._createUniqueChild(solverIBC, 'KWire'       , 'DataArray_t', 0)
        Internal._createUniqueChild(solverIBC, 'DiameterWire', 'DataArray_t', 0)
        Internal._createUniqueChild(solverIBC, 'CtWire'      , 'DataArray_t', 0)
    if 'TBLE' in ibctypes or 'TBLE_FULL' in ibctypes:
        Internal._createUniqueChild(solverIBC, 'isTBLE'        , 'DataArray_t', 'True')
        Internal._createUniqueChild(solverIBC, 'alphaGrad'     , 'DataArray_t', 0)
        Internal._createUniqueChild(solverIBC, 'NbPtsLinelits' , 'DataArray_t', 0)

    return None

def _tInitialize__(t, tinit=None, model='NSTurbulent', isWireModel=False):
    if tinit is None: I._initConst(t, loc='centers')
    else: t = Pmpi.extractMesh(tinit, t, mode='accurate')

    if model != "Euler": C._initVars(t, 'centers:ViscosityEddy', 0.)
    if isWireModel:
        vars_wm = ['Density','VelocityX','VelocityY','VelocityZ','Temperature']
        if model == 'NSTurbulent':vars_wm.append('TurbulentSANuTilde')
        for z in Internal.getZones(t):
            for v_local in vars_wm:
                C._initVars(z,'{centers:'+v_local+'_WM}=0.')
    return None

def initializeIBM(t, tc, tb, tinit=None, tbCurvi=None, dimPb=3, twoFronts=False, tbFilament=None, cleanCellN=True):
    """Initialize the computational and connectivity trees for IBM pre-processing."""

    isFilamentOnly, isWireModel = D_IBM.localWMMFlags__(tb, tbFilament)

    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('initializeIBM: GoverningEquations is missing in input geometry tree.')
    model = Internal.getValue(model)

    ibctypes = Internal.getNodesFromName(tb, 'ibctype')
    if ibctypes is None: raise ValueError('initializeIBM: ibctype is missing in input geometry tree.')
    ibctypes = list(set(Internal.getValue(ibc) for ibc in ibctypes))

    if model != 'Euler' and t is not None:
        _recomputeDistForViscousWall__(t, tb, tbCurvi=tbCurvi, dimPb=dimPb, tbFilament=tbFilament)
        # if not cleanCellN: C._initVars(t, '{centers:TurbulentDistanceWallBC}={centers:TurbulentDistance}')

    tc2 = Internal.copyTree(tc) if twoFronts or isWireModel else None

    if isWireModel:
        Internal._rmNodesByName(tc2, 'IBCD*')
        Internal._rmNodesByName(tc, '2_IBCD*')
        D_IBM._transformTc2(tc2)
        tc2 = Internal.rmNodesByName(tc2, 'ID*')
        NewIBCD= 141
        _changeNameIBCD__(tc2,NewIBCD)
        tc  = Internal.merge([tc,tc2])
        tc2 = None

    _tcInitialize__(tc, tc2=tc2, ibctypes=ibctypes, isWireModel=isWireModel)
    if t is not None: _tInitialize__(t, tinit=tinit, model=model, isWireModel=isWireModel)

    return t, tc, tc2

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## BLANKING
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#==============================================================================
#
#==============================================================================
def _blankClosestTargetCells(t, cellNName='cellN', depth=3):
    """Blank target points."""
    for z in Internal.getZones(t):
        connector._blankClosestTargetCells(z, depth, cellNName,
                                           Internal.__GridCoordinates__,
                                           Internal.__FlowSolutionNodes__,
                                           Internal.__FlowSolutionCenters__)
    return None

#==============================================================================
# Remove fully blanked grids considering cellNNIBC and cellNChim
#==============================================================================
def _removeBlankedGrids(t, loc='centers'):
    """Delete zones when fully blanked considering cellNNIBC and cellNChim."""
    vari = 'cellNIBC'
    varc = 'cellNChim'
    flag = 'flag'
    if loc == 'centers':
        vari = 'centers:'+vari
        varc = 'centers:'+varc
        flag = 'centers:'+flag

    C._initVars(t,'{%s}=abs(1.-{%s}*{%s})<0.5'%(flag,vari,varc))

    #poschim = C.isNamePresent(t,loc+varc)
    #posibc = C.isNamePresent(t,loc+vari)
    #if poschim != -1 and posibc != -1: C._initVars(t,'{%s}=abs(1.-{%s}*{%s})<0.5'%(flag,vari,varc))
    #elif poschim != -1 and posibc==-1: flag=varc
    #elif poschim == -1 and posibc!=-1: flag=vari
    #else: return None
    for z in Internal.getZones(t):
        if C.getMaxValue(z,flag) < 0.5:
            (parent,noz) = Internal.getParentOfNode(t, z)
            del parent[2][noz]
        else:
            #if poschim != -1 and posibc != -1: C._rmVars(z,[flag])
            C._rmVars(z,[flag])
    return None

#==============================================================================
# masquage par les corps IBC
# gridType = single or composite - composite means that an off body grid exists
# IN: tb: immersed bodies
# IN: loc: "centers" or "nodes"
# IN: dim: 2 or 3
#==============================================================================
def blankByIBCBodies(t, tb, loc, dim, cellNName='cellN'):
    """Blank by immersed bodies."""
    tp = Internal.copyRef(t)
    _blankByIBCBodies(tp, tb, loc, dim, cellNName=cellNName)

    return tp

def _blankByIBCBodies(t, tb, loc, dim, cellNName='cellN'):
    """Blank by immersed bodies."""
    DIM = dim
    blankalgo='tri'
    #blankalgo='xray'
    if DIM == 2: blankalgo = 'xray'

    bodies = []
    for b in Internal.getBases(tb):
        wallsl = Internal.getNodesFromType1(b, 'Zone_t')
        #soldef = Internal.getNodeFromName(wallsl,'.Solver#define')
        bodies.append(wallsl)
        # if wallsl != []:
        #     try:
        #         wallsl = C.convertArray2Tetra(wallsl)
        #         wallsl = T.join(wallsl)
        #         wallsl = G.close(wallsl)
        #         Internal.addChild(wallsl,soldef)
        #         bodies.append([wallsl])
        #         # if DIM == 3:
        #         #     try: P.exteriorFaces(wallsl)
        #         #     except: pass
        #         #     bodies.append([wallsl])
        #         # else: bodies.append([wallsl])
        #     except:
        #         wallsl = C.convertArray2Tetra(wallsl)
        #         bodies.append(wallsl)

    nbodies = len(bodies)
    if nbodies == 0:
        print("Warning: blankByIBCBodies: no body defined.")
        return t

    print('Info: blankByIBCBodies: blanking mesh by %d bodies'%nbodies)
    if loc == 'centers': typeb = 'center_in'
    else: typeb = 'node_in'
    nbases = len(Internal.getBases(t))

    bodiesInv=[]
    for body in bodies:
        inv = Internal.getNodeFromName(body,'inv')
        if inv is not None: inv = Internal.getValue(inv)
        else: inv = 0
        if inv == 1:
            bodies.remove(body)
            bodiesInv.append(body)

    if blankalgo == 'xray' or DIM == 2:
        BM = numpy.ones((nbases,nbodies),dtype=Internal.E_NpyInt)
        dh_min = G_IBM_Height.getMinimumCartesianSpacing(t)
        XRAYDIM1 = 2000; XRAYDIM2 = XRAYDIM1
        if dh_min > 0.:
            bb = G.bbox(tb)
            Lxref = bb[3]-bb[0]
            Lyref = bb[4]-bb[1]
            XRAYDIM1 = max(XRAYDIM1,int(Lxref/(0.15*dh_min)))
            XRAYDIM2 = max(XRAYDIM2,int(Lyref/(0.15*dh_min)))
        if DIM == 2: XRAYDIM2 = 2

        if loc == 'centers':
            tc = C.node2Center(t)
            for body in bodiesInv:
                print('Info: blankByIBCBodies: reverse blanking for body.')
                tc = X.blankCells(tc, [body], BM, blankingType='node_in', XRaydim1=XRAYDIM1, XRaydim2=XRAYDIM2, dim=DIM, cellNName=cellNName)
                C._initVars(tc,'{%s}=1.-{%s}'%(cellNName,cellNName)) # ecoulement interne

            for body in bodies:
                tc = X.blankCells(tc, [body], BM, blankingType='node_in', XRaydim1=XRAYDIM1, XRaydim2=XRAYDIM2, dim=DIM, cellNName=cellNName)

            C._cpVars(tc,'%s'%cellNName,t,'centers:%s'%cellNName)
        else:
            X._blankCells(t, bodies, BM, blankingType=typeb, delta=TOLDIST, XRaydim1=XRAYDIM1, XRaydim2=XRAYDIM2, dim=DIM, cellNName=cellNName)
    else:
        BM2 = numpy.ones((nbases,1),dtype=Internal.E_NpyInt)
        for body in bodiesInv:
            print('Info: blankByIBCBodies: reverse blanking for body.')
            X._blankCellsTri(t, [body], BM2, blankingType=typeb, cellNName=cellNName)
            if loc == 'centers': C._initVars(t,'{centers:%s}=1-{centers:%s}'%(cellNName,cellNName)) # ecoulement interne
            else:                C._initVars(t,'{%s}=1-{%s}'%(cellNName,cellNName)) # ecoulement interne
        for body in bodies:
            X._blankCellsTri(t, [body], BM2, blankingType=typeb, cellNName=cellNName)

    return None

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## FRONT
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#=============================================================================
# Returns the front defining the image points
# IN: frontvar:
# IN: dim: 2 or 3
# IN: frontType: 0 (equal distance), 1 (minimum distance), 42 (modelling height)
#=============================================================================
def getIBMFrontType1__(tc, frontvar, dim):
    if dim == 2:
        z0 = Internal.getNodeFromType2(tc, 'Zone_t')
        if z0 is not None:
            zmean = C.getValue(z0, 'CoordinateZ', 0)
            dz = 2*zmean
        else: zmean=0.5; dz = 1.
    else: dz = 0.
    front = []
    for z in Internal.getZones(tc):
        if C.getMinValue(z,frontvar)<0.2 and C.getMaxValue(z,frontvar)>0.8:
            X._maximizeBlankedCells(z, depth=1, dir=1, loc='nodes', cellNName='cellNChim')
            C._initVars(z,'{cellNChim}=minimum(1.,{cellNChim})')
            f = P.frontFaces(z, frontvar)
            if Internal.getZoneDim(f)[1]>0:
                Internal._rmNodesByName(f,'ID_*')
                front.append(f)
    C._initVars(front,'{tag}=({cellNChim}>0.5)*({cellNChim}<1.5)')
    front = P.selectCells2(front, 'tag', strict=1)
    Internal._rmNodesByName(front,Internal.__FlowSolutionNodes__)
    Internal._rmNodesByName(front,Internal.__FlowSolutionCenters__)
    if dim == 2:
        T._addkplane(front)
        T._translate(front,(0,0,-zmean))
        T._contract(front, (0,0,0), (1,0,0), (0,1,0), 0.9*dz)
    return front

def getIBMFrontType0__(tc, frontvar, dim, isFront2=False, frontType=0, SHIFTB=0.):
    import Converter.Mpi as Cmpi

    if dim == 2:
        z0 = Internal.getNodeFromType2(tc, 'Zone_t')
        zmean = C.getValue(z0, 'CoordinateZ', 0)
        dz = 2*zmean
    else: dz = 0.

    SHIFTD = 1.+SHIFTF
    # SHIFTD *= SHIFTB
    front = []
    tf = Internal.copyRef(tc)
    C._initVars(tf,'{%s}={%s}-2.*({%s}>1.5)'%(frontvar,frontvar,frontvar))
    if isFront2:
        if dim == 2:
            front = []
            # Creation du corps 2D pour le preprocessing IBC
            for z in Internal.getZones(tf):
                epsilon_dist = abs(C.getValue(z,'CoordinateX',1)-C.getValue(z,'CoordinateX',0))
                dmin = math.sqrt(2)*5*epsilon_dist
                if frontType == 42:
                    dmin = max(dmin, SHIFTB+math.sqrt(2)*3*epsilon_dist) # where shiftb = hmod
                tcl = T.addkplane(z)
                T._translate(tcl,(0,0,-zmean))
                T._contract(tcl, (0,0,0), (1,0,0), (0,1,0), dz)
                front.append(P.isoSurfMC(tcl,'TurbulentDistance',dmin*SHIFTD))
                del tcl
            front = C.newPyTree(['Base']+front)
        else:
            front = []
            # Creation du corps 2D pour le preprocessing IBC
            for z in Internal.getZones(tf):
                epsilon_dist = abs(C.getValue(z,'CoordinateX',1)-C.getValue(z,'CoordinateX',0))
                dmin = math.sqrt(3)*5*epsilon_dist
                if frontType == 42:
                    dmin = max(dmin, SHIFTB+math.sqrt(3)*3*epsilon_dist) # where shiftb = hmod
                front.append(P.isoSurfMC(z,'TurbulentDistance',dmin*SHIFTD))
            front = C.newPyTree(['Base']+front)

    else:
        for z in Internal.getZones(tf):
            if C.getMinValue(z,frontvar)==0. and C.getMaxValue(z,frontvar)==1.:
                f = P.frontFaces(z, frontvar)
                if Internal.getZoneDim(f)[1]>0:
                    Internal._rmNodesByName(f,'ID_*')
                    front.append(f)
        if dim == 2:
            dmin = C.getMaxValue(front, 'TurbulentDistance')

            if Cmpi.KCOMM is not None:
                dmin = numpy.array([dmin], dtype=numpy.float64)
                dmin_max = numpy.zeros(1, dtype=numpy.float64)
                Cmpi.KCOMM.Allreduce(dmin, dmin_max, Cmpi.MAX)
                dmin = dmin_max[0]

            tcl = T.addkplane(tc)
            T._translate(tcl,(0,0,-zmean))
            T._contract(tcl, (0,0,0), (1,0,0), (0,1,0), dz)
            front = P.isoSurfMC(tcl,'TurbulentDistance',dmin*SHIFTD)
            del tcl
        else:
            dmin = C.getMaxValue(front, 'TurbulentDistance')

            if Cmpi.KCOMM is not None:
                dmin = numpy.array([dmin], dtype=numpy.float64)
                dmin_max = numpy.zeros(1, dtype=numpy.float64)
                Cmpi.KCOMM.Allreduce(dmin, dmin_max, Cmpi.MAX)
                dmin = dmin_max[0]

            front = P.isoSurfMC(tc, 'TurbulentDistance', dmin*SHIFTD)
    return front

def getIBMFront(tc, frontvar, dim, frontType, isFront2=False, SHIFTB=0.):
    """Returns the front used to define the location of the IBM image points."""

    # if frontType == 1 or frontType == 2 : front = getIBMFrontType1__(tc,frontvar,dim)
    if frontType == 0:
        front = getIBMFrontType0__(tc,frontvar,dim,isFront2,frontType,SHIFTB)
    else:
        front = getIBMFrontType1__(tc, frontvar, dim)
    front = C.deleteEmptyZones(front)
    Internal._rmNodesFromName(front, "ID_*")
    if frontType == 0: return front

    dxmin = 1.e12
    if frontType>0:
        front = Internal.getZones(front)
        dxmax = 1.e-12
        dht = [[]]*len(front)
        nof = 0
        for f in front:
            subf = T.subzone(f, [0], type='elements')
            dx = C.getMaxValue(subf,"CoordinateX")-C.getMinValue(subf,"CoordinateX")
            dy = C.getMaxValue(subf,"CoordinateY")-C.getMinValue(subf,"CoordinateY")
            if dim == 3:
                dz = C.getMaxValue(subf,"CoordinateZ")-C.getMinValue(subf,"CoordinateZ")
                dht[nof] = max(dx,dy,dz)
            else:
                dht[nof] = max(dx,dy)

            if dht[nof] < dxmin and dht[nof] > 1.e-12:
                dxmin = dht[nof]
            if dht[nof] > dxmax:
                dxmax = dht[nof]
            nof += 1

        nlevels = int(math.log(dxmax/dxmin)/math.log(2)+1)

        dictOfLevels={}
        for l in range(nlevels): dictOfLevels[l]=[]
        for nof in range(len(dht)):
            if dht[nof]>1.e-12:
                nolvl = int(math.log(dht[nof]/dxmin)/math.log(2))
                dictOfLevels[nolvl].append(front[nof])

        front=[]
        for nol in range(nlevels):
            if len(dictOfLevels[nol])>0:
                front.append(T.join(dictOfLevels[nol]))

    return front

#=============================================================================
# Gather front
# Si un front est calcule par morceau sur chaque proc, ramene le meme front
# sur tous les procs
#=============================================================================
def gatherFront(front):
    """Gather front for mpi."""
    import Converter.Mpi as Cmpi
    zones = Internal.getNodesFromType1(front, 'Zone_t')
    for z in zones: z[0] += '_'+str(Cmpi.rank)

    if Cmpi.KCOMM is not None:
        #allFront = Cmpi.KCOMM.allgather(front)
        #front = []
        #for f in allFront: front += f
        front = Cmpi.allgatherZones(front)
        return front
    else: return front

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## IBM INFO
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#=============================================================================
# For FastLBM with IBMs. As FastLBM/IBC will evolve making this function decrepit
# no additional information will be currently added [May 12, 2023].
# Get total number of IBM points & puts in param_int. Need to FastLBM/IBC that
# applies the IBMs on the distribution functions.
#=============================================================================
def _extractIBMInfo_param(t, tc):
    """Extract IBM informations from tc."""
    XPC   ={}
    Zones = []
    for z in Internal.getZones(tc):
        allIBCD = Internal.getNodesFromName(z, "IBCD_*")
        for IBCD in allIBCD:
            znames = Internal.getValue(IBCD)
            Zones.append(znames)
            xPC = Internal.getNodesFromName(IBCD,"CoordinateX_PC")[0][1]

            if znames in XPC:
                a = numpy.concatenate((XPC[znames][0],xPC))
                XPC[znames] = [a]
            else:
                XPC[znames] = [xPC]
    Zones = list(set(Zones))
    for zname in Zones:
        xPC = XPC[zname];
        size = xPC[0].shape[0]
        z          = Internal.getNodeFromName ( t, zname)
        o          = Internal.getNodeFromName1( z, '.Solver#ownData')
        param_int  = Internal.getNodeFromName1( o, 'Parameter_int')
        param_real = Internal.getNodeFromName1( o, 'Parameter_real')
        param_int[1][LBM_IBC_NUM] = size
    return None

#=============================================================================
# Extraction des pts IBM: retourne un arbre avec les coordonnees des
# pts IBM a corriger, paroi, miroirs
#=============================================================================
def extractIBMInfo(tc_in, IBCNames="IBCD_*", fileout=None):
    """Extract IBM informations in a pyTree."""
    if isinstance(tc_in, str): tc = Cmpi.convertFile2PyTree(tc_in, proc=Cmpi.rank)
    else: tc = tc_in

    XPC={}; YPC={}; ZPC={}
    XPW={}; YPW={}; ZPW={}
    XPI={}; YPI={}; ZPI={}
    Zones = []
    for z in Internal.getZones(tc):
        allIBCD = Internal.getNodesFromName(z, IBCNames)
        for IBCD in allIBCD:
            znames = Internal.getValue(IBCD)
            Zones.append(znames)
            xPC = Internal.getNodesFromName(IBCD,"CoordinateX_PC")[0][1]
            yPC = Internal.getNodesFromName(IBCD,"CoordinateY_PC")[0][1]
            zPC = Internal.getNodesFromName(IBCD,"CoordinateZ_PC")[0][1]
            xPI = Internal.getNodesFromName(IBCD,"CoordinateX_PI")[0][1]
            yPI = Internal.getNodesFromName(IBCD,"CoordinateY_PI")[0][1]
            zPI = Internal.getNodesFromName(IBCD,"CoordinateZ_PI")[0][1]
            xPW = Internal.getNodesFromName(IBCD,"CoordinateX_PW")[0][1]
            yPW = Internal.getNodesFromName(IBCD,"CoordinateY_PW")[0][1]
            zPW = Internal.getNodesFromName(IBCD,"CoordinateZ_PW")[0][1]
            if znames in XPW:
                a = numpy.concatenate((XPW[znames][0],xPW))
                XPW[znames] = [a]
                a = numpy.concatenate((YPW[znames][0],yPW))
                YPW[znames] = [a]
                a = numpy.concatenate((ZPW[znames][0],zPW))
                ZPW[znames] = [a]
                a = numpy.concatenate((XPI[znames][0],xPI))
                XPI[znames] = [a]
                a = numpy.concatenate((YPI[znames][0],yPI))
                YPI[znames] = [a]
                a = numpy.concatenate((ZPI[znames][0],zPI))
                ZPI[znames] = [a]
                a = numpy.concatenate((XPC[znames][0],xPC))
                XPC[znames] = [a]
                a = numpy.concatenate((YPC[znames][0],yPC))
                YPC[znames] = [a]
                a = numpy.concatenate((ZPC[znames][0],zPC))
                ZPC[znames] = [a]

            else:
                XPW[znames] = [xPW]
                YPW[znames] = [yPW]
                ZPW[znames] = [zPW]
                XPI[znames] = [xPI]
                YPI[znames] = [yPI]
                ZPI[znames] = [zPI]
                XPC[znames] = [xPC]
                YPC[znames] = [yPC]
                ZPC[znames] = [zPC]
    Zones = list(set(Zones))
    corrected = []; wall = []; interp = []
    t = C.newPyTree(['IBM','Wall','Image'])
    for zname in Zones:
        xPC = XPC[zname]; yPC = YPC[zname]; zPC = ZPC[zname]
        size = xPC[0].shape[0]
        zone = G.cart((0,0,0),(1,1,1),(size,1,1))
        zone[0] = 'correctedPts_'+zname

        XPC0 = Internal.getNodeFromName2(zone, 'CoordinateX')
        XPC0[1] = xPC[0]
        YPC0 = Internal.getNodeFromName2(zone, 'CoordinateY')
        YPC0[1] = yPC[0]
        ZPC0 = Internal.getNodeFromName2(zone, 'CoordinateZ')
        ZPC0[1] = zPC[0]
        corrected.append(zone)
        #
        xPI = XPI[zname]; yPI = YPI[zname]; zPI = ZPI[zname]
        size = xPI[0].shape[0]
        zone = G.cart((0,0,0),(1,1,1),(size,1,1))
        zone[0] = 'interpPts_'+zname

        XPI0 = Internal.getNodeFromName2(zone, 'CoordinateX')
        XPI0[1] = xPI[0]
        YPI0 = Internal.getNodeFromName2(zone, 'CoordinateY')
        YPI0[1] = yPI[0]
        ZPI0 = Internal.getNodeFromName2(zone,'CoordinateZ')
        ZPI0[1] = zPI[0]
        interp.append(zone)

        xPW = XPW[zname]; yPW = YPW[zname]; zPW = ZPW[zname]
        size = xPW[0].shape[0]
        zone = G.cart((0,0,0),(1,1,1),(size,1,1))
        zone[0] = 'wallPts_'+zname

        XPW0 = Internal.getNodeFromName2(zone, 'CoordinateX')
        XPW0[1] = xPW[0]
        YPW0 = Internal.getNodeFromName2(zone, 'CoordinateY')
        YPW0[1] = yPW[0]
        ZPW0 = Internal.getNodeFromName2(zone, 'CoordinateZ')
        ZPW0[1] = zPW[0]
        wall.append(zone)

    t[2][1][2] = corrected; t[2][2][2] = wall; t[2][3][2] = interp
    t = C.convertArray2Node(t)

    # Avoid that two procs write the same information
    if Cmpi.size > 1:
        for z in Internal.getZones(t):
            if int(z[0][-1]) != Cmpi.rank:
                z[0] = z[0]+"%"+str(Cmpi.rank)

    if isinstance(fileout, str): Cmpi.convertPyTree2File(t, fileout)

    return t

# =============================================================================
# Compute the IBM points: wall, image, & target
# IN: loc: "centers" or "nodes"
# IN: hi:NOT SURE BUT USED WHEN tb IS NOT PROVIDED.
# IN: he:NOT SURE BUT USED WHEN tb IS NOT PROVIDED.
# IN: tb: immersed bodies
# IN: tFront: IBM front for first set of image points (always needed)
# IN: frontType: type of front= 0,1,42
# IN: cellNName: name of cellN used to get IBM point.
# IN: IBCType: type of IBM method: -1 points inside solid; 1 point in fluid
# IN: depth : depth of overlaps & some IBM points
# IN: Reynolds: Reynolds Number. Need as this will be used to get the modelisation height for frontType 42
# IN: yplus: target y+. Need as this will be used to get the modelisation height for frontType 42
# IN: Lref: reference length: Need as this will be used to get the modelisation height for frontType 42
# IN: hmod: [AJ] not sure. It is not used. see [BC].
# IN: isLBM: is it an LBM run?
# IN: isWireModel: is the Wire Mesh Model used?
# IN: isOrthoFirst: Apply orthogonal projection first and level set normal decscent second. Needed for IBC filaments only.
# IN: twoFront: Needed to change the output of the projection write files.
# =============================================================================
def getAllIBMPoints(t, loc='nodes', hi=0., he=0., tb=None, tfront=None, frontType=0, cellNName='cellN',
                    IBCType=1, depth=2, Reynolds=6.e6, yplus=100., Lref=1.,
                    hmod=0.1, isLBM=False, isWireModel=False, isOrthoFirst=False, check=False, twoFronts=False):
    """Returns the dictionary of IBM points."""
    if IBCType == -1: signOfDistCorrected = -1
    else: signOfDistCorrected=1 # signe de la distance aux points corriges

    allCorrectedPts = []; allWallPts = []; allInterpPts = []
    #allInterpPts_layer2 = []; # [AJ] Keep for now
    #-------------------------------------------
    # 1. Get the list of IBC corrected pts
    #-------------------------------------------
    listOfSnearsLoc = []
    listOfModelisationHeightsLoc = []
    if loc == 'nodes':
        for z in Internal.getZones(t):
            an = C.getFields(Internal.__GridCoordinates__, z)[0]
            ac1 = C.getField(cellNName, z)[0]
            ac1[0] = 'cellN'
            ac2 = C.getField('TurbulentDistance', z)[0]
            ac3 = C.getField('gradxTurbulentDistance', z)[0]
            ac4 = C.getField('gradyTurbulentDistance', z)[0]
            ac5 = C.getField('gradzTurbulentDistance', z)[0]
            an = Converter.addVars([an,ac1,ac2,ac3,ac4,ac5])
            ah = C.getField('hi', z)[0]
            if ah != []: an = Converter.addVars([an,ah])
            ah = C.getField('he', z)[0]
            if ah != []: an = Converter.addVars([an,ah])
            correctedPts = Connector.getInterpolatedPoints__(an)
            allCorrectedPts.append(correctedPts)
            xt = C.getField('CoordinateX', z)[0][1][0]
            snearl = xt[1]-xt[0]
            if twoFronts: snearl *= 2
            listOfSnearsLoc.append(snearl)

            if frontType == 42:
                hmod = G_IBM_Height.computeModelisationHeight(Re=Reynolds, yplus=yplus, L=Lref)
                if twoFronts: hmod *= 2
                listOfModelisationHeightsLoc.append(hmod)
            else:
                listOfModelisationHeightsLoc.append(0.)
    else:
        for z in Internal.getZones(t):
            an = C.getFields(Internal.__GridCoordinates__,z)[0]
            an = Converter.node2Center(an)
            ac1 = C.getField('centers:'+cellNName,z)[0]
            ac1[0] = 'cellN'
            ac2 = C.getField('centers:TurbulentDistance',z)[0]
            ac3 = C.getField('centers:gradxTurbulentDistance',z)[0]
            ac4 = C.getField('centers:gradyTurbulentDistance',z)[0]
            ac5 = C.getField('centers:gradzTurbulentDistance',z)[0]
            an = Converter.addVars([an,ac1,ac2,ac3,ac4,ac5])
            ah = C.getField('centers:hi',z)[0]
            if ah != []: an = Converter.addVars([an,ah])
            ah = C.getField('centers:he',z)[0]
            if ah != []: an = Converter.addVars([an,ah])
            correctedPts = Connector.getInterpolatedPoints__(an)
            allCorrectedPts.append(correctedPts)
            xt = C.getField('CoordinateX',z)[0][1][0]
            snearl = xt[1]-xt[0]
            if twoFronts: snearl *= 2
            listOfSnearsLoc.append(snearl)

            if frontType == 42:
                hmod = G_IBM_Height.computeModelisationHeight(Re=Reynolds, yplus=yplus, L=Lref)
                if twoFronts: hmod *= 2
                listOfModelisationHeightsLoc.append(hmod)
            else:
                listOfModelisationHeightsLoc.append(0.)
    #-------------------------------------------
    # 2. Get the list of IBC wall and interp pts
    #-------------------------------------------
    indcell = Converter.extractVars(allCorrectedPts,['indcell'])
    if tb is None or tfront is None: # constant hi, he
        for nozc in range(len(allCorrectedPts)):
            poshe = KCore.isNamePresent(allCorrectedPts[nozc],'he')
            if poshe == -1: allCorrectedPts[nozc] = Converter.initVars(allCorrectedPts[nozc],'he',he)
            poshi = KCore.isNamePresent(allCorrectedPts[nozc],'hi')
            if poshi == -1: allCorrectedPts[nozc] = Converter.initVars(allCorrectedPts[nozc],'hi',hi)
        res = connector.getIBMPtsBasic(allCorrectedPts, varsn, 'TurbulentDistance')
    else:
        dictOfBodiesByIBCType={}
        for s in Internal.getZones(tb):
            sdd = Internal.getNodeFromName1(s,".Solver#define")
            if sdd is not None: # check consistency of ibc type with flow equations
                ibctype = Internal.getNodeFromName1(sdd, "ibctype")
                if ibctype is not None:
                    ibctype = Internal.getValue(ibctype)
                else:
                    if IBCType == -1: ibctype = 'slip'
                    else: ibctype = 'Musker'
            else: # type of IBC not found: Euler -> slip, other : Musker
                if IBCType == -1: ibctype = 'slip'
                else: ibctype = 'Musker'
	    #Over ride as LBM only supports no slip at the moment
            if isLBM:
                ibctype = 'noslip'
            famName = Internal.getNodeFromType1(s, 'FamilyName_t')
            if famName is not None:
                famName = Internal.getValue(famName)

            ibctypeI = TypesOfIBC[ibctype]
            if famName is not None: ibctype2 = str(ibctypeI)+"#"+famName
            else: ibctype2 = str(ibctypeI)
            if ibctype2 not in dictOfBodiesByIBCType: dictOfBodiesByIBCType[ibctype2]=[s]
            else: dictOfBodiesByIBCType[ibctype2]+=[s]

        # Regroupement des corps par type de BC - optimise les projections ensuite
        bodies = []; listOfIBCTypes=[]
        for itype in dictOfBodiesByIBCType:
            s = dictOfBodiesByIBCType.get(itype)
            body = C.getFields(Internal.__GridCoordinates__,s)
            body = Converter.convertArray2Tetra(body)
            body = Transform.join(body)
            bodies.append(body)
            listOfIBCTypes.append(itype)

        if frontType == 0:
            dmin = C.getMaxValue(tfront, 'TurbulentDistance')
            allCorrectedPts = Converter.initVars(allCorrectedPts,'dist',dmin)
            res = connector.getIBMPtsWithoutFront(allCorrectedPts, bodies, varsn, 'dist', signOfDistCorrected)
        else:
            front = C.getFields(Internal.__GridCoordinates__,tfront)

            front = Converter.convertArray2Tetra(front)
            allCorrectedPts = Converter.extractVars(allCorrectedPts,['CoordinateX','CoordinateY','CoordinateZ']+varsn)
            res = connector.getIBMPtsWithFront(allCorrectedPts, listOfSnearsLoc, listOfModelisationHeightsLoc, bodies,
                                               front, varsn, signOfDistCorrected, depth, int(isWireModel), int(isOrthoFirst))
    allWallPts = res[0]
    allWallPts = Converter.extractVars(allWallPts,['CoordinateX','CoordinateY','CoordinateZ'])

    allInterpPts = res[1]
    allInterpPts = Converter.extractVars(allInterpPts,['CoordinateX','CoordinateY','CoordinateZ'])
    allInterpPts = Converter.addVars([allInterpPts,indcell])
    allCorrectedPts = Converter.extractVars(allCorrectedPts,['CoordinateX','CoordinateY','CoordinateZ'])

    dictOfInterpPtsByIBCType={}
    dictOfCorrectedPtsByIBCType={}
    dictOfWallPtsByIBCType={}
    nzonesR = len(allInterpPts)

    ## Ouput the IBM points that have a type 3 and type 4 projection
    if len(res)>3:
        nameZone = ['IBM', 'Wall', 'Image']
        tLocal3   = C.newPyTree(nameZone)
        tLocal4   = C.newPyTree(nameZone)
        isWrite3  = 0
        isWrite4  = 0
        allProjectPts = res[3]
        allProjectPts = Converter.extractVars(allProjectPts,['ProjectionType'])
        outputProjection3 = [[],[],[],[],[],[],[],[],[]]
        outputProjection4 = [[],[],[],[],[],[],[],[],[]]
        for noz in range(nzonesR):
            arrayLocal = allProjectPts[noz][1][0]
            type_3     = numpy.count_nonzero(arrayLocal==3)
            type_4     = numpy.count_nonzero(arrayLocal==4)
            if type_3>0: _prepOutputProject__(outputProjection3, 3, arrayLocal, allCorrectedPts[noz][1], allWallPts[noz][1], allInterpPts[noz][1])
            if type_4>0: _prepOutputProject__(outputProjection4, 4, arrayLocal, allCorrectedPts[noz][1], allWallPts[noz][1], allInterpPts[noz][1])

        if outputProjection3[0] and check:
            tLocal3  = _writeOutputProject__(outputProjection3, tLocal3)
            isWrite3 = 1
        if outputProjection4[0] and check:
            tLocal4  = _writeOutputProject__(outputProjection4, tLocal4)
            isWrite4 = 1

        if check:
            suffixLocal = '.cgns'
            if twoFronts: suffixLocal = '_Front2.cgns'
            if cellNName=='cellNFilWMM':
                suffixLocal = '_WMM_Front1.cgns'
                if isWireModel: suffixLocal = '_WMM_Front2.cgns'
            if Cmpi.allreduce(isWrite3,op=Cmpi.MAX)>0: Cmpi.convertPyTree2File(tLocal3, 'projection3'+suffixLocal)
            if Cmpi.allreduce(isWrite4,op=Cmpi.MAX)>0: Cmpi.convertPyTree2File(tLocal4, 'projection4'+suffixLocal)

        del tLocal3
        del tLocal4

    if len(res)>=3:
        allIndicesByIBCType = res[2]
        for noz in range(nzonesR):
            indicesByTypeForZone = res[2][noz]
            nbTypes = len(indicesByTypeForZone)
            for nob in range(nbTypes):
                ibcTypeL = listOfIBCTypes[nob]
                indicesByTypeL = indicesByTypeForZone[nob]
                if indicesByTypeL.shape[0] > 0:
                    correctedPtsL = Transform.subzone(allCorrectedPts[noz], indicesByTypeL)
                    interpPtsL = Transform.subzone(allInterpPts[noz], indicesByTypeL)
                    wallPtsL = Transform.subzone(allWallPts[noz], indicesByTypeL)
                else:
                    correctedPtsL=[]; interpPtsL = []; wallPtsL = []

                if noz == 0:
                    dictOfCorrectedPtsByIBCType[ibcTypeL] = [correctedPtsL]
                    dictOfWallPtsByIBCType[ibcTypeL] = [wallPtsL]
                    dictOfInterpPtsByIBCType[ibcTypeL] = [interpPtsL]
                else:
                    dictOfCorrectedPtsByIBCType[ibcTypeL] += [correctedPtsL]
                    dictOfWallPtsByIBCType[ibcTypeL] += [wallPtsL]
                    dictOfInterpPtsByIBCType[ibcTypeL] += [interpPtsL]
    else:
        if IBCType == -1: ibcTypeL = 0
        else: ibcTypeL = 3
        for noz in range(nzonesR):
            if noz == 0:
                dictOfCorrectedPtsByIBCType[ibcTypeL] = [allCorrectedPts[noz]]
                dictOfWallPtsByIBCType[ibcTypeL] = [allWallPts[noz]]
                dictOfInterpPtsByIBCType[ibcTypeL] = [allInterpPts[noz]]
            else:
                dictOfCorrectedPtsByIBCType[ibcTypeL] += [allCorrectedPts[noz]]
                dictOfWallPtsByIBCType[ibcTypeL] += [allWallPts[noz]]
                dictOfInterpPtsByIBCType[ibcTypeL] += [allInterpPts[noz]]

    return dictOfCorrectedPtsByIBCType, dictOfWallPtsByIBCType, dictOfInterpPtsByIBCType

#=============================================================================
# Extraction des infos pour adaptation de front
# [BC] to complete
#=============================================================================
def createWallAdapt(tc):
    """Create a cloud of target points from y+ information."""
    listOfZones = []
    DictOfZones = {}

    for zc in Internal.getZones(tc):
        subRegions = Internal.getNodesFromType1(zc, 'ZoneSubRegion_t')
        for s in subRegions:
            sname = s[0][0:4]
            zname = s[0].split('_')[-1]
            # zname = zname.split('X')[0]
            if sname == 'IBCD':
                if zname not in DictOfZones:
                    DictOfZones[zname] = {"CoordinateX_PC":[], "CoordinateY_PC":[], "CoordinateZ_PC":[], "yplus":[], "PointListDonor":[]}

                xPC = Internal.getNodeFromName1(s,"CoordinateX_PC")[1]
                yPC = Internal.getNodeFromName1(s,"CoordinateY_PC")[1]
                zPC = Internal.getNodeFromName1(s,"CoordinateZ_PC")[1]
                DictOfZones[zname]["CoordinateX_PC"].append(xPC)
                DictOfZones[zname]["CoordinateY_PC"].append(yPC)
                DictOfZones[zname]["CoordinateZ_PC"].append(zPC)

                yplus = Internal.getNodeFromName1(s, "yplus")[1]
                DictOfZones[zname]["yplus"].append(yplus)

                pointListDnr = Internal.getNodeFromName1(s, "PointListDonor")[1]
                DictOfZones[zname]["PointListDonor"].append(pointListDnr)

    for zname in DictOfZones:
        xcNP = numpy.concatenate(DictOfZones[zname]["CoordinateX_PC"])
        ycNP = numpy.concatenate(DictOfZones[zname]["CoordinateY_PC"])
        zcNP = numpy.concatenate(DictOfZones[zname]["CoordinateZ_PC"])
        yplusNP = numpy.concatenate(DictOfZones[zname]["yplus"])
        pointListNP = numpy.concatenate(DictOfZones[zname]["PointListDonor"])

        # Creation d une seule zone
        zsize = numpy.empty((1,3), Internal.E_NpyInt, order='F')
        zsize[0,0] = xcNP.shape[0]; zsize[0,1] = 0; zsize[0,2] = 0
        z = Internal.newZone(name=zname,zsize=zsize,ztype='Unstructured')
        gc = Internal.newGridCoordinates(parent=z)
        coordx = ['CoordinateX',xcNP,[],'DataArray_t']
        coordy = ['CoordinateY',ycNP,[],'DataArray_t']
        coordz = ['CoordinateZ',zcNP,[],'DataArray_t']
        gc[2] = [coordx,coordy,coordz]
        n = Internal.createChild(z, 'GridElements', 'Elements_t', [2,0])
        Internal.createChild(n, 'ElementRange', 'IndexRange_t', [1,0])
        Internal.createChild(n, 'ElementConnectivity', 'DataArray_t', None)
        gridLocation = 'Vertex'
        FSN = Internal.newFlowSolution(name=Internal.__FlowSolutionNodes__,
                                       gridLocation=gridLocation, parent=z)

        FSN[2].append(["yplus",yplusNP, [],'DataArray_t'])
        FSN[2].append(["PointListDonor",pointListNP, [],'DataArray_t'])
        listOfZones.append(z)

    t = C.newPyTree(['Base', listOfZones])
    return t

#=============================================================================
# Creation des zones 1D de nom 'Zone#IBCD_*' pour pouvoir etre identifie en post
# retourne tw avec une zone 1D par IBCD de depart
# Sert au post traitement paroi en parallele et instationnaire
#=============================================================================
def createIBMWZones(tc, variables=[]):
    """Extract IBM informations from tc before post-processing."""
    import Converter.Mpi as Cmpi
    tw = C.newPyTree(['IBM_WALL'])
    for z in Internal.getZones(tc):
        ZSR = Internal.getNodesFromType2(z,'ZoneSubRegion_t')
        procNode = Internal.getNodeFromName(z,'proc')
        proc = -1
        if procNode is not None:
            proc = Internal.getValue(procNode)

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
            if proc != -1: Cmpi._setProc(zw,proc)
            tw[2][1][2].append(zw)
    return tw

#=============================================================================
# distance signee en fonction du masquage Chimere et IBC
#=============================================================================
def _signDistance(t):
    """Modify distance to be signed based on cellNIBC and cellNChim (outside fluid domain)."""
    C._initVars(t,'{centers:TurbulentDistance}=-1.*({centers:cellNIBC}*{centers:cellNChim}<1.)*{centers:TurbulentDistance}+({centers:cellNIBC}*{centers:cellNChim}>0.)*{centers:TurbulentDistance}')
    return None

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## TO BE DELETED?
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# isosurface of max dist of the first interpolable cells
def getIBMFrontType0Old__(tc, frontvar, dim, isFront2=False, frontType=0, SHIFTB=0.):
    import Converter.Mpi as Cmpi

    if dim == 2:
        z0 = Internal.getNodeFromType2(tc, 'Zone_t')
        zmean = C.getValue(z0, 'CoordinateZ', 0)
        dz = 2*zmean
    else: dz = 0.

    SHIFTD = 1.+SHIFTF
    # SHIFTD *= SHIFTB
    front = []
    tf = Internal.copyRef(tc)
    C._initVars(tf,'{%s}={%s}-2.*({%s}>1.5)'%(frontvar,frontvar,frontvar))
    if isFront2:
        for z in Internal.getZones(tf):
            epsilon_dist = abs(C.getValue(z,'CoordinateX',1)-C.getValue(z,'CoordinateX',0))
            if epsilon_dist < SHIFTB:
                if C.getMinValue(z,frontvar)==0. and C.getMaxValue(z,frontvar)==1.:
                    f = P.frontFaces(z, frontvar)
                    if Internal.getZoneDim(f)[1]>0:
                        Internal._rmNodesByName(f,'ID_*')
                        front.append(f)

        if dim == 2:
            if frontType != 1:
                dmin = C.getMaxValue(front, 'TurbulentDistance')

                if Cmpi.KCOMM is not None:
                    dmin = numpy.array([dmin], dtype=numpy.float64)
                    dmin_max = numpy.zeros(1, dtype=numpy.float64)
                    Cmpi.KCOMM.Allreduce(dmin, dmin_max, Cmpi.MAX)
                    dmin = dmin_max[0]
            else:
                dmin = 2*SHIFTB

            if frontType == 42: dmin += SHIFTB
            front = []
            # Creation du corps 2D pour le preprocessing IBC
            for z in Internal.getZones(tf):
                epsilon_dist = abs(C.getValue(z,'CoordinateX',1)-C.getValue(z,'CoordinateX',0))
                if epsilon_dist < SHIFTB:
                    tcl = T.addkplane(z)
                    T._translate(tcl,(0,0,-zmean))
                    T._contract(tcl, (0,0,0), (1,0,0), (0,1,0), dz)
                    front.append(P.isoSurfMC(tcl,'TurbulentDistance',dmin*SHIFTD))
                    del tcl
            front = C.newPyTree(['Base']+front)
        else:
            dmin = C.getMaxValue(front, 'TurbulentDistance')

            if Cmpi.KCOMM is not None:
                dmin = numpy.array([dmin], dtype=numpy.float64)
                dmin_max = numpy.zeros(1, dtype=numpy.float64)
                Cmpi.KCOMM.Allreduce(dmin, dmin_max, Cmpi.MAX)
                dmin = dmin_max[0]

            if frontType == 42: dmin += SHIFTB
            front = []
            # Creation du corps 2D pour le preprocessing IBC
            for z in Internal.getZones(tf):
                epsilon_dist = abs(C.getValue(z,'CoordinateX',1)-C.getValue(z,'CoordinateX',0))
                if epsilon_dist < SHIFTB:
                    front.append(P.isoSurfMC(z,'TurbulentDistance',dmin*SHIFTD))
            front = C.newPyTree(['Base']+front)

    else:
        for z in Internal.getZones(tf):
            if C.getMinValue(z,frontvar)==0. and C.getMaxValue(z,frontvar)==1.:
                f = P.frontFaces(z, frontvar)
                if Internal.getZoneDim(f)[1]>0:
                    Internal._rmNodesByName(f,'ID_*')
                    front.append(f)
        if dim == 2:
            dmin = C.getMaxValue(front, 'TurbulentDistance')

            if Cmpi.KCOMM is not None:
                dmin = numpy.array([dmin], dtype=numpy.float64)
                dmin_max = numpy.zeros(1, dtype=numpy.float64)
                Cmpi.KCOMM.Allreduce(dmin, dmin_max, Cmpi.MAX)
                dmin = dmin_max[0]

            if frontType == 42: dmin += SHIFTB
            # Creation du corps 2D pour le preprocessing IBC
            tcl = T.addkplane(tc)
            T._translate(tcl,(0,0,-zmean))
            T._contract(tcl, (0,0,0), (1,0,0), (0,1,0), dz)
            front = P.isoSurfMC(tcl,'TurbulentDistance',dmin*SHIFTD)
            del tcl
        else:
            dmin = C.getMaxValue(front, 'TurbulentDistance')

            if Cmpi.KCOMM is not None:
                dmin = numpy.array([dmin], dtype=numpy.float64)
                dmin_max = numpy.zeros(1, dtype=numpy.float64)
                Cmpi.KCOMM.Allreduce(dmin, dmin_max, Cmpi.MAX)
                dmin = dmin_max[0]

            if frontType == 42: dmin += SHIFTB
            front = P.isoSurfMC(tc, 'TurbulentDistance', dmin*SHIFTD)
    return front

#=============================================================================
# Performs smoothing after _pushBackImageFront2
# Only for 2D
#=============================================================================
def _smoothImageFront(t, tc, dimPb=2):
    """Add target points to smooth the front."""
    for z in Internal.getZones(t):
        cellNFront  = Internal.getNodeFromName(z,'cellNFront')
        cellNIBC    = Internal.getNodeFromName(z,'cellNIBC')
        cellNIBCDnr = Internal.getNodeFromName(z,'cellNIBCDnr')
        dim         = Internal.getZoneDim(z)
        k           = 0
        if cellNIBC != []:
            cellNIBC = cellNIBC[1]
            sizeTot = cellNIBC.shape[0]*cellNIBC.shape[1]*cellNIBC.shape[2]
            sizeOne =  int(numpy.sum(cellNIBC))
            # if sizeOne < sizeTot:
            cellNFront  = cellNFront[1]
            cellNIBCDnr = cellNIBCDnr[1]
            for i in range(1, int(dim[1])-2):
                for j in range(1, int(dim[2])-2):
                    if cellNIBC[i,j,k] == 1 and cellNIBC[i,j-1,k] == 2:
                        if cellNIBC[i-1,j+1,k] == 2:
                            cellNFront[i,j,k]    = 0
                            cellNIBC[i,j,k]      = 2
                            cellNIBCDnr[i,j,k]   = 2

                for j in range(int(dim[2])-3, 0, -1):
                    if cellNIBC[i,j,k] == 1 and cellNIBC[i,j+1,k] == 2:
                        if cellNIBC[i-1,j-1,k] == 2:
                            cellNFront[i,j,k]    = 0
                            cellNIBC[i,j,k]      = 2
                            cellNIBCDnr[i,j,k]   = 2

            for i in range(int(dim[1])-3, 0, -1):
                for j in range(1, int(dim[2])-2):
                    if cellNIBC[i,j,k] == 1 and cellNIBC[i,j-1,k] == 2:
                        if cellNIBC[i+1,j+1,k] == 2:
                            cellNFront[i,j,k]    = 0
                            cellNIBC[i,j,k]      = 2
                            cellNIBCDnr[i,j,k]   = 2

                for j in range(int(dim[2])-3, 0, -1):
                    if cellNIBC[i,j,k] == 1 and cellNIBC[i,j+1,k] == 2:
                        if cellNIBC[i+1,j-1,k] == 2:
                            cellNFront[i,j,k]    = 0
                            cellNIBC[i,j,k]      = 2
                            cellNIBCDnr[i,j,k]   = 2

    C._cpVars(t,'centers:cellNIBC',tc,'cellNIBC')
    C._cpVars(t,'centers:cellNIBC',t,'centers:cellN')
    C._cpVars(t,'centers:cellN',tc,'cellN')

    return None

#=============================================================================
# Compute curvature parameter "K" from geom in 2D
#=============================================================================
def _computeKcurvParameter(tc, tb):
    """Compute curvature parameter for simple 2D geometry."""

    ################################
    ## EXTRADOS ##
    ################################

    tb_extrados = P.selectCells(tb, '{CoordinateY}>0.')
    x_extrados  = Internal.getNodeFromName(tb_extrados, 'CoordinateX')[1]
    y_extrados  = Internal.getNodeFromName(tb_extrados, 'CoordinateY')[1]

    inds_extrados = x_extrados.argsort()
    x_extrados    = x_extrados[inds_extrados]
    y_extrados    = y_extrados[inds_extrados]

    # evite les divisions par 0
    nPts = len(x_extrados)
    listOfPos = []
    for i in range(1,nPts):
        if x_extrados[i] == x_extrados[i-1]:
            listOfPos.append(i)

    for pos in listOfPos:
        x_extrados = numpy.delete(x_extrados, pos)
        y_extrados = numpy.delete(y_extrados, pos)

    nPts = len(x_extrados)
    y1_extrados = numpy.zeros(nPts, dtype=numpy.float64)
    y2_extrados = numpy.zeros(nPts, dtype=numpy.float64)

    # derivee premiere y1
    # schema decentre pour CL
    y1_extrados[0] = (y_extrados[1]-y_extrados[0])/((x_extrados[1]-x_extrados[0]))
    y1_extrados[-1] = (y_extrados[-2]-y_extrados[-1])/((x_extrados[-2]-x_extrados[-1]))

    # schema centre
    for i in range(1, nPts-1):
        y1_extrados[i] = (y_extrados[i+1]-y_extrados[i-1])/((x_extrados[i+1]-x_extrados[i-1]))

    # derivee seconde y2
    # schema decentre pour CL
    y2_extrados[0] = (y1_extrados[1]-y1_extrados[0])/((x_extrados[1]-x_extrados[0]))
    y2_extrados[-1] = (y1_extrados[-2]-y1_extrados[-1])/((x_extrados[-2]-x_extrados[-1]))

    # schema centre
    for i in range(1, nPts-1):
        y2_extrados[i] = (y1_extrados[i+1]-y1_extrados[i-1])/((x_extrados[i+1]-x_extrados[i-1]))

    k_extrados = y2_extrados/(numpy.power(1 + numpy.power(y1_extrados, 2), 1.5))
    ka_extrados = numpy.abs(y2_extrados)/(numpy.power(1 + numpy.power(y1_extrados, 2), 1.5))

    ################################
    ## INTRADOS ##
    ################################

    tb_intrados = P.selectCells(tb, '{CoordinateY}<0.')
    x_intrados = Internal.getNodeFromName(tb_intrados, 'CoordinateX')[1]
    y_intrados = numpy.abs(Internal.getNodeFromName(tb_intrados, 'CoordinateY')[1]) # abs pour avoir vraie valeur k

    inds_intrados = x_intrados.argsort()
    x_intrados = x_intrados[inds_intrados]
    y_intrados = y_intrados[inds_intrados]

    # evite les divisions par 0
    nPts = len(x_intrados)
    listOfPos = []
    for i in range(1,nPts):
        if x_intrados[i] == x_intrados[i-1]:
            listOfPos.append(i)

    for pos in listOfPos:
        x_intrados = numpy.delete(x_intrados, pos)
        y_intrados = numpy.delete(y_intrados, pos)

    nPts = len(x_intrados)
    y1_intrados = numpy.zeros(nPts)
    y2_intrados = numpy.zeros(nPts)

    # derivee premiere y1
    # schema decentre pour CL
    y1_intrados[0] = (y_intrados[1]-y_intrados[0])/((x_intrados[1]-x_intrados[0]))
    y1_intrados[-1] = (y_intrados[-2]-y_intrados[-1])/((x_intrados[-2]-x_intrados[-1]))

    # schema centre
    for i in range(1, nPts-1):
        y1_intrados[i] = (y_intrados[i]-y_intrados[i-1])/((x_intrados[i]-x_intrados[i-1]))

    # derivee seconde y2
    # schema decentre pour CL
    y2_intrados[0] = (y1_intrados[1]-y1_intrados[0])/((x_intrados[1]-x_intrados[0]))
    y2_intrados[-1] = (y1_intrados[-2]-y1_intrados[-1])/((x_intrados[-2]-x_intrados[-1]))

    # schema centre
    for i in range(1, nPts-1):
        y2_intrados[i] = (y1_intrados[i]-y1_intrados[i-1])/((x_intrados[i]-x_intrados[i-1]))

    k_intrados = y2_intrados/(numpy.power(1 + numpy.power(y1_intrados, 2), 1.5))
    ka_intrados = numpy.abs(y2_intrados)/(numpy.power(1 + numpy.power(y1_intrados, 2), 1.5))

    ################################
    ## MaJ tc ##
    ################################
    for z in Internal.getZones(tc):
        for zsr in Internal.getNodesFromType1(z, 'ZoneSubRegion_t'):
            nameSubRegion = zsr[0]
            if nameSubRegion[0:4]=='IBCD':
                ibctypeCR = nameSubRegion.split('_')[1]
                if ibctypeCR=='100':
                    KCurv = Internal.getNodeFromName(zsr, XOD.__KCURV__)[1]
                    coordX = Internal.getNodeFromName(zsr, 'CoordinateX_PW')[1]
                    coordY = Internal.getNodeFromName(zsr, 'CoordinateY_PW')[1]
                    nIBC = numpy.shape(coordX)[0]

                    for i in range(nIBC):
                        if coordY[i] > 0:
                            j = 0
                            while(coordX[i] > x_extrados[j] and j < numpy.shape(x_extrados)[0]-1):
                                j += 1
                            KCurv[i] = k_extrados[j]
                        else:
                            j = 0
                            while(coordX[i] > x_intrados[j] and j < numpy.shape(x_intrados)[0]-1):
                                j += 1
                            KCurv[i] = k_intrados[j]
                        KCurv[i] = min(KCurv[i], 100.)
                        KCurv[i] = max(KCurv[i], -100.)

                    Internal.getNodeFromName(zsr, XOD.__KCURV__)[1] = KCurv

    return None

#=============================================================================
# Performs specific treatment for frontType=2
# Modify/add target points in matching boundary regions
#=============================================================================
def _pushBackImageFront2(t, tc, tbb, interpDataType=1):
    """Modify/add target points in matching boundary regions for frontType=2."""
    intersectionsDict = X.getIntersectingDomains(tbb, method='AABB', taabb=tbb)
    C._initVars(t,'{centers:cellNFront_origin}={centers:cellNFront}')
    C._initVars(t,'{centers:cellNIBC_origin}={centers:cellNIBC}')
    C._cpVars(t,'centers:cellNFront',tc,'cellNFront')
    C._cpVars(t,'centers:cellNFront_origin',tc,'cellNFront_origin')
    C._cpVars(t,'centers:cellNIBC_origin',tc,'cellNIBC_origin')

    C._initVars(t,'{centers:cellNFront2}=1.-({centers:cellNFront}<1.)*(abs({centers:cellNChim})>1.)')

    for z in Internal.getZones(t):
        cellNFront = Internal.getNodeFromName2(z,'cellNFront2')
        if cellNFront != []:
            cellNFront = cellNFront[1]
            sizeTot = cellNFront.shape[0]*cellNFront.shape[1]*cellNFront.shape[2]
            sizeOne =  int(numpy.sum(cellNFront))
            if sizeOne < sizeTot:
                X._setHoleInterpolatedPoints(z, depth=1, dir=0, loc='centers',cellNName='cellNFront2',addGC=False)
                res = X.getInterpolatedPoints(z,loc='centers', cellNName='cellNFront2') #indices,X,Y,Z
                if res is not None:
                    indicesI = res[0]
                    XI = res[1]; YI = res[2]; ZI = res[3]
                    allInterpFields=[]
                    for zc in Internal.getZones(tc):
                        if zc[0] in intersectionsDict[z[0]]:
                            if interpDataType==1: HOOKADT = C.createHook(zc, 'adt')
                            else: HOOKADT = None
                            fields = X.transferFields(zc, XI, YI, ZI, hook=HOOKADT, variables=['cellNFront_origin','cellNIBC_origin'], interpDataType=interpDataType)
                            allInterpFields.append(fields)
                            if interpDataType == 1: C.freeHook(HOOKADT)
                    if allInterpFields != []:
                        C._filterPartialFields(z, allInterpFields, indicesI, loc='centers', startFrom=0, filterName='donorVol', verbose=False)
                        # C._initVars(z,'{centers:cellNFront}=({centers:cellNFront}>0.5)') #ancienne version
                        C._initVars(z,'{centers:cellNFront}={centers:cellNFront}*({centers:cellNFront_origin}>0.5)') # Modification du Front uniquement lorsque celui-ci est repousse
                        # i.e. if cellNFront_origin == 0 and cellNFront == 1 => cellNfront = 0

                        C._initVars(z,'{centers:cellNIBC}={centers:cellNIBC}*(1.-({centers:cellNChim}==1.)*({centers:cellNIBC_origin}>1.5)*({centers:cellNIBC_origin}<2.5)) \
                            + 2.*({centers:cellNChim}==1.)*({centers:cellNIBC_origin}>1.5)*({centers:cellNIBC_origin}<2.5)')
                        # i.e. if cellNChim == 1 and cellNIBC_origin == 2 => cellNIBC = 2

                        C._initVars(z,'{centers:cellNIBCDnr}={centers:cellNIBCDnr}*(1.-({centers:cellNChim}==1.)*({centers:cellNIBC_origin}>1.5)*({centers:cellNIBC_origin}<2.5)) \
                            + 2.*({centers:cellNChim}==1.)*({centers:cellNIBC_origin}>1.5)*({centers:cellNIBC_origin}<2.5)')
                        # le front donneur doit etre egalement maj sur l'exemple du champ cellNIBC

        C._rmVars(z,['centers:cellNFront2'])
        C._rmVars(z,['centers:cellNFront_origin'])
        C._rmVars(z,['centers:cellNIBC_origin'])

    C._cpVars(t,'centers:cellNIBC',tc,'cellNIBC')
    C._cpVars(t,'centers:cellNIBC',t,'centers:cellN')
    C._cpVars(t,'centers:cellN',tc,'cellN')

    return None

#=============================================================================
# setInterpDataIBM, modifie tc -> [BC] : to ve replaced by new version?
# IN: tbb: bboundingbox tree
# IN: tb: immersed bodies
# IN: typeI: type of Connectivity to do (ID=Chimere; IBCD=IBM)
# IN: dim: dimension
# IN: dictofADT: dictionary of alternative directivity tree (?)
# IN: front: front of upon which IBM image points are placed
# IN: frontType: type of front= 0,1,42
# IN: depth : depth of overlaps & some IBM points
# IN: IBCType: type of IBM method: -1 points inside solid; 1 point in fluid
# IN: interDataType: type of interpolation beteween zones: 0=Cartesian; 1=Chimere
# IN: Reynolds: Reynolds Number. Need as this will be used to get the modelisation height for frontType 42
# IN: yplus: target y+. Need as this will be used to get the modelisation height for frontType 42
# IN: Lref: reference length: Need as this will be used to get the modelisation height for frontType 42
# IN: isLBM: is the case LBM?
#=============================================================================
def doInterp(t, tc, tbb, tb=None, typeI='ID', dim=3, dictOfADT=None, front=None, frontType=0, depth=2, IBCType=1, interpDataType=1, Reynolds=6.e6, yplus=100., Lref=1., isLBM=False, check=False):
    ReferenceState = Internal.getNodeFromType2(t, 'ReferenceState_t')
    model = Internal.getNodeFromName(t, 'GoverningEquations')
    if model is not None: model = Internal.getValue(model)
    else:                 model = "Euler"
    if typeI == 'ID':
        # toutes les zones sont interpolables en Chimere
        intersectionsDict = X.getIntersectingDomains(tbb, method='AABB', taabb=tbb)
        rcvZones = []
        for zrcv in Internal.getZones(t):
            if C.getMaxValue(zrcv,'centers:cellN')==2.:
                zrcvname = zrcv[0]; rcvZones.append(zrcv)
        nozr = 0
        nbZonesChim = len(rcvZones)
        for nozr in range(nbZonesChim):
            zrcv = rcvZones[nozr]
            zrcvname = zrcv[0]
            nozr += 1; hook0 = []
            nobOfDnrBases = []; nobOfDnrZones=[]; dnrZones=[]
            for nobd in range(len(tc[2])):
                if tc[2][nobd][3] == 'CGNSBase_t':
                    for nozd in range(len(tc[2][nobd][2])):
                        zdnr = tc[2][nobd][2][nozd]
                        if zdnr[3] == 'Zone_t':
                            zdnrname = zdnr[0]
                            if zdnrname in intersectionsDict[zrcvname]:
                                nobOfDnrBases.append(nobd)
                                nobOfDnrZones.append(nozd)
                                dnrZones.append(zdnr)
                                if interpDataType==1 and dictOfADT is not None:
                                    hook0.append(dictOfADT[zdnrname])
            if interpDataType == 0: hook0 = None

            X._setInterpData(zrcv, dnrZones, nature=1,penalty=1,loc='centers',storage='inverse',sameName=1,\
                             interpDataType=interpDataType, itype='chimera')
            for nod in range(len(dnrZones)):
                nobd = nobOfDnrBases[nod]
                nozd = nobOfDnrZones[nod]
                tc[2][nobd][2][nozd] = dnrZones[nod]

    elif typeI == 'IBCD':
        # detection des zones IBC
        zonesRIBC = []
        for zrcv in Internal.getZones(t):
            if C.getMaxValue(zrcv,'centers:cellNIBC')==2.:
                zrcvname = zrcv[0]; zonesRIBC.append(zrcv)

        if zonesRIBC == []: return tc

        res = getAllIBMPoints(zonesRIBC, loc='centers', tb=tb, tfront=front, frontType=frontType, \
                              cellNName='cellNIBC', depth=depth, IBCType=IBCType, Reynolds=Reynolds, yplus=yplus, Lref=Lref, isLBM=isLBM, check=check)
        nbZonesIBC = len(zonesRIBC)
        dictOfADT = {}
        dictOfCorrectedPtsByIBCType = res[0]
        dictOfWallPtsByIBCType = res[1]
        dictOfInterpPtsByIBCType = res[2]
        for ibcTypeL in  dictOfCorrectedPtsByIBCType:
            allCorrectedPts = dictOfCorrectedPtsByIBCType[ibcTypeL]
            allWallPts = dictOfWallPtsByIBCType[ibcTypeL]
            allInterpPts = dictOfInterpPtsByIBCType[ibcTypeL]
            for nozr in range(nbZonesIBC):
                if allCorrectedPts[nozr] != []:
                    interpPtsBB=Generator.BB(allInterpPts[nozr])
                    zrcv = zonesRIBC[nozr]
                    zrcvname = zrcv[0]
                    nobOfDnrBases = []; nobOfDnrZones=[]; dnrZones=[]
                    if interpDataType == 1: hook0 = []
                    else: hook0 = None
                    for nobd in range(len(tc[2])):
                        if tc[2][nobd][3] == 'CGNSBase_t':
                            for nozd in range(len(tc[2][nobd][2])):
                                zdnr = tc[2][nobd][2][nozd]
                                if zdnr[3] == 'Zone_t':
                                    zdnrname = zdnr[0]
                                    zbb = tbb[2][nobd][2][nozd]
                                    bba = C.getFields(Internal.__GridCoordinates__,zbb)[0]
                                    if Generator.bboxIntersection(interpPtsBB,bba,isBB=True) == 1:
                                        if interpDataType == 1:
                                            if zdnrname not in dictOfADT:
                                                HOOKADT = C.createHook(zdnr, 'adt')
                                                dictOfADT[zdnrname] = HOOKADT
                                            hook0.append(dictOfADT[zdnrname])

                                        dnrZones.append(zdnr)
                                        nobOfDnrBases.append(nobd)
                                        nobOfDnrZones.append(nozd)

                    XOD._setIBCDataForZone__(zrcv, dnrZones, allCorrectedPts[nozr], allWallPts[nozr], allInterpPts[nozr], \
                                             nature=1, penalty=1, loc='centers', storage='inverse',
                                             interpDataType=interpDataType, hook=hook0, dim=dim, \
                                             ReferenceState=ReferenceState, bcType=ibcTypeL,model=model)
                    nozr += 1
                    for nod in range(len(dnrZones)):
                        nobd = nobOfDnrBases[nod]
                        nozd = nobOfDnrZones[nod]
                        tc[2][nobd][2][nozd] = dnrZones[nod]

        if dictOfADT is not None:
            for dnrname in dictOfADT: C.freeHook(dictOfADT[dnrname])

    return tc

def prepareIBMData_legacy(t, tbody, DEPTH=2, loc='centers', frontType=1, interpDataType=0, smoothing=False, yplus=100., Lref=1., wallAdapt=None, blankingF42=False, isLBM=False,LBMQ=False,isPrintDebug=False):
    tb =  Internal.copyRef(tbody)

    # tb: fournit model et dimension
    dimPb = Internal.getNodeFromName(tb,'EquationDimension')
    if dimPb is None: raise ValueError('prepareIBMData: EquationDimension is missing in input body tree.')
    dimPb = Internal.getValue(dimPb)

    # type de traitement paroi: pts interieurs ou externes
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('prepareIBMData: GoverningEquations is missing in input body tree.')
    # model: Euler, NSLaminar, NSTurbulent, LBMLaminar
    model = Internal.getValue(model)

    if model == 'Euler': IBCType =-1
    elif model == 'LBMLaminar':
        IBCType = -1
        if LBMQ: IBCType =1
    else: IBCType = 1 # Points cibles externes
    if loc == 'nodes':
        raise NotImplemented("prepareIBMData: prepareIBMData at nodes not yet implemented.")

    #------------------------
    # Ghost cells (overlaps)
    #------------------------
    X._applyBCOverlaps(t, depth=DEPTH,loc='centers',val=2, cellNName='cellN')
    C._initVars(t,'{centers:cellNChim}={centers:cellN}')

    #------------------------
    # Blanking IBM
    #------------------------
    C._initVars(t,'centers:cellN',1.)
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

    t = blankByIBCBodies(t, tb, 'centers', dimPb)
    C._initVars(t,'{centers:cellNIBC}={centers:cellN}')

    #-----------------------------------------
    # calcul de la normale et distance signee
    #-----------------------------------------
    COMPDIST = False # distance deja calculee ou non
    if Internal.getNodeFromName(t, 'TurbulentDistance') is None: COMPDIST=True
    if COMPDIST:
        print('Computing distance field...')
        DTW._distance2Walls(t, tb, loc='centers', type='ortho', signed=0)
    else: pass
    _signDistance(t)

    #-----------------------------------------
    # Pts IBC
    #-----------------------------------------
    C._initVars(t,'{centers:cellN}={centers:cellNIBC}')
    # determination des pts IBC
    Reynolds = Internal.getNodeFromName(tb, 'Reynolds')
    if Reynolds is not None: Reynolds = Internal.getValue(Reynolds)
    else: Reynolds = 6.e6
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
        # determination des pts IBC en fonction de la distance a la paroi
        # cree d'abord un front de type 1 pour assurer que l'adaptation ne descende pas en dessous de la limite...
        C._initVars(t,'{centers:cellNMin}={centers:cellNIBC}')
        if IBCType == -1: X._setHoleInterpolatedPoints(t,depth=-DEPTH,dir=0,loc='centers',cellNName='cellNMin',addGC=False)
        elif IBCType == 1: X._setHoleInterpolatedPoints(t,depth=1,dir=1,loc='centers',cellNName='cellNMin',addGC=False) # pour les gradients
        X._setHoleInterpolatedPoints(t,depth=DEPTH,dir=0,loc='centers',cellNName='cellNMin',addGC=False)

        for z in Internal.getZones(t):
            h = abs(C.getValue(z,'CoordinateX',0)-C.getValue(z,'CoordinateX',1))
            if yplus > 0.:
                height = G_IBM_Height.computeModelisationHeight(Re=Reynolds, yplus=yplus, L=Lref)
            else:
                height = G_IBM_Height.computeBestModelisationHeight(Re=Reynolds, h=h) # meilleur compromis entre hauteur entre le snear et la hauteur de modelisation
                yplus  = G_IBM_Height.computeYplus(Re=Reynolds, height=height, L=Lref)
            C._initVars(z,'{centers:cellN}=({centers:TurbulentDistance}>%20.16g)+(2*({centers:TurbulentDistance}<=%20.16g)*({centers:TurbulentDistance}>0))'%(height,height))

        # Si wallAdapt, utilisation de la solution precedente pour ne garder que les pts cibles tq y+PC <= y+ref : rapproche le front de la paroi (utile proche bord d'attaque)
        # Attention, le fichier d'adaptation doit etre un nuage de points dont les coordonnees correspondent aux points cibles (Cf modification dans extractIBMWallFields avec tb=None)
        if wallAdapt is not None:
            C._initVars(t,'{centers:yplus}=100000.')
            w = C.convertFile2PyTree(wallAdapt)
            total = len(Internal.getZones(t))
            cpt = 1
            for z in Internal.getZones(t):
                print("%d / %d"%(cpt, total))
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

        # Securite finale, on aura au min deux rangees de points cibles
        C._initVars(t, '{centers:cellN} = maximum({centers:cellN}, {centers:cellNMin})')
        C._rmVars(t,['centers:yplus', 'centers:cellNMin'])

        # permet de ne garder que les deux rangees superieures (prises en compte par le solveur), mais complique l'adaptation suivante et la visualisation
        if blankingF42: X._maximizeBlankedCells(t, depth=2, addGC=False)

    _removeBlankedGrids(t, loc='centers')

    print('Nb of Cartesian grids=%d for rank %d.'%(len(Internal.getZones(t)), Cmpi.rank))
    npts = 0
    for i in Internal.getZones(t):
        dims = Internal.getZoneDim(i)
        npts += dims[1]*dims[2]*dims[3]
    print('Final number of points=%5.4f millions for rank %d.'%(npts/1000000., Cmpi.rank))

    C._initVars(t,'{centers:cellNIBC}={centers:cellN}')

    #------------------------------------------------------------------------
    # Nature des points en fonction de leur nature Chimere et leur nature IBC
    #------------------------------------------------------------------------
    # -3 : agit comme un point masque - non donneur pour le type de point
    #  3  : agit comme donneur uniquement
    # updateNatureForIBM: modifie cellNChim, cellNFront, cellNIBM
    # cellNChim=-3, si cellNIBC=0 (masque)
    if IBCType == 1: # Points corriges IBM externes
        C._initVars(t,'{centers:cellNFront}=logical_and({centers:cellNIBC}>0.5, {centers:cellNIBC}<1.5)')
        for z in Internal.getZones(t):
            connector._updateNatureForIBM(z, IBCType,
                                          Internal.__GridCoordinates__,
                                          Internal.__FlowSolutionNodes__,
                                          Internal.__FlowSolutionCenters__)

    else: # EN 2 PARTIES : NECESSITE LE TRANSFERT DU FRONT PAR INTERPOLATION, QUI EST CALCULEE APRES
        print('Euler: on repousse le front un peu plus loin.')
        C._initVars(t,'{centers:dummy}={centers:cellN}') # sauvegarde
        C._initVars(t,'{centers:cellN}=({centers:cellNIBC}>0.5)*({centers:cellNIBC}<1.5)')
        X._setHoleInterpolatedPoints(t,depth=1,dir=1,loc='centers',cellNName='cellN',addGC=False)
        C._initVars(t,'{centers:cellNFront}=logical_and({centers:cellN}>0.5, {centers:cellN}<1.5)')
        C._cpVars(t,'centers:dummy',t,'centers:cellN')
        C._rmVars(t, ['centers:dummy'])
        for z in Internal.getZones(t):
            connector._updateNatureForIBM(z, IBCType,
                                          Internal.__GridCoordinates__,
                                          Internal.__FlowSolutionNodes__,
                                          Internal.__FlowSolutionCenters__)
    #------------------------------------------------------------------------
    # setInterpData - Chimere
    C._initVars(t,'{centers:cellN}=maximum(0.,{centers:cellNChim})')# vaut -3, 0, 1, 2 initialement

    # maillage donneur: on MET les pts IBC comme donneurs
    tc = C.node2Center(t)
    FSN = Internal.getNodesFromName(tc, Internal.__FlowSolutionNodes__)
    Internal._rmNodesByName(FSN,'cellNFront')
    Internal._rmNodesByName(FSN,'cellNIBC')
    Internal._rmNodesByName(FSN, "TurbulentDistance")

    tbb = G.BB(tc)

    # Creation du dictionnaire des ADT pour les raccords
    if interpDataType == 1:
        dictOfADT = {}
        for zdnr in Internal.getZones(tc):
            zdnrname = zdnr[0]
            if zdnrname not in dictOfADT:
                HOOKADT = C.createHook(zdnr, 'adt')
                dictOfADT[zdnrname] = HOOKADT
    else: dictOfADT = None
    print('Interpolations Chimere.')
    tc = doInterp(t, tc, tbb, tb=None, typeI='ID', dim=dimPb,
                  interpDataType=interpDataType, dictOfADT=dictOfADT)
    if dictOfADT is not None:
        for dnrname in dictOfADT: C.freeHook(dictOfADT[dnrname])

    # setIBCData - IBC
    C._initVars(t,'{centers:cellNIBCDnr}=minimum(2.,abs({centers:cellNIBC}))')
    C._initVars(t,'{centers:cellNIBC}=maximum(0.,{centers:cellNIBC})')# vaut -3, 0, 1, 2, 3 initialement
    C._initVars(t,'{centers:cellNIBC}={centers:cellNIBC}*({centers:cellNIBC}<2.5)')
    C._cpVars(t,'centers:cellNIBC',t,'centers:cellN')
    C._cpVars(t,'centers:cellN',tc,'cellN')

    #-----------------------------------------------
    # Transfert du cellNFront
    C._cpVars(t,'centers:cellNFront',tc,'cellNFront')

    for zc in Internal.getZones(tc):
        cellNFront = Internal.getNodeFromName2(zc,'cellNFront')
        if cellNFront != []:
            cellNFront = cellNFront[1]
            sizeTot = cellNFront.shape[0]*cellNFront.shape[1]*cellNFront.shape[2]
            sizeOne =  int(numpy.sum(cellNFront))
            if sizeOne < sizeTot:
                XOD._setInterpTransfers(t,zc,variables=['cellNFront'],cellNVariable='cellNFront',compact=0)

    if frontType == 2:
        _pushBackImageFront2(t, tc, tbb, interpDataType=interpDataType)
        if smoothing and dimPb==2: _smoothImageFront(t, tc); _pushBackImageFront2(t, tc, tbb, interpDataType=interpDataType)

        C._cpVars(t,'centers:cellNFront',tc,'cellNFront')

        for zc in Internal.getZones(tc):
            cellNFront = Internal.getNodeFromName2(zc,'cellNFront')
            if cellNFront != []:
                cellNFront = cellNFront[1]
                sizeTot = cellNFront.shape[0]*cellNFront.shape[1]*cellNFront.shape[2]
                sizeOne =  int(numpy.sum(cellNFront))
                if sizeOne < sizeTot:
                    XOD._setInterpTransfers(t,zc,variables=['cellNFront'],cellNVariable='cellNFront',compact=0)

    ## Fin traitement specifique, vaut 0 ou 1 apres la ligne suivante
    # C._cpVars(t,'centers:cellNFront',tc,'cellNFront')
    C._rmVars(t,['centers:cellNFront'])
    C._cpVars(t,'centers:TurbulentDistance',tc,'TurbulentDistance')

    print('Minimum distance: %f.'%C.getMinValue(t,'centers:TurbulentDistance'))
    P._computeGrad2(t, 'centers:TurbulentDistance', withCellN=False)
    print('Building the IBM front.')
    front = getIBMFront(tc, 'cellNFront', dimPb, frontType)
    if isPrintDebug:
        C.convertPyTree2File(front, 'IB_front.cgns')
    # C.convertPyTree2File(front, 'front.cgns')
    print('Interpolations IBM')
    tc = doInterp(t, tc, tbb, tb=tb,typeI='IBCD', dim=dimPb, dictOfADT=None, front=front, frontType=frontType, depth=DEPTH, IBCType=IBCType, interpDataType=interpDataType, Reynolds=Reynolds, yplus=yplus, Lref=Lref, isLBM=isLBM)

    # cleaning...
    Internal._rmNodesByName(tc, Internal.__FlowSolutionNodes__)
    if not Internal.getNodeFromName(t,'TimeMotion'):
        Internal._rmNodesByName(tc, Internal.__GridCoordinates__)
    C._initVars(t,'{centers:cellN}=minimum({centers:cellNChim}*{centers:cellNIBCDnr},2.)')
    varsRM = ['centers:gradxTurbulentDistance','centers:gradyTurbulentDistance','centers:gradzTurbulentDistance','centers:cellNFront','centers:cellNIBCDnr']
    varsRM += ['centers:cellNChim','centers:cellNIBC']
    C._rmVars(t, varsRM)
    C._rmVars(tc,['cellNChim','cellNIBC','TurbulentDistance'])
    #----------
    # SORTIE
    #----------
    return t, tc

#=============================================================================
# Version local explicit [GJ]
#=============================================================================
def doInterp2(t, tc, tbb, tb=None, typeI='ID', dim=3, dictOfADT=None, front=None, frontType=0, depth=2, IBCType=1, interpDataType=1, Reynolds=6.e6, yplus=100., Lref=1., check=False):
    ReferenceState = Internal.getNodeFromType2(t, 'ReferenceState_t')
    model = Internal.getNodeFromName(t, 'GoverningEquations')
    if model is not None: model = Internal.getValue(model)
    else:                 model = "Euler"

    bases  = Internal.getNodesFromType1(t     , 'CGNSBase_t')
    dimmm  = Internal.getNodeFromName2(bases[0], 'EquationDimension')
    dimPb   = Internal.getValue(dimmm)
    dxmax = 0.0


    zones = Internal.getZones(t)

    dico_dx = {}
    dico_dy = {}
    dico_dz = {}

    for z in zones:
        nodes = Internal.getNodesFromName(z, 'GridCoordinates')
        coordx = nodes[0][2][0][1]
        coordy = nodes[0][2][1][1]
        coordz = nodes[0][2][2][1]

        dxx  = abs(coordx[1,0,0]   - coordx[0,0,0])
        dyy  = abs(coordy[0,1,0]   - coordy[0,0,0])
        dzz  = abs(coordz[0,0,1]   - coordz[0,0,0])

        dico_dx[z[0]] = dxx
        dico_dy[z[0]] = dyy
        if dimPb == 2:dico_dz[z[0]] = 1
        else : dico_dz[z[0]] = dzz

        if dimPb == 2:dzz=max(dxx,dyy)

        dx = min(dxx,dyy,dzz)
        if dx > dxmax:dxmax=dx

    niveaux_temps = {}
    cx = {}

    for z in zones:
        nodes = Internal.getNodesFromName(z, 'GridCoordinates')
        coordx = nodes[0][2][0][1]
        coordy = nodes[0][2][1][1]
        coordz = nodes[0][2][2][1]

        dxx  = abs(coordx[1,0,0]   - coordx[0,0,0])
        dyy  = abs(coordy[0,1,0]   - coordy[0,0,0])
        dzz  = abs(coordz[0,0,1]   - coordz[0,0,0])

        if dimPb == 2:dzz=max(dxx,dyy)


        dx = min(dxx,dyy,dzz)

        #cx[z[0]]= coordx[1,0,0]

        N = math.log(dxmax/dx)/math.log(2.0)
        N = round(N) - 2
        if N < 0:
            niveaux_temps[z[0]] = 2**0
        else :
            niveaux_temps[z[0]] = 2**N
        ##if (N < 6):N=0
        ##else:N=1
        #if (cx[z[0]] < 0.98): niveaux_temps[z[0]] = 2**N
        #else:  niveaux_temps[z[0]] = 1
        #niveaux_temps['cart.264'] = 2
        #niveaux_temps['cart.294'] = 2

        #niveaux_temps[z[0]] = 2**N

        print(niveaux_temps[z[0]])
        #print(round(dxmax/dx))


    if typeI == 'ID':
        # toutes les zones sont interpolables en Chimere
        intersectionsDict = X.getIntersectingDomains(tbb, method='AABB', taabb=tbb)
        rcvZones = []
        for zrcv in Internal.getZones(t):
            if C.getMaxValue(zrcv,'centers:cellN')==2.:
                zrcvname = zrcv[0]; rcvZones.append(zrcv)
        nozr = 0
        nbZonesChim = len(rcvZones)
        for nozr in range(nbZonesChim):
            zrcv = rcvZones[nozr]
            zrcvname = zrcv[0]
            nozr += 1; hook0 = []
            nobOfDnrBases = []; nobOfDnrZones=[]; dnrZones=[]
            for nobd in range(len(tc[2])):
                if tc[2][nobd][3] == 'CGNSBase_t':
                    for nozd in range(len(tc[2][nobd][2])):
                        zdnr = tc[2][nobd][2][nozd]
                        if zdnr[3] == 'Zone_t':
                            zdnrname = zdnr[0]
                            if zdnrname in intersectionsDict[zrcvname]:
                                nobOfDnrBases.append(nobd)
                                nobOfDnrZones.append(nozd)
                                dnrZones.append(zdnr)
                                if interpDataType==1 and dictOfADT is not None:
                                    hook0.append(dictOfADT[zdnrname])
            if interpDataType == 0: hook0 = None

            X._setInterpData(zrcv, dnrZones, nature=1,penalty=1,loc='centers',storage='inverse',sameName=1,\
                             interpDataType=interpDataType, itype='chimera')


            levelrcv = niveaux_temps[zrcv[0]]


            for nod in range(len(dnrZones)):
                dim__ = Internal.getZoneDim(dnrZones[nod])
                prange = numpy.zeros(6,dtype=Internal.E_NpyInt)
                prangedonor = numpy.zeros(6,dtype=Internal.E_NpyInt)
                profondeur=numpy.zeros(1,dtype=Internal.E_NpyInt)
                dirD=numpy.zeros(1,dtype=Internal.E_NpyInt)
                dirR=numpy.zeros(1,dtype=Internal.E_NpyInt)

                plist = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][2][1]
                plistdnr = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][3][1]
                coeff = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][4][1]
                typ = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][5][1]

                leveldnr = niveaux_temps[dnrZones[nod][0]]

                nobd = nobOfDnrBases[nod]
                nozd = nobOfDnrZones[nod]
                tc[2][nobd][2][nozd] = dnrZones[nod]

                prangebis=numpy.reshape(prange,6)

                info = dnrZones[nod][2][len(dnrZones[nod][2])-1]
                info[2].append(['PointRange', prangebis , [], 'IndexArray_t'])

                transfo=numpy.zeros(3,dtype=Internal.E_NpyInt)#XOD.getTransfo(dnrZones[nod],zrcv)

                connector.indiceToCoord2(plist,prangedonor,transfo,profondeur,dirD,typ,dirR,plist.size,dim__[1],dim__[2],dim__[3])


                #connector.correctCoeffList(plist, coeff, typ, plist.size , dim__[1] , dim__[2] , dim__[3])

                NMratio = numpy.zeros(3,dtype=Internal.E_NpyInt)
                NMratio[0]=1
                NMratio[1]=1
                NMratio[2]=1

                info[2].append(['PointRangeDonor', prangedonor , [], 'IndexArray_t'])
                info[2].append(['DirReceveur', dirR , [], 'IndexArray_t'])
                info[2].append(['DirDonneur', dirD , [], 'IndexArray_t'])
                info[2].append(['Transform', transfo , [], 'IndexArray_t'])
                info[2].append(['Profondeur', profondeur , [], 'IndexArray_t'])
                info[2].append(['PointPivot', transfo , [], 'IndexArray_t'])
                info[2].append(['NMratio', NMratio , [], 'IndexArray_t'])
                info[2].append(['LevelZRcv', levelrcv , [], 'IndexArray_t'])
                info[2].append(['LevelZDnr', leveldnr , [], 'IndexArray_t'])



    elif typeI == 'IBCD':
        # detection des zones IBC
        zonesRIBC = []
        for zrcv in Internal.getZones(t):
            if C.getMaxValue(zrcv,'centers:cellNIBC')==2.:
                zrcvname = zrcv[0]; zonesRIBC.append(zrcv)

        if zonesRIBC == []: return tc

        res = getAllIBMPoints(zonesRIBC, loc='centers',tb=tb, tfront=front, frontType=frontType, \
                              cellNName='cellNIBC', depth=depth, IBCType=IBCType, Reynolds=Reynolds, yplus=yplus, Lref=Lref, check=check)
        nbZonesIBC = len(zonesRIBC)
        dictOfADT = {}
        dictOfCorrectedPtsByIBCType = res[0]
        dictOfWallPtsByIBCType = res[1]
        dictOfInterpPtsByIBCType = res[2]
        for ibcTypeL in  dictOfCorrectedPtsByIBCType:
            allCorrectedPts = dictOfCorrectedPtsByIBCType[ibcTypeL]
            allWallPts = dictOfWallPtsByIBCType[ibcTypeL]
            allInterpPts = dictOfInterpPtsByIBCType[ibcTypeL]
            for nozr in range(nbZonesIBC):
                if allCorrectedPts[nozr] != []:
                    interpPtsBB=Generator.BB(allInterpPts[nozr])
                    zrcv = zonesRIBC[nozr]
                    zrcvname = zrcv[0]
                    nobOfDnrBases = []; nobOfDnrZones=[]; dnrZones=[]
                    if interpDataType == 1: hook0 = []
                    else: hook0 = None
                    for nobd in range(len(tc[2])):
                        if tc[2][nobd][3] == 'CGNSBase_t':
                            for nozd in range(len(tc[2][nobd][2])):
                                zdnr = tc[2][nobd][2][nozd]
                                if zdnr[3] == 'Zone_t':
                                    zdnrname = zdnr[0]
                                    zbb = tbb[2][nobd][2][nozd]
                                    bba = C.getFields(Internal.__GridCoordinates__,zbb)[0]
                                    if Generator.bboxIntersection(interpPtsBB,bba,isBB=True) == 1:
                                        if interpDataType == 1:
                                            if zdnrname not in dictOfADT:
                                                HOOKADT = C.createHook(zdnr, 'adt')
                                                dictOfADT[zdnrname] = HOOKADT
                                            hook0.append(dictOfADT[zdnrname])

                                        dnrZones.append(zdnr)
                                        nobOfDnrBases.append(nobd)
                                        nobOfDnrZones.append(nozd)

                    XOD._setIBCDataForZone__(zrcv, dnrZones, allCorrectedPts[nozr], allWallPts[nozr], allInterpPts[nozr], \
                                             nature=1, penalty=1, loc='centers', storage='inverse',
                                             interpDataType=interpDataType, hook=hook0, dim=dim, \
                                             ReferenceState=ReferenceState, bcType=ibcTypeL,model=model)

                    nozr += 1

                    levelrcv = niveaux_temps[zrcv[0]]

                    for nod in range(len(dnrZones)):

                        dim__ = Internal.getZoneDim(dnrZones[nod])
                        prange = numpy.zeros(6,dtype=Internal.E_NpyInt)
                        prangedonor = numpy.zeros(6,dtype=Internal.E_NpyInt)
                        profondeur=numpy.zeros(1,dtype=Internal.E_NpyInt)
                        dirD=numpy.zeros(1,dtype=Internal.E_NpyInt)
                        dirR=numpy.zeros(1,dtype=Internal.E_NpyInt)

                        plist = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][2][1]
                        plistdnr = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][3][1]
                        coeff = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][4][1]
                        typ = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][5][1]

                        leveldnr = niveaux_temps[dnrZones[nod][0]]

                        nobd = nobOfDnrBases[nod]
                        nozd = nobOfDnrZones[nod]

                        tc[2][nobd][2][nozd] = dnrZones[nod]


                        prangebis=numpy.reshape(prange,6)

                        info = dnrZones[nod][2][len(dnrZones[nod][2])-1]
                        info[2].append(['PointRange', prangebis , [], 'IndexArray_t'])

                        transfo=numpy.zeros(3,dtype=Internal.E_NpyInt)#XOD.getTransfo(dnrZones[nod],zrcv)

                        connector.indiceToCoord2(plist,prangedonor,transfo,profondeur,dirD,typ,dirR,plist.size,dim__[1],dim__[2],dim__[3])

                        #connector.correctCoeffList(plist, coeff, typ, plist.size , dim__[1] , dim__[2] , dim__[3])

                        NMratio = numpy.zeros(3,dtype=Internal.E_NpyInt)
                        NMratio[0]=1
                        NMratio[1]=1
                        NMratio[2]=1

                        info[2].append(['PointRangeDonor', prangedonor , [], 'IndexArray_t'])
                        info[2].append(['DirReceveur', dirR , [], 'IndexArray_t'])
                        info[2].append(['DirDonneur', dirD , [], 'IndexArray_t'])
                        info[2].append(['Transform', transfo , [], 'IndexArray_t'])
                        info[2].append(['Profondeur', profondeur , [], 'IndexArray_t'])
                        info[2].append(['PointPivot', transfo , [], 'IndexArray_t'])
                        info[2].append(['NMratio', NMratio , [], 'IndexArray_t'])
                        info[2].append(['LevelZRcv', levelrcv , [], 'IndexArray_t'])
                        info[2].append(['LevelZDnr', leveldnr , [], 'IndexArray_t'])


        if dictOfADT is not None:
            for dnrname in dictOfADT: C.freeHook(dictOfADT[dnrname])

    return tc

def prepareIBMData2(t, tbody, DEPTH=2, loc='centers', frontType=1, inv=False, interpDataType=1):
    tb =  Internal.copyRef(tbody)

    # tb: fournit model et dimension
    dimPb = Internal.getNodeFromName(tb,'EquationDimension')
    if dimPb is None: raise ValueError('prepareIBMData: EquationDimension is missing in input body tree.')
    dimPb = Internal.getValue(dimPb)

    # type de traitement paroi: pts interieurs ou externes
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('prepareIBMData: GoverningEquations is missing in input body tree.')
    # model: Euler, NSLaminar, NSTurbulent
    model = Internal.getValue(model)

    if model == 'Euler': IBCType =-1
    else: IBCType = 1 # Points cibles externes
    if loc == 'nodes':
        raise NotImplemented("prepareIBMData: prepareIBMData at nodes not yet implemented.")

    #------------------------
    # Ghost cells (overlaps)
    #------------------------
    X._applyBCOverlaps(t, depth=DEPTH,loc='centers',val=2, cellNName='cellN')
    C._initVars(t,'{centers:cellNChim}={centers:cellN}')

    #------------------------
    # Blanking IBM
    #------------------------
    C._initVars(t,'centers:cellN',1.)
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

    t = blankByIBCBodies(t,tb,'centers',dimPb)
    if not inv: C._initVars(t,'{centers:cellNIBC}={centers:cellN}')
    if inv: C._initVars(t,'{centers:cellNIBC}=1-{centers:cellN}') # ecoulement interne



    #-----------------------------------------
    # calcul de la normale et distance signee
    #-----------------------------------------
    COMPDIST = False # distance deja calculee ou non
    if Internal.getNodeFromName(t, 'TurbulentDistance') is None: COMPDIST=True
    if COMPDIST:
        print('Computing distance field...')
        DTW._distance2Walls(t,tb,loc='centers',type='ortho',signed=0)
    else: pass
    _signDistance(t)

    #-----------------------------------------
    # Pts IBC
    #-----------------------------------------
    C._initVars(t,'{centers:cellN}={centers:cellNIBC}')
    # determination des pts IBC
    if IBCType == -1: X._setHoleInterpolatedPoints(t,depth=-DEPTH,dir=0,loc='centers',cellNName='cellN',addGC=False)
    elif IBCType == 1:
        X._setHoleInterpolatedPoints(t,depth=1,dir=1,loc='centers',cellNName='cellN',addGC=False) # pour les gradients
        if frontType < 2:
            X._setHoleInterpolatedPoints(t,depth=DEPTH,dir=0,loc='centers',cellNName='cellN',addGC=False)
        else:
            DEPTHL=DEPTH+1
            X._setHoleInterpolatedPoints(t,depth=DEPTHL,dir=0, loc='centers',cellNName='cellN',addGC=False)
            #cree des pts extrapoles supplementaires
            # _blankClosestTargetCells(t,cellNName='cellN', depth=DEPTHL)
    else:
        raise ValueError('prepareIBMData: not valid IBCType. Check model.')
    _removeBlankedGrids(t, loc='centers')
    print('Nb of Cartesian grids=%d.'%len(Internal.getZones(t)))
    npts = 0
    for i in Internal.getZones(t):
        dims = Internal.getZoneDim(i)
        npts += dims[1]*dims[2]*dims[3]
    print('Final number of points=%5.4f millions.'%(npts/1000000.))

    C._initVars(t,'{centers:cellNIBC}={centers:cellN}')

    #------------------------------------------------------------------------
    # Nature des points en fonction de leur nature Chimere et leur nature IBC
    #------------------------------------------------------------------------
    # -3 : agit comme un point masque - non donneur pour le type de point
    #  3  : agit comme donneur uniquement
    # updateNatureForIBM: modifie cellNChim, cellNFront, cellNIBM
    # cellNChim=-3, si cellNIBC=0 (masque)
    if IBCType == 1: # Points corriges IBM externes
        C._initVars(t,'{centers:cellNFront}=logical_and({centers:cellNIBC}>0.5, {centers:cellNIBC}<1.5)')
        for z in Internal.getZones(t):
            connector._updateNatureForIBM(z, IBCType,
                                          Internal.__GridCoordinates__,
                                          Internal.__FlowSolutionNodes__,
                                          Internal.__FlowSolutionCenters__)

    else: # EN 2 PARTIES : NECESSITE LE TRANSFERT DU FRONT PAR INTERPOLATION, QUI EST CALCULEE APRES
        print('Euler: on repousse le front un peu plus loin.')
        C._initVars(t,'{centers:dummy}={centers:cellN}') # sauvegarde
        C._initVars(t,'{centers:cellN}=({centers:cellNIBC}>0.5)*({centers:cellNIBC}<1.5)')
        X._setHoleInterpolatedPoints(t,depth=1,dir=1,loc='centers',cellNName='cellN',addGC=False)
        C._initVars(t,'{centers:cellNFront}=logical_and({centers:cellN}>0.5, {centers:cellN}<1.5)')
        C._cpVars(t,'centers:dummy',t,'centers:cellN')
        C._rmVars(t, ['centers:dummy'])
        for z in Internal.getZones(t):
            connector._updateNatureForIBM(z, IBCType,
                                          Internal.__GridCoordinates__,
                                          Internal.__FlowSolutionNodes__,
                                          Internal.__FlowSolutionCenters__)
    #------------------------------------------------------------------------
    # setInterpData - Chimere
    C._initVars(t,'{centers:cellN}=maximum(0.,{centers:cellNChim})')# vaut -3, 0, 1, 2 initialement

    # maillage donneur: on MET les pts IBC comme donneurs
    tc = C.node2Center(t)
    FSN = Internal.getNodesFromName(tc, Internal.__FlowSolutionNodes__)
    Internal._rmNodesByName(FSN,'cellNFront')
    Internal._rmNodesByName(FSN,'cellNIBC')
    Internal._rmNodesByName(FSN, "TurbulentDistance")

    tbb = G.BB(tc)

    # Creation du dictionnaire des ADT pour les raccords
    if interpDataType == 1:
        dictOfADT = {}
        for zdnr in Internal.getZones(tc):
            zdnrname = zdnr[0]
            if zdnrname not in dictOfADT:
                HOOKADT = C.createHook(zdnr, 'adt')
                dictOfADT[zdnrname] = HOOKADT
    else: dictOfADT = None
    print('Interpolations Chimere.')
    tc = doInterp2(t, tc, tbb, tb=None, typeI='ID', dim=dimPb,
                   interpDataType=interpDataType, dictOfADT=dictOfADT)
    if dictOfADT is not None:
        for dnrname in dictOfADT: C.freeHook(dictOfADT[dnrname])

    # setIBCData - IBC
    C._initVars(t,'{centers:cellNIBCDnr}=minimum(2.,abs({centers:cellNIBC}))')
    C._initVars(t,'{centers:cellNIBC}=maximum(0.,{centers:cellNIBC})')# vaut -3, 0, 1, 2, 3 initialement
    C._initVars(t,'{centers:cellNIBC}={centers:cellNIBC}*({centers:cellNIBC}<2.5)')
    C._cpVars(t,'centers:cellNIBC',t,'centers:cellN')
    C._cpVars(t,'centers:cellN',tc,'cellN')

    #-----------------------------------------------
    # Transfert du cellNFront
    C._cpVars(t,'centers:cellNFront',tc,'cellNFront')

    for zc in Internal.getZones(tc):
        cellNFront = Internal.getNodeFromName2(zc,'cellNFront')
        if cellNFront != []:
            cellNFront = cellNFront[1]
            sizeTot = cellNFront.shape[0]*cellNFront.shape[1]*cellNFront.shape[2]
            sizeOne =  int(numpy.sum(cellNFront))
            if sizeOne < sizeTot:
                XOD._setInterpTransfers(t,zc,variables=['cellNFront'],cellNVariable='cellNFront',compact=0)

    if frontType==2 or frontType==3: _pushBackImageFront2(t, tc, tbb, interpDataType=interpDataType)

    ## Fin traitement specifique, vaut 0 ou 1 apres la ligne suivante
    C._cpVars(t,'centers:cellNFront',tc,'cellNFront')
    C._rmVars(t,['centers:cellNFront'])
    C._cpVars(t,'centers:TurbulentDistance',tc,'TurbulentDistance')

    print('Minimum distance: %f.'%C.getMinValue(t,'centers:TurbulentDistance'))
    P._computeGrad2(t, 'centers:TurbulentDistance', withCellN=False)
    print('Building the IBM front.')
    front = getIBMFront(tc, 'cellNFront', dimPb, frontType)
    print('Interpolations IBM')
    tc = doInterp2(t,tc,tbb, tb=tb,typeI='IBCD',dim=dimPb, dictOfADT=None, front=front, frontType=frontType, depth=DEPTH, IBCType=IBCType, interpDataType=interpDataType)

    # cleaning...
    Internal._rmNodesByName(tc, Internal.__FlowSolutionNodes__)
    Internal._rmNodesByName(tc, Internal.__GridCoordinates__)
    C._initVars(t,'{centers:cellN}=minimum({centers:cellNChim}*{centers:cellNIBCDnr},2.)')
    varsRM = ['centers:gradxTurbulentDistance','centers:gradyTurbulentDistance','centers:gradzTurbulentDistance','centers:cellNFront','centers:cellNIBCDnr']
    varsRM += ['centers:cellNChim','centers:cellNIBC']
    C._rmVars(t, varsRM)
    C._rmVars(tc,['cellNChim','cellNIBC','TurbulentDistance'])
    #----------
    # SORTIE
    #----------
    return t, tc

#=============================================================================
# never called ?
#=============================================================================
def doInterp3(t, tc, tbb, tb=None, typeI='ID', dim=3, dictOfADT=None, frontType=0, depth=2, IBCType=1, interpDataType=1, Reynolds=6.e6, yplus=100., Lref=1., check=False):


    ReferenceState = Internal.getNodeFromType2(t, 'ReferenceState_t')
    model = Internal.getNodeFromName(t, 'GoverningEquations')
    if model is not None: model = Internal.getValue(model)
    else:                 model = "Euler"

    bases  = Internal.getNodesFromType1(t     , 'CGNSBase_t')
    dimmm  = Internal.getNodeFromName2(bases[0], 'EquationDimension')
    dimPb   = Internal.getValue(dimmm)
    dxmax = 0.0


    zones = Internal.getZones(t)

    dico_dx = {}
    dico_dy = {}
    dico_dz = {}

    for z in zones:
        nodes = Internal.getNodesFromName(z, 'GridCoordinates')
        coordx = nodes[0][2][0][1]
        coordy = nodes[0][2][1][1]
        coordz = nodes[0][2][2][1]

        dxx  = abs(coordx[1,0,0]   - coordx[0,0,0])
        dyy  = abs(coordy[0,1,0]   - coordy[0,0,0])
        dzz  = abs(coordz[0,0,1]   - coordz[0,0,0])

        dico_dx[z[0]] = dxx
        dico_dy[z[0]] = dyy
        if dimPb == 2: dico_dz[z[0]] = 1
        else: dico_dz[z[0]] = dzz

        if dimPb == 2:dzz=max(dxx,dyy)

        dx = min(dxx,dyy,dzz)
        if dx > dxmax:dxmax=dx

    niveaux_temps = {}
    cx = {}

    for z in zones:
        nodes = Internal.getNodesFromName(z, 'GridCoordinates')
        coordx = nodes[0][2][0][1]
        coordy = nodes[0][2][1][1]
        coordz = nodes[0][2][2][1]

        dxx  = abs(coordx[1,0,0]   - coordx[0,0,0])
        dyy  = abs(coordy[0,1,0]   - coordy[0,0,0])
        dzz  = abs(coordz[0,0,1]   - coordz[0,0,0])

        if dimPb == 2:dzz=max(dxx,dyy)

        dx = min(dxx,dyy,dzz)

        #cx[z[0]]= coordx[1,0,0]

        N = math.log(dxmax/dx)/math.log(2.0)
        N = round(N)
        if N < 6: N=0
        else: N=1
        #if (cx[z[0]] < 0.98): niveaux_temps[z[0]] = 2**N
        #else:  niveaux_temps[z[0]] = 1
        #niveaux_temps['cart.264'] = 2
        #niveaux_temps['cart.294'] = 2

        niveaux_temps[z[0]] = 2**N


        print(niveaux_temps[z[0]])
        #print(round(dxmax/dx))


    if typeI == 'ID':
        # toutes les zones sont interpolables en Chimere
        intersectionsDict = X.getIntersectingDomains(tbb, method='AABB', taabb=tbb)


        rcvZones = []
        for zrcv in Internal.getZones(t):
            if C.getMaxValue(zrcv,'centers:cellN')==2.:
                zrcvname = zrcv[0]; rcvZones.append(zrcv)

        #dico={}
        #for zrcv in Internal.getZones(t):
               # listofjoins = Internal.getNodesFromType2(zrcv, 'GridConnectivity_t')
               # if listofjoins is not None:
                #    prange_list=[]
                #    dir_list=[]
                #    for join in listofjoins:
                #        prange_ = Internal.getNodeFromName1(join,'PointRange')[1]
                #dirR = CGC.getDirBorderStruct__(prange_,dimPb)
                #dir_list.append(dirR)
                #print 'prange_= ', prange_
                #       for i in range(3):
                #           if prange_[i,1] == prange_[i,0] and prange_[i,1] != 1:
                #               prange_[i,1] =  prange_[i,1]-1
                #               prange_[i,0] =  prange_[i,0]-1
                #           elif prange_[i,1] != prange_[i,0] and prange_[i,1] != 1 :
                #              prange_[i,1] =  prange_[i,1]-1
                #      prange_=numpy.reshape(prange_,6)
                #      prange_list.append(prange_)
                #  dico[zrcv[0]]=prange_list
                #dico[zrcv[0]]=dir_list
                # print prange_, zrcv[0]


        nozr = 0
        nbZonesChim = len(rcvZones)
        for nozr in range(nbZonesChim):
            zrcv = rcvZones[nozr]
            dim_ = Internal.getZoneDim(zrcv)
            zrcvname = zrcv[0]
            nozr += 1; hook0 = []
            nobOfDnrBases = []; nobOfDnrZones=[]; dnrZones=[]
            for nobd in range(len(tc[2])):
                if tc[2][nobd][3] == 'CGNSBase_t':
                    for nozd in range(len(tc[2][nobd][2])):
                        zdnr = tc[2][nobd][2][nozd]
                        if zdnr[3] == 'Zone_t':
                            zdnrname = zdnr[0]
                            if zdnrname in intersectionsDict[zrcvname]:
                                nobOfDnrBases.append(nobd)
                                nobOfDnrZones.append(nozd)
                                dnrZones.append(zdnr)
                                hook0.append(dictOfADT[zdnrname])

            dnrZones = X.setInterpData(zrcv,dnrZones,nature=1,penalty=1,loc='centers',storage='inverse',sameName=1,\
                                       hook=hook0, itype='chimera')


            levelrcv = niveaux_temps[zrcv[0]]

            for nod in range(len(dnrZones)):

                dim__ = Internal.getZoneDim(dnrZones[nod])
                prange = numpy.zeros(6,dtype=Internal.E_NpyInt)
                prangedonor = numpy.zeros(6,dtype=Internal.E_NpyInt)
                profondeur=numpy.zeros(1,dtype=Internal.E_NpyInt)
                dirD=numpy.zeros(1,dtype=Internal.E_NpyInt)
                dirR=numpy.zeros(1,dtype=Internal.E_NpyInt)

                plist = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][2][1]
                plistdnr = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][3][1]
                coeff = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][4][1]
                typ = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][5][1]

                leveldnr = niveaux_temps[dnrZones[nod][0]]

                nobd = nobOfDnrBases[nod]
                nozd = nobOfDnrZones[nod]

                tc[2][nobd][2][nozd] = dnrZones[nod]


                prangebis=numpy.reshape(prange,6)

                info = dnrZones[nod][2][len(dnrZones[nod][2])-1]
                info[2].append(['PointRange', prangebis , [], 'IndexArray_t'])

                transfo = XOD.getTransfo(dnrZones[nod],zrcv)

                connector.indiceToCoord2(plist,prangedonor,transfo,profondeur,dirD,typ,dirR,plist.size,dim__[1],dim__[2],dim__[3])


                #connector.correctCoeffList(plist, coeff, typ, plist.size , dim__[1] , dim__[2] , dim__[3])

                NMratio = numpy.zeros(3,dtype=Internal.E_NpyInt)
                NMratio[0]=1
                NMratio[1]=1
                NMratio[2]=1

                info[2].append(['PointRangeDonor', prangedonor , [], 'IndexArray_t'])
                info[2].append(['DirReceveur', dirR , [], 'IndexArray_t'])
                info[2].append(['DirDonneur', dirD , [], 'IndexArray_t'])
                info[2].append(['Transform', transfo , [], 'IndexArray_t'])
                info[2].append(['Profondeur', profondeur , [], 'IndexArray_t'])
                info[2].append(['PointPivot', transfo , [], 'IndexArray_t'])
                info[2].append(['NMratio', NMratio , [], 'IndexArray_t'])
                info[2].append(['LevelZRcv', levelrcv , [], 'IndexArray_t'])
                info[2].append(['LevelZDnr', leveldnr , [], 'IndexArray_t'])



    elif typeI == 'IBCD':
        # detection des zones IBC
        zonesRIBC = []
        for zrcv in Internal.getZones(t):
            if C.getMaxValue(zrcv,'centers:cellNIBC')==2.:
                zrcvname = zrcv[0]; zonesRIBC.append(zrcv)

        if zonesRIBC == []: return tc

        print('Building the IBM front.')

        front = getIBMFront(tc, 'cellNFront', dim, frontType)
        # Sortie du front pour debug
        C.convertPyTree2File(front, 'front.cgns')

        res = getAllIBMPoints(zonesRIBC, loc='centers',tb=tb, tfront=front, frontType=frontType, \
                              cellNName='cellNIBC', depth=depth, IBCType=IBCType, Reynolds=Reynolds, yplus=yplus, Lref=Lref, check=check)
        nbZonesIBC = len(zonesRIBC)
        dictOfADT = {}
        dictOfCorrectedPtsByIBCType = res[0]
        dictOfWallPtsByIBCType = res[1]
        dictOfInterpPtsByIBCType = res[2]
        for ibcTypeL in dictOfCorrectedPtsByIBCType:
            allCorrectedPts = dictOfCorrectedPtsByIBCType[ibcTypeL]
            allWallPts = dictOfWallPtsByIBCType[ibcTypeL]
            allInterpPts = dictOfInterpPtsByIBCType[ibcTypeL]
            for nozr in range(nbZonesIBC):
                if allCorrectedPts[nozr] != []:
                    interpPtsBB=Generator.BB(allInterpPts[nozr])

                    zrcv = zonesRIBC[nozr]
                    zrcvname = zrcv[0]
                    nobOfDnrBases = []; nobOfDnrZones=[]; dnrZones=[]; hook0 = []
                    for nobd in range(len(tc[2])):
                        if tc[2][nobd][3] == 'CGNSBase_t':
                            for nozd in range(len(tc[2][nobd][2])):
                                zdnr = tc[2][nobd][2][nozd]
                                if zdnr[3] == 'Zone_t':
                                    zdnrname = zdnr[0]
                                    zbb = tbb[2][nobd][2][nozd]
                                    bba = C.getFields(Internal.__GridCoordinates__,zbb)[0]
                                    if Generator.bboxIntersection(interpPtsBB,bba,isBB=True) == 1:
                                        if zdnrname not in dictOfADT:
                                            HOOKADT = C.createHook(zdnr, 'adt')
                                            dictOfADT[zdnrname] = HOOKADT
                                        dnrZones.append(zdnr)
                                        hook0.append(dictOfADT[zdnrname])
                                        nobOfDnrBases.append(nobd)
                                        nobOfDnrZones.append(nozd)

                    XOD._setIBCDataForZone__(zrcv, dnrZones, allCorrectedPts[nozr], allWallPts[nozr], allInterpPts[nozr], \
                                             loc='centers', storage='inverse',  hook=hook0, dim=dim, ReferenceState=ReferenceState, bcType=ibcTypeL,model=model)
                    nozr += 1

                    levelrcv = niveaux_temps[zrcv[0]]

                    for nod in range(len(dnrZones)):

                        dim__ = Internal.getZoneDim(dnrZones[nod])
                        prange = numpy.zeros(6,dtype=Internal.E_NpyInt)
                        prangedonor = numpy.zeros(6,dtype=Internal.E_NpyInt)
                        profondeur=numpy.zeros(1,dtype=Internal.E_NpyInt)
                        dirD=numpy.zeros(1,dtype=Internal.E_NpyInt)
                        dirR=numpy.zeros(1,dtype=Internal.E_NpyInt)

                        plist = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][2][1]
                        plistdnr = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][3][1]
                        coeff = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][4][1]
                        typ = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][5][1]

                        leveldnr = niveaux_temps[dnrZones[nod][0]]

                        nobd = nobOfDnrBases[nod]
                        nozd = nobOfDnrZones[nod]

                        tc[2][nobd][2][nozd] = dnrZones[nod]


                        prangebis=numpy.reshape(prange,6)

                        info = dnrZones[nod][2][len(dnrZones[nod][2])-1]
                        info[2].append(['PointRange', prangebis , [], 'IndexArray_t'])

                        transfo=XOD.getTransfo(dnrZones[nod],zrcv)

                        connector.indiceToCoord2(plist,prangedonor,transfo,profondeur,dirD,typ,dirR,plist.size,dim__[1],dim__[2],dim__[3])

                        #connector.correctCoeffList(plist, coeff, typ, plist.size , dim__[1], dim__[2], dim__[3])

                        NMratio = numpy.zeros(3,dtype=Internal.E_NpyInt)
                        NMratio[0]=1
                        NMratio[1]=1
                        NMratio[2]=1

                        info[2].append(['PointRangeDonor', prangedonor , [], 'IndexArray_t'])
                        info[2].append(['DirReceveur', dirR , [], 'IndexArray_t'])
                        info[2].append(['DirDonneur', dirD , [], 'IndexArray_t'])
                        info[2].append(['Transform', transfo , [], 'IndexArray_t'])
                        info[2].append(['Profondeur', profondeur , [], 'IndexArray_t'])
                        info[2].append(['PointPivot', transfo , [], 'IndexArray_t'])
                        info[2].append(['NMratio', NMratio , [], 'IndexArray_t'])
                        info[2].append(['LevelZRcv', levelrcv , [], 'IndexArray_t'])
                        info[2].append(['LevelZDnr', leveldnr , [], 'IndexArray_t'])

                        print('LEVELS= ', levelrcv, leveldnr)


        for dnrname in dictOfADT.keys(): C.freeHook(dictOfADT[dnrname])

        for dnrname in dictOfADT: C.freeHook(dictOfADT[dnrname])

    return tc


def _setInjOutlet__(tc, tb):
    DicInj  = {}; DicOutP = {}
    for z in Internal.getZones(tb):
        sol = Internal.getNodeFromName(z, '.Solver#define')
        if sol:
            ibctype = Internal.getValue(Internal.getNodeFromName(sol, 'ibctype'))
            if ibctype == 'outpress':
                famName = Internal.getValue(Internal.getNodeFromName(z, 'FamilyName'))
                if Internal.getNodeFromName(sol, 'pStatic'):
                    pStatic  = Internal.getValue(Internal.getNodeFromName(sol, 'pStatic'));
                    isDensity = Internal.getValue(Internal.getNodeFromName(sol, 'isDensityConstant'));
                    DicOutP[famName] = str(pStatic)+'|'+str(isDensity)
            if ibctype == 'inj':
                famName = Internal.getValue(Internal.getNodeFromName(z, 'FamilyName'));
                if Internal.getNodeFromName(sol,'StagnationPressure'):
                    pStag    = Internal.getValue(Internal.getNodeFromName(sol, 'StagnationPressure'));
                    hStag    = Internal.getValue(Internal.getNodeFromName(sol, 'StagnationEnthalpy'));
                    dirx     = Internal.getValue(Internal.getNodeFromName(sol, 'dirx'));
                    diry     = Internal.getValue(Internal.getNodeFromName(sol, 'diry'));
                    dirz     = Internal.getValue(Internal.getNodeFromName(sol, 'dirz'));
                    DicInj[famName] = str(pStag)+'|'+str(hStag)+'|'+str(dirx)+'|'+str(diry)+'|'+str(dirz)

    for dic in DicOutP:
        PStatic   = float(DicOutP[dic].split('|')[0])
        isDensity = float(DicOutP[dic].split('|')[-1])
        isDensityConstant=False
        if isDensity>0: isDensityConstant=True
        D_IBM._initOutflow(tc, dic, PStatic, InterpolPlane=None, PressureVar=0, isDensityConstant=isDensityConstant)

    for dic in DicInj:
        PTot = float(DicInj[dic].split('|')[0])
        HTot = float(DicInj[dic].split('|')[1])
        dirx = float(DicInj[dic].split('|')[2])
        diry = float(DicInj[dic].split('|')[3])
        dirz = float(DicInj[dic].split('|')[4])
        D_IBM._initInj(tc, dic, PTot, HTot, injDir=[dirx,diry,dirz], InterpolPlane=None, PressureVar=0, EnthalpyVar=0)

    return None

## Separate IBM points that have a type 3 and type 4 projection
def _prepOutputProject__(outputProjection, typeValue, arrayLocal, allCorrectedPts, allWallPts, allInterpPts):
    indicesSave = numpy.argwhere(arrayLocal==typeValue)
    for i in indicesSave:
        outputProjection[0].append(allCorrectedPts[0][i][0])
        outputProjection[1].append(allCorrectedPts[1][i][0])
        outputProjection[2].append(allCorrectedPts[2][i][0])

        outputProjection[3].append(allWallPts[0][i][0])
        outputProjection[4].append(allWallPts[1][i][0])
        outputProjection[5].append(allWallPts[2][i][0])

        outputProjection[6].append(allInterpPts[0][i][0])
        outputProjection[7].append(allInterpPts[1][i][0])
        outputProjection[8].append(allInterpPts[2][i][0])
    return None

## Write CGNS file with the IBM points that have a type 3 and type 4 projection
def _writeOutputProject__(outputProjection, tLocal):
    nameZone = ['IBM', 'Wall', 'Image']
    for i in range(3):
        size     = len(outputProjection[i*3])
        coordxPC = numpy.array(outputProjection[i*3  ])
        coordyPC = numpy.array(outputProjection[i*3+1])
        coordzPC = numpy.array(outputProjection[i*3+2])
        zone = G.cart((0,0,0),(1,1,1),(size,1,1))
        zone[0] = nameZone[i]+"%"+str(Cmpi.rank)
        cont = Internal.getNodeFromName1(zone, Internal.__GridCoordinates__)
        XPC0 = Internal.getNodeFromName1(cont, 'CoordinateX')
        XPC0[1] = coordxPC
        YPC0 = Internal.getNodeFromName1(cont, 'CoordinateY')
        YPC0[1] = coordyPC
        ZPC0 = Internal.getNodeFromName1(cont, 'CoordinateZ')
        ZPC0[1] = coordzPC
        Internal.addChild(Internal.getNodeFromName(tLocal,nameZone[i]), zone, pos=-1) # at the end
    tLocal = C.convertArray2Node(tLocal)
    Cmpi._setProc(tLocal, Cmpi.rank)
    return tLocal



## construction info pour raccord nearmatch conservatif
def _buildConservativeFlux(t, tc, verbose=0):


    bases = Internal.getNodesFromType1(t, 'CGNSBase_t')
    for b in bases:
        model = "Euler"
        a = Internal.getNodeFromName2(b, 'GoverningEquations')
        if a is not None: model = Internal.getValue(a)
        neq = 5
        if model == 'nsspalart' or model =='NSTurbulent': neq = 6

    ## determine dx=dy for each zone & store per zone
    levelZone={}
    deltaZ={}
    zones = Internal.getZones(t)
    level0 = Internal.getNodeFromName(zones[0],'niveaux_temps')
    hmin_loc=1e30
    for z in zones:

        coordz        = Internal.getNodeFromName(z,'CoordinateZ')[1]
        dz            = coordz[0,0,1]-coordz[0,0,0]
        deltaZ[ z[0] ]= dz
        level = Internal.getNodeFromName(z,'niveaux_temps')

        if level is not None:
            levelZone[ z[0] ] = level[1][0]
        else:
            h              = abs(C.getValue(z,'CoordinateX',0)-C.getValue(z,'CoordinateX',1))
            #print('h_loc=', h, z[0], flush=True)
            levelZone[z[0]]= h
            if h < hmin_loc : hmin_loc = h

    #print('hmin_loc=', hmin_loc, flush=True)
    hmin = hmin_loc
    if Cmpi.size > 1 :
        hmin_loc = Cmpi.allgather(hmin_loc)
        hmin=1e30
        for h in hmin_loc:
            if h < hmin : hmin = h

    #print("hminGlob", hmin, flush=True)

    ## go from dx to dx/dx_min
    if level0 is None:
        Nlevels=1
        for i in levelZone:
            levelZone[i]= int( math.log( int(levelZone[i]/hmin + 0.00000001)  , 2) )
            if levelZone[i] +1  > Nlevels : Nlevels = int(levelZone[i]) +1

    #construction arbre skeleton global (tcs) pour calcul graph
    if Cmpi.size > 1:
        tcs_local = Cmpi.convert2SkeletonTree(tc)
        tcs       = Cmpi.allgatherTree(tcs_local)
        procDict  = Cmpi.getProcDict(tcs)
        graph     = Cmpi.computeGraph(tcs, type='ID', reduction=True, procDict=procDict)
        rank      = Cmpi.rank
        ## partage des info level delta_z en mpi
        levelZone = Cmpi.allgatherDict2(levelZone)
        deltaZ    = Cmpi.allgatherDict2(deltaZ)
    else:
        graph=None
        rank = 0
        procDict={}
        for z in Internal.getZones(t):
            procDict[ z[0] ]=0

    #construction list pour envoie fenetre zone distante
    ratio={}
    datas = {}
    dimR_loc={}
    for z in Internal.getZones(tc):
        dimR_loc[ z[0] ] = Internal.getZoneDim(z) # taille en centre

    if Cmpi.size > 1: dimR_loc = Cmpi.allgatherDict2(dimR_loc)

    for z in Internal.getZones(tc):

        zd_t  = Internal.getNodeFromName(t, z[0])
        levelD = levelZone[z[0]]
        subRegions =  Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        for s in subRegions:
            zRname = Internal.getValue(s)
            levelR = levelZone[zRname]
            if levelR > levelD:

                #print("raccord conservatif: levelRD", levelR , levelD,'zRD', zRname, z[0])
                proc   = procDict[zRname]
                dimD   = Internal.getZoneDim(z)  # taille en centre
                dimR   = dimR_loc[ zRname ]
                dimPb  = dimR[4]

                ratio_k=1
                if dimPb !=2:
                    dz_R = deltaZ[ zRname ]
                    dz_D = deltaZ[ z[0] ]
                    if dz_R/dz_D < 0.9999 or dz_R/dz_D > 1.0001: ratio_k=2

                ratio[z[0]]=[2,2,ratio_k]

                sh =[ dimR[1],dimR[2],dimR[3] ]
                shD=[ dimD[1],dimD[2],dimD[3] ]
                #print("dimR", dimR, "dimD", dimD)
                ptList = Internal.getNodeFromName1(s, 'PointList')[1]
                ptListD= Internal.getNodeFromName1(s, 'PointListDonor')[1]
                lmin = numpy.amin(ptListD)
                kmin =  lmin//(sh[0]*sh[1])
                jmin =  (lmin -kmin*sh[0]*sh[1])//sh[0]
                imin =  lmin -kmin*sh[0]*sh[1] -jmin*sh[0]
                lmax = numpy.amax(ptListD)
                kmax =  lmax//(sh[0]*sh[1])
                jmax =  (lmax -kmax*sh[0]*sh[1])//sh[0]
                imax =  lmax -kmax*sh[0]*sh[1] -jmax*sh[0]
                #imax,imin,...: adresse C


                win     = numpy.empty( (6,6), Internal.E_NpyInt)
                winD    = numpy.empty( (6,6), Internal.E_NpyInt)
                win[0,:]=100000 ;win[2,:]=100000 ;win[4,:]=100000
                win[1,:]=-1 ;win[3,:]=-1 ;win[5,:]=-1


                s1 = max( dimR[1],dimR[2])
                s2 = max( dimR[1],dimR[3])
                c0=0;c1=0;c2=0;c3=0;c4=0;c5=0
                count = numpy.zeros( 6, Internal.E_NpyInt)
                lmin  = numpy.zeros( 6, Internal.E_NpyInt)

                for l in range( numpy.size(ptListD)):
                   #i,j,k receveur
                    k =  ptListD[l]//(sh[0]*sh[1])
                    j = (ptListD[l] -k*sh[0]*sh[1])//sh[0]
                    i =  ptListD[l] -k*sh[0]*sh[1] -j*sh[0]
                    if dimPb==2:
                        if i==1 and j > 1 and j < dimR[2]-2:  #flux en imin
                            idir=0
                            if j < win[2,idir]: win[2,idir]=j; lmin[idir]= ptList[l]
                            if j > win[3,idir]: win[3,idir]=j
                            count[idir]+=1
                        elif i==dimR[1]-2 and j > 1 and j < dimR[2]-2: #flux en imax
                            idir=1
                            if j < win[2,idir]: win[2,idir]=j; lmin[idir]= ptList[l]
                            if j > win[3,idir]: win[3,idir]=j
                            count[idir]+=1
                        elif j==1 and i > 1 and i < dimR[1]-2: #flux en jmin
                            idir=2
                            if i < win[0,idir]: win[0,idir]=i; lmin[idir]= ptList[l]
                            if i > win[1,idir]: win[1,idir]=i
                            count[idir]+=1
                        elif j== dimR[2]-2 and i > 1 and i < dimR[1]-2: #flux en jmax
                            idir=3
                            if i < win[0,idir]: win[0,idir]=i; lmin[idir]= ptList[l]
                            if i > win[1,idir]: win[1,idir]=i
                            count[idir]+=1
                    else:  #Pb 3D
                       #print('k,j,i:', k,j,i)
                        if i==1 and j > 1 and j < dimR[2]-2 and k > 1 and k < dimR[3]-2:  #flux en imin

                            idir=0
                            if j < win[2,idir]: win[2,idir]=j
                            if j > win[3,idir]: win[3,idir]=j
                            if k < win[4,idir]: win[4,idir]=k
                            if k > win[5,idir]: win[5,idir]=k
                            if win[4,idir]==k and win[2,idir]==j: lmin[idir]= ptList[l]
                            count[idir]+=1

                        elif i==dimR[1]-2 and j > 1 and j < dimR[2]-2 and k > 1 and k < dimR[3]-2:  #flux en imin

                            idir=1
                            if j < win[2,idir]: win[2,idir]=j
                            if j > win[3,idir]: win[3,idir]=j
                            if k < win[4,idir]: win[4,idir]=k
                            if k > win[5,idir]: win[5,idir]=k
                            if win[4,idir]==k and win[2,idir]==j: lmin[idir]= ptList[l]
                            count[idir]+=1

                        elif j ==1 and i > 1 and i < dimR[1]-2 and k > 1 and k < dimR[3]-2    :
                            idir=2
                            if i < win[0,idir]: win[0,idir]=i
                            if i > win[1,idir]: win[1,idir]=i
                            if k < win[4,idir]: win[4,idir]=k
                            if k > win[5,idir]: win[5,idir]=k
                            if win[4,idir]==k and win[0,idir]==i: lmin[idir]= ptList[l]
                            count[idir]+=1

                        elif j == dimR[2]-2 and i > 1 and i < dimR[1]-2 and k > 1 and k < dimR[3]-2    :
                            idir=3
                            if i < win[0,idir]: win[0,idir]=i
                            if i > win[1,idir]: win[1,idir]=i
                            if k < win[4,idir]: win[4,idir]=k
                            if k > win[5,idir]: win[5,idir]=k
                            if win[4,idir]==k and win[0,idir]==i: lmin[idir]= ptList[l]
                            count[idir]+=1

                        elif k==1 and j > 1 and j < dimR[2]-2 and i > 1 and i < dimR[1]-2: #flux en kmin
                            idir=4
                            if i < win[0,idir]: win[0,idir]=i
                            if i > win[1,idir]: win[1,idir]=i
                            if j < win[2,idir]: win[2,idir]=j
                            if j > win[3,idir]: win[3,idir]=j
                            if win[2,idir]==j and win[0,idir]==i: lmin[idir]= ptList[l]
                            count[idir]+=1

                        elif k==dimR[3]-2 and j > 1 and j < dimR[2]-2 and i > 1 and i < dimR[1]-2: #flux en kmin
                            idir=5
                            if i < win[0,idir]: win[0,idir]=i
                            if i > win[1,idir]: win[1,idir]=i
                            if j < win[2,idir]: win[2,idir]=j
                            if j > win[3,idir]: win[3,idir]=j
                            if win[2,idir]==j and win[0,idir]==i: lmin[idir]= ptList[l]
                            count[idir]+=1

                #Adressage Fast
                #receveur
                win[0:6,:]-=1

                if dimPb==2:
                    win[4,:]=1
                    win[5,:]=1

                idirs=[]
                for i in range(6):
                    if count[i] !=0: idirs.append(i)

                #print('idirs', idirs)

                for idir in idirs:

                    #ijkminD,...: adresse C
                    kminD =  lmin[idir]//(shD[0]*shD[1])
                    jminD =  (lmin[idir] -kminD*shD[0]*shD[1])//shD[0]
                    iminD =  lmin[idir] -kminD*shD[0]*shD[1] -jminD*shD[0]
                    k1D = kminD-1
                    #k2D = kminD-2+(win[5,idir]-win[4,idir]+1)*2
                    k2D = kminD-2+(win[5,idir]-win[4,idir]+1)*ratio_k
                    if dimPb==2: k1D=1; k2D=1

                    if idir < 2:
                        sz=(win[3,idir]-win[2,idir]+1)*(win[5,idir]-win[4,idir]+1)
                        i1 =dimR[1]-3
                        i1D=iminD-1
                        name='imax'
                        if idir==0:
                            i1 =1
                            i1D= iminD+1
                            name='imin'
                        win[0:2 ,idir]=i1
                        winD[0:2,idir]=i1D
                        winD[2,idir]  =jminD-1
                        winD[3,idir]  =jminD-2+(win[3,idir]-win[2,idir]+1)*2
                        winD[4,idir]  =k1D
                        winD[5,idir]  =k2D
                    elif idir < 4:
                        sz=(win[1,idir]-win[0,idir]+1)*(win[5,idir]-win[4,idir]+1)
                        j1 =dimR[2]-3
                        j1D=jminD-1
                        name='jmax'
                        if idir==2:
                            j1 =1
                            j1D= jminD+1
                            name='jmin'
                        win[2:4 ,idir]=j1
                        winD[2:4,idir]=j1D
                        winD[0,idir]  =iminD-1
                        winD[1,idir]  =iminD-2+(win[1,idir]-win[0,idir]+1)*2
                        winD[4,idir]  =k1D
                        winD[5,idir]  =k2D
                    else:
                        sz=(win[1,idir]-win[0,idir]+1)*(win[3,idir]-win[2,idir]+1)
                        k1 =dimR[3]-3
                        k1D=kminD-1
                        name='kmax'
                        if idir==4:
                            k1 =1
                            k1D= kminD+1
                            name='kmin'
                        win[4:6 ,idir]=k1
                        winD[4:6,idir]=k1D
                        winD[0,idir]  =iminD-1
                        winD[1,idir]  =iminD-2+(win[1,idir]-win[0,idir]+1)*2
                        winD[2,idir]  =jminD-1
                        winD[3,idir]  =jminD-2+(win[3,idir]-win[2,idir]+1)*2

                    #print("verif count",sz , count[idir], ratio_k )
                    if sz== count[idir]:
                        if verbose==0:
                            #print('min ', imin,jmin,kmin, 'max ', imax, jmax,kmax)
                            #print('minD', iminD,jminD,kminD, k1D, k2D)
                            name1="#Flux_"+zRname+'_'+name
                            print("raccord conservatif: zD=", z[0], name1, 'win:', win[:,idir], 'winD:', winD[:,idir], 'taille win:', sz//2, 'min', imin,jmin,kmin)

                        name4 = 'Flux_'+zd_t[0]+'_'+name

                        #creation BC sur grille grossiere
                        if proc == rank:
                            zr_t= Internal.getNodeFromName(t, zRname)
                            #on nome le type BC avec name4 pour la retrouver apres et eviter le renommage des node par cassiopee
                            C._addBC2Zone(zr_t, name4, name4, wrange=win[:,idir])
                            bcs   = Internal.getNodesFromType2(zr_t, 'BC_t')
                            for bc in bcs:
                                btype = Internal.getValue(bc)
                                if btype == name4:
                                    Internal.setValue(bc,'BCFluxOctreeC')
                                    Prop = Internal.getNodeFromName(bc,'.Solver#Property')
                                    if Prop is None:
                                        Internal.createUniqueChild(bc,'.Solver#Property','UserDefinedData_t')
                                    Prop = Internal.getNodeFromName(bc,'.Solver#Property')

                                    ratioNM =  numpy.ones(3, numpy.float64)
                                    ratioNM[0]= ratio[zd_t[0]][0]
                                    ratioNM[1]= ratio[zd_t[0]][1]
                                    ratioNM[2]= ratio[zd_t[0]][2]
                                    Internal.createUniqueChild(Prop, 'ratioNM', 'DataArray_t', value=ratioNM)

                                    ptrange = Internal.getNodesFromType1(bc, 'IndexRange_t')
                                    rg  = ptrange[0][1]
                                    sz  = max(1, rg[0,1]-rg[0,0]+1) * max(1, rg[1,1]-rg[1,0]+1) * max(1, rg[2,1]-rg[2,0]+1)
                                    tab =  numpy.ones(sz*neq, numpy.float64)
                                    Internal.createUniqueChild(Prop, 'FluxFaces', 'DataArray_t', value=tab)

                        else:
                            if proc not in datas: datas[proc] = [ [ zRname, name4, win[:,idir], ratioNM ] ]
                            else: datas[proc] += [ [ zRname, name4, win[:,idir], ratioNM ] ]


                        #creation BC sur grille fine
                        if name[2]=='i':
                            name2= name[0:2]+'ax'
                        else:
                            name2= name[0:2]+'in'
                        name4 = 'Flux_'+ zRname +'_'+name2
                        C._addBC2Zone(zd_t, name4, name4, wrange=winD[:,idir])
                        bcs   = Internal.getNodesFromType2(zd_t, 'BC_t')
                        for bc in bcs:
                            btype = Internal.getValue(bc)
                            if btype == name4:
                                Internal.setValue(bc,'BCFluxOctreeF')
                                Prop = Internal.getNodeFromName(bc,'.Solver#Property')
                                if Prop is None:
                                    Internal.createUniqueChild(bc,'.Solver#Property','UserDefinedData_t')
                                Prop = Internal.getNodeFromName(bc,'.Solver#Property')

                                ratioNM =  numpy.ones(3, numpy.float64)
                                ratioNM[0]= ratio[zd_t[0]][0]
                                ratioNM[1]= ratio[zd_t[0]][1]
                                ratioNM[2]= ratio[zd_t[0]][2]
                                Internal.createUniqueChild(Prop, 'ratioNM', 'DataArray_t', value=ratioNM)

                                ptrange = Internal.getNodesFromType1(bc, 'IndexRange_t')
                                rg  = ptrange[0][1]
                                sz  = max(1, rg[0,1]-rg[0,0]+1) * max(1, rg[1,1]-rg[1,0]+1) * max(1, rg[2,1]-rg[2,0]+1)
                                tab =  numpy.ones(sz*neq, numpy.float64)
                                Internal.createUniqueChild(Prop, 'FluxFaces', 'DataArray_t', value=tab)

                    else:  print("Error: build flux conservative octree"+name, sz, count[idir],  win[:,idir])

    if Cmpi.size > 1:
        # Envoie des BC suivant le graph
        rcvDatas = Cmpi.sendRecv(datas, graph)

        # Remise des champs interpoles dans l'arbre receveur
        for i in rcvDatas:
            #print(rank, 'recoit de',i, '->', len(rcvDatas[i]), flush=True)
            for n in rcvDatas[i]:
                rcvName = n[0]
                #print('reception', Cmpi.rank, rcvName, flush=True)
                ptlistD = n[1]
                zr_t  = Internal.getNodeFromName(t,rcvName)
                C._addBC2Zone(zr_t, n[1], n[1], wrange=n[2])
                bcs   = Internal.getNodesFromType2(zr_t, 'BC_t')
                for bc in bcs:
                    btype = Internal.getValue(bc)
                    if btype == n[1]:
                        Internal.setValue(bc,'BCFluxOctreeC')
                        Prop = Internal.getNodeFromName(bc,'.Solver#Property')
                        if Prop is None:
                            Internal.createUniqueChild(bc,'.Solver#Property','UserDefinedData_t')
                        Prop = Internal.getNodeFromName(bc,'.Solver#Property')

                        ratioNM =  numpy.ones(3, numpy.float64)
                        ratioNM[0]= n[3][0]
                        ratioNM[1]= n[3][1]
                        ratioNM[2]= n[3][2]
                        Internal.createUniqueChild(Prop, 'ratioNM', 'DataArray_t', value=ratioNM)

                        ptrange = Internal.getNodesFromType1(bc, 'IndexRange_t')
                        rg  = ptrange[0][1]
                        sz  = max(1, rg[0,1]-rg[0,0]+1) * max(1, rg[1,1]-rg[1,0]+1) * max(1, rg[2,1]-rg[2,0]+1)
                        tab =  numpy.ones(sz*neq, numpy.float64)
                        Internal.createUniqueChild(Prop, 'FluxFaces', 'DataArray_t', value=tab)
