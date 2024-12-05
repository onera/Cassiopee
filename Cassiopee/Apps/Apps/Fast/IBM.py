# Class for FastS "IBM" prepare and compute
import FastS.ToolboxChimera as TBX
import Converter.PyTree as C
import Geom.PyTree as D
import Generator.PyTree as G
import Transform.PyTree as T
import Post.PyTree as P
import Converter.Internal as Internal
import Converter.GhostCells as Ghost
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

from Geom.IBM import setSnear, _setSnear, setDfar, _setDfar, snearFactor, _snearFactor, \
setFluidInside, _setFluidInside, setIBCType, _setIBCType, changeIBCType, \
initOutflow, initInj, _initOutflow, _initInj, transformTc2, _addOneOverLocally

import Post.IBM as P_IBM
import Connector.IBM as X_IBM
import Generator.IBM as G_IBM
import Generator.IBMmodelHeight as G_IBM_Height

KCOMM = Cmpi.KCOMM

RENAMEIBCNODES=False # renommage des ibcd*

__IBCNameServer__={}


## ============================= June 2024 =================================
## No further developments in this file. All IBM developments need to be
## performed in the IBM.py files in each Cassiopee module (e.g.$CASSIOPEE/Cassiopee/Connector/Connector/IBM.py)
## Any developments in this file will be REJECTED to be merged (pull request) in the official/upstream ONERA Cassiopee.
## Please address any questions to: Christophe Benoit, Benjamin Constant, Antoine Jost, and Stephanie Peron.

def getIBCDName(proposedName):
    global __IBCNameServer__
    (ibcname,__IBCNameServer__) = C.getUniqueName(proposedName, __IBCNameServer__)
    return ibcname

def _changeNameIBCD__(tc2, NewIBCD):
    ZsubR=Internal.getNodesByType(tc2, 'ZoneSubRegion_t')
    for z in ZsubR:
        zsplit = z[0].split('_')
        if zsplit[0] == 'IBCD' and zsplit[1]=='140':
            zsplit[1] = str(NewIBCD)
            znew = '_'.join(zsplit)
            Internal.setName(z, znew)
    return None   

def dist2wallNearBody__(t, tb, type='ortho', signed=0, dim=3, loc='centers', model='NSLaminar'):
    if model == 'NSLaminar':
        list_final_zones=[]
        for z in Internal.getZones(t):
            list_final_zones.append(z[0])

        tBB =G.BB(t)
        tbBB=G.BB(tb)

        ##PRT1 - Zones flagged by the intersection of the bounding boxes of t and tb
        interDict = X.getIntersectingDomains(tBB, tbBB)    
        zt       = []
        zt_names = []
        for i in interDict:
            if interDict[i]:
                zt.append(Internal.getNodeByName(t,i))
                zt_names.append(i)

        if zt_names:
            DTW._distance2Walls(zt, tb, type=type, signed=signed, dim=dim, loc=loc)
        del zt

        ###PRT2 - Zones flagged by the intersection of the bounding boxes of t and scaled version of tb
        list_additional_zones = getZonesScaleUpDown__(tbBB,tBB,zt_names,dim=dim)

        ##PRT2
        if list_additional_zones:        
            zt=[]
            for i in list_additional_zones:
                zt.append(Internal.getNodeByName(t,i))
            DTW._distance2Walls(zt, tb, type=type, signed=signed, dim=dim, loc=loc)

    else:
        DTW._distance2Walls(t, tb, type=type, signed=signed, dim=dim, loc=loc)    

    return None


def getZonesScaleUpDown__(tbBB, tBB, zt_names, diff_percent=0.15, sweep_num=4, scaleDirection=0, dim=2):
    mean_tb=[]
    for bb in Internal.getZones(tbBB):
        minval_tb = C.getMinValue(bb, ['CoordinateX', 'CoordinateY','CoordinateZ']); 
        maxval_tb = C.getMaxValue(bb, ['CoordinateX', 'CoordinateY','CoordinateZ']); 
        mean_tb.append(getMean__(maxval_tb,minval_tb))

    diff_percentz = diff_percent
    if dim == 2: diff_percentz=0

    list_additional_zones=[]
    list_additional_zonesCGNSFile=[]
    for i in range(1, sweep_num+1):
        if scaleDirection >= 0:            
            tbBB_scale = T.scale(tbBB, factor=(1.0+i*diff_percent,1.0+i*diff_percent,1.0+i*diff_percentz))
            add2listAdditionalZones__(list_additional_zones,tbBB_scale,tBB,mean_tb,zt_names)

        if scaleDirection<=0:
            tbBB_scale = T.scale(tbBB, factor=(1.0-i*diff_percent,1.0-i*diff_percent,1.0-i*diff_percentz))
            add2listAdditionalZones__(list_additional_zones,tbBB_scale,tBB,mean_tb,zt_names)

    return list_additional_zones


def getMean__(max_local, min_local):
    mean_local=[]
    for i in range(len(max_local)):
        mean_local.append((max_local[i]+min_local[i])/2)
    return mean_local


def add2listAdditionalZones__(list_additional_zones, tbBB_scale, tBB,mean_tb, zt_names):
    count=0
    for bb in Internal.getZones(tbBB_scale):
        minval_tbscale = C.getMinValue(bb, ['CoordinateX', 'CoordinateY','CoordinateZ']); 
        maxval_tbscale = C.getMaxValue(bb, ['CoordinateX', 'CoordinateY','CoordinateZ']);
        mean_tbscale   = getMean__(maxval_tbscale,minval_tbscale)
        T._translate(bb, (mean_tb[count][0]-mean_tbscale[0],mean_tb[count][1]-mean_tbscale[1],mean_tb[count][2]-mean_tbscale[2]))

        interDict_scale = X.getIntersectingDomains(tBB, bb)
        for i in interDict_scale:
            if interDict_scale[i] and i not in list_additional_zones and i not in zt_names:
                list_additional_zones.append(i)
        count += 1        
    return None


#=============================================================================
# Post - General
# IN: t_case: geometry file name or tree
# IN: t_in: result file name or tree
# IN: tc_in: connectivity file name or tree
# OUT: t_out ou None: output file name - values at nodes
# OUT: wall_out ou None: output file name - wall values
#==============================================================================
def post(t_case, t_in, tc_in, t_out, wall_out):
    """Extracts the target point values from tc to the geometry surface mesh tb.  Move the solution from the cell centers to the nodes in t..
        Usage: post(t_case, t_in, tc_in, t_out, wall_out)"""
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
# IN: famZones: list of family names of surface zones on which the solution is projected
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
    if famZones: ts[2][1][2] = zw
    else: ts[2][1][2] = zw[2][1][2]
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

def _projectMeshSize(t, NPas=10, span=1, dictNz=None, isCartesianExtrude=False):
    """Predicts the final size of the mesh when extruding 2D to 3D in the z-direction.
    Usage: loads(t, NPas, span, dictNz, isCartesianExtrude)"""
    NP             = Cmpi.size
    rank           = Cmpi.rank
    NPTS           = numpy.zeros(NP, dtype=Internal.E_NpyInt)
    NCELLS         = numpy.zeros(NP, dtype=Internal.E_NpyInt)
    NPTS_noghost   = numpy.zeros(NP, dtype=Internal.E_NpyInt)
    NCELLS_noghost = numpy.zeros(NP, dtype=Internal.E_NpyInt)    
    if isinstance(t, str):
        h = Filter.Handle(t)
        t = h.loadFromProc(loadVariables=False)
        h._loadVariables(t, var=['CoordinateX'])


    for z in Internal.getZones(t):
        name_zone = z[0]
        if not isCartesianExtrude:
            if dictNz is not None: Nk = int(dictNz[name_zone])
            else: Nk = NPas-1
        else:
            h = abs(C.getValue(z,'CoordinateX',0)-C.getValue(z,'CoordinateX',1))
            NPas_local = int(round(span/h)/4)
            if NPas_local < 2:
                print("WARNING:: MPI rank %d || Zone %s has Nz=%d and is being clipped to Nz=2"%(rank,z[0],NPas_local))
                NPas_local = 2
            Nk = NPas_local
        Nk += 1
        NPTS[rank]           += C.getNPts(z)/2*(Nk+4)
        NCELLS[rank]         += C.getNCells(z)*(Nk+3)
        NPTS_noghost[rank]   += C.getNPts(C.rmGhostCells(z, z, 2))*Nk
        NCELLS_noghost[rank] += C.getNCells(C.rmGhostCells(z, z, 2))*(Nk-1)
    NPTS             = Cmpi.allreduce(NPTS  ,op=Cmpi.SUM)
    NCELLS           = Cmpi.allreduce(NCELLS,op=Cmpi.SUM)
    NPTS_noghost     = Cmpi.allreduce(NPTS_noghost  ,op=Cmpi.SUM)
    NCELLS_noghost   = Cmpi.allreduce(NCELLS_noghost,op=Cmpi.SUM)   
    if rank ==0:
        print('Projected mesh size with ghost: {} million points & {} million cells'.format(numpy.sum(NPTS)/1.e6,numpy.sum(NCELLS)/1.e6))
        print('Projected mesh size without ghost: {} million points & {} million cells'.format(numpy.sum(NPTS_noghost)/1.e6,numpy.sum(NCELLS_noghost)/1.e6))
    return None


#================================================================================
# extrude mesh and case for z or theta periodic cases
#================================================================================
def extrudeCartesian(t,tb, check=False, extrusion="cart", dz=0.01, NPas=10, span=1 , Ntranche=1,
                     dictNz=None, ific=2,isCartesianExtrude=False,isAutoPeriodic=False):
    """Extrudes a 2D IBM grid and geoemtry. The extraction is done in the z-direction..
        Usage: extrudeCartesian(t, tb, check, extrusion, dz, NPas, span, Ntranche, dictNz, ific, isCartesianExtrude)"""

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
            Nk_min = 1000000000
            for z in Internal.getZones(tree):
                name_zone = z[0]
                if not isCartesianExtrude:
                    if dictNz is not None:
                        Nk[z[0]]     = int(dictNz[name_zone])
                        dz_loc[z[0]] = span/float(Ntranche*Nk[z[0]])
                    else:
                        Nk[z[0]]     = NPas-1
                        dz_loc[z[0]] = dz
                else:
                    if tree == t:
                        h = abs(C.getValue(z,'CoordinateX',0)-C.getValue(z,'CoordinateX',1))
                        NPas_local = int(round(span/h))
                        if NPas_local<4:
                            print("WARNING:: Zone %s has Nz=%d and is being clipped to Nz=4"%(z[0],NPas_local))
                            NPas_local=4

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

        ##Adding ghost cells in K direction
        for z in Internal.getZones(tree):
            Nk[z[0]] += 2*ific-1 #+c  # -1, car un plan de maillage deja ajoute dans le cas 2D

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

            if  len(sh_2d) == 1: #NGON
                if  extrusion == "cart":
                    for k in range( sh[1]):
                        zz[:,k] = (k-ific)* dz_loc[z[0]]
                else:
                    theta0 = numpy.arctan(zz_2d/yy_2d)
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
        c += 1
        ng = 0
        #Modif Ivan 11/10/24
        #ng = 2
        for z in Internal.getZones(tree):    
            zdim = Internal.getValue(z)

            # Modif rind cellule GH en kmin et kmax
            rind = Internal.getNodeFromName1(z, "Rind")[1]
            rind[4] = ific; rind[5] = ific
            # Modif range BC
            BCs = Internal.getNodesFromType(z, "BC_t")
            for bc in BCs:
                ptrg = Internal.getNodeFromName(bc, "PointRange")[1]
                ptrg[2,0] = 3 
                ptrg[2,1] = Nk[z[0]]
            if not isAutoPeriodic:
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
                    datap = numpy.empty( (3,2) , dtype=Internal.E_NpyInt)
                    datap[0,0]=1+ng;datap[1,0]=1+ng;datap[2,0]=ktg
                    datap[0,1]=zdim[0,0]-ng;datap[1,1]=zdim[1,0]-ng;datap[2,1]= ktg
                    Internal.createUniqueChild(tmp1 ,'PointRange', 'IndexRange_t',datap)
                    datap = numpy.empty( (3,2) , dtype=Internal.E_NpyInt)
                    datap[0,0]=1+ng;datap[1,0]=1+ng;datap[2,0]=ktgD
                    datap[0,1]=zdim[0,0]-ng;datap[1,1]=zdim[1,0]-ng;datap[2,1]= ktgD
                    Internal.createUniqueChild(tmp1 ,'PointRangeDonor', 'IndexRange_t',datap)
                    datap = numpy.empty( 3 , dtype=Internal.E_NpyInt)
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
        if isAutoPeriodic:    
            for node in Internal.getNodesFromName(t,'EquationDimension'): Internal.setValue(node,3)
            for z in Internal.getZones(t):
                C._addBC2Zone(z, z[0]+'periodKmin', 'BCautoperiod', 'kmin')
                C._addBC2Zone(z, z[0]+'periodKmax', 'BCautoperiod', 'kmax')
            BCs = Internal.getNodesFromType(t, "BC_t")
            for bc in BCs:
                if Internal.getValue(bc)=='BCautoperiod':
                    ptrg = Internal.getNodeFromName(bc, "PointRange")[1]
                    ptrg[0,0] = 3 
                    ptrg[0,1] = ptrg[0,1]-2

                    ptrg[1,0] = 3 
                    ptrg[1,1] = ptrg[1,1]-2


    return t, tb

#================================================================================
# Chimere: modif ghost planSym pour interp chimere 
# IN: t
# IN: depth= nombre de couche
# IN: plan= plan cible de la symmetrie ('X': x=cte= plan sym, ...)
#================================================================================
def _modifGhostSymmetryPlan(t, depth=2,dim=3,Plan='Y'):

    for z in Internal.getZones(t):
        bc = Internal.getNodeFromName1(z,'ZoneBC')
        bcs=[]
        if bc is not None: bcs = Internal.getNodesFromType(bc,'BC_t')
        for bc in bcs:
            bctyp = Internal.getValue(bc)
            if bctyp=='BCSymmetryPlane': 
                ptrg = Internal.getNodeFromName1(bc,'PointRange')
                dimZ = Internal.getZoneDim(z)
                #print(z[0], dimZ,'range', ptrg[1])

                idir = Ghost.getDirection__(dim, [ptrg])

                #print('idir=',idir)
                var = 'Coordinate'+Plan
                coord= Internal.getNodeFromName2(z,var)[1]

                if idir==0: 
                    jmin =  ptrg[1][1,0]-3; jmax =  ptrg[1][1,1]+2
                    kmin =  ptrg[1][2,0]-3; kmax =  ptrg[1][2,1]+2
                    for k in range(kmin,kmax):
                        for j in range(jmin,jmax):
                            for d in range(depth):
                                coord[depth-1-d,j,k]=  2*coord[depth-d,j,k]-coord[depth+1-d,j,k]
                if idir==1: 
                    i1 = ptrg[1][0,1]-depth-1
                    jmin =  ptrg[1][1,0]-3; jmax =  ptrg[1][1,1]+2
                    kmin =  ptrg[1][2,0]-3; kmax =  ptrg[1][2,1]+2
                    for k in range(kmin,kmax):
                        for j in range(jmin,jmax):
                            for d in range(depth):
                                coord[i1+1+d,j,k]=   2.*coord[i1+d  ,j,k]-coord[i1-1+d,j,k]
                if idir==2: 
                    imin =  ptrg[1][0,0]-3; imax =  ptrg[1][0,1]+2
                    kmin =  ptrg[1][2,0]-3; kmax =  ptrg[1][2,1]+2
                    for k in range(kmin,kmax):
                        for i in range(imin,imax):
                            for d in range(depth):
                                coord[i,depth-1-d,k]=  2*coord[i,depth-d,k]-coord[i,depth+1-d,k]
                if idir==3: 
                    js = ptrg[1][1,1]-depth-1
                    imin =  ptrg[1][0,0]-3; imax =  ptrg[1][0,1]+2
                    kmin =  ptrg[1][2,0]-3; kmax =  ptrg[1][2,1]+2
                    for k in range(kmin,kmax):
                        for i in range(imin,imax):
                            for d in range(depth):
                                coord[i,js+1+d,k]=   2.*coord[i,js+d,k]-coord[i,js-1+d,k]
                if idir==4: 
                    imin =  ptrg[1][0,0]-3; imax =  ptrg[1][0,1]+2
                    jmin =  ptrg[1][1,0]-3; jmax =  ptrg[1][1,1]+2
                    for j in range(jmin,jmax):
                        for i in range(imin,imax):
                            for d in range(depth):
                                coord[i,j,depth-1-d]=  2*coord[i,j,depth-d]-coord[i,j,depth+1-d]

                if idir==5: 
                    ks = ptrg[1][2,1]-depth-1
                    imin =  ptrg[1][0,0]-3; imax =  ptrg[1][0,1]+2
                    jmin =  ptrg[1][1,0]-3; jmax =  ptrg[1][1,1]+2
                    for j in range(jmin,jmax):
                        for i in range(imin,imax):
                            for d in range(depth):
                                coord[i,j,ks+1+d]=   2.*coord[i,j,ks+d]-coord[i,j,ks-1+d]
    return None

#================================================================================
# Chimere: modif cellN planSym pour interp chimere 
# IN: t
# IN: depth= nombre de couche
# IN: plan= plan cible de la symmetrie ('X': x=cte= plan sym, ...)
#================================================================================
def _modifcellNSymmetryPlan(t, depth=2,dim=3):

    for z in Internal.getZones(t):
        bc = Internal.getNodeFromName1(z,'ZoneBC')
        bcs=[]
        if bc is not None: bcs = Internal.getNodesFromType(bc,'BC_t')
        for bc in bcs:
            bctyp = Internal.getValue(bc)
            if bctyp=='BCSymmetryPlane': 
                ptrg = Internal.getNodeFromName1(bc,'PointRange')
                dimZ = Internal.getZoneDim(z)
                #print(z[0], dimZ,'range', ptrg[1])

                idir = Ghost.getDirection__(dim, [ptrg])

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

#================================================================================
# IBM prepare hybride cart + curvi
# IN: t_octree arbre cartesien   (from prepare1)
# IN: tc_octree arbre cartesien connectivite   (from prepare1)
# IN: t_curvi arbre curviligne   (with bc and connectivity but without ghost)
#================================================================================
def setInterpData_Hybride(t_octree, tc_octree, t_curvi, blankingMatrix=None, blankingBodies=None, extrusion=None, interpDataType=1, order=2, verbose=2, filter_zone=None):
    overlap = '14'
    InterpOrder = order
    Nproc =Cmpi.size
    rank =Cmpi.rank

    dim = 3
    node = Internal.getNodeFromName(t_curvi, 'EquationDimension')
    if node is not None: dim = Internal.getValue(node)

    Model='NSTurbulent'
    node= Internal.getNodeFromName(t_curvi, 'GoverningEquation')
    if node is not None: Model = Internal.getValue(node)

    Cmpi.barrier()
    #ajout 2 ghost maillag curvi
    if Nproc !=1:
        C._deleteGridConnectivity__(t_curvi, type='BCMatch')
        bases_1=[];bases_2=[]; bases_3=[];bases_4=[]
        for b in Internal.getBases(t_curvi):
            #if b[0]=='ROTOR' or b[0]=='STATOR': bases_1.append(b)
            if   b[0]=='Rx': bases_1.append(b)
            elif b[0]=='Fuselage' or b[0]=='Aile': bases_2.append(b)
            elif b[0]=='Stator' or 'STATOR' in b[0]: bases_3.append(b)
            elif b[0]=='Rotor'  or 'ROTOR'  in b[0]: bases_4.append(b)
            #else: bases_4.append(b)
            #if 'I1' in b[0] or 'I2' in b[0] or 'I3' in b[0] or 'I4' in b[0] or 'I5' in b[0]: bases_3.append(b)
            #if 'Proche' in b[0] or 'Far' in b[0]: bases_3.append(b)

        bout=[]
        c= 0
        #ATTENTION: ordre des bases doit respecter l'ordre employe pour creer la matrice BM, sinon blanking foireux
        #for b in [bases_3, bases_1, bases_2]:
        for b in [ bases_1, bases_2, bases_3, bases_4]:
            #print('traitement base:', b[0][0], flush=True)
            ##Cmpi._addGhostCells(tmp, tmp, 2, adaptBCs=1, fillCorner=0, modified=None)
            #Cmpi._addGhostCells(b, b, 2, adaptBCs=1, fillCorner=0, modified=None)
            ##Cmpi._addGhostCells(b, b, 2, adaptBCs=1, fillCorner=0, modified=[['centers:TurbulentDistance']])
            tmp = C.newPyTree(b)
            Cmpi._addBXZones(tmp, depth=3, NoVar=True)
            Cmpi.barrier()
            tmp = X.connectMatch(tmp, dim=dim)
            Internal._addGhostCells(tmp, tmp, 2, adaptBCs=1, fillCorner=0)
            Cmpi.barrier()
            Cmpi._rmBXZones(tmp)
            baseLoc = Internal.getBases(tmp)
            '''
          baseLoc = Internal.getBases(b)
          '''
            if len(baseLoc)!=0: bout= bout + baseLoc

        t_curvi = C.newPyTree(bout)

    else: Internal._addGhostCells(t_curvi, t_curvi, 2, adaptBCs=1, fillCorner=0)


    #print(" addGhost OK", rank, flush=True)
    Cmpi.barrier()

    #Ajout 1 plan k pour cas bidim
    if dim==2:
        t_curvi = T.addkplane(t_curvi)
        t_curvi = T.contract(t_curvi, (0,0,0), (1,0,0), (0,1,0), 0.01)
        t_curvi = T.makeDirect(t_curvi)

    #creation "vraie" cellule ghost pour BC plansymetry pour faciliter interp Chimere
    _modifGhostSymmetryPlan(t_curvi, depth=2, dim=dim)

    #calcul distance paroi si necessaire
    ret = C.isNamePresent(t_curvi, 'centers:TurbulentDistance')
    #print(" avt Dist: rank=", rank, 'ret=', ret, Model, flush=True)
    Cmpi.barrier()
    #for z in Internal.getZones(t_curvi):
    #   sol= Internal.getNodeFromName(z,'FlowSolution#Centers')
    #   dist = Internal.getNodeFromName(sol,'TurbulentDistance')

    #ret=1
    #if Model == 'NSTurbulent' and ret != 1: # Wall distance not present
    if Model == 'NSTurbulent':
        import Dist2Walls.PyTree as DTW
        walls = C.extractBCOfType(t_curvi, 'BCWall')
        wall_gather = Cmpi.allgather(walls)
        #print(" apr wall gather: rank=", rank, flush=True)
        if rank==0: C.convertPyTree2File(walls,'wall.cgns')
        Cmpi.barrier()
        walls = []
        for i in wall_gather: walls += i
        if ret != 1:
            #print(" Calcul Dist pourquoi: rank=", rank, 'ret=', ret, Model, flush=True)
            if walls != []: # Wall distance not present:
                DTW._distance2Walls(t_curvi, walls, loc='centers', type='ortho')
            else:
                C._initVars(t_curvi, 'centers:TurbulentDistance', 1000000.)

    #print('apr Dist', flush=True)
    Cmpi.barrier()

    #init cellN pour interp chimere sur grille curvi
    C._initVars(t_curvi,"centers:cellN",1.)
    t_curvi = TBX.modifyBCOverlapsForGhostMesh(t_curvi,2)

    for z in Internal.getZones(t_curvi):
        apply=False
        if  filter_zone is None: apply=True
        elif filter_zone in z[0]: apply=True
        if apply:
            C._fillEmptyBCWith(z,'inactive','BCExtrapolated',dim=dim)
            C._rmBCOfType(z,'BCOverlap')
            C._fillEmptyBCWith(z,'trou','BCOverlap',dim=dim)
            X._applyBCOverlaps(z,loc='centers',depth=4, val=2)
            X._applyBCOverlaps(z,loc='centers',depth=2, val=0)

    # blanking des zones curvi si necessaire 
    if blankingMatrix is not None:
        if dim==2: 
            print('revoir blanking en 2d')
            import sys; sys.exit()

        bases = Internal.getBases(t_curvi); nbases  = len(bases)

        for base in bases: print("verif base order", base[0])

        X._blankCells(t_final, blankingBodies, blankingMatrix, cellNName='cellN', XRaydim1=300, XRaydim2=300, dim=3, delta=0.)
        '''
        BM = numpy.zeros((nbases, 1),dtype=numpy.int32)
        c=0
        for body in blankingBodies:           
          for b in range(numpy.shape(blankingMatrix)[0]): BM[b,0]= blankingMatrix[b,c]
          #X._blankCellsTri(t_curvi, [body], blankingMatrix, cellNName='cellN', blankingType='center_in')
          print('bm', BM[:,0], c)
          X._blankCellsTri(t_curvi, [body], BM, cellNName='cellN', blankingType='cell_intersect')
          c+=1
        '''
        #print('apres blanking',flush=True)

        X._setHoleInterpolatedPoints(t_curvi, depth=2, loc='centers', addGC=False, cellNName='cellN', dir=2)

        #print('apres HoleInterpolated',flush=True)
        _modifcellNSymmetryPlan(t_curvi, depth=2, dim=dim)
        X._modifcellNBC(t_curvi, ['BCWall'], depth=2, dim=dim)

    #print('apr planSym', flush=True)
    Cmpi.barrier()

    #Passage en coordonne cylindrique si necessaire a l'interp
    if extrusion == 'cyl':
        T._cart2Cyl(t_octree , (0,0,0),(1,0,0))
        T._cart2Cyl(tc_octree, (0,0,0),(1,0,0))
        T._cart2Cyl(t_curvi   , (0,0,0),(1,0,0))


    #calcul raccord match sur zone curvi
    tc_curvi = C.node2Center(t_curvi)
    if Internal.getNodeFromType(t_curvi, "GridConnectivity1to1_t") is not None:
        if Nproc !=1:
            Xmpi._setInterpData(t_curvi, tc_curvi, nature=1, loc='centers', storage='inverse', sameName=1, dim=dim, itype='abutting')
        else:
            X._setInterpData(t_curvi, tc_curvi, nature=1, loc='centers', storage='inverse', sameName=1, dim=dim, itype='abutting')

    #print('apr abbuting', flush=True)
    Cmpi.barrier()
    #sauvegarde cellN
    C._initVars(t_octree,"{centers:cellN#Init}={centers:cellN}")
    C._initVars(t_curvi ,"{centers:cellN#Init}={centers:cellN}")

    #
    #tag les points interpoles cartesion de type wall et overlap a l'aide du pointListD du tc cartesian
    #
    for var in ['wall','racChimer']:
        for z in Internal.getZones(t_octree):
            sol            = Internal.getNodeFromName(z,'FlowSolution#Centers')
            cellN          = Internal.getNodeFromName(sol,'cellN')[1]
            sh             = numpy.shape(cellN)

            for k in range(sh[2]):
                for j in range(sh[1]):
                    for i in range(sh[0]):
                        if  cellN[i,j,k] != 0:  cellN[i,j,k] =1

        if var == 'wall': itype ="3"
        else: itype = overlap
        #creation nouveau tc pour garder que les raccords cibles (overlap, wall,...)
        tc_filtre = Internal.copyRef(tc_octree)
        for zc in Internal.getZones(tc_filtre):
            for zsr in Internal.getNodesFromType(zc, "ZoneSubRegion_t"):
                zsrname = Internal.getName(zsr)
                zsrname = zsrname.split('_')
                if zsrname[0]!='IBCD' or zsrname[1] != itype:
                    Internal._rmNodesFromName(zc,zsr[0])

        #construction arbre skeleton global (tcs) pour calcul graph          
        tcs_local = Cmpi.convert2SkeletonTree(tc_filtre)
        tcs       = Cmpi.allgatherTree(tcs_local)
        graph     = Cmpi.computeGraph(tcs, type='IBCD', reduction=True)
        procDict  = Cmpi.getProcDict(tcs)

        #construction list pour envoie pointListD  zone distante ou modif direct cellN si zone local
        datas = {}
        for zc in Internal.getZones(tc_filtre):
            for zsr in Internal.getNodesFromType(zc, "ZoneSubRegion_t"):
                zrname = Internal.getValue(zsr)
                ptlistD= Internal.getNodeFromName(zsr,'PointListDonor')[1]
                proc   = procDict[zrname]
                if proc == rank:
                    zloc   = Internal.getNodeFromName(t_octree,zrname)
                    sol    = Internal.getNodeFromName(zloc,'FlowSolution#Centers')
                    #print("shape", zrname, zloc[0], numpy.shape(ptlistD),  numpy.size(ptlistD), ptlistD[0], flush=True)
                    cellN  = Internal.getNodeFromName(sol,'cellN')[1]
                    sh     = numpy.shape(cellN)
                    ni= sh[0]; ninj = sh[0]*sh[1]
                    for l0 in range(numpy.size(ptlistD)):
                        l = ptlistD[ l0]
                        k = l//ninj
                        j = (l-k*ninj)//ni
                        i =  l-k*ninj - j*ni
                        cellN[ i,j,k ]=2
                else:
                        #print("shapeD", numpy.shape(ptlistD),  numpy.size(ptlistD),flush=True)
                    if proc not in datas: datas[proc] = [[zrname, ptlistD]]
                    else: datas[proc] += [[zrname, ptlistD]]
                    #print ("dataS", proc,  datas, flush=True)
        # Envoie des numpys suivant le graph
        rcvDatas = Cmpi.sendRecv(datas, graph)

        # Remise des champs interpoles dans l'arbre receveur
        for i in rcvDatas:
            #print(rank, 'recoit de',i, '->', len(rcvDatas[i]), flush=True)
            for n in rcvDatas[i]:
                rcvName = n[0]
                #print('reception', Cmpi.rank, rcvName, flush=True)
                ptlistD = n[1]
                zloc    = Internal.getNodeFromName(t_octree,rcvName)
                #if zloc is not None: print('la zone', zloc[0],'est OK, size PtlistD=', numpy.size(ptlistD) , flush=True)
                #else: print('aille', flush=True)

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

        if var == 'wall': C._initVars(t_octree,"{centers:cellN#wall}={centers:cellN}")
        else:             C._initVars(t_octree,"{centers:cellN#racChim}={centers:cellN}")
    #
    # fin Flag cellN=2 pour les points interpoles cartesien associes a la condition overlap
    #

    #on purge le tc cartesien des raccord overlap IBCD. SEra pris en compte par raccord ID plus loin
    Internal._rmNodesFromName(tc_octree,'IBCD_'+overlap+'_*')
    #Recopie cellN dans tc  octree en vue interp chimere
    for z in Internal.getZones(t_octree):
        zc = Internal.getNodeFromName(tc_octree, z[0])
        C._cpVars( z, 'centers:cellN', zc, 'cellN')

    #print('apr manip cellN', flush=True)
    Cmpi.barrier()

    # 
    # double wall operation
    # 
    # tag for double wall static
    for n in ['Rotor', 'Stator']:
        base = Internal.getNodeFromName1(t_curvi, n)
        for z in Internal.getNodesFromName(base, '*contact*'):
            bcs = Internal.getNodesFromType(z, 'BC_t')
            for bc in bcs:
                bcname=Internal.getValue(bc)
                if bcname != 'BCExtrapolated': C._tagWithFamily(bc, 'WALL_{}'.format(n))

    X._doubleWall(t_curvi,tc_curvi,'WALL_Rotor','WALL_ROTORBLADE',ghostCells=True,check=False)
    #Internal._rmNodesByName(t_curvi,'*inactive*')
    Cmpi.barrier()
    #if rank == 0: print("DB : {} -> {}".format('ROTORBLADES', 'ROTOR'), flush=True)
    Cmpi.barrier()

    base1 = Internal.getNodeFromName1(t_curvi, 'Rotor')
    base1c = Internal.getNodeFromName1(tc_curvi, 'Rotor')
    for i in range(14):
        base2 = Internal.getNodeFromName1(t_curvi, 'ROTOR%d'%(i+1))
        base2c = Internal.getNodeFromName1(tc_curvi, 'ROTOR%d'%(i+1))
        tp = C.newPyTree([base1, base2])
        tpc = C.newPyTree([base1c, base2c])
        #if rank == 0: 
        #   C.convertPyTree2File(tp,'bug'+str(rank)+'.cgns')
        #   C.convertPyTree2File(tpc,'bugC'+str(rank)+'.cgns')
        Cmpi.barrier()
        X._doubleWall(tp,tpc,'WALL_ROTORBLADE','WALL_Rotor',ghostCells=True,check=False) #Souci 
        #Internal._rmNodesByName(tp,'*inactive*')
        Cmpi.barrier()
        #if rank == 0: print("DB : {} -> {}".format('ROTOR', 'ROTORBLADE%d'%(i+1)), flush=True)
        Cmpi.barrier()

    X._doubleWall(t_curvi,tc_curvi,'WALL_Stator','WALL_STATORBLADE',ghostCells=True,check=False)
    #Internal._rmNodesByName(t_curvi,'*inactive*')
    Cmpi.barrier()
    #if rank == 0: print("DB : {} -> {}".format('STATORBLADES', 'STATOR'), flush=True)
    Cmpi.barrier()

    base1 = Internal.getNodeFromName1(t_curvi, 'Stator')
    base1c = Internal.getNodeFromName1(tc_curvi, 'Stator')
    for i in range(12):
        base2 = Internal.getNodeFromName1(t_curvi, 'STATOR%d'%(i+1))
        base2c = Internal.getNodeFromName1(tc_curvi, 'STATOR%d'%(i+1))
        tp = C.newPyTree([base1, base2])
        tpc = C.newPyTree([base1c, base2c])
        X._doubleWall(tp,tpc,'WALL_STATORBLADE','WALL_Stator',ghostCells=True,check=False)
        #Internal._rmNodesByName(tp,'*inactive*')
        Cmpi.barrier()
        #if rank == 0: print("DB : {} -> {}".format('STATOR', 'STATORBLADE%d'%(i+1)), flush=True)
        Cmpi.barrier()

    # tag for double wall dynamic
    for n in ['Rotor', 'Stator']:
        base = Internal.getNodeFromName1(t_curvi, n)
        for z in Internal.getNodesFromName(base, '*contactR*' if n == 'Stator' else '*contactS*'):
            bcs = Internal.getNodesFromType(z, 'BC_t')
            for bc in bcs:
                C._tagWithFamily(bc, 'WALL_{}_DYNAMIC'.format(n))

    #C.convertPyTree2File(t_curvi, 'CHECK/t_test3_'+str(rank)+'.cgns')



    if Nproc==1:
        tBB2 = G.BB(t_curvi)
        tBB  = G.BB(t_octree)
        interDict = X.getIntersectingDomains(tBB, t2=tBB2, method='AABB', taabb=tBB, taabb2=tBB2)
        zonelist=[]
        for z in Internal.getZones(t_octree):
            if C.getMaxValue(z,'centers:cellN')==2:
                dnrZones = []    
                for zdname in interDict[z[0]]:
                    zd = Internal.getNodeFromName(tc_curvi,zdname)
                    dnrZones.append(zd)
                X._setInterpData(z,dnrZones, nature=0,penalty=1,loc='centers',storage='inverse',sameName=1,\
                                 interpDataType=interpDataType, itype='chimera', order=InterpOrder,verbose=verbose)
                z = X.getOversetInfo(z, dnrZones, loc='center',type='extrapolated')
                zonelist.append(z)

        interDict_curvi = X.getIntersectingDomains(tBB2, t2=tBB, method='AABB', taabb=tBB2, taabb2=tBB)
        print(" Interp curvi from Cartesian")
        for z in Internal.getZones(t_curvi):
            if C.getMaxValue(z,'centers:cellN')==2:
                dnrZones = []    
                for zdname in interDict_curvi[z[0]]:
                    zd = Internal.getNodeFromName(tc_octree,zdname)
                    dnrZones.append(zd)

                X._setInterpData(z,dnrZones,nature=0,penalty=1,loc='centers',storage='inverse',sameName=1,\
                                 interpDataType=interpDataType, itype='chimera', order=InterpOrder,verbose=verbose)

                # to check orphans
                #z=X.getOversetInfo(z, dnrZones, loc='centers',type='extrapolated')
                #z = X.getOversetInfo(z, dnrZones, loc='centers',type='extrapolated')
                #z = X.getOversetInfo(z, dnrZones, loc='centers',type='orphan')
                #zonelist.append(z)
                #C.convertPyTree2File(z,"rcv2.cgns")

        t  = C.mergeTrees(t_octree ,t_curvi )
        tc = C.mergeTrees(tc_octree,tc_curvi)

    else:
        t  = C.mergeTrees(t_octree ,t_curvi )
        tc = C.mergeTrees(tc_octree,tc_curvi)

        '''
       #interp chimere curvi classique
       C._initVars(t,"{centers:cellN}={centers:cellN#Init}")
       for base in Internal.getBases(t):
         if base[0]=='CARTESIAN':
           for z in Internal.getZones(base):
             sol            = Internal.getNodeFromName(z,'FlowSolution#Centers')
             cellN          = Internal.getNodeFromName(sol,'cellN')[1]
             sh             = numpy.shape(cellN)
             for k in range(sh[2]):
               for j in range(sh[1]):
                 for i in range(sh[0]):
                     if  cellN[i,j,k] == 2:  cellN[i,j,k] =0

       #Recopie cellN dans tc  octree en vue interp chimere
       for z in Internal.getZones(t):
          zc = Internal.getNodeFromName(tc, z[0])
          C._cpVars( z, 'centers:cellN', zc, 'cellN')
       Xmpi._setInterpData(t, tc, nature=1, loc='centers', storage='inverse', sameName=1, sameBase=0, dim=dim, itype='chimera', order=InterpOrder,verbose=verbose)
       
       for base in Internal.getBases(t):
         if base[0]=='CARTESIAN':
           C._initVars(base,"{centers:cellN}={centers:cellN#racChim}")
         else:
           for z in Internal.getZones(base):
             sol            = Internal.getNodeFromName(z,'FlowSolution#Centers')
             cellN          = Internal.getNodeFromName(sol,'cellN')[1]
             sh             = numpy.shape(cellN)
             for k in range(sh[2]):
               for j in range(sh[1]):
                 for i in range(sh[0]):
                     if  cellN[i,j,k] != 0:  cellN[i,j,k] =1
       #Recopie cellN dans tc  octree en vue interp chimere
       for z in Internal.getZones(t):
          zc = Internal.getNodeFromName(tc, z[0])
          C._cpVars( z, 'centers:cellN', zc, 'cellN')
       #interp chimere cartIBM/curvi: on impose nature 0 pour eviter extrap couche limite
       '''
        Xmpi._setInterpData(t, tc, nature=0, loc='centers', storage='inverse', sameName=1, sameBase=0, dim=dim, itype='chimera', order=InterpOrder,verbose=verbose)


    #reaffection du vrai cellN 
    C._initVars(t, '{centers:cellN}={centers:cellN#Init}')

    if extrusion == 'cyl':
        T._cyl2Cart(t,   (0,0,0),(1,0,0))
        T._cyl2Cart(tc,  (0,0,0),(1,0,0))

    for z in Internal.getZones(t):
        for bc in  Internal.getNodesFromType(z,'BC_t'):
            if 'inactive' in bc[0] and Internal.getValue(bc) == 'BCExtrapolated':
                Internal._rmNodesByName(z, bc[0])

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

    _checkNcellsNptsPerProc(ts,NP)
    return None


#====================================================================================
# Check number of points and cells per zone & in total
#====================================================================================
def _checkNcellsNptsPerProc(ts, NP, isAtCenter=False):
    NPTS   = numpy.zeros(NP, dtype=Internal.E_NpyInt)
    NCELLS = numpy.zeros(NP, dtype=Internal.E_NpyInt)
    ##Done by zone for flexibility
    for z in Internal.getZones(ts):
        proc_num        = Cmpi.getProc(z)        
        NPTS[proc_num] += C.getNPts(z)
        if not isAtCenter:
            NCELLS[proc_num] += C.getNCells(z)
        else:
            NCELLS[proc_num]  = NPTS[proc_num]

    NPTS     = Cmpi.allreduce(NPTS  ,op=Cmpi.SUM)
    NCELLS   = Cmpi.allreduce(NCELLS,op=Cmpi.SUM)
    NptsTot  = numpy.sum(NPTS)
    NcellsTot= numpy.sum(NCELLS)
    ncells_percent=[]
    if Cmpi.rank == 0:
        for i in range(NP):
            ncells_percent.append(NCELLS[i]/NcellsTot*100.)
            if isAtCenter: print('Rank {} has {} cells & {} % of cells'.format(i,int(NCELLS[i]),round(ncells_percent[i],2)))
            else:          print('Rank {} has {} points & {} cells & {} % of cells'.format(i,int(NPTS[i]),int(NCELLS[i]),round(ncells_percent[i],2)))

        if isAtCenter:   print('All points: {} million cells'.format(NcellsTot/1.e6))
        else:            print('All points: {} million points & {} million cells'.format(NptsTot/1.e6,NcellsTot/1.e6))
        print('Range of % of cells: {} - {}'.format(round(min(ncells_percent),2),round(max(ncells_percent),2)))

    return None


#====================================================================================
##This interpolation complements the function _initialize_t in the class IBM.
##If the meshes are large the function in the class IBM can fail due to memory limits
##and the whole mesh generation fails - i.e. have to rerun to get a t and tc.
##This complementatry function does the solution extraction a posterior guaranteeing that
##a mesh generation is performed and if needed it can use more MPI tasks to perform the interpolation
##yielding a result.
#====================================================================================
def interpolateSolutionCoarse2Fine(tCoarse, tFine, NPprep, NPinterp):
    """Interpolate the solution from a coarse IBM mesh to fine IBM mesh. Useful when doing grid convergence studies. This can speed up the convergence rate of the solution on a finer mesh.
        Usage: interpolateSolutionCoarse2Fine(tCoarse,tFine,NPprep,NPinterp)"""
    if isinstance(tCoarse, str):
        h = Filter.Handle(tCoarse)
        if NPinterp == NPprep:
            tCoarse = h.loadFromProc()
        else:
            tCoarse = h.loadAndDistribute()

    if isinstance(tFine, str):
        h         = Filter.Handle(tFine)
        if NPinterp == NPprep:
            tFine = h.loadFromProc()
        else:
            tFine = h.loadAndDistribute()

    C._rmVars(tCoarse,['centers:cellN','centers:TurbulentDistance'])
    if C.isNamePresent(tCoarse, 'centers:MomentumX')==-1 and C.isNamePresent(tFine, 'centers:MomentumX')>-1:
        C._rmVars(tFine,['centers:MomentumX','centers:MomentumY','centers:MomentumZ',
                             'centers:EnergyStagnationDensity','centers:TurbulentSANuTildeDensity'])
    ## Extract mesh from COARSE to CURRENT
    tFine = Pmpi.extractMesh(tCoarse, tFine)
    tFine = C.node2Center(tFine,'FlowSolution')
    Internal._rmNodesByName(tFine,'FlowSolution')
    return tFine



# IN: maillage surfacique + reference State + snears
#================================================================================
# IBM prepare
# IN: t_case: fichier ou arbre body
# OUT: t_out, tc_out: fichier ou arbres de sorties
# snears         : liste des snear, mais c'est mieux si ils sont dans t_case
# dfar, dfarList : liste des dfars, mais c'est mieux si ils sont dans t_case
# tbox           : arbre de raffinement
# check          : si true, fait des sorties
# frontType=1,2,3: type de front
# expand=1,2,3   : sorte d'expand des niveaux (1:classque,2:minimum,3:deux niveaux)
# tinit          : arbre de champ d'avant pour la reprise
# smoothing      : smooth the front during the front 2 specific treatment in the cases of local refinements
# balancing      : balance the entire distribution after the octree generation, useful for symetries
# distrib        : new distribution at the end of prepare1
#===================================================================================================================

class IBM_Input:
    def __init__(self):
        self.IBCType                   = 1
        self.Lref                      = 1.
        self.NP                        = 0
        self.NPas                      = 10
        self.Ntranche                  = 1
        self.balancing                 = False
        self.blankingF42               = False
        self.cartesian                 = True
        self.check                     = False
        self.check_snear               = False
        self.cleanCellN                = True
        self.generateCartesianMeshOnly = False
        self.tbOneOver                 = None
        self.correctionMultiCorpsF42   = False
        self.dfars                     = 10.
        self.dfarDir                   = 0
        self.dictNz                    = None
        self.distrib                   = True
        self.dz                        = 0.01 
        self.expand                    = 3
        self.nature                    = 1
        self.ext                       = 2
        self.extrusion                 = None
        self.format                    = 'single'
        self.frontType                 = 1
        self.height_in                 = -1.0
        self.ific                      = 2
        self.initWithBBox              = -1.
        self.interpDataType            = 0
        self.isCartesianExtrude        = False
        self.isWireModel               = False
        self.optimized                 = 1
        self.order                     = 2
        self.recomputeDist             = False
        self.redistribute              = False
        self.smoothing                 = False
        self.snears                    = 0.01
        self.snearsf                   = None
        self.span                      = 1 
        self.t_in                      = None
        self.tbox                      = None
        self.tinit                     = None
        self.to                        = None
        self.twoFronts                 = False
        self.vmin                      = 21
        self.wallAdapt                 = None
        self.yplus                     = 100.
        self.yplusAdapt                = 100.

class IBM(Common):
    """Prepare IBM for FastS"""
    def __init__(self, format=None, numb=None, numz=None):
        Common.__init__(self, format, numb, numz)
        self.__version__ = "0.0"
        self.authors = ["defi@onera.fr"]
        self.input_var = IBM_Input()

        self.DEPTH                 = None
        self.RefState              = None
        self.Reynolds              = None
        self.cptBody               = None
        self.dimPb                 = None
        self.isOrthoProjectFirst   = None
        self.model                 = None
        self.nbZonesIBC            = None    
        self.procDict              = None
        self.rank                  = None
        self.refstate              = None
        self.res                   = None   
        self.res2                  = None  
        self.tbFilament            = None
        self.tbFilamentWMM         = None
        self.filamentBases         = None
        self.filamentBasesWMM      = None
        self.isFilamentOnly        = None
        self.tbbc                  = None
        self.tbsave                = None
        self.zonesRIBC             = None
        self.ibctypes              = None

    #================================================================================
    # Selection/determination of tb for closed solid & filament
    #================================================================================    
    def _determineClosedSolidFilament__(self,tb):
        ## General case where only a closeSolid is in tb
        ## or tb only has a filament

        self.filamentBases    = []
        len_tb = len(Internal.getBases(tb))
        for b in Internal.getBases(tb):
            if "IBCFil" in b[0]:self.filamentBases.append(b[0])

        self.isFilamentOnly=False
        if len(self.filamentBases) == len_tb:           
            self.isFilamentOnly=True
        self.isOrthoProjectFirst = self.isFilamentOnly
        ## if tb has both a closed solid and filaments

        self.tbFilament = Internal.getBases(tb)
        if not self.isFilamentOnly:
            tbFilament = []
            for b in self.filamentBases:
                node_local = Internal.getNodeFromNameAndType(tb, b, 'CGNSBase_t')
                tbFilament.append(node_local)
                Internal._rmNode(tb,node_local)     
                self.isOrthoProjectFirst = True
            self.tbFilament = C.newPyTree(tbFilament);        
        return None


    #================================================================================
    # Cartesian grid generator
    #================================================================================
    def generateCartesian__(self, tb, ext_local=3, to=None):    
        # list of dfars
        if self.input_var.dfars == []:
            zones = Internal.getZones(tb)
            self.input_var.dfars = [10.]*len(zones)
            for c, z in enumerate(zones):
                n = Internal.getNodeFromName2(z, 'dfar')
                if n is not None: self.input_var.dfars[c] = Internal.getValue(n)*1.

         # refinementSurfFile: surface meshes describing refinement zones
        if self.input_var.tbox is not None:
            if isinstance(self.input_var.tbox, str): self.input_var.tbox = C.convertFile2PyTree(self.input_var.tbox)
            else: self.input_var.tbox = self.input_var.tbox
            if self.input_var.snearsf is None:
                self.input_var.snearsf = []
                zones = Internal.getZones(self.input_var.tbox)
                for z in zones:
                    sn = Internal.getNodeFromName2(z, 'snear')
                    if sn is not None: self.input_var.snearsf.append(Internal.getValue(sn))
                    else: self.input_var.snearsf.append(1.)

        fileout = None
        if self.input_var.check: fileout = 'octree.cgns'
        # Octree identical on all procs
        test.printMem('>>> Octree unstruct [start]')
        if self.input_var.to is not None:
            if isinstance(self.input_var.to, str):
                o = C.convertFile2PyTree(self.input_var.to)
                o = Internal.getZones(o)[0]
            else:
                o = Internal.getZones(self.input_var.to)[0]
            parento = None
        else:
            o = G_IBM.buildOctree(tb, snears=self.input_var.snears, snearFactor=1., dfars=self.input_var.dfars,
                                  to=self.input_var.to, tbox=self.input_var.tbox, snearsf=self.input_var.snearsf,
                                  dimPb=self.dimPb, vmin=self.input_var.vmin,
                                  expand=self.input_var.expand, dfarDir=self.input_var.dfarDir)

        if self.rank==0 and self.input_var.check: C.convertPyTree2File(o, fileout)
        # build parent octree 3 levels higher
        # returns a list of 4 octants of the parent octree in 2D and 8 in 3D
        parento = G_IBM.buildParentOctrees__(o, tb, snears=self.input_var.snears, snearFactor=4., dfars=self.input_var.dfars, to=self.input_var.to,
                                             tbox=self.input_var.tbox, snearsf=self.input_var.snearsf,
                                             dimPb=self.dimPb, vmin=self.input_var.vmin)
        test.printMem(">>> Octree unstruct [end]")

        if self.input_var.check_snear: exit()

        # Split octree
        test.printMem(">>> Octree unstruct split [start]")
        bb = G.bbox(o)
        self.NP = Cmpi.size
        if self.NP == 1: p = Internal.copyRef(o) # keep reference
        else: p = T.splitNParts(o, N=self.NP, recoverBC=False)[self.rank]
        del o
        test.printMem(">>> Octree unstruct split [end]")

        # fill vmin + merge in parallel
        test.printMem(">>> Octree struct [start]")
        res = G_IBM.octree2StructLoc__(p, vmin=self.input_var.vmin, ext=-1, optimized=0, parento=parento, sizeMax=1000000,tbOneOver=self.input_var.tbOneOver)
        del p
        if parento is not None:
            for po in parento: del po
        t = C.newPyTree(['CARTESIAN', res])
        zones = Internal.getZones(t)
        for z in zones: z[0] = z[0]+'X%d'%self.rank
        Cmpi._setProc(t, self.rank)

        C._addState(t, 'EquationDimension', self.dimPb)
        test.printMem(">>> Octree struct [end]")

        # Add xzones for ext       
        test.printMem(">>> extended cart grids [start]")
        tbb = Cmpi.createBBoxTree(t)
        interDict = X.getIntersectingDomains(tbb)
        graph = Cmpi.computeGraph(tbb, type='bbox', intersectionsDict=interDict, reduction=False)
        del tbb
        Cmpi._addXZones(t, graph, variables=[], cartesian=True)

        if self.input_var.generateCartesianMeshOnly:
            Cmpi.convertPyTree2File(t,'CartMesh.cgns')
            if Cmpi.rank == 0: print("Wrote CartMesh.cgns...Exiting Prep...Bye")
            exit()

        #Turn Cartesian grid into a rectilinear grid
        test.printMem(">>> cart grids --> rectilinear grids [start]")        
        listDone = []
        if self.input_var.tbOneOver:
            tbb = G.BB(t)
            if self.dimPb==2:
                T._addkplane(tbb)
                T._contract(tbb, (0,0,0), (1,0,0), (0,1,0), 0.01)

            ##RECTILINEAR REGION
            ##Select regions that need to be coarsened
            tbbB            = G.BB(self.input_var.tbOneOver)                
            interDict_scale = X.getIntersectingDomains(tbbB, tbb)
            ##Avoid a zone to be coarsened twice
            for i in interDict_scale:
                (b,btmp) = Internal.getParentOfNode(self.input_var.tbOneOver,Internal.getNodeByName(self.input_var.tbOneOver,i))
                b        = Internal.getNodeByName(b,".Solver#define")
                oneoverX = int(Internal.getNodeByName(b, 'dirx')[1])
                oneoverY = int(Internal.getNodeByName(b, 'diry')[1])
                oneoverZ = int(Internal.getNodeByName(b, 'dirz')[1])
                for z in interDict_scale[i]:
                    if z not in listDone:
                        zLocal = Internal.getNodeFromName(t,z)
                        T._oneovern(zLocal, (oneoverX,oneoverY,oneoverZ));
                        listDone.append(z)
        test.printMem(">>> cart grids --> rectilinear grids [end]")        

        zones  = Internal.getZones(t)
        coords = C.getFields(Internal.__GridCoordinates__, zones, api=2)

        #print("extension LOCAL=", ext_local)
        #coords, rinds = Generator.extendCartGrids(coords, ext=ext_local, optimized=1, extBnd=0)
        coords, rinds = Generator.extendCartGrids(coords, ext=ext_local, optimized=self.input_var.optimized, extBnd=0)

        C.setFields(coords, zones, 'nodes')
        for noz in range(len(zones)):
            #print("Rind", rinds[noz], "No zone=", noz)
            Internal.newRind(value=rinds[noz], parent=zones[noz])
        Cmpi._rmXZones(t)
        coords = None; zones = None
        test.printMem(">>> extended cart grids (after rmXZones) [end]")        

        G_IBM._addBCOverlaps(t, bbox=bb)
        G_IBM._addExternalBCs(t, bbox=bb, dimPb=self.dimPb)

        if self.dimPb == 2:
            self.input_var.dz = 0.01
            T._addkplane(t)
            T._contract(t, (0,0,0), (1,0,0), (0,1,0), self.input_var.dz)

        # ReferenceState
        C._addState(t, state=self.refstate)
        C._addState(t, 'GoverningEquations', self.model)
        C._addState(t, 'EquationDimension', self.dimPb)
        return t


    def _distance2wallCalc__(self,t,tb):
        test.printMem(">>> Wall distance [start]")        
        self.cptBody                =0
        FSC = Internal.getNodeFromType(t,"FlowSolution_t")
        if (FSC is None or Internal.getNodeFromName(FSC,'TurbulentDistance') is None) and self.input_var.extrusion is None:
            C._initVars(t,'{centers:TurbulentDistance}=1e06')
            if self.dimPb == 2:
                tb2 = C.initVars(tb, 'CoordinateZ', self.input_var.dz*0.5)
                self.tbsave = tb2

                dist2wallNearBody__(t, tb2, type='ortho', signed=0, dim=self.dimPb, loc='centers',model=self.model)

                if self.filamentBases:
                    C._initVars(t,'{centers:TurbulentDistanceSolid}={centers:TurbulentDistance}')

                    tb2             = C.initVars(self.tbFilament, 'CoordinateZ', self.input_var.dz*0.5)                   
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
                        self.tbsave   = Internal.merge([self.tbsave,tbFilamentnoWMM])
                        C._initVars(t,'{centers:TurbulentDistance}=1e06')
                        dist2wallNearBody__(t, tbFilamentnoWMM, type='ortho', signed=0, dim=self.dimPb, loc='centers',model=self.model)
                        C._initVars(t,'{centers:TurbulentDistanceFilament}={centers:TurbulentDistance}')
                        C._initVars(t,'{centers:TurbulentDistance}=minimum({centers:TurbulentDistanceSolid},{centers:TurbulentDistanceFilament})')


                    if tbFilamentWMM:
                        C._initVars(t,'{centers:TurbulentDistance}=1e06')
                        dist2wallNearBody__(t, tbFilamentWMM, type='ortho', signed=0, dim=self.dimPb, loc='centers',model=self.model)
                        C._initVars(t,'{centers:TurbulentDistanceFilamentWMM}={centers:TurbulentDistance}')
                        C._initVars(t,'{centers:TurbulentDistance}=minimum({centers:TurbulentDistanceSolid},{centers:TurbulentDistanceFilamentWMM})')
                        if tbFilamentnoWMM:
                            C._initVars(t,'{centers:TurbulentDistance}=minimum({centers:TurbulentDistance},{centers:TurbulentDistanceFilamentWMM})')

                    C._initVars(t,"{centers:TurbulentDistanceSolid}=({centers:TurbulentDistanceSolid}>1e03)*0+({centers:TurbulentDistanceSolid}<1e03)*{centers:TurbulentDistanceSolid}")
                    C._initVars(t,"{centers:TurbulentDistanceFilament}=({centers:TurbulentDistanceFilament}>1e03)*0+({centers:TurbulentDistanceFilament}<1e03)*{centers:TurbulentDistanceFilament}")
                    C._initVars(t,"{centers:TurbulentDistanceFilamentWMM}=({centers:TurbulentDistanceFilamentWMM}>1e03)*0+({centers:TurbulentDistanceFilamentWMM}<1e03)*{centers:TurbulentDistanceFilamentWMM}")

            else:
                dist2wallNearBody__(t, tb, type='ortho', signed=0, dim=self.dimPb, loc='centers',model=self.model)

                self.tbsave = tb
                if self.filamentBases:
                    C._initVars(t,'{centers:TurbulentDistanceSolid}={centers:TurbulentDistance}')

                    tbFilamentnoWMM = []
                    tbFilamentWMM   = []
                    for z in Internal.getZones(self.tbFilament):
                        soldef  = Internal.getNodeFromName(z,'.Solver#define')
                        ibctype = Internal.getNodeFromName(soldef,'ibctype')
                        if Internal.getValue(ibctype) != "wiremodel":
                            tbFilamentnoWMM.append(z)
                        else:
                            tbFilamentWMM.append(z)

                    if tbFilamentnoWMM:
                        tbFilamentnoWMM = C.newPyTree(['BasenoWMM', tbFilamentnoWMM])
                        self.tbsave   = Internal.merge([self.tbsave,tbFilamentnoWMM])
                        C._initVars(t,'{centers:TurbulentDistance}=1e06')
                        dist2wallNearBody__(t, tbFilamentnoWMM, type='ortho', signed=0, dim=self.dimPb, loc='centers',model=self.model)
                        C._initVars(t,'{centers:TurbulentDistanceFilament}={centers:TurbulentDistance}')
                        C._initVars(t,'{centers:TurbulentDistance}=minimum({centers:TurbulentDistanceSolid},{centers:TurbulentDistanceFilament})')


                    if tbFilamentWMM:
                        C._initVars(t,'{centers:TurbulentDistance}=1e06')
                        dist2wallNearBody__(t, tbFilamentWMM, type='ortho', signed=0, dim=self.dimPb, loc='centers',model=self.model)
                        C._initVars(t,'{centers:TurbulentDistanceFilamentWMM}={centers:TurbulentDistance}')
                        C._initVars(t,'{centers:TurbulentDistance}=minimum({centers:TurbulentDistanceSolid},{centers:TurbulentDistanceFilamentWMM})')

                        if tbFilamentnoWMM:
                            C._initVars(t,'{centers:TurbulentDistance}=minimum({centers:TurbulentDistance},{centers:TurbulentDistanceFilamentWMM})')

                    C._initVars(t,"{centers:TurbulentDistanceSolid}=({centers:TurbulentDistanceSolid}>1e03)*0+({centers:TurbulentDistanceSolid}<1e03)*{centers:TurbulentDistanceSolid}")
                    C._initVars(t,"{centers:TurbulentDistanceFilament}=({centers:TurbulentDistanceFilament}>1e03)*0+({centers:TurbulentDistanceFilament}<1e03)*{centers:TurbulentDistanceFilament}")
                    C._initVars(t,"{centers:TurbulentDistanceFilamentWMM}=({centers:TurbulentDistanceFilamentWMM}>1e03)*0+({centers:TurbulentDistanceFilamentWMM}<1e03)*{centers:TurbulentDistanceFilamentWMM}")

            C._initVars(t,"{centers:TurbulentDistance}=({centers:TurbulentDistance}>1e03)*0+({centers:TurbulentDistance}<1e03)*{centers:TurbulentDistance}")


            # Compute turbulentdistance wrt each body that is not a sym plan
            if self.input_var.correctionMultiCorpsF42 and self.input_var.frontType==42:
                test.printMem(">>> Individual wall distance [start]")
                # Keep track of the general turbulentDistance
                C._initVars(t,'{centers:TurbulentDistance_ori}={centers:TurbulentDistance}')

                if self.input_var.yplus > 0.:
                    shiftDist = G_IBM_Height.computeModelisationHeight(Re=self.Reynolds, yplus=self.input_var.yplus, L=self.input_var.Lref)
                else:
                    self.input_var.snears    = Internal.getNodesFromName(tb, 'snear')
                    h         = max(self.input_var.snears, key=lambda x: x[1])[1]
                    shiftDist = G_IBM_Height.computeBestModelisationHeight(Re=self.Reynolds, h=h) # meilleur compromis entre hauteur entre le snear et la hauteur de modelisation

                if self.input_var.height_in>0.:
                    if shiftDist>self.input_var.height_in: shiftDist=self.input_var.height_in

                for z in Internal.getZones(t):
                    self.cptBody = 1
                    if self.dimPb == 3: tb2 = tb
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
                                DTW._distance2Walls(z, body, type='ortho', signed=0, dim=self.dimPb, loc='centers')
                                C._initVars(z,'{centers:TurbulentDistance_body%i={centers:TurbulentDistance}'%self.cptBody)
                            else:
                                C._initVars(z,'{centers:TurbulentDistance_body%i=1000'%self.cptBody)
                            self.cptBody += 1
                    if self.dimPb == 3: del tb2

                C._initVars(t,'{centers:TurbulentDistance}={centers:TurbulentDistance_ori}')
                C._rmVars(t,['centers:TurbulentDistance_ori'])

            #
            if self.dimPb == 2 and self.input_var.cleanCellN == False: C._initVars(t, '{centers:TurbulentDistanceAllBC}={centers:TurbulentDistance}')

        else:
            C._initVars(t, '{centers:TurbulentDistance}={centers:TurbulentDistanceAllBC}')
        test.printMem(">>> Wall distance [end]")
        return None


    def blanking__(self, t, tb):
        test.printMem(">>> Blanking [start]")

        snear_min = 10.e10
        for z in Internal.getZones(tb):
            sdd = Internal.getNodeFromName1(z, ".Solver#define")
            if sdd is not None:
                snearl = Internal.getNodeFromName1(sdd, "snear")
                if snearl is not None:
                    snearl = Internal.getValue(snearl)
                if snearl is not None: snear_min = min(snear_min, snearl)
        snear_min = Cmpi.allreduce(snear_min, op=Cmpi.MIN)

        if self.input_var.extrusion is None:
            if not self.isFilamentOnly:
                t = X_IBM.blankByIBCBodies(t, tb, 'centers', self.dimPb)
            if self.dimPb == 2 and self.input_var.cleanCellN == False:
                C._initVars(t, '{centers:cellNIBC_blank}={centers:cellN}')
        else:
            C._initVars(t, '{centers:cellN}={centers:cellNIBC_blank}')

        C._initVars(t, '{centers:cellNIBC}={centers:cellN}')

        if not self.isFilamentOnly: X_IBM._signDistance(t)

        if self.input_var.extrusion is not None:
            C._initVars(t,'{centers:cellN}={centers:cellNIBC_blank}')
        else:
            C._initVars(t,'{centers:cellN}={centers:cellNIBC}')

        if (self.filamentBases or self.isFilamentOnly) and self.input_var.extrusion is None:
            C._initVars(t,'{centers:cellNFil}={centers:cellN}')
            C._initVars(t,'{centers:cellNFilWMM}={centers:cellN}')

        # determination des pts IBC
        if self.input_var.extrusion is None:
            if self.input_var.frontType != 42:
                if self.input_var.IBCType == -1: X._setHoleInterpolatedPoints(t,depth=-self.DEPTH,dir=0,loc='centers',cellNName='cellN',addGC=False)
                elif self.input_var.IBCType == 1:
                    X._setHoleInterpolatedPoints(t, depth=1, dir=1, loc='centers', cellNName='cellN', addGC=False) # pour les gradients
                    if self.input_var.frontType < 2:
                        X._setHoleInterpolatedPoints(t,depth=self.DEPTH,dir=0,loc='centers',cellNName='cellN',addGC=False)
                    else:
                        DEPTHL = self.DEPTH+1
                        X._setHoleInterpolatedPoints(t,depth=DEPTHL,dir=0,loc='centers',cellNName='cellN',addGC=False)
                        #cree des pts extrapoles supplementaires
                        # _blankClosestTargetCells(t,cellNName='cellN', depth=DEPTHL)
                else:
                    raise ValueError('prepareIBMData: not valid IBCType. Check model.')

            else:
                # F42: tracking of IB points using distance information
                # the classical algorithm (front 1) is first used to ensure a minimum of two rows of target points around the geometry
                C._initVars(t,'{centers:cellNMin}={centers:cellNIBC}')
                if self.input_var.IBCType == -1: X._setHoleInterpolatedPoints(t,depth=-self.DEPTH,dir=0,loc='centers',cellNName='cellNMin',addGC=False)
                elif self.input_var.IBCType == 1: X._setHoleInterpolatedPoints(t,depth=1,dir=1,loc='centers',cellNName='cellNMin',addGC=False) # pour les gradients
                X._setHoleInterpolatedPoints(t,depth=self.DEPTH,dir=0,loc='centers',cellNName='cellNMin',addGC=False)

                for z in Internal.getZones(t):
                    h = abs(C.getValue(z,'CoordinateX',0)-C.getValue(z,'CoordinateX',1))
                    if self.input_var.yplus > 0.:
                        height = G_IBM_Height.computeModelisationHeight(Re=self.Reynolds, yplus=self.input_var.yplus, L=self.input_var.Lref)
                    else:
                        height = G_IBM_Height.computeBestModelisationHeight(Re=self.Reynolds, h=h) # meilleur compromis entre hauteur entre le snear et la hauteur de modelisation
                        self.input_var.yplus  = G_IBM_Height.computeYplus(Re=self.Reynolds, height=height, L=self.input_var.Lref)
                    if self.input_var.height_in>0.:
                        if height>self.input_var.height_in:
                            height=self.input_var.height_in
                            #print("Snear min (SM) = %g || Wall Modeling Height (WMH) = %g || WMH/SM = %g"%(snear_min,height,height/snear_min))
                    C._initVars(z,'{centers:cellN}=({centers:TurbulentDistance}>%20.16g)+(2*({centers:TurbulentDistance}<=%20.16g)*({centers:TurbulentDistance}>0))'%(height,height))

                    if self.input_var.correctionMultiCorpsF42:
                        # Prevent different body modeling from overlapping -> good projection of image points in the wall normal direction

                        epsilon_dist = 2*(abs(C.getValue(z,'CoordinateX',1)-C.getValue(z,'CoordinateX',0)))
                        max_dist     = 2*0.1*self.input_var.Lref

                        # Try to find the best route between two adjacent bodies by finding optimal iso distances
                        def correctionMultiCorps(cellN, cellNF):
                            if cellN == 2 and cellNF == 2: return 1
                            return cellN

                        def findIsoFront(cellNFront, Dist_1, Dist_2):
                            if Dist_1 < max_dist and Dist_2 < max_dist:
                                if abs(Dist_1-Dist_2) < epsilon_dist: return 2
                            return max(cellNFront,1)

                        for i in range(1, self.cptBody):
                            for j in range(1, self.cptBody):
                                if j != i:
                                    C._initVars(z,'centers:cellNFrontIso', findIsoFront, ['centers:cellNFrontIso', 'centers:TurbulentDistance_body%i'%i, 'centers:TurbulentDistance_body%i'%j])

                        C._initVars(z,'centers:cellN', correctionMultiCorps, ['centers:cellN', 'centers:cellNFrontIso'])

                        for i in range(1,self.cptBody):
                            C._rmVars(z,['centers:cellN_body%i'%i, 'centers:TurbulentDistance_body%i'%i])

                if self.input_var.wallAdapt is not None:
                    # Use previous computation to adapt the positioning of IB points around the geometry (impose y+PC <= y+ref)
                    # Warning: the wallAdapt file has to be obtained with TIBM.createWallAdapt(tc)
                    C._initVars(t,'{centers:yplus}=100000.')
                    w = C.convertFile2PyTree(self.input_var.wallAdapt)
                    total = len(Internal.getZones(t))
                    cpt = 1
                    for z in Internal.getZones(t):
                        print("{} / {}".format(cpt, total))
                        cellN = Internal.getNodeFromName(z,'cellN')[1]
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

                    C._initVars(t,'{centers:cellN}=({centers:cellN}>0) * ( (({centers:cellN}) * ({centers:yplus}<=%20.16g)) + ({centers:yplus}>%20.16g) )'%(self.input_var.yplus,self.input_var.yplus))

                # final security gate, we ensure that we have at least to layers of target points
                C._initVars(t, '{centers:cellN} = maximum({centers:cellN}, {centers:cellNMin})')
                C._rmVars(t,['centers:yplus', 'centers:cellNMin'])

                # propagate max yplus between procs
                yplus                = numpy.array([float(self.input_var.yplus)])
                self.input_var.yplus = Cmpi.allreduce(yplus, op=Cmpi.MAX)[0]

                # Only keep the layer of target points useful for solver iterations, particularly useful in 3D
                if self.input_var.blankingF42: X._maximizeBlankedCells(t, depth=2, addGC=False)

            if self.dimPb == 2 and self.input_var.cleanCellN == False: C._initVars(t, '{centers:cellNIBC_hole}={centers:cellN}')

        else:  # extrusion
            C._initVars(t, '{centers:cellN}={centers:cellNIBC_hole}')

        if (self.filamentBases or self.isFilamentOnly) and self.input_var.extrusion is None:
            if self.isFilamentOnly:
                C._initVars(t,'{centers:TurbulentDistanceFilament}={centers:TurbulentDistance}')
                maxy=C.getMaxValue(tb, ['CoordinateY']);
                miny=C.getMinValue(tb, ['CoordinateY']);
            if self.filamentBases:
                tbFilamentWMM   = []
                for z in Internal.getZones(self.tbFilament):
                    soldef  = Internal.getNodeFromName(z,'.Solver#define')
                    ibctype = Internal.getNodeFromName(soldef,'ibctype')
                    if Internal.getValue(ibctype) == "wiremodel":
                        tbFilamentWMM.append(z)
                maxy=C.getMaxValue(tbFilamentWMM, ['CoordinateY']);
                miny=C.getMinValue(tbFilamentWMM, ['CoordinateY']);

                if self.dimPb == 3:
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
                    for k in range(sh[2]):
                        for j in range(sh[1]):
                            for i in range(sh[0]):
                                valy=0.5*(ycord[i,j,k]+ycord[i,j+1,k])
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
                    if self.dimPb == 3:
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
            if self.filamentBases: C._rmVars(t,['centers:TurbulentDistanceSolid'])

        if self.input_var.extrusion is None:
            if not self.isFilamentOnly:X_IBM._removeBlankedGrids(t, loc='centers')


        localDistrib = 0
        if not Internal.getZones(t):
            localDistrib=1
            print("Proc %d  has 0 zones"%Cmpi.rank)

        localDistrib = Cmpi.allreduce(localDistrib  ,op=Cmpi.SUM)
        if localDistrib>0 or (self.input_var.tbOneOver and Cmpi.size>1):
            import Distributor2.Mpi as D2mpi
            if Cmpi.rank==0 and localDistrib>0:print('>>> Distribute (atleast one Proc has 0 Zones [start])')
            tcs    = Cmpi.allgatherTree(Cmpi.convert2SkeletonTree(t))
            stats  = D2._distribute(tcs, self.NP, algorithm='graph')
            D2._copyDistribution(t , tcs)
            D2mpi._redispatch(t)
            del tcs
            if Cmpi.rank==0 and localDistrib>0:print('>>> Distribute (atleast one Proc has 0 Zones [end])')
            if self.input_var.tbOneOver and Cmpi.size>1: self.input_var.redistribute = True

        test.printMem(">>> Blanking [end]")

        return t


    def setInterpDataAndSetInterpTransfer__(self, t):

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
        self.tbbc = Cmpi.createBBoxTree(tc)
        interDict = X.getIntersectingDomains(self.tbbc)
        graph     = Cmpi.computeGraph(self.tbbc, type='bbox', intersectionsDict=interDict, reduction=False)
        Cmpi._addXZones(tc, graph, variables=['cellN'], cartesian=self.input_var.cartesian)
        test.printMem(">>> Interpdata [after addXZones]")

        self.procDict = Cmpi.getProcDict(tc)
        datas = {}


        if self.input_var.order != 2:
            #modif Ivan pour limiter recouvrement en ordre 5
            ##  etape1 : interp ordre2 nature 1 entre grille de meme niveau
            ##  etape2 : modif cellN=0 pour les points interpole dans etape 1
            ##  etape3 : interp ordre5 nature 0 entre grille de niveau N (Receveur) et N-1, N+1 (Donneur)

            ## determine dx=dy for each zone & store per zone
            levelZone_loc={}
            hmin = 1.e30
            for z in Internal.getZones(t):
                h = abs(C.getValue(z,'CoordinateX',0)-C.getValue(z,'CoordinateX',1))
                #print("dx=",h)
                levelZone_loc[z[0]]=h
                if h < hmin : hmin = h

            ## go from dx to dx/dx_min
            Nlevels=1
            for i in levelZone_loc:
                levelZone_loc[i]= math.log( int(levelZone_loc[i]/hmin + 0.00000001)  , 2)
                if levelZone_loc[i] +1  > Nlevels : Nlevels = int(levelZone_loc[i]) +1

            ## partage des info level en mpi
            levelZone = Cmpi.allgather(levelZone_loc) 

            #
            #etape1: calcul interpolation entre grille de meme niveau (ordre2, nature1, pas d'extrap)a
            #
            for level in range(Nlevels):
                print("Interp level=",level)
                #filtre les zones par niveau de resolution
                zones=[]
                for z in Internal.getZones(t):
                    if levelZone[z[0]]==level:
                        zones.append(z)
                        print("level=", level,"zone In", z[0])

                for zrcv in zones:
                    zrname = zrcv[0]

                    dnrZones = []
                    for zdname in interDict[zrname]:
                        zd = Internal.getNodeFromName2(tc, zdname)
                        if levelZone[zd[0]]==level: dnrZones.append(zd)
                    if dnrZones:
                        X._setInterpData(zrcv, dnrZones, nature=1, penalty=1, extrap=0,loc='centers', storage='inverse',
                                         sameName=1, interpDataType=self.input_var.interpDataType, order=2, itype='chimera')
                    #
                    #etape 2 : modif CellN=0 pour les points interpoles a l'etape ci dessus et on vire info orphelin
                    #
                    for zd in dnrZones:
                        zdname = zd[0]
                        destProc = self.procDict[zdname]

                        IDs = []
                        for i in zd[2]:
                            modif=0
                            if i[0][0:2] == 'ID':
                                Internal.rmNodesByName(i,'OrphanPointList')

                                if Internal.getValue(i)==zrname: 
                                    IDs.append(i)
                                    PtlistD = Internal.getNodeFromName(i, "PointListDonor")[1]

                                    for zR in zones:
                                        if zR[0]==zrname: break
                                    print("zone receveuse", zrname, 'zD=', zdname)
                                    cellN_R = Internal.getNodeFromName(zR, "cellN")[1]
                                    sh_R    = numpy.shape(cellN_R)

                                    if len(sh_R)==2:
                                        for l in range(numpy.size(PtlistD)):
                                            jR = (PtlistD[l])//sh_R[0]
                                            iR =  PtlistD[l] -jR*sh_R[0]
                                            cellN_R[iR,jR]=0
                                    else:
                                        for l in range(numpy.size(PtlistD)):
                                            kR =  PtlistD[l]//(sh_R[0]*sh_R[1])
                                            jR = (PtlistD[l] -kR*sh_R[0]*sh_R[1])//sh_R[0]
                                            iR =  PtlistD[l] -kR*sh_R[0]*sh_R[1] -jR*sh_R[0]
                                            cellN_R[iR,jR,kR]=0

                        if IDs != []:
                            if destProc == self.rank:
                                zD = Internal.getNodeFromName2(tc, zdname)
                                zD[2] += IDs
                            else:
                                if destProc not in datas: datas[destProc] = [[zdname,IDs]]
                                else: datas[destProc].append([zdname,IDs])
                        else:
                            if destProc not in datas: datas[destProc] = []
            '''
           #t4 = X.getOversetInfo(t, tc, loc='center',type='extrapolated')
           #t = X.getOversetInfo(t, tc, loc='center',type='orphan')
           #C.convertPyTree2File(t,'t_prep.cgns')
           #C.convertPyTree2File(tc,'tc2_prep.cgns')
           #stop
           '''

            #
            #etape 3: calcul interpolation entre grille de niveau N (Receveur) et N+1, N-1 (donneurs)
            #
            for level in range(Nlevels):
                print("Interp level=",level)
                #filtre les zones par niveau de resolution
                zones=[]
                for z in Internal.getZones(t):
                    if levelZone[z[0]]==level:
                        zones.append(z)
                        print("level=", level,"zone In", z[0])

                for zrcv in zones:
                    zrname = zrcv[0]
                    print("ZoneR=", zrname)
                    dnrZones = []
                    for zdname in interDict[zrname]:
                        zd = Internal.getNodeFromName2(tc, zdname)
                        if levelZone[zd[0]]==level+1 or  levelZone[zd[0]]==level-1 : 
                            dnrZones.append(zd)
                            print("ZoneD=", zd[0], levelZone[zd[0]])
                    if dnrZones:
                        X._setInterpData(zrcv, dnrZones, nature=self.input_var.nature, penalty=1, loc='centers', storage='inverse',
                                         sameName=1, interpDataType=self.input_var.interpDataType, order=self.input_var.order, itype='chimera')
                    for zd in dnrZones:
                        zdname = zd[0]
                        destProc = self.procDict[zdname]

                        IDs = []
                        for i in zd[2]:
                            if i[0][0:2] == 'ID':
                                if Internal.getValue(i)==zrname: 
                                    IDs.append(i)
                        if IDs != []:
                            if destProc == self.rank:
                                zD = Internal.getNodeFromName2(tc, zdname)
                                zD[2] += IDs
                            else:
                                if destProc not in datas: datas[destProc] = [[zdname,IDs]]
                                else: datas[destProc].append([zdname,IDs])
                        else:
                            if destProc not in datas: datas[destProc] = []
            #Fin modif Ivan
        else:
            for zrcv in Internal.getZones(t):
                zrname = zrcv[0]
                dnrZones = []
                for zdname in interDict[zrname]:
                    zd = Internal.getNodeFromName2(tc, zdname)
                    dnrZones.append(zd)
                if dnrZones:
                    X._setInterpData(zrcv, dnrZones, nature=self.input_var.nature, penalty=1, loc='centers', storage='inverse',
                                  sameName=1, interpDataType=self.input_var.interpDataType, order=self.input_var.order, itype='chimera')
                for zd in dnrZones:
                    zdname = zd[0]
                    destProc = self.procDict[zdname]

                    IDs = []
                    for i in zd[2]:
                        if i[0][0:2] == 'ID':
                            if Internal.getValue(i) == zrname: IDs.append(i)

                    if IDs != []:
                        if destProc == self.rank:
                            zD = Internal.getNodeFromName2(tc, zdname)
                            zD[2] += IDs
                        else:
                            if destProc not in datas: datas[destProc] = [[zdname,IDs]]
                            else: datas[destProc].append([zdname,IDs])
                    else:
                        if destProc not in datas: datas[destProc] = []
        #Fin Interp classique


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

        #tmp = Internal.copyTree(t)   
        #t4 = X.getOversetInfo(tmp, tc, loc='center',type='extrapolated')
        #t5 = X.getOversetInfo(tmp, tc, loc='center',type='orphan')
        #C.convertPyTree2File(tc,'tc3_old.cgns')
        #C.convertPyTree2File(t4,'extrapolated.cgns')
        #C.convertPyTree2File(t5,'orphan.cgns')
        #stop


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

        if self.input_var.twoFronts:
            C._cpVars(t,'centers:cellNFront_2',tc,'cellNFront_2')
            Xmpi._setInterpTransfers(t,tc,variables=['cellNFront_2'], cellNVariable='cellNFront_2', compact=0)
        return tc


    def _specialFront2__(self, t, tc):
        test.printMem(">>> pushBackImageFront2 [start]")

        # bboxDict needed for optimised AddXZones (i.e. "layers" not None)
        # Return a dict with the zones of t as keys and their specific bboxes as key values
        bboxDict  = Cmpi.createBboxDict(t)
        interDict = X.getIntersectingDomains(self.tbbc)
        graph     = Cmpi.computeGraph(self.tbbc, type='bbox', intersectionsDict=interDict, reduction=False)

        # if subr, the tree subregions are kept during the exchange
        # if layers not None, only communicate the desired number of layers
        Cmpi._addLXZones(tc, graph, variables=['cellNIBC','cellNChim','cellNFront'], cartesian=self.input_var.cartesian,
                         interDict=interDict,bboxDict=bboxDict, layers=4, subr=False)
        Cmpi._addLXZones(t, graph, variables=['centers:cellNIBC', 'centers:cellNChim', 'centers:cellNFront'], cartesian=self.input_var.cartesian,
                         interDict=interDict, bboxDict=bboxDict, layers=4, subr=False)

        # Zones of tc are modified after addXZones, new tbbc, interDict and intersectionDict
        tbbcx             = G.BB(tc)
        interDict         = X.getIntersectingDomains(tbbcx)
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
                                fields = X.transferFields(zc, XI, YI, ZI, hook=None, variables=['cellNFront_origin','cellNIBC_origin'],
                                                          interpDataType=self.input_var.interpDataType, nature=1)
                                allInterpFields.append(fields)
                        if allInterpFields!=[]:
                            C._filterPartialFields(z, allInterpFields, indicesI, loc='centers', startFrom=0, filterName='donorVol',
                                                   verbose=False)

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

        # Smooth the front in case of a local refinement - only work in 2D
        if self.input_var.smoothing and self.dimPb == 2: X_IBM._smoothImageFront(t, tc)

        C._cpVars(t,'centers:cellNFront',tc,'cellNFront')

        Xmpi._setInterpTransfers(t,tc,variables=['cellNFront'], cellNVariable='cellNFront', compact=0)
        test.printMem(">>> pushBackImageFront2 [end]")

        return None


    def _buildFront__(self, t, tc, tb):
        test.printMem(">>> Building IBM front [start]")
        front = X_IBM.getIBMFront(tc, 'cellNFront', dim=self.dimPb, frontType=self.input_var.frontType)
        front = X_IBM.gatherFront(front)

        if self.input_var.twoFronts:
            front2 = X_IBM.getIBMFront(tc, 'cellNFront_2', dim=self.dimPb, frontType=self.input_var.frontType)
            front2 = X_IBM.gatherFront(front2)
        if self.isWireModel:
            frontWMM = X_IBM.getIBMFront(tc, 'cellNFrontFilWMM', dim=self.dimPb, frontType=self.input_var.frontType)
            frontWMM = X_IBM.gatherFront(frontWMM)

        if self.input_var.check and self.rank == 0:
            C.convertPyTree2File(front, 'front.cgns')
            if self.input_var.twoFronts: C.convertPyTree2File(front2, 'front2.cgns')
            if self.isWireModel: C.convertPyTree2File(frontWMM, 'frontWMM.cgns')

        self.zonesRIBC = []
        for zrcv in Internal.getZones(t):
            if C.getMaxValue(zrcv, 'centers:cellNIBC')==2.:
                zrcvname = zrcv[0]; self.zonesRIBC.append(zrcv)

        self.nbZonesIBC = len(self.zonesRIBC)
        if self.nbZonesIBC == 0:
            self.res = [{},{},{}]
            if self.input_var.twoFronts or self.isWireModel: self.res2 = [{},{},{}]
        else:
            if self.tbFilament:
                if self.dimPb ==2:
                    tb_filament_local = T.addkplane(self.tbFilament)
                    T._contract(tb_filament_local, (0,0,0), (1,0,0), (0,1,0), self.input_var.dz*0.5)
                    tb_local = Internal.merge([tb,tb_filament_local])
                else:
                    tb_filament_local = self.tbFilament
                    tb_local          = Internal.merge([tb,self.tbFilament])

            self.res = X_IBM.getAllIBMPoints(self.zonesRIBC, loc='centers',tb=tb_local, tfront=front, frontType=self.input_var.frontType,
                                             cellNName='cellNIBC', depth=self.DEPTH, IBCType=self.input_var.IBCType, Reynolds=self.Reynolds,
                                             yplus=self.input_var.yplus, Lref=self.input_var.Lref, isOrthoFirst=self.isOrthoProjectFirst)
            if self.input_var.twoFronts:
                self.res2 = X_IBM.getAllIBMPoints(self.zonesRIBC, loc='centers',tb=tb, tfront=front2, frontType=self.input_var.frontType,
                                                  cellNName='cellNIBC', depth=self.DEPTH, IBCType=self.input_var.IBCType, Reynolds=self.Reynolds,
                                                  yplus=self.input_var.yplus, Lref=self.input_var.Lref)
            if self.isWireModel:
                tb_filament_localWMM=Internal.copyTree(tb_filament_local)
                for z in Internal.getZones(tb_filament_local):
                    ibctype = Internal.getNodeFromName(Internal.getNodeFromName(z,'.Solver#define'),'ibctype')
                    if Internal.getValue(ibctype) != "wiremodel":
                        zlocal=Internal.getNodeFromName(tb_filament_localWMM,z[0])
                        Internal._rmNode(tb_filament_localWMM,zlocal)
                self.res2 = X_IBM.getAllIBMPoints(self.zonesRIBC, loc='centers',tb=tb_filament_localWMM, tfront=frontWMM, frontType=self.input_var.frontType,
                                                  cellNName='cellNFilWMM', depth=self.DEPTH, IBCType=self.input_var.IBCType, Reynolds=self.Reynolds,
                                                  yplus=self.input_var.yplus, Lref=self.input_var.Lref,isWireModel=self.isWireModel,
                                                  isOrthoFirst=self.isOrthoProjectFirst)

                self.restmp = X_IBM.getAllIBMPoints(self.zonesRIBC, loc='centers',tb=tb_filament_localWMM, tfront=frontWMM, frontType=self.input_var.frontType,
                                                    cellNName='cellNFilWMM', depth=self.DEPTH, IBCType=self.input_var.IBCType, Reynolds=self.Reynolds,
                                                    yplus=self.input_var.yplus, Lref=self.input_var.Lref, isOrthoFirst=self.isOrthoProjectFirst)

                for j in range(3):
                    ##delete in res
                    item_del=[]
                    for ii in self.res[j]:
                        if "140" in ii:
                            item_del.append(ii)
                    for ii in item_del:
                        del self.res[j][ii]

                    ##detele in res2
                    item_del=[]
                    for ii in self.res2[j]:
                        if "140" not in ii:
                            item_del.append(ii)

                    for ii in item_del:
                        del self.res2[j][ii]

                    ##add to res
                    for ii in self.restmp[j]:
                        self.res[j][ii] = self.restmp[j][ii]

        del tb


        # cleaning
        C._rmVars(tc,['cellNChim','cellNIBC','TurbulentDistance','cellNFront'])

        # dans t, il faut cellNChim et cellNIBCDnr pour recalculer le cellN a la fin
        varsRM = ['centers:gradxTurbulentDistance','centers:gradyTurbulentDistance','centers:gradzTurbulentDistance','centers:cellNFront','centers:cellNIBC']
        C._rmVars(t, varsRM)
        front = None
        if self.input_var.twoFronts: front2 = None

        if self.isWireModel:
            C._rmVars(tc,['cellNFrontFilWMM'])
            C._rmVars(t,['cellNFrontFilWMM'])
            frontWMM = None

        test.printMem(">>> Building IBM front [end]")
        return None


    def _ibcInterpolation__(self, t, tc):            

        test.printMem(">>> Interpolating IBM [start]")

        graph = {}; datas = {}

        # graph d'intersection des pts images de ce proc et des zones de tbbc
        zones  = Internal.getZones(self.tbbc)
        allBBs = []
        dictOfCorrectedPtsByIBCType = self.res[0]
        dictOfWallPtsByIBCType      = self.res[1]
        dictOfInterpPtsByIBCType    = self.res[2]
        interDictIBM={}

        if self.input_var.twoFronts or self.isWireModel:
            dictOfCorrectedPtsByIBCType2 = self.res2[0]
            dictOfWallPtsByIBCType2      = self.res2[1]
            dictOfInterpPtsByIBCType2    = self.res2[2]
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
                for nozr in range(self.nbZonesIBC):
                    if allCorrectedPts[nozr] != []:
                        zrname = self.zonesRIBC[nozr][0]
                        interpPtsBB = Generator.BB(allInterpPts[nozr])
                        for z in zones:
                            bba = C.getFields('GridCoordinates', z)[0]
                            if Generator.bboxIntersection(interpPtsBB, bba, isBB=True):
                                zname = z[0]
                                popp  = Cmpi.getProc(z)
                                if self.NP > 1:
                                    Distributed.updateGraph__(graph, popp, self.rank, zname)

                                if zrname not in interDictIBM: interDictIBM[zrname]=[zname]
                                else:
                                    if zname not in interDictIBM[zrname]: interDictIBM[zrname].append(zname)
            if self.input_var.twoFronts or self.isWireModel:
                for ibcTypeL in dictOfCorrectedPtsByIBCType2:
                    allCorrectedPts2 = dictOfCorrectedPtsByIBCType2[ibcTypeL]
                    allWallPts2      = dictOfWallPtsByIBCType2[ibcTypeL]
                    allInterpPts2    = dictOfInterpPtsByIBCType2[ibcTypeL]
                    for nozr in range(self.nbZonesIBC):
                        if allCorrectedPts2[nozr] != []:
                            zrname = self.zonesRIBC[nozr][0]
                            interpPtsBB2 = Generator.BB(allInterpPts2[nozr])
                            for z in zones:
                                bba = C.getFields('GridCoordinates', z)[0]
                                if Generator.bboxIntersection(interpPtsBB2,bba,isBB=True):
                                    zname = z[0]
                                    popp  = Cmpi.getProc(z)
                                    if self.NP > 1:
                                        Distributed.updateGraph__(graph, popp, self.rank, zname)
                                    if zrname not in interDictIBM2: interDictIBM2[zrname]=[zname]
                                    else:
                                        if zname not in interDictIBM2[zrname]: interDictIBM2[zrname].append(zname)
        else: graph={}

        if KCOMM is not None: allGraph = KCOMM.allgather(graph)
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
        Cmpi._addXZones(tc, graph, variables=['cellN'], cartesian=self.input_var.cartesian, subr=False)
        test.printMem(">>> Interpolating IBM [after addXZones]")

        ReferenceState = Internal.getNodeFromType2(t, 'ReferenceState_t')
        self.nbZonesIBC = len(self.zonesRIBC)

        for i in range(Cmpi.size): datas[i] = [] # force

        model = Internal.getNodeFromName(t, 'GoverningEquations')
        if model is not None: model = Internal.getValue(model)
        else:                 model = "Euler"

        if dictOfCorrectedPtsByIBCType!={}:
            for ibcTypeL in dictOfCorrectedPtsByIBCType:
                allCorrectedPts = dictOfCorrectedPtsByIBCType[ibcTypeL]
                allWallPts = dictOfWallPtsByIBCType[ibcTypeL]
                allInterpPts = dictOfInterpPtsByIBCType[ibcTypeL]
                for nozr in range(self.nbZonesIBC):
                    if allCorrectedPts[nozr] != []:
                        if ibcTypeL=="140#filament":
                            nlen  = numpy.shape(allInterpPts[nozr][1])[1]
                            save2file = numpy.zeros((nlen,2),dtype=float)
                            save2file[:,0]=allInterpPts[nozr][1][0][:]
                            save2file[:,1]=allInterpPts[nozr][1][1][:]
                            numpy.savetxt('allInterpPts.txt', save2file, delimiter=',')   # X is an array
                        zrcv = self.zonesRIBC[nozr]
                        zrname = zrcv[0]
                        dnrZones = []
                        for zdname in interDictIBM[zrname]:
                            zd = Internal.getNodeFromName2(tc, zdname)
                            if zd is None: print('!!!Zone None', zrname, zdname)
                            else: dnrZones.append(zd)

                        XOD._setIBCDataForZone__(zrcv, dnrZones, allCorrectedPts[nozr], allWallPts[nozr], allInterpPts[nozr],
                                                 nature=1, penalty=1, loc='centers', storage='inverse', dim=self.dimPb,
                                                 interpDataType=self.input_var.interpDataType, ReferenceState=ReferenceState, bcType=ibcTypeL,model=model)

                        nozr += 1
                        for zd in dnrZones:
                            zdname = zd[0]
                            destProc = self.procDict[zdname]

                            #allIDs = Internal.getNodesFromName(zd, 'IBCD*')
                            #IDs = []
                            #for zsr in allIDs:
                            #    if Internal.getValue(zsr)==zrname: IDs.append(zsr)

                            IDs = []
                            for i in zd[2]:
                                if i[0][0:4] == 'IBCD':
                                    if Internal.getValue(i)==zrname: IDs.append(i)

                            if IDs != []:
                                if destProc == self.rank:
                                    zD = Internal.getNodeFromName2(tc,zdname)
                                    zD[2] += IDs
                                else:
                                    if destProc not in datas: datas[destProc]=[[zdname,IDs]]
                                    else: datas[destProc].append([zdname,IDs])
                            else:
                                if destProc not in datas: datas[destProc] = []

        if dictOfCorrectedPtsByIBCType2 != {}:
            for ibcTypeL in dictOfCorrectedPtsByIBCType2:
                allCorrectedPts2 = dictOfCorrectedPtsByIBCType2[ibcTypeL]
                allWallPts2      = dictOfWallPtsByIBCType2[ibcTypeL]
                allInterpPts2    = dictOfInterpPtsByIBCType2[ibcTypeL]
                for nozr in range(self.nbZonesIBC):
                    if allCorrectedPts2[nozr] != []:
                        if ibcTypeL=="140#filament":
                            nlen  = numpy.shape(allInterpPts2[nozr][1])[1]
                            save2file = numpy.zeros((nlen,2),dtype=float)
                            save2file[:,0]=allInterpPts2[nozr][1][0][:]
                            save2file[:,1]=allInterpPts2[nozr][1][1][:]
                            numpy.savetxt('allInterpPts2.txt', save2file, delimiter=',')   # X is an array
                        zrcv     = self.zonesRIBC[nozr]
                        zrname   = zrcv[0]
                        dnrZones = []
                        for zdname in interDictIBM2[zrname]:
                            zd = Internal.getNodeFromName2(tc, zdname)
                            #if zd is not None: dnrZones.append(zd)
                            if zd is None: print('!!!Zone None', zrname, zdname)
                            else: dnrZones.append(zd)
                        XOD._setIBCDataForZone2__(zrcv, dnrZones, allCorrectedPts2[nozr], allWallPts2[nozr], None, allInterpPts2[nozr],
                                                  nature=1, penalty=1, loc='centers', storage='inverse', dim=self.dimPb,
                                                  interpDataType=self.input_var.interpDataType, ReferenceState=ReferenceState, bcType=ibcTypeL)

                        nozr += 1
                        for zd in dnrZones:
                            zdname = zd[0]
                            destProc = self.procDict[zdname]

                            IDs = []
                            for i in zd[2]:
                                if i[0][0:6] == '2_IBCD':
                                    if Internal.getValue(i)==zrname: IDs.append(i)

                            if IDs != []:
                                if destProc == self.rank:
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
        dictOfWallPtsByIBCType      = None
        dictOfInterpPtsByIBCType    = None
        interDictIBM = None
        if self.input_var.twoFronts or self.isWireModel:
            dictOfCorrectedPtsByIBCType2 = None
            dictOfWallPtsByIBCType2      = None
            dictOfInterpPtsByIBCType2    = None
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

        return None


    def _printCheckIBMInfo__(self,tc):        
        test.printMem(">>> Saving IBM infos [start]")
        tibm = X_IBM.extractIBMInfo(tc)

        # Avoid that two procs write the same information
        for z in Internal.getZones(tibm):
            if int(z[0][-1]) != self.rank:
                # Internal._rmNodesByName(tibm, z[0])
                z[0] = z[0]+"%{}".format(self.rank)
        Cmpi.convertPyTree2File(tibm, 'IBMInfo.cgns')

        if self.input_var.twoFronts:
            tibm2 = X_IBM.extractIBMInfo2(tc)

            # Avoid that two procs write the same information
            for z in Internal.getZones(tibm2):
                if int(z[0][-1]) != self.rank:
                    # Internal._rmNodesByName(tibm, z[0])
                    z[0] = z[0]+"%{}".format(self.rank)

            Cmpi.convertPyTree2File(tibm2, 'IBMInfo2.cgns')

        test.printMem(">>> Saving IBM infos [end]")
        del tibm
        if self.input_var.twoFronts: del tibm2
        return None


    def _recomputeDistRANS__(self,t,tb):       
        if 'outpress' in self.ibctypes or 'inj' in self.ibctypes or 'slip' in self.ibctypes or 'wallmodel' in self.ibctypes or 'overlap' in self.ibctypes:
            test.printMem(">>> wall distance for viscous wall only - RANS [start]")     
            for z in Internal.getZones(tb):
                ibc = Internal.getNodeFromName(z,'ibctype')
                if Internal.getValue(ibc)=='outpress' or Internal.getValue(ibc)=='inj' or Internal.getValue(ibc)=='slip' or Internal.getValue(ibc)=='wallmodel' or Internal.getValue(ibc)=='overlap':
                    Internal._rmNode(tb,z)

            if self.dimPb == 2:
                DTW._distance2Walls(t,self.tbsave,type='ortho', signed=0, dim=self.dimPb, loc='centers')
            else:
                DTW._distance2Walls(t,self.tbsave,type='ortho', signed=0, dim=self.dimPb, loc='centers')
            C._initVars(t, '{centers:TurbulentDistance}={centers:TurbulentDistance}*({centers:cellN}>0.)+(-1.)*{centers:TurbulentDistance}*({centers:cellN}<1.)')        
            test.printMem(">>> wall distance for viscous wall only - RANS [end]")  
        return None


    def _tInitialize__(self,t,tb=None):        
        if self.input_var.tinit is None: I._initConst(t, loc='centers')
        else:
            t = Pmpi.extractMesh(self.input_var.tinit, t, mode='accurate')
        if self.model != "Euler": C._initVars(t, 'centers:ViscosityEddy', 0.)

        if self.isWireModel:
            vars_wm = ['Density','VelocityX','VelocityY','VelocityZ','Temperature']
            if self.model == 'NSTurbulent':vars_wm.append('TurbulentSANuTilde')     
            for z in Internal.getZones(t):
                for v_local in vars_wm:
                    C._initVars(z,'{centers:'+v_local+'_WM}=0.')

        # Init with BBox
        if self.input_var.initWithBBox>0.:
            print('initialisation par bounding box')
            bodybb = C.newPyTree(['Base'])
            for base in Internal.getBases(tb):
                bbox = G.bbox(base)
                bodybbz = D.box(tuple(bbox[:3]),tuple(bbox[3:]), N=2, ntype='STRUCT')
                Internal._append(bodybb,bodybbz,'Base')
            T._scale(bodybb, factor=(self.input_var.initWithBBox,self.input_var.initWithBBox,self.input_var.initWithBBox))
            tbb = G.BB(t)
            interDict = X.getIntersectingDomains(tbb,bodybb,taabb=tbb,taabb2=bodybb)
            for zone in Internal.getZones(t):
                zname = Internal.getName(zone)
                if interDict[zname] != []:
                    C._initVars(zone, 'centers:MomentumX', 0.)
                    C._initVars(zone, 'centers:MomentumY', 0.)
                    C._initVars(zone, 'centers:MomentumZ', 0.)
        return None


    # Prepare
    def prepare(self, t_case, t_out, tc_out):        
        self.rank  = Cmpi.rank
        self.DEPTH = 2
        self.isWireModel=False
        for z in Internal.getZones(t_case):
            soldef  = Internal.getNodeFromName(z,'.Solver#define')
            if soldef is not None:
                ibctype = Internal.getNodeFromName(soldef,'ibctype')
                if ibctype is not None:
                    if Internal.getValue(ibctype) == "wiremodel":
                        self.isWireModel=True
                        break
        if self.input_var.extrusion is not None: self.input_var.cartesian = False

        if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
        else: tb = Internal.copyTree(t_case)

        if self.input_var.t_in is not None:
            refState=Internal.getNodeFromName(tb,'ReferenceState')
            flowEqn =Internal.getNodeFromName(tb,'FlowEquationSet')
            for b in Internal.getBases(self.input_var.t_in):
                Internal.addChild(b, refState, pos=0)
                Internal.addChild(b, flowEqn , pos=0)
        ## ================================================
        ## ============== Prelim. Filament ================
        ## ================================================
        ##OUT - filamentBases        :list of names of open geometries
        ##OUT - isFilamentOnly       :boolean if there is only a filament in tb
        ##OUT - isOrthoProjectFirst  :bolean to do orthonormal projection first
        ##OUT - tb(in place)         :tb of solid geometries only
        ##OUT - tbFilament           :tb of filament geometries only
        self._determineClosedSolidFilament__(tb)            

        ## ================================================
        ## ============ Get gen. info on test case ========
        ## ================================================
        # reference state
        self.refstate = C.getState(tb)
        Reynolds = Internal.getNodeFromName(tb, 'Reynolds')
        if Reynolds is not None:
            self.Reynolds = Internal.getValue(Reynolds)
            if self.Reynolds < 1.e5: self.input_var.frontType = 1
        else:
            self.Reynolds = 1.e6


        # dimension du pb
        dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
        self.dimPb = Internal.getValue(dimPb)

        model = Internal.getNodeFromName(tb, 'GoverningEquations')
        if model is None: raise ValueError('GoverningEquations is missing in input geometry tree.')
        self.model = Internal.getValue(model)    # model : Euler, NSLaminar, NSTurbulent

        # check Euler non consistant avec Musker
        if self.model == 'Euler':
            ibctype = Internal.getNodesFromName(tb, 'ibctype')
            if ibctype is not None:
                if 'Musker' in ibctype or 'Log' in ibctype:
                    raise ValueError("In tb: governing equations (Euler) not consistent with ibc type (%s)"%(ibctype))

        #Get the types of IBCs in tb
        ibctypes = set()
        for node in Internal.getNodesFromName(tb,'ibctype'):
            ibctypes.add(Internal.getValue(node))
        self.ibctypes = list(ibctypes)
        ## ================================================
        ## ===== Generate automatic Cartesian mesh ========
        ## ================================================
        ##IN  - tb                    :closed geometry tree (for filaments use tbbox)
        ##OUT - t                     :cartesian mesh
        if self.dimPb == 2 and self.input_var.cleanCellN == False: C._initVars(tb, 'CoordinateZ', 0.) # forced

        if self.input_var.t_in is None:
            t = self.generateCartesian__(Internal.merge([tb,self.tbFilament]),ext_local=self.input_var.ext+1)
        else:
            self.NP = Cmpi.size
            t = self.input_var.t_in

        # Balancing
        if self.input_var.balancing:
            test.printMem(">>> balancing [start]")
            import Distributor2.Mpi as D2mpi
            ts     = Cmpi.allgatherTree(Cmpi.convert2SkeletonTree(t))
            stats  = D2._distribute(ts, self.NP, algorithm='graph')
            D2._copyDistribution(t , ts)
            D2mpi._redispatch(t)
            del ts
            test.printMem(">>> balancing [end]")

        if self.input_var.extrusion == 'cyl':
            T._cart2Cyl(t, (0,0,0),(1,0,0))
            T._cart2Cyl(tb, (0,0,0),(1,0,0))

        ## ================================================
        ## ============== Distance to walls  ==============
        ## ================================================
        ##IN  - t                     :t tree (cartesian mesh or input mesh)
        ##IN  - tb                    :closed geometry tree (for filaments use tbbox)
        ##IN  - tbFilament            :filament geometry
        ##IN  - filamentBases         :list of filament bases
        ##OUT - cptBody               :Number of zones/bodies for F42 and multibody
        ##OUT - tbsave                :geometry (closed & filament) tree shifted in the z direction for 2D cases
        self._distance2wallCalc__(t,tb)

        if self.dimPb == 2:
            # Creation du corps 2D pour le preprocessing IBC
            T._addkplane(tb)
            T._contract(tb, (0,0,0), (1,0,0), (0,1,0), self.input_var.dz)

        X._applyBCOverlaps(t, depth=self.DEPTH, loc='centers', val=2, cellNName='cellN')
        C._initVars(t,'{centers:cellNChim}={centers:cellN}')
        C._initVars(t, 'centers:cellN', 1.)

        ## ================================================
        ## ============== Geometry blanking ===============
        ## ================================================
        ##IN  - t                     :t tree (cartesian mesh or input mesh)
        ##IN  - tb                    :closed geometry tree 
        ##OUT - t                     :t tree with blanked done
        t = self.blanking__(t,tb)

        C._initVars(t, '{centers:cellNIBC}={centers:cellN}')

        if self.input_var.IBCType == -1:
            #print('Points IBC interieurs: on repousse le front un peu plus loin.')
            C._initVars(t,'{centers:cellNDummy}=({centers:cellNIBC}>0.5)*({centers:cellNIBC}<1.5)')
            X._setHoleInterpolatedPoints(t,depth=1,dir=1,loc='centers',cellNName='cellNDummy',addGC=False)
            C._initVars(t,'{centers:cellNFront}=logical_and({centers:cellNDummy}>0.5, {centers:cellNDummy}<1.5)')
            C._rmVars(t, ['centers:cellNDummy'])
            for z in Internal.getZones(t):
                connector._updateNatureForIBM(z, self.input_var.IBCType,
                                              Internal.__GridCoordinates__,
                                              Internal.__FlowSolutionNodes__,
                                              Internal.__FlowSolutionCenters__)
        else:
            C._initVars(t,'{centers:cellNFront}=logical_and({centers:cellNIBC}>0.5, {centers:cellNIBC}<1.5)')
            if self.isWireModel:
                C._initVars(t,'{centers:cellNFrontFilWMM}={centers:cellNFilWMM}*({centers:cellNFilWMM}>0.5)+1*({centers:cellNFilWMM}<0.5)')
                C._initVars(t,'{centers:cellNFrontFilWMM}=logical_and({centers:cellNFrontFilWMM}>0.5, {centers:cellNFrontFilWMM}<1.5)')

            for z in Internal.getZones(t):
                if self.input_var.twoFronts:
                    epsilon_dist = abs(C.getValue(z,'CoordinateX',1)-C.getValue(z,'CoordinateX',0))
                    dmin = math.sqrt(3)*4*epsilon_dist
                    if self.input_var.frontType == 42:
                        SHIFTB = G_IBM_Height.computeModelisationHeight(Re=self.Reynolds, yplus=self.input_var.yplus, L=self.input_var.Lref)
                        dmin = max(dmin, SHIFTB+math.sqrt(3)*2*epsilon_dist) # where shiftb = hmod
                    C._initVars(z,'{centers:cellNIBC_2}=({centers:TurbulentDistance}>%20.16g)+(2*({centers:TurbulentDistance}<=%20.16g)*({centers:TurbulentDistance}>0))'%(dmin,dmin))
                    C._initVars(z,'{centers:cellNFront_2}=logical_and({centers:cellNIBC_2}>0.5, {centers:cellNIBC_2}<1.5)')

                connector._updateNatureForIBM(z, self.input_var.IBCType,
                                              Internal.__GridCoordinates__,
                                              Internal.__FlowSolutionNodes__,
                                              Internal.__FlowSolutionCenters__)

        ##Ghost kmin et kmax donneuse potentiel
        if self.input_var.extrusion is not None:
            listvars_local =['cellNChim','cellNIBC']
            for z in Internal.getZones(t):
                sol            = Internal.getNodeFromName(z,'FlowSolution#Centers')
                for var in listvars_local:
                    cellN          = Internal.getNodeFromName(sol,var)[1]
                    sh             = numpy.shape(cellN)
                    for k in [0,1, sh[2]-2, sh[2]-1]:
                        for j in range(sh[1]):
                            for i in range(sh[0]):
                                if  cellN[i,j,k] != 0:  cellN[i,j,k] =1

        ## ================================================
        ## ========== Interpdata & InterpTransfer =========
        ## ================================================
        ##IN  - t                     :t tree (cartesian mesh or input mesh)
        ##OUT - tc                    :front 42 if yplus is not provided or if the wallAdapt approach is used
        ##OUT - procDict              :Dictionary of procs
        ##OUT - tbbc                  :pytree of the bounding box of tc        
        tc=self.setInterpDataAndSetInterpTransfer__(t)

        ## ================================================
        ## ======= Specific treatment for front 2 =========
        ## ================================================
        if self.input_var.frontType == 2: self._specialFront2__(t, tc)

        C._rmVars(t,['centers:cellNFront'])
        if self.input_var.twoFronts:C._rmVars(t,['centers:cellNFront_2', 'centers:cellNIBC_2'])

        C._cpVars(t,'centers:TurbulentDistance',tc,'TurbulentDistance')

        print('Minimum distance: %f.'%C.getMinValue(t,'centers:TurbulentDistance'))
        P._computeGrad2(t, 'centers:TurbulentDistance', ghostCells=True, withCellN=False)

        ## ================================================
        ## ============== Building Front ==================
        ## ================================================
        ##IN  - t   :t tree (cartesian mesh or input mesh)
        ##IN  - tc  :tc tree
        ##IN  - tb  :tb tree
        ##OUT - res :IBM point for the 1st front
        ##OUT - res2:IBM point for the 2nd front
        self._buildFront__(t, tc, tb)


        ## ================================================
        ## ============== IBC Interpolation ===============
        ## ================================================
        self._ibcInterpolation__(t, tc) 

        C._initVars(t,'{centers:cellN}=minimum({centers:cellNChim}*{centers:cellNIBCDnr},2.)')
        varsRM = ['centers:cellNChim','centers:cellNIBCDnr']
        if self.model == 'Euler': varsRM += ['centers:TurbulentDistance']
        C._rmVars(t, varsRM)

        #-----------------------------------------
        # Computes distance field for Musker only
        #-----------------------------------------
        ## This was added in Revision 4265 - Comment: "Apps: IBM extrude cart et cylindrique"
        ## Need to understand why it was added before the recompute for RANS & was not included in the recompute of dist2wall for RANS.
        ## Left here for now. Not efficient but acceptable (for now).
        if self.model != 'Euler' and self.input_var.recomputeDist and (self.input_var.extrusion!='cyl' and self.input_var.extrusion !='cart'):
            if 'outpress' in self.ibctypes or 'inj' in self.ibctypes or 'slip' in self.ibctypes or 'wallmodel' in self.ibctypes or 'overlap' in self.ibctypes:
                test.printMem(">>> wall distance for viscous wall only [start]")
                for z in Internal.getZones(tb):
                    ibc = Internal.getNodeFromName(z,'ibctype')
                    if Internal.getValue(ibc)=='outpress' or Internal.getValue(ibc)=='inj' or Internal.getValue(ibc)=='slip' or Internal.getValue(ibc)=='wallmodel' or Internal.getValue(ibc)=='overlap':
                        Internal._rmNode(tb,z)
                if self.dimPb == 2:
                    DTW._distance2Walls(t,self.tbsave,type='ortho', signed=0, dim=self.dimPb, loc='centers')
                else:
                    DTW._distance2Walls(t,self.tbsave,type='ortho', signed=0, dim=self.dimPb, loc='centers')
                test.printMem(">>> wall distance for viscous wall only [end]")

                if self.dimPb == 2 and self.input_var.cleanCellN == False: C._initVars(t, '{centers:TurbulentDistanceWallBC}={centers:TurbulentDistance}')
        else:
            for z in Internal.getZones(t):
                dist = Internal.getNodeFromName2(z,'TurbulentDistanceWallBC')
                if dist is not None:  C._initVars(t, '{centers:TurbulentDistance}={centers:TurbulentDistanceWallBC}')

        ## ================================================
        ## =========== Save IBM/IBC info ==================
        ## ================================================
        if self.input_var.check: self._printCheckIBMInfo__(tc)

        if self.input_var.extrusion == 'cyl':
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
                            for l in range(numpy.size(r)):
                                yy  = r[l]*numpy.cos( theta[l] )
                                zz  = r[l]*numpy.sin( theta[l] )
                                r[l]= yy; theta[l] = zz

        ## ================================================
        ## =========== Redistribute - Final ===============
        ## ================================================
        # distribution par defaut (sur NP)
        #note:distrib does not work as tbbc does not have ID
        #     can be deleted (in the future)

        # Perform the final distribution
        if self.input_var.distrib:
            stats = D2._distribute(self.tbbc, self.NP, algorithm='graph', useCom='ID')
            D2._copyDistribution(tc, self.tbbc)
            D2._copyDistribution(t, self.tbbc)
        self.tbbc = None

        if self.input_var.redistribute:
            import Distributor2.Mpi as D2mpi
            tcs    = Cmpi.allgatherTree(Cmpi.convert2SkeletonTree(tc))
            stats  = D2._distribute(tcs, self.NP, algorithm='graph')
            D2._copyDistribution(tc, tcs)
            D2._copyDistribution(t , tcs)
            D2mpi._redispatch(tc)
            D2mpi._redispatch(t)
            self._checkNcellsNptsPerProc(tc, isAtCenter=True)

        ## ================================================
        ## === Recompute Distance for Wall (RANS Only) ====
        ## ================================================
        if self.model == 'NSTurbulent':self._recomputeDistRANS__(t, tb)

        ## ================================================
        ## ======== Remaning & Saving tc tree  ============
        ## ================================================
        if self.input_var.twoFronts or self.isWireModel:
            tc2 = Internal.copyTree(tc)
            tc2 = Internal.rmNodesByName(tc2, 'IBCD*')
            tc  = Internal.rmNodesByName(tc, '2_IBCD*')

            if self.isWireModel:
                tc2 = Internal.rmNodesByName(tc2, 'ID*')
                tc2 = Internal.rmNodesByName(tc2, 'gradxPressure')
                tc2 = Internal.rmNodesByName(tc2, 'gradyPressure')
                tc2 = Internal.rmNodesByName(tc2, 'gradzPressure')

                tc2    = transformTc2(tc2)
                NewIBCD= 141
                _changeNameIBCD__(tc2,NewIBCD)
                tc=Internal.merge([tc,tc2])

                del tc2

        if RENAMEIBCNODES:
            for zc in Internal.getZones(tc):
                for ibcd in Internal.getNodesFromName1(zc,'IBCD_*'):            
                    proposedName = Internal.getName(ibcd)[0:6]+'_X%d'%(self.rank)
                    ibcd[0] = getIBCDName(proposedName)

            if self.input_var.twoFronts:
                for zc in Internal.getZones(tc2):
                    for ibcd in Internal.getNodesFromName1(zc,'2_IBCD_*'):            
                        proposedName = Internal.getName(ibcd)[0:8]+'_X%d'%(self.rank)
                        ibcd[0] = getIBCDName(proposedName)

        ##Adding a userdefined node to the tc tree for the IBC conditions that are provided
        ##to FastS solver to reduce the number of input arguments and to make a clear distinction
        ##of the solver parameters and the IBC parameters
        base      = Internal.getBases(tc)[0]
        Internal._createUniqueChild(base, '.Solver#IBCdefine', 'UserDefinedData_t')
        solverIBC = Internal.getNodeFromName1(base, '.Solver#IBCdefine')

        Internal._createUniqueChild(solverIBC, 'Reref', 'DataArray_t', -1)
        Internal._createUniqueChild(solverIBC, 'Lref' , 'DataArray_t',  1)

        Internal._createUniqueChild(solverIBC, 'isgradP'    , 'DataArray_t', 'False')
        Internal._createUniqueChild(solverIBC, 'isWireModel', 'DataArray_t', 'False')
        Internal._createUniqueChild(solverIBC, 'isTBLE'     , 'DataArray_t', 'False')

        ##note: here alphagrad is the corrected nomenclature for alghagradp found in param_solver.h (B.Constant confirmed)
        ##      changed some other variables names to be consistent with other options/coding "guidelines"
        if 'Mafzal' in self.ibctypes:
            Internal._createUniqueChild(solverIBC, 'isgradP'   , 'DataArray_t', 'True')
            Internal._createUniqueChild(solverIBC, 'mafzalMode', 'DataArray_t', 0)
            Internal._createUniqueChild(solverIBC, 'alphaGrad' , 'DataArray_t', 0)            
        if self.input_var.twoFronts:
            Internal._createUniqueChild(solverIBC, 'isgradP'    , 'DataArray_t', 'True')
            Internal._createUniqueChild(solverIBC, 'alphaGrad'  , 'DataArray_t', 0)
        if self.isWireModel:
            Internal._createUniqueChild(solverIBC, 'isWireModel' , 'DataArray_t', 'True')
            Internal._createUniqueChild(solverIBC, 'DeltaVWire'  , 'DataArray_t', 0)
            Internal._createUniqueChild(solverIBC, 'KWire'       , 'DataArray_t', 0)
            Internal._createUniqueChild(solverIBC, 'DiameterWire', 'DataArray_t', 0)
            Internal._createUniqueChild(solverIBC, 'CtWire'      , 'DataArray_t', 0) 
        if 'TBLE' in self.ibctypes or 'TBLE_FULL' in self.ibctypes:
            Internal._createUniqueChild(solverIBC, 'isTBLE'        , 'DataArray_t', 'True')
            Internal._createUniqueChild(solverIBC, 'alphaGrad'     , 'DataArray_t', 0)
            Internal._createUniqueChild(solverIBC, 'NbPtsLinelits' , 'DataArray_t', 0)            

        if isinstance(tc_out, str):
            tcp = Compressor.compressCartesian(tc)
            Cmpi.convertPyTree2File(tcp, tc_out, ignoreProcNodes=True)

            if self.input_var.twoFronts:
                tc2  = transformTc2(tc2)
                tcp2 = Compressor.compressCartesian(tc2)
                Cmpi.convertPyTree2File(tcp2, 'tc2.cgns', ignoreProcNodes=True)
                del tc2

        ## ================================================
        ## ========= Initialization of t tree  ============
        ## ================================================
        self._tInitialize__(t)

        ## ================================================
        ## ========= t tree cleaning & saving  ============
        ## ================================================
        #Nettoyage arbre
        if self.input_var.extrusion is not None:
            vars = ['centers:TurbulentDistanceAllBC','centers:TurbulentDistanceWallBC', 'centers:cellNIBC_hole']
            C._rmVars(t, vars)

        if not self.input_var.redistribute:
            self._checkNcellsNptsPerProc(tc,isAtCenter=True)
        # Save t         
        if isinstance(t_out, str):
            tp = Compressor.compressCartesian(t)
            Cmpi.convertPyTree2File(tp, t_out, ignoreProcNodes=True)

        if Cmpi.size > 1: Cmpi.barrier()
        return t, tc


    ##"STAND ALONE" FUNCTIONS
    # distribute
    def _distribute(self, t_in, tc_in,algorithm='graph', tc2_in=None):
        return _distribute(t_in, tc_in, self.NP, algorithm='graph', tc2_in=None)


    # check nulber of points and cells per zone and i total
    def _checkNcellsNptsPerProc(self, ts, isAtCenter=False):
        return _checkNcellsNptsPerProc(ts,self.NP,isAtCenter=isAtCenter)


    # post-processing: extrait la solution aux noeuds + le champs sur les surfaces
    def post(self, t_case, t_in, tc_in, t_out, wall_out):
        return post(t_case, t_in, tc_in, t_out, wall_out)


    # post-processing: extrait les efforts sur les surfaces
    def loads(self, t_case, tc_in=None, wall_out=None, alpha=0., beta=0., Sref=None, famZones=[]):
        return loads(t_case, tc_in=tc_in, wall_out=wall_out, alpha=alpha, beta=beta, Sref=Sref, famZones=famZones)


    # post-processing: extrait les efforts sur les surfaces
    def extrudeCartesian(self, t, tb):
        return extrudeCartesian(t,tb,check=self.input_var.check, extrusion=self.input_var.extrusion,
                                dz=self.input_var.dz, NPas=self.input_var.NPas, span=self.input_var.span,
                                Ntranche=self.input_var.Ntranche,
                                dictNz=self.input_var.dictNz, ific=self.input_var.ific,
                                isCartesianExtrude=self.input_var.isCartesianExtrude)


    def setInterpData_Hybride(self,t_3d, tc_3d, t_curvi, interpDataType=1):
        return setInterpData_Hybride(t_3d, tc_3d, t_curvi, extrusion=self.input_var.extrusion, interpDataType=interpDataType)


## IMPORTANT NOTE !!
## FUNCTIONS MIGRATED TO $CASSIOPEE/Cassiopee/Post/Post/IBM.py
## The functions below will become deprecated after Jan. 1 2023
#====================================================================================
def extractIBMInfo(tc_in, t_out='IBMInfo.cgns'):
    tibm = P_IBM.extractIBMInfo(tc_in, t_out=t_out)
    return tibm

def _loads0(ts, Sref=None, Pref=None, Qref=None, alpha=0., beta=0., dimPb=3, verbose=False):
    return P_IBM.loads0(ts, Sref=Sref, Pref=Pref, Qref=Qref, alpha=alpha, beta=beta, dimPb=dimPb, verbose=verbose)

def loads0(ts, Sref=None, alpha=0., beta=0., dimPb=3, verbose=False):
    return P_IBM.loads0(ts, Sref=Sref, alpha=alpha, beta=beta, dimPb=dimPb, verbose=verbose)

def extractPressureHO(tc):
    tp = P_IBM.extractPressureHO(tc)
    return tp

def extractPressureHO2(tc):
    tp = P_IBM.extractPressureHO2(tc)
    return tp

def extractConvectiveTerms(tc):
    return P_IBM.extractConvectiveTerms(tc)

def unsteadyLoads(tb, Sref=None, Pref=None, Qref=None, alpha=0., beta=0.):
    return P_IBM.unsteadyLoads(tb, Sref=Sref, Pref=Pref, Qref=Qref, alpha=alpha, beta=beta)

def _unsteadyLoads(tb, Sref=None, Pref=None, Qref=None, alpha=0., beta=0.):
    return P_IBM._unsteadyLoads(tb, Sref=Sref, Pref=Pref, Qref=Qref, alpha=alpha, beta=beta)

def _prepareSkinReconstruction(ts, tc):
    return P_IBM._prepareSkinReconstruction(ts, tc)

def _computeSkinVariables(ts, tc, tl, graphWPOST,ProjectOldVersion=False):
    return P_IBM._computeSkinVariables(ts, tc, tl, graphWPOST, ProjectOldVersion=ProjectOldVersion)

def _modifIBCD(tc):
    raise NotImplementedError("_modifyIBCD is obsolete. Use _initOutflow and _initInj functions.")


## IMPORTANT NOTE !!
## FUNCTIONS MIGRATED TO $CASSIOPEE/Cassiopee/Generator/Generator/IBMmodelHeight.py
## The functions below will become deprecated after Jan. 1 2023
#====================================================================================
def compute_Cf(Re, Cf_law='ANSYS'):
    val=G_IBM_Height.compute_Cf(Re, Cf_law=Cf_law)
    return val

def computeYplusOpt(Re=None, tb=None, Lref=1., q=1.2, snear=None, Cf_law='ANSYS'):
    val=G_IBM_Height.computeYplusOpt(Re=Re, tb=tb, Lref=Lref, q=q, snear=snear, Cf_law=Cf_law)
    return val

def computeSnearOpt(Re=None, tb=None, Lref=1., q=1.2, yplus=300., Cf_law='ANSYS'):
    val=G_IBM_Height.computeSnearOpt(Re=Re, tb=tb, Lref=Lref, q=q, yplus=yplus, Cf_law=Cf_law)
    return val

## IMPORTANT NOTE !!
## FUNCTIONS MIGRATED TO $CASSIOPEE/Cassiopee/Post/Post/IBM.py
## The functions below will become deprecated after April. 1 2023
#================================================================================
# IBM prepare - seq
def prepare0(t_case, t_out, tc_out, snears=0.01, dfars=10.,
             tbox=None, snearsf=None, yplus=100.,
             vmin=21, check=False, format='single', frontType=1, recomputeDist=False,
             expand=3, tinit=None, initWithBBox=-1., wallAdapt=None, dfarDir=0, check_snear=False):
    prep_local=IBM()
    prep_local.input_var.snears                 =snears
    prep_local.input_var.dfars                  =dfars
    prep_local.input_var.tbox                   =tbox
    prep_local.input_var.snearsf                =snearsf
    prep_local.input_var.yplus                  =yplus
    prep_local.input_var.vmin                   =vmin
    prep_local.input_var.check                  =check
    prep_local.input_var.check_snear            =check_snear
    prep_local.input_var.format                 =format
    prep_local.input_var.frontType              =frontType
    prep_local.input_var.recomputeDist          =recomputeDist
    prep_local.input_var.expand                 =expand
    prep_local.input_var.tinit                  =tinit
    prep_local.input_var.initWithBBox           =initWithBBox
    prep_local.input_var.wallAdapt              =wallAdapt
    prep_local.input_var.dfarDir                =dfarDir

    t,tc = prep_local.prepare(t_case, t_out, tc_out)

    return t, tc

#================================================================================
# IBM prepare - parallel
def prepare1(t_case, t_out, tc_out, t_in=None, to=None, snears=0.01, dfars=10.,
             tbox=None, snearsf=None, yplus=100., Lref=1.,
             vmin=21, check=False, format='single', interpDataType=0, order=2, ext=2, nature=1,optimized=1,
             frontType=1, extrusion=None, smoothing=False, balancing=False, recomputeDist=False,
             distrib=True, expand=3, tinit=None, initWithBBox=-1., wallAdapt=None, yplusAdapt=100., dfarDir=0, 
             correctionMultiCorpsF42=False, blankingF42=False, twoFronts=False, redistribute=False, IBCType=1,
             height_in=-1.0,isFilamentOnly=False, cleanCellN=True, check_snear=False, generateCartesianMeshOnly=False, tbOneOver=None):
    prep_local=IBM()
    prep_local.input_var.t_in                   =t_in
    prep_local.input_var.to                     =to
    prep_local.input_var.snears                 =snears
    prep_local.input_var.dfars                  =dfars
    prep_local.input_var.tbox                   =tbox
    prep_local.input_var.snearsf                =snearsf
    prep_local.input_var.yplus                  =yplus
    prep_local.input_var.Lref                   =Lref
    prep_local.input_var.vmin                   =vmin
    prep_local.input_var.check                  =check
    prep_local.input_var.check_snear            =check_snear
    prep_local.input_var.format                 =format
    prep_local.input_var.interpDataType         =interpDataType
    prep_local.input_var.nature                 =nature
    prep_local.input_var.optimized              =optimized
    prep_local.input_var.order                  =order
    prep_local.input_var.ext                    =ext
    prep_local.input_var.frontType              =frontType
    prep_local.input_var.extrusion              =extrusion
    prep_local.input_var.smoothing              =smoothing
    prep_local.input_var.balancing              =balancing
    prep_local.input_var.recomputeDist          =recomputeDist
    prep_local.input_var.distrib                =distrib
    prep_local.input_var.expand                 =expand
    prep_local.input_var.tinit                  =tinit
    prep_local.input_var.initWithBBox           =initWithBBox
    prep_local.input_var.wallAdapt              =wallAdapt
    prep_local.input_var.yplusAdapt             =yplusAdapt
    prep_local.input_var.dfarDir                =dfarDir
    prep_local.input_var.correctionMultiCorpsF42=correctionMultiCorpsF42
    prep_local.input_var.blankingF42            =blankingF42
    prep_local.input_var.twoFronts              =twoFronts
    prep_local.input_var.redistribute           =redistribute
    prep_local.input_var.height_in              =height_in
    prep_local.input_var.IBCType                =IBCType
    prep_local.input_var.cleanCellN             =cleanCellN
    prep_local.input_var.generateCartesianMeshOnly  = generateCartesianMeshOnly
    prep_local.input_var.tbOneOver              = tbOneOver

    t,tc = prep_local.prepare(t_case, t_out, tc_out)
    return t, tc

