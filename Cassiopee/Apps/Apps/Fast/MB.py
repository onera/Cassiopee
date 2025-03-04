# Class for FastS "Multiblock" prepare and compute

# IN: maillage volumique + BCs + raccords + reference State
# optionel: solution initiale
# Si 2D, maillage 1 seul plan (XY)

import FastC.PyTree as FastC
import Converter.PyTree as C
import Transform.PyTree as T
import Converter.Internal as Internal
import Connector.PyTree as X

#================================================================================
# Multibloc prepare (avec split)
# NP is the target number of processors
#================================================================================
def prepare(t_case, t_out, tc_out, NP=0, format='single', removeGC=True):
    import Converter.Mpi as Cmpi
    rank = Cmpi.rank; size = Cmpi.size
    ret = None
    # sequential prep
    if rank == 0: ret = prepare0(t_case, t_out, tc_out, NP=NP, format=format, removeGC=removeGC)
    return ret

#================================================================================
# Multibloc prepare - seq
#================================================================================
def prepare0(t_case, t_out, tc_out, NP=0, format='single', removeGC=True):

    if isinstance(t_case, str): t = C.convertFile2PyTree(t_case)
    else: t = t_case

    if NP > 0: import Distributor2.PyTree as D2

    # Get dim
    dim = 3
    node = Internal.getNodeFromName(t, 'EquationDimension')
    if node is not None: dim = Internal.getValue(node)

    # Split size
    if NP > 0:
        perioInfo = Internal.getPeriodicInfo(Internal.getZones(t)[0])
        t = T.splitSize(t, R=NP, type=2, minPtsPerDir=9)
        t = X.connectMatch(t, dim=dim)
        for p in perioInfo:
            if p[0] != []:
                [xc,yc,zc,vx,vy,vz,angle] = p[0]
                t = X.connectMatchPeriodic(t, rotationCenter=[xc,yc,zc], rotationAngle=[vx*angle,vy*angle,vz*angle], tol=1.e-6, dim=dim)
            if p[1] != []:
                [tx,ty,tz] = p[1]
                t = X.connectMatchPeriodic(t, translation=[tx,ty,tz], tol=1.e-6, dim=dim)
        stats = D2._distribute(t, NP, useCom='match')
    else:
        Internal._rmNodesByName(t, 'proc')

    # Solution initiale
    eqs = Internal.getNodeFromType(t, 'GoverningEquations_t')
    Model = 'NSTurbulent'
    if eqs is not None: Model = Internal.getValue(eqs)

    ret = C.isNamePresent(t, 'centers:Density')
    if ret != 1: # Density not present, init from ref state
        state = Internal.getNodeFromType(t, 'ReferenceState_t')
        if state is None:
            raise ValueError('Reference state is missing in input cgns.')
        vars = ['Density', 'MomentumX', 'MomentumY', 'MomentumZ',
                'EnergyStagnationDensity']
        for v in vars:
            node = Internal.getNodeFromName(state, v)
            if node is not None:
                val = float(node[1][0])
                C._initVars(t, 'centers:'+v, val)
            else:
                raise ValueError(v + ' is missing in ReferenceState.')

        if Model == 'NSTurbulent':
            vars = ['TurbulentSANuTildeDensity']
            for v in vars:
                node = Internal.getNodeFromName(state, v)
                if node is not None:
                    val = float(node[1][0])
                    C._initVars(t, 'centers:'+v, val)

    ret = C.isNamePresent(t, 'centers:TurbulentDistance')
    if Model == 'NSTurbulent' and ret != 1: # Wall distance not present
        import Dist2Walls.PyTree as DTW
        walls = C.extractBCOfType(t, 'BCWall')
        if walls != []:
            DTW._distance2Walls(t, walls, loc='centers', type='ortho')
        else: C._initVars(t, 'centers:TurbulentDistance', 1000000.)

    # Ajout des ghost-cells
    C.addState2Node__(t, 'EquationDimension', dim)
    Internal._addGhostCells(t, t, 2, adaptBCs=1, fillCorner=0)
    if dim == 2:
        T._addkplane(t)
        T._contract(t, (0,0,0), (1,0,0), (0,1,0), 0.01)
        T._makeDirect(t)

    # Raccords
    tc = C.node2Center(t)
    tc = X.setInterpData(t, tc, nature=1, loc='centers', storage='inverse',
                         sameName=1, dim=dim)
    C._rmVars(tc, 'FlowSolution')
    if removeGC: C._rmVars(tc, 'GridCoordinates')

    if isinstance(tc_out, str):
        #Fast.save(tc, tc_out, split=format, NP=-NP)
        C.convertPyTree2File(tc, tc_out) # pour eviter l'ecriture parallele
    if isinstance(t_out, str):
        #Fast.save(t, t_out, split=format, NP=-NP)
        C.convertPyTree2File(t, t_out) # pour eviter l'ecriture parallele
    return t, tc

#================================================================================
# Multibloc prepare - parallel - en cours
# NP is the target number of processors
#================================================================================
def prepare1(t_case, t_out, tc_out, NP=0, format='single', removeGC=True):

    import Distributor2.PyTree as D2
    import Converter.Mpi as Cmpi
    import Converter.Filter as Filter
    import Connector.PyTree as X

    h = Filter.Handle(t_case)
    t = h.loadAndSplit(NParts=max(NP, Cmpi.size))

    # Get dim
    dim = 3
    node = Internal.getNodeFromName(t, 'EquationDimension')
    if node is not None: dim = Internal.getValue(node)

    # Ajout des ghost-cells
    C.addState2Node__(t, 'EquationDimension', dim)
    # Je suis oblige d'enlever les raccords de t et de les reconstruire a cause que les raccords ne sont
    # pas dans les BXZones
    C._deleteGridConnectivity__(t)
    Cmpi._addBXZones(t, depth=3)
    t = X.connectMatch(t, dim=dim)
    Internal._addGhostCells(t, t, 2, adaptBCs=1, fillCorner=0)
    Cmpi._rmBXZones(t)
    #mergeWindows(a)

    # tc local
    tc = C.node2Center(t)

    # Arbre des bbox sur t
    #tbb = Cmpi.createBBoxTree(t)
    #interDict = X.getIntersectingDomains(tbb)
    #graph = Cmpi.computeGraph(tbb, type='bbox', intersectionsDict=interDict, reduction=False)
    #del tbb
    graph = Cmpi.computeGraph(t, type='match', reduction=True)

    Cmpi._addXZones(t, graph, variables=[], noCoordinates=True)
    Cmpi._addXZones(tc, graph, variables=[], noCoordinates=True)
    X._setInterpData(t, tc, nature=1, loc='centers', storage='inverse',
                     sameName=1, dim=dim, itype='abutting')
    Cmpi._rmXZones(t)
    Cmpi._rmXZones(tc)

    if dim == 2:
        T._addkplane(t)
        T._contract(t, (0,0,0), (1,0,0), (0,1,0), 0.01)
        T._makeDirect(t)
    Internal._rmNodesByName(tc, Internal.__FlowSolutionNodes__)
    if removeGC: Internal._rmNodesByName(tc, Internal.__GridCoordinates__)

    # Solution initiale
    eqs = Internal.getNodeFromType(t, 'GoverningEquations_t')
    Model = 'NSTurbulent'
    if eqs is not None: Model = Internal.getValue(eqs)

    ret = C.isNamePresent(t, 'centers:Density')
    if ret != 1: # Density not present, init from ref state
        state = Internal.getNodeFromType(t, 'ReferenceState_t')
        if state is None:
            raise ValueError('Reference state is missing in input cgns.')
        vars = ['Density', 'MomentumX', 'MomentumY', 'MomentumZ', 'EnergyStagnationDensity']
        for v in vars:
            node = Internal.getNodeFromName(state, v)
            if node is not None:
                val = float(node[1][0])
                C._initVars(t, 'centers:'+v, val)
            else:
                raise ValueError(v + ' is missing in ReferenceState.')

        if Model == 'NSTurbulent':
            vars = ['TurbulentSANuTildeDensity']
            for v in vars:
                node = Internal.getNodeFromName(state, v)
                if node is not None:
                    val = float(node[1][0])
                    C._initVars(t, 'centers:'+v, val)

    ret = C.isNamePresent(t, 'centers:TurbulentDistance')
    if Model == 'NSTurbulent' and ret != 1: # Wall distance not present
        import Dist2Walls.PyTree as DTW
        walls = C.extractBCOfType(t, 'BCWall')
        ret = Cmpi.allgather(walls)
        walls = []
        for i in ret: walls += i
        if walls != []:
            DTW._distance2Walls(t, walls, loc='centers', type='ortho')
        else: C._initVars(t, 'centers:TurbulentDistance', 1000000.)

    Cmpi.convertPyTree2File(t, t_out)
    Cmpi.convertPyTree2File(tc, tc_out)

    return t, tc

#==========================================================================
# Post
#===========================================================================
def post(t_in, t_out, wall_out, format='single'):
    import Converter.Mpi as Cmpi
    rank = Cmpi.rank; size = Cmpi.size
    # sequential post
    ret = None
    if rank == 0: ret = post0(t_in, t_out, wall_out, format)
    return ret

#====================================================================================
def post0(t_in, t_out, wall_out, format='single'):
    import Transform.PyTree as T
    from math import sqrt

    # Use filter load here!
    if isinstance(t_in, str): a = FastC.loadFile(t_in)
    else: a = t_in

    #=============================
    # Supprime les champs inutiles
    #=============================
    vars = ['centers:Density_M1', 'centers:VelocityX_M1', 'centers:VelocityY_M1', 'centers:VelocityZ_M1', 'centers:Temperature_M1', 'centers:Density_P1', 'centers:VelocityX_P1', 'centers:VelocityY_P1', 'centers:VelocityZ_P1', 'centers:Temperature_P1','centers:TurbulentSANuTilde']
    C._rmVars(a, vars)

    dimPb = Internal.getNodeFromName(a, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)
    if dimPb == 2: a = T.subzone(a,(1,1,1),(-1,-1,1))

    #=================================
    # Model
    #=================================
    model = Internal.getNodeFromName(a, 'GoverningEquations')
    model = Internal.getValue(model)
    if model is None: raise ValueError('GoverningEquations is missing in input file.')

    #==============================
    # Sortie champs aux noeuds
    #==============================
    vars = ['centers:Density', 'centers:VelocityX', 'centers:VelocityY', 'centers:VelocityZ', 'centers:Temperature']
    if model != 'Euler': vars += ['centers:ViscosityEddy', 'centers:TurbulentSANuTilde']
    for v in vars: a = C.center2Node(a, v)
    Internal._rmNodesByName(a,'FlowSolution#Centers')

    # Supprime les ghost cells
    Internal._rmGhostCells(a, a, 2, adaptBCs=1)
    if isinstance(t_out, str): C.convertPyTree2File(a, t_out)

    #=================================
    # Extraction Kp sur les surfaces
    #=================================
    [RoInf, RouInf, RovInf, RowInf, RoeInf, PInf, TInf, cvInf, MInf,
     ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
     Mus, Cs, Ts, Pr] = C.getState(a)
    gam1cv = (Gamma-1.)*cvInf
    RoUInf2I = 1./(RouInf*RouInf+RovInf*RovInf+RowInf*RowInf)
    C._initVars(a,'{Pressure}={Density}*{Temperature}*%f'%gam1cv)
    C._initVars(a,'{Kp}=-2*%f*({Pressure}-%f)*%f'%(RoInf,PInf,RoUInf2I))

    if model != 'Euler':
        betas = Mus*(Ts+Cs)/(Ts**(3./2.))
        C._initVars(a,'{ViscosityMolecular} = %20.16g*sqrt({Temperature})/(1.+%20.16g/{Temperature})'%(betas,Cs))
        C._initVars(a,'{mutsmu}=({ViscosityEddy})/({ViscosityMolecular})-1.')

    w = C.extractBCOfType(a,'BCWall')
    Internal._rmNodesByName(w, '.Solver#Param')
    Internal._rmNodesByName(w, '.Solver#ownData')
    if isinstance(wall_out, str): C.convertPyTree2File(w, wall_out)
    return a, w

#====================================================================================
def calc_cp_cf(t,h,listzones,isRANS=False,wallType='BCWall',mode=None):
    h._loadZones(t, znp=listzones)

    ##Removing the zones of no interest
    bases = Internal.getBases(t)
    for b in bases:
        zones = Internal.getZones(b)
        for z in zones:
            if b[0]+'/'+z[0] not in listzones:
                Internal._rmNode(t,z)

    ##Getting and calculating ref. values
    [RoInf, RouInf, RovInf, RowInf, RoeInf, PInf, TInf, cvInf, MInf,
     ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
     Mus, Cs, Ts, Pr] = C.getState(t)

    RoUInf2I = 1./(RouInf*RouInf+RovInf*RovInf+RowInf*RowInf)

    w=FastC.get_wall_values(t,isRANS=isRANS,wallType=wallType,mode=mode)
    w=FastC.get_skin_friction(w,RoInf,PInf,RoUInf2I)

    return w
