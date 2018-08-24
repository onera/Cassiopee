# Class for FastS "Multiblock" prepare and compute

# IN: maillage volumique + BCs + raccords + reference State
# optionel: solution initiale 
# Si 2D, maillage 1 seul plan (XY)

import Fast.PyTree as Fast
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.Internal as Internal
import Connector.PyTree as X
from Apps.App import App

# Redistribue un fichier in place sans com pour l'instant
def distribute(t_in, NP):
    import Filter
    t  = Filter.convertFile2SkeletonTree(t_in, maxDepth=2,maxSize=6)
    t, stats = D2.distribute(t, NP, algorithm='graph', useCom=None)
    nodes = Internal.getNodesFromName(t, 'Proc')
    for n in nodes:
        p = Internal.getPath(n)
        Filter.writeNodesFromPaths(t_in, p, n)
    return t

#================================================================================
# Multibloc prepare
# NP is the target number of processors
#================================================================================ 
def prepare(t_case, t_out, tc_out, NP=0, format='single'):
    import Converter.Mpi as Cmpi
    rank = Cmpi.rank; size = Cmpi.size
    ret = None
    if rank == 0: ret = prepare0(t_case, t_out, tc_out, NP, format)
    #prepare1(t_case, t_out, tc_out, NP, format)
    return ret

#================================================================================
# Multibloc prepare - seq
# NP is the target number of processors
#================================================================================
def prepare0(t_case, t_out, tc_out, NP=0, format='single'):

    if isinstance(t_case, str): t = C.convertFile2PyTree(t_case)
    else: t = t_case 
        
    if NP > 0: import Distributor2.PyTree as D2

    # Get dim
    dim = 3
    node = Internal.getNodeFromName(t, 'EquationDimension')
    if node is not None: dim = Internal.getValue(node)

    # Split size
    if NP > 0:
        t = T.splitSize(t, R=NP, type=2, minPtsPerDir=9)
        t = X.connectMatch(t, dim=dim)
        t, stats = D2.distribute(t, NP, useCom='match')

    # Solution initiale
    eqs = Internal.getNodeFromType(t, 'GoverningEquations_t')
    Model = 'NSTurbulent'
    if eqs is not None: Model = Internal.getValue(eqs)

    ret = C.isNamePresent(t, 'centers:Density')
    if ret != 1: # Density not present, init from ref state
        state = Internal.getNodeFromType(t, 'ReferenceState_t')
        if state is None:
            raise ValueError, 'Reference state is missing in input cgns.'
        vars = ['Density', 'MomentumX', 'MomentumY', 'MomentumZ',
        'EnergyStagnationDensity']
        for v in vars:
            node = Internal.getNodeFromName(state, v)
            if node is not None:
                val = float(node[1][0])
                C._initVars(t, 'centers:'+v, val)
            else:
                raise ValueError, v + ' is missing in ReferenceState.'

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
            t = DTW.distance2Walls(t, walls, loc='centers', type='ortho')
        else: C._initVars(t, 'centers:TurbulentDistance', 1000000.)

    # Ajout des ghost-cells
    C.addState2Node__(t, 'EquationDimension', dim)
    t = Internal.addGhostCells(t, t, 2, adaptBCs=1, fillCorner=0)
    if dim == 2: 
        t = T.addkplane(t)
        t = T.contract(t, (0,0,0), (1,0,0), (0,1,0), 0.01)
        t = T.makeDirect(t)

    # Raccords
    tc = C.node2Center(t)
    tc = X.setInterpData(t, tc, nature=1, loc='centers', storage='inverse', 
                         sameName=1, dim=dim)
    C._rmVars(tc, 'FlowSolution')
    
    if isinstance(tc_out, str):
        Fast.save(tc, tc_out, split=format, NP=-NP)
    if isinstance(t_out, str):
        Fast.save(t, t_out, split=format, NP=-NP)
    return t, tc

#================================================================================
# Multibloc prepare - parallel
# NP is the target number of processors
#================================================================================
def prepare1(t_case, t_out, tc_out, NP=0, format='single'):
    import Distributor2.PyTree as D2
    import Converter.Mpi as Cmpi

    t = Cmpi.convertFile2SkeletonTree(t_case)
    allCells = C.getNCells(t)
    graph = Cmpi.computeGraph(t, type='match')
    #t, stats = D2.distribute(t, Cmpi.size, useCom='match')
    t, stats = D2.distribute(t, Cmpi.size, useCom=None, algorithm='graph')
    t = Cmpi.readZones(t, t_case, rank=Cmpi.rank)
    Cmpi._convert2PartialTree(t)
    #if Cmpi.rank == 0: Internal.printTree(t)

    # Get dim
    dim = 3
    node = Internal.getNodeFromName(t, 'EquationDimension')
    if node is not None: dim = Internal.getValue(node)

    # Charge du proc
    ncells = C.getNCells(t)
    print Cmpi.rank, ncells, allCells
    #Cmpi.convertPyTree2File(t, 'test.cgns')
    return

    # Split size
    if NP > 0:
        t = T.splitSize(t, R=NP, type=2, minPtsPerDir=9)
        t = X.connectMatch(t, dim=dim)
        t, stats = D2.distribute(t, NP, useCom='match')
        # redispatch

    # Solution initiale
    eqs = Internal.getNodeFromType(t, 'GoverningEquations_t')
    Model = 'NSTurbulent'
    if eqs is not None: Model = Internal.getValue(eqs)

    ret = C.isNamePresent(t, 'centers:Density')
    if ret != 1: # Density not present, init from ref state
        state = Internal.getNodeFromType(t, 'ReferenceState_t')
        if state is None:
            raise ValueError, 'Reference state is missing in input cgns.'
        vars = ['Density', 'MomentumX', 'MomentumY', 'MomentumZ',
        'EnergyStagnationDensity']
        for v in vars:
            node = Internal.getNodeFromName(state, v)
            if node is not None:
                val = float(node[1][0])
                C._initVars(t, 'centers:'+v, val)
            else:
                raise ValueError, v + ' is missing in ReferenceState.'

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
            t = DTW.distance2Walls(t, walls, loc='centers', type='ortho')
        else: C._initVars(t, 'centers:TurbulentDistance', 1000000.)

    # Ajout des ghost-cells
    C.addState2Node__(t, 'EquationDimension', dim)
    t = Internal.addGhostCells(t, t, 2, adaptBCs=1, fillCorner=0)
    if dim == 2: 
        t = T.addkplane(t)
        t = T.contract(t, (0,0,0), (1,0,0), (0,1,0), 0.01)
        t = T.makeDirect(t)

    # Raccords
    tc = C.node2Center(t)
    tc = X.setInterpData(t, tc, nature=1, loc='centers', storage='inverse', 
                         sameName=1, dim=dim)
    C._rmVars(tc, 'FlowSolution')
    
    if isinstance(tc_out, str):
        Fast.save(tc, tc_out, split=format, NP=-NP)
    if isinstance(t_out, str):
        Fast.save(t, t_out, split=format, NP=-NP)
    return t, tc

#=====================================================================================
# NP is the currently running number of processors
# IN: file names
#======================================================================================
def compute(t_in, tc_in, t_out,
            numb, numz,
            NIT, 
            NP=0, format='single'):

    if NP > 0:
        import Converter.Mpi as Cmpi
        import FastS.Mpi as FastS
        rank = Cmpi.rank; size = Cmpi.size
    else:
        import FastS.PyTree as FastS
        rank = 0; size = 1

    if NP != 0 and NP != size:
        print 'Warning: you are running not on the prepared number of processors %d != %d', NP, size

    t,tc,ts,graph = Fast.load(t_in, tc_in, split=format, NP=NP)

    it0 = 0; time0 = 0.
    first = Internal.getNodeFromName1(t, 'Iteration')
    if first is not None: it0 = Internal.getValue(first)
    first = Internal.getNodeFromName1(t, 'Time')
    if first is not None: time0 = Internal.getValue(first)

    # Numerics
    Fast._setNum2Zones(t, numz); Fast._setNum2Base(t, numb)

    (t, tc, metrics) = FastS.warmup(t, tc, graph)

    time_step = Internal.getNodeFromName(t, 'time_step')
    time_step = Internal.getValue(time_step)
    for it in xrange(NIT):
        FastS._compute(t, metrics, it, tc, graph)
        if it%100 == 0:
            if rank == 0: print '- %d / %d - %f'%(it+it0, NIT+it0, time0)
        #FastS.display_temporal_criteria(t, metrics, it, format='double')
        time0 += time_step

    # time stamp
    Internal.createUniqueChild(t, 'Iteration', 'DataArray_t', value=it0+NIT)
    Internal.createUniqueChild(t, 'Time', 'DataArray_t', value=time0)
    Fast.save(t, t_out, split=format, NP=NP)
    return t

#====================================================================================
def Post(t_in, t_out, surf_out):
    import Transform.PyTree as T
    import Converter.Internal as Internal
    import Post.PyTree as P
    import os
    from math import sqrt

    # Use filter load here!
    a = C.convertFile2PyTree(t_in)

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
    if model is None: raise ValueError, 'GoverningEquations is missing in input file.'

    #==============================
    # Sortie champs aux noeuds
    #==============================
    vars = ['centers:Density', 'centers:VelocityX', 'centers:VelocityY', 'centers:VelocityZ', 'centers:Temperature']
    if model != 'Euler': vars += ['centers:ViscosityEddy', 'centers:TurbulentSANuTilde']
    for v in vars: a = C.center2Node(a, v)
    Internal._rmNodesByName(a,'FlowSolution#Centers')

    # Supprime les ghost cells
    a = Internal.rmGhostCells(a, a, 2, adaptBCs=1)
    C.convertPyTree2File(a, t_out)

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
    C.convertPyTree2File(w, surf_out)

#====================================================================================
class MB(App):
    """Prepation et caculs avec le module FastS."""
    def __init__(self):
        App.__init__(self)
        self.__version__ = "0.0"
        self.authors = ["ash@onera.fr"]
        self.requires(['NP', 'format', 'numb', 'numz'])
        # default values
        self.set(NP=0)
        self.set(format='single')
        
    # Prepare n'utilise qu'un proc pour l'instant
    def prepare(self, t_case, t_out, tc_out):
        NP = self.data['NP']
        if NP == 0: print 'Preparing for a sequential computation.'
        else: print 'Preparing for a computation on %d processors.'%NP
        ret = prepare(t_case, t_out, tc_out, NP, self.data['format'])
        return ret

    # Compute peut etre lance en sequentiel ou en parallele
    def compute(self, t_in, tc_in, t_out, nit):
        numb = self.data['numb']
        numz = self.data['numz']
        compute(t_in, tc_in, t_out,
                numb, numz,
                nit, 
                NP=self.data['NP'], 
                format=self.data['format'])
