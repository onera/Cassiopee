# Class for FastS "IBM" prepare and compute

# IN: maillage surfacique + reference State + snears

#================================================================================
# Redistribue un fichier in place sans com pour l'instant
# Change les noeuds procs seulement
#================================================================================
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
# IBM prepare
# NP is the target number of processors
#================================================================================ 
def prepare(t_case, t_out, tc_out, NP=0, format='single'):
    import Converter.Mpi as Cmpi
    rank = Cmpi.rank; size = Cmpi.size
    ret = None
    # sequential prep
    if rank == 0: ret = prepare0(t_case, t_out, tc_out, NP, format)
    # parallel prep
    #prepare1(t_case, t_out, tc_out, NP, format)
    return ret

#================================================================================
# Multibloc prepare - seq
#================================================================================
def prepare0(t_case, t_out, tc_out, vmin, dfar, NP=0, format='single'):

    if isinstance(t_case, str): t = C.convertFile2PyTree(t_case)
    else: t = t_case 

    #-------------------------------------------------------
    # Refinement surfaces in the fluid
    #-------------------------------------------------------
    # snearsf: list of spacing required in the surfaces
    # refinementSurfFile: surface meshes describing refinement zones
    refinementSurfFile = 'refinementBody.cgns' 
    try: tbox = C.convertFile2PyTree(refinementSurfFile)
    except: tbox=None # no refinement surface

    # if check=True: extract some files for checks
    check = False

    #--------------------------------------------------------
    # Get Reference State and model from body pyTree
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError, 'GoverningEquations is missing in input cgns.'
    # model : Euler, NSLaminar, NSTurbulent
    model = Internal.getValue(model)

    # reference state
    refstate = C.getState(tb)
    # dimension du pb
    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)
    if dimPb == 2: C._initVars(tb, 'CoordinateZ', 0.) # forced

    #--------------------------------------------------------
    # Generates the full Cartesian mesh
    t = IBM.generateIBMMesh(tb,vmin,snears,dfar,DEPTH=2,NP=NP,tbox=tbox,snearsf=snearsf,check=check)

    #------------------------------------------------------
    # distribute the mesh over NP processors
    if NP > 0:
        print 'distribution over %d processors'%NP
        t,stats = D2.distribute(t, NP)
        if check == True: print stats

    #------------------------------------------------
    # Add reference state to the pyTree and init fields
    # Add viscosity if model is not Euler
    if model != "Euler": C._initVars(t, 'centers:ViscosityEddy', 0.)
    C._addState(t, state=refstate)
    C._addState(t, 'GoverningEquations', model)
    C._addState(t, 'EquationDimension', dimPb)
    if check: C.convertPyTree2File(t, 'mesh1.cgns')

    #----------------------------------------
    # Computes distance field
    #----------------------------------------
    if dimPb == 2:
        z0 = Internal.getNodeFromType2(t,"Zone_t")
        bb = G.bbox(z0); dz = bb[5]-bb[2]
        tb2 = C.initVars(tb,'CoordinateZ',dz*0.5)
        t = DTW.distance2Walls(t,tb2,type='ortho',signed=0, dim=dimPb,loc='centers')
    else:
        t = DTW.distance2Walls(t,tb,type='ortho',signed=0, dim=dimPb,loc='centers')
    
    #----------------------------------------
    # Create IBM info
    #----------------------------------------
    t,tc = IBM.prepareIBMData(t, tb, frontType=1)

    # arbre donneur
    D2._copyDistribution(tc,t)
    Fast.save(tc, 'tc.cgns', split='single', NP=-NP)

    #----------------------------------------
    # Extraction des coordonnees des pts IBM
    #----------------------------------------
    if check:
        tibm = IBM.extractIBMInfo(tc)
        C.convertPyTree2File(tibm, 'IBMInfo.cgns')
        del tibm

    # arbre de calcul
    del tc
    I._initConst(t, loc='centers')
    Fast.save(t, 't.cgns', split='single', NP=-NP)

def post():
    import Converter.PyTree as C
    import Post.PyTree as P
    import Transform.PyTree as T
    import Connector.PyTree as X
    import Connector.ToolboxIBM as IBM
    import Converter.Internal as Internal
    import os
    from math import *

    t = C.convertFile2PyTree('restart.cgns')
    tb = C.convertFile2PyTree("case.cgns")

    #=============================
    # Supprime les champs inutiles
    #=============================
    vars = ['centers:Density_M1', 'centers:VelocityX_M1', 'centers:VelocityY_M1', 'centers:VelocityZ_M1', 'centers:Temperature_M1', 'centers:Density_P1', 'centers:VelocityX_P1', 'centers:VelocityY_P1', 'centers:VelocityZ_P1', 'centers:Temperature_P1','centers:TurbulentDistance']
    C._rmVars(t, vars)

    #=============================
    # Arbre de connectivite
    #=============================
    tc = C.convertFile2PyTree('tc.cgns')
    Internal._rmNodesByName(tc, 'GridCoordinates')

    #==========================================================
    # Extraction Cp, Cf, ... sur les surfaces par interpolation
    #==========================================================
    tb = C.convertArray2Tetra(tb)

    #--------------------------------------------------------
    # Get Reference State and model from body pyTree
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError, 'GoverningEquations is missing in input cgns.'
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

    for z in Internal.getNodesFromType2(t,"Zone_t"):
        zc = Internal.getNodeFromName(tc,z[0])
        for v in varsIBC: C._cpVars(z, 'centers:'+v, zc, v)

    X._setInterpTransfers(t, tc, variables=vars,
                          variablesIBC=varsIBC, bcType=bcType,
                          varType=varType, storage=1,
                          Gamma=Gamma, Cv=cvInf, MuS=Mus,
                          Cs=Cs, Ts=Ts)
    zw = IBM.extractIBMWallFields(tc, tb=tb)

    RoUInf2I = 1./(RouInf*RouInf+RovInf*RovInf+RowInf*RowInf)
    C._initVars(zw,'{Cp}=2*%f*({Pressure}-%f)*%f'%(RoInf,PInf,RoUInf2I))
    if model != 'Euler':C._initVars(zw,'{Cf}=2*{Density}*{utau}**2*%f'%RoUInf2I)

    Internal._rmNodesByName(zw, '.Solver#Param')
    Internal._rmNodesByName(zw, '.Solver#ownData')

    C.convertPyTree2File(zw, 'wall.cgns')

    #===============================
    # En 2D, extrait un seul plan k
    #================================
    if dimPb == 2:
        t = T.subzone(t, (1,1,1), (-1,-1,1))
        C._initVars(tb, 'CoordinateZ', 0.) # forced

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
    vars = ['centers:Density','centers:VelocityX', 'centers:Temperature','centers:ViscosityEddy', 
    'centers:TurbulentSANuTilde','centers:ViscosityMolecular', 'centers:mutsmu', 'centers:cellN']
    for v in vars: t = C.center2Node(t, v)
    Internal._rmNodesByName(t, 'FlowSolution#Centers')
    C.convertPyTree2File(t, 'out.cgns')
