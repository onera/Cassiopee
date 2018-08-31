# Class for FastS "IBM" prepare and compute

import Fast.PyTree as Fast
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.Internal as Internal
import Connector.PyTree as X
import Connector.ToolboxIBM as TIBM
import Dist2Walls.PyTree as DTW
import Distributor2.PyTree as D2
import Initiator.PyTree as I
from Apps.Fast.Common import Common

# IN: maillage surfacique + reference State + snears

#================================================================================
# IBM prepare
# NP is the target number of processors
#================================================================================ 
def prepare(t_case, t_out, tc_out, dfar=10., vmin=21, check=False, NP=0, format='single'):
    import Converter.Mpi as Cmpi
    rank = Cmpi.rank; size = Cmpi.size
    ret = None
    # sequential prep
    if rank == 0: ret = prepare0(t_case, t_out, tc_out, dfar, vmin, check, NP, format)
    return ret

#================================================================================
# IBM prepare - seq
#================================================================================
def prepare0(t_case, t_out, tc_out, dfar=10., vmin=21, check=False, NP=0, format='single'):

    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = t_case 

    snears = 5.e-3

    #-------------------------------------------------------
    # Refinement surfaces in the fluid
    #-------------------------------------------------------
    # snearsf: list of spacing required in the refinement surfaces
    snearsf = None    
    tbox = None
    # refinementSurfFile: surface meshes describing refinement zones
    #refinementSurfFile = 'refinementBody.cgns'
    #try: tbox = C.convertFile2PyTree(refinementSurfFile)
    #except: tbox=None # no refinement surface

    #--------------------------------------------------------
    # Get Reference State and model from body pyTree
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError, 'GoverningEquations is missing in input cgns.'
    # model: Euler, NSLaminar, NSTurbulent
    model = Internal.getValue(model)

    # reference state
    refstate = C.getState(tb)
    # dimension du pb
    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)
    if dimPb == 2: C._initVars(tb, 'CoordinateZ', 0.) # forced

    #--------------------------------------------------------
    # Generates the full Cartesian mesh
    t = TIBM.generateIBMMesh(tb, vmin, snears, dfar, DEPTH=2, NP=NP,
                             tbox=tbox, snearsf=snearsf, check=check)

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
        z0 = Internal.getZones(t)
        bb = G.bbox(z0); dz = bb[5]-bb[2]
        tb2 = C.initVars(tb,'CoordinateZ',dz*0.5)
        t = DTW.distance2Walls(t,tb2,type='ortho',signed=0, dim=dimPb,loc='centers')
    else:
        t = DTW.distance2Walls(t,tb,type='ortho',signed=0, dim=dimPb,loc='centers')
    
    #----------------------------------------
    # Create IBM info
    #----------------------------------------
    t,tc = TIBM.prepareIBMData(t, tb, frontType=1)

    # arbre donneur
    D2._copyDistribution(tc,t)
    Fast.save(tc, tc_out, split=format, NP=-NP)

    #----------------------------------------
    # Extraction des coordonnees des pts IBM
    #----------------------------------------
    if check:
        tibm = TIBM.extractIBMInfo(tc)
        C.convertPyTree2File(tibm, 'IBMInfo.cgns')
        del tibm

    # arbre de calcul
    del tc
    I._initConst(t, loc='centers')
    Fast.save(t, t_out, split=format, NP=-NP)
    return t

#=============================================================================
# Post
#==============================================================================
def post(t_case, t_in, tc_in, t_out, wall_out, NP=0, format='single'):
    import Post.PyTree as P
    from math import *

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
    zw = TIBM.extractIBMWallFields(tc, tb=tb)

    RoUInf2I = 1./(RouInf*RouInf+RovInf*RovInf+RowInf*RowInf)
    C._initVars(zw,'{Cp}=2*%f*({Pressure}-%f)*%f'%(RoInf,PInf,RoUInf2I))
    if model != 'Euler':C._initVars(zw,'{Cf}=2*{Density}*{utau}**2*%f'%RoUInf2I)

    Internal._rmNodesByName(zw, '.Solver#Param')
    Internal._rmNodesByName(zw, '.Solver#ownData')

    if isinstance(wall_out, str): C.convertPyTree2File(zw, wall_out)

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
    if isinstance(t_out, str): C.convertPyTree2File(t, t_out)

    return t, zw

#====================================================================================
class IBM(Common):
    """Preparation et caculs avec le module FastS."""
    def __init__(self, NP=None, format=None, numb=None, numz=None):
        Common.__init__(self, NP, format, numb, numz)
        self.__version__ = "0.0"
        self.authors = ["ash@onera.fr"]
        
    # Prepare : n'utilise qu'un proc pour l'instant
    def prepare(self, t_case, t_out, tc_out, dfar=10., vmin=21, check=False):
        NP = self.data['NP']
        if NP == 0: print 'Preparing for a sequential computation.'
        else: print 'Preparing for a computation on %d processors.'%NP
        ret = prepare(t_case, t_out, tc_out, dfar, vmin, check, NP, self.data['format'])
        return ret

    # post-processing: extrait la solution aux noeuds + le champs sur les surfaces
    def post(self, t_case, t_in, tc_in, t_out, wall_out):
        return post(t_case, t_in, tc_in, t_out, wall_out, self.data['NP'], self.data['format'])
