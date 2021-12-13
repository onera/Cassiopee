"""Compute other complex variables than computeVariables."""

from . import PyTree as P
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import Transform.PyTree as T
import Geom.PyTree as D

#==============================================================================
# Vorticite en centres
#==============================================================================
def computeVorticity(t):
    t2 = Internal.copyRef(t)
    presvx = C.isNamePresent(t2, 'VelocityX')
    if presvx == -1: 
        presvxc = C.isNamePresent(t2, 'centers:VelocityX') 
        if presvxc > -1: t2 = C.center2Node(t2, ['centers:VelocityX','centers:VelocityY','centers:VelocityZ'])
        else: t2 = P.computeVariables(t2, ['VelocityX', 'VelocityY', 'VelocityZ'])

    t2 = P.computeCurl(t2, ['VelocityX', 'VelocityY', 'VelocityZ'])
    if presvx == -1: t2 = C.rmVars(t2, ['VelocityX','VelocityY','VelocityZ'])
    Internal._renameNode(t2, 'rotx', 'VorticityX')
    Internal._renameNode(t2, 'roty', 'VorticityY')
    Internal._renameNode(t2, 'rotz', 'VorticityZ')
    return t2

#==============================================================================
# Norme de la vorticite en centres
#==============================================================================
def computeVorticityMagnitude(t):
    t2 = Internal.copyRef(t)
    presvx = C.isNamePresent(t2, 'centers:VorticityX')
    if presvx==-1: t2 = computeVorticity(t)
    t2 = C.magnitude(t2, ['centers:VorticityX', 'centers:VorticityY', 'centers:VorticityZ'])
    if presvx==-1: t2 = C.rmVars(t2, ['centers:VorticityX','centers:VorticityY','centers:VorticityZ'])
    Internal._renameNode(t2, 'magnitudeVorticityXVorticityYVorticityZ', 'VorticityMagnitude')
    return t2

#==============================================================================
# Critere Q en centres
#==============================================================================
def computeQCriterion(t):
    vars0 = ['centers:gradxVelocityX', 'centers:gradyVelocityX', 'centers:gradzVelocityX',
             'centers:gradxVelocityY', 'centers:gradyVelocityY', 'centers:gradzVelocityY',
             'centers:gradxVelocityZ', 'centers:gradyVelocityZ', 'centers:gradzVelocityZ']
    t2 = Internal.copyRef(t)
    presvx = C.isNamePresent(t2, 'VelocityX') 
    if presvx == -1: 
        presvxc = C.isNamePresent(t2, 'centers:VelocityX') 
        if presvxc > -1: t2 = C.center2Node(t2,['centers:VelocityX','centers:VelocityY','centers:VelocityZ'])
        else: t2 = P.computeVariables(t2, ['VelocityX', 'VelocityY', 'VelocityZ'])

    presgx = C.isNamePresent(t2, 'centers:gradxVelocityX')
    if presgx == -1:
        presgxn = C.isNamePresent(t2, 'gradxVelocityX')
        if presgxn > -1: t2 = C.center2Node(t2, vars0)
        else:
            t2 = P.computeGrad(t2, 'VelocityX')
            t2 = P.computeGrad(t2, 'VelocityY')
            t2 = P.computeGrad(t2, 'VelocityZ')

    if presvx == -1: t2 = C.rmVars(t2, ['VelocityX','VelocityY','VelocityZ'])
    t2 = C.initVars(t2, '{centers:QCriterion} = -0.5*({centers:gradxVelocityX}*{centers:gradxVelocityX}+{centers:gradyVelocityY}*{centers:gradyVelocityY}+{centers:gradzVelocityZ}*{centers:gradzVelocityZ}+2*{centers:gradyVelocityX}*{centers:gradxVelocityY}+2*{centers:gradzVelocityX}*{centers:gradxVelocityZ}+2*{centers:gradzVelocityY}*{centers:gradyVelocityZ})')
    if presgx == -1: t2 = C.rmVars(t2, vars0)
    return t2

#==============================================================================
# Tenseur des contraintes en centres
#==============================================================================
def computeShearStress(t, gamma=1.4, rgp=287.053,
                       Cs=110.4, mus=1.76e-5, Ts=273.15):
    
    vars1 = ['centers:gradxVelocityX', 'centers:gradyVelocityX', 'centers:gradzVelocityX',
             'centers:gradxVelocityY', 'centers:gradyVelocityY', 'centers:gradzVelocityY',
             'centers:gradxVelocityZ', 'centers:gradyVelocityZ', 'centers:gradzVelocityZ']

    t2 = Internal.copyRef(t)

    # Molecular Viscosity
    presmuc = C.isNamePresent(t2, 'centers:ViscosityMolecular')
    if presmuc == -1: 
        presmu = C.isNamePresent(t2, 'ViscosityMolecular')
        if presmu == -1: 
            t2 = P.computeVariables(t2, ['ViscosityMolecular'],
                                    gamma=gamma, rgp=rgp, Cs=Cs, mus=mus, Ts=Ts)
        t2 = C.node2Center(t2, 'ViscosityMolecular')
        if presmu == -1: t2 = C.rmVars(t2, ['ViscosityMolecular'])

    # Gradient of velocity
    presgx = C.isNamePresent(t2, 'centers:gradxVelocityX')
    if presgx == -1: 
        presgxn = C.isNamePresent(t2, 'gradxVelocityX')
        if presgxn > -1: t2 = C.center2Node(t2,vars1)
        else: 
            # No gradient found
            presvx = C.isNamePresent(t2, 'VelocityX')
            if presvx > -1:
                # Using values at nodes
                t2 = P.computeGrad(t2, 'VelocityX')
                t2 = P.computeGrad(t2, 'VelocityY')
                t2 = P.computeGrad(t2, 'VelocityZ')
            else:
                # Using value at cell center
                presvxc = C.isNamePresent(t2, 'centers:VelocityX')
                if presvxc == -1:
                    # No value at cell center -> creation of them
                    t2 = P.computeVariables(t2, ['centers:VelocityX', 'centers:VelocityY', 'centers:VelocityZ'])
                t2 = P.computeGrad(t2, 'centers:VelocityX')
                t2 = P.computeGrad(t2, 'centers:VelocityY')
                t2 = P.computeGrad(t2, 'centers:VelocityZ')
                # If no value intially -> delete the created values at cell center
                if presvxc == -1: t2 = C.rmVars(t2, ['centers:VelocityX','centers:VelocityY','centers:VelocityZ'])


    t2 = C.initVars(t2, '{centers:ShearStressXX}=-2./3.*{centers:ViscosityMolecular}*({centers:gradxVelocityX}+{centers:gradyVelocityY}+{centers:gradzVelocityZ})+2.*{centers:ViscosityMolecular}*{centers:gradxVelocityX}')
    t2 = C.initVars(t2, '{centers:ShearStressYY}=-2./3.*{centers:ViscosityMolecular}*({centers:gradxVelocityX}+{centers:gradyVelocityY}+{centers:gradzVelocityZ})+2.*{centers:ViscosityMolecular}*{centers:gradyVelocityY}')
    t2 = C.initVars(t2, '{centers:ShearStressZZ}=-2./3.*{centers:ViscosityMolecular}*({centers:gradxVelocityX}+{centers:gradyVelocityY}+{centers:gradzVelocityZ})+2.*{centers:ViscosityMolecular}*{centers:gradzVelocityZ}')
    t2 = C.initVars(t2, '{centers:ShearStressXY}={centers:ViscosityMolecular}*({centers:gradyVelocityX}+{centers:gradxVelocityY})')
    t2 = C.initVars(t2, '{centers:ShearStressXZ}={centers:ViscosityMolecular}*({centers:gradzVelocityX}+{centers:gradxVelocityZ})')
    t2 = C.initVars(t2, '{centers:ShearStressYZ}={centers:ViscosityMolecular}*({centers:gradzVelocityY}+{centers:gradyVelocityZ})')

    if presmuc==-1: t2 = C.rmVars(t2, ['centers:ViscosityMolecular'])
    if presgx ==-1: t2 = C.rmVars(t2,vars1)
    return t2
    

#-------------------------------------------------------------------------------
# INPUT: t: tree of skin/wall borders (velocity gradients must be defined yet)
#-------------------------------------------------------------------------------
def _computeWallShearStress(t):
    dimPb = Internal.getNodeFromName(t, 'EquationDimension')
    if dimPb is None: dimPb = 3
    else: dimPb = Internal.getValue(dimPb)

    if dimPb == 2: 
        vars1 = ['gradxVelocityX', 'gradyVelocityX', 'gradxVelocityY', 'gradyVelocityY'] 
    else: 
        vars1 = ['gradxVelocityX', 'gradyVelocityX', 'gradzVelocityX',\
        'gradxVelocityY', 'gradyVelocityY', 'gradzVelocityY',\
        'gradxVelocityZ', 'gradyVelocityZ', 'gradzVelocityZ']

    presgx = C.isNamePresent(t,vars1[0])
    if presgx == 1: loc = 'nodes'
    else:
        presgx = C.isNamePresent(t, 'centers:'+vars1[0])
        if presgx == 1: loc = 'centers'
        else:
           raise ValueError('gradxVelocity is required in tree.')
        for nov in range(len(vars1)): 
            vars1[nov]='centers:'+vars1[nov]

    [RoInf, RouInf, RovInf, RowInf, RoeInf, PInf, TInf, cvInf, MInf, ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf, Mus, Cs, Ts, Pr] = C.getState(t)
    gam1cv = (Gamma-1.)*cvInf
    RoUInf2I = 1./(RouInf*RouInf+RovInf*RovInf+RowInf*RowInf)
    betas = Mus*(Ts+Cs)/(Ts**(3./2.))

    # mu=betas*sqrt(temp)/(1+Cs/temp)
    if loc == 'centers':
        varsRM = ['centers:dummy']
        presmuc = C.isNamePresent(t, 'centers:ViscosityMolecular')
        if presmuc == -1:
            varsRM+=['centers:ViscosityMolecular'] # a detruire a la fin
            # temperature
            prestemp = C.isNamePresent(t,'centers:Temperature')
            if prestemp == -1:
                varsRM += ['centers:Temperature']
                P._computeVariables(t,['centers:Temperature'])

            C._initVars(t,'{centers:ViscosityMolecular}=%5.12f*sqrt({centers:Temperature})/(1.+%5.12f/{centers:Temperature})'%(betas,Cs))
        # ShearStress tensor
        if dimPb == 3:
            C._initVars(t,'{centers:dummy}={centers:gradxVelocityX}+{centers:gradyVelocityY}+{centers:gradzVelocityZ}')# du/dx+dv/dy+dw/dz
            C._initVars(t, '{centers:ShearStressXX}={centers:ViscosityMolecular}*(-2./3.*{centers:dummy}+2.*{centers:gradxVelocityX})')
            C._initVars(t, '{centers:ShearStressYY}={centers:ViscosityMolecular}*(-2./3.*{centers:dummy}+2.*{centers:gradyVelocityY})')
            C._initVars(t, '{centers:ShearStressZZ}={centers:ViscosityMolecular}*(-2./3.*{centers:dummy}+2.*{centers:gradzVelocityZ})')
            C._initVars(t, '{centers:ShearStressXY}={centers:ViscosityMolecular}*({centers:gradyVelocityX}+{centers:gradxVelocityY})')
            C._initVars(t, '{centers:ShearStressXZ}={centers:ViscosityMolecular}*({centers:gradxVelocityZ}+{centers:gradzVelocityX})')
            C._initVars(t, '{centers:ShearStressYZ}={centers:ViscosityMolecular}*({centers:gradzVelocityY}+{centers:gradyVelocityZ})')
        else:
            C._initVars(t,'{centers:dummy}={centers:gradxVelocityX}+{centers:gradyVelocityY}')# du/dx+dv/dy+dw/dz
            C._initVars(t, '{centers:ShearStressXX}={centers:ViscosityMolecular}*(-2./3.*{centers:dummy}+2.*{centers:gradxVelocityX})')
            C._initVars(t, '{centers:ShearStressYY}={centers:ViscosityMolecular}*(-2./3.*{centers:dummy}+2.*{centers:gradyVelocityY})')
            C._initVars(t, '{centers:ShearStressXY}={centers:ViscosityMolecular}*({centers:gradyVelocityX}+{centers:gradxVelocityY})')

        if presgx == -1: varsRM+=vars1
            
    else:
        varsRM = ['dummy']
        presmuc = C.isNamePresent(t,'ViscosityMolecular')
        if presmuc == -1:
            varsRM += ['ViscosityMolecular'] # a detruire a la fin
            # temperature
            prestemp = C.isNamePresent(t,'Temperature')
            if prestemp == -1:
                varsRM+=['Temperature']
                P._computeVariables(t,['Temperature'])
            C._initVars(t,'{ViscosityMolecular}=%5.12f*sqrt({Temperature})/(1.+%5.12f/{Temperature})'%(betas,Cs))
        # ShearStress tensor
        if dimPb == 3:
            C._initVars(t,'{dummy}={gradxVelocityX}+{gradyVelocityY}+{gradzVelocityZ}')# du/dx+dv/dy+dw/dz
            C._initVars(t, '{ShearStressXX}={ViscosityMolecular}*(-2./3.*{dummy}+2.*{gradxVelocityX})')
            C._initVars(t, '{ShearStressYY}={ViscosityMolecular}*(-2./3.*{dummy}+2.*{gradyVelocityY})')
            C._initVars(t, '{ShearStressZZ}={ViscosityMolecular}*(-2./3.*{dummy}+2.*{gradzVelocityZ})')
            C._initVars(t, '{ShearStressXY}={ViscosityMolecular}*({gradyVelocityX}+{gradxVelocityY})')
            C._initVars(t, '{ShearStressXZ}={ViscosityMolecular}*({gradxVelocityZ}+{gradzVelocityX})')
            C._initVars(t, '{ShearStressYZ}={ViscosityMolecular}*({gradzVelocityY}+{gradyVelocityZ})')
        else:
            C._initVars(t,'{dummy}={gradxVelocityX}+{gradyVelocityY}')# du/dx+dv/dy+dw/dz
            C._initVars(t, '{ShearStressXX}={ViscosityMolecular}*(-2./3.*{dummy}+2.*{gradxVelocityX})')
            C._initVars(t, '{ShearStressYY}={ViscosityMolecular}*(-2./3.*{dummy}+2.*{gradyVelocityY})')
            C._initVars(t, '{ShearStressXY}={ViscosityMolecular}*({gradyVelocityX}+{gradxVelocityY})')

        if presgx  ==-1: varsRM+=vars1
    # end of shear stress tensor
    C._rmVars(t,varsRM)
    return None

#==============================================================================
# Skin friction: input t must be a 2D mesh (skin)
#==============================================================================
def computeSkinFriction(t, centers=0, tangent=0):
    try: import Generator.PyTree as G
    except: raise ImportError("computeExtraVariable: skinFriction requires Generator module.")
    t2 = Internal.copyRef(t)
    dimPb = Internal.getNodeFromName(t,'EquationDimension')
    if dimPb is None: dimPb = 3
    else: dimPb = Internal.getValue(dimPb)

    # variable shearstress at cell centers
    centersShear = C.isNamePresent(t2, 'centers:ShearStressXX')
    if centersShear == 1: # centers
        if dimPb == 2:
            C._initVars(t2, "centers:ShearStressXZ", 0.)
            C._initVars(t2, "centers:ShearStressYZ", 0.)
            C._initVars(t2, "centers:ShearStressZZ", 0.)

        # normal vector for faces
        pres1 = C.isNamePresent(t2, 'centers:sx')
        if pres1 == -1: t2 = G.getNormalMap(t2)
        else:
            pres1n = C.isNamePresent(t2, 'sx')
            if pres1n > -1: t2 = C.node2Center(t2,['sx','sy','sz'])                
        t2 = C.normalize(t2, ['centers:sx', 'centers:sy', 'centers:sz'])

        pres2 = C.isNamePresent(t2, 'centers:SkinFrictionX')
        if pres2 == -1:
            pres2n = C.isNamePresent(t2, 'SkinFrictionX')
            if pres2n > -1: t2 = C.node2Center(t2,['SkinFrictionX','SkinFrictionY','SkinFrictionZ'])     
            else:
                t2 = C.initVars(t2, '{centers:SkinFrictionX}={centers:ShearStressXX}*{centers:sx}+{centers:ShearStressXY}*{centers:sy}+{centers:ShearStressXZ}*{centers:sz}')
                t2 = C.initVars(t2, '{centers:SkinFrictionY}={centers:ShearStressXY}*{centers:sx}+{centers:ShearStressYY}*{centers:sy}+{centers:ShearStressYZ}*{centers:sz}')
                t2 = C.initVars(t2, '{centers:SkinFrictionZ}={centers:ShearStressXZ}*{centers:sx}+{centers:ShearStressYZ}*{centers:sy}+{centers:ShearStressZZ}*{centers:sz}')

        if tangent == 1:
            pres3 = C.isNamePresent(t2, 'centers:SkinFrictionTangentialX')
            if pres3 == -1:
                pres3n = C.isNamePresent(t2,'SkinFrictionTangentialX')
                if pres3n > -1:
                    t2 = C.node2Center(t2,['SkinFrictionTangentialX','SkinFrictionTangentialY','SkinFrictionTangentialZ'])   
                else:
                    t2 = C.initVars(t2, '{centers:SkinFrictionTangentialX}={centers:SkinFrictionX} - {centers:sx}*({centers:SkinFrictionX}*{centers:sx}+{centers:SkinFrictionY}*{centers:sy}+{centers:SkinFrictionZ}*{centers:sz})')
                    t2 = C.initVars(t2, '{centers:SkinFrictionTangentialY}={centers:SkinFrictionY} - {centers:sy}*({centers:SkinFrictionX}*{centers:sx}+{centers:SkinFrictionY}*{centers:sy}+{centers:SkinFrictionZ}*{centers:sz})')
                    t2 = C.initVars(t2, '{centers:SkinFrictionTangentialZ}={centers:SkinFrictionZ} - {centers:sz}*({centers:SkinFrictionX}*{centers:sx}+{centers:SkinFrictionY}*{centers:sy}+{centers:SkinFrictionZ}*{centers:sz})')

            if pres2 == -1: t2 = C.rmVars(t2, ['centers:SkinFrictionX', 'centers:SkinFrictionY', 'centers:SkinFrictionZ'])
        if pres1 == -1: t2 = C.rmVars(t2, ['centers:sx', 'centers:sy', 'centers:sz'])
        if centers == 1: return t2
        else:
            if tangent == 1:
                pv = ['centers:SkinFrictionTangentialX', 'centers:SkinFrictionTangentialY', 'centers:SkinFrictionTangentialZ']
                t2 = C.center2Node(t2, pv)
                t2 = C.rmVars(t2, pv)
            else:
                pv = ['centers:SkinFrictionX', 'centers:SkinFrictionY', 'centers:SkinFrictionZ']
                t2 = C.center2Node(t2, pv)
                t2 = C.rmVars(t2, pv)
            return t2

    else: # node
        if dimPb == 2:
            C._initVars(t2, 'ShearStressXZ', 0.)
            C._initVars(t2, 'ShearStressYZ', 0.)
            C._initVars(t2, 'ShearStressZZ', 0.)

        pres3 = C.isNamePresent(t2, 'sx')
        if pres3 == -1:
            pres1 = C.isNamePresent(t, 'centers:sx')
            if pres1 == -1: t2 = G.getNormalMap(t2)
            t2 = C.center2Node(t2, ['centers:sx', 'centers:sy', 'centers:sz'])
            if pres1 == -1: t2 = C.rmVars(t2, ['centers:sx', 'centers:sy', 'centers:sz'])
        t2 = C.normalize(t2, ['sx', 'sy', 'sz'])

        # Array must contain ShearStress
        pres2 = C.isNamePresent(t2, 'SkinFrictionX')
        if pres2 == -1:
            pres2n = C.isNamePresent(t2, 'centers:SkinFrictionX')
            if pres2n > -1:
                t2 = C.center2Node(t2, ['centers:SkinFrictionX','centers:SkinFrictionY','centers:SkinFrictionZ'])
            else:
                t2 = C.initVars(t2, '{SkinFrictionX}={ShearStressXX}*{sx}+{ShearStressXY}*{sy}+{ShearStressXZ}*{sz}')
                t2 = C.initVars(t2, '{SkinFrictionY}={ShearStressXY}*{sx}+{ShearStressYY}*{sy}+{ShearStressYZ}*{sz}')
                t2 = C.initVars(t2, '{SkinFrictionZ}={ShearStressXZ}*{sx}+{ShearStressYZ}*{sy}+{ShearStressZZ}*{sz}')
        
        if tangent == 1:
            t2 = C.initVars(t2, '{SkinFrictionTangentialX}={SkinFrictionX} - {sx}*({SkinFrictionX}*{sx}+{SkinFrictionY}*{sy}+{SkinFrictionZ}*{sz})')
            t2 = C.initVars(t2, '{SkinFrictionTangentialY}={SkinFrictionY} - {sy}*({SkinFrictionX}*{sx}+{SkinFrictionY}*{sy}+{SkinFrictionZ}*{sz})')
            t2 = C.initVars(t2, '{SkinFrictionTangentialZ}={SkinFrictionZ} - {sz}*({SkinFrictionX}*{sx}+{SkinFrictionY}*{sy}+{SkinFrictionZ}*{sz})')
            if pres2 == -1:
                t2 = C.rmVars(t2, ['SkinFrictionX', 'SkinFrictionY', 'SkinFrictionZ'])
        if pres3 == -1: t2 = C.rmVars(t2, ['sx', 'sy', 'sz'])
        if centers == 1:
            if tangent == 1:
                pv = ['SkinFrictionTangentialX', 'SkinFrictionTangentialY', 'SkinFrictionTangentialZ']
                t2 = C.node2Center(t2, pv)
                t2 = C.rmVars(t2, pv)
                return t2
            else:
                pv = ['SkinFrictionX', 'SkinFrictionY', 'SkinFrictionZ']
                t2 = C.node2Center(t2, pv)
                t2 = C.rmVars(t2, pv)
                return t2
        else: return t2


#=========================================================================================
# Volumic extractions
#=========================================================================================

# Extract tree
# Simplifie l'arbre. Ne recupere que certaines variables.
# Recupere le refState
def extractTree(t, vars=['centers:Density', 'centers:VelocityX','centers:VelocityY','centers:VelocityZ','centers:Temperature', 'centers:TurbulentSANuTilde']):
    """Create a mirror tree with less vars."""
    tp = C.extractVars(t, vars, keepOldNodes=False)
    # Recupere le RefState
    refState = Internal.getNodeFromName2(t, 'ReferenceState')
    if refState is not None:
        for b in Internal.getBases(tp):
            b[2].append(refState)
    return tp

# IN: centers:Velocity
# OUT: centers:Vorticity
# Set ghostCells to True if t has ghost cells
def _computeVorticity2(t, ghostCells=False):
    """Compute vorticity for velocity in centers."""
    P._computeGrad2(t, 'centers:VelocityX', ghostCells)
    C._initVars(t, '{centers:VorticityX}=0.')
    C._initVars(t, '{centers:VorticityY}={centers:gradzVelocityX}')
    C._initVars(t, '{centers:VorticityZ}=-{centers:gradyVelocityX}')
    C._rmVars(t, ['centers:gradxVelocityX', 'centers:gradyVelocityX', 'centers:gradzVelocityX'])
    P._computeGrad2(t, 'centers:VelocityY', ghostCells)
    C._initVars(t, '{centers:VorticityX}={centers:VorticityX}-{centers:gradzVelocityY}')
    C._initVars(t, '{centers:VorticityZ}={centers:VorticityZ}+{centers:gradxVelocityY}')
    C._rmVars(t, ['centers:gradxVelocityY', 'centers:gradyVelocityY', 'centers:gradzVelocityY'])
    P._computeGrad2(t, 'centers:VelocityZ', ghostCells)
    C._initVars(t, '{centers:VorticityX}={centers:VorticityX}+{centers:gradyVelocityZ}')
    C._initVars(t, '{centers:VorticityY}={centers:VorticityY}-{centers:gradxVelocityZ}')
    C._rmVars(t, ['centers:gradxVelocityZ', 'centers:gradyVelocityZ', 'centers:gradzVelocityZ'])
    return None

# IN: centers:Velocity 
# OUT: centers:VorticityMagnitude
def _computeVorticityMagnitude2(t, ghostCells=False):
    """Compute vorticity magnitude for velocity in centers."""    
    _computeVorticity2(t, ghostCells)
    C._magnitude(t, ['centers:VorticityX', 'centers:VorticityY', 'centers:VorticityZ'])
    C._rmVars(t, ['centers:VorticityX', 'centers:VorticityY', 'centers:VorticityZ'])
    Internal._renameNode(t, 'magnitudeVorticityXVorticityYVorticityZ', 'VorticityMagnitude')
    return None

# IN: centers:Velocity 
# OUT: centers:QCriterion
def _computeQCriterion2(t, ghostCells=False):
    """Compute Q criterion for velocity in centers."""
    P._computeGrad2(t, 'centers:VelocityX', ghostCells)
    C._initVars(t, '{centers:QCriterion}=-0.5*{centers:gradxVelocityX}*{centers:gradxVelocityX}')
    C._rmVars(t, ['centers:gradxVelocityX'])
    P._computeGrad2(t, 'centers:VelocityY', ghostCells)
    C._initVars(t, '{centers:QCriterion}={centers:QCriterion}-0.5*{centers:gradyVelocityY}*{centers:gradyVelocityY}-{centers:gradyVelocityX}*{centers:gradxVelocityY}')
    C._rmVars(t, ['centers:gradyVelocityY', 'centers:gradyVelocityX', 'centers:gradxVelocityY'])
    P._computeGrad2(t, 'centers:VelocityZ', ghostCells)
    C._initVars(t, '{centers:QCriterion}={centers:QCriterion}-0.5*{centers:gradzVelocityZ}*{centers:gradzVelocityZ}-{centers:gradzVelocityX}*{centers:gradxVelocityZ}-{centers:gradzVelocityY}*{centers:gradyVelocityZ}')
    C._rmVars(t, ['centers:gradxVelocityZ', 'centers:gradzVelocityZ', 'centers:gradzVelocityX', 'centers:gradzVelocityY', 'centers:gradzVelocityX', 'centers:gradyVelocityZ'])
    return None
    
# Calcul de la pression dans le volume
# IN: centers:Density
# IN: centers:Temperature
# IN: cv: refState (r a partir de cv)
# OUT: centers:Pressure
# P = ro r T
def _extractPressure(t):
    """Extract Pressure."""
    P._computeVariables2(t, ['centers:Pressure'])
    return None

# Calcul module de la vitesse 
# IN: centers:VelocityX, centers:VelocityY, centers:VelocityZ
# OUT: centers:VelocityMagnitude
# v = sqrt(vx**2+vy**2+vz**2)
def _extractVelocityMagnitude(t):
    """Extract velocity magnitude."""
    P._computeVariables2(t, ['centers:VelocityMagnitude'])

# Mach volumique
# IN: centers:Temperature
# IN: centers:VelocityX, centers:VelocityY, centers:VelocityZ
# IN: cv: refState (r a partir de cv)
# OUT: centers:Mach
# M = u/sqrt(gamma p/ro) - p = ro r T
def _extractMach(t):
    """Extract Mach."""
    P._computeVariables2(t, ['centers:Mach'])

# Mu volumique
# mu = sutherland(T)
# IN: centers:Temperature
# IN: Cs, Mus, Ts from refstate
# OUT: centers:ViscosityMolecular
def _extractViscosityMolecular(t):
    """Extract Viscosity molecular."""
    P._computeVariables2(t, ['centers:ViscosityMolecular'])

# mut from spallart
# IN: centers:TurulentSANuTilde
# IN: centers:Density
# OUT: centers:ViscosityEddy
# kappa = ro * nutilde / mu
# mut = ro * nutilde * kappa^3 / (kappa^3 + 7.1^3)
def _extractViscosityEddy(t):
    """Extract eddy viscosity."""
    C._initVars(t, '{centers:Kappa} = {centers:Density} * {centers:TurbulentSANuTilde} / {centers:ViscosityMolecular}')
    C._initVars(t, '{centers:ViscosityEddy} = {centers:Density}*{centers:TurbulentSANuTilde} * {centers:Kappa}**3 / ( {centers:Kappa}**3 + 7.1**3 )')
    C._rmVars(t, 'centers:Kappa')
    return None

# mut sur mu
# IN: centers:ViscosityEddy (mut)
# IN: centers:ViscosityMolecular (mu)
# OUT: centers:MutSurMu
# mutSurmu = mut / mu
def _extractMutSurMu(t):
    """Extract Mut over mu."""
    C._initVars(t, '{centers:MutSurMu} = {centers:ViscosityEddy}/{centers:ViscosityMolecular}')

#======================================================
# Volumic extractions with gradient
# input : avec ghost cells
#======================================================

# IN: centers:VelocityX, centers:VelocityY, centers:VelocityZ
# OUT: centers:VorticityMagnitude
def _extractVorticityMagnitude(t):
    """Extract vorticity magnitude."""
    _computeVorticityMagnitude2(t, ghostCells=True)
    return None

# IN: centers:VelocityX, centers:VelocityY, centers:VelocityZ
# OUT: centers:VorticityX, centers:VorticityY, centers:VorticityZ
def _extractVorticity(t):
    """Extract vorticity."""
    _computeVorticity2(t, ghostCells=True)
    return None

# IN: centers:VelocityX, centers:VelocityY, centers:VelocityZ
# OUT: centers:QCriterion
def _extractQCriterion(t):
    """Extract Q criterion."""
    _computeQCriterion2(t, ghostCells=True)
    return None

#===========================================================
# Surfacic extractions 
# input : effort tree (teff)
#===========================================================

# calcul les efforts aerodynamiques
# IN: centers:MomentumX,centers:MomentumY,centers:MomentumZ (flux de quantite de mouvement)
# OUT: centers:Fx, centers:Fy, centers:Fz
# Attention, il y a -Pinf dans le flux et ils sont deja multiplies par 0.5*roinf*uinf^2
# De plus, ils sont deja integres sur la surface de la cellule
# Pour pouvoir ensuite l'integre avec P.integ, on multiplie ici par 1./s
def _extractForceLoads(teff):
    """Extract forces."""
    G._getNormalMap(teff)
    C._magnitude(teff, ['centers:sx', 'centers:sy', 'centers:sz'])
    C._initVars(teff, '{centers:Fx}={centers:MomentumX}/{centers:sMagnitude}')
    C._initVars(teff, '{centers:Fy}={centers:MomentumY}/{centers:sMagnitude}')
    C._initVars(teff, '{centers:Fz}={centers:MomentumZ}/{centers:sMagnitude}')
    C._rmVars(teff, ['centers:sx','centers:sy','centers:sz','sMagnitude'])
    return None

# extract shear stress
# IN: centers:ViscosityMolecular (mu)
# IN: centers:gradxVelocityX, centers:gradxVelocityY, ...
# OUT: centers:ShearStressXX,YY,ZZ,XY,XZ,YZ (tenseur symetrique)
# tau = -2/3 mu div u I + 2 mu D
# Attention : dans teff, c'est ViscosityMolecular+ViscosityEddy qui est dans la variable ViscosityMolecular
# ce qui est correct ici
def _extractShearStress(teff):
    """Extract shearStress."""
    C._initVars(teff,'{centers:divu}={centers:gradxVelocityX}+{centers:gradyVelocityY}+{centers:gradzVelocityZ}') # du/dx+dv/dy+dw/dz
    C._initVars(teff, '{centers:ShearStressXX}={centers:ViscosityMolecular}*(-2./3.*{centers:divu}+2.*{centers:gradxVelocityX})')
    C._initVars(teff, '{centers:ShearStressYY}={centers:ViscosityMolecular}*(-2./3.*{centers:divu}+2.*{centers:gradyVelocityY})')
    C._initVars(teff, '{centers:ShearStressZZ}={centers:ViscosityMolecular}*(-2./3.*{centers:divu}+2.*{centers:gradzVelocityZ})')
    C._initVars(teff, '{centers:ShearStressXY}={centers:ViscosityMolecular}*({centers:gradyVelocityX}+{centers:gradxVelocityY})')
    C._initVars(teff, '{centers:ShearStressXZ}={centers:ViscosityMolecular}*({centers:gradxVelocityZ}+{centers:gradzVelocityX})')
    C._initVars(teff, '{centers:ShearStressYZ}={centers:ViscosityMolecular}*({centers:gradzVelocityY}+{centers:gradyVelocityZ})')
    C._rmVars(teff, 'centers:divu')
    return None

# Extract tangential friction vector
# IN: centers:shearStressXX,...
# OUT: centers:frictionX, centers:frictionY, centers:frictionZ
# taut = (n. tau.n) n - tau.n
def _extractFrictionVector(teff):
    """Extract tangiential firction vector."""
    G._getNormalMap(teff)
    C._normalize(teff, ['centers:sx', 'centers:sy','centers:sz'])
    C._initVars(teff, '{centers:frictionX}={centers:ShearStressXX}*{centers:sx}+{centers:ShearStressXY}*{centers:sy}+{centers:ShearStressXZ}*{centers:sz}')
    C._initVars(teff, '{centers:frictionY}={centers:ShearStressXY}*{centers:sx}+{centers:ShearStressYY}*{centers:sy}+{centers:ShearStressYZ}*{centers:sz}')
    C._initVars(teff, '{centers:frictionZ}={centers:ShearStressXZ}*{centers:sx}+{centers:ShearStressYZ}*{centers:sy}+{centers:ShearStressZZ}*{centers:sz}')
    C._initVars(teff, '{centers:scal} = {centers:frictionX}*{centers:sx}+{centers:frictionY}*{centers:sy}+{centers:frictionZ}*{centers:sz}')
    C._initVars(teff, '{centers:frictionX}={centers:scal}*{centers:sx}-{centers:frictionX}')
    C._initVars(teff, '{centers:frictionY}={centers:scal}*{centers:sy}-{centers:frictionY}')
    C._initVars(teff, '{centers:frictionZ}={centers:scal}*{centers:sz}-{centers:frictionZ}')
    C._rmVars(teff, ['centers:scal', 'centers:sx', 'centers:sy', 'centers:sz'])
    return None

# tauw
# IN: centers:frictionX, centers:frictionY, centers:frictionZ
# OUT: centers: frictionMagnitude
# tauw = ||taut||
def _extractFrictionMagnitude(teff):
    """Extract friction magnitude."""
    C._magnitude(teff, ['centers:frictionX', 'centers:frictionY', 'centers:frictionZ'])
    return None

# IN: centers:frictionMagnitude
# IN: centers:Density2 - density
# OUT: centers:utau
# utau = sqrt(tauw/ro)
def _extractUTau(teff):
    """Extract utau."""
    C._initVars(teff, '{centers:utau}=sqrt({centers:frictionMagnitude}/{centers:Density2})')
    return None

# IN: centers:frictionMagnitude
# IN: norm : (1/2*roinf*uinf**2)
# OUT: centers:Cf
# cf = tauw / (1/2*roinf*uinf**2)
def _extractCf(teff, norm):
    """Extract Cf."""
    C._initVars(teff, '{centers:Cf} = {centers:frictionMagnitude} / %20.16g'%norm)
    return None

# IN: centers:Pressure
# IN: pinf: infinite pressure
# IN: norm : (1/2*roinf*uinf**2)
# OUT: centers:Cp
# cp = (p-pinf) / (1/2*roinf*uinf**2)
def _extractCp(teff, pinf, norm):
    """Extract Cp."""
    C._initVars(teff, '{centers:Cp} = ({centers:Pressure}-%20.16g) / %20.16g'%(pinf,norm))
    return None

# integration des pressions Integ( p.n ds )
# IN: centers:Pressure
def integPressure(teff):
    """Integ p.n.ds"""
    ret = P.integNorm(teff, 'centers:Pressure')
    return ret[0]

# integration coefficient de pression Integ( Cp.n ds )
def integCp(teff):
    """Integ Cp.n.ds"""
    ret = P.integNorm(teff, 'centers:Cp')
    return ret[0]

# integration Integ ( tau.n ds )
# IN: centers:ShearStress
def integTaun(teff):
    """Integ tau.n.ds"""
    G._getNormalMap(teff)
    C._normalize(teff, ['centers:sx', 'centers:sy','centers:sz'])
    C._initVars(teff, '{centers:taunx} = {centers:ShearStressXX}*{centers:sx}+{centers:ShearStressXY}*{centers:sy}+{centers:ShearStressXZ}*{centers:sz}')
    C._initVars(teff, '{centers:tauny} = {centers:ShearStressXY}*{centers:sx}+{centers:ShearStressYY}*{centers:sy}+{centers:ShearStressYZ}*{centers:sz}')
    C._initVars(teff, '{centers:taunz} = {centers:ShearStressXZ}*{centers:sx}+{centers:ShearStressYZ}*{centers:sy}+{centers:ShearStressZZ}*{centers:sz}')
    retx = P.integ(teff, 'centers:taunx')
    rety = P.integ(teff, 'centers:tauny')
    retz = P.integ(teff, 'centers:taunz')
    C._rmVars(teff, ['centers:sx', 'centers:sy', 'centers:sz', 'centers:taunx', 'centers:tauny', 'centers:taunz'])
    return [-retx[0],-rety[0],-retz[0]]

# Integration des Cf Integ( Cf. ds )
def integCf(teff):
    """Integ Cf.ds"""
    ret = P.integ(teff, 'centers:Cf')
    return ret[0]

# Attention : deja -pinf dans les Fx,Fy,Fz
# et deja division par 0.5*roinf*uinf^2
# retourne donc des charges adimensionnees
def integLoads(teff):
    """Integ F.ds"""
    retx = P.integ(teff, 'centers:Fx')
    rety = P.integ(teff, 'centers:Fy')
    retz = P.integ(teff, 'centers:Fz')
    return [retx[0],rety[0],retz[0]]

#=============================================================
# Profile extractions 
# input : effort tree (teff) + volume tree t
#=============================================================

# Extrait des profiles proche paroi en i,j,k (deux grandeurs constantes)
# IN: zonePath: chemin dela zone pour extraire
# IN: i,j de la ligne
def extractProfile(t, zonePath, i=-1, j=-1, k=-1):
    """Extract given profile."""
    z = Internal.getNodeFromPath(t, zonePath)
    #zp = T.subzone(z, (i,j,3), (i,j,-1))
    if k == -1: # suivant k
        zp = T.subzone(z, (i,j,3), (i,j,-1))
    elif j == -1:
        zp = T.subzone(z, (i,3,k), (i,-1,k))
    else:
        zp = T.subzone(z, (3,j,k), (-1,j,k))

    Internal._rmNodesFromName(zp, 'ZoneBC')
    zp = C.node2Center(zp)
    L = D.getLength(zp)
    D._getCurvilinearAbscissa(zp)
    C._initVars(zp, '{yw}={s}*%20.16g'%L)
    return zp

# Calcul y+ et u+ sur le profil
# IN: teff: with utau computed
def _extractYPlus(zp, teff):
    """Extract y+ and u+."""
    # Extract point wall de zp
    Pwall = T.subzone(zp, (1,1,1), (1,1,1))
    # localise Pwall sur teff
    zeff = Internal.getZones(teff)[0]
    hook = C.createHook(zeff, function='nodes')
    nodes, dist = C.nearestNodes(hook, Pwall)
    # Get tauw
    n = nodes[0]
    utau = C.getValue(zeff, 'utau', n)
    mu = C.getValue(zeff, 'ViscosityMolecular', n)
    ro = C.getValue(zeff, 'Density2', n)
    print('INFO: found wall point ro=%g, utau=%g, mu=%g, index=%d\n'%(ro, utau, mu, n))
    C._initVars(zp, '{yPlus}={yw} * %20.16g * %20.16g / %20.16g'%(ro,utau,mu))
    C._initVars(zp, '{yPlus}=log10({yPlus})')
    C._initVars(zp, '{uPlus}={VelocityMagnitude}/%20.16g'%utau)
    return None
