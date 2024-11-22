"""Compute other complex variables than computeVariables."""

from . import PyTree as P
from . import Mpi as Pmpi
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import Transform.PyTree as T
import Geom.PyTree as D
import numpy

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
    if presvx == -1: t2 = computeVorticity(t)
    t2 = C.magnitude(t2, ['centers:VorticityX', 'centers:VorticityY', 'centers:VorticityZ'])
    if presvx == -1: t2 = C.rmVars(t2, ['centers:VorticityX','centers:VorticityY','centers:VorticityZ'])
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