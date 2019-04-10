# Calcul des autres variables de computeVariables

from . import Post as P
import Converter as C
import KCore

#==============================================================================
# Vorticite en centres
#==============================================================================
def computeVorticity(array):
    presvx = KCore.isNamePresent(array, 'VelocityX') 
    if presvx == -1: 
        ret = P.computeVariables(array, ['VelocityX', 'VelocityY', 'VelocityZ'])
        a = C.addVars([array, ret])
        rot = P.computeCurl(a, ['VelocityX', 'VelocityY', 'VelocityZ'])
    else: rot = P.computeCurl(array, ['VelocityX', 'VelocityY', 'VelocityZ'])

    if isinstance(rot[0], list):
        for r in rot:
            r[0] = 'VorticityX,VorticityY,VorticityZ'
    else: rot[0] = 'VorticityX,VorticityY,VorticityZ'
    return rot

#==============================================================================
# Norme de la vorticite en centres
#==============================================================================
def computeVorticityMagnitude(array):
    presvx = KCore.isNamePresent(array, 'centers:VorticityX') 
    if presvx == -1:
        ret = P.computeVariables(array, ['VelocityX', 'VelocityY', 'VelocityZ'])
        a = C.addVars([array, ret])
        rot = P.computeCurl(a, ['VelocityX', 'VelocityY', 'VelocityZ'])
    else: rot = P.computeCurl(array, ['VelocityX', 'VelocityY', 'VelocityZ'])
        
    if isinstance(rot[0], list):
        for r in rot:
            r[0] = 'VorticityX,VorticityY,VorticityZ'
    else: rot[0] = 'VorticityX,VorticityY,VorticityZ'
    rot = C.magnitude(rot, ['VorticityX', 'VorticityY', 'VorticityZ'])
    return rot

#==============================================================================
# Calcul du critere Q aux centres
#==============================================================================
def computeQCriterion(array):
    vars0 = ['gradxVelocityX', 'gradyVelocityX', 'gradzVelocityX',
             'gradxVelocityY', 'gradyVelocityY', 'gradzVelocityY',
             'gradxVelocityZ', 'gradyVelocityZ', 'gradzVelocityZ']
    # Compute velocity
    presvx = KCore.isNamePresent(array, 'VelocityX') 
    if presvx==-1: 
        ret = P.computeVariables(array, ['VelocityX', 'VelocityY', 'VelocityZ'])
        a = C.addVars([array, ret])
    else: a = array
    # compute gradient of velocity
    presgx = KCore.isNamePresent(a,'gradxVelocityX')
    if presgx == -1:
        gradU = P.computeGrad(a, 'VelocityX')
        gradV = P.computeGrad(a, 'VelocityY')
        gradW = P.computeGrad(a, 'VelocityZ')
        grads = C.addVars([gradU, gradV, gradW])
    else: grads = C.extractVars(a,vars0)
    
    grads = C.initVars(grads, '{QCriterion}=-0.5*({gradxVelocityX}*{gradxVelocityX}+{gradyVelocityY}*{gradyVelocityY}+{gradzVelocityZ}*{gradzVelocityZ}+2*{gradyVelocityX}*{gradxVelocityY}+2*{gradzVelocityX}*{gradxVelocityZ}+2*{gradzVelocityY}*{gradyVelocityZ})')
    ret = C.extractVars(grads, ['QCriterion'])
    return ret

#==============================================================================
# Tenseur des contraintes en centres
#==============================================================================
def computeShearStress(array, gamma=1.4, rgp=287.053,
                       Cs=110.4, mus=1.76e-5, Ts=273.15):
    vars0 = ['gradxVelocityX', 'gradyVelocityX', 'gradzVelocityX',
             'gradxVelocityY', 'gradyVelocityY', 'gradzVelocityY',
             'gradxVelocityZ', 'gradyVelocityZ', 'gradzVelocityZ']
    vars1 = ['ViscosityMolecular']+vars0
    presvx = KCore.isNamePresent(array, 'VelocityX')
    if presvx == -1: 
        ret = P.computeVariables(array, ['VelocityX', 'VelocityY', 'VelocityZ'])
        a = C.addVars([array, ret])
    else: a = array

    presmu = KCore.isNamePresent(array, 'ViscosityMolecular')
    if presmu == -1: mu = P.computeVariables(array, ['ViscosityMolecular'],
                                             gamma=gamma, rgp=rgp, Cs=Cs, mus=mus, Ts=Ts)
    mu = C.node2Center(mu)
    
    presgx = KCore.isNamePresent(a,'gradxVelocityX')
    if presgx == -1:
        gradU = P.computeGrad(a, 'VelocityX')
        gradV = P.computeGrad(a, 'VelocityY')
        gradW = P.computeGrad(a, 'VelocityZ')
        tau = C.addVars([gradU, gradV, gradW, mu])
    else: 
        grads =  C.extractVars(a,vars0)
        tau = C.addVars([grads, mu])

    tau = C.initVars(tau, '{ShearStressXX}=-2./3.*{ViscosityMolecular}*({gradxVelocityX}+{gradyVelocityY}+{gradzVelocityZ})+{ViscosityMolecular}*2*{gradxVelocityX}')
    tau = C.initVars(tau, '{ShearStressXY}={ViscosityMolecular}*({gradyVelocityX}+{gradxVelocityY})')
    tau = C.initVars(tau, '{ShearStressXZ}={ViscosityMolecular}*({gradzVelocityX}+{gradxVelocityZ})')
    tau = C.initVars(tau, '{ShearStressYX}={ViscosityMolecular}*({gradxVelocityY}+{gradyVelocityX})')
    tau = C.initVars(tau, '{ShearStressYY}=-2./3.*{ViscosityMolecular}*({gradxVelocityX}+{gradyVelocityY}+{gradzVelocityZ})+{ViscosityMolecular}*2*{gradyVelocityY}')
    tau = C.initVars(tau, '{ShearStressYZ}={ViscosityMolecular}*({gradzVelocityY}+{gradyVelocityZ})')
    tau = C.initVars(tau, '{ShearStressZX}={ViscosityMolecular}*({gradxVelocityZ}+{gradzVelocityX})')
    tau = C.initVars(tau, '{ShearStressZY}={ViscosityMolecular}*({gradyVelocityZ}+{gradzVelocityY})')
    tau = C.initVars(tau, '{ShearStressZZ}=-2./3.*{ViscosityMolecular}*({gradxVelocityX}+{gradyVelocityY}+{gradzVelocityZ})+{ViscosityMolecular}*2*{gradzVelocityZ}')
    ret = C.extractVars(tau,
                        ['ShearStressXX', 'ShearStressXY', 'ShearStressXZ',
                         'ShearStressYX', 'ShearStressYY', 'ShearStressYZ',
                         'ShearStressZX', 'ShearStressZY', 'ShearStressZZ'])
    return ret

#==============================================================================
# Skin friction: array must be a surface array
#==============================================================================
def computeSkinFriction(array, tangent=0):
    pres1 = KCore.isNamePresent(array, 'sx')
    if pres1 == -1: 
        try: import Generator as G
        except: raise ImportError("computeExtraVariable: skinFriction requires Generator module.")
        n = G.getNormalMap(array)
        n = C.center2Node(n)
    else: n = C.extractVars(array,['sx','sy','sz'])
    n = C.normalize(n, ['sx', 'sy', 'sz'])
    a = C.addVars([array, n])
    # Array must contain ShearStress
    a = C.initVars(a, '{SkinFrictionX}={ShearStressXX}*{sx}+{ShearStressXY}*{sy}+{ShearStressXZ}*{sz}')
    a = C.initVars(a, '{SkinFrictionY}={ShearStressYX}*{sx}+{ShearStressYY}*{sy}+{ShearStressYZ}*{sz}')
    a = C.initVars(a, '{SkinFrictionZ}={ShearStressZX}*{sx}+{ShearStressZY}*{sy}+{ShearStressZZ}*{sz}')
    a = C.extractVars(a, ['SkinFrictionX', 'SkinFrictionY', 'SkinFrictionZ'])
    if tangent == 1:
        a = C.addVars([a, n])
        a = C.initVars(a, '{SkinFrictionTangentialX}={SkinFrictionX} - {sx}*({SkinFrictionX}*{sx}+{SkinFrictionY}*{sy}+{SkinFrictionZ}*{sz})')
        a = C.initVars(a, '{SkinFrictionTangentialY}={SkinFrictionY} - {sy}*({SkinFrictionX}*{sx}+{SkinFrictionY}*{sy}+{SkinFrictionZ}*{sz})')
        a = C.initVars(a, '{SkinFrictionTangentialZ}={SkinFrictionZ} - {sz}*({SkinFrictionX}*{sx}+{SkinFrictionY}*{sy}+{SkinFrictionZ}*{sz})')
        a = C.extractVars(a, ['SkinFrictionTangentialX', 'SkinFrictionTangentialY', 'SkinFrictionTangentialZ'])
    return a
