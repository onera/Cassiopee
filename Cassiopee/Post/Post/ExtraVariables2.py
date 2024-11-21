"""Compute other derived variables from primitive variables."""

from . import PyTree as P
from . import Mpi as Pmpi
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import Transform.PyTree as T
import Geom.PyTree as D
import numpy

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
def computeVorticity2(t, ghostCells=False):
    """Compute vorticity from velocity in centers."""
    tp = Internal.copyRef(t)
    _computeVorticity2(tp, ghostCells=ghostCells)
    return tp

def _computeVorticity2(t, ghostCells=False):
    """Compute vorticity from velocity in centers."""
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
def computeVorticityMagnitude2(t, ghostCells=False):
    """Compute vorticity magnitude from velocity in centers."""
    tp = Internal.copyRef(t)
    _computeVorticityMagnitude2(tp, ghostCells=ghostCells)
    return tp

def _computeVorticityMagnitude2(t, ghostCells=False):
    """Compute vorticity magnitude from velocity in centers."""    
    _computeVorticity2(t, ghostCells)
    C._magnitude(t, ['centers:VorticityX', 'centers:VorticityY', 'centers:VorticityZ'])
    C._rmVars(t, ['centers:VorticityX', 'centers:VorticityY', 'centers:VorticityZ'])
    Internal._renameNode(t, 'magnitudeVorticityXVorticityYVorticityZ', 'VorticityMagnitude')
    return None

# IN: centers:Velocity
# OUT: centers:QCriterion
def computeQCriterion2(t, ghostCells=False):
    """Compute Q criterion from velocity in centers."""
    tp = Internal.copyRef(t)
    _computeQCriterion2(tp, ghostCells=ghostCells)
    return tp

def _computeQCriterion2(t, ghostCells=False):
    """Compute Q criterion from velocity in centers."""
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

# IN: centers:Velocity 
# OUT: centers:lambda2
def computeLambda2(t, ghostCells=False):
    """Compute lambda2 criterion from velocity in centers."""
    tp = Internal.copyRef(t)
    _computeLambda2(tp, ghostCells=ghostCells)
    return tp

def _computeLambda2(t, ghostCells=False):
    """Compute lambda2 criterion for velocity in centers."""
    P._computeGrad2(t, 'centers:VelocityX', ghostCells)
    P._computeGrad2(t, 'centers:VelocityY', ghostCells)
    P._computeGrad2(t, 'centers:VelocityZ', ghostCells)

    # Rate of strain tensor
    C._initVars(t, 'centers:s11 = {centers:gradxVelocityX}')   
    C._initVars(t, 'centers:s12 = 0.5*({centers:gradyVelocityX}+{centers:gradxVelocityY})')
    C._initVars(t, 'centers:s13 = 0.5*({centers:gradzVelocityX}+{centers:gradxVelocityZ})')
    C._initVars(t, 'centers:s22 = {centers:gradyVelocityY}')
    C._initVars(t, 'centers:s23 = 0.5*({centers:gradzVelocityY}+{centers:gradyVelocityZ})')
    C._initVars(t, 'centers:s33 = {centers:gradzVelocityZ}')

    # Vorticity tensor
    C._initVars(t, 'centers:o12 = 0.5*({centers:gradyVelocityX}-{centers:gradxVelocityY})')
    C._initVars(t, 'centers:o13 = 0.5*({centers:gradzVelocityX}-{centers:gradxVelocityZ})')
    C._initVars(t, 'centers:o23 = 0.5*({centers:gradzVelocityY}-{centers:gradyVelocityZ})')

    # SikSkj + OikOkj
    C._initVars(t, 'centers:a11 = {centers:s11}**2+{centers:s12}**2+{centers:s13}**2 - ({centers:o12}**2+{centers:o13}**2)')
    C._initVars(t, 'centers:a22 = {centers:s12}**2+{centers:s22}**2+{centers:s23}**2 - ({centers:o12}**2+{centers:o23}**2)')
    C._initVars(t, 'centers:a33 = {centers:s13}**2+{centers:s23}**2+{centers:s33}**2 - ({centers:o13}**2+{centers:o23}**2)')
    C._initVars(t, 'centers:a12 = {centers:s11}*{centers:s12}+{centers:s12}*{centers:s22}+{centers:s13}*{centers:s23} - {centers:o13}*{centers:o23}')
    C._initVars(t, 'centers:a13 = {centers:s11}*{centers:s13}+{centers:s12}*{centers:s23}+{centers:s13}*{centers:s33} - {centers:o12}*{centers:o23}')
    C._initVars(t, 'centers:a23 = {centers:s12}*{centers:s13}+{centers:s22}*{centers:s23}+{centers:s23}*{centers:s33} - {centers:o12}*{centers:o13}')
    C._initVars(t, 'centers:Lambda2 = 0.')

    zones = Internal.getZones(t)
    for z in zones:
        nc = C.getNCells(z)
        AAA = numpy.empty( (nc,3,3), dtype=numpy.float64)
        p = Internal.getNodeFromName(z, 'a11')[1].ravel('k')
        AAA[:,0,0] = p[:]
        p = Internal.getNodeFromName(z, 'a22')[1].ravel('k')
        AAA[:,1,1] = p[:]
        p = Internal.getNodeFromName(z, 'a33')[1].ravel('k')
        AAA[:,2,2] = p[:]
        p = Internal.getNodeFromName(z, 'a12')[1].ravel('k')
        AAA[:,0,1] = p[:]
        AAA[:,1,0] = p[:]
        p = Internal.getNodeFromName(z, 'a13')[1].ravel('k')
        AAA[:,0,2] = p[:]
        AAA[:,2,0] = p[:]
        p = Internal.getNodeFromName(z, 'a23')[1].ravel('k')
        AAA[:,1,2] = p[:]
        AAA[:,2,1] = p[:]

        # Criterion : second eigenvalue (real and sorted) < 0
        lambda2 = numpy.linalg.eigvals(AAA)
        AAA = None
        lambda2 = numpy.real(lambda2)
        lambda2 = numpy.sort(lambda2, axis=1)

        s = Internal.getNodeFromName2(z, 'a13')[1].shape
        lambda2out = numpy.empty( (s), dtype=numpy.float64, order='F')
        lambda2out.ravel('k')[:] = lambda2[:,1]
        p = Internal.getNodeFromName2(z, 'Lambda2')
        p[1] = lambda2out

    C._rmVars(t, ['centers:gradxVelocityX', 'centers:gradyVelocityX', 'centers:gradzVelocityX'])
    C._rmVars(t, ['centers:gradxVelocityY', 'centers:gradyVelocityY', 'centers:gradzVelocityY'])
    C._rmVars(t, ['centers:gradxVelocityZ', 'centers:gradyVelocityZ', 'centers:gradzVelocityZ'])
    C._rmVars(t, ['centers:a11', 'centers:a22', 'centers:a33', 'centers:a12', 'centers:a13', 'centers:a23'])
    C._rmVars(t, ['centers:a11', 'centers:a22', 'centers:a33', 'centers:a12', 'centers:a13', 'centers:a23'])
    C._rmVars(t, ['centers:s11', 'centers:s22', 'centers:s33', 'centers:s12', 'centers:s13', 'centers:s23'])
    C._rmVars(t, ['centers:o12', 'centers:o13', 'centers:o23'])
    return None

def computeLogGradField2(t, name, ghostCells=False):
    """Compute log(grad field) for field in centers."""
    tp = Internal.copyRef(t)
    _computeLogGradField2(tp, name, ghostCells=ghostCells)
    return tp

# IN: centers:Field
# OUT: centers:logGradField (log10)
def _computeLogGradField2(t, name, ghostCells=False):
    """Compute log(grad field) for field in centers."""
    sp = name.split(':')
    if len(sp) == 2: name = sp[1]
    P._computeGrad2(t, 'centers:'+name, ghostCells)
    C._magnitude(t, ['centers:gradx'+name, 'centers:grady'+name, 'centers:gradz'+name])
    C._rmVars(t, ['centers:gradx'+name, 'centers:grady'+name, 'centers:gradz'+name])
    for z in Internal.getZones(t):
        n = Internal.getNodeFromName(z, 'gradx'+name+'grady'+name+'gradz'+name+'Magnitude')
        p = numpy.clip(n[1], 1.e-12, None)
        n[1] = numpy.log10(p)
        n[0] = 'LogGrad'+name
    return None

# Calcul de la pression dans le volume
# IN: centers:Density
# IN: centers:Temperature
# IN: cv: refState (r a partir de cv)
# OUT: centers:Pressure
# P = ro r T
def extractPressure(t):
    """Extract Pressure."""
    tp = Internal.copyRef(t)
    _extractPressure(tp)
    return tp

def _extractPressure(t):
    """Extract Pressure."""
    P._computeVariables2(t, ['centers:Pressure'])
    return None

# Calcul module de la vitesse 
# IN: centers:VelocityX, centers:VelocityY, centers:VelocityZ
# OUT: centers:VelocityMagnitude
# v = sqrt(vx**2+vy**2+vz**2)
def extractVelocityMagnitude(t):
    """Extract velocity magnitude."""
    tp = Internal.copyRef(t)
    _extractVelocityMagnitude(tp)
    return tp

def _extractVelocityMagnitude(t):
    """Extract velocity magnitude."""
    P._computeVariables2(t, ['centers:VelocityMagnitude'])
    return None

# Mach volumique
# IN: centers:Temperature
# IN: centers:VelocityX, centers:VelocityY, centers:VelocityZ
# IN: cv: refState (r a partir de cv)
# OUT: centers:Mach
# M = u/sqrt(gamma p/ro) - p = ro r T
def extractMach(t):
    """Extract Mach."""
    tp = Internal.copyRef(t)
    _extractMach(tp)
    return tp

def _extractMach(t):
    """Extract Mach."""
    P._computeVariables2(t, ['centers:Mach'])
    return None

# Mu volumique
# mu = sutherland(T)
# IN: centers:Temperature
# IN: Cs, Mus, Ts from refstate
# OUT: centers:ViscosityMolecular
def extractViscosityMolecular(t):
    """Extract Viscosity molecular."""
    tp = Internal.copyRef(t)
    _extractViscosityMolecular(tp)
    return tp

def _extractViscosityMolecular(t):
    """Extract Viscosity molecular."""
    P._computeVariables2(t, ['centers:ViscosityMolecular'])
    return None

# mut from spalart
# IN: centers:TurbulentSANuTilde
# IN: centers:Density
# IN: centers:ViscosityMolecular
# OUT: centers:ViscosityEddy
# kappa = ro * nutilde / mu
# mut = ro * nutilde * kappa^3 / (kappa^3 + 7.1^3)
def extractViscosityEddy(t):
    """Extract eddy viscosity."""
    tp = Internal.copyRef(t)
    _extractViscosityEddy(tp)
    return tp

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
def extractMutSurMu(t):
    """Extract Mut over mu."""
    tp = Internal.copyRef(t)
    _extractMutSurMu(tp)
    return tp

def _extractMutSurMu(t):
    """Extract Mut over mu."""
    C._initVars(t, '{centers:MutSurMu} = {centers:ViscosityEddy}/{centers:ViscosityMolecular}')
    return None

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
    C._initVars(teff, '{centers:FxA}={centers:MomentumX}/{centers:sMagnitude}')
    C._initVars(teff, '{centers:FyA}={centers:MomentumY}/{centers:sMagnitude}')
    C._initVars(teff, '{centers:FzA}={centers:MomentumZ}/{centers:sMagnitude}')
    C._rmVars(teff, ['centers:sx','centers:sy','centers:sz','sMagnitude'])
    return None

# extract shear stress
# IN: centers:ViscosityMolecular (mu)
# IN: centers:gradxVelocityX, centers:gradxVelocityY, ...
# OUT: centers:ShearStressXX,YY,ZZ,XY,XZ,YZ (tenseur symetrique)
# tau = -2/3 mu div u I + 2 mu D
# Attention : dans teff, c'est ViscosityMolecular+ViscosityEddy qui est dans la variable ViscosityMolecular
# ce qui est correct ici car sur les parois, viscosityEddy = 0
def extractShearStress(teff):
    """Extract shearStress."""
    tp = Internal.copyRef(teff)
    _extractShearStress(tp)
    return tp

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

# Extract tau.n
# IN: centers:ShearStress
# taun = tau.n
def extractTaun(teff):
    """Extract tau.n."""
    tp = Internal.copyRef(teff)
    _extractTaun(tp)
    return tp

def _extractTaun(teff):
    """Extract tau.n."""
    G._getNormalMap(teff)
    C._normalize(teff, ['centers:sx', 'centers:sy', 'centers:sz'])
    C._initVars(teff, '{centers:taunx} = {centers:ShearStressXX}*{centers:sx}+{centers:ShearStressXY}*{centers:sy}+{centers:ShearStressXZ}*{centers:sz}')
    C._initVars(teff, '{centers:tauny} = {centers:ShearStressXY}*{centers:sx}+{centers:ShearStressYY}*{centers:sy}+{centers:ShearStressYZ}*{centers:sz}')
    C._initVars(teff, '{centers:taunz} = {centers:ShearStressXZ}*{centers:sx}+{centers:ShearStressYZ}*{centers:sy}+{centers:ShearStressZZ}*{centers:sz}')
    return None

# Extract p.n
# IN: centers:Pressure
# pn = p.n
def extractPn(teff):
    """Extract p.n."""
    tp = Internal.copyRef(teff)
    _extractPn(tp)
    return tp

def _extractPn(teff):
    """Extract p.n."""
    G._getNormalMap(teff)
    C._normalize(teff, ['centers:sx', 'centers:sy', 'centers:sz'])
    C._initVars(teff, '{centers:pnx} = {centers:Pressure}*{centers:sx}')
    C._initVars(teff, '{centers:pny} = {centers:Pressure}*{centers:sy}')
    C._initVars(teff, '{centers:pnz} = {centers:Pressure}*{centers:sz}')
    return None

# Extract Force
# IN: centers:Pressure
# IN: centers:ShearStress
# IN: si withPInf = None: F = -p.n + tau.n
#     sinon F = -(p-pinf).n + tau.n
def extractForce(teff, withPInf=None):
    """Extract forces."""
    tp = Internal.copyRef(teff)
    _extractForce(tp, withPInf=withPInf)
    return tp

def _extractForce(teff, withPInf=None):
    """Extract forces."""
    G._getNormalMap(teff)
    C._normalize(teff, ['centers:sx', 'centers:sy', 'centers:sz'])
    if withPInf is None:
        C._initVars(teff, '{centers:Fx} = {centers:ShearStressXX}*{centers:sx}+{centers:ShearStressXY}*{centers:sy}+{centers:ShearStressXZ}*{centers:sz}-{centers:Pressure}*{centers:sx}')
        C._initVars(teff, '{centers:Fy} = {centers:ShearStressXY}*{centers:sx}+{centers:ShearStressYY}*{centers:sy}+{centers:ShearStressYZ}*{centers:sz}-{centers:Pressure}*{centers:sy}')
        C._initVars(teff, '{centers:Fz} = {centers:ShearStressXZ}*{centers:sx}+{centers:ShearStressYZ}*{centers:sy}+{centers:ShearStressZZ}*{centers:sz}-{centers:Pressure}*{centers:sz}')
    else:
        C._initVars(teff, '{centers:Fx} = ({centers:Pressure}-%20.16g)*{centers:sx}'%withPInf)
        C._initVars(teff, '{centers:Fy} = ({centers:Pressure}-%20.16g)*{centers:sy}'%withPInf)
        C._initVars(teff, '{centers:Fz} = ({centers:Pressure}-%20.16g)*{centers:sz}'%withPInf)
        C._initVars(teff, '{centers:Fx} = {centers:ShearStressXX}*{centers:sx}+{centers:ShearStressXY}*{centers:sy}+{centers:ShearStressXZ}*{centers:sz}-{centers:Fx}')
        C._initVars(teff, '{centers:Fy} = {centers:ShearStressXY}*{centers:sx}+{centers:ShearStressYY}*{centers:sy}+{centers:ShearStressYZ}*{centers:sz}-{centers:Fy}')
        C._initVars(teff, '{centers:Fz} = {centers:ShearStressXZ}*{centers:sx}+{centers:ShearStressYZ}*{centers:sy}+{centers:ShearStressZZ}*{centers:sz}-{centers:Fz}')
    return None

# Extract tangential friction vector
# IN: centers:shearStressXX,...
# OUT: centers:frictionX, centers:frictionY, centers:frictionZ
# taut = tau.n - (n. tau.n) n
def extractFrictionVector(teff):
    """Extract tangential friction vector."""
    tp = Internal.copyRef(teff)
    _extractFrictionVector(tp)
    return tp

def _extractFrictionVector(teff):
    """Extract tangential friction vector."""
    G._getNormalMap(teff)
    C._normalize(teff, ['centers:sx', 'centers:sy', 'centers:sz'])
    C._initVars(teff, '{centers:frictionX}={centers:ShearStressXX}*{centers:sx}+{centers:ShearStressXY}*{centers:sy}+{centers:ShearStressXZ}*{centers:sz}')
    C._initVars(teff, '{centers:frictionY}={centers:ShearStressXY}*{centers:sx}+{centers:ShearStressYY}*{centers:sy}+{centers:ShearStressYZ}*{centers:sz}')
    C._initVars(teff, '{centers:frictionZ}={centers:ShearStressXZ}*{centers:sx}+{centers:ShearStressYZ}*{centers:sy}+{centers:ShearStressZZ}*{centers:sz}')
    C._initVars(teff, '{centers:scal} = {centers:frictionX}*{centers:sx}+{centers:frictionY}*{centers:sy}+{centers:frictionZ}*{centers:sz}')
    C._initVars(teff, '{centers:frictionX}={centers:frictionX}-{centers:scal}*{centers:sx}')
    C._initVars(teff, '{centers:frictionY}={centers:frictionY}-{centers:scal}*{centers:sy}')
    C._initVars(teff, '{centers:frictionZ}={centers:frictionZ}-{centers:scal}*{centers:sz}')
    C._rmVars(teff, ['centers:scal', 'centers:sx', 'centers:sy', 'centers:sz'])
    return None

# tauw
# IN: centers:frictionX, centers:frictionY, centers:frictionZ
# OUT: centers: frictionMagnitude
# tauw = ||taut||
def extractFrictionMagnitude(teff):
    """Extract friction magnitude."""
    tp = Internal.copyRef(teff)
    _extractFrictionMagnitude(tp)
    return tp

def _extractFrictionMagnitude(teff):
    """Extract friction magnitude."""
    _extractFrictionVector(teff)
    C._magnitude(teff, ['centers:frictionX', 'centers:frictionY', 'centers:frictionZ'])
    return None

# IN: centers:frictionMagnitude
# IN: centers:Density2 - density
# OUT: centers:utau
# utau = sqrt(tauw/ro)
def extractUTau(teff):
    """Extract utau."""
    tp = Internal.copyRef(teff)
    _extractUTau(tp)
    return tp

def _extractUTau(teff):
    """Extract utau."""
    C._initVars(teff, '{centers:utau}=sqrt({centers:frictionMagnitude}/{centers:Density2})')
    return None

# IN: centers:frictionMagnitude
# IN: norm : (1/2*roinf*uinf**2)
# OUT: centers:Cf
# cf = tauw / (1/2*roinf*uinf**2)
def extractCf(teff, pinf, norm):
    """Extract friction coefficient."""
    tp = Internal.copyref(teff)
    _extractCf(tp, pinf, norm)
    return tp

def _extractCf(teff, norm):
    """Extract friction coefficient."""
    C._initVars(teff, '{centers:Cf} = {centers:frictionMagnitude} / %20.16g'%norm)
    return None

# IN: centers:Pressure
# IN: pinf: infinite pressure
# IN: norm : (1/2*roinf*uinf**2)
# OUT: centers:Cp
# cp = (p-pinf) / (1/2*roinf*uinf**2)
def extractCp(teff, pinf, norm):
    """Extract pressure coefficient."""
    tp = Internal.copyref(teff)
    _extractCp(tp, pinf, norm)
    return tp

def _extractCp(teff, pinf, norm):
    """Extract pressure coefficient."""
    C._initVars(teff, '{centers:Cp} = ({centers:Pressure}-%20.16g) / %20.16g'%(pinf,norm))
    return None

# integration des pressions Integ( p.n ds )
# IN: centers:Pressure
def integPressure(teff):
    """Integ p.n.ds"""
    ret = Pmpi.integNorm(teff, 'centers:Pressure')
    return ret[0]

# integration coefficient de pression Integ( Cp.n ds )
# IN: centers:Cp
def integCp(teff):
    """Integ Cp.n.ds"""
    ret = Pmpi.integNorm(teff, 'centers:Cp')
    return ret[0]

# integration coefficient de pression Integ( CM x Cp.n ds )
# IN: centers:Cp
def integMomentCp(teff, center):
    """Integ CM x Cp.n.ds"""
    ret = Pmpi.integMomentNorm(teff, center, 'centers:Cp')
    return ret

# integration Integ ( tau.n ds )
# IN: centers:ShearStress
def integTaun(teff):
    """Integ tau.n.ds"""
    G._getNormalMap(teff)
    C._normalize(teff, ['centers:sx', 'centers:sy', 'centers:sz'])
    C._initVars(teff, '{centers:taunx} = {centers:ShearStressXX}*{centers:sx}+{centers:ShearStressXY}*{centers:sy}+{centers:ShearStressXZ}*{centers:sz}')
    C._initVars(teff, '{centers:tauny} = {centers:ShearStressXY}*{centers:sx}+{centers:ShearStressYY}*{centers:sy}+{centers:ShearStressYZ}*{centers:sz}')
    C._initVars(teff, '{centers:taunz} = {centers:ShearStressXZ}*{centers:sx}+{centers:ShearStressYZ}*{centers:sy}+{centers:ShearStressZZ}*{centers:sz}')
    retx = Pmpi.integ(teff, 'centers:taunx')
    rety = Pmpi.integ(teff, 'centers:tauny')
    retz = Pmpi.integ(teff, 'centers:taunz')
    C._rmVars(teff, ['centers:sx', 'centers:sy', 'centers:sz', 'centers:taunx', 'centers:tauny', 'centers:taunz'])
    return [retx[0],rety[0],retz[0]]

# integration Integ ( CM x tau.n ds )
# IN: centers:ShearStress
def integMomentTaun(teff, center):
    """Integ CM x tau.n.ds"""
    G._getNormalMap(teff)
    C._normalize(teff, ['centers:sx', 'centers:sy', 'centers:sz'])
    C._initVars(teff, '{centers:taunx} = {centers:ShearStressXX}*{centers:sx}+{centers:ShearStressXY}*{centers:sy}+{centers:ShearStressXZ}*{centers:sz}')
    C._initVars(teff, '{centers:tauny} = {centers:ShearStressXY}*{centers:sx}+{centers:ShearStressYY}*{centers:sy}+{centers:ShearStressYZ}*{centers:sz}')
    C._initVars(teff, '{centers:taunz} = {centers:ShearStressXZ}*{centers:sx}+{centers:ShearStressYZ}*{centers:sy}+{centers:ShearStressZZ}*{centers:sz}')
    ret = Pmpi.integMoment(teff, center, ['centers:taunx', 'centers:tauny', 'centers:taunz'])
    C._rmVars(teff, ['centers:sx', 'centers:sy', 'centers:sz', 'centers:taunx', 'centers:tauny', 'centers:taunz'])
    return ret


# Integration des Cf Integ( Cf. ds )
# IN: centers:Cf
def integCf(teff):
    """Integ Cf.ds"""
    ret = Pmpi.integ(teff, 'centers:Cf')
    return ret[0]

# Attention : deja -pinf dans les Fx,Fy,Fz
# et deja division par 0.5*roinf*uinf^2
# retourne donc des charges adimensionnees
def integLoads(teff):
    """Integ F.ds"""
    retx = Pmpi.integ(teff, 'centers:FxA')
    rety = Pmpi.integ(teff, 'centers:FyA')
    retz = Pmpi.integ(teff, 'centers:FzA')
    return [retx[0],rety[0],retz[0]]

#=============================================================
# Profile extractions 
# input : effort tree (teff) + volume tree t
#=============================================================

# Extrait des profiles proche paroi en i,j,k (deux grandeurs constantes)
# IN: zonePath: chemin de la zone a extraire
# IN: i,j,k de la ligne constante
# OUT: zp: ligne avec yw (distance a la paroi suivant la ligne)
def extractProfile(t, zonePath, i=-1, j=-1, k=-1):
    """Extract given profile."""
    z = Internal.getNodeFromPath(t, zonePath)
    if z is None: return None
    #zp = T.subzone(z, (i,j,3), (i,j,-1))
    if k == -1: # suivant k
        zp = T.subzone(z, (i,j,3), (i,j,-3)) # suppress Ghostcells
    elif j == -1:
        zp = T.subzone(z, (i,3,k), (i,-3,k))
    else:
        zp = T.subzone(z, (3,j,k), (-3,j,k))

    Internal._rmNodesFromName(zp, 'ZoneBC')
    zp = C.node2Center(zp)
    L = D.getLength(zp)
    D._getCurvilinearAbscissa(zp)
    C._initVars(zp, '{yw}={s}*%20.16g'%L)
    return zp

# Calcul la norme de la vitesse tangentielle
# IN: zp: ligne with VelocityX, VelocityY, VelocityZ
# OUT: VelocityTangential: norm of tangential velocity
def extractVelocityTangential(zp, teff):
    """Extract tangential speed."""
    zpp = Internal.copyRef(zp)
    _extractVelocityTangential(zpp, teff)
    return zpp

def _extractVelocityTangential(zp, teff):
    """Extract tangential speed."""
    # calcul la normale sur teff
    G._getNormalMap(teff)
    teff = C.center2Node(teff, ['centers:sx','centers:sy','centers:sz'])
    C._normalize(teff, ['sx', 'sy','sz'])
    C._rmVars(teff, ['centers:sx', 'centers:sy', 'centers:sz'])
    # Extract point wall de zp
    Pwall = T.subzone(zp, (1,1,1), (1,1,1))
    # localise Pwall sur teff
    zeff = Internal.getZones(teff)[0]
    hook = C.createHook(zeff, function='nodes')
    nodes, dist = C.nearestNodes(hook, Pwall)
    n = nodes[0]
    sx = C.getValue(zeff, 'sx', n) # must be normalized
    sy = C.getValue(zeff, 'sy', n)
    sz = C.getValue(zeff, 'sz', n)
    print('INFO: extractVelocityTangential: found wall point sx=%g, sy=%g, sz=%g, index=%d\n'%(sx, sy, sz, n))
    C._initVars(zp, '{VelocityTangential1}={VelocityX} * %20.16g+{VelocityY} * %20.16g +{VelocityZ} * %20.16g'%(sx,sy,sz))
    C._initVars(zp, '{VelocityTangential}=sqrt(({VelocityX}-{VelocityTangential1}*%20.16g)**2+({VelocityY}-{VelocityTangential1}*%20.16g)**2+({VelocityZ}-{VelocityTangential1}*%20.16g)**2)'%(sx,sy,sz))
    C._rmVars(zp, ['VelocityTangential1','sx','sy','sz'])
    return None

# Calcul y+ et u+ sur le profil
# on suppose que sur la ligne i=1 est la paroi (sinon reorder)
# IN: zp: ligne de maillage with VelocityTangential, yw
# IN: teff: with utau, ViscosityMolecular computed
def extractYPlus(zp, teff):
    """Extract y+ and u+."""
    zpp = Internal.copyRef(zp)
    _extractYPlus(zpp, teff)
    return zpp

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
    print('INFO: extractYPlus: found wall point ro=%g, utau=%g, mu=%g, index=%d\n'%(ro, utau, mu, n))
    C._initVars(zp, '{yPlus}={yw} * %20.16g * %20.16g / %20.16g'%(ro,utau,mu))
    C._initVars(zp, '{yPlus}=log10({yPlus})')
    C._initVars(zp, '{uPlus}={VelocityTangential}/%20.16g'%utau)
    return None

# Calcul delta (hauteur couche limite), teta (moment thickness), Rteta
# on suppose que sur la ligne i=1 est la paroi (sinon reorder)
# IN: zp: ligne de maillage with VelocityTangential
# IN: Uinf,Roinf,MuInf etat infini ou au moins tout en haut de la ligne
# OUT: delta,theta,Rtheta
def extractRtheta(zp, Uinf, Roinf, Muinf):
    """Compute delta, teta and Rteta."""
    # delta
    j = 1
    while (C.getValue(zp, 'VelocityTangential', (j,1,1)) < 0.99*Uinf): j += 1
    dx = C.getValue(zp, 'CoordinateX', (j,1,1)) - C.getValue(zp, 'CoordinateX', (1,1,1))
    dy = C.getValue(zp, 'CoordinateY', (j,1,1)) - C.getValue(zp, 'CoordinateY', (1,1,1))
    dz = C.getValue(zp, 'CoordinateZ', (j,1,1)) - C.getValue(zp, 'CoordinateZ', (1,1,1))
    delta = numpy.sqrt(dx*dx+dy*dy+dz*dz)

    # subzone la ligne pour faire l'integrale de teta
    zpp = T.subzone(zp, (1,1,1),(j+1,-1,-1))
    ret = P.integ(zpp, 'dTheta')
    theta = ret[0]
    Rtheta = Uinf * theta * Roinf / Muinf

    C._initVars(zpp, '{dDelta} = ( 1. - {VelocityTangential} / %20.16g ) * ( {Density} / %20.16g)'%(Uinf,Roinf) )
    ret2 = P.integ(zpp, 'dDelta')
    delta = ret2[0]
    ShapeFactor = delta / theta    

    return delta,theta,Rtheta,ShapeFactor
