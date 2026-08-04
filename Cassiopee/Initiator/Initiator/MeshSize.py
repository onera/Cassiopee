# Return first mesh size for y+=1
import math

# Mesh size as pointwise (attached turbulent boundary layer)
def meshSize1(UInf, RoInf, MuInf, LInf, yplus=1.):
    """Return first cell height to match y+."""
    ReInf = RoInf*UInf*LInf/MuInf
    Cf = 0.026/ReInf**(1./7.)
    tauw = Cf*RoInf*(UInf**2)/2.
    utau = math.sqrt(tauw/RoInf)
    hp = yplus*MuInf/(utau*RoInf)
    print('INFO: tauw: %20.16g'%tauw)
    print('INFO: utau: %20.16g'%utau)
    return hp

# From the Reynolds number
def meshSize2(UInf, RoInf, ReInf, LInf, yplus=1.):
    """Return first cell height to match y+."""
    MuInf = RoInf*UInf*LInf/ReInf
    return meshSize1(UInf, RoInf, MuInf, LInf, yplus)

# Marco - attached turbulent boundary layer
# with correction of e/c = thickness to chord ratio
def meshSize3(UInf, RoInf, ReInf, LInf, esurc=0.012, yplus=1.):
    """Return first cell height to match y+."""
    xsurL = 0.5
    MuInf = RoInf*UInf*LInf/ReInf
    correction = math.exp(4.5*pow(esurc, 1.3))
    Cf = 0.058*math.pow(ReInf*xsurL, -0.2) * correction
    utau = math.sqrt(RoInf*Cf*UInf*UInf*0.5)
    hp = MuInf * yplus / RoInf / utau
    print('INFO: utau: %20.16g'%utau)
    return hp

# Marco - laminar boundary layer
def meshSize4(UInf, RoInf, ReInf, LInf, esurc=0.012, yplus=1.):
    """Return first cell height to match y+."""
    xsurL = 0.5
    MuInf = RoInf*UInf*LInf/ReInf
    correction = math.exp(4.5*pow(esurc, 1.3))
    Cf = 0.664/math.sqrt(ReInf*xsurL)* correction
    utau = math.sqrt(RoInf*Cf*UInf*UInf*0.5)
    hp = MuInf * yplus / RoInf / utau
    print('INFO: utau: %20.16g'%utau)
    return hp

def meshSize(UInf, RoInf, ReInf, LInf, esurc=0.012, yplus=1., algo='Turbulent'):
    """Return the height of first wall cell to match a certain y+."""
    if algo == 'Turbulent':
        return meshSize2(UInf, RoInf, ReInf, LInf, yplus)
    elif algo == 'TurbulentCorr':
        return meshSize3(UInf, RoInf, ReInf, LInf, esurc, yplus)
    elif algo == 'LaminarCorr':
        return meshSize4(UInf, RoInf, ReInf, LInf, esurc, yplus)
    else:
        raise ValueError('meshSize: unknown algo.')

def boundaryLayerHeight(ReInf, algo='Turbulent'):
    """Return adimensioned boundary layer height."""
    if algo == 'Laminar':
        delta = 0.75*5.0*ReInf**(-0.5)
    elif algo == 'Turbulent':
        delta = 0.75*0.37*ReInf**(-1./5.)
    return delta

if __name__ == '__main__':
    print(meshSize1(1., 1.225, 0.000017894, 1.0))
    print(meshSize2(1., 1.225, 68458., 1.0))
    print(meshSize3(1., 1.225, 68458., 1.0))
    print(meshSize4(1., 1.225, 68458., 1.0))
