"""Initialization of grid solutions.
"""
__version__ = '4.0'
__author__ = "Stephanie Peron, Christophe Benoit, Gaelle Jeanfaivre, Pascal Raud"

from . import initiator

from .MeshSize import meshSize
from .Adim import adim1, adim2, adim3, dim1, dim2, dim3, dim4

def _applyGaussianAL(a, listOfLoads, listOfALPositions, listOfRotMat,
                    localEpsX, localEpsY, localEpsZ, NbBlades, NbPointsAL,
                    TruncVarLoads, TruncVarVelos):
    if isinstance(a[0], list):
        initiator.applyGaussianAL(a,
                listOfLoads,listOfALPositions,listOfRotMat,
                localEpsX, localEpsY, localEpsZ, NbBlades, NbPointsAL, TruncVarLoads, TruncVarVelos)
    else:
        initiator.applyGaussianAL(
            [a],listOfLoads,listOfALPositions,listOfRotMat,
            localEpsX, localEpsY, localEpsZ, NbBlades, NbPointsAL, TruncVarLoads, TruncVarVelos)
    return None

def initConst(a, adim='adim1', MInf=None, alphaZ=0., alphaY=0., ReInf=1.e8):
    """Init a by a constant field.
    Usage: initConst(a, MInf, alphaZ, alphaY, ReInf)"""
    try:
        import Converter as C
        import KCore.Adim as Adim
    except ImportError:
        raise ImportError("initConst: requires Converter module.")

    if MInf is None: raise ValueError("initConst: MInf must be defined.")
    if adim == 'adim1': cons = Adim.adim1(MInf, alphaZ, alphaY, ReInf)
    else: cons = Adim.adim2(MInf, alphaZ, alphaY, ReInf)

    a = C.initVars(a, 'ro', cons[0])
    a = C.initVars(a, 'rou', cons[1])
    a = C.initVars(a, 'rov', cons[2])
    a = C.initVars(a, 'row', cons[3])
    a = C.initVars(a, 'roE', cons[4])
    return a

def initLamb(a, position=(0.,0.), Gamma=2., MInf=0.5):
    """Init a with a Lamb vortex of
    intensity Gamma and position (x0,y0).
    Usage: initLamb(a, (x0,y0), Gamma, MInf)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(initiator.initLamb(i, position, Gamma, MInf))
        return b
    else:
        return initiator.initLamb(a, position, Gamma, MInf)

def initVisbal(a, position=(0.,0.), Gamma=2., MInf=0.5):
    """Init a with a Visbal vortex of intensity Gamma and position (x0,y0).
    Usage: initVisbal(array, (x0,y0), Gamma, MInf)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(initiator.initVisbal(i, position, Gamma, MInf))
        return b
    else:
        return initiator.initVisbal(a, position, Gamma, MInf)

def initYee(a, position=(0.,0.), Gamma=2., MInf=0.5):
    """Init a with a Yee vortex of intensity Gamma and position (x0,y0).
    Usage: initYee(a, (x0,y0), Gamma, Minf)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(initiator.initYee(i, position, Gamma, MInf))
        return b
    else:
        return initiator.initYee(a, position, Gamma, MInf)

def initScully(a, position=(0.,0.), Gamma=2.,
               coreRadius=1., MInf=0.5, model=0):
    """Init a with a Scully vortex
    of intensity Gamma, core radius coreRadius and position (x0,y0).
    Usage:
    initScully(a, (x0,y0), Gamma, coreRadius, MInf, model)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(initiator.initScully(i, position, Gamma, coreRadius,
                                          MInf, model))
        return b
    else:
        return initiator.initScully(a, position, Gamma,
                                    coreRadius, MInf, model)

def overlayField(a1, a2, MInf=0.5):
    """Overlay the field of a1 and a2.
    Usage: overlayField(a1, a2, MInf)"""
    if isinstance(a1[0], list):
        b = []; c = 0
        for i in a1:
            b.append(initiator.overlayField(i, a2[c], MInf))
            c += 1
        return b
    else:
        return initiator.overlayField(a1, a2, MInf)

def initWissocq(a, position=(0.5,0.5), Gamma=0.07, MInf=0.5):
    """Init a with Wissocq's vortex of
    intensity Gamma and position (x0,y0).
    Usage: initWissocq(a, (x0,y0), Gamma, MInf)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(initiator.initWissocq(i, position, Gamma, MInf))
        return b
    else:
        return initiator.initWissocq(a, position, Gamma, MInf)

def cons2Prim(a, Gamma=1.4, Rgas=287.053):
    """Compute primitive variables from conservative variables"""
    a = C.initVars(t, '{VelocityX} = {MomentumX}/{Density}')
    a = C.initVars(t, '{VelocityY} = {MomentumY}/{Density}')
    a = C.initVars(t, '{VelocityZ} = {MomentumZ}/{Density}')
    K = (Gamma - 1.)/Rgas
    a = C.initVars(t, '{Temperature} = ({EnergyStagnationDensity}/{Density} - 0.5*({VelocityX}**2+{VelocityY}**2+{VelocityZ}**2))*%20.16g'%(K))
    return a
