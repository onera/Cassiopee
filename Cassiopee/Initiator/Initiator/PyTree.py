"""Initialization of grid solutions.
"""
from . import Initiator
__version__ = Initiator.__version__

from .MeshSize import meshSize
from .Adim import adim1, adim2, adim3, dim1, dim2, dim3, dim4

try:
    import Converter
    import Converter.PyTree as C
    import Converter.Internal as Internal
except ImportError:
    raise ImportError("Initiator.PyTree: requires Converter.PyTree module.")

def _applyGaussianAL(t, listOfLoads, listOfALPositions, listOfRotMat,
                     localEpsX, localEpsY, localEpsZ, NbBlades, NbPointsAL,
                     TruncVarLoads, TruncVarVelos):
    writeDim = False
    fa = C.getFields(Internal.__FlowSolutionCenters__, t, api=3)
    if fa != []:
        Initiator._applyGaussianAL(fa, listOfLoads, listOfALPositions, listOfRotMat,
                                   localEpsX, localEpsY, localEpsZ, NbBlades, NbPointsAL, TruncVarLoads, TruncVarVelos)
    return None

def initConst(t, adim='adim1', MInf=None, alphaZ=0., alphaY=0., ReInf=1.e8,
              loc='nodes'):
    """Init the pyTree by the reference state if it is defined in t, else by input parameters.
    Usage: initConst(t, adim, MInf, alphaZ, alphaY, ReInf, loc)"""
    tp = Internal.copyRef(t)
    _initConst(tp, adim, MInf, alphaZ, alphaY, ReInf, loc)
    return tp

def _initConst(t, adim='adim1', MInf=None, alphaZ=0., alphaY=0., ReInf=1.e8,
               loc='nodes'):
    """Init the pyTree by the reference state if it is defined in t, else by input parameters."""
    if MInf is None: # recuperation de reference state
        eq = Internal.getNodeFromName(t, 'GoverningEquations')
        state = Internal.getNodeFromName(t, 'ReferenceState')
        if state is None: raise ValueError("initConst: no reference state and no argument.")
        vars0 = ['Density', 'MomentumX', 'MomentumY', 'MomentumZ', 'EnergyStagnationDensity']
        if eq is not None and Internal.getValue(eq) == 'NSTurbulent':
            vars0 += ['TurbulentSANuTildeDensity', 'TurbulentEnergyKineticDensity', 'TurbulentDissipationDensity']
        for v in vars0:
            node = Internal.getNodeFromName(state, v)
            if node is not None:
                val = Internal.getValue(node)
                C._initVars(t, loc+':'+v, val)
    else: # recuperation des arguments
        nodes = Internal.getZones(t)
        for z in nodes:
            a = C.getFields(Internal.__GridCoordinates__, z, api=3)[0]
            if loc == 'nodes':
                a = Initiator.initConst(a, adim, MInf, alphaZ, alphaY, ReInf)
                z = C.setFields([a], z, 'nodes')
            elif loc == 'centers':
                a = Converter.node2Center(a)
                a = Initiator.initConst(a, adim, MInf, alphaZ, alphaY, ReInf)
                a = Converter.rmVars(a, ['CoordinateX', 'CoordinateY', 'CoordinateZ'])
                z = C.setFields([a], z, 'centers')
            else:
                raise ValueError("initConst: wrong location: %s."%loc)
    return None

# Initialise a partir d'un ReferenceState les grandeurs conservatives
# La vitesse est proportionelle a la distance
def _initDist(t, adim='adim1', loc='nodes'):
    eq = Internal.getNodeFromName(t, 'GoverningEquations')
    state = Internal.getNodeFromName(t, 'ReferenceState')
    if state is None: raise ValueError("initConst: no reference state and no argument.")
    vars0 = ['Density', 'MomentumX', 'MomentumY', 'MomentumZ',
             'EnergyStagnationDensity']
    if eq is not None and Internal.getValue(eq) == 'NSTurbulent':
        vars0 += ['TurbulentSANuTildeDensity', 'TurbulentEnergyKineticDensity', 'TurbulentDissipationDensity']
    for v in vars0:
        node = Internal.getNodeFromName(state, v)
        if node is not None:
            val = Internal.getValue(node)
            C._initVars(t, loc+':'+v, val)
    maxd = C.getMaxValue(t, '%s:TurbulentDistance'%loc)
    mind = C.getMinValue(t, '%s:TurbulentDistance'%loc)

    C._initVars(t, '{%s:MomentumX}={%s:MomentumX}*{%s:TurbulentDistance}/%g'%(loc,loc,loc,maxd))
    C._initVars(t, '{%s:MomentumY}={%s:MomentumY}*{%s:TurbulentDistance}/%g'%(loc,loc,loc,maxd))
    C._initVars(t, '{%s:MomentumZ}={%s:MomentumZ}*{%s:TurbulentDistance}/%g'%(loc,loc,loc,maxd))
    return None

def initLamb(t, position=(0.,0.), Gamma=2., MInf=0.5, loc='nodes'):
    """Init the pyTree with a Lamb vortex of intensity Gamma and position (x0,y0).
    Usage: initLamb(t, (x0,y0), Gamma, MInf)"""
    tp = Internal.copyRef(t)
    _initLamb(tp, position, Gamma, MInf, loc)
    return tp

def _initLamb(t, position=(0.,0.), Gamma=2., MInf=0.5, loc='nodes'):
    """Init the pyTree with a Lamb vortex of intensity Gamma and position (x0,y0)."""
    nodes = Internal.getZones(t)
    for z in nodes:
        coordn = C.getFields(Internal.__GridCoordinates__, z, api=3)
        if coordn == []:
            print('Warning: initLamb: zone '+z[0]+' has no coordinates. Skipped...')
            continue
        coordn = coordn[0]
        if loc == 'nodes':
            a = C.getFields(Internal.__FlowSolutionNodes__, z, api=3)[0]
            if a == []: a = coordn
            else: Converter._addVars([a, coordn])
            a = Initiator.initLamb(a, position, Gamma, MInf)
            z = C.setFields([a], z, 'nodes')
        else:
            coordc = Converter.node2Center(coordn)
            ac = C.getFields(Internal.__FlowSolutionCenters__, z, api=3)[0]
            if ac == []: ac = coordc
            else: Converter._addVars([ac, coordc])
            ac = Initiator.initLamb(ac, position, Gamma, MInf)
            ac = Converter.rmVars(ac, ['x', 'y', 'z'])
            z = C.setFields([ac], z, 'centers')
    return None

def initWissocq(t, position=(0.5,0.5), Gamma=0.07, MInf=0.5, loc='nodes'):
    """Init the pyTree with Wissocq's vortex of
    intensity Gamma and position (x0,y0).
    Usage: initWissocq(t, (x0,y0), Gamma, MInf)"""
    tp = Internal.copyRef(t)
    _initWissocq(tp, position, Gamma, MInf, loc)
    return tp

def _initWissocq(t, position=(0.5,0.5), Gamma=0.07, MInf=0.5, loc='nodes'):
    nodes = Internal.getZones(t)
    for z in nodes:
        coordn = C.getFields(Internal.__GridCoordinates__, z, api=3)
        if coordn == []:
            print ('Warning: initWissocq: zone '+z[0]+' has no coordinates. Skipped...')
            continue
        coordn = coordn[0]
        if loc == 'nodes':
            a = C.getFields(Internal.__FlowSolutionNodes__, z, api=3)[0]
            if a == []: a = coordn
            else: Converter._addVars([a, coordn])
            a = Initiator.initWissocq(a, position, Gamma, MInf)
            z = C.setFields([a], z, 'nodes')
        else:
            coordc = Converter.node2Center(coordn)
            ac = C.getFields(Internal.__FlowSolutionCenters__, z, api=3)[0]
            if ac == []: ac = coordc
            else: Converter._addVars([ac, coordc])
            ac = Initiator.initWissocq(ac, position, Gamma, MInf)
            ac = Converter.rmVars(ac, ['x', 'y', 'z'])
            z = C.setFields([ac], z, 'centers')
    return None


def initVisbal(t, position=(0.,0.), Gamma=2., MInf=0.5, loc='nodes'):
    """Init the array defining a grid with a Visbal vortex of intensity Gamma and position (x0,y0).
    Returns the array of the grid + cfd field in centers
    Usage: initVisbal(array, (x0,y0), Gamma, MInf)"""
    tp = Internal.copyRef(t)
    _initVisbal(tp, position, Gamma, MInf, loc)
    return tp

def _initVisbal(t, position=(0.,0.), Gamma=2., MInf=0.5, loc='nodes'):
    """Init the array defining a grid with a Visbal vortex of intensity Gamma and position (x0,y0)."""
    nodes = Internal.getZones(t)
    for z in nodes:
        coordn = C.getFields(Internal.__GridCoordinates__, z, api=3)
        if coordn == []:
            print ('Warning: initVisbal: zone '+z[0]+' has no coordinates. Skipped...')
            continue
        coordn = coordn[0]
        if loc == 'nodes':
            a = C.getFields(Internal.__FlowSolutionNodes__, z, api=3)[0]
            if a == []: a = coordn
            else: Converter._addVars([a, coordn])
            a = Initiator.initVisbal(a, position, Gamma, MInf)
            z = C.setFields([a], z, 'nodes')
        else:
            coordc = Converter.node2Center(coordn)
            ac = C.getFields(Internal.__FlowSolutionCenters__, z, api=3)[0]
            if ac == []: ac = coordc
            else: Converter._addVars([ac, coordc])
            ac = Initiator.initVisbal(ac, position, Gamma, MInf)
            ac = Converter.rmVars(ac, ['x', 'y', 'z'])
            z = C.setFields([ac], z, 'centers')
    return None

def initYee(t, position=(0.,0.), Gamma=2., MInf=0.5, loc='nodes'):
    """Init the array defining a grid with a Yee vortex of intensity Gamma and position (x0,y0).
    Returns the array of the grid + cfd field in centers
    Usage: initYee(array, (x0,y0), Gamma, MInf)"""
    tp = Internal.copyRef(t)
    _initYee(tp, position, Gamma, MInf, loc)
    return tp

def _initYee(t, position=(0.,0.), Gamma=2., MInf=0.5, loc='nodes'):
    """Init the array defining a grid with a Yee vortex of intensity Gamma and position (x0,y0)."""
    nodes = Internal.getZones(t)
    for z in nodes:
        coordn = C.getFields(Internal.__GridCoordinates__, z, api=3)
        if coordn == []:
            print ('Warning: initYee: zone '+z[0]+' has no coordinates. Skipped...')
            continue
        coordn = coordn[0]
        if loc == 'nodes':
            a = C.getFields(Internal.__FlowSolutionNodes__, z, api=3)[0]
            if a == []: a = coordn
            else: Converter._addVars([a, coordn])
            a = Initiator.initYee(a, position, Gamma, MInf)
            z = C.setFields([a], z, 'nodes')
        else:
            coordc = Converter.node2Center(coordn)
            ac = C.getFields(Internal.__FlowSolutionCenters__, z, api=3)[0]
            if ac == []: ac = coordc
            else: Converter._addVars([ac, coordc])
            ac = Initiator.initYee(ac, position, Gamma, MInf)
            ac = Converter.rmVars(ac, ['x', 'y', 'z'])
            z = C.setFields([ac], z, 'centers')
    return None

def initScully(t, position=(0.,0.), Gamma=2.,
               coreRadius=1., MInf=0.5, model=0, loc='nodes'):
    """Init the array defining a block field with a Scully vortex
    of intensity Gamma, core radius coreRadius and position (x0,y0).
    Usage: initScully(array, (x0,y0), Gamma, coreRadius, MInf, model)"""
    tp = Internal.copyRef(t)
    _initScully(tp, position, Gamma, coreRadius, MInf, model, loc)
    return tp

def _initScully(t, position=(0.,0.), Gamma=2.,
                coreRadius=1., MInf=0.5, model=0, loc='nodes'):
    """Init the array defining a block field with a Scully vortex
    of intensity Gamma, core radius coreRadius and position (x0,y0)."""
    nodes = Internal.getZones(t)
    for z in nodes:
        coordn = C.getFields(Internal.__GridCoordinates__, z, api=3)
        if coordn == []:
            print ('Warning: initScully: zone '+z[0]+' has no coordinates. Skipped...')
            continue
        coordn = coordn[0]
        if loc == 'nodes':
            a = C.getFields(Internal.__FlowSolutionNodes__, z, api=3)[0]
            if a == []: a = coordn
            else: Converter._addVars([a, coordn])
            a = Initiator.initScully(a, position, Gamma, coreRadius, MInf,
                                     model)
            C.setFields([a], z, 'nodes')
        else:
            coordc = Converter.node2Center(coordn)
            ac = C.getFields(Internal.__FlowSolutionCenters__, z, api=3)[0]
            if ac == []: ac = coordc
            else: Converter._addVars([ac, coordc])
            ac = Initiator.initScully(ac, position, Gamma, coreRadius, MInf,
                                      model)
            ac = Converter.rmVars(ac, ['x', 'y', 'z'])
            C.setFields([ac], z, 'centers')
    return None

def overlayField(t1, t2, MInf=0.5, loc='nodes'):
    """Overlay the field of zone1 and zone2 in a unique zone.
    Usage: overlayField(z1, z2, MInf, loc)"""
    tp = Internal.copyRef(t1)
    _overlayField(tp, t2, MInf, loc)
    return tp

def _overlayField(t1, t2, MInf=0.5, loc='nodes'):
    """Overlay the field of zone1 and zone2 in a unique zone."""
    nodes = Internal.getZones(t1)
    nodes2 = Internal.getZones(t2)
    for c, z1 in enumerate(nodes):
        if loc == 'centers':
            a1 = C.getAllFields(z1, 'centers', api=1)[0]
            x1 = C.getFields(Internal.__GridCoordinates__, z1, api=1)[0]
            x1 = Converter.node2Center(x1)
            a1 = Converter.addVars([x1, a1])
            z2 = nodes2[c]
            a2 = C.getAllFields(z2, 'centers', api=1)[0]
            x2 = C.getFields(Internal.__GridCoordinates__, z2, api=1)[0]
            x2 = Converter.node2Center(x2)
            a2 = Converter.addVars([x2, a2])
            ret = Initiator.overlayField(a1, a2, MInf)
            C.setFields([ret], z1, 'centers')
        else:
            a1 = C.getAllFields(z1, 'nodes', api=1)[0]
            z2 = nodes2[c]
            a2 = C.getAllFields(z2, 'nodes', api=1)[0]
            ret = Initiator.overlayField(a1, a2, MInf)
            C.setFields([ret], z1, 'nodes')
    return None

# passage variables conservatives en variables primitives (ro,u,T)
def cons2Prim(t, Gamma=1.4, Rgas=287.053, loc='centers'):
    """Compute primitive variables from conservative variables"""
    tp = Internal.copyRef(t)
    _cons2Prim(t, Gamma=Gamma, Rgas=Rgas, loc=loc)
    return tp

def _cons2Prim(t, Gamma=1.4, Rgas=287.053, loc='centers'):
    """Compute primitive variables from conservative variables"""
    if loc == 'centers':
        C._initVars(t, '{centers:VelocityX} = {centers:MomentumX}/{centers:Density}')
        C._rmVars(t, 'centers:MomentumX')
        C._initVars(t, '{centers:VelocityY} = {centers:MomentumY}/{centers:Density}')
        C._rmVars(t, 'centers:MomentumY')
        C._initVars(t, '{centers:VelocityZ} = {centers:MomentumZ}/{centers:Density}')
        C._rmVars(t, 'centers:MomentumZ')
        K = (Gamma - 1.)/Rgas
        C._initVars(t, '{centers:Temperature} = ({centers:EnergyStagnationDensity}/{centers:Density} - 0.5*({centers:VelocityX}**2+{centers:VelocityY}**2+{centers:VelocityZ}**2))*%20.16g'%(K))
        C._rmVars(t, 'centers:EnergyStagnationDensity')
    else:
        C._initVars(t, '{VelocityX} = {MomentumX}/{Density}')
        C._rmVars(t, 'MomentumX')
        C._initVars(t, '{VelocityY} = {MomentumY}/{Density}')
        C._rmVars(t, 'MomentumY')
        C._initVars(t, '{VelocityZ} = {MomentumZ}/{Density}')
        C._rmVars(t, 'MomentumZ')
        K = (Gamma - 1.)/Rgas
        C._initVars(t, '{Temperature} = ({EnergyStagnationDensity}/{Density} - 0.5*({VelocityX}**2+{VelocityY}**2+{VelocityZ}**2))*%20.16g'%(K))
        C._rmVars(t, 'centers:EnergyStagnationDensity')
    return None
