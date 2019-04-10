# 
# Python Interface to initialize pyTrees solutions
#
from . import Initiator
__version__ = Initiator.__version__

try: range = xrange
except: pass

try:
    import Converter
    import Converter.PyTree as C
    import Converter.Internal as Internal
except:
    raise ImportError("Initiator.PyTree: requires Converter.PyTree module.")

def initConst(t, adim='adim1', MInf=None, alphaZ=0., alphaY=0., ReInf=1.e8, 
              loc='nodes'):
    """Init the pyTree by the reference state if it is defined in t, else by
    input parameters.
    Usage: initConst(t, adim, MInf, alphaZ, alphaY, ReInf, loc)"""
    tp = Internal.copyRef(t)
    _initConst(tp, adim, MInf, alphaZ, alphaY, ReInf, loc)
    return tp

def _initConst(t, adim='adim1', MInf=None, alphaZ=0., alphaY=0., ReInf=1.e8, 
               loc='nodes'):
    if MInf is None: # recuperation de reference state
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
    else: # recuperation des arguments
        nodes = Internal.getZones(t)
        for z in nodes:
            a = C.getFields(Internal.__GridCoordinates__, z)[0]
            if loc == 'nodes':
                a = Initiator.initConst(a, adim, MInf, alphaZ, alphaY, ReInf)
                z = C.setFields([a], z, 'nodes')
            elif loc == 'centers':
                a = Converter.node2Center(a)
                a = Initiator.initConst(a, adim, MInf, alphaZ, alphaY, ReInf)
                a = Converter.rmVars(a,
                                     ['CoordinateX', 'CoordinateY', 'CoordinateZ'])
                z = C.setFields([a], z, 'centers')
            else:
                raise ValueError("initConst: wrong location: %s."%loc)
    return None

def initLamb(t, position=(0.,0.), Gamma=2., MInf=0.5, loc='nodes'):
    """Init the pyTree with a Lamb vortex of
    intensity Gamma and position (x0,y0).
    Usage: initLamb(t, (x0,y0), Gamma, MInf)"""
    tp = Internal.copyRef(t)
    _initLamb(tp, position, Gamma, MInf, loc)
    return tp

def _initLamb(t, position=(0.,0.), Gamma=2., MInf=0.5, loc='nodes'):
    nodes = Internal.getZones(t)
    for z in nodes:
        coordn = C.getFields(Internal.__GridCoordinates__, z)
        if coordn == []:
            print ('Warning: initLamb: zone '+z[0]+' has no coordinates. Skipped...')
            continue
        coordn = coordn[0]
        if loc == 'nodes':
            a = C.getFields(Internal.__FlowSolutionNodes__, z)[0]
            if a == []: a = coordn
            else: Converter._addVars([a, coordn])
            a = Initiator.initLamb(a, position, Gamma, MInf)
            z = C.setFields([a], z, 'nodes')
        else:
            coordc = Converter.node2Center(coordn)
            ac = C.getFields(Internal.__FlowSolutionCenters__, z)[0]
            if ac == []: ac = coordc
            else: Converter._addVars([ac, coordc])
            ac = Initiator.initLamb(ac, position, Gamma, MInf)
            ac = Converter.rmVars(ac, ['x', 'y', 'z'])
            z = C.setFields([ac], z, 'centers')
    return None

def initVisbal(t, position=(0.,0.), Gamma=2., MInf=0.5, loc='nodes'):
    """Init the array defining a grid with a Visbal vortex of
    intensity Gamma and position (x0,y0).
    Returns the array of the grid + cfd field in centers
    Usage: initVisbal(array, (x0,y0), Gamma, MInf)"""
    tp = Internal.copyRef(t)
    _initVisbal(tp, position, Gamma, MInf, loc)
    return tp

def _initVisbal(t, position=(0.,0.), Gamma=2., MInf=0.5, loc='nodes'):    
    nodes = Internal.getZones(t)
    for z in nodes:
        coordn = C.getFields(Internal.__GridCoordinates__, z)
        if coordn == []:
            print ('Warning: initVisbal: zone '+z[0]+' has no coordinates. Skipped...')
            continue
        coordn = coordn[0]
        if loc == 'nodes':
            a = C.getFields(Internal.__FlowSolutionNodes__, z)[0]
            if a == []: a = coordn
            else: Converter._addVars([a, coordn])
            a = Initiator.initVisbal(a, position, Gamma, MInf)
            z = C.setFields([a], z, 'nodes')
        else:
            coordc = Converter.node2Center(coordn)
            ac = C.getFields(Internal.__FlowSolutionCenters__, z)[0]
            if ac == []: ac = coordc
            else: Converter._addVars([ac, coordc])
            ac = Initiator.initVisbal(ac, position, Gamma, MInf)
            ac = Converter.rmVars(ac, ['x', 'y', 'z'])
            z = C.setFields([ac], z, 'centers')
    return None

def initYee(t, position=(0.,0.), Gamma=2., MInf=0.5, loc='nodes'):
    """Init the array defining a grid with a Yee vortex of
    intensity Gamma and position (x0,y0).
    Returns the array of the grid + cfd field in centers
    Usage: initYee(array, (x0,y0), Gamma, MInf)"""
    tp = Internal.copyRef(t)
    _initYee(tp, position, Gamma, MInf, loc)
    return tp

def _initYee(t, position=(0.,0.), Gamma=2., MInf=0.5, loc='nodes'):
    nodes = Internal.getZones(t)
    for z in nodes:
        coordn = C.getFields(Internal.__GridCoordinates__, z)
        if coordn == []:
            print ('Warning: initYee: zone '+z[0]+' has no coordinates. Skipped...')
            continue
        coordn = coordn[0]
        if loc == 'nodes':
            a = C.getFields(Internal.__FlowSolutionNodes__, z)[0]
            if a == []: a = coordn
            else: Converter._addVars([a, coordn])
            a = Initiator.initYee(a, position, Gamma, MInf)
            z = C.setFields([a], z, 'nodes')
        else:
            coordc = Converter.node2Center(coordn)
            ac = C.getFields(Internal.__FlowSolutionCenters__, z)[0]
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
    nodes = Internal.getZones(t)
    for z in nodes:
        coordn = C.getFields(Internal.__GridCoordinates__, z)
        if (coordn == []):
            print ('Warning: initScully: zone '+z[0]+' has no coordinates. Skipped...')
            continue
        coordn = coordn[0]
        if loc == 'nodes':
            a = C.getFields(Internal.__FlowSolutionNodes__, z)[0]
            if a == []: a = coordn
            else: Converter._addVars([a, coordn])
            a = Initiator.initScully(a, position, Gamma, coreRadius, MInf,
                                     model)
            C.setFields([a], z, 'nodes')
        else:
            coordc = Converter.node2Center(coordn)
            ac = C.getFields(Internal.__FlowSolutionCenters__, z)[0]
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
    nodes = Internal.getZones(t1)
    nodes2 = Internal.getZones(t2)
    for c in range(len(nodes)):
        z1 = nodes[c]
        if loc == 'centers':
            a1 = C.getAllFields(z1, 'centers')[0]
            x1 = C.getFields(Internal.__GridCoordinates__, z1)[0]
            x1 = Converter.node2Center(x1)
            a1 = Converter.addVars([x1, a1])
            z2 = nodes2[c]
            a2 = C.getAllFields(z2, 'centers')[0]
            x2 = C.getFields(Internal.__GridCoordinates__, z2)[0]
            x2 = Converter.node2Center(x2)
            a2 = Converter.addVars([x2, a2])
            ret = Initiator.overlayField(a1, a2, MInf)
            C.setFields([ret], z1, 'centers')
        else:
            a1 = C.getAllFields(z1, 'nodes')[0]
            z2 = nodes2[c]
            a2 = C.getAllFields(z2, 'nodes')[0]
            ret = Initiator.overlayField(a1, a2, MInf)
            C.setFields([ret], z1, 'nodes')
    return None
