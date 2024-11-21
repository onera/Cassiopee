"""Module of computation of distance to walls.
"""
from . import dist2walls
try: import Converter as C
except ImportError: raise ImportError("Dist2Walls: requires Converter modules.")
__version__ = '4.0'
__author__ = "Stephanie Peron, Christophe Benoit, Pascal Raud, Sam Landier"

try: range = xrange
except: pass

# Types de solver pour Eikonal
fmm=0; fim=1; fim_old=2 # Temporaire

def signDistance__(zones, distances, bodies, loc, dimPb):
    import Connector as X
    cellns = []
    if loc == 'nodes':  # blankingType=0: node_in
        cellns = C.initVars(zones, 'cellN', 1.)
        cellns = C.extractVars(cellns, ['cellN'])
        for b in bodies:
            cellns = X.blankCells(zones, cellns, [b], blankingType=0,
                                  dim=dimPb)
    else:  # blankingType=2: center_in
        cellns = C.initVars(distances, 'cellN', 1.)
        cellns = C.extractVars(cellns, ['cellN'])
        for b in bodies:
            cellns = X.blankCells(zones, cellns, [b], blankingType=2,
                                  dim=dimPb)
    # calcul de distances
    distances = C.addVars([distances, cellns])
    distances = C.initVars(distances, '{TurbulentDistance}={TurbulentDistance}*({cellN}>0.)+(-1.)*{TurbulentDistance}*({cellN}<1.)')
    distances = C.extractVars(distances, ['TurbulentDistance'])
    return distances


# ==============================================================================
# Calcul de la distance a des surfaces pour une liste de zones, a partir d'une
# liste de surfaces (bodies) et d'une liste de champs celln associee a bodies
# ==============================================================================
def distance2Walls(zones, bodies, flags=None, cellnbodies=[], type='ortho',
                   loc='centers', signed=0, dim=3, isIBM_F1=False, dTarget=1000.):
    """Compute distance to walls.
       Usage: distance2Walls(zones, bodies, cellnbodies, type, loc, signed, dim)"""

    # firewalls
    if len(zones) == 0: return
    if bodies == []:
        print('Warning: distance2Walls: no body defined, no distance computed.')
        return zones

    onezone = 0
    if not isinstance(zones[0], list):
        onezone = 1
        zones = [zones]

    if loc != 'centers' and loc != 'nodes':
        raise ValueError("distance2Walls: loc must be centers or nodes.")
    if type != 'ortho' and type != 'mininterf' and type != 'mininterf_ortho' and type != 'ortho_local':
        raise ValueError("distance2Walls: type must be ortho, mininterf, mininterf_ortho or ortho_local.")

    # Recuperation du cellN en noeuds ou centres selon loc
    bodies0 = []  # argument dans la fonction c associe a bodies et cellN
    if cellnbodies == []:
        bodies0 = C.initVars(bodies, 'cellN', 1.)
    elif len(cellnbodies) != len(bodies):
        print('Warning: distance2Walls: cellN not defined for some bodies: cellN set to 1 for invalid body zones.')
        bodies0 = C.initVars(bodies, 'cellN', 1.)
    else:
        bodies0 = bodies
        for c in range(len(bodies0)):
            if cellnbodies[c] == []:
                bodies0[c] = C.initVars(bodies0[c], 'cellN', 1.)
            elif bodies0[c][1].shape[1] == cellnbodies[c][1].shape[1]:
                bodies0[c] = C.addVars([bodies0[c], cellnbodies[c]])
            else:
                print('Warning: distance2Walls: bodies and celln must be of same dimensions: cellN set to 1 for invalid body zones.')
                bodies0[c] = C.initVars(bodies0[c], 'cellN', 1.)
    # conversion en triangles de bodies (important pour les champs aux centres
    # sur maillages en centres)
    bodies0 = C.convertArray2Tetra(bodies0, split='withBarycenters')

    # calcul de la distance a la paroi localisee aux centres ou aux noeuds
    dist = []
    if loc == 'nodes':
        if type == 'ortho' or type == 'ortho_local':
            if flags is not None:
                for noz in range(len(zones)):
                    if flags[noz] != []:
                        zones[noz] = C.addVars([zones[noz], flags[noz]])
            isminortho = 0
            if type == "ortho_local": isminortho = 1
            dist = dist2walls.distance2WallsOrtho(zones, bodies0, isminortho, int(isIBM_F1), dTarget)
        else:
            isminortho = 0
            if type == 'mininterf_ortho': isminortho = 1
            dist = dist2walls.distance2Walls(zones, bodies0, isminortho)

    elif loc == 'centers':
        zonesc = C.node2Center(zones)
        if type == 'ortho' or type == 'ortho_local':
            if flags is not None:
                for noz in range(len(zonesc)):
                    if flags[noz] != []:
                        zonesc[noz] = C.addVars([zonesc[noz], flags[noz]])
            isminortho = 0
            if type == "ortho_local": isminortho = 1
            dist = dist2walls.distance2WallsOrtho(zonesc, bodies0, isminortho, int(isIBM_F1), dTarget)
        else:
            isminortho = 0
            if type == 'mininterf_ortho': isminortho = 1
            dist = dist2walls.distance2Walls(zonesc, bodies0, isminortho)

    # distance signee
    if signed == 1:
        dist = signDistance__(zones, dist, bodies0, loc, dim)

    if onezone == 1: return dist[0]
    else: return dist

# ========================================================================================
# Solve eikonal on a Cartesian grid
# ========================================================================================
def eikonal(a, algo=fim_old):
    """Solve Eikonal equation on a Cartesian grid.
       Usage: eikonal(a)"""
    return dist2walls.eikonal(a, algo)
