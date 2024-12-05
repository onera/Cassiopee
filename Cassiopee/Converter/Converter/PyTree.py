from __future__ import print_function
try: range = xrange
except: pass

from . import Converter
from . import Internal
import KCore
import numpy, fnmatch, re, os.path
__version__ = Converter.__version__

# Variables globales
__ZoneNameServer__ = {}
__BCNameServer__ = {}
__BaseNameServer__ = {}

# Force R4 periodic nodes
FORCER4PERIODIC = True

#==============================================================================
# -- Gestion du nommage unique --
#==============================================================================

# -- getZoneName
# Retourne un nom de zone (unique par rapport au __ZoneNameServer__)
# IN: proposedName: nom propose
# OUT: nouveau nom unique
def getZoneName(proposedName):
    global __ZoneNameServer__
    (name, __ZoneNameServer__) = getUniqueName(proposedName, __ZoneNameServer__)
    return name

# Enregistre les noms de zones de t dans le __ZoneNameServer__
def registerZoneNames(t):
    global __ZoneNameServer__
    nodes = Internal.getZones(t)
    for i in nodes: __ZoneNameServer__[i[0]] = 0

# -- getBCName
# Retourne un nom de BC unique, stocke les noms de BC_t, GridConnectivity_t
# IN: proposedName: nom propose
def getBCName(proposedName):
    global __BCNameServer__
    (name, __BCNameServer__) = getUniqueName(proposedName, __BCNameServer__)
    return name

# Retourne le denier nom de BC propose par le serveur pour le nom proposedName
def getLastBCName(proposedName):
    return getLastName(proposedName, __BCNameServer__)

# Enregistre les noms des BCs de t dans le __BCNameServer__
def registerBCNames(t):
    global __BCNameServer__
    zones = Internal.getZones(t)
    for z in zones:
        zoneBCs = Internal.getNodesFromType1(z, 'ZoneBC_t')
        for zbc in zoneBCs:
            nodes = Internal.getNodesFromType1(zbc, 'BC_t')
            for i in nodes: __BCNameServer__[i[0]] = 0
        zoneGCs = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
        for zgc in zoneGCs:
            nodes = Internal.getNodesFromType1(zgc, 'GridConnectivity1to1_t')
            for i in nodes: __BCNameServer__[i[0]] = 0
            nodes = Internal.getNodesFromType1(zgc, 'GridConnectivity_t')
            for i in nodes: __BCNameServer__[i[0]] = 0

# -- getBaseName
# Retourne un nom de base (unique)
# IN: proposedName: nom propose
# OUT: nouveau nom unique
def getBaseName(proposedName):
    global __BaseNameServer__
    (name, __BaseNameServer__) = getUniqueName(proposedName, __BaseNameServer__)
    return name

# Enregistre les noms de base de t dans le __BaseNameServer__
def registerBaseNames(t):
    global __BaseNameServer__
    nodes = Internal.getBases(t)
    for i in nodes:
        __BaseNameServer__[i[0]] = 0

# Enregistre les Zone names, les Base names, les BC names
def registerAllNames(t):
    registerZoneNames(t); registerBaseNames(t); registerBCNames(t)

# Clear all names
def clearAllNames():
    global __ZoneNameServer__
    __ZoneNameServer__ = {}
    global __BCNameServer__
    __BCNameServer__ = {}
    global __BaseNameServer__
    __BaseNameServer__ = {}
    return None

# Retourne proposedName ou proposedName.count
def getUniqueName(proposedName, server):
    namespl = proposedName.rsplit('.', 1)
    if len(namespl) == 2:
        try: c = int(namespl[1]); name = namespl[0]
        except: name = proposedName
    else: name = proposedName
    if name not in server:
        server[name] = 0
        return (name, server)
    else:
        c = server[name]; ret = 1
        while ret == 1:
            name2 = '%s.%d'%(name,c)
            if name2 not in server: ret = 0
            else: ret = 1
            c += 1
        server[name2] = 0
        server[name] = c
        return (name2, server)

# Retourne le dernier nom fourni par le serveur
def getLastName(proposedName, server):
    namespl = proposedName.rsplit('.', 1)
    if len(namespl) == 2:
        try: c = int(namespl[1]); name = namespl[0]
        except: name = proposedName
    else: name = proposedName
    if name not in server:
        return None
    else:
        c = server[name]
        if c == 0: return name
        else: return name+'.'+str(c-1)

#==============================================================================
# -- pyTree informations --
#==============================================================================

# -- printTree
# print a tree to file or screen (obsolete, use Internal)
def printTree(t, file=None, stdOut=None, editor=False):
    Internal.printTree(t, file, stdOut, editor)

# -- getNPts
# Retourne le nombre de pts dans t
def getNPts(t):
    """Return the number of points in t.
    Usage: getNPts(t)"""
    zones = Internal.getZones(t)
    npts = 0
    for z in zones:
        dim = Internal.getZoneDim(z)
        if dim[0] == 'Structured': npts += dim[1]*dim[2]*dim[3]
        elif dim[0] == 'Unstructured': npts += dim[1]
    return npts

# -- getNCells
# Retourne le nombre de cellules dans t
def getNCells(t):
    """Return the number of cells in t.
    Usage: getNCells(t)"""
    zones = Internal.getZones(t)
    ncells = 0
    for z in zones:
        dim = Internal.getZoneDim(z)
        if dim[0] == 'Structured':
            ni1 = max(dim[1]-1,1); nj1 = max(dim[2]-1,1); nk1 = max(dim[3]-1,1)
            ncells += ni1*nj1*nk1
        elif dim[0] == 'Unstructured': ncells += dim[2]
    return ncells

# -- convertPyTree2ZoneNames
# Return the list of path of zones of a python tree as a list of strings
def convertPyTree2ZoneNames(t):
    """Return the list of zone names of a py tree.
    Usage: convertPyTree2ZoneNames(t)"""
    l = []
    toptree = Internal.isTopTree(t)
    if toptree:
        bases = Internal.getBases(t)
        for b in bases:
            baseName = b[0]
            zones = Internal.getNodesFromType1(b, 'Zone_t')
            for z in zones:
                zoneName = z[0]
                l.append('%s/%s'%(baseName,zoneName))
    return l

# -- getStdNodesFromName. Applique uniquement sur une zone.
# Retourne les noeuds des champs de nom 'name' dans les conteneurs standards
# Accepte les variables 'centers:var', les noms de containers (GridCoordinates)
# ou les noms Container/name
# Si pas trouve, retourne []
def getStdNodesFromName(z, name):
    loc = '*'
    v = name.split(':',1)
    if len(v) > 1:
        if v[0] == 'centers' or v[0] == 'nodes': var = v[1]; loc = v[0]
        else: var = name
    else: var = name
    v = name.split('/')
    if len(v) > 1: var = v[1]; loc = v[0]
    result = []
    if loc == 'nodes' or loc == '*':
        node = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
        if node is not None:
            if name == Internal.__GridCoordinates__: result.append(node)
            for j in node[2]:
                if j[0] == var: result.append(j); break
        node = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
        if node is not None:
            if name == Internal.__FlowSolutionNodes__: result.append(node)
            for j in node[2]:
                if j[0] == var: result.append(j); break
    if loc == 'centers' or loc == '*':
        node = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
        if node is not None:
            if name == Internal.__FlowSolutionCenters__: result.append(node)
            for j in node[2]:
                if j[0] == var: result.append(j); break
    if loc != 'nodes' and loc != 'centers': # given container name
        node = Internal.getNodeFromName1(z, loc)
        if node is not None:
            for j in node[2]:
                if j[0] == var: result.append(j); break
    return result

# -- getVarNames
# Retourne une liste des noms de variables presents pour chaque zone
# avec leur localisation (mode 0)
# localisation=centers: ou nom du container:
# Ex: pour t contenant 2 zones:
# [ ['CoordinateX','CoordinateY'], ['centers:F'] ]
# Si excludeXYZ=True, les coordonnees ne font pas partie des champs retournes
# Si loc='both', retourne les variables en centres et en noeuds
# Sinon loc='centers' ou loc='nodes'
# Si mode=1: retourne l'union des variables de toutes les zones
# Si mode=2: retourne les variables communes a toutes les zones
def getVarNames(t, excludeXYZ=False, loc='both', mode=0):
    allvars = []
    nodes = Internal.getZones(t)
    for z in nodes:
        varlist = []
        if not excludeXYZ:
            nodesGC = Internal.getNodesFromName1(z, Internal.__GridCoordinates__)
            for i in nodesGC:
                nodesVar = Internal.getNodesFromType1(i, 'DataArray_t')
                for j in nodesVar: varlist.append(j[0])

        if loc == 'nodes' or loc == 'both':
            nodesSol = Internal.getNodesFromName1(z, Internal.__FlowSolutionNodes__)
            location = ''
            for i in nodesSol:
                nodesVar = Internal.getNodesFromType1(i, 'DataArray_t')
                for j in nodesVar: varlist.append(location+j[0])

        if loc == 'centers' or loc == 'both':
            nodesSol = Internal.getNodesFromName1(z, Internal.__FlowSolutionCenters__)
            location = 'centers:'
            for i in nodesSol:
                nodesVar = Internal.getNodesFromType1(i, 'DataArray_t')
                for j in nodesVar: varlist.append(location+j[0])
        allvars.append(varlist)

        if mode == 1: # all vars everywhere
            out = []
            for zvars in allvars:
                if out == []: out = zvars
                else:
                    for v in zvars:
                        if v not in out: out.append(v)
            allvars = [out]

        if mode == 2: # common vars only
            out = []
            d = {}
            for zvars in allvars:
                for v in zvars:
                    if v in d: d[v] += 1
                    else: d[v] = 1
            out = []
            for k in d:
                if d[k] == len(allvars): out.append(k)
            allvars = [out]
    return allvars

# -- isNamePresent
# isNamePresent dans un arbre t
# Retourne -1: la variable n'est presente dans aucune zone de t
# Retourne 0: la variable est presente dans au moins une zone, mais pas toutes.
# Retourne 1: la variables est presente dans toutes les zones
def isNamePresent(t, varname):
    v = varname.split(':',1)
    if len(v) > 1 and v[0] == 'nodes': varname = v[1]
    zvars = getVarNames(t)
    if len(zvars) == 0: return -1
    one = 0; n = 0
    for z in zvars:
        found = 0
        for v in z:
            if v == varname: found = 1; one = 1; n += 1; break
        if found == 0:
            if one == 1: return 0
    if one == 0: return -1
    elif n != len(zvars): return 0
    else: return 1

# -- getNobNozOfZone
# IN: a: zone
# IN: t: le top tree de a
# OUT: (nob, noz): no de la base et no de la zone a dans toptree
# If zone can not be found, return (-1, -1)
# Must be faster than getParentOfNode
def getNobNozOfZone(a, t):
    """Return the (nob, noz) of a in t.
     Usage: getNobNozOfZone(a, t)"""
    nob = 0; noz = 0
    for b in t[2]:
        noz = 0
        for z in b[2]:
            if id(z) == id(a): return (nob, noz)
            noz += 1
        nob += 1
    return (-1, -1)

# -- getNobOfBase
# IN: b: base
# IN: t: le top tree de a
# OUT: nob: no de la base dans le top tree
# If base can not be found, return -1
# Must be faster than getParentOfNode
def getNobOfBase(base, t):
    """Return the nob of a base in t.
    Usage: getNobOfBase(base, t)"""
    nob = 0
    for b in t[2]:
        if id(b) == id(base): return nob
        nob += 1
    return -1

# -- GetZoneNames
# Retourne une liste contenant le nom des zones de l'arbre, trie ainsi:
# En premier, les zones structurees
# En second, les zones non structurees
def getZoneNames(t, prefixByBase=True):
    bases = Internal.getBases(t)
    names = []
    if bases == [] or not prefixByBase: # Juste zone name
        zones = Internal.getZones(t)
        for z in zones:
            type = Internal.getZoneType(z)
            if type == 1: names.append(z[0])
        for z in zones:
            type = Internal.getZoneType(z)
            if type == 2: names.append(z[0])
    else: # Base/Zone name
        for b in bases:
            baseName = b[0]
            zones = Internal.getNodesFromType1(b, 'Zone_t')
            for z in zones:
                type = Internal.getZoneType(z)
                if type == 1: names.append(baseName+Internal.SEP1+z[0])
        for b in bases:
            baseName = b[0]
            zones = Internal.getNodesFromType1(b, 'Zone_t')
            for z in zones:
                type = Internal.getZoneType(z)
                if type == 2: names.append(baseName+Internal.SEP1+z[0])
    return names

# -- getConnectedZones
# IN: a: zone, tree, base, liste de zones
# IN: toptree: le toptree de a
# OUT: la liste des zones connectees a a par:
# - des matchs
# - des near matchs
# Attention: le nom des zones doit etre unique dans l'arbre!
def getConnectedZones(a, topTree, type='all'):
    """Return the list of zones connected to a through match or nearMatch.
    Usage: getConnectedZones(a, topTree, type)"""
    zones = Internal.getZones(a)
    out = set()
    match = 0
    if type == 'all' or type == 'BCMatch': match = 1
    nearMatch = 0
    if type == 'all' or type == 'BCMatch': nearMatch = 1

    if match == 1:
        for z in zones:
            nodes = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
            for n in nodes:
                donor = Internal.getValue(n)
                out.add(donor)

    if nearMatch == 1:
        for z in zones:
            nodes = Internal.getNodesFromType2(z, 'GridConnectivity_t')
            for n in nodes:
                donor = Internal.getValue(n)
                out.add(donor)

    # Retrouve les zones
    outz = []
    zones = Internal.getZones(topTree)
    for i in out:
        for z in zones:
            if z[0] == i: outz.append(z); break
    return outz

#==============================================================================
# -- Get / Set field values --
#==============================================================================

# -- getValue
# Retourne la valeur d'un champ standard pour un indice donne
def getValue(t, var, ind):
    """Return the values for a point of index ind or (i,j,k)...
    Usage: getValue(t, var, ind)"""
    if isinstance(var, list):
        out = []
        for i in var:
            ret = getValue__(t, i, ind)
            if isinstance(ret, list): out += ret
            else: out.append(ret)
        return out
    else: return getValue__(t, var, ind)

def getValue__(t, var, ind):
    result = []
    zones = Internal.getZones(t)
    if zones == []: raise ValueError("getValue: not a zone node.")
    z = zones[0]

    # localisation
    u = var.split(':')
    if len(u) >= 2: loc = u[0]; var = u[1]
    else: loc = '*'; var = u[0]

    # dim
    dim = Internal.getZoneDim(z); cellDim = dim[4]
    if dim[0] == 'Structured': # output im,jm,km
        ni = dim[1]; nj = dim[2]; ni1 = max(ni-1,1); nj1 = max(nj-1,1)
        ninj = ni*nj; ni1nj1 = ni1*nj1
        if isinstance(ind, tuple):
            if len(ind) == 3:
                im = ind[0]-1; jm = ind[1]-1; km = ind[2]-1
            elif len(ind) == 2:
                im = ind[0]-1; jm = ind[1]-1; km = 0
            elif len(ind) == 1:
                im = ind[0]-1; jm = 0; km = 0
            else:
                raise ValueError("getValue: too much values in tuple.")
        else:
            if loc == 'nodes' or loc == '*':
                km = ind // (ninj)
                jm = (ind - km*ninj) // ni
                im = ind - jm*ni - km*ninj
            else:
                km = ind // (ni1nj1)
                jm = (ind - km*ni1nj1) // ni1
                im = ind - jm*ni1 - km*ni1nj1

        # GridCoordinates
        v = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
        if v is not None:
            for i in v[2]:
                if ((i[0] == var or var == Internal.__GridCoordinates__)
                    and (loc == 'nodes' or loc == '*')
                    and i[3] == 'DataArray_t'):
                    if cellDim == 3: result.append(i[1][im,jm,km])
                    elif cellDim == 2: result.append(i[1][im,jm])
                    else: result.append(i[1][im])

        # FlowSolutionNodes
        v = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
        if v is not None:
            for i in v[2]:
                if ((i[0] == var or var == Internal.__FlowSolutionNodes__)
                    and (loc == 'nodes' or loc == '*')
                    and i[3] == 'DataArray_t'):
                    if cellDim == 3: result.append(i[1][im,jm,km])
                    elif cellDim == 2: result.append(i[1][im,jm])
                    else: result.append(i[1][im])

        # FlowSolutionCenters
        v = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
        if v is not None:
            for i in v[2]:
                if ((i[0] == var or var == Internal.__FlowSolutionCenters__)
                    and (loc == 'centers' or loc == '*')
                    and i[3] == 'DataArray_t'):
                    if cellDim == 3: result.append(i[1][im,jm,km])
                    elif cellDim == 2: result.append(i[1][im,jm])
                    else: result.append(i[1][im])

    else: # output index
        # GridCoordinates
        v = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
        if v is not None:
            for i in v[2]:
                if ((i[0] == var or var == Internal.__GridCoordinates__)
                    and (loc == 'nodes' or loc == '*')
                    and i[3] == 'DataArray_t'):
                    result.append(i[1][ind])

        # FlowSolutionNodes
        v = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
        if v is not None:
            for i in v[2]:
                if ((i[0] == var or var == Internal.__FlowSolutionNodes__)
                    and (loc == 'nodes' or loc == '*')
                    and i[3] == 'DataArray_t'):
                    result.append(i[1][ind])

        # FlowSolutionCenters
        v = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
        if v is not None:
            for i in v[2]:
                if ((i[0] == var or var == Internal.__FlowSolutionCenters__)
                    and (loc == 'centers' or loc == '*')
                    and i[3] == 'DataArray_t'):
                    result.append(i[1][ind])

    if len(result) == 1: return result[0]
    else: return result

# -- setValue
# Set a value in a zone field pour l'indice donne
# Attention: pas de copie ici, c'est directement le champ qui est modifie
def setValue(t, var, ind, val):
    """Set the values in an array for a point of index ind or (i,j,k)...
    Usage: setValue(t, var, ind, value)"""
    nodes = Internal.getZones(t)
    if nodes == []: raise ValueError("setValue: not a zone node.")
    z = nodes[0]

    # localisation
    u = var.split(':')
    if len(u) >= 2: loc = u[0]; var = u[1]
    else: loc = '*'; var = u[0]

    c = 0
    if not isinstance(val, list): val = [val]

    # dim
    dim = Internal.getZoneDim(z); cellDim = dim[4]
    if dim[0] == 'Structured': # output im,jm,km
        ni = dim[1]; nj = dim[2]; ni1 = max(ni-1,1); nj1 = max(nj-1,1)
        ninj = ni*nj; ni1nj1 = ni1*nj1
        if isinstance(ind, tuple):
            ll = len(ind)
            if ll == 3:
                im = ind[0]-1; jm = ind[1]-1; km = ind[2]-1
            elif ll == 2:
                im = ind[0]-1; jm = ind[1]-1; km = 0
            elif ll == 1:
                im = ind[0]-1; jm = 0; km = 0
            else:
                raise ValueError("setValue: too much values in tuple.")
        else:
            if loc == 'nodes' or loc == '*':
                km = ind // (ninj)
                jm = (ind - km*ninj) // ni
                im = ind - jm*ni - km*ninj
            else:
                km = ind // (ni1nj1)
                jm = (ind - km*ni1nj1) // ni1
                im = ind - jm*ni1 - km*ni1nj1

        # GridCoordinates
        v = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
        if v is not None:
            for i in v[2]:
                if ((i[0] == var or var == Internal.__GridCoordinates__)
                and (loc == 'nodes' or loc == '*')):
                    if cellDim == 3:
                        i[1][im,jm,km] = val[c]; c += 1
                    elif cellDim == 2:
                        i[1][im,jm] = val[c]; c += 1
                    else: i[1][im] = val[c]; c += 1

        # FlowSolutionNodes
        v = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
        if v is not None:
            for i in v[2]:
                if ((i[0] == var or var == Internal.__FlowSolutionNodes__)
                and (loc == 'nodes' or loc == '*')):
                    if cellDim == 3:
                        i[1][im,jm,km] = val[c]; c += 1
                    elif cellDim == 2:
                        i[1][im,jm] = val[c]; c += 1
                    else: i[1][im] = val[c]; c += 1

        # FlowSolutionCenters
        v = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
        if v is not None:
            for i in v[2]:
                if ((i[0] == var or var == Internal.__FlowSolutionCenters__)
                and (loc == 'centers' or loc == '*')):
                    if cellDim == 3:
                        i[1][im,jm,km] = val[c]; c += 1
                    elif cellDim == 2:
                        i[1][im,jm] = val[c]; c += 1
                    else: i[1][im] = val[c]; c += 1

    else:
        # GridCoordinates
        v = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
        if v is not None:
            for i in v[2]:
                if ((i[0] == var or var == Internal.__GridCoordinates__)
                and (loc == 'nodes' or loc == '*')):
                    i[1][ind] = val[c]; c += 1

        # FlowSolutionNodes
        v = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
        if v is not None:
            for i in v[2]:
                if ((i[0] == var or var == Internal.__FlowSolutionNodes__)
                and (loc == 'nodes' or loc == '*')):
                    i[1][ind] = val[c]; c += 1

        # FlowSolutionCenters
        v = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
        if v is not None:
            for i in v[2]:
                if ((i[0] == var or var == Internal.__FlowSolutionCenters__)
                and (loc == 'centers' or loc == '*')):
                    i[1][ind] = val[c]; c += 1
    return None

def _getIndexField(t):
    return _TZGC1(t, 'nodes', True, Converter.getIndexField)

def getIndexField(t):
    """Return the index field."""
    return TZGC1(t, 'nodes', True, Converter.getIndexField)

#==============================================================================
# -- Create PyTree --
#==============================================================================

# -- newPyTree
def newPyTree(args=[]):
    """Create a new PyTree.
    Usage: newPyTree(['Base', z1, z2,  ...])"""
    t = Internal.createRootNode()
    t[2].append(Internal.createCGNSVersionNode())
    l = len(args)
    base = None; cellDim = 3; name = 'Base'; noforce = []
    for i in args:
        if isinstance(i, str): # a baseName
            name = i
            base = Internal.createBaseNode(name, cellDim); t[2].append(base)
        elif isinstance(i, int): # base dim explicitely given
            if base is not None: base[1][0] = i; noforce.append(base[0])
        else:
            if len(i) == 4: # maybe a standard node
                if i[3] == 'Zone_t':
                    if base is None:
                        base = Internal.createBaseNode(name, cellDim); t[2].append(base)
                    base[2].append(i)
                elif i[3] == 'CGNSBase_t':
                    base = i; name = i[0]; t[2].append(base)
                else: # liste de zones
                    for z in i:
                        if isinstance(z, list):
                            if len(z) == 4 and z[3] == 'Zone_t':
                                if base is None:
                                    base = Internal.createBaseNode(name, cellDim); t[2].append(base)
                                base[2].append(z)
            else: # a list of zones?
                for z in i:
                    if isinstance(z, list):
                        if len(z) == 4 and z[3] == 'Zone_t':
                            if base is None:
                                base = Internal.createBaseNode(name, cellDim); t[2].append(base)
                            base[2].append(z)
    #for b in Internal.getBases(t):
    #  if b[0] not in noforce:
    #    Internal._correctBaseZonesDim(b, splitBases=False) # force the zones cellDim if cellDim not specified
    return t

# -- addBase2PyTree
# IN: cellDim=1,2,3 (dimension des cellules des zones de cette base)
def addBase2PyTree(t, baseName, cellDim=3):
    """Add a base name to a pyTree.
    Usage: addBase2PyTree(t, baseName, cellDim)"""
    if not isinstance(baseName, str):
        raise TypeError("addBase2PyTree: baseName must be a string.")
    a = Internal.copyRef(t)
    if a == []:
        a = Internal.createRootNode()
        a[2].append(Internal.createCGNSVersionNode())
        base = Internal.createBaseNode(baseName, cellDim)
        a[2].append(base)
    else: _addBase2PyTree(a, baseName, cellDim)
    return a

def _addBase2PyTree(a, baseName, cellDim=3):
    bases = Internal.getBases(a)
    found = False
    for i in bases:
        if i[0] == baseName: found = True
    if not found:
        base = Internal.createBaseNode(baseName, cellDim)
        a[2].append(base)
    return None

#==============================================================================
# -- delete certain nodes --
#==============================================================================

# -- deleteFlowSolutions__
# Enleve les noeuds FlowSolutionNodes, ou FlowSolutionCenters ou les 2
# loc='nodes', 'centers', 'both' respectivement
def deleteFlowSolutions__(a, loc='both'):
    b = Internal.copyRef(a)
    _deleteFlowSolutions__(b, loc)
    return b

def _deleteFlowSolutions__(a, loc='both'):
    if loc == 'centers' or loc == 'both':
        centers = Internal.getNodesFromName3(a, Internal.__FlowSolutionCenters__)
        if centers != []: _rmNodes(a, Internal.__FlowSolutionCenters__)
    if loc == 'nodes' or loc == 'both':
        nodes = Internal.getNodesFromName3(a, Internal.__FlowSolutionNodes__)
        if nodes != []: _rmNodes(a, Internal.__FlowSolutionNodes__)
    return None

# -- deleteGridConnectivity__
# Enleve les noeuds GridConnectivity
# type='None' par defaut: enleve le noeud complet
# type='BCMatch': enleve les BCMatch.
#   si kind='self' detruit uniquement les match sur la meme zone
#   si kind='other', n'enleve pas les raccords coincidents sur la meme zone
# type='BCOverlap' : enleve uniquement les BCOverlap
#   si kind='self' detruit les BCOverlap en autoattach (donorZoneName=zoneName)
#   si kind='other' detruit les BCOverlap pas autoattach mais definies par des zones donneuses
#   si kind='families': detruit les BCOverlap non autoattach definies par des familles de zones donneuses
# Par defaut, kind='all', detruit alors les 2 (self + other)
def deleteGridConnectivity__(a, type='None', kind='all'):
    if type == 'None':
        cn = Internal.getNodesFromName3(a, 'ZoneGridConnectivity')
        if cn != []: a = rmNodes(a, 'ZoneGridConnectivity')
    elif type == 'BCMatch':
        if kind == 'self': _deleteSelfBCMatch__(a)
        elif kind == 'other': _deleteOtherBCMatch__(a)
        else: a = rmBCOfType(a, 'BCMatch')
    elif type == 'BCOverlap':
        if kind == 'other': _deleteBCOverlapWithDonorZone__(a)
        elif kind == 'self': _deleteBCOverlapWithoutDonorZone__(a)
        else: a = rmBCOfType(a, 'BCOverlap')
    return a

def _deleteGridConnectivity__(a, type='None', kind='all'):
    if type == 'None':
        cn = Internal.getNodesFromName3(a, 'ZoneGridConnectivity')
        if cn != []: _rmNodes(a, 'ZoneGridConnectivity')
    elif type == 'BCMatch':
        if kind == 'self': _deleteSelfBCMatch__(a)
        elif kind == 'other': _deleteOtherBCMatch__(a)
        else: _rmBCOfType(a, 'BCMatch')
    elif type == 'BCOverlap':
        if kind == 'other': _deleteBCOverlapWithDonorZone__(a, removeDnrZones=True, removeDnrFam=False)
        elif kind == 'self': _deleteBCOverlapWithoutDonorZone__(a)
        else: _rmBCOfType(a, 'BCOverlap')
    return None

# enleve les BCMatch sur une meme zone
def _deleteSelfBCMatch__(a):
    zones = Internal.getZones(a)
    for z in zones:
        zoneName = z[0]
        nodes = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
        for i in nodes:
            donorName = i[1]
            if donorName == zoneName:
                (parent, d) = Internal.getParentOfNode(z, i)
                del parent[2][d]
    return None

# enleve les BCMatch entre deux zones differentes
def _deleteOtherBCMatch__(a):
    zones = Internal.getZones(a)
    for z in zones:
        zoneName = z[0]
        nodes = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
        for i in nodes:
            donorName = i[1]
            if donorName != zoneName:
                (parent, d) = Internal.getParentOfNode(z, i)
                del parent[2][d]
    return None

# enleve les BCOverlap de type autoattach
# identification par valeur du noeud = nom de la zone
# si les donneurs sont definis par familles de zone, les overlaps ne sont pas detruits
def _deleteBCOverlapWithoutDonorZone__(a):
    zones = Internal.getZones(a)
    for z in zones:
        zoneName = z[0]
        nodes = Internal.getNodesFromType2(z, 'GridConnectivity_t')
        for i in nodes:
            donorName = Internal.getValue(i)
            r = Internal.getNodeFromType1(i, 'GridConnectivityType_t')
            if r is not None:
                val = Internal.getValue(r)
                if val == 'Overset' and zoneName == donorName:
                    (parent, d) = Internal.getParentOfNode(z, i)
                    del parent[2][d]
    return None

# enleve les BCOverlap avec domaine attache
# si les donneurs sont des familles de zones alors ne sont pas detruits
def _deleteBCOverlapWithDonorZone__(a, removeDnrZones=True, removeDnrFam=True):
    zones = Internal.getZones(a)
    families = Internal.getNodesFromType2(a,'Family_t')
    for z in zones:
        zoneName = z[0]
        nodes = Internal.getNodesFromType2(z, 'GridConnectivity_t')
        for i in nodes:
            donorName = Internal.getValue(i)
            r = Internal.getNodeFromType1(i, 'GridConnectivityType_t')
            if r is not None:
                val = Internal.getValue(r)
                if val == 'Overset' and zoneName != donorName: # no auto attach
                    removeGC = False
                    typeOfDnr = 0 # 1: zone, 2: family
                    # check if donorName is a zone name
                    if Internal.getNodeFromName(zones,donorName) is not None: typeOfDnr = 1
                    else: # check if it is referenced to a family
                        if Internal.getNodeFromName(families,donorName) is not None: typeOfDnr=2

                    if (removeDnrZones and typeOfDnr==1) or (removeDnrFam and typeOfDnr==2) or (typeOfDnr==0):
                        (parent, d) = Internal.getParentOfNode(z, i)
                        del parent[2][d]
    return None

# -- deleteZoneBC
# Enleve les noeuds ZoneBC d'un type donne
# si type != 'None' enleve seulement le type de BC specifie
def deleteZoneBC__(a, type='None'):
    if type == 'None':
        bc = Internal.getNodesFromName3(a, 'ZoneBC')
        if bc != []: a = rmNodes(a, 'ZoneBC')
    else: a = rmBCOfType(a, type)
    return a

def _deleteZoneBC__(a, type='None'):
    if type == 'None':
        bc = Internal.getNodesFromName3(a, 'ZoneBC')
        if bc != []: _rmNodes(a, 'ZoneBC')
    else: _rmBCOfType(a, type)
    return None

# -- deleteAllBCAndSolutions__
# Enleve les noeuds ZoneBC, ZoneGridConnectivity, FlowSolution,
# FlowSolutionCenters
def deleteAllBCAndSolutions__(a):
    a = deleteGridConnectivity__(a)
    a = deleteZoneBC__(a)
    a = deleteFlowSolutions__(a)
    return a

def _deleteAllBCAndSolutions__(a):
    _deleteGridConnectivity__(a)
    _deleteZoneBC__(a)
    _deleteFlowSolutions__(a)
    return None

def deleteChimeraInfo__(a):
    ap = Internal.copyRef(a)
    _deleteChimeraInfo__(ap)
    return ap

def _deleteChimeraInfo__(a):
    Internal._rmNodesByName(a, 'OversetHoles')
    Internal._rmNodesByName(a, 'ID_*')
    return None

# -- delete Nodes specific to solvers --
def _deleteSolverNodes__(a):
    Internal._rmNodesByName(a, ':elsA#Hybrid') # elsA
    Internal._rmNodesByName(a, '.Solver#ownData') # Fast
    return None

# -- deleteEmptyZones
# Supprime les zones ayant un nombre de noeuds ou d'elements
# nul (non structure) ou un ni, nj, nk de 0 (structure)
# IN: t: tree, base, liste de zones
# OUT: isomorphe a l'entree moins les zones ayant aucun element
def deleteEmptyZones(t):
    """Delete zones with null number of points or elements."""
    tp = Internal.copyRef(t)
    _deleteEmptyZones(tp)
    return tp

def _deleteEmptyZones(t):
    zones = Internal.getZones(t)
    for z in zones:
        dim = Internal.getZoneDim(z)
        removeZ = False
        if dim[0] == 'Structured':
            if dim[1]*dim[2]*dim[3] == 0: removeZ = True
        else:
            if dim[3] == 'NODE':
                if dim[1] == 0: removeZ = True
            else:
                if dim[1]*dim[2] == 0: removeZ = True

        if removeZ:
            (p, c) = Internal.getParentOfNode(t, z)
            if id(p) == id(t): del p[c] # this is a list
            else: del p[2][c]
    return None

# -- rmNodes
def rmNodes(z, name):
    """Remove nodes name from z.
    Usage: rmNodes(z, name)"""
    zc = Internal.copyRef(z)
    _rmNodes(zc, name)
    return zc

def _rmNodes(z, name):
    zn = Internal.getZones(z)
    for i in zn:
        if isinstance(name, list):
            for v in name:
                nodes = Internal.getNodesFromName2(i, v)
                for j in nodes:
                    (parent, d) = Internal.getParentOfNode(i, j)
                    del parent[2][d]
        else:
            nodes = Internal.getNodesFromName2(i, name)
            for j in nodes:
                (parent, d) = Internal.getParentOfNode(i, j)
                del parent[2][d]
    return None

# Upgrade tree (applique apres lecture)
def _upgradeTree(t, uncompress=True, oldcompress=False):
    #Internal._adaptTypes(t)
    _relaxCGNSProfile__(t)
    Internal._correctPyTree(t, level=10) # force CGNS names
    Internal._correctPyTree(t, level=2) # force unique name
    Internal._correctPyTree(t, level=7) # create familyNames
    registerAllNames(t)
    if uncompress:
        try:
            import Compressor.PyTree as Compressor
            if oldcompress:
                Compressor._uncompressCartesian_old(t)
                Compressor._uncompressAll_old(t)
            else:
                Compressor._uncompressCartesian(t)
                Compressor._uncompressAll(t)
        except: pass
    return None

# Hack pour les arrays en centres avec sentinelle - 1.79769e+308
# copie sur les champs - plus utilise depuis que l'on sort directement
# les champs en centres
def hackCenters(a):
    varString = a[0].split(',')
    np = a[1]
    nfield = np.shape[0]; size = np.shape[1]
    ncenterFields = 0
    varStringN = ''; varStringC = ''
    for n in range(nfield):
        if np[n, size-1] == -1.79769e+308: # sentinelle for centers
            ncenterFields += 1; varStringC += varString[n]+','
        else: varStringN += varString[n]+','
    if len(varStringC)>0: varStringC = varStringC[:-1]
    if len(varStringN)>0: varStringN = varStringN[:-1]
    if ncenterFields == 0: return a, []
    if isinstance(a[3], str):
        npts = Converter.getNPts(a); nelts = Converter.getNCells(a)
        a1 = [varStringN,numpy.zeros((nfield-ncenterFields,npts), numpy.float64),a[2],a[3]]
        a2 = [varStringC,numpy.zeros((ncenterFields,nelts), numpy.float64),a[2],a[3]+'*']
    else:
        ni = a[2], nj = a[3]; nk = a[4]; npts = ni*nj*nk
        ni1 = max(ni-1,1); nj1 = max(nj-1,1); nk1 = max(nk-1,1); nelts = ni1*nj1*nk1
        a1 = Converter.array(varStringN,ni,nj,nk); a2 = Converter.array(varStringC,ni1,nj1,nk1)

    nn = 0; nc = 0
    for n in range(nfield):
        if np[n, size-1] == -1.79769e+308: # sentinelle for centers
            a2[1][nc,0:nelts] = np[n,0:nelts]; nc += 1
        else:
            a1[1][nn,0:npts] = np[n,0:npts]; nn += 1
    return a1, a2

#==============================================================================
# -- File / pyTree conversions --
#==============================================================================

# -- convertFile2PyTree
def convertFile2PyTree(fileName, format=None, nptsCurve=20, nptsLine=2,
                       density=-1., skeletonData=None, dataShape=None,
                       links=None, skipTypes=None, uncompress=True,
                       hmax=0.0, hausd=1., grow=0.0, mergeTol=-1, occAlgo=4,
                       oldcompress=False, readMode=0, api=1):
    """Read a file and return a pyTree containing file data.
    Usage: convertFile2PyTree(fileName, format, options)"""
    if format is None:
        format = Converter.convertExt2Format__(fileName); autoTry = True
    else: autoTry = False
    exists = os.path.exists(fileName)
    if not exists: raise IOError("convertFile2PyTree: file %s not found."%fileName)

    if format == 'bin_cgns' or format == 'unknown':
        format = Converter.checkFileType(fileName)

    if format == 'bin_cgns' or format == 'bin_adf' or format == 'bin_hdf':
        try:
            t = Converter.converter.convertFile2PyTree(fileName, format, skeletonData, dataShape, links, skipTypes, readMode)
            t = Internal.createRootNode(children=t[2])
            _upgradeTree(t, uncompress, oldcompress)
            CAD = Internal.getNodeFromName1(t, 'CAD')
            if CAD is not None: # reload CAD
                file = Internal.getNodeFromName1(CAD, 'file')
                if file is not None: file = Internal.getValue(file)
                fmt = Internal.getNodeFromName1(CAD, 'format')
                if fmt is not None: fmt = Internal.getValue(fmt)
                if file is not None and fmt is not None:
                    import OCC.PyTree as OCC
                    import CPlot.Tk as CTK
                    hook = OCC.readCAD(file, fmt)
                    CTK.CADHOOK = hook
            return t
        except:
            if format == 'bin_cgns' or format == 'bin_adf':
                try:
                    t = Converter.converter.convertFile2PyTree(fileName, 'bin_hdf', skeletonData, dataShape, links, skipTypes, readMode)
                    t = Internal.createRootNode(children=t[2])
                    _upgradeTree(t, uncompress, oldcompress)
                    return t
                except: pass
            else: # adf par defaut
                try:
                    t = Converter.converter.convertFile2PyTree(fileName, 'bin_adf', skeletonData, dataShape, links, skipTypes, readMode)
                    t = Internal.createRootNode(children=t[2])
                    _upgradeTree(t)
                    return t
                except: pass
    elif format == 'unknown':
        try:
            t = Converter.converter.convertFile2PyTree(fileName, 'bin_adf', skeletonData, dataShape, links, skipTypes, readMode)
            t = Internal.createRootNode(children=t[2])
            _upgradeTree(t)
            return t
        except: pass
        try:
            t = Converter.converter.convertFile2PyTree(fileName, 'bin_hdf', skeletonData, dataShape, links, skipTypes, readMode)
            t = Internal.createRootNode(children=t[2])
            _upgradeTree(t)
            return t
        except: pass

    if occAlgo == 4 and (format == 'fmt_iges' or format == 'fmt_step'): # new cassiopee cad system
        import OCC.PyTree as OCC
        import CPlot.Tk as CTK
        hook = OCC.readCAD(fileName, format)
        if hmax == 0.:
                # auto setting
            (hmin,hmax,hausd) = OCC.occ.analyseEdges(hook)
        CTK.CADHOOK = hook
        t = OCC.meshAll(hook, hmax, hmax, hausd) # constant hmax
        _upgradeTree(t)
        return t

    if format == 'bin_pickle':
        try: import cPickle as pickle
        except: import pickle
        print('Reading %s (bin_pickle)...'%fileName, end='')
        try:
            file = open(fileName, 'rb')
            oldData = False
            if oldData: a = pickle.load(file, encoding='latin1')
            else: a = pickle.load(file)
            file.close()
        except:
            raise TypeError("convertFile2PyTree: file %s can not be read."%fileName)
        else:
            # modify i8/i4 types
            if Internal.E_NpyInt == numpy.int32: Internal._adaptTypes(a, convertI82I4=True, convertI42I8=False)
            elif Internal.E_NpyInt == numpy.int64: Internal._adaptTypes(a, convertI82I4=False, convertI42I8=True)
            print('done.')

        if Internal.isTopTree(a): # top tree
            registerAllNames(a)
            return a # OK
        ret = Internal.isStdNode(a)
        if ret != -2: # standard node
            t, ntype = Internal.node2PyTree(a)
            registerAllNames(t)
            return t # OK
        # sinon, c'est un arrays (traite dans la suite)

    zn = []; bcfaces = []; centerArrays = []; bcfields = []
    if format != 'bin_pickle': # autres formats
        if autoTry: format = None
        a = Converter.convertFile2Arrays(fileName, format, nptsCurve, nptsLine,
                                         density, zoneNames=zn, BCFaces=bcfaces,
                                         BCFields=bcfields,
                                         hmax=hmax, hausd=hausd, grow=grow,
                                         mergeTol=mergeTol, occAlgo=occAlgo,
                                         centerArrays=centerArrays, api=api)
    t = newPyTree([])
    base = 1; c = 0
    isBaseOfDim = [0 for _ in range(4)]
    unsEltsDim = {'NODE': 0, 'BAR': 1, 'TRI': 2, 'QUAD': 2}

    for c, i in enumerate(a):
        #a1, a2 = hackCenters(i)
        a1 = i
        if centerArrays != []: a2 = centerArrays[c]
        else: a2 = []

        if len(i) == 5: # Structure
            if i[3] == 1 and i[4] == 1: dim = 1
            elif i[4] == 1: dim = 2
            else: dim = 3
        else: # non structure
            dim = unsEltsDim.get(i[3].split(',')[0], 3) # BE, ME, (HOBE,NGON? TODO)

        suffix = '' if dim == 3 else str(dim)
        if not isBaseOfDim[dim]:
            t = addBase2PyTree(t, 'Base{}'.format(suffix), dim)
            isBaseOfDim[dim] = 1; base += 1
        z = Internal.createZoneNode(getZoneName('Zone'), a1, a2,
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)
        t[2][isBaseOfDim[dim]][2].append(z)

        if len(zn) > c: z[0] = zn[c]
        c += 1
    Internal._correctPyTree(t, level=10) # force CGNS names
    Internal._correctPyTree(t, level=2) # force unique name
    _addBCFaces(t, bcfaces)
    #_addBCFields(t, bcfields)
    Internal._correctPyTree(t, level=7) # create familyNames
    registerAllNames(t)
    return t

# Suppress invalid links
def checkLinks__(links, t):
    out = []
    for l in links:
        a = l[2]; b = l[3]
        if len(a) > 8 and a[0:8] == 'CGNSTree': l[2] = a.replace('CGNSTree', '')
        elif a[0] != '/': l[2] = '/'+a
        if len(b) > 8 and b[0:8] == 'CGNSTree': l[3] = b.replace('CGNSTree', '')
        elif b[0] != '/': l[3] = '/'+b
        # verifie que le noeud existe dans l'arbre (supprime)
        #if Internal.getNodeFromPath(t, b) is not None:
        # verifie que le noeud parent existe dans l'arbre
        #if Internal.getNodeFromPath(t, Internal.getPathAncestor(b)) is None:
        #    print("Warning: link %s is skipped, path ancestor not found."%b)
        #else:
        out.append(l)
    #for c, i in enumerate(links): print(out[c], links[c])
    return out

# Force periodic nodes to be R4 (old CGNS norm)
def _forceR4PeriodicNodes__(t):
    periodicNodeNames = ["RotationCenter", "RotationAngle", "RotationRateVector",
                         "Translation"]
    for name in periodicNodeNames:
        nodes = Internal.getNodesFromName(t, name)
        Internal._adaptValueType(nodes, numpy.float32)
    return None

# Force type of specific nodes as prescribed by v4 CGNS norm
def _forceCGNSProfile__(t):
    # CGNS version (R4)
    n = Internal.getNodeFromType1(t, "CGNSVersion_t")
    if n is not None: n[1] = n[1].astype(numpy.float32)
    #if FORCER4PERIODIC: _forceR4PeriodicNodes__(t)

    # CGNSBase_t, Elements_t, Rind_t, BaseIterativeData_t, Ordinal_t,
    # ConvergenceHistory_t (I4)
    nodes = Internal.getNodesFromType1(t, "CGNSBase_t")
    if nodes is not None:
        for n in nodes: n[1] = n[1].astype(numpy.int32)
    nodeTypes2 = ["Elements_t", "BaseIterativeData_t", "ConvergenceHistory_t"]
    for ntype in nodeTypes2:
        nodes = Internal.getNodesFromType2(t, ntype)
        Internal._adaptValueType(nodes, numpy.int32)
    zones = Internal.getZones(t)
    for z in zones:
        #nodes = Internal.getNodesFromType1(z, "Rind_t")
        #Internal._adaptValueType(nodes, numpy.int32)
        nodes = Internal.getNodesFromType(z, "Ordinal_t")
        Internal._adaptValueType(nodes, numpy.int32)
    return None

# Relax type of specific nodes as prescribed by v4 CGNS norm
def _relaxCGNSProfile__(t):
    if Internal.E_NpyInt == numpy.int32: return None
    # Rind_t cannot be forced to I4 internally
    zones = Internal.getZones(t)
    for z in zones:
        nodes = Internal.getNodesFromType1(z, "Rind_t")
        Internal._adaptValueType(nodes, Internal.E_NpyInt)
    return None

# -- convertPyTree2File
def convertPyTree2File(t, fileName, format=None, isize=4, rsize=8,
                       endian='big', colormap=0, dataFormat='%.9e ', links=[]):
    """Write a pyTree to a file.
    Usage: convertPyTree2File(t, fileName, format, options)"""
    if t == []: print('Warning: convertPyTree2File: nothing to write.'); return
    if links is not None: links = checkLinks__(links, t)
    if format is None:
        format = Converter.convertExt2Format__(fileName)
        if format == 'unknown': format = 'bin_cgns'
    if format == 'bin_cgns' or format == 'bin_adf' or format == 'bin_hdf':
        tp, ntype = Internal.node2PyTree(t)
        Internal._adaptZoneNamesForSlash(tp)
        Internal._correctBaseZonesDim(tp, splitBases=False)
        _forceCGNSProfile__(tp)
        Converter.converter.convertPyTree2File(tp[2], fileName, format, links)
    elif format == 'bin_pickle':
        try: import cPickle as pickle
        except: import pickle
        file = open(fileName, 'wb')
        print('Writing '+fileName+'...', end='')
        pickle.dump(t, file, protocol=pickle.HIGHEST_PROTOCOL); file.close()
        print('done.')
    else:
        tp = fillMissingVariables(t) # force all zones to have the same variables
        a = center2Node(tp, Internal.__FlowSolutionCenters__); tp = None
        a = getAllFields(a, 'nodes', api=3)
        a = Internal.clearList(a)
        zoneNames = getZoneNames(t, prefixByBase=False)
        BCFaces = getBCFaces(t, nameType=1)
        Converter.convertArrays2File(a, fileName, format, isize, rsize, endian,
                                     colormap, dataFormat, zoneNames, BCFaces)

# Fonction utilisee dans PPart
def convertFile2PartialPyTreeFromPath(fileName, Filter, comm=None,
                                      format=None, nptsCurve=20, nptsLine=2,
                                      density=-1., skeletonData=None, readMode=0):
    """Convert a file to pyTree.
    Usage: convertFile2PartialPyTree(fileName, format, options)"""
    if format is None:
        format = Converter.convertExt2Format__(fileName); autoTry = True
    else: autoTry = False
    try: file = open(fileName, 'r')
    except: raise IOError("convertFile2PartialPyTreeFromPath: file %s not found."%fileName)
    file.close()
    t = Converter.converter.convertFile2PartialPyTree(fileName, format, skeletonData, comm, Filter, readMode)
    return t

# Fonction utilisee dans PPart
def convertPyTree2FilePartial(t, fileName, comm, Filter, ParallelHDF=False,
                              format=None, links=[]):
    """Convert a pyTree to a file.
    Usage: convertPyTree2File(t, fileName, format, options)"""
    import collections
    # > GardeFou
    if t == []: print('Warning: convertPyTree2File: nothing to write.'); return
    format = 'bin_hdf'

    if not ParallelHDF:
            # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            # > Write Tree Data except Data in Filter
        SkeletonTree = Internal.copyRef(t)
        for path in Filter:
            #print(path)
            Node = Internal.getNodeFromPath(SkeletonTree, path)
            Node[1] = None
        # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Si MPI Mode Off (HDF Not Parallel)
        if comm.Get_rank() == 0:
            convertPyTree2File(SkeletonTree, fileName, format)
            # > Fill up Dimension
            skeletonData = None
            Converter.converter.convertPyTree2FilePartial(t, fileName, format, skeletonData, comm, Filter)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Wait for Skeleton write
        comm.barrier()

        # > Set skeletonData to Not None
        skeletonData = []

        # > Cette maniere de faire provoque de la non reproductibilite ...
        # Converter.converter.convertPyTree2FilePartial(t, fileName, format, skeletonData, comm, Filter)

        FilterSort = collections.OrderedDict(sorted(Filter.items()))  # Because sometimes proc have not the same order in key and HDF get in trouble !

        # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        # > On serialize
        for lock in range(comm.Get_size()):
            if lock == comm.Get_rank():
                Converter.converter.convertPyTree2FilePartial(t, fileName, format, skeletonData, comm, FilterSort)
            comm.barrier()
        # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    else:  # > HDF Parallel
        # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Write Tree Data except Data in Filter
        SkeletonTree = Internal.copyRef(t)
        for path in Filter:
            Node = Internal.getNodeFromPath(SkeletonTree, path)
            Node[1] = None
        # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        if comm.Get_rank() == 0:
            convertPyTree2File(SkeletonTree, fileName, format, links=links)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        # > On wait l'ecriture Skelette ...
        comm.barrier()

        # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Write data in filter in file (With creation of DataSpace)
        skeletonData = None  # Skeleton Data is inefective (Normaly)
        FilterSort = collections.OrderedDict(sorted(Filter.items()))  # Because sometimes proc have not the same order in key and HDF get in trouble !
        Converter.converter.convertPyTree2FilePartial(t, fileName, format, skeletonData, comm, FilterSort)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#==============================================================================
# -- Array / PyTree conversions --
#==============================================================================

# -- convertArrays2ZoneNode
# Cree une nouvelle zone a partir d'arrays
# IN: A: liste d'arrays. Contient un array ou 2 arrays (dans ce cas, le
# premier array est en noeuds et le deuxieme en centres).
# OUT: noeud zone
def convertArrays2ZoneNode(zoneName, A):
    """Convert arrays to a zone node.
    Usage: convertArrays2ZoneNode(zoneName, A)"""
    if len(A) == 1:
        z = Internal.createZoneNode(getZoneName(zoneName), A[0], [],
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)
    elif len(A) == 2:
        z = Internal.createZoneNode(getZoneName(zoneName), A[0], A[1],
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)
    else:
        raise TypeError("convertArrays2ZoneNode: A must be [a] or [node, center] when creating %s."%zoneName)
    return z

# -- getField
# Retourne des arrays correspondant a un nom de DataArray_t dans t
# Retourne une liste d'arrays par zone. Si une zone ne possede pas ce champ,
# retourne [] pour cette zone.
# Ex: getField('CoordinateX', t), retourne le champ CoordinateX pour toutes
# les zones de t
# Ex: getField('centers:Density', z), retourne le champ de Density pour la
# zone z.
def getField(name, t, api=1):
    zones = Internal.getZones(t)
    arrays = []
    loc = 'nodes'
    spl = name.split(':',1)
    if len(spl) > 1:
        if spl[0] == 'centers': loc = 'centers'; name = spl[1]
        elif spl[0] == 'nodes': name = spl[1]

    for z in zones:
        dim = Internal.getZoneDim(z)
        if dim[0] == 'Structured':
        #np = dim[1]*dim[2]*dim[3]
            connects = []
        else:
            #np = dim[1]
            connects = Internal.getElementNodes(z)

        info = z[2]; a = None
        if loc == 'nodes': # regarde les containeurs GridCoordinates et Nodes
            for i in info:
                if i[0] == Internal.__GridCoordinates__ or i[0] == Internal.__FlowSolutionNodes__:
                    info2 = i[2]
                    for j in info2:
                        if j[0] == name and j[3] == 'DataArray_t':
                            if api==1: r = Internal.convertDataNode2Array(j, dim, connects, 0)
                            elif api==2: r = Internal.convertDataNode2Array2(j, dim, connects, 0)
                            elif api==3: r = Internal.convertDataNode2Array3(j, dim, connects, 0)
                            else: raise ValueError('getField: unknow api.')
                            a = r[1]
                            if a is not None: break

        elif loc == 'centers':
            for i in info:
                if i[0] == Internal.__FlowSolutionCenters__:
                    info2 = i[2]
                    for j in info2:
                        if j[0] == name and j[3] == 'DataArray_t':
                            if api==1: r = Internal.convertDataNode2Array(j, dim, connects, 1)
                            elif api==2: r = Internal.convertDataNode2Array2(j, dim, connects, 1)
                            elif api==3: r = Internal.convertDataNode2Array3(j, dim, connects, 1)
                            else: raise ValueError('getField: unknow api.')
                            a = r[1]
                            if a is not None: break

        if a is not None: arrays.append(a)
        else: arrays.append([])
    return arrays

# -- getFields
# Retourne les arrays a partir d'un nom de conteneur a partir de n'importe quel
# noeud d'un arbre
# Retourne une liste de arrays
# IN: t: n'importe quel noeud de l'arbre
# IN: containerName: GridCoordinates, FlowSolution, FlowSolution#Centers (conteneur)
# ou liste de conteneurs homogenes en localisation
# IN: vars: optionel, liste des variables a recuperer (sans centers: car deja specifie dans container)
# IN: api=1, sortie array (avec copie), api=2, sortie array2 sans copie, api=3, sortie array3 sans copie
# OUT: arrays: solution (un par zone)
# OUT: peut contenir des arrays vides ([])
# Attention: il faut envoyer que des containeurs homogenes en localisation
def getFields(containerName, t, vars=None, api=1):
    zones = Internal.getZones(t)
    arrays = []
    if containerName == 'nodes': names = [Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__]
    elif containerName == 'centers': names = [Internal.__FlowSolutionCenters__]
    elif containerName == 'coords': names = [Internal.__GridCoordinates__]
    elif isinstance(containerName, list): names = containerName
    else: names = [containerName]
    for z in zones:
        dim = Internal.getZoneDim(z)
        if dim[0] == 'Structured':
            #np = dim[1]*dim[2]*dim[3]
            connects = []
        else:
            #np = dim[1]
            connects = Internal.getElementNodes(z)

        info = z[2]; out = []; loc = 0
        for i in info:
            if i[0] in names:
                locf = Internal.getNodeFromType1(i, 'GridLocation_t')
                if locf is not None and Internal.getValue(locf) == 'CellCenter': loc = 1
                for j in i[2]:
                    if j[3] == 'DataArray_t':
                        # Firewall for Coordinates in Fields containers
                        #if j[0] == 'CoordinateX' or j[0] == 'CoordinateY' or j[0] == 'CoordinateZ':
                        #  if i[0] == Internal.__GridCoordinates__: out.append(j)
                        #else: out.append(j)
                        if vars is None: out.append(j)
                        else:
                            if j[0] in vars: out.append(j)

        if out != []:
            if api==1: array = Internal.convertDataNodes2Array(out, dim, connects, loc)
            elif api==2: array = Internal.convertDataNodes2Array2(out, dim, connects, loc)
            elif api==3: array = Internal.convertDataNodes2Array3(out, dim, connects, loc)
            else: raise ValueError('getFields: unknow api.')
        else: array = []
        arrays.append(array)
    return arrays

# -- getAllFields
# Retourne les arrays correspondant a une localisation
# si loc='nodes', retourne les champs de
# __GridCoordinates__ + __FlowSolutionNodes__
# si loc='centers', retourne les champs de __FlowSolutionCenters__
# OUT: peut contenir des arrays vides ([])
def getAllFields(t, loc, api=1):
    zones = Internal.getZones(t)
    result = []
    for z in zones:
        if loc == 'nodes':
            f = getFields([Internal.__GridCoordinates__,Internal.__FlowSolutionNodes__], z, api=api)[0]
        else:
            f = getFields(Internal.__FlowSolutionCenters__, z, api=api)[0]
        result.append(f)
    return result

def filterPartialFields(a, arrays, listIndices, loc='nodes', startFrom=0, filterName='', verbose=True):
    ap = Internal.copyRef(a)
    _filterPartialFields(ap, arrays, listIndices, loc=loc, startFrom=startFrom, filterName=filterName, verbose=verbose)
    return ap

def _filterPartialFields(a, arrays, listIndices, loc='nodes', startFrom=0, filterName='', verbose=True):
    if loc == 'nodes': locI=0
    elif loc == 'centers': locI=1
    else: raise ValueError("_filterPartialFields: value of loc argument is not valid.")

    typeN = Internal.typeOfNode(a)
    if typeN != 1: raise ValueError("_filterPartialFields: 1st arg must be a zone.")

    if filterName == '': raise ValueError("_filterPartialFields: filter name must be provided.")
    iverbose = int(verbose)
    Converter.converter.filterPartialFields(a, arrays, listIndices, locI, startFrom, filterName,
                                            Internal.__GridCoordinates__,
                                            Internal.__FlowSolutionNodes__,
                                            Internal.__FlowSolutionCenters__, iverbose)
    return None

# -- setPartialFields
# Set field values in zones at given indices
# IN: t: pyTree a modifier
# IN: arrays: liste des arrays des champs a modifier.
# Chaque array correspond a une zone de t a modifier. Il doit
# contenir le meme nbre de champs que chaque zone de t.
# IN: listIndices: liste des numpys des indices a modifier pour chaque zone
# sous forme de numpys d'entiers.
# NB: un array de la liste arrays peut etre vide (=[])
def setPartialFields(t, arrays, listIndices, loc='nodes', startFrom=0):
    """Set some field values for given indices."""
    tp = Internal.copyRef(t)
    _setPartialFields(tp, arrays, listIndices, loc, startFrom)
    return tp

def _setPartialFields(t, arrays, listIndices, loc='nodes', startFrom=0):
    if loc == 'nodes': locI = 0
    else: locI = 1
    nodes = Internal.getZones(t)
    nzones = len(nodes)
    for c in range(nzones):
        a = arrays[c]; indices = listIndices[c]
        z = nodes[c] # zone
        Converter.converter.setPartialFieldsPT(z, a, indices, locI,
                                               Internal.__GridCoordinates__,
                                               Internal.__FlowSolutionNodes__,
                                               Internal.__FlowSolutionCenters__,
                                               startFrom)
    return None

# --setPartialFields1
# Set field values in zones at given indices
# IN: t: pyTree a modifier
# IN: listFields: liste pour chaque zone de t des numpys des champs.
# Chaque zone  est representee par une liste de numpys (1 par champ).
# IN: listIndices: liste des numpys des indices a modifier pour chaque zone
# sous forme de numpys d'entiers.
# NB: un numpy de la liste des numpys peut etre vide (=[])
def setPartialFields1(t, listFields, listIndices, loc='nodes', startFrom=0):
    tp = Internal.copyRef(t)
    _setPartialFields1(tp, listFields, listIndices, loc,startFrom)
    return tp

def _setPartialFields1(t, listFields, listIndices, loc='nodes', startFrom=0):
    if loc == 'nodes': locI = 0
    else: locI = 1
    zones = Internal.getZones(t)
    for c, z in enumerate(zones):
        Converter.converter._setPartialFields(z, listFields[c], listIndices[c], locI,
                                              Internal.__GridCoordinates__,
                                              Internal.__FlowSolutionNodes__,
                                              Internal.__FlowSolutionCenters__,
                                              startFrom)
    return None

def _updatePartialFields(t, arrays, listIndices, loc='nodes', startFrom=0):
    if loc == 'nodes': locI = 0
    else: locI = 1
    zones = Internal.getZones(t)
    for c, z in enumerate(zones):
        a = arrays[c]; indices = listIndices[c]
        Converter.converter.updatePartialFieldsPT(z, a, indices, locI,
                                                  Internal.__GridCoordinates__,
                                                  Internal.__FlowSolutionNodes__,
                                                  Internal.__FlowSolutionCenters__,
                                                  startFrom)
    return None

# -- setFields
# Set fields a partir de n'importe quel noeud de l'arbre et
# d'une liste d'arrays dans les conteneurs standards:
# si loc=nodes, conteneurs=__GridCoordinates__, __FlowSolutionNodes__
# si loc=centers, conteneurs=__FlowSolutionCenters__
# IN: arrays: liste d'array qui seront mis dans les noeuds zones
# IN: doit contenir des arrays vides si certains champs sont vides.
# IN: t: noeud designant une base, un arbre ou une zone
# IN: loc='nodes' ou 'centers': localisation des arrays
# IN: writeDim: True: force l'ecriture de la dimension/connect du noeud zone
def setFields(arrays, t, loc, writeDim=True):
    # Recherche de tous les noeuds zones a partir de t
    nodes = Internal.getZones(t)

    # Verification de la coherence
    if not isinstance(arrays[0], list):
        raise TypeError("setFields: arrays should be a list of array.")
    if len(nodes) != len(arrays):
        raise ValueError("setFields: more zones in tree than in arrays.")

    for c, a in enumerate(arrays):
        z = nodes[c] # zone
        info = z[2]
        if writeDim and loc == 'nodes' and a != []:
            d = Internal.array2PyTreeDim(a)
            cellDim = d.shape[0]
        else:
            cellDim = z[1].shape[0]

        # Remplace les noeuds contenant les variables
        if a == []: vars = []
        else:
            vars = a[0].split(",")
            # un array * ne peut pas etre mis en nodes
            if loc == 'nodes' and len(a) == 4:
                elt = a[3]
                if elt[len(elt)-1] == '*':
                    print('Warning: setFields: %s array is not set.'%elt)
                    vars = []
        p = 0
        for v in vars:
            renamed = 0 # si le nom est change en nom CGNS = 1
            if v in Internal.name2CGNS:
                variable = Internal.name2CGNS[v]
                if variable != v: renamed = 1
            else: variable = v
            if (variable == 'CoordinateX' or variable == 'CoordinateY'
                or variable == 'CoordinateZ') and loc == 'nodes':
                coordNode = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
                if coordNode is None:
                    info = [Internal.__GridCoordinates__, None, [], 'GridCoordinates_t']
                    z[2].append(info)
                else: info = coordNode
                l = Internal.getNodesFromName(info, variable)
                if l != []:
                    l[0][1] = Internal.createDataNode(variable, a, p, cellDim)[1]
                else:
                    node = Internal.createDataNode(variable, a, p, cellDim)
                    info[2].append(node)
                    if renamed == 1: Internal._rmNodesByName(info, v)

            else: # FlowSolution
                if loc == 'nodes':
                    flowNode = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
                else:
                    flowNode = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
                if flowNode is None:
                    if loc == 'nodes':
                        info = [Internal.__FlowSolutionNodes__, None, [], 'FlowSolution_t']
                    else:
                        info = [Internal.__FlowSolutionCenters__, None, [], 'FlowSolution_t']
                        Internal.createChild(info, 'GridLocation', 'GridLocation_t', value='CellCenter')
                    z[2].append(info)
                else:
                    info = flowNode
                l = Internal.getNodesFromName(info, variable)
                if loc == 'nodes':
                    node = Internal.createDataNode(variable, a, p, cellDim)
                else:
                    node = Internal.createDataNode(variable, a, p, cellDim)
                if l != []: l[0][1] = node[1]
                else:
                    info[2].append(node)
                    if renamed == 1: Internal._rmNodesByName(info, v)

            p += 1

        # update les dimensions si necessaire
        if writeDim and loc == 'nodes' and vars != []:
            z[1] = Internal.array2PyTreeDim(a)
            if len(a) == 5: # structure
                typeNodes = Internal.getNodesFromType1(z, 'ZoneType_t')
                val = Internal.getValue(typeNodes[0])
                if val != 'Structured':
                    # Supprimer GridElementsNode
                    GENodes = Internal.getElementNodes(z)
                    for n in GENodes:
                        p, r = Internal.getParentOfNode(z, n)
                        del p[2][r]
                    Internal._setValue(typeNodes[0], 'Structured')

            elif len(a) == 4: # non structure
                typeNodes = Internal.getNodesFromType1(z, 'ZoneType_t')
                val = Internal.getValue(typeNodes[0])
                if val == 'Structured':
                    Internal._setValue(typeNodes[0], 'Unstructured')
                    if isinstance(a[2], list): # Array2/Array3
                        Internal.setElementConnectivity2(z, a)
                    else: # Array1
                        Internal.setElementConnectivity(z, a)
                else:
                    if isinstance(a[2], list): # Array2/Array3
                        Internal.setElementConnectivity2(z, a)
                    else: # Array1
                        Internal.setElementConnectivity(z, a)
    return t

# -- getNumpyArrays
# Retourne une reference sur le tableau numpy correspondant a un champ
# IN: t: noeud zone, base, tree
# IN: name: nom du champ a extraire (pas de conteneur)
# OUT: liste des numpy arrays (references)
def getNumpyArrays(t, name):
    zones = Internal.getZones(t)
    result = []
    for z in zones:
        sol = Internal.getNodeFromName2(z, name)
        if sol is not None: result.append(sol[1])
        else: result.append([])
    return result

# -- setNumpyArrays
# Remplace une reference sur le tableau numpy correspondant a un champ
# IN: t: noeud zone, base, tree
# IN: name: nom du champ a extraire
# IN: liste des numpy arrays (references)
# OUT: t modifie
def setNumpyArrays(t, name, arrays):
    zones = Internal.getZones(t)
    if len(arrays) != len(zones):
        print('Error: setNumpyArrays: not enough arrays.'); return
    i = 0
    for z in zones:
        sol = Internal.getNodeFromName2(z, name)
        if sol is not None: sol[1] = arrays[i]
        i += 1
    return

# -- ownNumpyArrays
# For numpys allocated with KCore empty (OWNDATA=False), reown it
def _ownNumpyArrays(t):
    if t == []: return None
    n = t[1]
    if isinstance(n, numpy.ndarray) and not n.flags['OWNDATA']:
        b = numpy.copy(n); t[1] = b
    for i in t[2]: _ownNumpyArrays(i)
    return None

# -- convertPyTree2Array
# Convert a python tree node value to an array if possible
def convertPyTree2Array(path, tree):
    """Convert a Python tree node to an array.
    Usage: convertPyTree2Array(path, tree)"""
    p = path.split('/')
    zone = p[0]+'/'+p[1]
    z = Internal.getNodeFromPath(tree, zone)
    if z is None:
        raise ValueError("convertPyTree2Array: zone %s not found."%zone)

    dim = Internal.getZoneDim(z)
    connects = Internal.getElementNodes(z)

    a = Internal.getNodeFromPath(tree, path)

    if a is None:
        raise ValueError("convertPyTree2Array: path %s not found."%path)

    if a[3] == 'DataArray_t':
        array = Internal.convertDataNode2Array(a, dim, connects)[1]
        return array

    elif a[3] == 'IndexRange_t':
        b = numpy.zeros((1,6))
        b.flat = a[1].flat
        return ['index', b, 6, 1, 1]

    elif a[3] == 'GridCoordinates_t':
        info = a[2]; out = []
        for i in info:
            if i[3] == 'DataArray_t':
                a = Internal.convertDataNode2Array(i, dim, connects, 0)[1]
                if a is not None: out.append(a)
        if out != []:
            array = Converter.addVars(out); return array
        else: return None
    elif a[3] == 'FlowSolution_t':
        loc = 0
        locf = Internal.getNodeFromType1(a, 'GridLocation_t')
        if locf is not None and Internal.getValue(locf) == 'CellCenter': loc = 1
        info = a[2]; out = []
        for i in info:
            if i[3] == 'DataArray_t':
                a = Internal.convertDataNode2Array(i, dim, connects, loc)[1]
                if a is not None: out.append(a)
        if out != []:
            array = Converter.addVars(out); return array
        else: return None

    elif a[3] == 'Zone_t':
        if dim[0] == 'Structured':
            np = dim[1]*dim[2]*dim[3]
        else:
            np = dim[1]
        info = a[2]; out = []; out2 = []
        for i in info:
            if i[3] == 'GridCoordinates_t':
                info2 = i[2]
                for j in info2:
                    if j[3] == 'DataArray_t':
                        a = Internal.convertDataNode2Array(j, dim, connects, 0)[1]
                        if a is not None:
                            if a[1].shape[1] == np: out.append(a)
                            else: out2.append(a)

            if i[3] == 'FlowSolution_t':
                loc = 0
                locf = Internal.getNodeFromType1(i, 'GridLocation_t')
                if locf is not None and Internal.getValue(locf) == 'CellCenter': loc = 1
                info2 = i[2]
                for j in info2:
                    if j[3] == 'DataArray_t':
                        a = Internal.convertDataNode2Array(j, dim, connects, loc)[1]
                        if a is not None:
                            if a[1].shape[1] == np: out.append(a)
                            else: out2.append(a)

        if out != [] and out2 != []:
            array = Converter.addVars(out)
            array2 = Converter.addVars(out2)
            print("Warning: convertPyTree2Array: only node field are in array.")
            return array # return array, array2
        elif out != []:
            array = Converter.addVars(out); return array
        elif out2 != []:
            array2 = Converter.addVars(out2); return array2
        else: return None

    else:
        raise ValueError("convertPyTree2Array: in given path, no data node was found.")

# -- makeAllNumpy
# Remplace tous les noeuds 'valeur' (int ou float) par des noeuds numpy de
# une valeur
def makeAllNumpy(t):
    """Make a pyTree all numpy compliant.
    Usage: makeAllNumpy(a)"""
    a = Internal.copyRef(t)
    _makeAllNumpy(a)
    return a

def _makeAllNumpy(a):
    if Internal.isTopTree(a):
        for c in a[1:]: makeAllNumpy__(c)
    else: makeAllNumpy__(a)
    return None

def makeAllNumpy__(node):
    for c in node[2]:
        if isinstance(c[1], str):
            c[1] = numpy.array(c[1], dtype='c')
        if isinstance(c[1], int):
            c[1] = numpy.array(c[1], dtype=Internal.E_NpyInt)
        if isinstance(c[1], float):
            c[1] = numpy.array(c[1], dtype=numpy.float64)
        makeAllNumpy__(c)

# -- makeNumpyStringString
# Remplace tous les noeuds numpy strings (tableau de lettres) par
# des noeuds string
def makeNumpyStringString(t):
    """Make a pyTree all numpy string to strings.
    Usage: makeNumpyStringString(a)"""
    a = Internal.copyRef(t)
    if Internal.isTopTree(a):
        for c in a[1:]: makeNumpyStringString__(c)
    else: makeNumpyStringString__(a)
    return a

def makeNumpyStringString__(node):
    for c in node[2]:
        c[1] = getValue(c)
        makeNumpyStringString__(c)

#==============================================================================
# -- Traitements generiques --
#==============================================================================

# -- TZGC

# Recupere les coords, applique _F dessus sans changement topologique
def __TZGCX(api, t, _F, *args):
    zones = Internal.getZones(t)
    for z in zones:
        fc = getFields(Internal.__GridCoordinates__, z, api=api)[0]
        if fc != []: _F(fc, *args)
    return None

def __TZGC1(t, _F, *args):
    return __TZGCX(1, t, _F, *args)
def __TZGC2(t, _F, *args):
    return __TZGCX(2, t, _F, *args)
def __TZGC3(t, _F, *args):
    return __TZGCX(3, t, _F, *args)

# Recupere les coords, applique F qui renvoie une copie, fait un setFields a locout
def _TZGCX(api, t, locout, writeDim, F, *args):
    zones = Internal.getZones(t)
    for z in zones:
        fc = getFields(Internal.__GridCoordinates__, z, api=api)[0]
        if fc != []:
            fcp = F(fc, *args) # copy
            setFields([fcp], z, locout, writeDim)
    return None

def _TZGC1(t, locout, writeDim, F, *args):
    return _TZGCX(1, t, locout, writeDim, F, *args)
def _TZGC2(t, locout, writeDim, F, *args):
    return _TZGCX(2, t, F, locout, writeDim, F, *args)
def _TZGC3(t, locout, writeDim, F, *args):
    return _TZGCX(3, t, F, locout, writeDim, F, *args)

def TZGCX(api, t, locout, writeDim, F, *args):
    tp = Internal.copyRef(t)
    _TZGCX(api, tp, locout, writeDim, F, *args)
    return tp

def TZGC1(t, locout, writeDim, F, *args):
    return TZGCX(1, t, locout, writeDim, F, *args)
def TZGC2(t, locout, writeDim, F, *args):
    return TZGCX(2, t, locout, writeDim, F, *args)
def TZGC3(t, locout, writeDim, F, *args):
    return TZGCX(3, t, locout, writeDim, F, *args)


# -- TZGF
# Traitement agissant sur le conteneur fieldName, effectue par zones.
# Prend les champs definis dans le conteneur fieldName, applique F, met
# le resultat dans un conteneur suivant locout.
def TZGF(t, locout, fieldName, F, *args):
    tp = Internal.copyRef(t)
    _TZGF(tp, locout, fieldName, F, *args)
    return tp

def _TZGF(t, locout, fieldName, F, *args):
    zones = Internal.getZones(t)
    for z in zones:
        fa = getFields(fieldName, z)[0]
        if fa != []:
            if F is not None: fa = F(fa, *args)
            setFields([fa], z, locout)
    return None

# -- TZA
# Traitement requerant les coord + tous les champs, effectue par zones.
# Prend les champs situes dans locin, effectue le traitement F sur ces champs,
# Remet le resultat dans locout.
# IN: locin: nodes, centers, both
# IN: locout: nodes, centers, both
# sont permis nodes / nodes, nodes / centers, centers / nodes,
# centers / centers, both / both (dans ce cas => nodes/nodes + centers/centers)
# IN: F: fonction de traitement pour les noeuds (possible None)
# IN: Fc: fonction de traitement pour les centres (possible None)
# IN: args: argument de la fonction F
# Si locin=both, args doit contenir les arguments de F pour les noeuds
# suivis des arguments de Fc pour les centres. Les champs en centres ne
# contiennent pas les coordonnees. La fonction F est appelee une
# fois pour les noeuds, Fc une fois pour les centres.
def TZA(t, locin, locout, F, Fc, *args):
    tp = Internal.copyRef(t)
    _TZA(tp, locin, locout, F, Fc, *args)
    return tp

# obsolete : use _TZA1, _TZA3 instead
def _TZA(t, locin, locout, F, Fc, *args):
    zones = Internal.getZones(t)
    for z in zones:
        if locin == 'nodes':
            fc = getFields(Internal.__GridCoordinates__, z)[0]
            fa = getFields(Internal.__FlowSolutionNodes__, z)[0]
            if fc != [] and fa != []:
                Converter._addVars([fc, fa]) # modifie fc
                fp = F(fc, *args)
                setFields([fp], z, locout)
            elif fa != []:
                fp = F(fa, *args)
                setFields([fp], z, locout)
            elif fc != []:
                fp = F(fc, *args)
                setFields([fp], z, locout)
        elif locin == 'centers':
            fa = getFields(Internal.__FlowSolutionCenters__, z)[0]
            if fa != []:
                fp = Fc(fa, *args)
                setFields([fp], z, locout)
        else: # both
            # Dans ce cas, on suppose que F ne change pas la localisation
            l = len(args)//2
            args1 = args[0:l]; args2 = args[l:]
            fc = getFields(Internal.__GridCoordinates__, z)[0]
            fa = getFields(Internal.__FlowSolutionNodes__, z)[0]
            fb = getFields(Internal.__FlowSolutionCenters__, z)[0]
            if fc != [] and fa != []:
                Converter._addVars([fc, fa]) # modifie fc
                fp = F(fc, *args1)
                setFields([fp], z, 'nodes')
            elif fa != []:
                fp = F(fa, *args1)
                setFields([fp], z, 'nodes')
            elif fc != []:
                fp = F(fc, *args1)
                setFields([fp], z, 'nodes')
            fa = None
            if fb != []:
                if Fc is not None: fb = Fc(fb, *args2)
                setFields([fb], z, 'centers')
    return None

# obsolete : use _TZAGC1, _TZAGC3 instead
def _TZAGC(t, locin, locout, writeDim, F, Fc, *args):
    zones = Internal.getZones(t)
    for z in zones:
        if locin == 'nodes':
            fc = getFields(Internal.__GridCoordinates__, z)[0]
            fa = getFields(Internal.__FlowSolutionNodes__, z)[0]
            if fc != [] and fa != []:
                Converter._addVars([fc, fa]) # modifie fc
                fp = F(fc, *args)
                setFields([fp], z, locout, writeDim)
            elif fa != []:
                fp = F(fa, *args)
                setFields([fp], z, locout, writeDim)
            elif fc != []:
                fp = Fc(fc, *args)
                setFields([fp], z, locout, writeDim)
        elif locin == 'centers':
            fc = getFields(Internal.__GridCoordinates__, z)[0]
            fa = getFields(Internal.__FlowSolutionCenters__, z)[0]
            if fc != [] and fa != []:
                posx = KCore.isNamePresent(fa,'CoordinateX')
                posy = KCore.isNamePresent(fa,'CoordinateY')
                posz = KCore.isNamePresent(fa,'CoordinateZ')
                if posx == -1 or posy == -1 or posz == -1:
                    f = Converter.node2Center(fc)
                    Converter._addVars([f, fa])
                else: f = fa

                fp = Fc(f, *args)
                st = fp[0].split(',')
                vars = []
                for i in st:
                    if i != 'CoordinateX' and i != 'CoordinateY' and i != 'CoordinateZ':
                        vars.append(i)
                if vars != []:
                    fp = Converter.extractVars(fp, vars)
                    setFields([fp], z, locout, writeDim)
            elif fc != []:
                fc2 = Converter.node2Center(fc)
                fp = Fc(fc2, *args)
                st = fp[0].split(',')
                vars = []
                for i in st:
                    if i != 'CoordinateX' and i != 'CoordinateY' and i != 'CoordinateZ':
                        vars.append(i)
                if vars != []:
                    fp = Converter.extractVars(fp, vars)
                    setFields([fp], z, locout, writeDim)
            elif fa != []:
                fp = Fc(fa, *args)
                setFields([fp], z, locout)
        else: # both
            l = len(args)//2
            args1 = args[0:l]; args2 = args[l:]
            fc = getFields(Internal.__GridCoordinates__, z)[0]
            fa = getFields(Internal.__FlowSolutionNodes__, z)[0]
            fb = getFields(Internal.__FlowSolutionCenters__, z)[0]
            if fc != [] and fa != []:
                f = Converter.addVars([fc, fa])
                fp = F(f, *args1)
                setFields([fp], z, 'nodes', writeDim)
            elif fa != []:
                fp = F(fa, *args1)
                setFields([fp], z, 'nodes', writeDim)
            elif fc != []:
                fp = F(fc, *args1)
                setFields([fp], z, 'nodes', writeDim)

            fa = None
            if fc != [] and fb != []:
                f = Converter.node2Center(fc)
                Converter._addVars([f, fb])
                fp = Fc(f, *args2)
                st = fp[0].split(',')
                vars = []
                for i in st:
                    if i != 'CoordinateX' and i != 'CoordinateY' and i != 'CoordinateZ':
                        vars.append(i)
                if vars != []:
                    fp = Converter.extractVars(fp, vars)
                    setFields([fp], z, 'centers', writeDim)
            elif fb != []:
                fp = Fc(fb, *args2)
                setFields([fp], z, 'centers', writeDim)
    return None

def TZAGC(t, locin, locout, writeDim, F, Fc, *args):
    tp = Internal.copyRef(t)
    _TZAGC(tp, locin, locout, writeDim, F, Fc, *args)
    return tp

# -- TZANC
# Traitement effectue pour tous les champs + coord. memes pour les centres.
# Dans ce cas, on passe le champ aux centres en noeuds, on applique F,
# et on repasse le champ en centres
# IN: t: arbre a traiter
# IN: locin: nodes, centers, both
# IN: locout: nodes, centers, both
# IN: F: fonction a appliquer pour les noeuds
# IN: Fc: fonction a appliquer pour les centres
# IN: args: arguments de F pour les noeuds + arguments de Fc pour les centres
# dans le cas both.
# La fonction F est appelee une fois pour les noeuds, Fc une fois pour les
# centres.
def TZANC(t, locin, locout, F, Fc, *args):
    tp = Internal.copyRef(t)
    _TZANC(tp, locin, locout, F, Fc, *args)
    return tp

def _TZANC(t, locin, locout, F, Fc, *args):
    zones = Internal.getZones(t)
    l = len(args)//2
    args1 = args[0:l]; args2 = args[l:]
    for z in zones:
        if locin == 'nodes':
            fc = getFields(Internal.__GridCoordinates__, z)[0]
            fa = getFields(Internal.__FlowSolutionNodes__, z)[0]
            if fc != [] and fa != []:
                f = Converter.addVars([fc, fa])
                fp = F(f, *args)
                setFields([fp], z, locout)
            elif fa != []:
                fp = F(fa, *args)
                setFields([fp], z, locout)
            elif fc != []:
                fp = F(fc, *args)
                setFields([fp], z, locout)
        elif locin == 'centers':
            zp = Internal.copyRef(z)
            _deleteFlowSolutions__(zp, 'nodes')
            zp = center2Node(zp, Internal.__FlowSolutionCenters__)
            fa = getFields(Internal.__FlowSolutionNodes__, zp)[0]
            if fa != []:
                fa = Fc(fa, *args)
                if locout == 'nodes': setFields([fa], z, 'nodes')
                else:
                    zp = node2Center(zp, Internal.__FlowSolutionNodes__)
                    fa = getFields(Internal.__FlowSolutionCenters__, zp)[0]
                    setFields([fa], z, 'centers')

        else: # both
            l = len(args)//2
            args1 = args[0:l]; args2 = args[l:]
            fc = getFields(Internal.__GridCoordinates__, z)[0]
            fa = getFields(Internal.__FlowSolutionNodes__, z)[0]
            zp = Internal.copyRef(z)
            _deleteFlowSolutions__(zp, 'nodes')
            zp = center2Node(zp, Internal.__FlowSolutionCenters__)
            fb = getFields(Internal.__FlowSolutionNodes__, zp)[0]
            if fc != [] and fa != []:
                f = Converter.addVars([fc, fa])
                fp = F(f, *args1)
                setFields([fp], z, 'nodes')
            elif fa != []:
                fp = F(fa, *args1)
                setFields([fp], z, 'nodes')
            elif fc != []:
                fp = F(fc, *args1)
                setFields([fp], z, 'nodes')

            if fb != []:
                if Fc is not None: fb = Fc(fb, *args2)
                setFields([fb], zp, 'nodes')
                zp = node2Center(zp, Internal.__FlowSolutionNodes__)
                fa = getFields(Internal.__FlowSolutionCenters__, zp)[0]
                setFields([fa], z, 'centers')
    return None

# -- TLAGC
# Traitement effectue PAR LOT pour tous les champs + coord. memes pour
# les centres.
# Dans ce cas, on reconstruit un maillage en centres qui sera associe
# au champ en centres.
# IN: t: arbre a traiter
# IN: F: fonction a appliquer
# IN: args: arguments de F pour les noeuds + arguments de F pour les centres
# OUT: le resultat du traitement par lot (entier, reels, arrays)
# Si sortie=arrays de Converter
#  1- c'est a l'appelant de retourner une zone/liste de zones/arbre
#  2- traitement des centres: l'appelant doit enlever les coordonnees des centres
def TLAGC(t, F, *args):
    tp = Internal.copyRef(t)
    zones = Internal.getZones(tp)
    l = len(args)//2
    args1 = args[0:l]; args2 = args[l:]
    allfc = []; allfa = []; allfb = []
    nzones = len(zones)
    if nzones == 0: return tp

    for z in zones:
        fc = getFields(Internal.__GridCoordinates__, z)[0]
        fa = getFields(Internal.__FlowSolutionNodes__, z)[0]
        fb = getFields(Internal.__FlowSolutionCenters__, z)[0]
        allfc.append(fc); allfa.append(fa); allfb.append(fb)

    # Application par lot des noeuds
    allf = []
    for i in range(nzones):
        if allfc[i] != [] and allfa[i] != []:
            allf.append(Converter.addVars([allfc[i], allfa[i]]))
        elif allfa[i] != []: allf.append(allfa[i])
        elif allfc[i] != []: allf.append(allfc[i])
    if allf != []: allf = F(allf, *args1)
    allfa = allf # ref

    # Application par lot des centres
    allf = []
    allfc = Converter.node2Center(allfc)
    for i in range(nzones):
        if allfc[i] != [] and allfb[i] != []:
            allf.append(Converter.addVars([allfc[i],allfb[i]]))
        elif allfb[i] != []: allf.append(allfb[i])

    if allf != []: allf = F(allf, *args2)
    allfb = allf # ref

    if allfa != [] and allfb != []: return [allfa,allfb]
    elif allfa != []: return [allfa]
    elif allfb != []: return [allfb]
    else: return []

# -- __TZC = generique
# Traitement uniquement sur les champs
def __TZCX(api, t, locin, writeDim, _F, *args):
    zones = Internal.getZones(t)
    for z in zones:
        if   locin == 'nodes':
            fa = getFields(Internal.__FlowSolutionNodes__, z, api=api)[0]
        elif locin == 'centers':
            fa = getFields(Internal.__FlowSolutionCenters__, z, api=api)[0]

        if fa != []:
            _F(fa, *args)
            setFields([fa], z, locin, writeDim)
    return None

def __TZC1(t, locin, writeDim, _F, *args):
    return __TZCX(1, t, locin, writeDim, _F, *args)
def __TZC2(t, locin, writeDim, _F, *args):
    return __TZCX(2, t, locin, writeDim, _F, *args)
def __TZC3(t, locin, writeDim, _F, *args):
    return __TZCX(3, t, locin, writeDim, _F, *args)


# -- __TZA = generique
# Recupere les champs de locin en shared
# applique _F (array in place)
def __TZAX(api, t, locin, _F, *args):
    zones = Internal.getZones(t)
    for z in zones:
        if locin == 'nodes':
            fc = getFields(Internal.__GridCoordinates__, z, api=api)[0]
            fa = getFields(Internal.__FlowSolutionNodes__, z, api=api)[0]
            if fc != [] and fa != []:
                if api == 1: Converter._addVars([fc, fa])
                else:
                    fc[0] = fc[0]+','+fa[0]
                    fc[1] = fc[1]+fa[1]
                _F(fc, *args)
            elif fc != []: _F(fc, *args)
            elif fa != []: _F(fa, *args)
        elif locin == 'centers':
            fa = getFields(Internal.__FlowSolutionCenters__, z, api=api)[0]
            if fa != []: _F(fa, *args)
        elif locin == 'both':
            fc = getFields(Internal.__GridCoordinates__, z, api=api)[0]
            fa = getFields(Internal.__FlowSolutionNodes__, z, api=api)[0]
            if fc != [] and fa != []:
                if api == 1: fc = Converter.addVars([fc, fa])
                else:
                    fc[0] = fc[0]+','+fa[0]
                    fc[1] = fc[1]+fa[1]
                _F(fc, *args)
            elif fc != []: _F(fc, *args)
            elif fa != []: _F(fa, *args)
            fa = getFields(Internal.__FlowSolutionCenters__, z, api=api)[0]
            if fa != []: _F(fa, *args)
    return None

def __TZA1(t, locin, _F, *args):
    return __TZAX(1, t, locin, _F, *args)
def __TZA2(t, locin, _F, *args):
    return __TZAX(2, t, locin, _F, *args)
def __TZA3(t, locin, _F, *args):
    return __TZAX(3, t, locin, _F, *args)

# -- _TZA = generique
# Recupere les champs locin en shared
# Applique F qui rend une copie
# Remet cet array1/2/3 dans t a locout
def _TZAX(api, t, locin, locout, writeDim, F, *args):
    zones = Internal.getZones(t)
    for z in zones:
        if locin == 'nodes':
            fc = getFields(Internal.__GridCoordinates__, z, api=api)[0]
            fa = getFields(Internal.__FlowSolutionNodes__, z, api=api)[0]
            ret = None
            if fc != [] and fa != []:
                if api == 1: Converter._addVars([fc, fa])
                else:
                    fc[0] = fc[0]+','+fa[0]
                    fc[1] = fc[1]+fa[1]
                ret = F(fc, *args)
            elif fc != []: ret = F(fc, *args)
            elif fa != []: ret = F(fa, *args)
            if ret is not None: setFields([ret], z, locout, writeDim)
        elif locin == 'centers':
            fa = getFields(Internal.__FlowSolutionCenters__, z, api=api)[0]
            if fa != []:
                ret = F(fa, *args)
                setFields([ret], z, locout, writeDim)
        else: # both
            # Dans ce cas, on suppose que F ne change pas la localisation
            # les nodes restent en nodes, les centres restent en centres par F
            fc = getFields(Internal.__GridCoordinates__, z, api=api)[0]
            fa = getFields(Internal.__FlowSolutionNodes__, z, api=api)[0]
            fb = getFields(Internal.__FlowSolutionCenters__, z, api=api)[0]
            if fc != [] and fa != []:
                if api == 1: Converter._addVars([fc, fa]) # modifie fc
                else:
                    fc[0] = fc[0]+','+fa[0]
                    fc[1] = fc[1]+fa[1]
                fp = F(fc, *args)
                setFields([fp], z, 'nodes', writeDim)
            elif fa != []:
                fp = F(fa, *args)
                setFields([fp], z, 'nodes', writeDim)
            elif fc != []:
                fp = F(fc, *args)
                setFields([fp], z, 'nodes', writeDim)
            fa = None
            if fb != []:
                if locout != 'nodes': fb = F(fb, *args)
                setFields([fb], z, 'centers', writeDim)
    return None

def _TZA1(t, locin, locout, writeDim, F, *args):
    return _TZAX(1, t, locin, locout, writeDim, F, *args)
def _TZA2(t, locin, locout, writeDim, F, *args):
    return _TZAX(2, t, locin, locout, writeDim, F, *args)
def _TZA3(t, locin, locout, writeDim, F, *args):
    return _TZAX(3, t, locin, locout, writeDim, F, *args)

# Fait une ref copie en +
def TZAX(api, t, locin, locout, writeDim, F, *args):
    tp = Internal.copyRef(t)
    _TZAX(api, tp, locin, locout, writeDim, F, *args)
    return tp

def TZA1(t, locin, locout, writeDim, F, *args):
    return TZAX(1, t, locin, locout, writeDim, F, *args)
def TZA2(t, locin, locout, writeDim, F, *args):
    return TZAX(2, t, locin, locout, writeDim, F, *args)
def TZA3(t, locin, locout, writeDim, F, *args):
    return TZAX(3, t, locin, locout, writeDim, F, *args)

# -- TZAGC generique
# Traitement effectue pour tous les champs + coord. memes pour les centres.
# Dans ce cas, on reconstruit un maillage en centres qui sera associe
# aux champs en centres.
# IN: t: arbre a traiter
# IN: locin: nodes, centers, both
# IN: locout: nodes, centers, both
# IN: F: fonction a appliquer pour les noeuds
# IN: Fc: fonction a appliquer pour les centres
# IN: args: arguments de F pour les noeuds + arguments de Fc pour les centres
# dans le cas both.
# La fonction F est appelee une fois pour les noeuds, Fc une fois pour les
# centres.
def _TZAGCX(api, t, locin, locout, writeDim, F, Fc, *args):
    zones = Internal.getZones(t)
    for z in zones:
        if locin == 'nodes':
            fc = getFields(Internal.__GridCoordinates__, z, api)[0]
            fa = getFields(Internal.__FlowSolutionNodes__, z, api)[0]
            if fc != [] and fa != []:
                if api == 1: Converter._addVars([fc, fa])
                else: # array2/3
                    fc[0] = fc[0]+','+fa[0]
                    fc[1] = fc[1]+fa[1]
                fp = F(fc, *args)
                setFields([fp], z, locout, writeDim)
            elif fa != []:
                fp = F(fa, *args)
                setFields([fp], z, locout, writeDim)
            elif fc != []:
                fp = Fc(fc, *args)
                setFields([fp], z, locout, writeDim)
        elif locin == 'centers':
            fc = getFields(Internal.__GridCoordinates__, z, api=api)[0]
            fa = getFields(Internal.__FlowSolutionCenters__, z, api=api)[0]
            if fc != [] and fa != []:
                posx = KCore.isNamePresent(fa, 'CoordinateX')
                posy = KCore.isNamePresent(fa, 'CoordinateY')
                posz = KCore.isNamePresent(fa, 'CoordinateZ')
                if posx == -1 or posy == -1 or posz == -1:
                    f = Converter.node2Center(fc)
                    if api == 1: Converter._addVars([f, fa])
                    else: # array2/3
                        f[0] = f[0]+','+fa[0]
                        f[1] = f[1]+fa[1]
                else: f = fa

                fp = Fc(f, *args)
                st = fp[0].split(',')
                vars = []
                for i in st:
                    if i != 'CoordinateX' and i != 'CoordinateY' and i != 'CoordinateZ':
                        vars.append(i)
                if vars != []:
                    if api == 1: fp2 = Converter.extractVars(fp, vars)
                    else:
                        fp2 = fp[:]; fp2[1] = []; vs = ''
                        for i in vars:
                            vs += i+','
                            p = KCore.isNamePresent(fp, i)
                            fp2[1].append(fp[1][p])
                        fp2[0] = vs[:-1]
                    setFields([fp2], z, locout, writeDim)
            elif fc != []:
                fc2 = Converter.node2Center(fc)
                fp = Fc(fc2, *args)
                st = fp[0].split(',')
                vars = []
                for i in st:
                    if i != 'CoordinateX' and i != 'CoordinateY' and i != 'CoordinateZ':
                        vars.append(i)
                if vars != []:
                    if api == 1: fp2 = Converter.extractVars(fp, vars)
                    else: # array2/3
                        fp2 = fp[:]; fp2[1] = []; vs = ''
                        for i in vars:
                            vs += i+','
                            p = KCore.isNamePresent(fp, i)
                            fp2[1].append(fp[1][p])
                        fp2[0] = vs[:-1]
                    setFields([fp2], z, locout, writeDim)
            elif fa != []:
                fp = Fc(fa, *args)
                setFields([fp], z, locout, writeDim)
        else: # both
            l = len(args)//2
            args1 = args[0:l]; args2 = args[l:]
            fc = getFields(Internal.__GridCoordinates__, z, api=api)[0]
            fa = getFields(Internal.__FlowSolutionNodes__, z, api=api)[0]
            fb = getFields(Internal.__FlowSolutionCenters__, z, api=api)[0]
            if fc != [] and fa != []:
                if api == 1: f = Converter.addVars([fc, fa])
                else:
                    f = fc[:]
                    f[0] = fc[0]+','+fa[0]
                    f[1] = fc[1]+fa[1]
                fp = F(f, *args1)
                setFields([fp], z, 'nodes', writeDim)
            elif fa != []:
                fp = F(fa, *args1)
                setFields([fp], z, 'nodes', writeDim)
            elif fc != []:
                fp = F(fc, *args1)
                setFields([fp], z, 'nodes', writeDim)
            fa = None
            if fc != [] and fb != []:
                f = Converter.node2Center(fc)
                if api == 1: Converter._addVars([f, fb])
                else: # array2/3
                    f[0] = f[0]+','+fb[0]
                    f[1] = f[1]+fb[1]
                fp = Fc(f, *args2)
                st = fp[0].split(',')
                vars = []
                for i in st:
                    if i != 'CoordinateX' and i != 'CoordinateY' and i != 'CoordinateZ':
                        vars.append(i)
                if vars != []:
                    if api == 1: fp2 = Converter.extractVars(fp, vars)
                    else:
                        fp2 = fp[:]; fp2[1] = []; vs = ''
                        for i in vars:
                            vs += i+','
                            p = KCore.isNamePresent(fp, i)
                            fp2[1].append(fp[1][p])
                        fp2[0] = vs[:-1]
                    setFields([fp2], z, 'centers', writeDim)
            elif fb != []:
                fp = Fc(fb, *args2)
                setFields([fp], z, 'centers', writeDim)
    return None

def _TZAGC1(t, locin, locout, writeDim, F, Fc, *args):
    return _TZAGCX(1, t, locin, locout, writeDim, F, Fc, *args)
def _TZAGC2(t, locin, locout, writeDim, F, Fc, *args):
    return _TZAGCX(2, t, locin, locout, writeDim, F, Fc, *args)
def _TZAGC3(t, locin, locout, writeDim, F, Fc, *args):
    return _TZAGCX(3, t, locin, locout, writeDim, F, Fc, *args)

# Fait une ref copie en +
def TZAGCX(api, t, locin, locout, writeDim, F, Fc, *args):
    tp = Internal.copyRef(t)
    _TZAGCX(api, tp, locin, locout, writeDim, F, Fc, *args)
    return tp

def TZAGC1(t, locin, locout, writeDim, F, Fc, *args):
    return TZAGCX(1, t, locin, locout, writeDim, F, Fc, *args)
def TZAGC2(t, locin, locout, writeDim, F, Fc, *args):
    return TZAGCX(2, t, locin, locout, writeDim, F, Fc, *args)
def TZAGC3(t, locin, locout, writeDim, F, Fc, *args):
    return TZAGCX(3, t, locin, locout, writeDim, F, Fc, *args)

#==============================================================================
# -- Fields / Vars management --
#==============================================================================

# -- addVars: ajoute des variables dans un pyTree
def addVars(t, vars):
    """Add variables to a pyTree.
    Usage: a = addVars(t, vars)"""
    tp = Internal.copyRef(t)
    _addVars(tp, vars)
    return tp

def _addVars(t, vars):
    zones = Internal.getZones(t)
    if isinstance(vars, list):
        for var in vars:
            loc = 'nodes'; v = var.split(':',1)
            if len(v) > 1:
                if v[0] == 'nodes' or v[0] == 'centers':
                    var = v[1]; loc = v[0]
            for z in zones:
                found = 0
                if loc == 'centers':
                    node = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
                    if node is not None:
                        nodeofname = Internal.getNodeFromName1(node, var)
                        if nodeofname is not None: found = 1
                if loc == 'nodes':
                    node = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
                    if node is not None:
                        nodeofname = Internal.getNodeFromName1(node, var)
                        if nodeofname is not None: found = 1
                    node = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
                    if node is not None:
                        nodeofname = Internal.getNodeFromName1(node, var)
                        if nodeofname is not None: found = 1
                if found == 0:
                    dim = Internal.getZoneDim(z)
                    _addAVar__(z, dim, '%s:%s'%(loc,var))
    else:
        loc = 'nodes'; v = vars.split(':',1)
        if len(v) > 1:
            if loc == 'nodes' or loc == 'centers': vars = v[1]; loc = v[0]
        for z in zones:
            found = 0
            if loc == 'centers':
                node = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
                if node is not None:
                    nodeofname = Internal.getNodeFromName1(node, vars)
                    if nodeofname is not None: found = 1
            if loc == 'nodes':
                node = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
                if node is not None:
                    nodeofname = Internal.getNodeFromName1(node, vars)
                    if nodeofname is not None: found = 1
                node = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
                if node is not None:
                    nodeofname = Internal.getNodeFromName1(node, vars)
                    if nodeofname is not None: found = 1
            if found == 0:
                dim = Internal.getZoneDim(z)
                _addAVar__(z, dim, '%s:%s'%(loc,vars))
    return None

def _addAVar__(z, dim, var):
    vref = getVarNames(z)
    if var in vref: return # variable deja existante
    loc = 'nodes'
    v = var.split(':',1)
    if len(v) > 1:
        if v[0] == 'nodes' or v[0] == 'centers': var = v[1]; loc = v[0]
    if dim[0] == 'Structured':
        if loc == 'nodes':
            ni = dim[1]; nj = dim[2]; nk = dim[3]
            if var == 'cellN' or var == 'cellNF':
                n = numpy.ones((ni*nj*nk), dtype=numpy.float64)
            else:
                n = numpy.zeros((ni*nj*nk), dtype=numpy.float64)
            array2 = [var, [n], ni, nj, nk]
        elif loc == 'centers':
            ni1 = max(dim[1]-1,1); nj1 = max(dim[2]-1,1); nk1 = max(dim[3]-1,1)
            if var == 'cellN' or var == 'cellNF':
                n = numpy.ones((ni1*nj1*nk1), dtype=numpy.float64)
            else:
                n = numpy.zeros((ni1*nj1*nk1), dtype=numpy.float64)
            array2 = [var, [n], ni1, nj1, nk1]
    else:
        if loc == 'nodes':
            npts = dim[1]
            if var == 'cellN' or var == 'cellNF':
                n = numpy.ones((npts), dtype=numpy.float64)
            else:
                n = numpy.zeros((npts), dtype=numpy.float64)
            array2 = [var, [n], npts, 1, 1]
        else:
            nelts = dim[2]
            if var == 'cellN' or var == 'cellNF':
                n = numpy.ones((nelts), dtype=numpy.float64)
            else:
                n = numpy.zeros((nelts), dtype=numpy.float64)
            array2 = [var, [n], nelts, 1, 1]
    setFields([array2], z, loc=loc, writeDim=False)
    return None

# -- initVars: initialise une variable
def initVars(t, varNameString, v1=[], v2=[], mode=0, isVectorized=False):
    """Init variables defined by varNameString.
    Usage: a = initVars(array, varNameString, val)
    or
    Usage: a = initVars(array, varNameString, F, [strings])"""
    tp = Internal.copyRef(t)
    _initVars(tp, varNameString, v1, v2, mode, isVectorized)
    return tp

def _initVars(t, varNameString, v1=[], v2=[], mode=0, isVectorized=False):
    loc = 'nodes'
    centerCoordNeeded = False
    # Check that the type of varNameString is correct for all types of initialisations
    isFuncOrConstInit = callable(v1) or isinstance(v1, (int, float, numpy.float64))

    if isFuncOrConstInit:
        # Initialisation(s) by constant or function
        varNames = []
        if not isinstance(varNameString, list): varNameString = [varNameString]
        for varName in varNameString:
            v = varName.replace('}', '').replace('{', '').strip().split(':',1)
            if len(v) > 1 and v[0] in ['nodes', 'centers']:
                loc = v[0]
                varNames.append(v[1])
            else:
                varNames.append(v[0])

        if callable(v1) and loc == 'centers':
            # Initialisation(s) by function of (a) cell-centered lhs variable(s)
            for i, farg in enumerate(v2):
                v = farg.replace('}', '').replace('{', '').strip().split(':',1)
                if len(v) > 1:
                    if v[0] == 'centers' and v[1].startswith('Coordinate'): centerCoordNeeded = True
                    v2[i] = v[1]

    else:
        # Initialisation by string
        s = varNameString.split('=')
        varName = s[0].replace('}', '').replace('{', '').strip()
        v = varName.split(':',1)
        if len(v) > 1 and v[0] == 'centers' and v1 == []:
            loc = v[0]
            rhsVars = re.findall("{.*?}", s[1])
            for var in rhsVars:
                v = var.replace('}', '').replace('{', '').strip().split(':',1)
                if len(v) > 1 and v[0] == 'centers' and v[1].startswith('Coordinate'): centerCoordNeeded = True

    #centerCoordNeeded = True # for DBX
    if centerCoordNeeded:
        if callable(v1):
            # Initialisation(s) by function # Reg
            _TZAGC3(t, loc, loc, False, Converter.initVars,
                    Converter.initVars, varNames, v1, v2, mode, isVectorized)
        elif isinstance(v1, (int, float, numpy.float64)):
            # Initialisation(s) by constant
            _TZAGC3(t, loc, loc, False, Converter.initVars, Converter.initVars,
                    varNames, v1, v2, mode)
        else:
            # Initialisation by string
            _TZAGC3(t, loc, loc, False, Converter.initVars, Converter.initVars,
                    varNameString, v1, v2, mode)
    else:
        if v1 == []:
            # Initialisation by string
            _addVars(t, varName)
            __TZA3(t, loc, Converter._initVars, varNameString, v1, v2, mode)
        else:
            # Initialisation(s) ...
            [_addVars(t, varName) for varName in varNameString]
            if callable(v1):
                # ... by function
                __TZA3(t, loc, Converter._initVars, varNames, v1, v2, mode, isVectorized)
            else:
                # ... by constant
                __TZA3(t, loc, Converter._initVars, varNames, v1, v2, mode)
    return None

# Merge BCDataSets
def _mergeBCDataSets(t):
    zones = Internal.getZones(t)
    if zones == []: zones = [t] # must be a BC node
    for z in zones:
        bcs = Internal.getNodesFromType2(z, 'BC_t')
        for b in bcs:
            Internal._mergeBCDataSets__(z,b)
    return None

# Nullify field in BCDataSet
def _nullifyBCDataSetVectors(t, bndType, loc='FaceCenter',
                             vectors=[['VelocityX','VelocityY','VelocityZ']]):
    locI = 0
    if loc == 'FaceCenter' or loc == 'CellCenter':
        locI = 1; FSCont = Internal.__FlowSolutionCenters__
    else:
        locI = 0; FSCont = Internal.__FlowSolutionNodes__

    zones = Internal.getZones(t)
    families = getFamilyBCNamesOfType(t, bndType)
    for z in zones:
        dimZ = Internal.getZoneDim(z)
        niZ = dimZ[1]; njZ = dimZ[2]; nkZ = dimZ[3]
        allbcs = Internal.getNodesFromType2(z, 'BC_t')
        bcs = Internal.getNodesFromValue(allbcs, bndType)
        bcs += getFamilyBCs(z,families)    # CW : plutot allbcs ??
        for bc in bcs:
            PR = Internal.getNodeFromName1(bc, 'PointRange')
            PL = Internal.getNodeFromName1(bc, 'PointList')
            np = 0
            if PR is not None:
                win =  Internal.range2Window(PR[1])
                imin = win[0]; imax = win[1]
                jmin = win[2]; jmax = win[3]
                kmin = win[4]; kmax = win[5]
                if locI == 0:
                    di = max(1,imax-imin+1)
                    dj = max(1,jmax-jmin+1)
                    dk = max(1,kmax-kmin+1)
                else:
                    di = max(1,imax-imin)
                    dj = max(1,jmax-jmin)
                    dk = max(1,kmax-kmin)
                np = di*dj*dk
            elif PL is not None:
                np = PL[1].size
            else:
                raise ValueError("nullifyVectorAtBCDataSet: no PointRange/PointList in BC.")

            datas = Internal.getBCDataSet(z,bc)
            if datas == []: # create the BCDataSet
                d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                                          gridLocation=loc, parent=bc)
                d = Internal.newBCData('BCDirichlet', parent=d)
                cont, noc = Internal.getParentOfNode(z, d)
                for vect in vectors:
                    vxname = vect[0]; fxInt = numpy.zeros((np),numpy.float64)
                    vyname = vect[1]; fyInt = numpy.zeros((np),numpy.float64)
                    vzname = vect[2]; fzInt = numpy.zeros((np),numpy.float64)
                    if PR is not None:
                        Converter.converter.nullifyVectorAtBCFaceStruct(z, fxInt, fyInt, fzInt,
                                                                        imin, imax, jmin, jmax, kmin, kmax, locI,
                                                                        vxname, vyname, vzname,
                                                                        Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__,Internal.__FlowSolutionCenters__)

                    elif PL is not None:
                        print("Warning: nullifyVectorAtBCDataSet: not implemented for PointList.")
                    Internal._createUniqueChild(d, vxname, 'DataArray_t', value=fxInt)
                    Internal._createUniqueChild(d, vyname, 'DataArray_t', value=fyInt)
                    Internal._createUniqueChild(d, vzname, 'DataArray_t', value=fzInt)

            else: # BCDataSet exists: add missing variables
                d = Internal.getNodeFromType2(bc,'BCData_t')
                for vect in vectors:
                    vxname = vect[0]; fxInt = numpy.zeros((np),numpy.float64)
                    vyname = vect[1]; fyInt = numpy.zeros((np),numpy.float64)
                    vzname = vect[2]; fzInt = numpy.zeros((np),numpy.float64)
                    Converter.converter.nullifyVectorAtBCFaceStruct(z, fxInt, fyInt, fzInt,
                                                                    imin, imax, jmin, jmax, kmin, kmax, locI,
                                                                    vxname, vyname, vzname,
                                                                    Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__,Internal.__FlowSolutionCenters__)
                    Internal._createUniqueChild(d, vxname, 'DataArray_t', value=fxInt)
                    Internal._createUniqueChild(d, vyname, 'DataArray_t', value=fyInt)
                    Internal._createUniqueChild(d, vzname, 'DataArray_t', value=fzInt)
    return None

# Create a BCDataSet for a bndType by Oth-order extrapolation from interior
# loc='FaceCenter' or 'Vertex'
# update = True: update the BCDataSet by extrapolation if field already exists
def _createBCDataSetOfType(t, bndType, loc='FaceCenter', update=True, vectors=[]):
    locI = 0
    if loc == 'FaceCenter': locI = 1; FSCont = Internal.__FlowSolutionCenters__
    else: locI = 0; FSCont = Internal.__FlowSolutionNodes__

    zones = Internal.getZones(t)
    families = getFamilyBCNamesOfType(t, bndType)
    for z in zones:
        dimZ = Internal.getZoneDim(z)
        niZ = dimZ[1]; njZ = dimZ[2]; nkZ = dimZ[3]
        allbcs = Internal.getNodesFromType2(z, 'BC_t')
        bcs = Internal.getNodesFromValue(allbcs, bndType)
        bcs += getFamilyBCs(z,families)
        FSNode = Internal.getNodeFromName1(z, FSCont)
        varnames=[]
        for fs in FSNode[2]:
            if fs[3]=='DataArray_t': varnames.append(fs[0])
        for bc in bcs:
            PR = Internal.getNodeFromName1(bc, 'PointRange')
            PL = Internal.getNodeFromName1(bc, 'PointList')
            np = 0
            if PR is not None:
                win =  Internal.range2Window(PR[1])
                imin = win[0]; imax = win[1]
                jmin = win[2]; jmax = win[3]
                kmin = win[4]; kmax = win[5]
                if locI == 0:
                    di = max(1,imax-imin+1)
                    dj = max(1,jmax-jmin+1)
                    dk = max(1,kmax-kmin+1)
                else:
                    di = max(1,imax-imin)
                    dj = max(1,jmax-jmin)
                    dk = max(1,kmax-kmin)
                np = di*dj*dk
            elif PL is not None:
                np = PL[1].size
            else:
                raise ValueError("createBCDataSetOfType: no PointRange/PointList in BC.")

            datas = Internal.getBCDataSet(z,bc)
            if datas == []: # create the BCDataSet
                d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                                          gridLocation=loc, parent=bc)
                d = Internal.newBCData('BCDirichlet', parent=d)
                cont, noc = Internal.getParentOfNode(z, d)
                for fs in FSNode[2]:
                    if fs[3]=='DataArray_t':
                        varname = fs[0]
                        fInt = numpy.zeros((np),numpy.float64)
                        if PR is not None:
                            Converter.converter.extrapInterior2BCFaceStruct(z, fInt, imin, imax, jmin, jmax, kmin, kmax, locI, varname,
                                                                            Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__,Internal.__FlowSolutionCenters__)

                        elif PL is not None:
                            print("Warning: createBCDataSetOfType: extrapolation from interior cells not implemented for PointList. Fields are initialized by O in BCDataSet.")
                        Internal._createUniqueChild(d, varname, 'DataArray_t', value=fInt)

            else: # BCDataSet exists: add missing variables
                d = Internal.getNodeFromType2(bc,'BCData_t')
                dataSetNames = []
                for dataSet in d[2]: dataSetNames.append(dataSet[0])
                for varname in varnames:
                    if varname not in dataSetNames or update is True:
                        fInt = numpy.zeros((np),numpy.float64)
                        Converter.converter.extrapInterior2BCFaceStruct(z, fInt, imin, imax, jmin, jmax, kmin, kmax, locI, varname,
                                                                        Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__,Internal.__FlowSolutionCenters__)
                        Internal._createUniqueChild(d, varname, 'DataArray_t', value=fInt)
    if bndType == 'BCWallInviscid':
        _nullifyBCDataSetVectors(t, bndType, loc=loc, vectors=[['VelocityX','VelocityY','VelocityZ']])
    elif bndType == 'BCSymmetryPlane':
        if vectors == []:
            _nullifyBCDataSetVectors(t, bndType, loc=loc, vectors=[['VelocityX','VelocityY','VelocityZ']])
        else:
            _nullifyBCDataSetVectors(t, bndType, loc=loc, vectors=[['VelocityX','VelocityY','VelocityZ']]+vectors)
    return None

# Apply init to all BCDataSets
def _initBCDataSet(t, varNameString, val1=[], val2=[]):
    zones = Internal.getZones(t)
    if zones == []: zones = [t] # must be a BC node
    for z in zones:
        bcs = Internal.getNodesFromType2(z, 'BC_t')
        for b in bcs:
            datas = Internal.getBCDataSet(z, b)

            fields = []; connects = []
            for d in datas: # build array
                np = d[1].size; ne = np-1
                dim = ['Unstructured', np, ne, 'NODE', 1]
                f = Internal.convertDataNode2Array(d, dim, connects)
                fields.append(f[1])
            if fields != []:
                fields = Converter.addVars(fields)
            else:
                np1 = Internal.getNodeFromName1(b, 'PointList')
                np2 = Internal.getNodeFromName1(b, 'PointRange')
                if np2 is not None:
                    win = Internal.range2Window(np2[1])
                    imin = win[0]; imax = win[1]
                    jmin = win[2]; jmax = win[3]
                    kmin = win[4]; kmax = win[5]
                    np = max(imax-imin,1)*max(jmax-jmin,1)*max(kmax-kmin,1)
                elif np1 is not None:
                    np = np1[1].size
                else: raise ValueError('initBCDataSet: no PointRange or PointList in BC.')
                fields = Converter.array('empty',np,1,1)
            fn = Converter.initVars(fields, varNameString, val1, val2)
            nofld = fn[1].shape[0]-1
            varName = varNameString.split('=')[0]
            varName = varName.replace('{', '')
            varName = varName.replace('}', '')
            varName = varName.replace('centers:', '')
            varName = varName.replace('nodes:', '')
            f = Converter.extractVars(fn, [varName])
            fieldFaceNode = Internal.createDataNode(varName, f, 0, cellDim=1)
            if datas != []:
                cont, c = Internal.getParentOfNode(z, datas[0])
                Internal._createUniqueChild(cont, varName, 'DataArray_t', value=fieldFaceNode[1])
            else:
                d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                                          gridLocation='FaceCenter', parent=b)
                d = Internal.newBCData('BCData', parent=d)
                Internal._createUniqueChild(d, varName, 'DataArray_t', value=fieldFaceNode[1])
    return None

# -- randomizeVar: ajoute du bruit sur une variable
def randomizeVar(t, var, deltaMin, deltaMax):
    """Randomize a field defined by var within a range [a-deltamin,a+deltamax]
    Usage: a = randomizeVar(a, var, deltaMin, deltaMax)"""
    tp = Internal.copyRef(t)
    _randomizeVar(tp, var, deltaMin, deltaMax)
    return tp

def _randomizeVar(t, var, deltaMin, deltaMax):
    loc = 'nodes'
    varname = var
    spl = var.split(':')
    if len(spl) != 1:
        if spl[0] == 'centers': loc = 'centers'
        varname = spl[1]
    _TZA(t, loc, loc, Converter.randomizeVar, Converter.randomizeVar, varname, deltaMin, deltaMax)
    return None

def _orderVariables(t, varsn=[], varsc=[]):
    nodes = Internal.getZones(t)
    if varsn==[]:
        for z in nodes:
            vars = getVarNames(z, excludeXYZ=True, loc='nodes')[0]
            for var in vars:
                found = 0
                for v in varsn:
                    if v == var: found = 1
                if found == 0: varsn.append(var)
    if varsc == []:
        for z in nodes:
            vars = getVarNames(z, excludeXYZ=True, loc='centers')[0]
            for var in vars:
                found = 0
                for v in varsc:
                    if v == var: found = 1
                if found == 0: varsc.append(var)

    #1/ Sort coordinates
    #2/ Sort FlowSolution nodes
    varx = ['CoordinateX', 'CoordinateY', 'CoordinateZ']
    for z in nodes:
        cont = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
        if cont is not None:
            children = cont[2]
            childrenNames = []
            for c in children: childrenNames.append(c[0])
            new = []
            for name in varx:
                try:
                    s1 = childrenNames.index(name)
                    new.append(children[s1])
                except: pass
            cont[2] = new

    for z in nodes:
        cont = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
        if cont is not None:
            children = cont[2]
            childrenNames = []
            for c in children: childrenNames.append(c[0])
            new = []
            for name in varsn:
                try:
                    s1 = childrenNames.index(name)
                    new.append(children[s1])
                except: pass
            pos = 0
            for name in childrenNames:
                try:
                    s1 = varsn.index(name)
                except: new.append(children[pos])
                pos += 1
            cont[2] = new

    #3. Sort flow solution centers
    for z in nodes:
        cont = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
        if cont is not None:
            children = cont[2]
            childrenNames = []
            for c in children: childrenNames.append('centers:'+c[0])
            new = []
            for name in varsc:
                try:
                    s1 = childrenNames.index(name)
                    new.append(children[s1])
                except: pass
            pos = 0
            for name in childrenNames:
                try:
                    s1 = varsc.index(name)
                except: new.append(children[pos])
                pos += 1
            cont[2] = new
    return None

#reorder containers: GC, then FSN, then FSC
def _orderContainers(t):
    for z in Internal.getZones(t):
        c = 0; i0 = -1; i1 = -1; i2 = -1
        for j in z[2]:
            if j[0] == Internal.__GridCoordinates__: i0 = c
            if j[0] == Internal.__FlowSolutionNodes__: i1 = c
            if j[0] == Internal.__FlowSolutionCenters__: i2 = c
            c += 1
        if i0 > i1 and i1 != -1: tmp = z[2][i1]; z[2][i1] = z[2][i0]; z[2][i0] = tmp; tmp = i1; i1 = i0; i0 = tmp
        if i0 > i2 and i2 != -1: tmp = z[2][i2]; z[2][i2] = z[2][i0]; z[2][i0] = tmp; tmp = i2; i2 = i0; i0 = tmp
        if i1 > i2 and i2 != -1: tmp = z[2][i2]; z[2][i2] = z[2][i1]; z[2][i1] = tmp; tmp = i2; i2 = i1; i1 = tmp
    return None

# -- fillMissingVariables: remplit les variables manquantes pour que toutes
# les zones aient les memes variables
def fillMissingVariables(t):
    """Fill FlowSolution nodes with variables, such that all the zones have
    the same variables."""
    tp = Internal.copyRef(t)
    _fillMissingVariables(tp)
    return tp

def _fillMissingVariables(t):
    # scan des variables
    varsn = []; varsc = []
    nodes = Internal.getZones(t)
    for z in nodes:
        vars = getVarNames(z, excludeXYZ=True, loc='nodes')[0]
        for var in vars:
            found = 0
            for v in varsn:
                if v == var: found = 1
            if found == 0: varsn.append(var)

    for z in nodes:
        vars = getVarNames(z, excludeXYZ=True, loc='centers')[0]
        for var in vars:
            found = 0
            for v in varsc:
                if v == var: found = 1
            if found == 0: varsc.append(var)

    # add vars
    _addVars(t, varsn+varsc)

    # Reorder vars for all zones
    _orderVariables(t, varsn, varsc)

    _orderContainers(t)

    return None

# -- cpVars
def cpVars(t1, var1, t2, var2):
    """Copy field variables."""
    tp2 = Internal.copyRef(t2)
    zones1 = Internal.getZones(t1)
    zones2 = Internal.getZones(tp2)
    if len(zones1) != len(zones2):
        raise TypeError("cpVars: zones in t1 and t2 must be coherents.")
    l = 0
    for z1 in zones1:
        z2 = zones2[l]
        _cpVars__(z1, var1, z2, var2)
        l += 1
    return tp2

# -- cpVars: in place in t2
def _cpVars(t1, var1, t2, var2):
    zones1 = Internal.getZones(t1)
    zones2 = Internal.getZones(t2)
    l = 0
    for z1 in zones1:
        z2 = zones2[l]
        _cpVars__(z1, var1, z2, var2)
        l += 1
    return None

# internal cpVars for zones - check supprimes car utilisation avec le solver
def _cpVars__(z1, var1, z2, var2):
    """Copy variables var1 from z1 to variables var2 of z2.
    Usage: cpVars(z1, var1, z2, var2)"""
    nodes1 = getStdNodesFromName(z1, var1)
    nodes2 = getStdNodesFromName(z2, var2)

    # Switch nodes
    if nodes2 != []: # var2 existe deja dans z2
        nodes2[0][1] = nodes1[0][1]
        return None
    # Create nodes
    s2 = var2.split(':')
    if len(s2) == 1 or s2[0] == 'nodes': loc2 = 'nodes'
    else: loc2 = 'centers'
    (r, c) = Internal.getParentOfNode(z1, nodes1[0])
    if r[0] == Internal.__GridCoordinates__:
        place2 = Internal.getNodeFromName1(z2, Internal.__GridCoordinates__)
    elif loc2 == 'nodes':
        place2 = Internal.getNodeFromName1(z2, Internal.__FlowSolutionNodes__)
    elif loc2 == 'centers':
        place2 = Internal.getNodeFromName1(z2, Internal.__FlowSolutionCenters__)

    if place2 is None: # Le containeur n'existe pas; create it
        if r[0] == Internal.__GridCoordinates__:
            z2[2].append([Internal.__GridCoordinates__, None, [], 'GridCoordinates_t'])
        elif loc2 == 'nodes':
            z2[2].append([Internal.__FlowSolutionNodes__, None, [], 'FlowSolution_t'])
        elif loc2 == 'centers':
            z2[2].append([Internal.__FlowSolutionCenters__, None, [], 'FlowSolution_t'])
            f = Internal.getNodeFromName1(z2, Internal.__FlowSolutionCenters__)
            Internal._createUniqueChild(f, 'GridLocation', 'GridLocation_t', 'CellCenter')
        h = z2[2][len(z2[2])-1]
    else:
        h = place2
    var2s = var2.split(':')
    if len(var2s) > 1: var2 = var2s[1]
    h[2].append([var2, nodes1[0][1], nodes1[0][2], nodes1[0][3]])
    return None

# -- extractVars
def extractVars(t, vars, keepOldNodes=True):
    """Only keep the given variables in t.
    Usage: extractVars(z, var)"""
    tp = Internal.copyRef(t)
    _extractVars(tp, vars, keepOldNodes)
    return tp

def _extractVars(t, vars, keepOldNodes=True):
    """Only keep the given variables in t."""
    if isinstance(vars, str): vars = [vars]
    zones = Internal.getZones(t)

    if keepOldNodes: # old legacy version
        for z in zones:
            varNames = getVarNames(z)[0]
            for v in varNames:
                if v not in vars: _rmVars(z, v)
    else:
        # rebuild a new zone with only vars and coordinates + ZoneBC + ZoneGridConnectivity
        for z in zones:
            zp = Internal.copyNode(z); zp[2] = []
            n = Internal.getNodeFromName1(z, 'ZoneType')
            if n is not None: zp[2].append(n)
            n = Internal.getNodeFromName1(z, 'ZoneRind')
            if n is not None: zp[2].append(n)
            n = Internal.getNodeFromName1(z, 'GridCoordinates')
            if n is not None: zp[2].append(n)
            if vars is not None:
                cont1 = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
                if cont1 is not None:
                    np1 = Internal.copyNode(cont1); np1[2] = []
                    zp[2].append(np1)
                    n = Internal.getNodeFromName1(cont1, 'GridLocation')
                    if n is not None: np1[2].append(n)
                else: np1 = None
                cont2 = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
                if cont2 is not None:
                    np2 = Internal.copyNode(cont2); np2[2] = []
                    zp[2].append(np2)
                    n = Internal.getNodeFromName1(cont2, 'GridLocation')
                    if n is not None: np2[2].append(n)
                else: np2 = None
                for v in vars:
                    v2 = v.split(':')
                    if len(v2) == 1: cont = cont2; np = np2
                    elif v2[0] == 'nodes': cont = cont2; v = v2[1]; np = np2
                    else: cont = cont1; v = v2[1]; np = np1
                    if cont is not None:
                        n = Internal.getNodeFromName1(cont, v)
                        if n is not None: np[2].append(n)
            n = Internal.getNodeFromName1(z, 'ZoneBC')
            if n is not None: zp[2].append(n)
            n = Internal.getNodeFromName1(z, 'ZoneGridConnectivity')
            if n is not None: zp[2].append(n)
            z[2] = zp[2]

    return None

# -- rmVars
def rmVars(z, var):
    """Remove variables var from t.
    Usage: rmVars(t, var)"""
    zc = Internal.copyRef(z)
    _rmVars(zc, var)
    return zc

def _rmVars(z, var):
    zn = Internal.getZones(z)
    for i in zn:
        if isinstance(var, list): # liste de vars
            for v in var:
                nodes = getStdNodesFromName(i, v)
                if nodes != []:
                    (parent, d) = Internal.getParentOfNode(i, nodes[0])
                    del parent[2][d]
        else: # une seule var
            nodes = getStdNodesFromName(i, var)
            if nodes != []:
                (parent, d) = Internal.getParentOfNode(i, nodes[0])
                del parent[2][d]
    return None

def rmBCDataVars(t,var):
    """Remove variables var from t in every BCDataSet.
       Usage rmBCDataVars(t,var)"""
    tc = Internal.copyRef(t)
    _rmBCDataVars(tc,var)
    return tc

def _rmBCDataVars(t,var):
    if not isinstance(var, list): vars = [var]
    else: vars = var

    isStd = Internal.isStdNode(t)
    if isStd >= 0:
        for c in t[isStd:]: rmBCDataVars__(c,vars)
    else:
        rmBCDataVars__(t,vars)

def rmBCDataVars__(t,vars):
    for v in vars:
        bcnodes = Internal.getNodesFromType(t, 'BC_t')
        for bcnode in bcnodes:
            bcdata = Internal.getNodeFromType2(bcnode, 'BCData_t')
            nodes  = Internal.getNodesFromName(bcdata, v)
            if nodes != []:
                (parent, d) = Internal.getParentOfNode(bcdata, nodes[0])
                del parent[2][d]
    return None

# -- normalize: normalise un jeu de variables
def normalize(t, vars):
    """Normalize the field defined by vars in the tree.
    Usage: normalize(t, vars)"""
    loc = ''; vars2 = []
    for v in vars:
        s = v.split(':')
        if len(s) == 2 and s[0] == 'centers': #centres
            if loc == '': loc = s[0]
            elif loc != s[0]: raise ValueError("normalize: vector components must have the same location (centers or nodes).")
            vars2.append(s[1])
        elif len(s) == 1 or (len(s) == 2 and s[0] == 'nodes'): #noeuds
            if loc == '': loc = 'nodes'
            elif loc == 'centers': raise ValueError("normalize: vector components must have the same location (centers or nodes).")
            if len(s) == 2: vars2.append(s[1])
            else: vars2.append(s[0])
        else: raise ValueError("normalize: invalid vector component.")
    return TZA3(t, loc, loc, False, Converter.normalize, vars2)

def _normalize(t, vars):
    loc = ''; vars2 = []
    for v in vars:
        s = v.split(':')
        if len(s) == 2 and s[0] == 'centers': #centres
            if loc == '': loc = s[0]
            elif loc != s[0]: raise ValueError("normalize: vector components must have the same location (centers or nodes).")
            vars2.append(s[1])
        elif len(s) == 1 or (len(s) == 2 and s[0] == 'nodes'): #noeuds
            if loc == '': loc = 'nodes'
            elif loc == 'centers': raise ValueError("normalize: vector components must have the same location (centers or nodes).")
            if len(s) == 2: vars2.append(s[1])
            else: vars2.append(s[0])
        else: raise ValueError("normalize: invalid vector component.")
    __TZA3(t, loc, Converter._normalize, vars2)
    return None

# -- magnitude: calcul la norme d'un jeu de variables
def magnitude(t, vars):
    """Compute the magnitude of the fields defined by vars in the tree.
    Usage: magnitude(t, vars)"""
    tp = Internal.copyRef(t)
    _magnitude(tp, vars)
    return tp

def _magnitude(t, vars):
    loc = ''
    vars2 = []
    for v in vars:
        s = v.split(':')
        if len(s) == 2: #centres
            if loc == '': loc = s[0]
            elif loc != s[0]: raise ValueError("magnitude: vector components must have the same location (centers or nodes).")
            vars2.append(s[1])
        elif len(s) == 1: #noeuds
            if loc == '': loc = 'nodes'
            elif loc == 'centers': raise ValueError("magnitude: vector components must have the same location (centers or nodes).")
        else:
            raise ValueError("magnitude: invalid vector component.")

    if loc == 'nodes':
        _TZA3(t, loc, loc, False, Converter.magnitude, vars)
    else:
        _TZA3(t, loc, loc, False, Converter.magnitude, vars2)
    return None

# -- normL0
def normL0(t, var):
    """Get the L0 norm of the field defined by varName in t.
    If celln exists in the array, the norm for blanked points is not computed.
    Usage: normL0(t, varName)"""
    A = getField(var, t, api=3)
    v = var.split(':')
    if len(v) > 1: var = v[1]
    return Converter.normL0(A, var)

# -- normL2
def normL2(t, var):
    """Get the L2 norm of the field defined by varName in t.
    If celln exists in the array, the norm for blanked points is not computed.
    Usage: normL0(t, varName)"""
    A = getField(var, t, api=3)
    v = var.split(':')
    if len(v) > 1: var = v[1]
    return Converter.normL2(A, var)

# -- getArgMin (only on nodes or centers separately)
def getArgMin(t, var):
    """Get field values where the variable defined by varName is minimum.
    Usage: getArgMin(t, var)"""
    v = var.split(':')
    centers = False
    if len(v) > 1:
        var = v[1]
        if v[0] == 'centers': centers = True
    if centers:
        A = getFields([Internal.__FlowSolutionCenters__], t, api=2)
    else:
        A = getFields([Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__], t, api=2)
    return Converter.getArgMin(A, var)

# -- getArgMax (only on nodes or centers separately)
def getArgMax(t, var):
    """Get field values where the variable defined by varName is maximum.
    Usage: getArgMax(t, var)"""
    v = var.split(':')
    centers = False
    if len(v) > 1:
        var = v[1]
        if v[0] == 'centers': centers = True
    if centers:
        A = getFields([Internal.__FlowSolutionCenters__], t, api=2)
    else:
        A = getFields([Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__], t, api=2)
    return Converter.getArgMax(A, var)

# -- getMinValue
def getMinValue(t, varName):
    """Get the minimum value of variable defined by var.
    Usage: getMinValue(t, var)"""
    if not isinstance(varName, list): varNames = [varName]
    else: varNames = varName
    out = []
    if varNames[0] == Internal.__GridCoordinates__: varNames = ['CoordinateX', 'CoordinateY', 'CoordinateZ']
    for v in varNames:
        A = getField(v, t, api=3)
        va = v.split(':')
        if len(va) > 1: v = va[1]
        minValue = 1.e9
        for i in A:
            if i != []:
                minValue = min(minValue, Converter.getMinValue(i, v))
        out.append(minValue)
    if len(out) == 1: return out[0]
    else: return out

# -- getMaxValue
def getMaxValue(t, varName):
    """Get the maximum value of variable defined by var.
    Usage: getMaxValue(t, var)"""
    if not isinstance(varName, list): varNames = [varName]
    else: varNames = varName
    out = []
    if varNames[0] == Internal.__GridCoordinates__: varNames = ['CoordinateX', 'CoordinateY', 'CoordinateZ']
    for v in varNames:
        A = getField(v, t, api=3)
        va = v.split(':')
        if len(va) > 1: v = va[1]
        maxValue = -1.e9
        for i in A:
            if i != []:
                maxValue = max(maxValue, Converter.getMaxValue(i, v))
        out.append(maxValue)
    if len(out) == 1: return out[0]
    else: return out

# -- getMeanValue
def getMeanValue(t, var):
    """Get the mean value of variable defined by var.
    Usage: getMeanValue(t, var)"""
    A = getField(var, t, api=3)
    v = var.split(':')
    if len(v) > 1: var = v[1]
    return Converter.getMeanValue(A, var)

# -- getMeanRangeValue
def getMeanRangeValue(t, var, rmin, rmax):
    """Get the mean value of variable defined by var in the sorted range between rmin and rmax.
    Usage: getMeanRangeValue(t, var, rmin, rmax)"""
    A = getField(var, t, api=2)
    v = var.split(':')
    if len(v) > 1: var = v[1]
    return Converter.getMeanRangeValue(A, var, rmin, rmax)

#==============================================================================
# -- Topology conversions --
#==============================================================================

# -- Convert a BAR array without branches, closed into an i-array
def convertBAR2Struct(t):
    """Convert a BAR array without branches, closed into an i-array.
    Usage: convertBAR2Struct(t)"""
    return TZA1(t, 'both', 'nodes', True, Converter.convertBAR2Struct)

def _convertBAR2Struct(t):
    _TZA1(t, 'both', 'nodes', True, Converter.convertBAR2Struct)
    return None

# -- convertArray2Tetra
# split='simple': pas d'ajout de points
# split='withBarycenters': ajout de points aux centres des elements et des faces
def convertArray2Tetra(t, split='simple'):
    """Convert a zone to an unstructured zone.
    Unstructured array can be triangular in 2D and tetraedrical in 3D.
    Usage: convertArray2Tetra(t, split)"""
    tp = Internal.copyRef(t)
    _convertArray2Tetra(tp, split)
    return tp

def _convertArray2Tetra(t, split='simple'):
    _deleteZoneBC__(t)
    _deleteGridConnectivity__(t)
    if split == 'simple':
        _TZANC(t, 'both', 'both', Converter.convertArray2Tetra,
               Converter.convertArray2Tetra, split, split)
    else:
        nodes = Internal.getZones(t)
        for z in nodes:
            fieldn = getAllFields(z, 'nodes')[0]
            fieldc = getAllFields(z, 'centers')[0]
            if fieldc == []:
                res = Converter.convertArray2Tetra1__(fieldn, split='withBarycenters')
                z = setFields([res], z, 'nodes', writeDim=True)
            else:
                res = Converter.convertArray2Tetra1__(fieldn, fieldc, split='withBarycenters')
                z = setFields([res[0]], z, 'nodes', writeDim=True)
                z = setFields([res[1]], z, 'centers', writeDim=False)
    return None

# -- convertArray2Hexa
def convertArray2Hexa(t):
    """Convert a structured zone to an unstructured quad/hexa zone.
    Convert an unstructured zone to a quad/hexa zone. If the original
    zone is a TRI,TETRA or PENTA zone, return a QUAD/HEXA/HEXA zone with
    degenerated edges.
    Usage: convertArray2Hexa(t)"""
    tp = Internal.copyRef(t)
    _convertArray2Hexa(tp)
    return tp

def _convertArray2Hexa(t):
    _deleteZoneBC__(t)
    _deleteGridConnectivity__(t)
    _TZA1(t, 'both', 'nodes', True, Converter.convertArray2Hexa)
    return None

# -- convertArray2NGon
def convertArray2NGon(t, recoverBC=True, api=1):
    """Convert a zone to a NGON zone.
    Usage: convertArray2NGon(t)"""
    tp = Internal.copyRef(t)
    _convertArray2NGon(tp, recoverBC=recoverBC, api=api)
    return tp

def _convertArray2NGon(t, recoverBC=True, api=1):
    zones = Internal.getZones(t)
    if recoverBC:
        gbcs = []
        for z in zones:
            dims = Internal.getZoneDim(z)
            if dims[0] == 'Unstructured' and dims[3] == 'NGON': gbcs.append([])
            else: gbcs.append(getBCs(z, extrapFlow=False))
    else: _deleteZoneBC__(t)
    _deleteGridConnectivity__(t)
    _TZA1(t, 'both', 'nodes', True, Converter.convertArray2NGon)
    #_TZAX(api, t, 'both', 'nodes', True, Converter.convertArray2NGon, api)
    Internal._fixNGon(t)

    # Recover BCs for NGon
    if recoverBC:
        zones = Internal.getZones(t); c = 0
        for z in zones:
            if gbcs[c] != [] and len(gbcs[c][0]) > 0: _recoverBCs(z, gbcs[c])
            c += 1
    return None

# -- convertPenta2Strand
def convertPenta2Strand(t):
    """Convert PENTA to STRAND."""
    tp = Internal.copyRef(t)
    _convertPenta2Strand(tp)
    return tp

def _convertPenta2Strand(t):
    """Convert PENTA to STRAND."""
    _deleteFlowSolutions__(t, 'centers')
    _deleteZoneBC__(t)
    _deleteGridConnectivity__(t)
    _TZA3(t, 'nodes', 'nodes', True, Converter.convertPenta2Strand)
    return None

def convertStrand2Penta(t):
    """Convert STRAND to PENTA."""
    tp = Internal.copyRef(t)
    _convertStrand2Penta(tp)
    return tp

# -- convertStrand2Penta
def _convertStrand2Penta(t):
    """Convert STRAND to PENTA."""
    _deleteFlowSolutions__(t, 'centers')
    _deleteZoneBC__(t)
    _deleteGridConnectivity__(t)
    _TZA3(t, 'nodes', 'nodes', True, Converter.convertStrand2Penta)
    return None

# -- convertArray2Node
def convertArray2Node(t):
    """Convert an array to a node array.
    Usage: convertArray2Node(t)"""
    tp = Internal.copyRef(t)
    _convertArray2Node(tp)
    return tp

def _convertArray2Node(t):
    _deleteFlowSolutions__(t, 'centers')
    _deleteZoneBC__(t)
    _deleteGridConnectivity__(t)
    zones = Internal.getZones(t)
    for z in zones:
        n = Internal.getNodeFromName1(z, 'ZoneType')
        Internal.setValue(n, 'Unstructured')
        x = Internal.getNodeFromName2(z, 'CoordinateX')
        z[1] = numpy.empty((1,3), Internal.E_NpyInt, order='F')
        z[1][0,0] = x[1].size; z[1][0,1] = 0; z[1][0,2] = 0
        Internal._rmNodesByType(z, 'Elements_t')
        n = Internal.createChild(z, 'GridElements', 'Elements_t', [2,0])
        Internal.createChild(n, 'ElementRange', 'IndexRange_t', [1,0])
        Internal.createChild(n, 'ElementConnectivity', 'DataArray_t', None)
    return None

# -- convertTri2Quad
def convertTri2Quad(z, alpha=30.):
    """Convert a TRI zone to a QUAD zone.
    Usage: convertTri2Quad(z, angle)"""
    a = getAllFields(z, 'nodes')
    b, c = Converter.convertTri2Quad(a, alpha)
    z1 = convertArrays2ZoneNode(getZoneName('quad'), [b])
    z2 = convertArrays2ZoneNode(getZoneName('tri'), [c])
    return z1, z2

# -- conformizeNGon
def conformizeNGon(a, tol=1.e-6):
    """Conformize topologically a NGON zone.
    Usage: conformizeNGon(a, tol)"""
    return TZGC1(a, 'nodes', True, Converter.conformizeNGon, tol)

def _conformizeNGon(a, tol=1.e-6):
    _TZGC1(a, 'nodes', True, Converter.conformizeNGon, tol)
    return None

# -- convertSurfaceNGon
def convertSurfaceNGon(a, rmEmptyNFaceElements=True):
    """Convert a surface NGon from (A: NGON=bars, NFACE=polygon)
    to (B: NGON=polygon, NFACE=NULL), or vice versa.
    Usage: convertSurfaceNGon(a)"""
    # add NFace node if necessary
    for z in Internal.getZones(a):
        nFace = Internal.getNodeFromName(z, 'NFaceElements')
        if nFace is None:
            nGon = Internal.getNodeFromName(z, 'NGonElements')
            offset = Internal.getNodeFromName(nGon, 'ElementStartOffset')
            api = 3 if offset is not None else 2
            rnGon = Internal.getNodeFromName(nGon, 'ElementRange')[1]

            nface = Internal.createNode('NFaceElements', 'Elements_t', parent=z,
                                        value=numpy.array([23,0],
                                        dtype=Internal.E_NpyInt, order='F'))

            value = numpy.array([rnGon[1]+1, rnGon[1]+1],
                                 dtype=Internal.E_NpyInt, order='F')
            Internal.createNode('ElementRange', 'IndexRange_t',
                                parent=nface, value=value)
            value = numpy.array([], dtype=Internal.E_NpyInt, order='F')
            Internal.createNode('ElementConnectivity', 'DataArray_t',
                                parent=nface, value=value)
            if api == 3:
                value = numpy.array([0], dtype=Internal.E_NpyInt, order='F')
            else: value =  numpy.array([], dtype=Internal.E_NpyInt, order='F')
            Internal.createNode('ElementStartOffset', 'DataArray_t',
                                parent=nface, value=value)

    a = TZGC3(a, 'nodes', True, Converter.convertSurfaceNGon)

    if rmEmptyNFaceElements:
        for z in Internal.getZones(a):
            nFace = Internal.getNodeFromName(z, 'NFaceElements')
            cnFace = Internal.getNodeFromName(nFace, 'ElementConnectivity')[1]
            if cnFace.size == 0: Internal._rmNodesByName(z, 'NFaceElements')

    return a

#=============================================================================
# -- Create BC(s) to a zone node --
#=============================================================================

# -- addBC2Zone
def addBC2Zone(a, bndName, bndType, wrange=[],
               zoneDonor=[], rangeDonor=[], trirac=[1,2,3],
               rotationCenter=[], rotationAngle=[], translation=[],
               faceList=None, pointList=None, elementList=None, elementRange=[], data=None,
               subzone=None, faceListDonor=None, elementListDonor=None,
               elementRangeDonor=None, tol=1.e-12, unitAngle=None):
    """Add a BC to a zone node.
    Usage: addBC2Zone(zone, bndName, bndType, wrange)"""
    ap = Internal.copyRef(a)
    _addBC2Zone(ap, bndName, bndType, wrange=wrange,
                zoneDonor=zoneDonor, rangeDonor=rangeDonor, trirac=trirac,
                rotationCenter=rotationCenter, rotationAngle=rotationAngle, translation=translation,
                faceList=faceList, pointList=pointList, elementList=elementList, elementRange=elementRange, data=data, subzone=subzone,
                faceListDonor=faceListDonor, elementListDonor=elementListDonor, elementRangeDonor=elementRangeDonor,
                tol=tol, unitAngle=unitAngle)
    return ap

def _addBC2Zone(a, bndName, bndType, wrange=[],
                zoneDonor=[], rangeDonor=[], trirac=[1,2,3],
                rotationCenter=[], rotationAngle=[], translation=[],
                faceList=None, pointList=None, elementList=None, elementRange=[], data=None,
                subzone=None, faceListDonor=None, elementListDonor=None,
                elementRangeDonor=None, tol=1.e-12, unitAngle=None):
    bndName = getBCName(bndName)
    zones = Internal.getZones(a)
    for z in zones:
        dims = Internal.getZoneDim(z)
        if dims[0] == 'Unstructured':
            eltType = dims[3]
            if eltType == 'NGON':
                if faceList is None and subzone is None:
                    raise TypeError("addBC2Zone: unstructured grids requires a faceList or a subzone.")
                _addBC2NGonZone__(z, bndName, bndType, faceList=faceList, data=data, subzone=subzone,
                                  zoneDonor=zoneDonor, faceListDonor=faceListDonor,
                                  rotationCenter=rotationCenter, rotationAngle=rotationAngle, translation=translation,
                                  tol=tol, unitAngle=unitAngle)
            else: # basic elements
                if elementList is None and elementRange == [] and subzone is None and faceList is None and pointList is None:
                    raise TypeError("addBC2Zone: unstructured grids requires an elementList, a elementRange or a subzone.")
                _addBC2UnstructZone__(z, bndName, bndType, elementList=elementList, elementRange=elementRange,
                                      faceList=faceList, pointList=pointList, data=data, subzone=subzone,
                                      zoneDonor=zoneDonor, elementListDonor=elementListDonor, elementRangeDonor=elementRangeDonor,
                                      faceListDonor=faceListDonor, rotationCenter=rotationCenter, rotationAngle=rotationAngle,
                                      translation=translation, unitAngle=unitAngle)
        else: # structured
            _addBC2StructZone__(z, bndName, bndType, wrange=wrange, faceList=faceList,
                                zoneDonor=zoneDonor, rangeDonor=rangeDonor, faceListDonor=faceListDonor,
                                trirac=trirac, rotationCenter=rotationCenter, rotationAngle=rotationAngle,
                                translation=translation, data=data, unitAngle=unitAngle)
    return None

# -- Ajout Infos peridiques pour les BC match periodiques
def _addPeriodicInfoInGC__(info, rotationCenter, rotationAngle, translation, unitAngle=None):
    if len(rotationCenter) != 0 or len(rotationAngle) != 0 or len(translation) != 0:
        if len(rotationCenter) == 0: rotCenter = (0,0,0)
        else: rotCenter = rotationCenter
        if len(rotationAngle) == 0: rotAngle = (0,0,0)
        else: rotAngle = rotationAngle
        if len(translation) == 0: trans = (0,0,0)
        else: trans = translation

        info[2].append(['GridConnectivityProperty', None, [], 'GridConnectivityProperty_t'])
        info = info[2][len(info[2])-1]
        info[2].append(['Periodic', None, [], 'Periodic_t'])
        info = info[2][0]
        if len(rotCenter) != 3:
            raise ValueError("_addPeriodicInfoInGC__: rotationCenter must be (Cx, Cy, Cz).")
        if len(rotAngle) != 3:
            raise ValueError("_addPeriodicInfoInGC__: rotationAngle must be (alpha, beta, gamma).")
        if len(trans) != 3:
            raise ValueError("_addPeriodicInfoInGC__: translation must be (tx, ty, tz).")

        v = numpy.zeros((3), numpy.float64)
        v[0] = rotCenter[0]; v[1] = rotCenter[1]; v[2] = rotCenter[2]
        info[2].append(['RotationCenter', v, [], 'DataArray_t'])

        v = numpy.zeros((3), numpy.float64)
        if unitAngle is None:
            unitAngle='Radian'
            RA1 = rotAngle[0]*Internal.__DEG2RAD__
            RA2 = rotAngle[1]*Internal.__DEG2RAD__
            RA3 = rotAngle[2]*Internal.__DEG2RAD__
        else:
            RA1 = rotAngle[0]
            RA2 = rotAngle[1]
            RA3 = rotAngle[2]
        v[0] = RA1; v[1] = RA2; v[2] = RA3
        info[2].append(['RotationAngle', v, [], 'DataArray_t'])

        if unitAngle is not None:
            Internal.newDimensionalUnits(angleUnit=unitAngle, parent=info[2][-1])

        v = numpy.zeros((3), numpy.float64)
        v[0] = trans[0]; v[1] = trans[1]; v[2] = trans[2]
        info[2].append(['Translation', v, [], 'DataArray_t'])
    return None

# typeZone=0 : struct, 1 : NGON, 2: unstructured BE
def _addFamilyOfStageGC__(z, bndName, bndType2, typeZone=0, faceList=None, elementList=None, elementRange=[], pointRange=[], zoneDonor=[]):
    zoneGC = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
    if zoneGC == []:
        z[2].append(['ZoneGridConnectivity', None, [], 'ZoneGridConnectivity_t'])
        zoneGC = z[2][len(z[2])-1]
    else:
        zoneGC = zoneGC[0]
    # Cree le noeud de GC
    if zoneDonor == []:
        # autoattach
        Internal._createChild(zoneGC, bndName, 'GridConnectivity_t', value=z[0])
    else:
        # donors donnes
        st = ""
        for i in zoneDonor:
            if isinstance(i, str):
                if st == "": st = i
                else: st = st+","+i
            else:
                if st == "": st = i[0]
                else: st = st+","+i[0]
        Internal._createChild(zoneGC, bndName, 'GridConnectivity_t', value=st)
    info = zoneGC[2][len(zoneGC[2])-1]

    if typeZone == 0: # STRUCTURED
        if pointRange == []: raise ValueError("_addFamilyOfStageGC__: pointRange is empty.")
        r = Internal.window2Range(pointRange)
        info[2].append(['PointRange', r, [], 'IndexRange_t'])

    elif typeZone == 1: # NGON
        if faceList is None: raise ValueError("_addFamilyOfStageGC__: faceList is empty.")

        if isinstance(faceList, numpy.ndarray): r = faceList
        else: r = numpy.array(faceList, dtype=Internal.E_NpyInt)
        r = r.reshape((1,r.size), order='F')
        info[2].append([Internal.__FACELIST__, r, [], 'IndexArray_t'])

    elif typeZone == 2: # UNS BE
        if elementList is not None:
            if isinstance(elementList, numpy.ndarray): r = elementList
            else: r = numpy.array(elementList, dtype=Internal.E_NpyInt)
            r = r.reshape((1,r.size), order='F')
            info[2].append([Internal.__ELEMENTLIST__, r, [], 'IndexArray_t'])
        elif elementRange != []:
            r = numpy.empty((1,2), dtype=Internal.E_NpyInt, order='F')
            r[0,0] = elementRange[0]
            r[0,1] = elementRange[1]
            info[2].append([Internal.__ELEMENTRANGE__, r, [], 'IndexRange_t'])
        elif faceList is not None:
            Internal._createChild(info, 'GridLocation', 'GridLocation_t', value='FaceCenter')
            if isinstance(faceList, numpy.ndarray): r = faceList
            else: r = numpy.array(faceList, dtype=Internal.E_NpyInt)
            r = r.reshape((1, r.size), order='F')
            info[2].append([Internal.__FACELIST__, r, [], 'IndexArray_t'])
        else:
            raise ValueError("_addFamilyOfStageGC__: elementList, elementRange and faceList are all empty.")

    else:
        raise ValueError("_addFamilyOfStageGC__: typeZone not valid.")

    Internal.createChild(info, 'GridConnectivityType', 'GridConnectivityType_t', 'Abutting')
    Internal.createChild(info, 'FamilyStage', 'FamilyName_t', bndType2)
    return None

# -- addBC2Zone pour les grilles structurees
def _addBC2StructZone__(z, bndName, bndType, wrange=[], faceList=[],
                        zoneDonor=[], rangeDonor=[], faceListDonor=[], trirac=[1,2,3],
                        rotationCenter=[], rotationAngle=[],
                        translation=[], data=None, unitAngle=None):
    # Type de condition aux limites definies par une famille
    # par ex: FamilySpecified:OVERLAP
    s = bndType.split(':')
    bndType1 = s[0]
    if len(s) > 1: bndType2 = s[1]
    else: bndType2 = ''

    typeR = -1 # 0 : PR, 1 PL
    if wrange != []: typeR=0
    elif len(faceList) != 0: typeR=1
    else: raise ValueError("addBC2Zone: match connectivity requires a range or face list.")

    # Range defini par une chaine ou non
    if typeR == 0:
        if isinstance(wrange, str):
            wrange = convertStringRange2Range__(wrange, z) #fenetre complete
        if isinstance(rangeDonor, str):
            if isinstance(zoneDonor, str):
                raise ValueError("addBC2Zone: donor range must be explicitly specified.")
            else:
                rangeDonor = convertStringRange2Range__(rangeDonor, zoneDonor)

    if bndType1 == 'BCMatch' or bndType1 == 'Abutting1to1':
        if typeR == 0 and rangeDonor == [] and trirac == []:
            raise ValueError("addBC2Zone: match connectivity requires a donor point range and a trirac.")
        elif typeR == 1 and len(faceListDonor) == 0:
            raise ValueError("addBC2Zone: match connectivity requires a donor face list.")

        if zoneDonor == []:
            raise ValueError("addBC2Zone: match connectivity requires a donor zone.")
        # Cree le noeud zoneGridConnectivity si besoin
        zoneGC = Internal.getNodeFromType1(z, 'ZoneGridConnectivity_t')
        if zoneGC is None:
            z[2].append(['ZoneGridConnectivity', None, [], 'ZoneGridConnectivity_t'])
            zoneGC = z[2][len(z[2])-1]
        # Cree le noeud de GC1-1
        if isinstance(zoneDonor, str): v = zoneDonor
        else: v = zoneDonor[0]
        Internal._createChild(zoneGC, bndName, 'GridConnectivity1to1_t', value=v)

        l = len(zoneGC[2])
        info = zoneGC[2][l-1]

        if typeR == 0:
            r = Internal.window2Range(wrange)
            info[2].append([Internal.__FACERANGE__, r, [], 'IndexRange_t'])
            o = Internal.window2Range(rangeDonor)
            info[2].append([Internal.__FACERANGE__+'Donor', o, [], 'IndexRange_t'])
            size = len(trirac)
            o = numpy.zeros((size), dtype=Internal.E_NpyInt)
            for i in range(size): o[i] = trirac[i]
            info[2].append(['Transform', o, [], '\"int[IndexDimension]\"'])

        elif typeR == 1:
            Internal._createChild(info, 'GridLocation', 'GridLocation_t', value='FaceCenter')
            if isinstance(faceList, numpy.ndarray): r = faceList
            else: r = numpy.array(faceList, dtype=Internal.E_NpyInt)
            r = r.reshape((1,r.size), order='F')
            info[2].append([Internal.__FACELIST__, r, [], 'IndexArray_t'])
            if isinstance(faceListDonor, numpy.ndarray): r = faceListDonor
            else: r = numpy.array(faceListDonor, dtype=Internal.E_NpyInt)
            r = r.reshape((1,r.size), order='F')
            info[2].append([Internal.__FACELIST__+'Donor', r, [], 'IndexArray_t'])
            Internal.createChild(info, 'GridConnectivityType', 'GridConnectivityType_t', 'Abutting1to1')

        # Ajout pour les BC match periodiques
        _addPeriodicInfoInGC__(info, rotationCenter, rotationAngle, translation, unitAngle=unitAngle)

    elif bndType1 == 'BCNearMatch':
        if typeR != 0:
            raise ValueError("addBC2Zone: nearmatch connectivity requires a point range.")

        if rangeDonor == [] or zoneDonor == [] or trirac == []:
            raise ValueError("addBC2Zone: nearmatch connectivity requires a donor and transform.")

        # Cree le noeud zoneGridConnectivity si besoin
        zoneGC = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
        if zoneGC == []:
            Internal._createChild(z, 'ZoneGridConnectivity', 'ZoneGridConnectivity_t', value=None)
            zoneGC = z[2][len(z[2])-1]
        else:
            zoneGC = zoneGC[0]
        if isinstance(zoneDonor, str): v = zoneDonor
        else: v = zoneDonor[0]
        Internal._createChild(zoneGC, bndName, 'GridConnectivity_t', value=v)

        l = len(zoneGC[2])

        info = zoneGC[2][l-1]
        r = Internal.window2Range(wrange)
        info[2].append(['PointRange', r, [], 'IndexRange_t'])
        Internal.createChild(info, 'GridConnectivityType', 'GridConnectivityType_t', 'Abutting')
        c = numpy.ones((3,1), dtype=Internal.E_NpyInt)
        info[2].append(['PointListDonor', c, [], 'IndexArray_t'])

        # Nearmatch: rajoute les noeuds PointRangeDonor et Transform en UserDefinedData
        info[2].append(['UserDefinedData', None, [], 'UserDefinedData_t'])
        l = len(info[2])-1; info = info[2][l]
        dd = Internal.window2Range(rangeDonor)
        info[2].append(['PointRangeDonor', dd, [], 'DataArray_t'])
        size = len(trirac); o = numpy.zeros((size), dtype=Internal.E_NpyInt)
        for i in range(size): o[i] = trirac[i]
        info[2].append(['Transform', o, [], 'DataArray_t'])
        ratio = getNMRatio__(wrange, rangeDonor, trirac)
        size = len(ratio); ro = numpy.zeros((size), numpy.float64)
        for i in range(size): ro[i] = ratio[i]
        info[2].append(['NMRatio', ro, [], 'DataArray_t'])

    elif bndType1 == 'BCOverlap':
        # Cree le noeud zoneGridConnectivity si besoin
        zoneGC = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
        if zoneGC == []:
            z[2].append(['ZoneGridConnectivity', None, [], 'ZoneGridConnectivity_t'])
            zoneGC = z[2][len(z[2])-1]
        else:
            zoneGC = zoneGC[0]

        # Cree le noeud de GC
        if zoneDonor == []:
            # autoattach
            v = z[0]

        else:
            # on specifie les donneurs sous forme de liste de noeuds zones, noms de zones, FamilySpecified:FAMILYZONES
            dnrZoneNames=[]
            for i in zoneDonor:
                if isinstance(i,str):
                    isp = i.split(":")
                    if len(isp)==1: dnrZoneNames.append(i)
                    else:
                        if isp[0]=='FamilySpecified': dnrZoneNames.append(isp[1])

                elif Internal.getType(i)=="Zone_t": dnrZoneNames.append(Internal.getName(i))
            if dnrZoneNames == []: raise ValueError("addBC2ZoneStruct: no donor zone for doubly defined overlap bc %s."%(zoneGC[0]))
            # donors donnes
            v = ",".join(dnrZoneNames)
            #print("addBC2ZoneStruct(overlap): liste des zones donneuses : ", v)

        Internal._createChild(zoneGC, bndName, 'GridConnectivity_t', value=v)

        info = zoneGC[2][len(zoneGC[2])-1]
        d = numpy.ones((3,1), dtype=Internal.E_NpyInt)
        c = numpy.ones((3,1), dtype=numpy.float64)
        r = Internal.window2Range(wrange)
        info[2].append(['PointRange', r, [], 'IndexRange_t'])
        Internal._createChild(info, 'GridConnectivityType', 'GridConnectivityType_t', value='Overset')
        info[2].append(['CellListDonor', d, [], 'IndexArray_t'])
        info[2].append(['InterpolantsDonor', c, [], 'DataArray_t'])
        if rangeDonor == 'doubly_defined':
            info[2].append(['UserDefinedData', None, [], 'UserDefinedData_t'])
            l = len(info[2])-1; info = info[2][l]
            dd = numpy.ones((1), dtype=Internal.E_NpyInt)
            info[2].append(['doubly_defined', dd, [], 'DataArray_t'])

    elif (bndType1 == 'FamilySpecified' and fnmatch.fnmatch(bndType2, 'BCStage*')) or (bndType1 == 'BCStage'):
        _addFamilyOfStageGC__(z, bndName, bndType2, typeZone=0, pointRange=wrange, zoneDonor=zoneDonor)

    else: # classical BC
        # Cree le noeud zoneBC si besoin
        zoneBC = Internal.getNodesFromType1(z, 'ZoneBC_t')
        if zoneBC == []:
            z[2].append(['ZoneBC', None, [], 'ZoneBC_t'])
            zoneBC = z[2][len(z[2])-1]
        else: zoneBC = zoneBC[0]
        # Cree le noeud de BC
        Internal._createChild(zoneBC, bndName, 'BC_t', value=bndType1)
        l = len(zoneBC[2])
        info = zoneBC[2][l-1]
        r = Internal.window2Range(wrange)
        info[2].append(['PointRange', r, [], 'IndexRange_t'])

        # Ajoute la famille si necessaire
        if bndType1 == 'FamilySpecified':
            Internal.createChild(info, 'FamilyName', 'FamilyName_t', bndType2)

        # Ajoute les Data si necessaire (donnees Dirichlet)
        if data is not None:
            node1 = Internal.createNode('State', 'DataArray_t', value=data)
            node2 = Internal.createNode('DirichletData', 'BCData_t', children=[node1])
            node3 = Internal.createNode('BCDataSet', 'BCDataSet_t', children=[node2])
            info[2].append(node3)
    return None

# -- addBC2Zone pour les zones NGon
def _addBC2NGonZone__(z, bndName, bndType, faceList, data, subzone,
                      zoneDonor, faceListDonor,
                      rotationCenter, rotationAngle, translation, tol, unitAngle=None):
    # Type de condition aux limites definies par une famille
    s = bndType.split(':')
    bndType1 = s[0]
    if len(s) > 1: bndType2 = s[1]
    else: bndType2 = ''

    # si subzone: on cree le faceList par identification
    if subzone is not None:
        hook = createHook(z, 'faceCenters')
        faceList = identifyElements(hook, subzone, tol)
        freeHook(hook)

    if bndType1 == 'BCMatch' or bndType1 == 'Abutting1to1':
        if (zoneDonor == [] or (faceListDonor is None and subzone is None)):
            raise ValueError("addBC2Zone: NGON match connectivity requires a donor and a faceListDonor or a subzone.")
        # cree le faceListDonor si subzone fourni
        if (subzone is not None and faceListDonor is None):
            hook = createHook(zoneDonor, 'faceCenters')
            faceListDonor = identifyElements(hook, subzone, tol)
            freeHook(hook)

        # Cree le noeud zoneGridConnectivity si besoin
        zoneGC = Internal.getNodeFromType1(z, 'ZoneGridConnectivity_t')
        if zoneGC is None:
            z[2].append(['ZoneGridConnectivity', None, [], 'ZoneGridConnectivity_t'])
            zoneGC = z[2][len(z[2])-1]

        if isinstance(zoneDonor, str): v = zoneDonor
        else: v = zoneDonor[0]
        Internal._createChild(zoneGC, bndName, 'GridConnectivity_t', value=v)

        l = len(zoneGC[2])
        info = zoneGC[2][l-1]
        Internal.createChild(info, 'GridLocation', 'GridLocation_t', 'FaceCenter')
        if isinstance(faceList, numpy.ndarray): r = faceList
        else: r = numpy.array(faceList, dtype=Internal.E_NpyInt)
        r = r.reshape((1, r.size), order='F')
        info[2].append([Internal.__FACELIST__, r, [], 'IndexArray_t'])
        if isinstance(faceListDonor, numpy.ndarray): r = faceListDonor
        else: r = numpy.array(faceListDonor, dtype=Internal.E_NpyInt)
        r = r.reshape((1, r.size), order='F')
        info[2].append([Internal.__FACELIST__+'Donor', r, [], 'IndexArray_t'])
        Internal.createChild(info, 'GridConnectivityType', 'GridConnectivityType_t', 'Abutting1to1')
        # Ajout pour les BC match periodiques
        _addPeriodicInfoInGC__(info, rotationCenter, rotationAngle, translation, unitAngle=unitAngle)

    elif bndType1 == 'BCNearMatch':
        print('Warning: addBC2Zone: BCNearMatch not valid for NGON zones.')

    # elif bndType1 == 'BCOverlap':
    #   # Cree le noeud zoneGridConnectivity si besoin
    #   zoneGC = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
    #   if zoneGC == []:
    #     z[2].append(['ZoneGridConnectivity', None, [], 'ZoneGridConnectivity_t'])
    #     zoneGC = z[2][len(z[2])-1]
    #   else:
    #     zoneGC = zoneGC[0]
    #   # Cree le noeud de GC
    #   if zoneDonor == []:
    #     # autoattach
    #     v = numpy.fromstring(z[0], 'c')
    #     zoneGC[2].append([bndName, v, [], 'GridConnectivity_t'])
    #   else:
    #     # donors donnes
    #     st = ""
    #     for i in zoneDonor:
    #       if isinstance(i, str):
    #         if (st == ""): st = i
    #         else: st = st+","+i
    #       else:
    #         if (st == ""): st = i[0]
    #         else: st = st+","+i[0]
    #     v = numpy.fromstring(st, 'c')
    #     zoneGC[2].append([bndName, v, [], 'GridConnectivity_t'])

        # l = len(zoneGC[2])
        # info = zoneGC[2][l-1]
        # v = numpy.fromstring('FaceCenter', 'c')
        # info[2].append(['GridLocation', v, [], 'GridLocation_t'])
        # if isinstance(faceList, numpy.ndarray): r = faceList
        # else: r = numpy.array(faceList, dtype=Internal.E_NpyInt)
        # info[2].append(['PointList', r, [], 'IndexArray_t'])
        # v = numpy.fromstring('Overset', 'c')
        # info[2].append(['GridConnectivityType', v, [], 'GridConnectivityType_t'])

        # d = numpy.ones((3,1), Internal.E_NpyInt)
        # c = numpy.ones((3,1), numpy.float64)
        # info[2].append(['PointListDonor', d, [], 'IndexArray_t'])
        # info[2].append(['InterpolantsDonor', c, [], 'DataArray_t'])


    elif (bndType1 == 'FamilySpecified' and fnmatch.fnmatch(bndType2, 'BCStage*')) or (bndType1 == 'BCStage'):
        _addFamilyOfStageGC__(z, bndName, bndType2, typeZone=1, faceList=faceList, zoneDonor=zoneDonor)

    else: # BC classiques
        # Cree le noeud zoneBC si besoin
        zoneBC = Internal.getNodesFromType1(z, 'ZoneBC_t')
        if zoneBC == []:
            z[2].append(['ZoneBC', None, [], 'ZoneBC_t'])
            zoneBC = z[2][len(z[2])-1]
        else: zoneBC = zoneBC[0]

        # Cree le noeud de BC
        Internal._createChild(zoneBC, bndName, 'BC_t', value=bndType1)
        l = len(zoneBC[2])
        info = zoneBC[2][l-1]
        Internal.createChild(info, 'GridLocation', 'GridLocation_t', 'FaceCenter')
        if isinstance(faceList, numpy.ndarray): r = faceList
        else: r = numpy.array(faceList, dtype=Internal.E_NpyInt)
        r = r.reshape((1,r.size), order='F')
        info[2].append([Internal.__FACELIST__, r, [], 'IndexArray_t'])

        # Ajoute la famille si necessaire
        if bndType1 == 'FamilySpecified':
            Internal.createChild(info, 'FamilyName', 'FamilyName_t', bndType2)

        # Ajoute les Data si necessaire (donnees Dirichlet)
        if data is not None:
            node1 = Internal.createNode('State', 'DataArray_t', value=data)
            node2 = Internal.createNode('DirichletData', 'BCData_t', children=[node1])
            node3 = Internal.createNode('BCDataSet', 'BCDataSet_t', children=[node2])
            info[2].append(node3)
    return None

# -- addBC2Zone pour les zones non structurees a elements basiques
def _addBC2UnstructZone__(z, bndName, bndType, elementList, elementRange,
                          faceList, pointList, data, subzone,
                          zoneDonor, elementListDonor, elementRangeDonor, faceListDonor,
                          rotationCenter, rotationAngle, translation, unitAngle=None):
    s = bndType.split(':')
    bndType1 = s[0]
    if len(s) > 1: bndType2 = s[1]
    else: bndType2 = ''

    # si subzone: on cree la connectivite BC, on passe en elementRange
    if subzone is not None and pointList is None:
        bcn = Internal.getNodeFromName1(z, subzone[0])
        if bcn is None:
            _mergeConnectivity(z, subzone, boundary=1)
        bcn = Internal.getNodeFromName1(z, subzone[0])
        bcnr = Internal.getNodeFromName1(bcn, 'ElementRange')
        elementRange = [bcnr[1][0], bcnr[1][1]]

    # si subzone + pointList=True, on identifie la subzone en pointList
    if subzone is not None and pointList == True:
        hook = createHook(z, function='nodes')
        pointList = identifyNodes(hook, subzone)
        freeHook(hook)

    if bndType1 == 'BCMatch' or bndType1 == 'Abutting1to1':
        if (zoneDonor == [] or
            faceListDonor is None and subzone is None and elementListDonor is None and elementRangeDonor is None):
            raise ValueError("addBC2Zone: unstructured match connectivity requires a donor zone and a faceListDonor or a subzone or an elementRangeDonor or an elementListDonor.")
        # si subzone fournie: cree le elementRangeDonor
        if subzone is not None:
            bcn = Internal.getNodeFromName1(zoneDonor, subzone[0])
            if bcn is None:
                _mergeConnectivity(zoneDonor, subzone, boundary=1)
            bcn = Internal.getNodeFromName1(zoneDonor, subzone[0])
            bcnr = Internal.getNodeFromName1(bcn, 'ElementRange')
            elementRangeDonor = [bcnr[1][0], bcnr[1][1]]

        # Cree le noeud zoneGridConnectivity si besoin
        zoneGC = Internal.createUniqueChild(z, 'ZoneGridConnectivity',
                                            'ZoneGridConnectivity_t')

        if isinstance(zoneDonor, str): v = zoneDonor
        else: v = zoneDonor[0]
        Internal._createChild(zoneGC, bndName, 'GridConnectivity_t', value=v)
        l = len(zoneGC[2])
        info = zoneGC[2][l-1]

        if elementList is not None:
            if isinstance(elementList, numpy.ndarray): r = elementList
            else: r = numpy.array(elementList, dtype=Internal.E_NpyInt)
            r = r.reshape((1,r.size), order='F')
            info[2].append([Internal.__ELEMENTLIST__, r, [], 'IndexArray_t'])
        elif elementRange != []:
            r = numpy.empty((1,2), dtype=Internal.E_NpyInt, order='F')
            r[0,0] = elementRange[0]
            r[0,1] = elementRange[1]
            info[2].append([Internal.__ELEMENTRANGE__, r, [], 'IndexRange_t'])
        elif faceList is not None:
            Internal._createChild(info, 'GridLocation', 'GridLocation_t', value='FaceCenter')
            if isinstance(faceList, numpy.ndarray): r = faceList
            else: r = numpy.array(faceList, dtype=Internal.E_NpyInt)
            r = r.reshape((1,r.size), order='F')
            info[2].append([Internal.__FACELIST__, r, [], 'IndexArray_t'])

        if elementListDonor is not None:
            if isinstance(elementListDonor, numpy.ndarray):
                r = elementListDonor
            else: r = numpy.array(elementListDonor, dtype=Internal.E_NpyInt)
            r = r.reshape((1,r.size))
            info[2].append([Internal.__ELEMENTLIST__+'Donor', r, [], 'IndexArray_t'])
        elif elementRangeDonor is not None:
            r = numpy.empty((1,2), dtype=Internal.E_NpyInt, order='F')
            r[0,0] = elementRangeDonor[0]
            r[0,1] = elementRangeDonor[1]
            info[2].append([Internal.__ELEMENTRANGE__+'Donor', r, [], 'IndexRange_t'])
        elif faceListDonor is not None:
            if isinstance(faceListDonor, numpy.ndarray): r = faceList
            else: r = numpy.array(faceListDonor, dtype=Internal.E_NpyInt)
            r = r.reshape((1,r.size), order='F')
            info[2].append([Internal.__FACELIST__+'Donor', r, [], 'IndexArray_t'])

        Internal.createChild(info, 'GridConnectivityType', 'GridConnectivityType_t', 'Abutting1to1')

        # Ajout pour les BC match periodiques
        _addPeriodicInfoInGC__(info, rotationCenter, rotationAngle, translation, unitAngle=unitAngle)

    elif bndType1 == 'BCOverlap':
        print('Warning: addBC2Zone: BCOverlap not valid for unstructured zones.')

    elif bndType1 == 'BCNearMatch':
        print('Warning: addBC2Zone: BCNearMatch not valid for unstructured zones.')

    elif (bndType1 == 'FamilySpecified' and fnmatch.fnmatch(bndType2, 'BCStage*')) or (bndType1 == 'BCStage'):
        _addFamilyOfStageGC__(z, bndName, bndType2, typeZone=2, elementRange=elementRange,
                              elementList=elementList, faceList=faceList, zoneDonor=zoneDonor)

    else: # BC classique
        # Cree le noeud zoneBC si besoin
        zoneBC = Internal.createUniqueChild(z, 'ZoneBC', 'ZoneBC_t')

        # Cree le noeud de BC
        info = Internal.createChild(zoneBC, bndName, 'BC_t', value=bndType1)
        if elementList is not None:
            if isinstance(elementList, numpy.ndarray): r = elementList
            else: r = numpy.array(elementList, dtype=Internal.E_NpyInt)
            r = r.reshape((1,r.size), order='F')
            Internal.createChild(info, Internal.__ELEMENTLIST__, 'IndexArray_t', value=r)
        elif elementRange != []:
            n = numpy.empty((1,2), dtype=Internal.E_NpyInt, order='F')
            n[0,0] = elementRange[0]
            n[0,1] = elementRange[1]
            Internal.createUniqueChild(info, Internal.__ELEMENTRANGE__,
                                       'IndexRange_t', value=n)
        elif faceList is not None:
            Internal.createUniqueChild(info, 'GridLocation', 'GridLocation_t',
                                       value='FaceCenter')
            if isinstance(faceList, numpy.ndarray): r = faceList
            else: r = numpy.array(faceList, dtype=Internal.E_NpyInt)
            r = r.reshape((1,r.size), order='F')
            info[2].append([Internal.__FACELIST__, r, [], 'IndexArray_t'])

        elif pointList is not None:
            Internal.createUniqueChild(info, 'GridLocation', 'GridLocation_t',
                                       value='Vertex')
            if isinstance(pointList, numpy.ndarray): r = pointList
            else: r = numpy.array(pointList, dtype=Internal.E_NpyInt)
            r = r.reshape((1,r.size), order='F')
            info[2].append([Internal.__POINTLIST__, r, [], 'IndexArray_t'])

        # Ajoute la famille si necessaire
        if bndType1 == 'FamilySpecified':
            Internal.createChild(info, 'FamilyName', 'FamilyName_t', bndType2)

        # Ajoute les Data si necessaire (donnees Dirichlet)
        if data is not None:
            node1 = Internal.createNode('State', 'DataArray_t', value=data)
            node2 = Internal.createNode('DirichletData', 'BCData_t', children=[node1])
            node3 = Internal.createNode('BCDataSet', 'BCDataSet_t', children=[node2])
            info[2].append(node3)
    return None

# -- converti une chaine en range (z Structure seulement)
# IN: 'imin', 'jmin', ...
# OUT: [imin,imax,jmin,jmax,kmin,kmax]
def convertStringRange2Range__(wrange, z):
    if wrange == 'doubly_defined': return wrange
    dim = Internal.getZoneDim(z)
    typeGrid = dim[0]
    if typeGrid == 'Unstructured':
        raise TypeError("addBC2Zone: range invalid for unstructured grids.")
    if wrange == 'imin': wrange = [1,1,1,dim[2],1,dim[3]]
    elif wrange == 'imax': wrange = [dim[1],dim[1],1,dim[2],1,dim[3]]
    elif wrange == 'jmin': wrange = [1,dim[1],1,1,1,dim[3]]
    elif wrange == 'jmax': wrange = [1,dim[1],dim[2],dim[2],1,dim[3]]
    elif wrange == 'kmin': wrange = [1,dim[1],1,dim[2],1,1]
    else: wrange = [1,dim[1],1,dim[2],dim[3],dim[3]]
    return wrange

# -- converti un BCRange en tableau de BCFaces
def convertBCRange2BCFace__(ni, nj, nk, range0):
    win = Internal.range2Window(range0)
    imin = win[0]; imax = win[1]
    jmin = win[2]; jmax = win[3]
    kmin = win[4]; kmax = win[5]

    dir = 0
    di = imax-imin; dj = jmax-jmin; dk = kmax-kmin
    if di == 0:
        di = 1
        if ni != 1: dir = 1
    if dj == 0:
        dj = 1
        if nj != 1: dir = 2
    if dk == 0:
        dk = 1
        if nk != 1: dir = 3

    if dir == 1: nfaces = dj*dk
    elif dir == 2: nfaces = di*dk
    else: nfaces = di*dj
    ni1 = max(1, ni-1); nj1 = max(1, nj-1); nk1 = max(1, nk-1)
    ninti = ni*nj1*nk1; nintj = ni1*nj*nk1
    bcfaces = numpy.empty((nfaces), dtype=Internal.E_NpyInt)
    ninj = ni*nj
    nof = 0
    if nk == 1: kmax += 1
    if nj == 1: jmax += 1
    if ni == 1: jmax += 1
    if dir == 1:
        for k in range(kmin,kmax):
            for j in range(jmin,jmax):
                indf = (imin-1) + (j-1)*ni+(k-1)*ni*nj1
                bcfaces[nof] = indf+1; nof+=1
    elif dir == 2:
        for k in range(kmin,kmax):
            for i in range(imin,imax):
                indf = (i-1) + (jmin-1)*ni1+(k-1)*ni1*nj
                bcfaces[nof] = indf+1+ninti; nof+=1
    else:
        for j in range(jmin,jmax):
            for i in range(imin,imax):
                indf = (i-1) + (j-1)*ni1+(kmin-1)*ni1*nj1
                bcfaces[nof] = indf+1+ninti+nintj; nof+=1
    return bcfaces

# -- Computes the nearmatch ratio between two near-matching grids (structured)
def getNMRatio__(win, winopp, trirac):
    i1 = win[0]; j1 = win[2]; k1 = win[4]
    i2 = win[1]; j2 = win[3]; k2 = win[5]
    i1o = i1; j1o = j1; k1o = k1; i2o = i2; j2o = j2; k2o = k2
    i1opp = winopp[0]; j1opp = winopp[2]; k1opp = winopp[4]
    i2opp = winopp[1]; j2opp = winopp[3]; k2opp = winopp[5]
    oi = trirac[0]
    if len(trirac) > 1: oj = trirac[1]
    else: oj = 2
    if len(trirac) > 2: ok = trirac[2]
    else: ok = 3
    if oi == 2 or oi == -2: j1o = i1; j2o = i2
    elif oi == 3 or oi == -3: k1o = i1; k2o = i2

    if oj == 1 or oj == -1: i1o = j1; i2o = j2
    elif oj == 3 or oj == -3: k1o = j1; k2o = j2

    if ok == 1 or ok == -1: i1o = k1; i2o = k2
    elif ok == 2 or ok == -2: j1o = k1; j2o = k2

    diopp = max(1,i2opp-i1opp); dio = max(1,i2o-i1o)
    djopp = max(1,j2opp-j1opp); djo = max(1,j2o-j1o)
    dkopp = max(1,k2opp-k1opp); dko = max(1,k2o-k1o)

    ri = float(diopp)/dio
    rj = float(djopp)/djo
    rk = float(dkopp)/dko
    res = [1.,1.,1.]

    if   oi == 2 or oi ==-2: res[0] = rj
    elif oi == 3 or oi ==-3: res[0] = rk
    elif oi == 1 or oi ==-1: res[0] = ri

    if   oj == 1 or oj==-1: res[1] = ri
    elif oj == 3 or oj==-3: res[1] = rk
    elif oj == 2 or oj==-2: res[1] = rj

    if   ok == 1 or ok ==-1: res[2] = ri
    elif ok == 2 or ok ==-2: res[2] = rj
    elif ok == 3 or ok ==-3: res[2] = rk
    return res

# -- recoverBCs
# Identifie des subzones comme BC Faces
# IN: liste de geometries de BC, liste de leur nom, liste de leur type
# OUT: a modifie avec des BCs ajoutees en BCFaces
def recoverBCs(a, T, tol=1.e-11):
    """Recover given BCs on a tree.
    Usage: recoverBCs(a, (BCs, BCNames, BCTypes), tol)"""
    tp = Internal.copyRef(a)
    _recoverBCs(tp, T, tol)
    return tp

def _recoverBCs(t, BCInfo, tol=1.e-11, removeBC=True):
    if removeBC: return _recoverBCs1(t, BCInfo, tol)
    else: return _recoverBCs2(t, BCInfo, tol)

# N'efface pas les matchs et bc deja existantes
def _recoverBCs2(t, BCInfo, tol):
    try: import Post.PyTree as P
    except: raise ImportError("_recoverBCs: requires Post module.")
    try: import Transform.PyTree as T
    except: raise ImportError("_recoverBCs: requires Transform module.")
    try: import Generator.PyTree as G
    except: raise ImportError("_recoverBCs: requires Generator module.")
    (BCs, BCNames, BCTypes) = BCInfo
    for z in Internal.getZones(t):
        indicesF = []
        zf = P.exteriorFaces(z, indices=indicesF)
        indicesF = indicesF[0]
        # BC classique
        bnds = Internal.getNodesFromType2(z, 'BC_t')
        # BC Match
        bnds += Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
        # BC Overlap/NearMatch/NoMatch
        bnds += Internal.getNodesFromType2(z, 'GridConnectivity_t')
        indicesBC = []
        for b in bnds:
            f = Internal.getNodeFromName1(b, 'PointList')
            indicesBC.append(f[1])

        undefBC = False
        if indicesBC != []:
            indicesBC = numpy.concatenate(indicesBC, axis=1)
            nfacesExt = indicesF.shape[0]
            nfacesDef = indicesBC.shape[1]
            if nfacesExt < nfacesDef:
                print('Warning: zone %s: number of faces defined by BCs is greater than the number of external faces. Try to reduce the matching tolerance.'%(z[0]))
            elif nfacesExt > nfacesDef:
                indicesBC = indicesBC.reshape( (indicesBC.size) )
                indicesE = Converter.converter.diffIndex(indicesF, indicesBC)
                undefBC = True
        else:
            undefBC = True
            indicesE = indicesF
        if undefBC:
            zf = T.subzone(z, indicesE, type='faces')
            hook = createHook(zf, 'elementCenters')
            for c in range(len(BCs)):
                if BCs[c] == []: raise ValueError("_recoverBCs: boundary is probably ill-defined.")
                for b in BCs[c]:
                    if G.bboxIntersection(zf, b):
                        # Break BC connectivity si necessaire
                        elts = Internal.getElementNodes(b)
                        size = 0
                        for e in elts:
                            erange = Internal.getNodeFromName1(e, 'ElementRange')[1]
                            size += erange[1]-erange[0]+1
                        n = len(elts)
                        if n == 1:
                            ids = identifyElements(hook, b, tol)
                        else:
                            bb = breakConnectivity(b)
                            ids = numpy.array([], dtype=Internal.E_NpyInt)
                            for bc in bb:
                                ids = numpy.concatenate([ids, identifyElements(hook, bc, tol)])

                        # Cree les BCs
                        ids0 = ids # keep ids for bcdata
                        ids = ids[ids > -1]
                        sizebc = ids.size
                        if sizebc > 0:
                            id2 = numpy.empty(sizebc, dtype=Internal.E_NpyInt)
                            id2[:] = indicesE[ids[:]-1]
                            _addBC2Zone(z, BCNames[c], BCTypes[c], faceList=id2)

                            # Recupere BCDataSets
                            fsc = Internal.getNodeFromName(b, Internal.__FlowSolutionCenters__)

                            if fsc is not None:
                                newNameOfBC = getLastBCName(BCNames[c])
                                bcz = Internal.getNodeFromNameAndType(z, newNameOfBC, 'BC_t')

                                ds = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                                                           gridLocation='FaceCenter', parent=bcz)
                                d = Internal.newBCData('NeumannData', parent=ds)

                                for node in Internal.getChildren(fsc):
                                    if Internal.isType(node, 'DataArray_t'):
                                        val0 = Internal.getValue(node)
                                        if isinstance(val0,numpy.ndarray):
                                            val0 = numpy.reshape(val0, val0.size, order='F')
                                        else:
                                            val0 = numpy.array([val0])
                                        val1 = val0[ids0>-1]
                                        Internal._createUniqueChild(d, node[0], 'DataArray_t', value=val1)
            freeHook(hook)
    return None

def _recoverBCs1(a, T, tol=1.e-11):
    """Recover given BCs on a tree.
    Usage: _recoverBCs(a, (BCs, BCNames, BCTypes), tol)"""
    try:import Post.PyTree as P
    except: raise ImportError("_recoverBCs: requires Post module.")
    _deleteZoneBC__(a)
    zones = Internal.getZones(a)
    (BCs, BCNames, BCTypes) = T
    for z in zones:
        indicesF = []
        try: f = P.exteriorFaces(z, indices=indicesF)
        except: continue
        indicesF = indicesF[0]
        hook = createHook(f, 'elementCenters')

        for c in range(len(BCs)):
            b = BCs[c]

            if b == []:
                raise ValueError("_recoverBCs: boundary is probably ill-defined.")
            # Break BC connectivity si necessaire
            elts = Internal.getElementNodes(b)
            size = 0
            for e in elts:
                erange = Internal.getNodeFromName1(e, 'ElementRange')[1]
                size += erange[1]-erange[0]+1
            n = len(elts)
            if n == 1:
                ids = identifyElements(hook, b, tol)
            else:
                bb = breakConnectivity(b)
                ids = numpy.array([], dtype=Internal.E_NpyInt)
                for bc in bb:
                    ids = numpy.concatenate([ids, identifyElements(hook, bc, tol)])

            # Cree les BCs
            ids0 = ids # keep ids for bcdata
            ids  = ids[ids > -1]
            sizebc = ids.size
            if sizebc > 0:
                id2 = numpy.empty(sizebc, dtype=Internal.E_NpyInt)
                id2[:] = indicesF[ids[:]-1]
                _addBC2Zone(z, BCNames[c], BCTypes[c], faceList=id2)

                # Recupere BCDataSets
                fsc = Internal.getNodeFromName(b, Internal.__FlowSolutionCenters__)

                if fsc is not None:
                    newNameOfBC = getLastBCName(BCNames[c])
                    bcz = Internal.getNodeFromNameAndType(z, newNameOfBC, 'BC_t')

                    ds = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                                             gridLocation='FaceCenter', parent=bcz)
                    d = Internal.newBCData('NeumannData', parent=ds)

                    for node in Internal.getChildren(fsc):
                        if Internal.isType(node, 'DataArray_t'):
                            val0 = Internal.getValue(node)
                            if isinstance(val0,numpy.ndarray):
                                val0 = numpy.reshape(val0, val0.size, order='F')
                            else:
                                val0 = numpy.array([val0])
                            val1 = val0[ids0>-1]
                            Internal._createUniqueChild(d, node[0], 'DataArray_t', value=val1)

        freeHook(hook)

    return None

# -- pushBC
# Recopie les BCs de z1 sur z2 (zones a zones)
# par identification geometrique. Seules les BCs correspondant a des faces
# identiques geometriquement seront reportees.
# IN: z1: any type (STRUCT,BE,NGON)
# IN: z2 tout sauf STRUCT
# IN: type='F': pour tout type: output en face list
#     type='BCC': seulement si z2 est BE: output en BCC (BC connect)
# IN: overwriteBC: efface les BC de t2 si existent
def pushBC(t1, t2, type='F', tol=1.e-12, overwriteBC=True):
    """Put BCs of t1 to t2.
    Usage: t2 = pushBC(t1, t2)"""
    zones1 = Internal.getZones(t1)
    t2p = Internal.copyRef(t2)
    zones2 = Internal.getZones(t2p)
    nz = min(len(zones1),len(zones2))
    for c in range(nz):
        z1 = zones1[c]; zp = zones2[c]
        BCs = Internal.getNodesFromType2(z1, 'BC_t')
        if len(BCs) != 0:
            if overwriteBC: _deleteZoneBC__(zp)
            dims = Internal.getZoneDim(zp)
            outBCC = 0
            if dims[0] == 'Unstructured' and dims[3] != 'NGON' and type != 'F':
                outBCC = 1
            if outBCC == 0: KDT = createHook(zp, function='faceCenters')

            for b in BCs:
                name = b[0]
                BCType = Internal.getValue(b)
                if BCType == 'FamilySpecified':
                    familyName = Internal.getNodeFromType1(b, 'FamilyName_t')
                    if familyName is not None: BCType = 'FamilySpecified:%s'%Internal.getValue(familyName)
                ext = extractBCOfName(z1, name)
                if outBCC == 0:
                    n = identifyElements(KDT, ext, tol)
                    n = n[n>0]
                    if n.size > 0: _addBC2Zone(zp, name, BCType, faceList=n)
                else: # output BCC
                    try:
                        import Transform.PyTree as T
                        ext = T.breakElements(ext[0])
                    except: pass
                    valid = []
                    for e in ext:
                        d = Internal.getZoneDim(e)
                        if d[0] == 'Structured': valid.append(convertArray2Hexa(e))
                        elif d[0] == 'Unstructured' and d[3] != 'NGON': valid.append(e)
                    if len(valid) == 1:
                        _addBC2Zone(zp, name, BCType, subzone=valid[0])
                    elif len(valid) > 1:
                        f = 0
                        for v in valid: _addBC2Zone(zp, name+str(f), BCType,
                                                    subzone=v); f += 1
            if outBCC == 0: freeHook(KDT)
    return t2p

def identifyBC(t, infos, tol=1.e-12):
    if infos == []: return t
    try: import Post.PyTree as P
    except: raise ImportError("identifyBC: requires Post.PyTree module.")
    allWins = []
    for info in infos: allWins.append(Converter.node2Center(info[3]))

    # Creation d'un hook global a partir de toutes les fenetres
    indirWins = [] # permet d'identifier la fenetre a laquelle se rapporte un pt du hook
    globalHook, indirWins = Converter.createGlobalHook(allWins, 'nodes', indir=1)
    # Identify and gather...
    tpp,typen = Internal.node2PyTree(t)
    for nob in range(len(tpp[2])):
        if tpp[2][nob][3] == 'CGNSBase_t':
            for noz in range(len(tpp[2][nob][2])):
                if tpp[2][nob][2][noz][3] == 'Zone_t':
                    z = tpp[2][nob][2][noz]
                    dimZ = Internal.getZoneDim(z)
                    if dimZ[0] == 'Structured':
                        niZ = dimZ[1]; njZ = dimZ[2]; nkZ = dimZ[3]
                        faces = P.exteriorFacesStructured(z)
                        dirf = 0
                        for face in faces:
                            dirf += 1
                            dimW = Internal.getZoneDim(face)
                            niw = dimW[1]; njw = dimW[2]
                            if niw == 1 and njw == 1: pass
                            else:
                                face = node2Center(face)
                                dimW = Internal.getZoneDim(face)
                                niw = dimW[1]; njw = dimW[2]
                                idHook = identifyNodes(globalHook, face, tol)
                                if max(idHook) > -1:
                                    ranges,nowins = Internal.gatherInStructPatch2D__(idHook,indirWins,niw,njw,dirf,niZ,njZ,nkZ)
                                    if ranges != []:
                                        nor = 0
                                        for r in ranges:
                                            noinfo = nowins[nor]
                                            info = infos[noinfo]
                                            bcname = info[0]
                                            bctype = info[1]
                                            bcdata = info[2]
                                            if bcdata is not None: bcdata = bcdata[1]
                                            win = info[3]
                                            famName = ''
                                            if bctype.split(':')[0] == 'FamilySpecified': famName=bctype.split(':')[1]

                                            # Passage en noeuds des indices
                                            istart=r[0]; iend=r[1]
                                            jstart=r[2]; jend=r[3]
                                            kstart=r[4]; kend=r[5]
                                            if istart != iend:
                                                if iend>1: iend=min(niZ,iend+1)
                                            if jstart != jend:
                                                if jend>1: jend=min(njZ,jend+1)
                                            if kstart != kend:
                                                if kend>1: kend=min(nkZ,kend+1)
                                            else: # cas 2D
                                                if nkZ == 2 and kend == 1: kend = 2
                                            nor += 1

                                            # addBC
                                            _addBC2Zone(z,bcname,bctype,[istart,iend,jstart,jend,kstart,kend],data=bcdata)
                                            tpp[2][nob][2][noz] = z
    Converter.freeHook(globalHook)
    tp = Internal.pyTree2Node(tpp, typen)
    return tp

# -- tagDefinedBC
# tag des points definis par une CL:
# tag=0 si pas de BC
#    =1 point frontiere
#    =2 point interieur
def tagDefinedBC(t):
    a = Internal.copyRef(t)
    # si a = arbre
    toptree = Internal.isTopTree(a)
    if toptree:
        bases = Internal.getBases(a)
        for b in bases:
            zones = Internal.getNodesFromType1(b, 'Zone_t')
            for z in zones:
                dim = Internal.getZoneDim(z)[4]
                p = Internal.getParentOfNode(b, z)
                if dim == 3: b[2][p[1]] = tagDefinedBCForZone3D__(z)
                else: b[2][p[1]] = tagDefinedBCForZone2D__(z)
    else:
        # si a = base
        base = Internal.getBases(a)
        if base != []:
            zones = Internal.getZones(a)
            for z in zones:
                dim = Internal.getZoneDim(z)[4]
                p = Internal.getParentOfNode(a, z)
                if dim == 3: a[2][p[1]] = tagDefinedBCForZone3D__(z)
                else: a[2][p[1]] = tagDefinedBCForZone2D__(z)
        else:
            # liste zones ou zone ?
            stdNode = Internal.isStdNode(a)
            if stdNode == 0 : # liste de zones
                zones = Internal.getZones(a)
                nzones = len(zones)
                for noz in range(nzones):
                    z = zones[noz]
                    dim = Internal.getZoneDim(z)[4]
                    if dim == 3: a[noz] = tagDefinedBCForZone3D__(z)
                    else: zones[noz] = tagDefinedBCForZone2D__(z)
            else:
                dim = Internal.getZoneDim(a)[4]
                if dim == 3: a = tagDefinedBCForZone3D__(a)
                else: a = tagDefinedBCForZone2D__(a)

    return a

def tagDefinedBCForZone2D__(z):
    dims = Internal.getZoneDim(z)
    if dims[0] == 'Unstructured': return z
    ni = dims[1]; nj = dims[2]; nk = dims[3]
    tag = Converter.array('definedBC',ni,nj,1)
    taga = tag[1][0]; wins = []
    # BC classique
    bnds = Internal.getNodesFromType2(z, 'BC_t')
    for bc in bnds:
        wrange = Internal.getNodeFromName1(bc, 'PointRange')
        r = wrange[1]
        wins.append(Internal.range2Window(r))
    # BCNearMatch
    bnds = Internal.getNodesFromType2(z, 'GridConnectivity_t')
    for bc in bnds:
        type = Internal.getNodeFromName1(bc, 'GridConnectivityType')
        if type is not None:
            val = Internal.getValue(type)
            if val == 'Abutting':
                wrange = Internal.getNodeFromName1(bc, 'PointRange')
                r = wrange[1]
                wins.append(Internal.range2Window(r))

    # BC match
    bnds = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
    for bc in bnds:
        gl = Internal.getNodeFromType1(bc,'GridLocation_t')
        isHybrid = False
        if gl is not None:
            gl = Internal.getValue(gl)
            if gl == 'FaceCenter' or gl == 'CellCenter': isHybrid=True
        if not isHybrid:
            wrange = Internal.getNodeFromName1(bc, 'PointRange')
            r = wrange[1]
            wins.append(Internal.range2Window(r))
    # BC Overlap
    bnds = Internal.getNodesFromType1(z, 'GridConnectivity_t')
    for bc in bnds:
        r = Internal.getNodesFromName1(bc, 'GridConnectivityType')
        if r is not None:
            val = Internal.getValue(r)
            if val == 'Overset':
                dd = 0 # doubly defined
                userDef = Internal.getNodesFromName(bc, 'UserDefinedData')
                if userDef != []:
                    if len(userDef[0]) == 4:
                        info = userDef[0][2][0]
                        if info[0] == 'doubly_defined': dd = 1
                if dd == 0:
                    wrange = Internal.getNodeFromName(bc, 'PointRange')
                    wins.append(Internal.range2Window(wrange[1]))

    tag = Converter.converter.tagDefinedBC(tag, wins, 2)
    z = setFields([tag], z, 'nodes')
    return z

def tagDefinedBCForZone3D__(z):
    dims = Internal.getZoneDim(z)
    if dims[0] == 'Unstructured': return z
    ni = dims[1]; nj = dims[2]; nk = dims[3]

    tag = Converter.array('definedBC', ni, nj, nk)
    wins = []
    # BC classique
    bnds = Internal.getNodesFromType2(z, 'BC_t')
    for bc in bnds:
        wrange = Internal.getNodeFromName1(bc, 'PointRange')
        r = wrange[1]
        wins.append(Internal.range2Window(r))
    # BCNearMatch
    bnds = Internal.getNodesFromType2(z, 'GridConnectivity_t')
    for bc in bnds:
        type = Internal.getNodeFromName1(bc, 'GridConnectivityType')
        if type is not None:
            val = Internal.getValue(type)
            if val == 'Abutting':
                wrange = Internal.getNodeFromName1(bc, 'PointRange')
                r = wrange[1]
                wins.append(Internal.range2Window(r))
    # BC match
    bnds = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
    for bc in bnds:
        gl = Internal.getNodeFromType1(bc, 'GridLocation_t')
        isHybrid = False
        if gl is not None:
            gl = Internal.getValue(gl)
            if gl == 'FaceCenter' or gl == 'CellCenter': isHybrid=True
        if not isHybrid:
            wrange = Internal.getNodeFromName1(bc, 'PointRange')
            r = wrange[1]
            wins.append(Internal.range2Window(r))
    # BC Overlap
    bnds = Internal.getNodesFromType2(z, 'GridConnectivity_t')
    for bc in bnds:
        r = Internal.getNodeFromName1(bc, 'GridConnectivityType')
        if r is not None:
            val = Internal.getValue(r)
            if val == 'Overset':
                dd = 0 # doubly defined
                userDef = Internal.getNodesFromName(bc, 'UserDefinedData')
                if userDef != []:
                    if len(userDef[0]) == 4:
                        info = userDef[0][2][0]
                        if info[0] == 'doubly_defined': dd = 1
                if dd == 0:
                    wrange = Internal.getNodeFromName1(bc, 'PointRange')
                    wins.append(Internal.range2Window(wrange[1]))
    tag = Converter.converter.tagDefinedBC(tag, wins, 3)
    z = setFields([tag], z, 'nodes')
    return z

# -- updateDefinedBCForWins__
# update the definedBC field for a list of windows wins defined by BC
# wins are [i1,i2,j1,j2,k1,k2]
def updateDefinedBCForWins__(z, wins, dim=3):
    tag = getField('definedBC', z)[0]
    for now in range(len(wins)):
        if len(wins[now]) < 6:
            for i in range(len(wins[now]),7): wins[now].append(1)
    tag = Converter.converter.tagDefinedBC(tag, wins, dim)
    z = setFields([tag], z, 'nodes')
    return z

#=============================================================================
# -- Remove BCs of type or name from a tree or a zone --
#=============================================================================

# -- rmBCOfType
def rmBCOfType(t, bndType):
    """Remove BCs of given type.
    Usage: rmBCOfType(t, bndType)"""
    tp = Internal.copyRef(t)
    _rmBCOfType(tp, bndType)
    return tp

def _rmBCOfType(t, bndType):
    names = bndType.split(':')
    if len(names) == 2: # Family
        _rmBCOfName(t, bndType)

    zones = Internal.getZones(t)
    if bndType == 'BCMatch':
        for z in zones:
            nodes = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
            for i in nodes:
                (parent, d) = Internal.getParentOfNode(z, i)
                del parent[2][d]
            nodes = Internal.getNodesFromType2(z, 'GridConnectivity_t')
            for i in nodes:
                r = Internal.getNodeFromType1(i, 'GridConnectivityType_t')
                if r is not None:
                    val = Internal.getValue(r)
                    if val == 'Abutting1to1':
                        (parent, d) = Internal.getParentOfNode(z, i)
                        del parent[2][d]
    elif bndType == 'BCNearMatch' or bndType == 'BCStage':
        for z in zones:
            nodes = Internal.getNodesFromType2(z, 'GridConnectivity_t')
            for i in nodes:
                r = Internal.getNodeFromType1(i, 'GridConnectivityType_t')
                if r is not None:
                    val = Internal.getValue(r)
                    if val == 'Abutting':
                        (parent, d) = Internal.getParentOfNode(z, i)
                        del parent[2][d]
    elif bndType == 'BCOverlap':
        for z in zones:
            nodes = Internal.getNodesFromType2(z, 'GridConnectivity_t')
            for i in nodes:
                r = Internal.getNodeFromType1(i, 'GridConnectivityType_t')
                if r is not None:
                    val = Internal.getValue(r)
                    if val == 'Overset':
                        (parent, d) = Internal.getParentOfNode(z, i)
                        del parent[2][d]
    else: # physical BC
        families = getFamilyBCNamesOfType(t, bndType)
        if bndType == 'BCWall':
            familiesI = getFamilyBCNamesOfType(t, 'BCWallInviscid')
            familiesV = getFamilyBCNamesOfType(t, 'BCWallViscous*')
        for z in zones:
            BCnodes = Internal.getNodesFromType1(z, 'ZoneBC_t')
            for n in BCnodes:
                nodes = Internal.getNodesFromValue(n, bndType)
                nodes += getFamilyBCs(z, families)
                if bndType == 'BCWall':
                    nodes += Internal.getNodesFromValue(n, 'BCWallInviscid')
                    nodes += Internal.getNodesFromValue(n, 'BCWallViscous*')
                    nodes += getFamilyBCs(z, familiesI)
                    nodes += getFamilyBCs(z, familiesV)
                for i in nodes:
                    (parent, d) = Internal.getParentOfNode(z, i)
                    del parent[2][d]
    return None

# -- rmBCOfName
def rmBCOfName(t, bndName):
    """Remove BCs of given name.
    Usage: rmBCOfName(t, bndName)"""
    tp = Internal.copyRef(t)
    _rmBCOfName(tp, bndName)
    return tp

def _rmBCOfName(t, bndName):
    names = bndName.split(':')
    if len(names) == 1: # real bnd name
        zones = Internal.getZones(t)
        for z in zones:
            nodes = Internal.getNodesFromName(z, bndName)
            for i in nodes:
                if Internal.getType(i) in ['BC_t', 'GridConnectivity1to1_t', 'GridConnectivity_t']:
                    (parent, d) = Internal.getParentOfNode(z, i)
                    del parent[2][d]
    else: # family specified BC
        zones = Internal.getZones(t)
        for z in zones:
            nodes = getFamilyBCs(z, names[1])
            for i in nodes:
                (parent, d) = Internal.getParentOfNode(z, i)
                del parent[2][d]
    return None

#==============================================================================
# -- Get empty BCs from a zone --
#==============================================================================

# -- getEmptyBC
# Return the list of empty BCs:
# return range (structured) or face list (unstructured) of undefined BCs
# for any zone in any bases
# if t=tree: returns [winsBase1, winsBase2,...],
# with winsBase=[winzone1Base1, winzone2Base1,...]
# if all the BCs are defined for a zone, [] is returned for the zone
# IN: dim: dimension of the pb
# IN: splitFactor: used only for unstructured grids, split the windows
# following split angle.
def getEmptyBC(a, dim=3, splitFactor=181.):
    """Return the range or facelist of unset boundary conditions."""
    # si a = arbre
    toptree = Internal.isTopTree(a)
    if toptree:
        winst = []
        bases = Internal.getBases(a)
        for b in bases:
            winsb = []
            zones = Internal.getNodesFromType1(b, 'Zone_t')
            for z in zones:
                winsz = getEmptyBCForZone__(z, dim, splitFactor)
                winsb.append(winsz)
            winst.append(winsb)
        return winst
    else:
        # si a = base
        base = Internal.getBases(a)
        if base != []:
            winsb = []
            zones = Internal.getNodesFromType1(base, 'Zone_t')
            for z in zones:
                winsz = getEmptyBCForZone__(z, dim, splitFactor)
                winsb.append(winsz)
            return winsb
        else:
            # liste zones ou zone ?
            stdNode = Internal.isStdNode(a)
            if stdNode == 0: # liste de zones
                wins = []
                for z in a[0:]:
                    winsz = getEmptyBCForZone__(z, dim, splitFactor)
                    wins.append(winsz)
                return wins
            else:
                return getEmptyBCForZone__(a, dim, splitFactor)
    return []

# -- Detect empty boundary conditions for zones
def getEmptyBCForZone__(z, pbDim, splitFactor):
    dims = Internal.getZoneDim(z)
    if dims[0] == 'Structured':
        return getEmptyBCForStructZone__(z, dims, pbDim, splitFactor)
    elif dims[3] == 'NGON':
        return getEmptyBCForNGonZone__(z, dims, pbDim, splitFactor)
    else: return getEmptyBCForBEZone__(z, dims, pbDim, splitFactor)

# -- Detect empty boundary conditions for struct zones
# Renvoie une liste de ranges des BC empty
def getEmptyBCForStructZone__(z, dims, pbDim, splitFactor):
    ni = dims[1]; nj = dims[2]; nk = dims[3]
    ni1 = ni-1; nj1 = nj-1; nk1 = nk-1
    nwins = []; wins = []

    # BC classique
    bnds = Internal.getNodesFromType2(z, 'BC_t')
    for bc in bnds:
        r = Internal.getNodeFromName1(bc, 'PointRange')
        if r is not None: wins.append(Internal.range2Window(r[1]))
    # BC match
    bnds = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
    for bc in bnds:
        r = Internal.getNodeFromName1(bc, 'PointRange')
        if r is not None: wins.append(Internal.range2Window(r[1]))
        r = Internal.getNodeFromName1(bc, 'PointList') # hybrid
        if r is not None:
            ret = Converter.converter.pointList2Ranges(r[1],ni,nj,nk)
            wins += ret

    # BC Overlap/NearMatch/NoMatch
    bnds = Internal.getNodesFromType2(z, 'GridConnectivity_t')
    for bc in bnds:
        r = Internal.getNodeFromName1(bc, 'PointRange')
        if r is not None: wins.append(Internal.range2Window(r[1]))

    # Parcours des faces
    directions = []
    if ni != 1: directions += [1,2]
    if nj != 1: directions += [3,4]
    if nk != 1 and pbDim == 3: directions += [5,6]
    for dir in directions:
        nwins = Converter.converter.detectEmptyBC(wins, ni, nj, nk, dir, nwins)
    return nwins

# -- Detect empty boundary conditions for unstruct zones
# Renvoie une liste des indices des faces des BC empty
def getEmptyBCForBEZone__(z, dims, pbDim, splitFactor):
    try: import Transform.PyTree as T; import Post.PyTree as P
    except: raise ImportError("getEmptyBC: requires Transform and Post modules.")
    indicesF = []
    f = P.exteriorFaces(z, indices=indicesF)
    indicesF = indicesF[0]
    _initVars(f, 'centers:__tag__', 1)

    # BC classique
    bnds = Internal.getNodesFromType2(z, 'BC_t')
    # BC Match
    bnds += Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
    # BC Overlap/NearMatch/NoMatch
    bnds += Internal.getNodesFromType2(z, 'GridConnectivity_t')

    zp = Internal.copyRef(z)
    _deleteZoneBC__(zp); _deleteFlowSolutions__(zp)

    defined = [] # BC deja definies
    for bc in bnds:
        flist = Internal.getNodeFromName1(bc, Internal.__FACELIST__)
        if flist is not None: defined.append(T.subzone(zp, flist[1], type='faces'))
        erange = Internal.getNodeFromName1(bc, Internal.__ELEMENTRANGE__)
        if erange is not None:
            r = erange[1]
            defined.append(selectOneConnectivity(zp, irange=[r[0,0],r[0,1]]))

    hook = createHook(f, 'elementCenters')
    if defined != []:
        tag = Internal.getNodeFromName2(f, '__tag__')[1]
        for d in defined:
            d = convertArray2NGon(d)
            id0 = identifyElements(hook, d)
            tag[id0[:]-1] = 0

    sel = P.selectCells2(f, 'centers:__tag__')
    if splitFactor >= 180.: sel = T.splitConnexity(sel)
    else: sel = T.splitSharpEdges(sel, alphaRef=splitFactor)

    id0 = []
    for s in sel:
        id1 = identifyElements(hook, s)
        #mask = (id1[:] >= 0) # enleve les elements non identifies
        id2 = numpy.empty(id1.size, dtype=Internal.E_NpyInt)
        id2[:] = indicesF[id1[:]-1]
        id0.append(id2)
    freeHook(hook)
    return id0

def getEmptyBCForNGonZone__(z, dims, pbDim, splitFactor):
    try: import Transform.PyTree as T; import Post.PyTree as P
    except: raise ImportError("getEmptyBC: requires Transform and Post modules.")
    indicesF=[]
    f = P.exteriorFaces(z, indices=indicesF)
    indicesF = indicesF[0]
    # BC classique
    bnds = Internal.getNodesFromType2(z, 'BC_t')
    # BC Match
    bnds += Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
    # BC Overlap/NearMatch/NoMatch
    bnds += Internal.getNodesFromType2(z, 'GridConnectivity_t')

    indicesBC = []
    for b in bnds:
        f = Internal.getNodeFromName1(b, 'PointList')
        indicesBC.append(f[1])

    undefBC = False
    if indicesBC != []:
        indicesBC = numpy.concatenate(indicesBC, axis=1)
        nfacesExt = indicesF.shape[0]
        nfacesDef = indicesBC.shape[1]
        if nfacesExt < nfacesDef:
            print('Warning: zone %s: number of faces defined by BCs is greater than the number of external faces. Try to reduce the matching tolerance.'%(z[0]))
        elif nfacesExt > nfacesDef:
            indicesBC = indicesBC.reshape( (indicesBC.size) )
            indicesE = Converter.converter.diffIndex(indicesF, indicesBC)
            #indicesE = numpy.setdiff1d(indicesF, indicesBC) #SP
            #indicesE = numpy.delete(indicesF,indicesBC[0,:])
            undefBC = True
    else:
        undefBC = True
        indicesE = indicesF
    if undefBC:
        ze = T.subzone(z, indicesE, type='faces')
        zea = getFields('GridCoordinates', ze)[0]
        id0 = T.transform.splitSharpEdgesList(zea, indicesE, splitFactor)
        return id0
    else: return []

#==============================================================================
# IN: z: zone subzone (structure)
# IN: w: window [imin,imax, jmin, jmax, kmin, kmax]
# Reorder unz subzone pour garder les normales exterieures
#==============================================================================
def _reorderSubzone__(z, w, T):
    dim = Internal.getZoneDim(z)
    ni = dim[1]; nj = dim[2]; nk = dim[3]
    imin = w[0]; imax = w[1]; jmin = w[2]; jmax = w[3]; kmin = w[4]; kmax = w[5]

    if imin == imax and imin > 1: T._reorder(z, (1,-2,3))
    elif jmin == jmax and jmin == 1: T._reorder(z, (-1,2,3))
    elif kmin == kmax and kmin > 1: T._reorder(z, (-1,2,3))
    return None

#==============================================================================
# -- Extract all BC of type or name --
#==============================================================================
# keep the bcdataset of bc from original zone in the subzone zbc
# extracted from the bc
def _keepBCDataSet(zbc, zorig, bc, extrapFlow=True):
    datas = Internal.getBCDataSet(zorig, bc)
    if datas == [] and not extrapFlow:
        Internal._rmNodesByName(zbc, Internal.__FlowSolutionCenters__)
    elif datas != []:
        if not extrapFlow:
            Internal._rmNodesByName(zbc, Internal.__FlowSolutionCenters__)
        f = Internal.createUniqueChild(zbc, Internal.__FlowSolutionCenters__,
                                        'FlowSolution_t')
        Internal.newGridLocation(value='CellCenter', parent=f)
        for d in datas: Internal.createUniqueChild(f, d[0], d[3], value=d[1])
    return None

# extract bc defined by zbc node as a subzone of z for BE zones
# with BCs defined by a PointList/FaceCenter refering to a connectivity
# of BC (TRI for TETRA for instance)
def getBC2__(zbc, z, T, res, extrapFlow=True):
    zdim = Internal.getZoneDim(z)
    if zdim[0] == 'Unstructured':
        ztype = zdim[3]
        if ztype != 'NGON':
            bndName = Internal.getName(zbc)
            bndType = Internal.getValue(zbc)
            gcl = Internal.getNodeFromType(zbc, 'GridLocation_t')
            PL = Internal.getNodeFromName(zbc, 'PointList')
            if PL is not None:
                if Internal.getValue(gcl) == 'FaceCenter':
                    PL = Internal.getValue(PL)
                    ermin = numpy.amin(PL)
                    ermax = numpy.amax(PL)
                    zp = Internal.copyRef(z)
                    connects = Internal.getNodesFromType(z,'Elements_t')
                    for cn in connects:
                        ERLoc = Internal.getNodeFromName(cn,'ElementRange')
                        ER_min = Internal.getValue(ERLoc)[0]
                        ER_max = Internal.getValue(ERLoc)[1]
                        if ermin >= ER_min and ermax <= ER_max:
                            for cn2 in connects:
                                if cn2[0] != cn[0]:
                                    Internal._rmNodesFromName(zp,cn2[0])
                            etype=Internal.getValue(cn)
                            if etype[1]!= 0: etype[1]=0
                            nfaces  = PL[0].shape[0]
                            facelist = numpy.zeros(nfaces, dtype=Internal.E_NpyInt)
                            for noet in range(nfaces): facelist[noet] = PL[0][noet]-ER_min
                            Internal._rmNodesFromType(zp,"ZoneBC_t")
                            Internal._rmNodesFromType(zp,"ZoneGridConnectivity_t")
                            Internal._rmNodesFromType(zp,"Family_t")
                            Internal._rmNodesFromType(zp,"FamilyName_t")
                            Internal._rmNodesFromName(zp,"ParentElement*")
                            zp=T.subzone(zp,facelist,type='elements')
                            zp[0] = z[0]+Internal.SEP1+zbc[0]
                            _keepBCDataSet(zp, z, zbc, extrapFlow=extrapFlow)
                            res.append(zp)
                            break
    return None

# Get the geometrical BC and append it to res
# IN: i: BC_t node or GridConnectivity_t node
# IN: z: zone owning BC
# IN: T: transform module
# IN: reorder: if True, reorder BC to get external normals
# IN: extrapFlow: if True, get the BCDataSet field
# IN: shift: if not 0, shift BC of index for structured grids only
def getBC__(i, z, T, res, reorder=True, extrapFlow=True, shift=0):

    connects = Internal.getNodesFromType1(z, "Elements_t")
    zdim = Internal.getZoneDim(z)
    if zdim[0] == 'Unstructured': ztype = zdim[3]

    # IndexRange (PointRange)
    r = Internal.getNodeFromType1(i, 'IndexRange_t')
    if r is not None and r[1].shape[0] > 1: # structure - suppose range in nodes
        wrange = r[1]
        w = Internal.range2Window(wrange)
        #zp = subzoneWithReorder__(z, w)
        imin = w[0]; imax = w[1]; jmin = w[2]; jmax = w[3]; kmin = w[4]; kmax = w[5]
        if imin == imax:
            if imin == 1: imin += shift; imax += shift
            else: imin -= shift; imax -= shift
        elif jmin == jmax:
            if jmin == 1: jmin += shift; jmax += shift
            else: jmin -= shift; jmax -= shift
        elif kmin == kmax:
            if kmin == 1: kmin += shift; kmax += shift
            else: kmin -= shift; kmax -= shift
        zp = T.subzone(z, (imin,jmin,kmin), (imax,jmax,kmax))
        zp[0] = z[0]+Internal.SEP1+i[0]
        # Get BCDataSet if any
        _keepBCDataSet(zp, z, i, extrapFlow=extrapFlow)

        if reorder: _reorderSubzone__(zp, w, T) # normales ext
        _deleteZoneBC__(zp)
        _deleteGridConnectivity__(zp)
        res.append(zp)
    elif r is not None and r[1].shape[0] == 1: # suppose BE + BCC
        r = r[1]
        # zp = selectOneConnectivity(z, irange=[r[0,0],r[0,1]])
        zp = selectConnectivity(z, irange=[r[0,0],r[0,1]])
        zp[0] = z[0]+Internal.SEP1+i[0]
        _deleteZoneBC__(zp)
        _deleteGridConnectivity__(zp)
        _deleteSolverNodes__(zp)

        # Get BCDataSet if any
        _keepBCDataSet(zp, z, i, extrapFlow=extrapFlow)
        res.append(zp)

    # IndexArray (PointList)
    if r is None: r = Internal.getNodeFromName(i, Internal.__FACELIST__)
    else: r = None

    if r is not None:
        loc = Internal.getNodeFromName1(i, 'GridLocation')
        if loc is not None:
            val = Internal.getValue(loc)
            if val == 'FaceCenter' or val == 'CellCenter': # Face list (BE ou NGON)
                faceList = r[1]
                rf = Internal.getElementRange(z, type='NGON')
                if rf is not None and rf[0] != 1: # decalage possible du NGON
                    faceList2 = numpy.copy(faceList)
                    faceList2[:] = faceList[:]-rf[0]+1
                    zp = T.subzone(z, faceList2, type='faces')
                else:
                    if len(connects)>1 and ztype != 'NGON':
                        return getBC2__(i, z, T, res, extrapFlow=extrapFlow)
                    else:
                        zp = T.subzone(z, faceList, type='faces') # BE

                zp[0] = z[0]+Internal.SEP1+i[0]
                _deleteZoneBC__(zp)
                _deleteGridConnectivity__(zp)
                _deleteSolverNodes__(zp)
                # Get BCDataSet if any
                _keepBCDataSet(zp, z, i, extrapFlow=extrapFlow)
                res.append(zp)
            elif val == 'Vertex': # vertex indices
                pointList = r[1]
                zp = T.subzone(z, pointList, type='nodes')
                zp[0] = z[0]+Internal.SEP1+i[0]
                _deleteZoneBC__(zp)
                _deleteGridConnectivity__(zp)
                _deleteSolverNodes__(zp)
                # Get BCDataSet if any
                _keepBCDataSet(zp, z, i, extrapFlow=extrapFlow)
                res.append(zp)
        else: # suppose FaceList
            faceList = r[1]
            rf = Internal.getElementRange(z, type='NGON')
            if rf is not None and rf[0] != 1:
                faceList2 = numpy.copy(faceList)
                faceList2[:] = faceList[:]-rf[0]+1
                zp = T.subzone(z, faceList2, type='faces')
            else: zp = T.subzone(z, faceList, type='faces')
            zp[0] = z[0]+Internal.SEP1+i[0]
            _deleteZoneBC__(zp)
            _deleteGridConnectivity__(zp)
            _deleteSolverNodes__(zp)
            # Get BCDataSet if any
            _keepBCDataSet(zp, z, i, extrapFlow=extrapFlow)
            res.append(zp)
    return None

# -- extractBCOfType
# Extract all BC of a given type
# Recognised bndType: classic CGNS + BCMatch + BCNearMatch + BCOverlap
# topTree: utile si t n'est pas un topTree, permet de trouver les familyBCs
def extractBCOfType(t, bndType, topTree=None, reorder=True, extrapFlow=True, shift=0):
    """Extract the grid coordinates of given BC type as zones."""
    try: import Transform.PyTree as T
    except: raise ImportError("extractBCOfType: requires Transform.PyTree module.")

    names = bndType.split(':')
    if len(names) == 2: # Family
        return extractBCOfName(t, bndType, reorder=reorder, extrapFlow=extrapFlow)

    zones = Internal.getZones(t)
    res = []
    if bndType == 'BCMatch':
        for z in zones:
                        # Cherche GridConnectivity1to1_t
            nodes = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
            for i in nodes: getBC__(i, z, T, res, reorder=reorder, extrapFlow=extrapFlow, shift=shift)
            # Cherche GridConnectivity_t + Abutting1to1
            nodes = Internal.getNodesFromType2(z, 'GridConnectivity_t')
            for i in nodes:
                type = Internal.getNodeFromName1(i, 'GridConnectivityType')
                if type is not None:
                    val = Internal.getValue(type)
                    if val == 'Abutting1to1': getBC__(i, z, T, res, reorder=reorder, extrapFlow=extrapFlow, shift=shift)
    elif bndType == 'BCNearMatch' or bndType == 'BCStage':
        for z in zones:
            nodes = Internal.getNodesFromType2(z, 'GridConnectivity_t')
            for i in nodes:
                type = Internal.getNodeFromName1(i, 'GridConnectivityType')
                if type is not None:
                    val = Internal.getValue(type)
                    if val == 'Abutting': getBC__(i, z, T, res, reorder=reorder, extrapFlow=extrapFlow, shift=shift)
    elif bndType == 'BCOverlap':
        for z in zones:
            nodes = Internal.getNodesFromType2(z, 'GridConnectivity_t')
            for i in nodes:
                r = Internal.getNodeFromName1(i, 'GridConnectivityType')
                if r is not None:
                    val = Internal.getValue(r)
                    if val == 'Overset': getBC__(i, z, T, res, reorder=reorder, extrapFlow=extrapFlow, shift=shift)
    else: # BC physiques
        if topTree is None: topTree = t
        families = getFamilyBCNamesOfType(topTree, bndType)
        if bndType == 'BCWall':
            families1 = getFamilyBCNamesOfType(topTree, 'BCWallInviscid')
            families2 = getFamilyBCNamesOfType(topTree, 'BCWallViscous*')

        for z in zones:
            nodes = []
            zbcs = Internal.getNodesFromType1(z, 'ZoneBC_t')
            for zbc in zbcs:
                for i in zbc[2]:
                    if Internal.isValue(i, bndType) and i[3] == 'BC_t': nodes.append(i)
            nodes += getFamilyBCs(z, families)
            if bndType == 'BCWall':
                for zbc in zbcs:
                    for i in zbc[2]:
                        if (Internal.isValue(i, 'BCWallInviscid') or Internal.isValue(i, 'BCWallViscous*')) and i[3] == 'BC_t':
                            nodes.append(i)
                nodes += getFamilyBCs(z, families1)
                nodes += getFamilyBCs(z, families2)
            for i in nodes: getBC__(i, z, T, res, reorder=reorder, extrapFlow=extrapFlow, shift=shift)
    return res

# -- extractBCOfName
def extractBCOfName(t, bndName, reorder=True, extrapFlow=True, shift=0):
    """Extract the grid coordinates of given BC name as zones.
    Usage: extractBCOfName(t, bndName)"""
    try: import Transform.PyTree as T
    except: raise ImportError("extractBCOfName: requires Transform.PyTree module.")
    names = bndName.split(':')
    res = []
    if len(names) == 1: # real bnd name
        for z in Internal.getZones(t):
        # Pas de niveau 2 car pas de wild card autorisee dans ce cas
            nodes = Internal.getNodesFromName(z, bndName)
            for i in nodes:
                if Internal.getType(i) in ['BC_t', 'GridConnectivity1to1_t', 'GridConnectivity_t']:
                    getBC__(i, z, T, res, reorder=reorder, extrapFlow=extrapFlow, shift=shift)
    else: # family specified BC
        for z in Internal.getZones(t):
            nodes = getFamilyBCs(z, names[1])
            for i in nodes: getBC__(i, z, T, res, reorder=reorder, extrapFlow=extrapFlow, shift=shift)
    return res

def extractBCOfSubRegionName__(t, zsrName, reorder=True, extrapFlow=True, shift=0):
    zones = []
    for z in Internal.getZones(t):
        zsr = Internal.getNodeFromName(z, zsrName)
        gl = Internal.getNodeFromType(zsr, 'GridLocation_t')
        bcname = Internal.getValue(Internal.getNodeFromName(zsr,'BCRegionName'))
        z_surf = extractBCOfName(z, bcname, reorder=reorder, extrapFlow=extrapFlow, shift=shift)[0]
        glname = 'Vertex'
        if Internal.getValue(gl) == 'FaceCenter': glname = 'CellCenter'

        Internal._rmNodesFromType(z_surf, 'FlowSolution_t')
        FS = Internal.newFlowSolution(name='FlowSolution#Centers', gridLocation=glname, parent=z_surf)
        for var in lvar: # BUGGED
            datan = Internal.getNodeFromName(zsr, var)
            FS[2].append(datan)

        Internal._rmNodesFromType(z_surf, "ZoneSubRegion_t")
        zones.append(z_surf)
        return zones

# -- getBCs
# Retourne la geometrie de toutes les BCs
# IN: t: une zone, une liste de zones, une base ou un arbre
# IN: extrapFlow: extrapolate flow solution if True
# IN: reorder: if True, extracted zones are reordered such that normals are oriented towards the interior of a
# OUT: BCs: liste des geometries de bcs
# OUT: BCNames: liste des noms des BCs
# OUT: BCTypes: liste des types des BCs
def getBCs(t, reorder=True, extrapFlow=True):
    """Return geometry, names and types of boundary conditions."""
    BCs = []; BCNames = []; BCTypes = []
    for z in Internal.getZones(t):
        nodes = Internal.getNodesFromType2(z, 'BC_t')
        for n in nodes:
            name = n[0]; typeBC = Internal.getValue(n)
            if typeBC == 'FamilySpecified':
                fname = Internal.getNodeFromType1(n, 'FamilyName_t')
                if fname is not None:
                    fname = Internal.getValue(fname)
                    typeBC = 'FamilySpecified:%s'%fname
            zBC = extractBCOfName(z, name, reorder, extrapFlow)
            if zBC == []:
                name2 = name.split(':')
                if len(name2)>1 and name2[0] == 'FamilySpecified':
                    raise ValueError("getBCs: BC of name FamilySpecified:* is not valid.")
            BCs.append(zBC); BCNames.append(name); BCTypes.append(typeBC)

        # Raccords Stage* definis par des familles
        nodes = Internal.getNodesFromType2(z, "ZoneGridConnectivity_t")
        for n in nodes:
            for gc in Internal.getNodesFromType1(n, "GridConnectivity_t"):
                name = Internal.getName(gc)
                fname = Internal.getNodeFromType1(gc, 'FamilyName_t')
                if fname is not None:
                    fname = Internal.getValue(fname)
                    zBC = extractBCOfName(z, name, reorder, extrapFlow)
                    typeGC = 'FamilySpecified:%s'%fname
                    BCs.append(zBC); BCNames.append(name); BCTypes.append(typeGC)
    return (BCs, BCNames, BCTypes)

# merge GridConnectivity nodes (with the same receptor/donor zones)
# Warning: different GridConnectivityProperty are not considered yet
def _mergeGCs(z):
    dictOfGCs={}
    nodes = Internal.getNodesFromType2(z,'GridConnectivity_t')
    for gc in nodes:
        gcname = Internal.getName(gc)
        zoppname = Internal.getValue(gc)
        gctype = Internal.getNodeFromType(gc,'GridConnectivityType_t')
        gctype = Internal.getValue(gctype)
        if gctype == 'Abutting1to1':
                    #GCP = Internal.getNodeFromType(gc,'GridConnectivityProperty_t')
            PL = Internal.getNodeFromName(gc,'PointList')[1]
            PLD = Internal.getNodeFromName(gc,'PointListDonor')[1]
            if zoppname not in dictOfGCs:
                dictOfGCs[zoppname]=[PL,PLD,gcname]
            else:
                PLN = dictOfGCs[zoppname][0]
                PLDN = dictOfGCs[zoppname][1]
                dictOfGCs[zoppname][0]=numpy.concatenate([PLN,PL],axis=1)
                dictOfGCs[zoppname][1]=numpy.concatenate([PLDN,PLD],axis=1)
                Internal._rmNodesFromName(z,gcname)

    for zoppname in dictOfGCs:
        PL = dictOfGCs[zoppname][0]
        PLD = dictOfGCs[zoppname][1]
        gcname= dictOfGCs[zoppname][2]
        gcnode = Internal.getNodeFromName(z, gcname)
        PLNode = Internal.getNodeFromName(gcnode, 'PointList')
        PLDNode = Internal.getNodeFromName(gcnode, 'PointListDonor')
        PLNode[1] = PL
        PLDNode[1] = PLD
    return None

# Merge BCs on an unstructured zone
def _mergeBCs(z):
    import Transform.PyTree as T
    bcs = getBCs(z)
    BCs=bcs[0]; BCNames=bcs[1]; BCTypes=bcs[2]
    # merge par type
    alltypes = {}
    for c, t in enumerate(BCTypes):
        if t not in alltypes: alltypes[t] = BCs[c]
        else: alltypes[t] += BCs[c]
    # rebuild
    for i in alltypes: alltypes[i] = T.join(alltypes[i])
    BCs=[]; BCNames=[]; BCTypes=[]
    for i in alltypes:
        BCs.append(alltypes[i])
        BCNames.append('Merged'+i)
        BCTypes.append(i)

    _recoverBCs(z, (BCs,BCNames,BCTypes))
    return None

def isXZone(zone):
    """Check whether a zone has been added by addXZones."""
    r = Internal.getNodeFromName1(zone, 'XZone')
    if r is None: return False
    else: return True

# Extract fields on all match connectivities
def extractAllBCMatch(t, varList=None):
    zones    = Internal.getZones(t)
    allMatch = {}

    for z in zones:
        if not isXZone(z):
            dim = Internal.getZoneDim(z)
            if dim[0] == 'Structured':
                gcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
            else:
                gcs = Internal.getNodesFromType2(z, 'GridConnectivity_t')

            for gc in gcs:
                zname  = Internal.getValue(gc)
                zdonor = Internal.getNodeFromName2(t, zname)

                # Extraction BCMatch pour la zone donneuse
                [indR,fldD]  = extractBCMatch(zdonor,gc,dim,varList)
                key          = z[0]+"/"+gc[0]
                if fldD is not None: allMatch[key] = [indR,fldD]

    return allMatch

def computeBCMatchField(z, allMatch, variables=None):

    # Fields array in receiver zone (z)
    # =================================
    fields = []
    dim    = Internal.getZoneDim(z)

    # Type de la zone
    # ================
    if dim[0] == 'Structured': zoneType=1
    else:
        zoneType = 2; eltName = dim[3]
        if eltName=='NGON': pass
        else: raise ValueError("computeBCMatchField: not yet implement for basic elements.")

    # Liste des variables
    # ====================
    if variables is not None:
        if not isinstance(variables, list): varList = [variables]
        else: varList = variables
    else:
        varList=[]
        FS = Internal.getNodeFromName1(z,Internal.__FlowSolutionCenters__)
        for fs in FS[2]:
            if Internal.getType(fs) == 'DataArray_t':
                varList.append(Internal.getName(fs))

    # Traitement pour maillage struture
    # =================================
    if zoneType == 1: # Structured mesh
        # Tableau des champs a extraire
        for var in varList:
                # on verifie qu'on cherche des variables aux centres
            spl = var.split(':')
            if len(spl) != 1:
                if spl[0] != 'centers':
                    raise TypeError("computeBCMatchField: expected variables at centers location.")
            else: var = 'centers:'+var

            fld = getField(var,z)[0]
            if fld != []: fields.append(fld)

        if fields != []:
            fields = Internal.convertDataNodes2Array2(fields, dim, connects=[], loc=1)

        # Fields array in zdonor (stored in allMatch)
        # ===========================================
        dim = Internal.getZoneDim(z)
        ni = dim[1]-1; nj = dim[2]-1; nk = dim[3]-1

        fld = None; indR = None

        # ============================= Traitement TNC  ===================================
        isTNC    = False
        allCount = {} # dictionnaire nb occurence des faces

        if isTNC:
            # Concatenation de tous les indices de la zone
            # l'objectif est de detecter les indices presents plusieurs fois et
            # de compter leur nombre d'occurence (-> allCount)
            indRzone = []
            for key in allMatch:
                if key.split("/")[0] == z[0]:
                    [indR,fldD] = allMatch[key]
                    if indRzone != []:
                        indRzone = numpy.concatenate((indRzone,indR))
                    else:
                        indRzone = indR

            # Test l'existence de doublons
            if indRzone != []:

                [indUniq, indIndir, indCount] = numpy.unique(indRzone,return_inverse=True,return_counts=True)

                nind1 = indRzone.size
                nind2 = indUniq.size

                if nind2 != nind1: # il y a des indices presents plusieurs fois

                    shift = 0

                    for key in allMatch:
                        if key.split("/")[0] == z[0]:
                            [indR,fldD] = allMatch[key]

                            ncount = numpy.zeros(indR.size, dtype=Internal.E_NpyInt)

                            for i in range(indR.size):
                                indx      = indIndir[i+shift]
                                ncount[i] = indCount[indx]

                            allCount[key] = ncount

                            shift += indR.size
        # ============================= Fin traitement TNC  ===============================

        for key in allMatch:
            if key.split("/")[0] == z[0]:
                [indR1,fldD] = allMatch[key]

                if key in allCount.keys():
                    ncount = allCount[key]
                    # print(key, ncount)
                else:
                    ncount = None

                if fields != []:
                    fld1 = Converter.converter.buildBCMatchFieldStruct(fields,indR1,fldD,ncount)

                    if fld is not None:
                        fld.append(fld1)
                        indR = numpy.concatenate((indR,indR1))
                    else:
                        fld  = [fld1]
                        indR = indR1


    # Traitement pour maillage NGON
    # ==============================
    else: # NGON
        varL = []
        for var in varList:
            spl = var.split(':')
            if len(spl) != 1: varL.append(spl[1])
            else: varL.append(spl[0])

        fld  = []; indR = None

        for key in allMatch:
            if key.split("/")[0] == z[0]:
                [indR1,fldD] = allMatch[key]

                fld1 = Converter.converter.buildBCMatchFieldNG(z,indR1,fldD, varL,
                                                               Internal.__GridCoordinates__,
                                                               Internal.__FlowSolutionNodes__,
                                                               Internal.__FlowSolutionCenters__)

                if indR is not None:
                    indR = numpy.concatenate((indR,indR1))
                else:
                    indR = indR1

                fld.append(fld1)

    return indR, fld


# ===================================================================================
# Extraction des champs sur les raccords de type match
# Le champs en centre est extrapole sur les centres des faces
# IN
# ===
# zdonor  : zone donneuse
# gc      : GridConnectivity_t de la zone receveuse > indique le raccord  extraire
# dimzR   : dimensions de la zone receveuse (pour construire le tab d'indirection)
# varList : liste des variables a extraire (champs en centres)
# OUT
# ===
# indFaceD : indices des faces de la frontiere dans la zone donneuse
# indFaceR : indices des faces de la frontiere dans la zone receveuse
# fldFace  : champ de la zone donneuse extrapole sur les faces frontieres
# ===================================================================================
def extractBCMatch(zdonor,gc,dimzR,variables=None):
        # On verifie que gc donne le raccord dans zdonor
        # ==============================================
        # print("================================================")
        # print("zdonor :", zdonor[0])
        # print("gc : ", gc[0])
        # if Internal.getValue(gc) != zdonor[0]:
                # raise ValueError("extractBCMatch: GridConnectivity doesn't match zdonor.")

    dim = Internal.getZoneDim(zdonor)

    # Type de la zone
    # ================
    if dim[0]=='Structured': zoneType=1
    else:
        zoneType = 2; eltName = dim[3]
        if eltName=='NGON': pass
        else: raise ValueError("extractBCMatch: not yet implement for basic elements.")

    fields = []

    # Liste des variables
    # ====================
    if variables is not None:
        if not isinstance(variables, list): varList = [variables]
        else: varList = variables
    else:
        varList = []
        FS = Internal.getNodeFromName1(zdonor,Internal.__FlowSolutionCenters__)
        for fs in FS[2]:
            if Internal.getType(fs) == 'DataArray_t':
                varList.append(Internal.getName(fs))

    # Traitement pour maillage struture
    # =================================
    if zoneType == 1: # Structured mesh
        # Tableau des champs a extraire
        for var in varList:
                # on verifie qu'on cherche des variables aux centres
            spl = var.split(':')
            if len(spl) != 1:
                if spl[0] != 'centers':
                    raise TypeError("extractBCMatch: expected variables at centers location.")
            else:
                var = 'centers:'+var
            fld = getField(var,zdonor)[0]

            if fld != []:
                fields.append(fld)

        if fields != [] :
            if zoneType==1: connects = []
            else: connects = Internal.getElementNodes(zdonor)

            fields = Internal.convertDataNodes2Array2(fields, dim, connects, loc=1)

    # else:
    #   if zoneType == 1: # Structured mesh

    #     fields = getAllFields(zdonor, 'centers')[0]
        if fields != []:
            # raise ValueError("extractBCMatch: Variable(s) not found:", variables)

            # Infos raccord
            # =============
            prr   = Internal.getNodeFromName1(gc,'PointRange')
            prd   = Internal.getNodeFromName1(gc,'PointRangeDonor')
            tri   = Internal.getNodeFromName1(gc,'Transform')
            tri   = Internal.getValue(tri)

            wr    = Internal.range2Window(prr[1])
            wd    = Internal.range2Window(prd[1])

            iminR = wr[0] ; imaxR = wr[1]
            jminR = wr[2] ; jmaxR = wr[3]
            kminR = wr[4] ; kmaxR = wr[5]

            iminD = wd[0] ; imaxD = wd[1]
            jminD = wd[2] ; jmaxD = wd[3]
            kminD = wd[4] ; kmaxD = wd[5]

            sizeR = (imaxR-iminR+1)*(jmaxR-jminR+1)*(kmaxR-kminR+1)
            sizeD = (imaxD-iminD+1)*(jmaxD-jminD+1)*(kmaxD-kminD+1)

            if sizeR != sizeD:
                fldD = None
                indR = None
                print("Warning: extractBCMatch: Not a coincident match: ", gc[0])
                return [indR,fldD]

            niR   = dimzR[1]-1
            njR   = dimzR[2]-1
            nkR   = dimzR[3]-1

            t1    = tri[0]
            t2    = tri[1]

            if len(tri) == 3: t3 = tri[2]
            else: t3 = 0

            [indR,fldD]  = Converter.converter.extractBCMatchStruct(fields,(iminD,jminD,kminD,imaxD,jmaxD,kmaxD),
                                                                           (iminR,jminR,kminR,imaxR,jmaxR,kmaxR),
                                                                           (niR,njR,nkR),(t1,t2,t3))
        else:
            fldD = None
            indR = None
            print('Warning: extractBCMatch: field not found ', variables)

    # Traitement pour maillage NGON
    # ==============================
    else: # NGON

        varL = []
        for var in varList:
            spl = var.split(':')
            if len(spl) !=1:
                varL.append(spl[1])
            else:
                varL.append(spl[0])

        indR = Internal.getNodeFromName1(gc, 'PointList')
        indD = Internal.getNodeFromName1(gc, 'PointListDonor')

        indR = indR[1][0]
        indD = indD[1][0]

        PE = Internal.getNodeFromName2(zdonor, 'ParentElements')
        if PE is None: Internal._adaptNFace2PE(zdonor, remove=False)

        fldD  = Converter.converter.extractBCMatchNG(zdonor, indD, varL,
                                                     Internal.__GridCoordinates__,
                                                     Internal.__FlowSolutionNodes__,
                                                     Internal.__FlowSolutionCenters__)

    # print("len(indR): ", len(indR))
    # print("len(fldD): ", len(fldD[1][0]) )

    return [indR,fldD]

# ===================================================================================
# ** ATTENTION : Code non fonctionnel pour le moment. **
# Extraction des champs sur les raccords de type no-match
# Le champs en centre est extrapole sur les centres des faces et pondere (via l'algo
# superMesh qui permet de decouper les faces)
# ===================================================================================
def extractAllBCMatchTNC(t,variables=None):
    zones       = Internal.getZones(t)
    allMatchTNC = {}

    # Variables a extraire
    # ====================
    if variables is not None:
        if not isinstance(variables, list): varList = [variables]
        else: varList = variables
    else:
        varList = []
        FS = Internal.getNodeFromName1(zones[0],Internal.__FlowSolutionCenters__)
        for fs in FS[2]:
            if Internal.getType(fs) == 'DataArray_t':
                varList.append(Internal.getName(fs))

    # Parcours des zones et raccords
    # ==============================
    for zoneA in zones:
        if not isXZone(zoneA):
            indRzA = []
            fldDzA = []
            dim    = Internal.getZoneDim(zoneA)

            if dim[0] != 'Structured':
                print("extractAllBCMatchTNC: not ready for unstructured grid.")
                return {}

            gcs  = Internal.getNodesFromType2(zoneA, 'GridConnectivity_t')     # TNC
            gcs += Internal.getNodesFromType2(zoneA, 'GridConnectivity1to1_t') # near-match

            for gcA in gcs:

                if (Internal.getNodeFromName1(gcA, '.Solver#Property') == None ):
                    continue
                else:
                    print("TNC or near-field match found: ", gcA[0])

                zname  = Internal.getValue(gcA)
                zoneB  = Internal.getNodeFromName(t,zname) # A MODIFIER POUR TNC

                gcBs   = Internal.getNodesFromType2(zoneB, 'GridConnectivity_t')
                gcBs  += Internal.getNodesFromType2(zoneB, 'GridConnectivity1to1_t')

                for gcB in gcBs:
                    zname = Internal.getValue(gcB)
                    if zname == zoneA[0]: break

                key    = zoneA[0]+"/"+gcA[0]

                [indR,fldD] = computeBCMatchTNC(zoneA,zoneB,gcA, gcB, varList)

                if indR != []:
                    allMatchTNC[key] = [indR,fldD]

    return allMatchTNC


def extractMatchOfName(t,bndName,reorder=True,extrapFlow=True):
    """Return match surfaces of given name as zones (in analogy with
    C.extractBCOfType for BCs).
    """
    try: import Transform.PyTree as T
    except: raise ImportError("extractMatchOfName: requires Transform module.")

    tp = list()
    for z in Internal.getZones(t):
        zgc = Internal.getNodeFromType1(z, 'ZoneGridConnectivity_t')
        if zgc is not None:
            for match in zgc[2]:
                if Internal.isName(match, bndName):
                    getBC__(match, z, T, tp, reorder=reorder, extrapFlow=extrapFlow)

    return tp



# ===================================================================================
# ** ATTENTION : Code non fonctionnel pour le moment. **
# Interpolation des champs sur les raccords de type no-match
# ===================================================================================
def computeBCMatchTNC(zoneA,zoneB,gcA,gcB, varList):

    try: import Transform.PyTree as T
    except: raise ImportError("computeBCMatchTNC: requires Transform module.")
    try: import Intersector.PyTree as XOR
    except: raise ImportError("computeBCMatchTNC: requires Intersector module.")

    # Type de la zone
    # ================
    dim = Internal.getZoneDim(zoneB)

    if dim[0] == 'Structured':
        zoneType = 1
    else:
        zoneType = 2; eltName = dim[3]
        if eltName == 'NGON':
            raise ValueError("computeBCMatchTNC: not yet implement for Ngon elements.")
        else:
            raise ValueError("computeBCMatchTNC: not yet implement for basic elements.")

    prr   = Internal.getNodeFromName1(gcA,'PointRange')
    wrA   = Internal.range2Window(prr[1])

    iminA = wrA[0] ; imaxA = wrA[1]
    jminA = wrA[2] ; jmaxA = wrA[3]
    kminA = wrA[4] ; kmaxA = wrA[5]

    prr   = Internal.getNodeFromName1(gcB,'PointRange')
    wrB   = Internal.range2Window(prr[1])

    iminB = wrB[0] ; imaxB = wrB[1]
    jminB = wrB[2] ; jmaxB = wrB[3]
    kminB = wrB[4] ; kmaxB = wrB[5]

    surfA0 = extractMatchOfName(zoneA, gcA[0], extrapFlow=True)
    clipB0 = extractMatchOfName(zoneB, gcB[0], extrapFlow=True)

    surfA = convertArray2NGon(surfA0, recoverBC=0)
    clipB = convertArray2NGon(clipB0, recoverBC=0)

    XOR._convertNGON2DToNGON3D(surfA) # convert to ngon format expected by XOR (faces/nodes)
    XOR._convertNGON2DToNGON3D(clipB) # convert to ngon format expected by XOR (faces/nodes)

    hook  = createHook(zoneA, function='faceCenters')
    indR  = identifyFaces(hook, surfA) # indice des faces dans le maillage vol.
    indR  = indR-1 # shift

    (ancA, ancB, weight, isMatch) = XOR.superMesh2(surfA, clipB, tol=-1.e-4, proj_on_first=True)

    if isMatch: return [[],[]]

    fields = []

    for var in varList:
        # on verifie qu'on cherche des variables aux centres
        spl = var.split(':')
        if len(spl) != 1:
            if spl[0] != 'centers':
                raise TypeError("computeBCMatchTNC: expected variables at centers location.")
        else:
            var = 'centers:'+var

        fld = getField(var, clipB0)[0]

        if fld != []:
            fields.append(fld)

        if fields != []:
            if zoneType == 1: connects = []
            else: connects = Internal.getElementNodes(zoneB)

            fields = Internal.convertDataNodes2Array2(fields, dim, connects, loc=1)

        fx =XOR.extractBCMatchTNC(ancA, ancB, weight, fields, iminA, jminA, kminA,
                                  imaxA, jmaxA, kmaxA)

        return [indR, fx]



# ===================================================================================
# Extract fields at face centers defining a BC
# If varList is None -> variables defined at cell centers are extracted
# If a variable is defined in a BCDataSet, it is used, else Oth order extrapolation elsewhere
# returns a list of var names, corresponding fields defined by numpy arrays and the numpy array of bc indices
def extractBCFields(z, varList=None):
    """Extract fields on BCs."""
    typeZ = Internal.typeOfNode(z)
    zp = Internal.copyRef(z)

    if typeZ == 1: pass
    else:
        zp = Internal.getZones(zp)[0]
        print('Warning: valid for only one zone. Zone %s is selected.'%(zp[0]))


    if varList is None:
        varList=[]
        FS = Internal.getNodeFromName1(zp,Internal.__FlowSolutionCenters__)
        for fs in FS[2]:
            if Internal.getType(fs)=='DataArray_t':
                varList.append(Internal.getName(fs))

    dimZone = Internal.getZoneDim(zp)
    PE = None
    if dimZone[0]=='Structured': zoneType=1
    else:
        zoneType = 2; eltName = dimZone[3]
        if eltName == 'NGON': pass
        else: raise ValueError("extractBCFields: not yet implement for basic elements.")

    bcnodes = Internal.getNodesFromType2(zp,'BC_t')
    allFields=[];allIndices=[]; allVars=[]
    for bc in bcnodes:
        fieldsL=[]; varsL=[]; indicesL=[]

        #1. extract face indices
        if zoneType == 2:
            indicesL = Internal.getNodeFromName1(bc, 'PointList')
            indicesL = Internal.getValue(indicesL)[0]
        else:
            PR = Internal.getNodeFromName1(bc,'PointRange')
            win = Internal.range2Window(PR[1])
            imin = win[0]; imax = win[1]
            jmin = win[2]; jmax = win[3]
            kmin = win[4]; kmax = win[5]
            ni = dimZone[1]; nj = dimZone[2]; nk = dimZone[3]
            indicesL = Converter.converter.range2PointList(imin, imax, jmin, jmax, kmin, kmax, ni, nj, nk)

        #2. extract fields from BCDataSet nodes
        bcdata = Internal.getNodesFromType1(bc, 'BCDataSet_t')
        if bcdata is not None:
            alldatanodes = Internal.getBCDataSet(z, bc)
            for datanode in alldatanodes:
                dataname = Internal.getName(datanode)
                if dataname in varList:
                    fieldsL.append(Internal.getValue(datanode))
                    varsL.append(dataname)

        #3. no BCDataSet or variable not in BCDataSet -> Oth order extrapolation
        if len(varsL) < len(varList):
            varsE=[]
            for var in varList:
                if var not in varsL: varsE.append(var)


            if zoneType==2:
                if eltName =='NGON':
                    PE = Internal.getNodeFromName2(zp, 'ParentElements')
                    if PE is None: Internal._adaptNFace2PE(zp, remove=False)
                else:
                    raise TypeError('extractBCFields: basic elements not yet implemented.')

            locI = 1# volume fields located at centers
            fieldsL += Converter.converter.extractBCFields(zp, indicesL, varsE, locI,
                                                           Internal.__GridCoordinates__,
                                                           Internal.__FlowSolutionNodes__,
                                                           Internal.__FlowSolutionCenters__)
            varsL += varsE

        #3. add to list of variables/fields
        allFields.append(fieldsL)
        allIndices.append(indicesL)
        allVars.append(varsL)
    return allVars, allFields, allIndices

# -- extractBCDataStruct__
# Extract BC characteristics of structured zones [Name,BCType,Data,coords]
# Name is the name of the BC
# Type is the type of BC
# Data is the first BCDataSet node attached to some BCs, None if not used
# coords: coordinates of 9 particular points defining the zone
# BCMatch and BCNearMatch are not taken into account...
def extractBCDataStruct__(z):
    import Transform as T
    infos = []
    dims = Internal.getZoneDim(z)
    ni = dims[1]; nj = dims[2]; nk = dims[3]

    # BCOverlap
    nodes = Internal.getNodesFromType2(z, 'GridConnectivity_t')
    for i in nodes:
        bcname = i[0]
        bctype = Internal.getValue(i)
        bcdata = Internal.getNodesFromType1(i, 'BCDataSet_t')
        if bcdata != []: bcdata = bcdata[0]
        else: bcdata = None
        fname = Internal.getNodeFromType1(i, 'FamilyName_t')
        if fname is not None: fname = Internal.getValue(fname)

        r = Internal.getNodeFromName1(i, 'GridConnectivityType')
        if r is not None:
            val = Internal.getValue(r)
            if val == 'Overset':
                pr = Internal.getNodeFromName1(i, 'PointRange')
                if pr is not None:
                    if fname is not None: bctype = 'FamilySpecified:'+fname
                    range0 = pr[1]
                    w = Internal.range2Window(range0)
                    imin = w[0]; imax = w[1]; jmin = w[2]; jmax = w[3]; kmin = w[4]; kmax = w[5]
                    coords = getFields(Internal.__GridCoordinates__,z)[0]
                    coords = T.subzone(coords, (imin,jmin,kmin), (imax,jmax,kmax))
                    info = [bcname, 'BCOverlap', bcdata, coords]
                    infos.append(info)

    # All classical BC
    nodes = Internal.getNodesFromType2(z, 'BC_t')
    for i in nodes:
        bcname = i[0]
        bctype = Internal.getValue(i)
        bcdata = Internal.getNodesFromType1(i, 'BCDataSet_t')
        if bcdata != []: bcdata = bcdata[0]
        else: bcdata = None
        fname = Internal.getNodeFromType1(i, 'FamilyName_t')
        if fname is not None: fname = Internal.getValue(fname)
        pr = Internal.getNodeFromName1(i, 'PointRange')
        if pr is not None:
            if fname is not None: bctype = 'FamilySpecified:'+fname
            range0 = pr[1]
            w = Internal.range2Window(range0)
            imin = w[0]; imax = w[1]; jmin = w[2]; jmax = w[3]; kmin = w[4]; kmax = w[5]
            coords = getFields(Internal.__GridCoordinates__,z)[0]
            coords = T.subzone(coords, (imin,jmin,kmin), (imax,jmax,kmax))
            info = [bcname, bctype, bcdata, coords]
            infos.append(info)
    return infos

def extractBCInfo(t):
    infos = []
    for z in Internal.getZones(t):
        dims = Internal.getZoneDim(z)
        if dims[0] == 'Structured': infos += extractBCDataStruct__(z)
    return infos

#=============================================================================
# -- Fill empty BCs avec une CL donnee --
#=============================================================================
def fillEmptyBCWith(t, bndName, bndType, dim=3):
    """Fill empty BCs with given type."""
    a = Internal.copyRef(t)
    _fillEmptyBCWith(a, bndName, bndType, dim)
    return a

def _fillEmptyBCWith(t, bndName, bndType, dim=3):
    for z in Internal.getZones(t):
        c = 1
        wins = getEmptyBC(z, dim)
        for w in wins:
            ztype = Internal.getZoneType(z)
            if ztype == 1: # structured
                _addBC2Zone(z, bndName+str(c), bndType, w); c += 1
            else:
                dims = Internal.getZoneDim(z)
                if dims[0] == 'Unstructured':
                    eltType = dims[3]
                    if eltType != 'MIXED':
                        _addBC2Zone(z, bndName+str(c), bndType, faceList=w); c += 1
                    else:
                        try:
                            import Transform.PyTree as T
                            zbc = T.subzone(z, w, type='faces')
                            _addBC2Zone(z, bndName+str(c), bndType, subzone=zbc); c += 1
                        except:
                            raise ImportError("_fillEmptyBCWith: requires Transform module for unstructured MIXED zones")

    return None

#==============================================================================
# -- Family management --
#==============================================================================

# -- tagWithFamily (familyZone)
def tagWithFamily(z, familyName, add=False):
    """Tag a zone node or a BC node with a familyName.
    Usage: tagWithFamily(z, familyName, add)"""
    a = Internal.copyRef(z)
    _tagWithFamily(a, familyName, add)
    return a

def _tagWithFamily__(a, familyName, add=False):
    if a[3] != 'Zone_t' and a[3] != 'BC_t':
        print('Warning: tagWithFamily: must be used on a Zone_t or BC_t node.')
    if not add:
        Internal._createUniqueChild(a, 'FamilyName', 'FamilyName_t', value=familyName)
    else:
        if Internal.getNodeFromType1(a, 'FamilyName_t') is None:
            Internal._createChild(a, 'FamilyName', 'FamilyName_t', value=familyName)
        else:
            Internal._createChild(a, 'AdditionalFamilyName', 'AdditionalFamilyName_t', value=familyName)
    return None

def _tagWithFamily(a, familyName, add=False):
    zones = Internal.getZones(a)
    if len(zones) > 0: # family of zones
        for z in zones: _tagWithFamily__(z, familyName, add)
    else:
        nodes = Internal.getNodesFromType(a, 'BC_t')
        for n in nodes: _tagWithFamily__(n, familyName, add)

# -- getFamilyZones (wildcard possible on familyName)
def getFamilyZones(t, familyName):
    """Return all zones that have this familyName.
    Usage: getFamilyZones(t, familyName)"""
    out = []
    if isinstance(familyName, str): families = [familyName]
    else: families = familyName

    for z in Internal.getZones(t):
        res = Internal.getNodesFromType1(z, 'FamilyName_t')
        res += Internal.getNodesFromType1(z, 'AdditionalFamilyName_t')
        for i in res:
            val = Internal.getValue(i)
            for f in families:
                if ('*' in f)|('?' in f)|('[' in f):
                    if fnmatch.fnmatch(val, f): out.append(z)
                else:
                    if val == f: out.append(z)
    return out

def getFamilyBCZones(t, familyBCName):
    """Return all zones that have this familyBCName.
    Usage: getFamilyBCZones(t, familyBCName)"""
    out = []
    if isinstance(familyBCName, str): families = [familyBCName]
    else: families = familyBCName

    for z in Internal.getZones(t):
        for bc in Internal.getNodesFromType2(z,'BC_t'):
            res = Internal.getNodesFromType1(bc, 'FamilyName_t')
            for i in res:
                val = Internal.getValue(i)
                for f in families:
                    if ('*' in f)|('?' in f)|('[' in f):
                        if fnmatch.fnmatch(val, f): out.append(z)
                    else:
                        if val == f: out.append(z)

        for gc in Internal.getNodesFromType1(z,"ZoneGridConnectivity_t"):
            res = Internal.getNodesFromType2(gc, 'FamilyName_t')
            for i in res:
                val = Internal.getValue(i)
                for f in families:
                    if ('*' in f)|('?' in f)|('[' in f):
                        if fnmatch.fnmatch(val, f): out.append(z)
                    else:
                        if val == f: out.append(z)
    return out

# -- getFamilyBCs (wildcards possible on familyName)
def getFamilyBCs(t, familyName):
    """Return all BC nodes that have this familyName.
    Usage: getFamilyBCs(t, familyName)"""
    out = []
    if isinstance(familyName, str): families = [familyName]
    else: families = familyName
    for z in Internal.getZones(t):
        nodes = Internal.getNodesFromType2(z, 'BC_t')
        nodes += Internal.getNodesFromType2(z, 'GridConnectivity_t')
        nodes += Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
        for n in nodes:
            res = Internal.getNodesFromType1(n, 'FamilyName_t')
            for i in res:
                val = Internal.getValue(i)
                val = val.strip()
                for f in families:
                    if ('*' in f)|('?' in f)|('[' in f):
                        if fnmatch.fnmatch(val, f): out.append(n)
                    else:
                        if val == f: out.append(n)
    return out

# -- getFamilyBCNamesOfType (wildcards possible on bndType)
def getFamilyBCNamesOfType(t, bndType=None):
    """Return the family BC names of a given type.
    Usage: names = getFamilyBCNamesOfType(t, 'BCWall')"""
    out = set()
    nodes = Internal.getNodesFromType2(t, 'Family_t')
    if bndType is None:
        for n in nodes:
            p = Internal.getNodeFromType1(n, 'FamilyBC_t')
            if p is not None: out.add(n[0])
    else:
        for n in nodes:
            p = Internal.getNodeFromType1(n, 'FamilyBC_t')
            if p is not None:
                if Internal.isValue(p, bndType): out.add(n[0])
    return list(out)

# -- getFamilyBCNamesDict
def getFamilyBCNamesDict(t):
    """Return the dictionary of familyBCs. d['CARTER'] = "BCWall"
    Usage: dict = getFamilyBCNamesDict(t)"""
    d = {}
    nodes = Internal.getNodesFromType2(t, 'Family_t')
    for n in nodes:
        btype = Internal.getNodeFromType1(n, 'FamilyBC_t')
        if btype is not None:
            d[n[0]] = Internal.getValue(btype)
    return d

# -- getFamilyZoneNames
def getFamilyZoneNames(t):
    """Return the family zone names in a tree or a base."""
    out = set()
    nodes = Internal.getNodesFromType2(t, 'Family_t')
    for n in nodes:
        ret = Internal.getNodeFromType(n, 'FamilyBC_t')
        if ret is None: out.add(n[0])
    return list(out)

# -- Add a Family node to a base node
def addFamily2Base(base, familyName, bndType=None):
    """Add a family node to a base node.
    Usage: addFamily2Base(base, familyName, bndType)"""
    a = Internal.copyRef(base)
    _addFamily2Base(a, familyName, bndType)
    return a

def _addFamily2Base(a, familyName, bndType=None):
    bases = Internal.getBases(a)
    for b in bases:
        res = Internal.getNodeFromName1(b, familyName)
        child = Internal.createUniqueChild(b, familyName, 'Family_t', None)

        if bndType == 'BCOverlap': # Cas special
            Internal.createUniqueChild(child, 'FamilyBC', 'FamilyBC_t', 'UserDefined')
            Internal.createUniqueChild(child, '.Solver#Overlap', 'UserDefinedData_t', None)
        elif bndType == 'Mask': # Cas special
            Internal.createUniqueChild(child, 'FamilyBC', 'FamilyBC_t', bndType)
            mask = Internal.createUniqueChild(child, '.Solver#Mask', 'DataArray_t', None)
            Internal.createUniqueChild(mask, 'type', 'DataArray_t', 'xray')
            Internal.createUniqueChild(mask, 'dim1', 'DataArray_t', 1000)
            Internal.createUniqueChild(mask, 'dim2', 'DataArray_t', 1000)
            Internal.createUniqueChild(mask, 'delta', 'DataArray_t', 0.1)
        elif bndType is not None:
            Internal.createUniqueChild(child, 'FamilyBC', 'FamilyBC_t', bndType)
    return None

#==============================================================================
# -- change of localisation --
#==============================================================================

# -- node2Center
def node2Center(t, var='', accurate=0):
    """Convert zone defined at nodes to a zone defined at centers.
    Usage: node2Center(t, varname)"""
    a = Internal.copyRef(t)
    if var == '': # all grid
        if Internal.isStdNode(a) == 0: la = a
        else: la = [a] # noeud standard => liste de noeuds standards
        for i in range(len(la)):
            zones = Internal.getZones(la[i])
            for z in zones:
                (p, pos) = Internal.getParentOfNode(la[i], z)
                fieldc = getFields(Internal.__FlowSolutionCenters__, z)
                _deleteFlowSolutions__(z, 'centers')
                _deleteZoneBC__(z)
                _deleteGridConnectivity__(z)
                dim = Internal.getZoneDim(z)[0]
                if dim == 'Unstructured':
                    f = fieldc[0]
                    if f != []: f[3] = f[3].replace('*', '')
                    try: import Transform
                    except: pass
                    else:
                        z = TZA1(z, 'nodes', 'nodes', True, Transform.dual, 0)
                else:
                    z = TZA1(z, 'nodes', 'nodes', True, Converter.node2Center, accurate)
                setFields(fieldc, z, 'nodes', writeDim=False)
                if p is None: la[i] = z
                else: p[2][pos] = z
        if Internal.isStdNode(a) == 0: a = la
        else: a = la[0] # liste de noeuds standards => noeud standard
        return a
    elif var == Internal.__GridCoordinates__: # coordinates
        fieldc = []
        fieldn = getFields(Internal.__GridCoordinates__, a, api=3)
        for i in fieldn:
            if i != []:
                b = Converter.node2Center(i, accurate)
                fieldc.append(b)
            else:
                fieldc.append([])
        setFields(fieldc, a, 'centers')
        return a
    elif var == Internal.__FlowSolutionNodes__: # all node fields
        fieldc = []
        fieldn = getFields(Internal.__FlowSolutionNodes__, a, api=3)
        for i in fieldn:
            if i != []:
                b = Converter.node2Center(i, accurate)
                fieldc.append(b)
            else:
                fieldc.append([])
        setFields(fieldc, a, 'centers', writeDim=False)
        return a
    else: # var lists
        if isinstance(var, list):
            for v in var: _node2Center__(a, v, accurate)
        else: _node2Center__(a, var, accurate)
        return a

# commentaire SP : mauvais fonctionnement : modifie a et retourne a
# est appele par etc tel quel...
def node2Center__(a, var, accurate=0):
    var, loc = Internal.fixVarName(var)
    if loc == 1: return a
    fieldn = getField(var, a, api=3)
    fieldc = Converter.node2Center(fieldn, accurate)
    setFields(fieldc, a, 'centers', writeDim=False)
    return a

def _node2Center__(a, var, accurate=0):
    var, loc = Internal.fixVarName(var)
    if loc == 1: return None
    fieldn = getField(var, a, api=3)
    fieldc = Converter.node2Center(fieldn, accurate)
    setFields(fieldc, a, 'centers', writeDim=False)
    return None

# Adapte les arrays pour les cas nk=1 => passe en nk=2 si la zone est nk=2
def _patchArrayForCenter2NodeNK1__(fields, a):
    zones = Internal.getZones(a)
    c = 0
    for f in fields:
        if f != []:
            z = zones[c]
            dim = Internal.getZoneDim(z)
            ni = dim[1]; nj = dim[2]; nk = dim[3]
            fp = f[1]
            if isinstance(fp, list): # array2/3
                nfld = len(fp); s = fp[0].size
            else: # array1
                nfld = fp.shape[0]; s = fp.shape[1]

            if dim[0] == 'Structured' and nk == 2 and s == ni*nj:
                if isinstance(fp, list): # array2/3
                    for n, fl in enumerate(fp):
                        b = numpy.empty( (ni*nj*2), dtype=numpy.float64 )
                        b[0:ni*nj] = fl[0:ni*nj]
                        b[ni*nj:2*ni*nj] = fl[0:ni*nj]
                        fp[n] = b
                else:
                    b = numpy.empty( (nfld, ni*nj*2), dtype=numpy.float64 )
                    b[:,0:ni*nj] = fp[:,0:ni*nj]
                    b[:,ni*nj:2*ni*nj] = fp[:,0:ni*nj]
                    f[1] = b; f[4] = 2
        c += 1

# -- center2Node
# Convert a zone defining centers to nodes or convert a field in a
# base/tree/zone located at centers to nodes
# if useGhost: first addGhostCells before center2Node
def center2Node(t, var=None, cellNType=0, useGhost=True):
    """Converts array defined at centers to an array defined at nodes.
    Usage: center2Node(t, var, cellNType)"""

# Preparation pour le topTree
#     istoptree = Internal.isTopTree(t)
#     if not istoptree and topTree != []:
#       zones = C.getConnectedZones(t, topTree=topTree)
#       tp = C.newPyTree(['Base'])
#       zonest = Internal.getZones(t)
#       tp[2][1][2] += zonest; tp[2][1][2] += zones
#       tp = C.center2Node(tp, Internal.__FlowSolutionCenters__)
#       zone = tp[2][1][2][0]

    ghost = Internal.getNodeFromName(t, 'ZoneRind')
    if var is None: # all grid
            # solution en centres
        res = Internal.getNodesFromName3(t, Internal.__FlowSolutionCenters__)
        fieldsc = []
        if res != []:
            if Internal.getNodesFromType1(res[0], 'DataArray_t') != []:
                if ghost is None and useGhost:
                    a = Internal.addGhostCells(t, t, 1, adaptBCs=0, modified=[Internal.__FlowSolutionCenters__])
                else: a = Internal.copyRef(t)
                fieldc = getFields(Internal.__FlowSolutionCenters__, a)
                fieldn = []
                listVar = []
                for i in fieldc:
                    if i != []:
                        b = Converter.center2Node(i, cellNType); fieldn.append(b)
                        for va in b[0].split(','):
                            if va not in listVar: listVar.append(va)
                    else: fieldn.append([])
                _patchArrayForCenter2NodeNK1__(fieldn, a)
                setFields(fieldn, a, 'nodes', writeDim=False)
                if ghost is None and useGhost:
                    a = Internal.rmGhostCells(a, a, 1, adaptBCs=0,
                                              modified=[listVar, Internal.__FlowSolutionCenters__])
                fieldsc = getFields(Internal.__FlowSolutionNodes__, a)

        # destruction
        t = deleteFlowSolutions__(t, 'centers')
        t = TZA1(t, 'nodes', 'nodes', True, Converter.center2Node, cellNType, None)
        if fieldsc != []: setFields(fieldsc, t, 'centers', writeDim=False)
        return t

    elif var == Internal.__FlowSolutionCenters__: # all center fields
        res = Internal.getNodesFromName3(t, Internal.__FlowSolutionCenters__)
        if res == []: return t
        if Internal.getNodesFromType1(res[0], 'DataArray_t') == []: return t
        if ghost is None and useGhost:
            a = Internal.addGhostCells(t, t, 1, adaptBCs=0,
                                       modified=[Internal.__FlowSolutionCenters__])
        else: a = Internal.copyRef(t)
        fieldc = getFields(Internal.__FlowSolutionCenters__, a, api=3)
        fieldn = []
        listVar = []
        for i in fieldc:
            if i != []:
                b = Converter.center2Node(i, cellNType); fieldn.append(b)
                for va in b[0].split(','):
                    if va not in listVar: listVar.append(va)
            else: fieldn.append([])
        _patchArrayForCenter2NodeNK1__(fieldn, a)
        setFields(fieldn, a, 'nodes', writeDim=False)
        if ghost is None and useGhost:
            a = Internal.rmGhostCells(a, a, 1, adaptBCs=0,
                                      modified=[listVar,Internal.__FlowSolutionCenters__])
        return a

    else: # var list
        if isinstance(var, list): vars = var
        else: vars = [var]
        if ghost is None and useGhost:
            a = Internal.addGhostCells(t, t, 1, adaptBCs=0, modified=vars)
        else: a = Internal.copyRef(t)
        for v in vars: _center2Node__(a, v, cellNType)
        # if there are center vars in the list, add equivalent node vars because
        # they have been created by center2Node
        ghost2 = Internal.getNodeFromName(a, 'ZoneRind')
        if ghost is None and ghost2 is not None and useGhost:
            var2 = vars[:]
            for v in vars:
                variable = v.split(':')
                if len(variable) == 2 and variable[0] == 'centers':
                    var2.append(variable[1])
            a = Internal.rmGhostCells(a, a, 1, adaptBCs=0, modified=var2)
        return a

# mauvais fonctionnement : modifie a et retourne a
# mais appele dans etc tel quel.
def center2Node__(a, var, cellNType):
    _center2Node__(a, var, cellNType)
    return a

def _center2Node__(a, var, cellNType):
    fieldc = getField(var, a, api=3)
    fieldn = []
    for i in fieldc:
        if i != []: i = Converter.center2Node(i, cellNType)
        fieldn.append(i)
    _patchArrayForCenter2NodeNK1__(fieldn, a)
    setFields(fieldn, a, 'nodes', writeDim=False)
    return None

# -- node2ExtCenters
# Convert a zone to an extended center zone
# If FlowSolutionCenters exist, they are also converted to extended centers
def node2ExtCenter(t, var=''):
    """Convert zones defined on nodes to zones defined on extended centers.
    Usage: node2ExtCenter(a)"""
    a = Internal.copyRef(t)
    if var == '': # on prend ts les champs
        if Internal.isStdNode(a) == 0: la = a
        else: la = [a] # noeud standard => liste de noeuds standards
        for i in range(len(la)):
            zones = Internal.getZones(la[i])
            for z in zones:
                (p, pos) = Internal.getParentOfNode(la[i], z)
                fieldn = getAllFields(z, loc='nodes')
                fieldc = getFields(Internal.__FlowSolutionCenters__, z)
                _deleteFlowSolutions__(z, 'centers')
                _deleteZoneBC__(z)
                _deleteGridConnectivity__(z)
                dim = Internal.getZoneDim(z)[0]
                if dim == 'Structured':
                    fieldn2 = Converter.node2ExtCenter(fieldn)
                    if fieldc != [[]]:
                        fieldc2 = Converter.center2ExtCenter(fieldc)
                        fieldn2 = Converter.addVars([fieldn2,fieldc2])
                    setFields(fieldn2, z, 'nodes', True)
                    if p is None: la[i] = z
                    else: p[2][pos] = z
        if Internal.isStdNode(a) == 0: a = la
        else: a = la[0] # liste de noeuds standards => noeud standard
        return a

    elif var == Internal.__GridCoordinates__:
        fieldn = getFields(Internal.__GridCoordinates__, a)
        fielde = Converter.node2ExtCenter(fieldn)
        setFields(fielde, a, 'nodes', writeDim=True)
        return a

    elif var == Internal.__FlowSolutionNodes__:
        fieldn = getFields(Internal.__FlowSolutionNodes__, a)
        fielde = Converter.node2ExtCenter(fieldn)
        setFields(fielde, a, 'nodes', writeDim=True)
        return a
    else:
        raise ValueError("node2ExtCenter: only for all fields, coordinates or flow solution located at nodes.")

#==============================================================================
# diff 2 pyTrees
#==============================================================================
def diffArrays(A, B, removeCoordinates=True):
    """Compute the difference of two pyTrees, topologically."""
    t1 = Internal.copyRef(A); t2 = Internal.copyRef(B)
    zones1 = Internal.getZones(t1)
    zones2 = Internal.getZones(t2)
    nz = len(zones1)
    if nz != len(zones2):
        raise ValueError("diffArrays: different number of zones (A=%d; B=%d)."%(nz,len(zones2)))
    for no in range(nz):
        # noeuds
        A1 = getAllFields(zones1[no], 'nodes', api=3); A1 = Internal.clearList(A1)
        A2 = getAllFields(zones2[no], 'nodes', api=3); A2 = Internal.clearList(A2)
        # elimination des solutions aux noeuds
        node = Internal.getNodesFromName1(zones1[no], Internal.__FlowSolutionNodes__)
        if node != []:
            (parent, d) = Internal.getParentOfNode(t1, node[0])
            if parent is not None: del parent[2][d]

        if A1 != [] and A2 != []:
            diff = Converter.diffArrays(A1, A2)
            setFields(diff, zones1[no], 'nodes')

        # centres
        A1 = getAllFields(zones1[no], 'centers', api=3); A1 = Internal.clearList(A1)
        A2 = getAllFields(zones2[no], 'centers', api=3); A2 = Internal.clearList(A2)
        node = Internal.getNodesFromName1(zones1[no], Internal.__FlowSolutionCenters__)
        if node != []:
            (parent, d) = Internal.getParentOfNode(t1, node[0])
            if parent is not None: del parent[2][d]

        if A1 != [] and A2 != []:
            diff = Converter.diffArrays(A1, A2)
            setFields(diff, zones1[no], 'centers')
    if removeCoordinates: t1 = rmNodes(t1, Internal.__GridCoordinates__)
    return t1

def diffArrayGeom(A, B, tol=1.e-10, removeCoordinates=True):
    """Compute the difference of two pyTrees, geometrically."""
    t1 = Internal.copyRef(A); t2 = Internal.copyRef(B)
    zones1 = Internal.getZones(t1)
    zones2 = Internal.getZones(t2)
    nz = len(zones1)
    if nz != len(zones2):
        raise ValueError("diffArrays: different number of zones"
                         "(A=%d; B=%d)."%(nz,len(zones2)))

    for no in range(nz):
        # geometric diff
        diffn, diffc = diffArrayGeom__(zones1[no], zones2[no], tol)
        if diffn is None: return None # one array is different on coordinates
        # remplacement des solutions aux noeuds par diffn
        A1 = getAllFields(zones1[no], 'nodes', api=3); A1 = Internal.clearList(A1)
        A2 = getAllFields(zones2[no], 'nodes', api=3); A2 = Internal.clearList(A2)
        node = Internal.getNodesFromName1(zones1[no], Internal.__FlowSolutionNodes__)
        if node != []:
            (parent, d) = Internal.getParentOfNode(t1, node[0])
            if parent is not None: del parent[2][d]
        if A1 != [] and A2 != []:
            setFields(diffn, zones1[no], 'nodes')

        # remplacement des solutions aux centres par diffc
        if diffc is None: continue # one array is different on coordinates
        A1 = getAllFields(zones1[no], 'centers', api=3); A1 = Internal.clearList(A1)
        A2 = getAllFields(zones2[no], 'centers', api=3); A2 = Internal.clearList(A2)
        node = Internal.getNodesFromName1(zones1[no], Internal.__FlowSolutionCenters__)
        if node != []:
            (parent, d) = Internal.getParentOfNode(t1, node[0])
            if parent is not None: del parent[2][d]
        if A1 != [] and A2 != []:
            setFields(diffc, zones1[no], 'centers')

    if removeCoordinates: t1 = rmNodes(t1, Internal.__GridCoordinates__)
    return t1

def diffArrayGeom__(z1, z2, tol):
    # compare nodes
    a1 = getFields(Internal.__GridCoordinates__, z1, api=3)[0]
    f1 = getFields(Internal.__FlowSolutionNodes__, z1, api=3)[0]
    if f1 != []: a1 = Converter.addVars2([a1,f1])
    a2 = getFields(Internal.__GridCoordinates__, z2, api=3)[0]
    f2 = getFields(Internal.__FlowSolutionNodes__, z2, api=3)[0]
    if f2 != []: a2 = Converter.addVars2([a2,f2])
    diffn = Converter.diffArrayGeom__(a1, a2)
    if diffn is None: return None, None
    # compare centers
    z1c = node2Center(z1)
    a1 = getFields(Internal.__GridCoordinates__, z1c, api=3)[0]
    f1 = getFields(Internal.__FlowSolutionNodes__, z1c, api=3)[0]
    if f1 != []: a1 = Converter.addVars2([a1,f1])
    z2c = node2Center(z2)
    a2 = getFields(Internal.__GridCoordinates__, z2c, api=3)[0]
    f2 = getFields(Internal.__FlowSolutionNodes__, z2c, api=3)[0]
    if f2 != []: a2 = Converter.addVars2([a2,f2])
    diffc = Converter.diffArrayGeom__(a1, a2)
    return diffn, diffc

# Check if all fields are finite (no NAN no INF)
def isFinite(a, var=None):
    """Return true if all fields in a have no NAN or INF values."""
    if var is not None:
        var = var.replace('centers:', '')
        var = var.replace('nodes:', '')
    zones = Internal.getZones(a)
    containers = [Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__, Internal.__FlowSolutionCenters__]
    ret = True
    for z in zones:
        for name in containers:
            c = Internal.getNodeFromName1(z, name)
            if c is not None:
                nodes = Internal.getNodesFromType1(c, 'DataArray_t')
                for n in nodes:
                    if var is None or n[0] == var:
                        array = n[1]
                        #b = numpy.isfinite(array)
                        #res = numpy.all(b)
                        array = array.ravel(order="K")
                        res = Converter.converter.isFinite(array)
                        if res > 0:
                            ret = False
                            print('Warning: NAN or INF value in %s (%s)'%(n[0],z[0]))
    return ret

def setNANValuesAt(a, var=None, value=0.):
    """Set value if field is NAN."""
    ap = Internal.copyTree(a)
    _setNANValuesAt(ap, var=var, value=value)
    return ap

def _setNANValuesAt(a, var=None, value=0.):
    """Set value if field is NAN."""
    zones = Internal.getZones(a)
    containers = [Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__, Internal.__FlowSolutionCenters__]
    for z in zones:
        for name in containers:
            c = Internal.getNodeFromName1(z, name)
            if c is not None:
                nodes = Internal.getNodesFromType1(c, 'DataArray_t')
                for n in nodes:
                    if var is None or n[0] == var:
                        array = n[1]
                        array = array.ravel(order="K")
                        Converter.converter.setNANValuesAt(array, value)
    return None

#==============================================================================
# - add specific nodes -
#==============================================================================
# -- addState
# Add a single state/value or a full reference state
def addState(t, state=None, value=None, adim='adim1',
             MInf=None, alphaZ=0., alphaY=0., ReInf=1.e8,
             UInf=None, TInf=None, PInf=None, RoInf=None, LInf=None,
             Mus=None, MutSMuInf=0.2, TurbLevelInf=1.e-4,
             EquationDimension=None, GoverningEquations=None, Mtip=None,Mu_multiplier=1):
    """Add single state value or a full reference state."""
    tp = Internal.copyRef(t)
    _addState(tp, state, value, adim,
              MInf, alphaZ, alphaY, ReInf, UInf, TInf, PInf, RoInf, LInf,
              Mus, MutSMuInf, TurbLevelInf, EquationDimension,
              GoverningEquations, Mtip,Mu_multiplier)
    return tp

def _addState(t, state=None, value=None, adim='adim1',
              MInf=0.5, alphaZ=0., alphaY=0., ReInf=1.e8,
              UInf=None, TInf=None, PInf=None, RoInf=None, LInf=None,
              Mus=None, MutSMuInf=0.2, TurbLevelInf=1.e-4,
              EquationDimension=None, GoverningEquations=None, Mtip=None, Mu_multiplier=1):
    """Add single state value or a full reference state."""
    ntype = Internal.typeOfNode(t)
    if state is not None and value is not None: # single state value
        addState2Node2__(t, ntype, state, value); return None

    import KCore.Adim
    if state is None: # compute state
        if adim == 'adim1':
            if MInf is None: raise ValueError("addState: MInf is missing.")
            state = KCore.Adim.adim1(MInf, alphaZ, alphaY, ReInf,
                                     MutSMuInf, TurbLevelInf, Mtip)
        elif adim == 'adim2' or adim == 'adim2funk':
            state = KCore.Adim.adim2(MInf, alphaZ, alphaY, ReInf,
                                     MutSMuInf, TurbLevelInf)
        elif adim == 'adim3':
            state = KCore.Adim.adim3(MInf, alphaZ, alphaY, ReInf, LInf,
                                     MutSMuInf, TurbLevelInf, Mtip)
        elif adim == 'dim1':
            if UInf is None: raise ValueError("addState: UInf is missing.")
            if TInf is None: raise ValueError("addState: TInf is missing.")
            if PInf is None: raise ValueError("addState: PInf is missing.")
            if LInf is None: raise ValueError("addState: LInf is missing.")
            state = KCore.Adim.dim1(UInf, TInf, PInf, LInf, alphaZ, alphaY,
                                    MutSMuInf, TurbLevelInf, Mtip, Mu_multiplier)
        elif adim == 'dim2':
            if UInf is None: raise ValueError("addState: UInf is missing.")
            if TInf is None: raise ValueError("addState: TInf is missing.")
            if RoInf is None: raise ValueError("addState: RoInf is missing.")
            if LInf is None: raise ValueError("addState: LInf is missing.")
            state = KCore.Adim.dim2(UInf, TInf, RoInf, LInf, alphaZ, alphaY,
                                    MutSMuInf, TurbLevelInf, Mtip, Mu_multiplier)
        elif adim == 'dim3':
            if UInf is None: raise ValueError("addState: UInf is missing.")
            if PInf is None: raise ValueError("addState: PInf is missing.")
            if RoInf is None: raise ValueError("addState: RoInf is missing.")
            if LInf is None: raise ValueError("addState: LInf is missing.")
            state = KCore.Adim.dim3(UInf, PInf, RoInf, LInf, alphaZ, alphaY,
                                    MutSMuInf, TurbLevelInf, Mtip, Mu_multiplier)
        elif adim == 'dim4':
            if UInf is None: raise ValueError("addState: UInf is missing.")
            if TInf is None: raise ValueError("addState: TInf is missing.")
            if PInf is None: raise ValueError("addState: PInf is missing.")
            if LInf is None: raise ValueError("addState: LInf is missing.")
            if Mus is None: raise ValueError("addState: Mus is missing.")
            state = KCore.Adim.dim4(UInf, TInf, PInf, LInf, alphaZ, alphaY,
                                    Mus, MutSMuInf, TurbLevelInf, Mtip)

    UInf   = state[1] / state[0]
    VInf   = state[2] / state[0]
    WInf   = state[3] / state[0]
    addState2Node2__(t, ntype, 'VelocityX', UInf)
    addState2Node2__(t, ntype, 'VelocityY', VInf)
    addState2Node2__(t, ntype, 'VelocityZ', WInf)
    addState2Node2__(t, ntype, 'Density', state[0])
    addState2Node2__(t, ntype, 'MomentumX', state[1])
    addState2Node2__(t, ntype, 'MomentumY', state[2])
    addState2Node2__(t, ntype, 'MomentumZ', state[3])
    addState2Node2__(t, ntype, 'EnergyStagnationDensity', state[4])
    addState2Node2__(t, ntype, 'Pressure', state[5])
    addState2Node2__(t, ntype, 'Temperature', state[6])
    addState2Node2__(t, ntype, 'Cv', state[7])
    addState2Node2__(t, ntype, 'Mach', state[8])
    addState2Node2__(t, ntype, 'Reynolds', state[9])
    addState2Node2__(t, ntype, 'Gamma', state[11])
    addState2Node2__(t, ntype, 'Rok', state[12])
    addState2Node2__(t, ntype, 'RoOmega', state[13])
    addState2Node2__(t, ntype, 'TurbulentSANuTildeDensity', state[14])
    addState2Node2__(t, ntype, 'Mus', state[15])
    addState2Node2__(t, ntype, 'Cs', state[16])
    addState2Node2__(t, ntype, 'Ts', state[17])
    addState2Node2__(t, ntype, 'Pr', state[18])

    if EquationDimension is not None:
        addState2Node2__(t, ntype, 'EquationDimension', EquationDimension)
    if GoverningEquations is not None:
        addState2Node2__(t, ntype, 'GoverningEquations', GoverningEquations)

    return None

# Ajoute un noeud state/value suivant le type de t (in place)
def addState2Node2__(t, ntype, state, value):
    if ntype == 1: # add to zone
        addState2Node__(t, state, value)
    elif ntype == 2: # add to all zones
        for i in t: addState2Node__(i, state, value)
    elif ntype == 3: # add to all bases
        bases = Internal.getBases(t)
        for b in bases: addState2Node__(b, state, value)
    elif ntype == 4: # add to base
        addState2Node__(t, state, value)
    elif ntype == 5: # add to all bases
        for b in t: addState2Node__(b, state, value)
    else: addState2Node__(t, state, value) # direct to node

# Ajoute un noeud state/value au noeud a (in place)
def addState2Node__(a, state, value):
    # Container: FlowEquationSet or ReferenceState
    if state == 'EquationDimension' or state == 'GoverningEquations' or state == 'TurbulenceModel':
        H = []
        for n in a[2]:
            if n[0] == 'FlowEquationSet': H = n; break

        if H == []:
            a[2].append(['FlowEquationSet', None, [], 'FlowEquationSet_t'])
            H = a[2][len(a[2])-1]
    else:
        H = []
        for n in a[2]:
            if n[0] == 'ReferenceState': H = n; break

        if H == []:
            a[2].append(['ReferenceState', None, [], 'ReferenceState_t'])
            H = a[2][len(a[2])-1]

    # EquationDimension
    if state == 'EquationDimension':
        nodes = Internal.getNodesFromName(H, state)
        if nodes == []:
            v = numpy.empty((1), dtype=Internal.E_NpyInt); v[0] = value
            H[2].append([state, v, [], '"int"']) # Better DataArray_t
        else:
            v = numpy.empty((1), dtype=Internal.E_NpyInt); v[0] = value
            nodes[0][1] = v

    # GoverningEquations
    elif state == 'GoverningEquations':
        Internal._createUniqueChild(H, state, 'GoverningEquations_t', value=value)

    # TurbulenceModel
    elif state == 'TurbulenceModel':
        Internal._createUniqueChild(H, state, 'TurbulenceModel_t', value=value)

    # Reference state
    else:
        Internal._createUniqueChild(H, state, 'DataArray_t', value=value)

    return a

# -- getState
# Retourne un vecteur identique a Adim.
# Doit etre l'exact inverse de addState.
def getState__(state, name):
    A = Internal.getNodeFromName1(state, name)
    if A is None: raise ValueError("getState: %s is missing in tree ReferenceState."%name)
    return Internal.getValue(A)

def getState(t, name=None):
    state = Internal.getNodeFromName(t, 'ReferenceState')
    if state is None: raise ValueError("getState: ReferenceState is missing in tree.")
    if name is not None:
        value = getState__(state, name)
        return value
    RoInf = getState__(state, 'Density')
    RouInf = getState__(state, 'MomentumX')
    RovInf = getState__(state, 'MomentumY')
    RowInf = getState__(state, 'MomentumZ')
    RoeInf = getState__(state, 'EnergyStagnationDensity')
    PInf = getState__(state, 'Pressure')
    TInf = getState__(state, 'Temperature')
    cvInf = getState__(state, 'Cv')
    MInf = getState__(state, 'Mach')
    ReInf = getState__(state, 'Reynolds')
    Gamma = getState__(state, 'Gamma')
    RokInf = getState__(state, 'Rok')
    RoomegaInf = getState__(state, 'RoOmega')
    RonutildeInf = getState__(state, 'TurbulentSANuTildeDensity')
    Mus = getState__(state, 'Mus')
    Cs = getState__(state, 'Cs')
    Ts = getState__(state, 'Ts')
    Pr = getState__(state, 'Pr')
    return [RoInf, RouInf, RovInf, RowInf, RoeInf, PInf, TInf, cvInf, MInf,
            ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
            Mus, Cs, Ts, Pr]

# -- addChimera2Base
# add a chimera user defined node to a base
# this node stores specific Cassiopee Chimera settings
def addChimera2Base(base, setting, value):
    """Add chimera setting as node in base."""
    basep = Internal.copyRef(base)
    _addChimera2Base(basep, setting, value)
    return basep

def _addChimera2Base(base, setting, value):
    node = Internal.getNodeFromName1(base, '.Solver#Chimera')
    if node is not None:
        chimera = node
    else:
        chimera = ['.Solver#Chimera', None, [], 'UserDefinedData_t']
        base[2].append(chimera)

    # Priority
    if setting == 'Priority':
        v = numpy.empty((1,1), dtype=Internal.E_NpyInt); v[0,0] = value
        node = Internal.getNodeFromName1(chimera, 'Priority')
        if node is not None:
            node[1] = v
        else:
            a = ['Priority', v, [], 'DataArray_t']
            chimera[2].append(a)

    # XRayTol
    elif setting == 'XRayTol':
        v = numpy.empty((1,1), numpy.float64); v[0,0] = value
        node = Internal.getNodeFromName1(chimera, 'XRayTol')
        if node is not None:
            node[1] = v
        else:
            a = ['XRayTol', v, [], 'DataArray_t']
            chimera[2].append(a)

    # XRayDelta
    elif setting == 'XRayDelta':
        v = numpy.empty((1,1), numpy.float64); v[0,0] = value
        node = Internal.getNodeFromName1(chimera, 'XRayDelta')
        if node is not None:
            node[1] = v
        else:
            a = ['XRayDelta', v, [], 'DataArray_t']
            chimera[2].append(a)

    # DoubleWallTol
    elif setting == 'DoubleWallTol':
        v = numpy.empty((1,1), numpy.float64); v[0,0] = value
        node = Internal.getNodeFromName1(chimera, 'DoubleWallTol')
        if node is not None:
            node[1] = v
        else:
            a = ['DoubleWallTol', v, [], 'DataArray_t']
            chimera[2].append(a)

    # relations
    elif setting == '+' or setting == '-' or setting == '0' or setting == 'N':
        Internal._createChild(chimera, setting, 'UserDefinedData_t', value=value)

    return None

#==============================================================================
# Retourne un arbre ou certaines bases ont ete eliminees
#==============================================================================
def getPrimaryTree(t):
    import re
    tp = newPyTree()
    excluded = ['CONTOURS', 'BODY#', 'SLICES', 'SURFACES']
    bases = t[1:]
    for b in bases:
        found = False
        if b[3] == 'CGNSBase_t':
            for i in excluded:
                exp = re.compile(i)
                if exp.search(b[0]) is not None:
                    found = True
                    break
        if not found: tp.append(b)
    return tp

#==============================================================================
# check if name is already in list
# return True or False
#==============================================================================
def checkNameInList(name, list):
    for i in list:
        if name == i: return True
    return False

#==============================================================================
# mergeTrees
# Merge 2 arbres. Les bases de t2 sont ajoutees a t1. Les noms des bases
# sont modifies si elles existent deja dans t1.
#==============================================================================
def mergeTrees(t1, t2):
    t1p = Internal.copyRef(t1)
    # Enregistre les noms des bases de t1
    t1BaseNames = []
    bases = Internal.getBases(t1)
    for b in bases: t1BaseNames.append(b[0])

    # traitement par base
    bases = Internal.getBases(t2)
    for b in bases:
        ret = checkNameInList(b[0], t1BaseNames)
        if not ret: t1p[2].append(b)
        else:
            c = 0
            while ret:
                ret = checkNameInList('%s.%d'%(b[0],c), t1BaseNames)
                c += 1
            b[0] = '%s.%d'%(b[0],c-1)
            t1p[2].append(b)

    # noeud extra base
    nodes = []
    for n in t2[2]:
        if n[3] != 'CGNSBase_t' and n[3] != 'CGNSLibraryVersion_t': nodes.append(n)

    for n in nodes: print(n[0])
    t1p[2] += nodes

    return t1p

#==============================================================================
# Retourne les faces des BCs pour les zones non structurees et pour
# les BCs definies par faces (BC physiques) une liste par zone
# [ ['nomBC', numpy(faces), 'nomBC', numpy(faces)...] ...]
# si nameType=0: nomBC est le nom de la BC
# si nameType=1: nomBC concatene le nom de la BC et son type (BCName_BCType
# ou BCName_FamilyName)
#==============================================================================
def getBCFaces(t, nameType=0):
    BCFaces = []
    zones = Internal.getZones(t)
    for z in zones:
        BCFZ = []
        BCs = Internal.getNodesFromType2(z, 'BC_t')
        for b in BCs:
            name = b[0]
            BCtype = Internal.getValue(b)
            p = Internal.getNodeFromType1(b, 'GridLocation_t')
            n = Internal.getNodeFromName1(b, 'PointList')
            if p is not None:
                loc = Internal.getValue(p)
                if n is not None and (loc == 'FaceCenter' or loc == 'CellCenter'):
                    if BCtype == 'FamilySpecified':
                        familyName = Internal.getNodeFromType1(b, 'FamilyName_t')
                        if familyName is not None:
                            ft = Internal.getValue(familyName)
                            if nameType == 1: BCFZ += [name+Internal.SEP2+ft, n[1]]
                            else: BCFZ += [ft, n[1]]
                        else: BCFZ += [name, n[1]]
                    else:
                        if nameType == 1: BCFZ += [name+Internal.SEP2+BCtype, n[1]]
                        else: BCFZ += [name, n[1]]
        BCFaces.append(BCFZ)
    return BCFaces

#==============================================================================
# Ajoute les BCFaces dans l'arbre
#==============================================================================
def addBCFaces(t, BCFaces):
    tp = Internal.copyRef(t)
    _addBCFaces(t, BCFaces)
    return tp

def _addBCFaces(t, BCFaces):
    if BCFaces == []: return None
    zones = Internal.getZones(t)
    nz = 0; c = 0
    for z in zones:
        myBC = BCFaces[c]
        l = len(myBC)
        for i in range(l//2):
            name = myBC[2*i]
            names = name.split(Internal.SEP2)
            if len(names) == 2:
                name = names[0]; bctype = names[1]
                if bctype not in Internal.KNOWNBCS:
                    name = myBC[2*i]; bctype = 'FamilySpecified:'+name
            else: name = names[0]; bctype = 'FamilySpecified:'+names[0]
            faces = myBC[2*i+1]
            _addBC2Zone(z, name, bctype, faceList=faces)
        c += 1
    return None

# Ajoute les BCFields dans l'arbre
#def addBCFields(t, BCFields):
#    tp = Internal.copyRef(t)
#    _addBCFields(t, BCFields)
#    return tp
#
#def _addBCFields(t, BCFields):
#    if BCFields == []: return None
#    return None

#==============================================================================
# -- Fonctions de preconditionnement (hook) --
#==============================================================================

# -- createHook
def createHook(a, function='None'):
    """Create a hook for a given function.
      Usage: hook = createHook(a, function)"""
    if function == 'extractMesh': # growOfEps pas pret pour Array2
        fields = getFields(Internal.__GridCoordinates__, a, api=1)
    else:
        fields = getFields(Internal.__GridCoordinates__, a, api=3)
    if function == 'extractMesh' or function == 'adt':
        return Converter.createHook(fields, function)
    else:
        if len(fields) == 1: return Converter.createHook(fields[0], function)
        else: return Converter.createHook(fields, function)

def createHookAdtCyl(a, center=(0,0,0), axis=(0,0,1), depth=0, thetaShift=0.):
    """Create a hook for a cylindrical adt."""
    fields = getFields(Internal.__GridCoordinates__, a, api=2)
    return Converter.createHookAdtCyl(fields, center, axis, depth, thetaShift)

# -- createGlobalHook
def createGlobalHook(a, function='None', indir=0):
    """Create a global hook for all zones in a.
    Usage: hook = createGlobalHook(a, function)"""
    fields = getFields(Internal.__GridCoordinates__, a, api=3)
    return Converter.createGlobalHook(fields, function, indir)

# -- freeHook
def freeHook(hook):
    """Free hook.
    Usage: freeHook(hook)"""
    Converter.freeHook(hook)

#==============================================================================
# -- Fonctions d'identification geometrique --
#==============================================================================

# -- identifyNodes: identifie les noeuds de a dans hook
def identifyNodes(hook, a, tol=1.e-11):
    """Find in a hook nearest points of nodes of a. return identified node indices.
    Usage: identifyNodes(hook, a)"""
    fields = getFields(Internal.__GridCoordinates__, a, api=3)
    if len(fields) == 1: return Converter.identifyNodes(hook, fields[0], tol)
    else: return Converter.identifyNodes(hook, fields, tol)

# -- identifyFaces: identifie les centres de faces de a dans hook
def identifyFaces(hook, a, tol=1.e-11):
    """Find in a hook nearest points of face centers of a. return identified face indices.
    Usage: identifyFaces(hook, a)"""
    fields = getFields(Internal.__GridCoordinates__, a, api=1)
    if len(fields) == 1: return Converter.identifyFaces(hook, fields[0], tol)
    else: return Converter.identifyFaces(hook, fields, tol)

# -- identifyElements: identifie le centre des elements de a dans hook
def identifyElements(hook, a, tol=1.e-11):
    """Find in a hook nearest points of element centers of a. return identified element indices.
    Usage: identifyElements(hook, a)"""
    fields = getFields(Internal.__GridCoordinates__, a, api=3)
    if len(fields) == 1: return Converter.identifyElements(hook, fields[0], tol)
    else: return Converter.identifyElements(hook, fields, tol)

# -- identifySolutions: recopie la solution de tDnr dans tRcv par identification
# des noeuds et des centres
def identifySolutions(tRcv, tDnr, hookN=None, hookC=None, vars=[], tol=1.e6):
    """Identify points in stored in a global hook to mesh points and set the solution if donor
    and receptor points are distant from tol.
    Usage: identifySolutions(tRcv, tDnr, hookN, hookC, vars, tol)"""
    tp = Internal.copyRef(tRcv)
    _identifySolutions(tp, tDnr, hookN, hookC, vars, tol)
    return tp

def _identifySolutions(tRcv, tDnr, hookN=None, hookC=None, vars=[], tol=1.e6):
    """Identify points in stored in a global hook to mesh points and set the solution if donor
    and receptor points are distant from tol.
    Usage: identifySolutions(tRcv, tDnr, hookN, hookC, vars, tol)"""
    varsC=[]; varsN=[]
    for v in vars:
        s = v.find('centers:')
        if s != -1: varsC.append(v)
        else: varsN.append(v)

    if varsC == [] and hookC is not None:
        varsC = getVarNames(tDnr, excludeXYZ=True, loc='centers')
        if varsC != []: varsC = varsC[0]
    if varsN == [] and hookN is not None:
        varsN = getVarNames(tDnr, excludeXYZ=True, loc='nodes')
        if varsN != []: varsN = varsN[0]

    if len(varsC) == 0 and len(varsN) == 0: return None
    if varsC != []:
        for nov in range(len(varsC)):
            vc = varsC[nov].split(':')[1]
            varsC[nov] = vc

    zones = Internal.getZones(tRcv)
    coordsR = getFields(Internal.__GridCoordinates__, zones, api=1)

    if varsN != [] and hookN is not None:
        fnodes = getFields(Internal.__FlowSolutionNodes__, tDnr, api=1)
        resn = Converter.identifySolutions(coordsR, fnodes, hookN, vars=varsN, tol=tol)
        setFields(resn, zones, 'nodes')
    if varsC != [] and hookC is not None:
        fcenters = getFields(Internal.__FlowSolutionCenters__, tDnr, api=1)
        centersR = Converter.node2Center(coordsR)
        resc = Converter.identifySolutions(centersR, fcenters, hookC, vars=varsC, tol=tol)
        setFields(resc, zones, 'centers')
    return None

# -- nearestNodes: identifie le noeud de a le plus proche d'un point de hook
def nearestNodes(hook, a):
    """Identify nearest nodes to a in hook. return identified face indices.
    Usage: nearestNodes(hook, a)"""
    fields = getFields(Internal.__GridCoordinates__, a, api=3)
    if len(fields) == 1: return Converter.nearestNodes(hook, fields[0])
    else: return Converter.nearestNodes(hook, fields)

# -- nearestFaces: identifie la face de a la plus proches d'un point de hook
def nearestFaces(hook, a):
    """Identify nearest face centers to a in hook. return identified face indices.
    Usage: nearestFaces(hook, a)"""
    fields = getFields(Internal.__GridCoordinates__, a, api=1)
    if len(fields) == 1: return Converter.nearestFaces(hook, fields[0])
    else: return Converter.nearestFaces(hook, fields)

# -- nearestElements: identifie le centre de l'elements de a le plus proche
# d'un point de hook
def nearestElements(hook, a):
    """Identify nearest element centers to a in hook. return identified face indices.
    Usage: nearestFaces(hook, a)"""
    fields = getFields(Internal.__GridCoordinates__, a, api=3)
    if len(fields) == 1: return Converter.nearestElements(hook, fields[0])
    else: return Converter.nearestElements(hook, fields)

# Create global index
def createGlobalIndex(a, start=0):
    """Create the global index field."""
    return TZA2(a, 'nodes', 'nodes', Converter.createGlobalIndex, start)

def _createGlobalIndex(a, start=0):
    """Create the global index field."""
    _initVars(a, 'globalIndex', 0)
    return __TZA2(a, 'nodes', Converter._createGlobalIndex, start)

# Recover field from global index
def recoverGlobalIndex(a, b):
    """Recover fields of b in a following the global index field."""
    fb = getFields(Internal.__FlowSolutionNodes__, b, api=2)[0]
    return TZA2(a, 'nodes', 'nodes', Converter.recoverGlobalIndex, fb)

def _recoverGlobalIndex(a, b):
    """Recover fields of b in a following the global index field."""
    variables = getVarNames(b)[0]
    _addVars(a, variables)
    fx = getFields(Internal.__GridCoordinates__, b, api=2)[0]
    fb = getFields(Internal.__FlowSolutionNodes__, b, api=2)[0]
    if fx != [] and fb != []:
        fb[0] = fb[0]+','+fx[0]
        fb[1] = fb[1]+fx[1]
    elif fb == []: fb = fx
    return __TZA2(a, 'nodes', Converter._recoverGlobalIndex, fb)

#==============================================================================
# -- Connectivity management --
#==============================================================================

# -- mergeConnectivity
# IN: z1: zone BE
# IN: z2: zone BE (to be merged in z1) avec un seul type d'elements
# si boundary==1, les noeuds de z2 sont identifies dans z1
# si boundary==0, les noeuds de z2 sont merges dans z1 et reidentifies
# IN: boundary: 0 (not a boundary zone), 1 (a boundary zone, add it as a
# boundary connectivity)
def mergeConnectivity(z1, z2=None, boundary=0):
    """Gather an additional zone connectivity in z1.
    Usage: mergeConnectivity(z1, z2)"""
    zout = Internal.copyRef(z1)
    _mergeConnectivity(zout, z2, boundary)
    if z2 is None: return zout[0]
    return zout

def _mergeConnectivity(z1, z2=None, boundary=0, shared=False):
    if z2 is None:
    # Checking that z1 is a list of zones
        zones = Internal.getZones(z1)
        if len(zones) > 1:
        # Merge all BE connectivities within z1
            return _mergeConnectivities(zones, boundary, shared)
        else: return z1

    # Analyse zone z2
    dims = Internal.getZoneDim(z2)
    neb = dims[2] # nbre d'elts de z2
    eltType, _ = Internal.eltName2EltNo(dims[3]) # type d'elements de z2

    # On cherche l'element max dans les connectivites de z1
    maxElt = 0
    connects = Internal.getNodesFromType(z1, 'Elements_t')
    for cn in connects:
        r = Internal.getNodeFromName1(cn, 'ElementRange')
        m = r[1][1]
        maxElt = max(maxElt, m)

    # connectivite ajoutee=volumique non shared
    if boundary == 0 and not shared:
        # on fusionne les coordonnees
        import Transform.PyTree as T
        zn1 = convertArray2Node(z1)
        zn2 = convertArray2Node(z2)
        zn = T.join(zn1, zn2)
        # reset Coordinates in z1 with merged zn coordinates
        cont1 = Internal.getNodeFromName1(z1, Internal.__GridCoordinates__)
        contn = Internal.getNodeFromName1(zn, Internal.__GridCoordinates__)
        for name in ['CoordinateX', 'CoordinateY', 'CoordinateZ']:
            p1 = Internal.getNodeFromName1(cont1, name)
            pn = Internal.getNodeFromName1(contn, name)
            p1[1] = pn[1]
        # Recupere le container noeud
        cont1 = Internal.getNodeFromName1(z1, Internal.__FlowSolutionNodes__)
        contn = Internal.getNodeFromName1(zn, Internal.__FlowSolutionNodes__)
        if contn is not None:
            for n in contn[2]:
                if n[3] == 'DataArray_t':
                    p1 = Internal.getNodeFromName1(cont1, n[0])
                    if p1 is not None: p1[1] = n[1]

        # Nouveau nbre de points dans z1
        np = Internal.getZoneDim(zn)[1]
        z1[1] = numpy.copy(z1[1])
        z1[1][0,0] = np

        # Reidentifie les connectivites de z1
        hook = createHook(zn, 'nodes')
        ids = identifyNodes(hook, z1)
        nodes = Internal.getNodesFromType1(z1, 'Elements_t')
        for n in nodes:
            node = Internal.getNodeFromName1(n, 'ElementConnectivity')
            oldc = node[1]
            newc = numpy.copy(oldc)
            newc[:] = ids[oldc[:]-1]
            node[1] = newc

        # Ajoute les FlowSolutions en centres
        cont1 = Internal.getNodeFromName1(z1, Internal.__FlowSolutionCenters__)
        cont2 = Internal.getNodeFromName1(z2, Internal.__FlowSolutionCenters__)
        if cont2 is not None:
            for n in cont2[2]:
                if n[3] == 'DataArray_t':
                    p1 = Internal.getNodeFromName1(cont1, n[0])
                    if p1 is not None: p1[1] = numpy.concatenate((p1[1],n[1]))

        # Ajoute les connectivites de z2
        z1[1][0,1] += neb # nouveau nbre de cellules

        # on identifie les noeuds connectivity de z2 dans zn
        ids = identifyNodes(hook, z2)
        # On ajoute les connectivites de z2
        nodes = Internal.getNodesFromType1(z2, 'Elements_t')
        for n in nodes:
            node = Internal.createUniqueChild(z1, n[0]+'-2', 'Elements_t', value=n[1])
            rn = Internal.getNodeFromType1(n, 'IndexRange_t')
            r = Internal.getValue(rn)
            r0 = r[0]+maxElt; r1 = r[1]+maxElt
            Internal.createUniqueChild(node, rn[0], 'IndexRange_t', value=[r0,r1])

            oldc = Internal.getNodeFromName1(n, 'ElementConnectivity')[1]
            newc = numpy.copy(oldc)
            newc[:] = ids[oldc[:]-1]
            Internal.createUniqueChild(node, 'ElementConnectivity', 'DataArray_t', value=newc)

        # decale les ranges des zoneBC de z2
        zbc1 = Internal.getNodeFromType1(z1, 'ZoneBC_t')
        zbc2 = Internal.getNodesFromType2(z2, 'BC_t')
        if zbc1 is None and zbc2 != []: zbc1 = Internal.newZoneBC(parent=z1)
        for bc in zbc2:
            rn = Internal.getNodeFromType1(bc, 'IndexRange_t')
            r = Internal.getValue(rn)
            r0 = r[0,0]+maxElt; r1 = r[0,1]+maxElt
            n = Internal.createUniqueChild(zbc1, bc[0], 'BC_t', value=bc[1])
            Internal.createUniqueChild(n, rn[0], rn[3], value=[[r0,r1]])

        # Tri des connectivites
        Internal._sortNodesInZone(z1)

    elif boundary == 0 and shared: # connectivite ajoutee=volumique shared
        # on cree des nouveaux noeuds connectivites dans z1
        elts = Internal.getNodesFromType2(z2, 'Elements_t')
        z1[1][0,1] += neb # update le nbre d'elements de z1
        nebb = 0
        for e in elts:
            if e[1][1] == 0: # volumique uniquement
                r = Internal.getNodeFromName1(e, 'ElementRange')[1]
                nbe2 = r[1]-r[0]+1
                e2 = Internal.createUniqueChild(z1, e[0], 'Elements_t', value=[eltType,0])
                Internal.createUniqueChild(e2, 'ElementRange', 'IndexRange_t',
                                           value=[maxElt+nebb+1,maxElt+nebb+nbe2])
                newc = Internal.getNodeFromName1(e, 'ElementConnectivity')[1]
                Internal.createUniqueChild(e2, 'ElementConnectivity', 'DataArray_t', value=newc)
                nebb += nbe2

        # Ajoute les FlowSolutions en centres
        cont1 = Internal.getNodeFromName1(z1, Internal.__FlowSolutionCenters__)
        cont2 = Internal.getNodeFromName1(z2, Internal.__FlowSolutionCenters__)
        if cont2 is not None:
            for n in cont2[2]:
                if n[3] == 'DataArray_t':
                    p1 = Internal.getNodeFromName1(cont1, n[0])
                    if p1 is not None: p1[1] = numpy.concatenate((p1[1],n[1]))

        # Tri des connectivites
        Internal._sortNodesInZone(z1)


    else: # connectivite ajoutee=boundary (subzone)
        # on identifie les noeuds de z2 dans z1
        hook = createHook(z1, 'nodes')
        ids = identifyNodes(hook, z2)
        z1[1] = numpy.copy(z1[1])
        z1[1][0,2] += neb

        # on cree un nouveau noeud connectivite dans z1 (avec le nom de la zone z2)
        nebb = neb
        node = Internal.createUniqueChild(z1, z2[0], 'Elements_t', value=[eltType,nebb])
        Internal.createUniqueChild(node, 'ElementRange', 'IndexRange_t',
                                    value=[maxElt+1,maxElt+neb])
        oldc = Internal.getNodeFromName2(z2, 'ElementConnectivity')[1]
        newc = numpy.copy(oldc)
        newc[:] = ids[oldc[:]-1]
        Internal.createUniqueChild(node, 'ElementConnectivity', 'DataArray_t', value=newc)
    return None

def _mergeConnectivities(z, boundary=0, shared=False):
    if len(z) == 1: return z
    nebs = []; eltTypes = []
    for zone in z:
        dims = Internal.getZoneDim(zone)
        nebs.append(dims[2]) # nbre d'elts de zone
        eltType, _ = Internal.eltNames2EltNos(dims[3]) # types d'elements
        eltTypes.extend(eltType)
    maxElts = numpy.cumsum(nebs)
    nnewElts = maxElts[-1] - maxElts[0] # nbre d'elts ajoutes a z[0]

    # connectivite ajoutee=volumique non shared
    if boundary == 0 and not shared:
        # on fusionne les coordonnees
        import Transform.PyTree as T
        zn = convertArray2Node(z)
        zn = T.join(zn)
        # reset Coordinates in z[0] with merged zn coordinates
        cont1 = Internal.getNodeFromName1(z[0], Internal.__GridCoordinates__)
        contn = Internal.getNodeFromName1(zn, Internal.__GridCoordinates__)
        for name in ['CoordinateX', 'CoordinateY', 'CoordinateZ']:
            p1 = Internal.getNodeFromName1(cont1, name)
            pn = Internal.getNodeFromName1(contn, name)
            p1[1] = pn[1]
        # Recupere le container noeud
        cont1 = Internal.getNodeFromName1(z[0], Internal.__FlowSolutionNodes__)
        contn = Internal.getNodeFromName1(zn, Internal.__FlowSolutionNodes__)
        if contn is not None:
            for n in contn[2]:
                if n[3] == 'DataArray_t':
                    p1 = Internal.getNodeFromName1(cont1, n[0])
                    if p1 is not None: p1[1] = n[1]

        # Nouveau nbre de points dans z[0]
        np = Internal.getZoneDim(zn)[1]
        z[0][1] = numpy.copy(z[0][1])
        z[0][1][0,0] = np

        # Reidentifie les connectivites de z[0]
        hook = createHook(zn, 'nodes')
        ids = identifyNodes(hook, z[0])
        nodes = Internal.getNodesFromType1(z[0], 'Elements_t')
        for n in nodes:
            node = Internal.getNodeFromName1(n, 'ElementConnectivity')
            oldc = node[1]
            newc = numpy.copy(oldc)
            newc[:] = ids[oldc[:]-1]
            node[1] = newc

        # Ajoute les FlowSolutions en centres
        cont1 = Internal.getNodeFromName1(z[0], Internal.__FlowSolutionCenters__)
        if cont1 is not None:
            for i in range(1, len(z)):
                cont2 = Internal.getNodeFromName1(z[i], Internal.__FlowSolutionCenters__)
                for n in cont2[2]:
                    if n[3] == 'DataArray_t':
                        p1 = Internal.getNodeFromName1(cont1, n[0])
                        if p1 is not None: p1[1] = numpy.concatenate((p1[1],n[1]))

        # Ajoute les connectivites de z[1] a z[-1]
        z[0][1][0,1] += nnewElts # nouveau nbre de cellules

        # on identifie les noeuds connectivity de z[1] a z[-1] dans zn
        zbc1 = Internal.getNodeFromType1(z[0], 'ZoneBC_t')
        for i in range(1, len(z)):
            ids = identifyNodes(hook, z[i])
            # On ajoute les connectivites de z[i]
            nodes = Internal.getNodesFromType1(z[i], 'Elements_t')
            for n in nodes:
                node = Internal.createUniqueChild(z[0], n[0]+'-'+str(i+1), 'Elements_t',
                                                  value=n[1])
                rn = Internal.getNodeFromType1(n, 'IndexRange_t')
                r = Internal.getValue(rn)
                r0 = r[0]+maxElts[i]-maxElts[0]; r1 = r[1]+maxElts[i]-maxElts[0]
                Internal.createUniqueChild(node, rn[0], 'IndexRange_t', value=[r0,r1])

                oldc = Internal.getNodeFromName1(n, 'ElementConnectivity')[1]
                newc = numpy.copy(oldc)
                newc[:] = ids[oldc[:]-1]
                Internal.createUniqueChild(node, 'ElementConnectivity', 'DataArray_t',
                                           value=newc)

            # decale les ranges des zoneBC de z[i]
            zbc2 = Internal.getNodesFromType2(z[i], 'BC_t')
            if zbc1 is None and zbc2 != []: zbc1 = Internal.newZoneBC(parent=z[0])
            for bc in zbc2:
                rn = Internal.getNodeFromType1(bc, 'IndexRange_t')
                r = Internal.getValue(rn)
                r0 = r[0,0]+maxElts[i]-maxElts[0]; r1 = r[0,1]+maxElts[i]-maxElts[0]
                n = Internal.createUniqueChild(zbc1, bc[0], 'BC_t', value=bc[1])
                Internal.createUniqueChild(n, rn[0], rn[3], value=[[r0,r1]])

        # Tri des connectivites
        Internal._sortNodesInZone(z[0])

    elif boundary == 0 and shared: # connectivite ajoutee=volumique shared
        # on cree des nouveaux noeuds connectivites dans z[0]
        z[0][1][0,1] += nnewElts # update le nbre d'elements de z[0]
        nebb = maxElts[0]
        for i in range(1, len(z)):
            elts = Internal.getNodesFromType2(z[i], 'Elements_t')
            for e in elts:
                if e[1][1] == 0: # volumique uniquement
                    r = Internal.getNodeFromName1(e, 'ElementRange')[1]
                    nbe2 = r[1]-r[0]+1
                    e2 = Internal.createUniqueChild(z[0], e[0], 'Elements_t',
                                                    value=[eltTypes[i],0])
                    Internal.createUniqueChild(e2, 'ElementRange', 'IndexRange_t',
                                               value=[nebb+1,nebb+nbe2])
                    newc = Internal.getNodeFromName1(e, 'ElementConnectivity')[1]
                    Internal.createUniqueChild(e2, 'ElementConnectivity', 'DataArray_t',
                                               value=newc)
                    nebb += nbe2

        # Ajoute les FlowSolutions en centres
        cont1 = Internal.getNodeFromName1(z[0], Internal.__FlowSolutionCenters__)
        if cont1 is not None:
            for i in range(1, len(z)):
                cont2 = Internal.getNodeFromName1(z[i], Internal.__FlowSolutionCenters__)
                for n in cont2[2]:
                    if n[3] == 'DataArray_t':
                        p1 = Internal.getNodeFromName1(cont1, n[0])
                        if p1 is not None: p1[1] = numpy.concatenate((p1[1],n[1]))

        # Tri des connectivites
        Internal._sortNodesInZone(z[0])

    else: # connectivite ajoutee=boundary (subzone)
        z[0][1] = numpy.copy(z[0][1])
        z[0][1][0,2] += nnewElts
        # on identifie les noeuds de z[i] dans z[0]
        hook = createHook(z[0], 'nodes')
        for i in range(1, len(z)):
            ids = identifyNodes(hook, z[i])
            # on cree un nouveau noeud connectivite dans z[0] (avec le nom de la zone z[i])
            node = Internal.createUniqueChild(z[0], z[i][0], 'Elements_t',
                                              value=[eltTypes[i],nebs[i]])
            Internal.createUniqueChild(node, 'ElementRange', 'IndexRange_t',
                                       value=[maxElts[i-1]+1,maxElts[i]])
            oldc = Internal.getNodeFromName2(z[i], 'ElementConnectivity')[1]
            newc = numpy.copy(oldc)
            newc[:] = ids[oldc[:]-1]
            Internal.createUniqueChild(node, 'ElementConnectivity', 'DataArray_t',
                                       value=newc)
    return None

#============================================
# Soit z une sous zone de zt
# Remap la connectivite de z sur celle de zt
# Share nodes field a la fin
#============================================
def _mapConnectOnZone(z, zt):
    hook = createHook(zt, 'nodes')
    ids = identifyNodes(hook, z)
    freeHook(hook)
    nodes = Internal.getNodesFromType(z, 'Elements_t')
    for n in nodes:
        node = Internal.getNodeFromName(n, 'ElementConnectivity')
        oldc = node[1]
        newc = numpy.copy(oldc)
        newc[:] = ids[oldc[:]-1]
        node[1] = newc
    gc = Internal.getNodeFromName1(zt, Internal.__GridCoordinates__)
    if gc is not None:
        fields = Internal.getNodesFromType1(gc, 'DataArray_t')
        for f in fields:
            Internal.getNodeFromName2(z, f[0])[1] = f[1]
    gc = Internal.getNodeFromName1(zt, Internal.__FlowSolutionNodes__)
    if gc is not None:
        fields = Internal.getNodesFromType1(gc, 'DataArray_t')
        for f in fields:
            Internal.getNodeFromName2(z, f[0])[1] = f[1]
    z[1][0,0] = zt[1][0,0]
    return None

# -- breakConnectivity
# break a multiple connectivity zone into single elements zones
# don't break boundary connectivity
# IN: t: t to break
def breakConnectivity(t):
    """Break a multi-element zone into single element zones.
    Usage: breakConnectivity(t)"""
    tp, typen = Internal.node2PyTree(t)
    bases = Internal.getBases(tp)
    for b in bases:
        c = 0; l = len(b[2])
        for c in range(l):
            z = b[2][c]
            if z[3] == 'Zone_t':
                            # compte les connectivites elements (hors boundary)
                connects = Internal.getElementNodes(z)
                N = len(connects)
                if N <= 1: break # une seule connectivite
                if N == 2:
                    type1 = connects[0][1][0]; type2 = connects[1][1][0]
                    if (type1 == 22 and type2 == 23) or (type1 == 23 and type2 == 22): # pur NGON
                        break

                iBE = []; iNGon = -1; iNFace = -1; i = 0
                for co in connects:
                    ctype = co[1][0]
                    if ctype == 22: iNGon = i
                    elif ctype == 23: iNFace = i
                    else: iBE.append(i)
                    i += 1

                N = len(iBE)
                # split les connectivites volumiques
                for p in range(N):
                    i = iBE[p]
                    zp = Internal.copyRef(z)
                    _deleteGridConnectivity__(zp)
                    GEl = Internal.getNodesFromType1(zp, 'Elements_t')
                    GE = Internal.getNodeFromName1(zp, connects[i][0])
                    eltType, nf = Internal.eltNo2EltName(GE[1][0])
                    # Nouveau nom de la zone
                    zp[0] = getZoneName(z[0]+'_'+eltType)
                    # Enleve toutes les connects volumiques a part la ieme
                    for GEj in GEl:
                        if GEj is not GE and GEj[1][1] == 0: Internal._rmNodesByName(zp, GEj[0])
                    # Renumerote la connectivite
                    r = Internal.getNodeFromName(GE, 'ElementRange'); r = r[1]
                    start = r[0]; end = r[1]
                    # modifie le range de la connectivite
                    r[0] = 1; r[1] = end-start+1
                    # modifie le nbre d'elements dans la zone
                    zp[1] = numpy.copy(z[1])
                    zp[1][0,1] = end-start+1
                    if i == 0: b[2][c] = zp
                    else: b[2] += [zp]
                    # On slice le champ en centres (shared)
                    cont = Internal.getNodeFromName1(zp, Internal.__FlowSolutionCenters__)
                    if cont is not None:
                        for v in cont[2]:
                            if v[3] == 'DataArray_t': v[1] = v[1][start-1:end] # view
                    # On modifie le range des connectivites BCs
                    Internal._sortNodesInZone(zp)

                if iNGon != -1 and iNFace != -1: # NGon additionnel
                    i1 = iNGon; i2 = iNFace
                    zp = Internal.copyRef(z)
                    _deleteFlowSolutions__(zp, 'centers')
                    _deleteGridConnectivity__(zp)
                    _deleteZoneBC__(zp)
                    GEl = Internal.getNodesFromType1(zp, 'Elements_t')
                    GE1 = Internal.getNodeFromName1(zp, connects[i1][0])
                    GE2 = Internal.getNodeFromName1(zp, connects[i2][0])
                    zp[0] = getZoneName(z[0]+'_NGON')
                    # Enleve toutes les connects a part la ieme
                    for GEj in GEl:
                        if GEj is not GE1 and GEj is not GE2: Internal._rmNodesByName(zp, GEj[0])
                    # Renumerote la connectivite
                    r = Internal.getNodeFromName(GE1, 'ElementRange'); r = r[1]
                    start = r[0]; end = r[1]
                    # modifie le range de la connectivite
                    r[0] = 1; r[1] = end-start+1; fin = end-start+1
                    r = Internal.getNodeFromName(GE2, 'ElementRange'); r = r[1]
                    start = r[0]; end = r[1]
                    # modifie le range de la connectivite
                    r[0] = fin; r[1] = fin+end-start+1
                    # modifie le nbre d'elements dans la zone
                    zp[1] = numpy.copy(z[1])
                    zp[1][0,1] = end-start+1
                    b[2] += [zp]
                    #zp = pushBC(z, zp, type='F')

    if typen == 1: typen = 2 # une zone renvoie une liste de zones
    return Internal.pyTree2Node(tp, typen)

#==============================================================================
# renumerote les ElementRanges pour qu'ils soient dans un ordre contigue
# y compris pour les BCCs
#==============================================================================
def _renumberElementConnectivity(t):
    zones = Internal.getZones(t)
    for z in zones:
        elts = Internal.getNodesFromType1(z, 'Elements_t')
        c = 1
        for e in elts:
            r = Internal.getNodeFromName1(e, 'ElementRange')
            r[1] = numpy.copy(r[1]); r = r[1]
            delta = r[1]-r[0]+1
            r[0] = c; r[1] = c+delta-1; c += delta
    return None

# -- selectOneConnectivity
# Retourne une nouvelle zone avec une seule de ses connectivites
# (passee en connectivite volumique)
# IN: name: nom de la ElementConnectivity a conserver
# or IN: number: no de la ElementConnectivity a conserver (first=0)
# or IN: irange=[min,max]
def selectOneConnectivity(z, name=None, number=None, irange=None):
    zp = Internal.copyRef(z)
    _selectOneConnectivity(zp, name=name, number=number, irange=irange)
    return zp

def _selectOneConnectivity(zp, name=None, number=None, irange=None):
    elts = Internal.getNodesFromType1(zp, 'Elements_t')

    if name is not None:
        for e in elts:
            if e[0] != name: Internal._rmNodesByName(zp, e[0])
            else: e[1] = numpy.copy(e[1]); e[1][1] = 0 # force volumique
    elif number is not None:
        c = 0
        for e in elts:
            if c != number: Internal._rmNodesByName(zp, e[0])
            else: e[1] = numpy.copy(e[1]); e[1][1] = 0 # force volumique
            c += 1
    elif irange is not None:
        for e in elts:
            r = Internal.getNodeFromName1(e, 'ElementRange')
            if r is not None:
                r = r[1]
                if r[0] > irange[0] or r[1] < irange[1]:
                    Internal._rmNodesByName(zp, e[0])
                else:
                    if r[0] == irange[0] and r[1] == irange[1]: # full
                        e[1] = numpy.copy(e[1]); e[1][1] = 0 # force volumique
                    else: # slice
                        (name, nnodes) = Internal.eltNo2EltName(e[1][0])
                        if name != 'NGON' and name != 'NFACE' and name != 'MIXED':
                            e[1] = numpy.copy(e[1]); e[1][1] = 0 # force volumique
                            r = Internal.getNodeFromName1(e, 'ElementRange')
                            r[1] = numpy.copy(r[1]); r[1][0] = 1; r[1][1] = irange[1]-irange[0]+1
                            c = Internal.getNodeFromName1(e, 'ElementConnectivity')
                            c[1] = c[1][nnodes*(irange[0]-1):nnodes*(irange[1])+1]
                        else: print ('Warning: selectOneConnectivity: slice impossible.')
            else: Internal._rmNodesByName(zp, e[0])
    _renumberElementConnectivity(zp)
    return None

# -- selectConnectivity
# Retourne une nouvelle zone avec une seule de ses connectivites
# (passee en connectivite volumique)
# IN: name: nom de la ElementConnectivity a conserver
# or IN: number: no de la ElementConnectivity a conserver (first=0)
# or IN: irange=[min,max]
def selectConnectivity(z, name=None, number=None, irange=None):
    zp = Internal.copyRef(z)
    elts = Internal.getNodesFromType1(zp, 'Elements_t')

    if name is not None:
        for e in elts:
            if e[0] != name: Internal._rmNodesByName(zp, e[0])
            else: e[1] = numpy.copy(e[1]); e[1][1] = 0 # force volumique
    elif number is not None:
        c = 0
        for e in elts:
            if c != number: Internal._rmNodesByName(zp, e[0])
            else: e[1] = numpy.copy(e[1]); e[1][1] = 0 # force volumique
            c += 1
    elif irange is not None:
        for e in elts:
            r = Internal.getNodeFromName1(e, 'ElementRange')
            if r is not None:
                r = r[1]
                if r[0] != irange[0] and r[1] != irange[1]:
                    Internal._rmNodesByName(zp, e[0])
                if r[0] == irange[0] and r[1] == irange[1]: # full
                    e[1] = numpy.copy(e[1]); e[1][1] = 0 # force volumique
                if r[0] == irange[0] and r[1] < irange[1]: # full
                    (name, nnodes) = Internal.eltNo2EltName(e[1][0])
                    if name != 'NGON' and name != 'NFACE' and name != 'MIXED':
                        e[1] = numpy.copy(e[1]); e[1][1] = 0 # force volumique
                        r0 = r[0]; r1 = r[1]
                        r = Internal.getNodeFromName1(e, 'ElementRange')
                        r[1] = numpy.copy(r[1]); r[1][0] = 1; r[1][1] = r1-r0+1
                        c = Internal.getNodeFromName1(e, 'ElementConnectivity')
                        r = r[1]
                        c[1] = c[1][nnodes*(r[0]-1):nnodes*(r[1])+1]
                    else: print('Warning: selectConnectivity: slice impossible.')
                if r[0] > irange[0] and r[1] == irange[1]: # full
                    (name, nnodes) = Internal.eltNo2EltName(e[1][0])
                    if name != 'NGON' and name != 'NFACE' and name != 'MIXED':
                        e[1] = numpy.copy(e[1]); e[1][1] = 0 # force volumique
                        r0 = r[0]; r1 = r[1]
                        r = Internal.getNodeFromName1(e, 'ElementRange')
                        r[1] = numpy.copy(r[1]); r[1][0] = 1; r[1][1] = r1-r0+1
                        c = Internal.getNodeFromName1(e, 'ElementConnectivity')
                        r = r[1]
                        c[1] = c[1][nnodes*(r[0]-1):nnodes*(r[1])+1]
                    else: print('Warning: selectConnectivity: slice impossible.')
            else: Internal._rmNodesByName(zp, e[0])
    _renumberElementConnectivity(zp)
    return zp

#=============================================================================
# Rm duplicated periodic zones in a according to the node 'TempPeriodicZone'
#=============================================================================
def removeDuplicatedPeriodicZones__(a):
    t = Internal.copyRef(a)
    _removeDuplicatedPeriodicZones__(t)
    return t

def _removeDuplicatedPeriodicZones__(a):
    for z in Internal.getZones(a):
        parent,d = Internal.getParentOfNode(a, z)
        isperiod = Internal.getNodeFromName1(z, 'TempPeriodicZone')
        if isperiod is not None: del parent[2][d]
    return None

#=============================================================================
# Extract periodic zone info and duplicate zone in same base
#=============================================================================
def addPeriodicZones__(a):
    t = Internal.copyRef(a)
    _addPeriodicZones__(t)
    return t

def _addPeriodicZones__(a):
    try: import Transform.PyTree as T
    except: raise ImportError("addPeriodicZones: requires Transform module.")
    atype = Internal.typeOfNode(a)
    if atype != 4:  # base
        print('Warning: addPeriodicZones: input node must be a CGNS basis.')
        print('Skipped.')
        return None

    zones = Internal.getNodesFromType1(a, 'Zone_t')
    zonesdup = []
    for z in zones:
        zname = Internal.getName(z)
        usd = Internal.getNodeFromName2(z, '.Solver#Param')
        periodicChimera = False
        if usd is not None:
            perdir = Internal.getNodeFromName1(usd, 'periodic_dir')
            if perdir is not None: periodicChimera=True

        # Periodic match info
        if not periodicChimera:
            gcmatch = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
            for gc in gcmatch:
                rotationData, translVect = Internal.getPeriodicInfo__(gc)

                if translVect == []: translated = False
                else: translated = True

                if rotationData == []: rotated = False
                else: rotated = True
                if rotated and translated:
                    print('Warning: addPeriodicZones: rotation and translation cannot be applied at the same time. %s periodic grid connectivity not taken into account.'%gc[0])
                elif not rotated and not translated: pass
                else:
                    zdonorname = Internal.getValue(gc)
                    zd = Internal.getNodeFromName2(a, zdonorname)
                    if zd is not None:
                        if rotated:
                            xc = rotationData[0]; yc = rotationData[1]; zc = rotationData[2]
                            vx = rotationData[3]; vy = rotationData[4]; vz = rotationData[5]
                            angle = rotationData[6]
                            if angle != 0.:
                                zddup = T.rotate(zd,(xc,yc,zc), (vx,vy,vz),-angle)
                                if angle >0.: signangle = -1
                                else: signangle = 1
                                rotationInfo = Internal.createNode('SignOfRotationAngle', 'UserDefinedData_t', value=signangle)
                                # creation du noeud temporaire le marquant comme periodique
                                zddup[0] = getZoneName(zddup[0]+'_dup')
                                Internal.createChild(zddup, 'TempPeriodicZone', 'UserDefinedData_t', value=zdonorname, children=[rotationInfo])
                                zonesdup.append(zddup)
                        elif translated:
                            tvx = translVect[0]; tvy = translVect[1]; tvz = translVect[2]
                            zddup = T.translate(zd,(-tvx,-tvy,-tvz))
                            # creation du noeud temporaire le marquant comme periodique
                            zddup[0] = getZoneName(zddup[0]+'_dup')
                            Internal.createChild(zddup, 'TempPeriodicZone', 'UserDefinedData_t', value=zdonorname, children=None)
                            zonesdup.append(zddup)

        # Chimere periodique: compatible avec elsA uniquement
        else:
            print("Info: Periodic Chimera for zone %s"%(z[0]))
            perdir = Internal.getNodeFromName1(usd,'periodic_dir')
            if perdir is not None:
                perdir = Internal.getValue(perdir)
                ang1 = Internal.getNodeFromName1(usd,'axis_ang_1'); ang1 = Internal.getValue(ang1)
                ang2 = Internal.getNodeFromName1(usd,'axis_ang_2'); ang2 = Internal.getValue(ang2)
                angle = float(ang2)/float(ang1)*360.
                xc = Internal.getNodeFromName1(usd,'axis_pnt_x'); xc = Internal.getValue(xc)
                yc = Internal.getNodeFromName1(usd,'axis_pnt_y'); yc = Internal.getValue(yc)
                zc = Internal.getNodeFromName1(usd,'axis_pnt_z'); zc = Internal.getValue(zc)
                vx = Internal.getNodeFromName1(usd,'axis_vct_x'); vx = Internal.getValue(vx)
                vy = Internal.getNodeFromName1(usd,'axis_vct_y'); vy = Internal.getValue(vy)
                vz = Internal.getNodeFromName1(usd,'axis_vct_z'); vz = Internal.getValue(vz)
                if angle != 0.:
                    if perdir == 1 or perdir == 3:
                        zdup = T.rotate(z,(xc,yc,zc), (vx,vy,vz),angle)
                        signangle = 1
                        rotationInfo = Internal.createNode('SignOfRotationAngle', 'UserDefinedData_t', value=signangle)
                        # creation du noeud temporaire le marquant comme periodique
                        #zdup[0] = getZoneName(zdup[0]+'_dup')
                        pref = zdup[0].split('_')
                        if len(pref) == 1: pref = pref[0]
                        else: pref = pref[0]+'_'+pref[1]
                        zdup[0] = getZoneName(pref+'_dup')
                        Internal.createChild(zdup, 'TempPeriodicZone', 'UserDefinedData_t', value=zname, children=[rotationInfo])
                        zonesdup.append(zdup)
                    if perdir == 2 or perdir == 3:
                        zdup = T.rotate(z,(xc,yc,zc), (vx,vy,vz),-angle)
                        signangle =-1
                        rotationInfo = Internal.createNode('SignOfRotationAngle', 'UserDefinedData_t', value=signangle)
                        # creation du noeud temporaire le marquant comme periodique
                        zdup[0] = getZoneName(zdup[0]+'_dup')
                        pref = zdup[0].split('_')
                        if len(pref) == 1: pref = pref[0]
                        else: pref = pref[0]+'_'+pref[1]
                        zdup[0] = getZoneName(pref+'_dup')

                        Internal.createChild(zdup, 'TempPeriodicZone', 'UserDefinedData_t', value=zname, children=[rotationInfo])
                        zonesdup.append(zdup)
    a[2] += zonesdup
    return None

#==============================================================================
def convertPyTree2FFD(zone, RefStat, FlowEq, nd):
    Converter.converter.convertPyTree2FFD(zone, RefStat, FlowEq, nd,
                                          Internal.__GridCoordinates__,
                                          Internal.__FlowSolutionNodes__,
                                          Internal.__FlowSolutionCenters__)
    return None

def signNGonFaces(t):
    """Sign NFACE connectivity in NGON zones."""
    tc = Internal.copyRef(t)
    _signNGonFaces(tc)
    return tc

def _signNGonFaces(t):
    """Sign NFACE connectivity in NGON zones."""
    __TZGC3(t, Converter._signNGonFaces)
    return None

def unsignNGonFaces(t):
    """Unsign NFACE connectivity in NGON zones."""
    tc = Internal.copyRef(t)
    _unsignNGonFaces(tc)
    return tc

def _unsignNGonFaces(t):
    """Unsign NFACE connectivity in NGON zones and create a node storing whether
    NFACE was initially signed.
    Usage: _unsignNGonFaces(a)"""
    zones = Internal.getZones(t)
    isSigned = -1
    for z in zones:
        n0 = Internal.getNodeFromName1(z, 'NFaceElements')
        if n0 is None: return None # not an NGon
        n = Internal.getNodeFromName1(n0, 'ElementStartOffset')
        if n is None: return None # not an NGon v4
        n = Internal.getNodeFromName1(n0, 'Signed')
        if n is not None: return None # function already called
        fc = getFields(Internal.__GridCoordinates__, z, api=3)[0]
        if isSigned != 0 and fc:
            isSigned = Converter._unsignNGonFaces(fc)
        if isSigned != -1:
            Internal._createUniqueChild(n0, 'Signed', 'DataArray_t', value=isSigned)
    return None

def resignNGonFaces(t):
    """Resign NFACE connectivity in NGON zones if the node `Signed` is present
    and has a value of 1 (meaning the NGonv4 was signed when first read)."""
    tc = Internal.copyRef(t)
    _resignNGonFaces(tc)
    return tc

def _resignNGonFaces(t):
    """Resign NFACE connectivity in NGON zones and delete node `Signed`."""
    z = Internal.getZones(t)[0]
    n0 = Internal.getNodeFromName1(z, 'NFaceElements')
    if n0 is None: return None # not an NGon
    n = Internal.getNodeFromName1(n0, 'ElementStartOffset')
    if n is None: return None # not an NGon v4
    n = Internal.getNodeFromName1(n0, 'Signed')
    if n is None: return None # NGon was not signed
    Internal._rmNodesFromName(n0, 'Signed')
    if Internal.getValue(n) == 1: _signNGonFaces(t)
    return None

def makeParentElements(t):
    tc = Internal.copyRef(t)
    _makeParentElements(tc)
    return tc

def _makeParentElements(t):
    zones = Internal.getZones(t)
    for z in zones:
        fc = getFields(Internal.__GridCoordinates__, z, api=3)[0]
        if fc != []:
            pe = Converter.makeParentElements(fc)
            Internal._createUniqueChild(z, 'ParentElements', 'DataArray_t', pe[0])
    return None

# Convert to low order mesh
def convertHO2LO(t, mode=0):
    """Convert a HO element mesh to linear mesh.
    Usage: convertHO2LO(t, mode)"""
    return TZGC2(t, 'nodes', True, Converter.convertHO2LO, mode)

# Convert to high order mesh
def convertLO2HO(t, mode=0, order=2):
    """Convert a LO element mesh to high order mesh.
    Usage: convertLO2HO(t, mode, order)"""
    return TZGC2(t, 'nodes', True, Converter.convertLO2HO, mode, order)

def convertMIXED2NGon(a, recoverBC=True, merged=False):
    """Convert a mixed-element monozone to an NGON.
    Usage: convertMIXED2NGon(a, recoverBC, merged)"""
    zones = Internal.getZones(a)
    # Extrait les BCs comme des zones
    AllBCs = []; AllBCNames = []; AllBCTypes = []
    for z in zones:
        dimZ = Internal.getZoneDim(z)
        if dimZ[3] != 'MIXED': pass
        zoneBC = Internal.getNodeFromType1(z, 'ZoneBC_t')
        bcs = Internal.getNodesFromType1(zoneBC, 'BC_t')
        for b in bcs:
            myZone = Internal.rmNodesFromType(z, 'Elements_t')
            ntype = Internal.getValue(b)
            # Trouve la connectivite avec le nom de la BC (!)
            split = b[0].split('-')
            namebc = split[0]
            bctype = split[-1]

            if len(split)>2:
                for nosuff in range(1,len(split)-1): namebc+='-%s'%(split[nosuff])

            ebc = Internal.getNodeFromName1(z, namebc)
            if ebc is not None and Internal.getValue(ebc)[0]==20:
                p = Internal.getNodeFromName1(ebc, 'ElementConnectivity')
                out = Converter.converter.convertMix2BE(p[1])
                if out[0] is not None:
                    p = out[0]
                    e = Internal.newElements(ebc[0]+'_BAR', etype='BAR', econnectivity=p, erange=[1,p.size/2], parent=myZone)
                if out[1] is not None:
                    p = out[1]
                    e = Internal.newElements(ebc[0]+'_TRI', etype='TRI', econnectivity=p, erange=[1,p.size/3], parent=myZone)
                if out[2] is not None:
                    p = out[2]
                    e = Internal.newElements(ebc[0]+'_QUAD', etype='QUAD', econnectivity=p, erange=[1,p.size/4], parent=myZone)
                Internal._rmNodesFromName(myZone, 'ZoneBC')
                Internal._updateElementRange(myZone)
                znamebc = myZone[0]
                myZone = breakConnectivity(myZone)

                for zbc in Internal.getZones(myZone):
                    zbc[0] = getZoneName(znamebc)
                    Internal.createNode("BCName","UserDefinedData_t",value=b[0], parent=zbc)
                    Internal.createNode("BCType","UserDefinedData_t",value=bctype, parent=zbc)
                    AllBCs.append(zbc)

    Internal._rmNodesFromType(a, 'ZoneBC_t')
    tb = newPyTree(["BCs"]); tb[2][1][2] = Internal.getZones(AllBCs)
    importG = False
    try:
        import Generator.PyTree as G
        tb = G.close(tb,tol=1e-6)
        importG = True
    except:
        pass

    # Convertit les connectivites MIXED volumique
    for z in Internal.getZones(a):
        nodes = Internal.getNodesFromType(z, 'Elements_t')
        for n in nodes:
            p = Internal.getNodeFromName1(n, 'ElementConnectivity')
            out = Converter.converter.convertMix2BE(p[1])
            # Replace et renumerotes
            if out[0] is not None:
                Internal._rmNode(z,p[0])
            if out[1] is not None:
                Internal._rmNode(z,p[0])
            if out[2] is not None:
                Internal._rmNode(z,p[0])
            if out[3] is not None:
                p = out[3]
                e = Internal.newElements(n[0]+'_TETRA', etype='TETRA', econnectivity=p, erange=[1,p.size/4], parent=z)
            if out[4] is not None:
                p = out[4]
                e = Internal.newElements(n[0]+'_PYRA', etype='PYRA', econnectivity=p, erange=[1,p.size/5], parent=z)
            if out[5] is not None:
                p = out[5]
                e = Internal.newElements(n[0]+'_PENTA', etype='PENTA', econnectivity=p, erange=[1,p.size/6], parent=z)
            if out[6] is not None:
                p = out[6]
                e = Internal.newElements(n[0]+'_HEXA', etype='HEXA', econnectivity=p, erange=[1,p.size/8], parent=z)

            Internal._rmNodesFromName(z, n[0])
        Internal._updateElementRange(z)

    # Join
    a = breakConnectivity(a)
    a = convertArray2NGon(a)

    if merged:
        try:
            import Transform.PyTree as T
            a = T.join(a)
        except:
            pass
    if recoverBC:
        Internal._rmNodesFromType(a,"ZoneBC_t")
        nzones = len(Internal.getZones(a))
        AllBCs = []; AllBCNames = []; AllBCTypes = []
        dictOfBCsPerBCName={}
        for zs in Internal.getZones(tb):
            bcname = Internal.getNodeFromName(zs, 'BCName')
            bcname = Internal.getValue(bcname)

            if bcname not in dictOfBCsPerBCName:
                dictOfBCsPerBCName[bcname]=[zs]
            else:
                #check connect
                zsref=dictOfBCsPerBCName[bcname][0]
                ecref = Internal.getNodeFromType(zsref,'Elements_t')
                ecref = Internal.getNodeFromName(ecref,'ElementConnectivity')
                ecref = Internal.getValue(ecref)
                # to compare with
                ec = Internal.getNodeFromType(zs,'Elements_t')
                ec = Internal.getNodeFromName(ec,'ElementConnectivity')
                ec = Internal.getValue(ec)
                if not numpy.array_equal(ecref, ec):
                    dictOfBCsPerBCName[bcname].append(zs)
                else:
                    print("WARNING: deux BCs identiques : %s : %s et %s !"%(bcname, zsref[0], zs[0]))

        for noz, z in enumerate(Internal.getZones(a)):
            AllBCs = []; AllBCNames = []; AllBCTypes = []
            for zsname in dictOfBCsPerBCName:
                for zs in dictOfBCsPerBCName[zsname]:
                    intersect = 1
                    if importG: intersect = G.bboxIntersection(z,zs)

                    if intersect==1:
                        bcname =Internal.getValue(Internal.getNodeFromName(zs,"BCName"))
                        AllBCNames.append(bcname)
                        AllBCTypes.append('FamilySpecified:'+bcname)
                        Internal._rmNodesFromType(z,"UserDefinedData_t")
                        AllBCs.append(zs)
            _recoverBCs(z, (AllBCs, AllBCNames, AllBCTypes))

        for base in Internal.getBases(a):
            for famname in Internal.getNodesFromType(base,'FamilyName_t'):
                name = Internal.getValue(famname)
                _addFamily2Base(base, name, bndType='BCWall')

    return a

# Ghost cells
def _addGhostCells(t, b, d, adaptBCs=1, modified=[], fillCorner=1):
    """Add ghost cells to pyTree."""
    Internal._addGhostCells(t, b, d, adaptBCs, modified, fillCorner)
    return None

def addGhostCells(t, b, d, adaptBCs=1, modified=[], fillCorner=1):
    """Add ghost cells to pyTree."""
    return  Internal.addGhostCells(t, b, d, adaptBCs, modified, fillCorner)

def rmGhostCells(t, b, d, adaptBCs=1, modified=[]):
    """Remove ghost cells to a pyTree."""
    return Internal.rmGhostCells(t, b, d, adaptBCs, modified)

def _rmGhostCells(t, b, d, adaptBCs=1, modified=[]):
    """Remove ghost cells to a pyTree."""
    Internal._rmGhostCells(t, b, d, adaptBCs, modified)
    return None

#==============================================================================
# - client/server -
#==============================================================================
def send(data, host, rank=0, port=15555):
    """Send data to server."""
    Converter.send(data, host, rank, port)

def createSockets(nprocs=1, port=15555):
    """Create sockets."""
    return Converter.createSockets(nprocs, port)

def listen(s):
    """Listen to one socket."""
    return Converter.listen(s)

