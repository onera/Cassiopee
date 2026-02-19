"""Read files with filters."""
from . import Internal
from . import PyTree
from . import Converter
from . import Mpi as Cmpi
from .Distributed import convert2PartialTree, _convert2PartialTree, convert2SkeletonTree, _convert2SkeletonTree, convertFile2SkeletonTree, \
    _readPyTreeFromPaths, readPyTreeFromPaths, _readZones, readNodesFromPaths, fixPaths__
from . import Distributed
import Compressor.PyTree as Compressor
import numpy

# Prend un fileName, si c'est toto/*, rend la liste des fichiers
def expand(fileName):
    s = fileName.rsplit('/', 1)
    if len(s) > 1 and (s[1] == '*' or s[1] == '*.cgns' or s[1] == '*.hdf'): # multiple
        import glob
        if s[1] == '*': files = glob.glob(s[0]+'/*.cgns')
        else: files = glob.glob(fileName)
    return files

# version collective/serialisee
def writeNodesFromPaths(fileName, paths, nodes, format=None, maxDepth=-1, mode=0):
    """Write nodes to file given their paths."""
    Cmpi.seq(Distributed.writeNodesFromPaths, fileName, paths, nodes, format, maxDepth, mode)
    return None

# version collective/serialisee
def writePyTreeFromPaths(t, fileName, paths, format=None, maxDepth=-1):
    """Write some nodes of the pyTree given their paths."""
    Cmpi.seq(Distributed.writePyTreeFromPaths, t, fileName, paths, format, maxDepth)
    return None

# version collective/serialisee
def deletePaths(fileName, paths, format=None):
    """Delete nodes in file given their paths."""
    Cmpi.seq(Distributed.deletePaths, fileName, paths, format)
    return None

#==============================================================================
# Lit seulement une partie des tableaux d'un fichier a partir de la definition d'un filtre.
# filter est un dictionnaire pour chaque path
# pour les grilles structurees : [[imin,jmin,kmin], [1,1,1], [imax,jmax,kmax], [1,1,1]]
# pour les grilles non structurees : [[istart], [1], [iend], [1]]
# Uniquement HDF
# Retourne un dictionnaire du numpy loade pour chaque path
#==============================================================================
def readNodesFromFilter(fileName, filter, format='bin_hdf', com=None):
    """Read nodes from file given a filter."""
    filter2 = {}
    for i in filter:
        b = fixPaths__([i])[0]
        val = filter[i]
        filter2[b] = val
    ret = Converter.converter.convertFile2PartialPyTree(fileName, format, None, com, filter2, 0)
    return ret

# Ecrit des tableaux ou des morceaux de tableau a certains endroits du fichier
# definis par filter
# t: pyTree avec les memes chemins
def writePyTreeFromFilter(t, fileName, filter, format='bin_hdf', com=None, skelData=None):
    """Write nodes to file given a filter."""
    filter2 = {}
    for i in filter:
        b = fixPaths__([i])[0]
        val = filter[i]
        filter2[b] = val
    # serialise (pas de com pour forcer le hdf seq)
    if com is None:
        Cmpi.seq(Converter.converter.convertPyTree2FilePartial, t, fileName, format, skelData, None, filter)
    # ok si parallel HDF (passer le com de mpi4py)
    else:
        Converter.converter.convertPyTree2FilePartial(t, fileName, format, skelData, com, filter)
    return None

#============================================================================
# Lecture des noms Base/Zones + dims
# Retourne un squelette (depth=3) + la liste des zones path names (znp)
# si baseNames: retourne les znp pour les bases donnees
# si familyZoneNames: retourne les znp des zones de famille donnee
# si BCType: retourne les znp des zones contenant une BC de type BCType
# ou de famille 'familySpecified:WALL'
# maxDepth: peut-etre mis a 2 ou 3 suivant l'utilisant que l'on veut faire du squelette
# readProcNode: si True, ajoute le procNode (pas lu si depth < 4)
#============================================================================
def readZoneHeaders(fileName, format=None, baseNames=None, familyZoneNames=None, BCType=None,
                    maxDepth=3, readProcNode=False, readGridElementRange=False):
    """Read zone headers."""
    a = convertFile2SkeletonTree(fileName, format, maxDepth=maxDepth, maxFloatSize=6)
    # filter by base names
    if baseNames is not None:
        if not isinstance(baseNames, list): baseNames = [baseNames]
        znp = []
        for b in baseNames:
            l = Internal.getNodeFromName1(a, b)
            if l is not None:
                zs = Internal.getZones(l)
                for z in zs: znp.append('/'+b+'/'+z[0])
    # filter by zone family names
    elif familyZoneNames is not None:
        if not isinstance(familyZoneNames, list): familyZoneNames = [familyZoneNames]
        znp = []
        bases = Internal.getBases(a)
        for b in bases:
            zones = Internal.getZones(b)
            for z in zones:
                f = Internal.getNodeFromType1(z, 'FamilyName_t')
                if f is not None:
                    for fz in familyZoneNames:
                        if Internal.getValue(f) == fz:
                            znp.append('/'+b[0]+'/'+z[0])
    elif BCType is not None:
        families = PyTree.getFamilyBCNamesOfType(a, BCType)
        s = BCType.split(':',1)
        if len(s) == 2: families.append(s[1])
        if BCType == 'BCWall':
            families1 = PyTree.getFamilyBCNamesOfType(a, 'BCWallInviscid')
            families2 = PyTree.getFamilyBCNamesOfType(a, 'BCWallViscous*')
            families += families1
            families += families2
        znp = []
        bases = Internal.getBases(a)
        for b in bases:
            zones = Internal.getZones(b)
            for z in zones:
                path = '/'+b[0]+'/'+z[0]
                _readPyTreeFromPaths(a, fileName, path+'/ZoneBC', format)
                nodes = Internal.getNodesFromValue(z, BCType)
                nodes += PyTree.getFamilyBCs(z, families)
                #_convert2SkeletonTree(z)
                if nodes != []: znp.append(path)
    else:
        znp = Internal.getZonePaths(a, pyCGNSLike=True)
    if readProcNode:
        if maxDepth < 2: raise ValueError('loadSkeleton: maxDepth must be >= 2 for use with readProcNode.')
        paths = []
        bases = Internal.getBases(a)
        for b in bases:
            zones = Internal.getZones(b)
            for z in zones:
                path = '/'+b[0]+'/'+z[0]+'/.Solver#Param/proc'
                paths.append(path)
            for z in zones:
                Internal._createUniqueChild(z, '.Solver#Param', 'UserDefinedData_t')
        _readPyTreeFromPaths(a, fileName, paths, format)
    if readGridElementRange:
        paths = []
        bases = Internal.getBases(a)
        for b in bases:
            zones = Internal.getZones(b)
            for z in zones:
                grelt = Internal.getNodesFromType1(z, 'Elements_t')
                for e in grelt:
                    path = '/'+b[0]+'/'+z[0]+'/'+e[0]+'/ElementRange'
                    paths.append(path)
        _readPyTreeFromPaths(a, fileName, paths, format)
    return a, znp

#========================================================================
# Load par containeurs
# Load les containeurs "cont" dans a pour les znp donnees
# cont='GridCoordinates', 'FlowSolution'
#========================================================================
def _loadContainer(a, fileName, znp, cont, format=None):
    if isinstance(cont, list): conts = cont
    else: conts = [cont]
    if isinstance(znp, list): znps = znp
    else: znps = [znp]
    for p in znps:
        paths = []
        for c in conts: paths.append('%s/%s'%(p,c))
        _readPyTreeFromPaths(a, fileName, paths, format)
    return None

# variablesN = ['GridCoordinates/CoordinateX',...]
# variablesC = ['FlowSolution#Centers/Density',...]
def _loadContainerPartial(a, fileName, znp, variablesN=[], variablesC=[], format=None):
    if isinstance(znp, list): znps = znp
    else: znps = [znp]
    for p in znps:
        f = {}
        spl = p.rsplit('/',1)
        zname = spl[1]
        bname = spl[0].rsplit('/',1)[1]
        # Get loc2Glob
        zone = Internal.getNodeFromPath(a, p)
        src = None
        if zone is not None:
            src, loc2glob = Internal.getLoc2Glob(zone)
        if src is None: continue

        j = loc2glob
        pname = src

        # Variables aux noeuds
        paths = []
        for v in variablesN: paths.append('/%s/%s/%s'%(bname, pname, v))

        dim = Internal.getZoneDim(zone)
        if dim[2] == 1 and dim[3] == 1:
            DataSpaceMMRY = [[0], [1], [j[1]-j[0]+1], [1]]
            DataSpaceFILE = [[j[0]-1], [1], [j[1]-j[0]+1], [1]]
            DataSpaceGLOB = [[0]]
        elif dim[3] == 1:
            DataSpaceMMRY = [[0,0], [1,1], [j[1]-j[0]+1,j[3]-j[2]+1], [1,1]]
            DataSpaceFILE = [[j[0]-1,j[2]-1], [1,1], [j[1]-j[0]+1,j[3]-j[2]+1], [1,1]]
            DataSpaceGLOB = [[0]]
        else:
            DataSpaceMMRY = [[0,0,0], [1,1,1], [j[1]-j[0]+1,j[3]-j[2]+1,j[5]-j[4]+1], [1,1,1]]
            DataSpaceFILE = [[j[0]-1,j[2]-1,j[4]-1], [1,1,1], [j[1]-j[0]+1,j[3]-j[2]+1,j[5]-j[4]+1], [1,1,1]]
            DataSpaceGLOB = [[0]]

        for pp in paths: f[pp] = DataSpaceMMRY+DataSpaceFILE+DataSpaceGLOB

        # Variables aux centres
        paths = []
        for v in variablesC: paths.append('/%s/%s/%s'%(bname, pname, v))

        if dim[2] == 1 and dim[3] == 1:
            DataSpaceMMRYC = [[0], [1], [max(j[1]-j[0],1)], [1]]
            DataSpaceFILEC = [[j[0]-1,j[2]-1], [1], [max(j[1]-j[0],1)], [1]]
            DataSpaceGLOBC = [[0]]
        elif dim[3] == 1:
            DataSpaceMMRYC = [[0,0], [1,1], [max(j[1]-j[0],1),max(j[3]-j[2],1)], [1,1]]
            DataSpaceFILEC = [[j[0]-1,j[2]-1], [1,1], [max(j[1]-j[0],1),max(j[3]-j[2],1)], [1,1]]
            DataSpaceGLOBC = [[0]]
        else:
            DataSpaceMMRYC = [[0,0,0], [1,1,1], [max(j[1]-j[0],1),max(j[3]-j[2],1),max(j[5]-j[4],1)], [1,1,1]]
            DataSpaceFILEC = [[j[0]-1,j[2]-1,j[4]-1], [1,1,1], [max(j[1]-j[0],1),max(j[3]-j[2],1),max(j[5]-j[4],1)], [1,1,1]]
            DataSpaceGLOBC = [[0]]

        for pp in paths: f[pp] = DataSpaceMMRYC+DataSpaceFILEC+DataSpaceGLOBC

        r = readNodesFromFilter(fileName, f, format)

        # Repositionne les chemins dans la zone
        for k in r:
            k2 = k.replace(pname, zname)
            spl = k2.rsplit('/',1)
            varName = spl[1]
            contName = spl[0].rsplit('/',1)[1]
            n = Internal.getNodeFromPath(a, k2)
            if n is None:
                parent = k2.rsplit('/',1)[0]
                parent = Internal.getNodeFromPath(a, parent)
                if parent is not None:
                    Internal.createChild(parent, varName, 'DataArray_t', value=r[k])
                    loc = Internal.getNodeFromType1(parent, 'GridLocation_t')
                    if loc is None and contName+'/'+varName in variablesC:
                        Internal.newGridLocation(value='CellCenter', parent=parent)
                else: print('Cannot set %s.'%k2)
            else: n[1] = r[k]
    return None

#=========================================================================
# Load connectivity
# Load les connectivites dans a pour les znps donnes
#=========================================================================
def _loadConnectivity(a, fileName, znp, format=None):
    if isinstance(znp, list): znps = znp
    else: znps = [znp]
    for p in znps:
        z = Internal.getNodeFromPath(a, p)
        elts = Internal.getNodesFromType1(z, 'Elements_t')
        paths = []
        for e in elts: paths.append(p+'/'+e[0])
        _readPyTreeFromPaths(a, fileName, paths, format)
    return None

# force le proc node des zones au processeur courant
def _enforceProcNode(a):
    zones = Internal.getZones(a)
    for z in zones:
        p = Internal.createUniqueChild(z, '.Solver#Param', 'UserDefinedData_t')
        Internal.createUniqueChild(p, 'proc', 'DataArray_t', value=Cmpi.rank)
    return None

#==========================================================================
# Load par variables
# Load les variables "var" pour les znp donnes
# var='Density', 'centers:Density', 'FlowSolution/Density' or list of vars
#==========================================================================
def _loadVariables(a, fileName, znp, var, format):
    if isinstance(var, list): vars = var
    else: vars = [var]
    if isinstance(znp, list): znps = znp
    else: znps = [znp]
    fvars = []; cont = []
    for v in vars:
        s = v.split(':',1)
        if v[0:10] == 'Coordinate':
            fvars.append(Internal.__GridCoordinates__+'/'+v); cont = Internal.__GridCoordinates__
        elif len(s) == 2 and s[0] == 'centers':
            fvars.append(Internal.__FlowSolutionCenters__+'/'+s[1]); cont = Internal.__FlowSolutionCenters__
        elif len(s) == 2 and s[0] == 'nodes':
            fvars.append(Internal.__FlowSolutionNodes__+'/'+s[1]); cont = Internal.__FlowSolutionNodes__
        else:
            s = v.split('/')
            if len(s) == 2:
                fvars.append(s[0]+'/'+s[1]); cont = s[0]
            else:
                fvars.append(Internal.__FlowSolutionNodes__+'/'+v); cont = Internal.__FlowSolutionNodes__

    for p in znps:
        paths = []
        for v in fvars: paths.append('%s/%s'%(p,v))
        _readPyTreeFromPaths(a, fileName, paths, format)
        # Ensure location in containers
        zp = Internal.getNodeFromPath(a, p)
        fp = Internal.getNodeFromPath(a, paths[0])
        Compressor._uncompressAll(zp)
        if zp is not None and fp is not None:
            c = Internal.getNodeFromName1(zp, cont)
            if c is not None and fp[1] is not None:
                if PyTree.getNPts(zp) != fp[1].size:
                    Internal._createUniqueChild(c, 'GridLocation', 'GridLocation_t', 'CellCenter')
    return None

#==================================================================================
# Fully load all path in a
#==================================================================================
def _loadZones(a, fileName, znp, format=None):
    if isinstance(znp, list): znps = znp
    else: znps = [znp]
    _readPyTreeFromPaths(a, fileName, znps, format)
    PyTree._upgradeZone(a, uncompress=True, upgradeNGon=True)


# Fully load zoneBC_t and GridConnectivity_t of znp
def _loadZoneBCs(a, fileName, znp, format=None):
    if isinstance(znp, list): znps = znp
    else: znps = [znp]
    paths = []
    for p in znps:
        n = Internal.getNodeFromPath(a, p)
        if n is not None:
            children = n[2]
            for i in children:
                if i[3] == 'ZoneBC_t': paths.append(p+'/'+i[0])
                elif i[3] == 'ZoneGridConnectivity_t': paths.append(p+'/'+i[0])
                if i[0] == '.Solver#define': paths.append(p+'/'+i[0])
                elif i[0] == '.Solver#ownData': paths.append(p+'/'+i[0])
    if paths != []: _readPyTreeFromPaths(a, fileName, paths, format)
    return None

# Load zone BC but without loading BCDataSet fields
def _loadZoneBCsWoData(a, fileName, znp, format=None):
    if isinstance(znp, list): znps = znp
    else: znps = [znp]
    paths = []
    for p in znps:
        n = Internal.getNodeFromPath(a, p)
        if n is not None:
            children = n[2]
            for i in children:
                if i[3] == 'GridConnectivity_t': paths.append(p+'/'+i[0])
    if paths != []: _readPyTreeFromPaths(a, fileName, paths, format)
    paths = []
    for p in znps:
        n = Internal.getNodeFromPath(a, p)
        if n is not None:
            children = n[2]
            for i in children:
                if i[3] == 'ZoneBC_t': paths.append(p+'/'+i[0])
    if paths != []: _readPyTreeFromPaths(a, fileName, paths, format, maxDepth=2, setOnlyValue=False)
    return None

# Load fully extra nodes of a tree (but no Zone or Base)
def _loadTreeExtras(a, fileName, format=None):
    """Load extra data in tree."""
    # Fully load nodes of level 1 except CGNSBase_t
    paths = []
    children = a[2]
    for i in children:
        if i[3] != 'CGNSBase_t': paths.append(i[0])

    # Fully load nodes of level 2 except Zone_t
    bases = Internal.getBases(a)
    for b in bases:
        children = b[2]
        for i in children:
            if i[3] != 'Zone_t': paths.append(b[0]+'/'+i[0])
    if paths != []: _readPyTreeFromPaths(a, fileName, paths, format)
    return None

# Load extra data of zones (cartesianData or solver#define)
# decompresse eventuellement
def _loadZoneExtras(a, fileName, znp, format=None):
    """Load extra data in zones."""
    if isinstance(znp, list): znps = znp
    else: znps = [znp]
    paths = []
    for p in znps:
        n = Internal.getNodeFromPath(a, p)

        # Level1
        # - DataArray_t)
        for i in n[2]:
            if i[3] == 'DataArray_t' and i[1] is None: paths.append(p+'/'+i[0])

        # Level2
        # UserDefinedData_t/DataArray_t
        for i in n[2]:
            if i[3] == 'UserDefinedData_t':
                for j in i[2]:
                    if j[3] == 'DataArray_t' and j[1] is None: paths.append(p+'/'+i[0]+'/'+j[0])

        # generique niveau 2
        #for i in n[2]:
        #  for j in i[2]:
        #    if j[3] == 'DataArray_t' and j[1] is None: paths.append(p+'/'+i[0]+'/'+j[0])

        # Level3
        # - FlowSolution_t/UserDefinedData_t/DataArray_t
        # - TimeMotion/TimeRigidMotion_t/DataArray_t
        for i in n[2]:
            if i[0] == 'TimeMotion':
                for j in i[2]:
                    if j[3] == 'TimeRigidMotion_t':
                        for k in j[2]:
                            if k[3] == 'DataArray_t' and k[1] is None: paths.append(p+'/'+i[0]+'/'+j[0]+'/'+k[0])

            if i[3] == 'FlowSolution_t':
                for j in i[2]:
                    if j[3] == 'UserDefinedData_t':
                        for k in j[2]:
                            if k[3] == 'DataArray_t' and k[1] is None: paths.append(p+'/'+i[0]+'/'+j[0]+'/'+k[0])

        # generique niveau 3
        #for i in n[2]:
        #  for j in i[2]:
        #    for k in j[2]:
        #      if k[3] == 'DataArray_t' and k[1] is None: paths.append(p+'/'+i[0]+'/'+j[0]+'/'+k[0])

    if paths != []: _readPyTreeFromPaths(a, fileName, paths, format)

    return None

# get variables: return a list
# this doesnt mean that all vars are defined in all zones!
def getVariables(fileName, znp, cont=None, format=None):
    """Return a list of variables contained in file."""
    if isinstance(znp, list): znps = znp
    else: znps = [znp]
    if cont is None: # check containers in files
        conts = set()
        nodes = readNodesFromPaths(fileName, znp, format, maxDepth=1, maxFloatSize=0)
        for n in nodes:
            for j in n[2]:
                if j[3] == 'FlowSolution_t': conts.add(j[0])
        cont = list(conts)
    elif cont == 'centers': cont = [Internal.__FlowSolutionCenters__]
    elif cont == 'nodes': cont = [Internal.__FlowSolutionNodes__]
    elif isinstance(cont, str): cont = [cont]
    paths = []
    for c in cont:
        for p in znps:
            paths.append(p+'/'+c)
    zvars = set()
    nodes = readNodesFromPaths(fileName, paths, format, maxDepth=1, maxFloatSize=0)
    for n in nodes:
        for j in n[2]:
            if j[3] == 'DataArray_t': zvars.add(n[0]+'/'+j[0])
    return list(zvars)

# Return the list of variables in BCDataSet
def getBCVariables(a, fileName, znp, cont=None, format=None):
    """Return a list of variables contained in BCDataSet."""
    if isinstance(znp, list): znps = znp
    else: znps = [znp]
    paths = []
    for p in znps:
        zp = Internal.getNodeFromPath(a, p)
        if zp is not None:
            zoneBC = Internal.getNodesFromType1(zp, 'ZoneBC_t')
            for zbc in zoneBC:
                paths.append(p+'/'+zbc[0])
    nodes = readNodesFromPaths(fileName, paths, format, maxDepth=-1, maxFloatSize=0)
    zvars = set()
    for n in nodes:
        p = Internal.getNodesFromType(n, 'BCData_t')
        for j in p:
            for k in j[2]:
                if k[3] == 'DataArray_t': zvars.add(k[0])
    return list(zvars)

# Load only zones that match a given bbox
def isInBBox(a, fileName, format, bbox, znp):
    """Load zones that lie in bbox."""
    xmin = bbox[0]; ymin = bbox[1]; zmin = bbox[2]
    xmax = bbox[3]; ymax = bbox[4]; zmax = bbox[5]
    if isinstance(znp, list): znps = znp
    else: znps = [znp]
    out = []
    for p in znps:
        z = Internal.getNodeFromPath(a, p)
        path = ['%s/GridCoordinates/CoordinateX'%p]
        _readPyTreeFromPaths(a, fileName, path, format)
        pt = Internal.getNodeFromName2(z, 'CoordinateX')
        vmin = numpy.min(pt[1])
        vmax = numpy.max(pt[1])
        pt[1] = None
        if vmax < xmin or vmin > xmax: out.append(False); continue
        path = ['%s/GridCoordinates/CoordinateY'%p]
        _readPyTreeFromPaths(a, fileName, path, format)
        pt = Internal.getNodeFromName2(z, 'CoordinateY')
        vmin = numpy.min(pt[1])
        vmax = numpy.max(pt[1])
        pt[1] = None
        if vmax < ymin or vmin > ymax: out.append(False); continue
        path = ['%s/GridCoordinates/CoordinateZ'%p]
        _readPyTreeFromPaths(a, fileName, path, format)
        pt = Internal.getNodeFromName2(z, 'CoordinateZ')
        vmin = numpy.min(pt[1])
        vmax = numpy.max(pt[1])
        pt[1] = None
        if vmax < zmin or vmin > zmax: out.append(False); continue
        out.append(True)
    if len(out) == 1: out = out[0]
    return out

#==========================================================
# a: must be a tree or a zone list coherent with znp
# znp: is the full path from top
#==========================================================
def writeZones(a, fileName, znp, format=None):
    """Write Zones in file."""
    if isinstance(znp, list): znps = znp
    else: znps = [znp]
    if len(a) < 4: zones = a # suppose single zone in a list
    else:
        if a[3] == 'CGNSTree_t':
            zones = []
            for p in znps: zones.append(Internal.getNodeFromPath(a, p))
        elif a[3] == 'Zone_t': zones = [a]
        else: zones = a

    if len(zones) != len(znps): raise ValueError('writeZones: data and znp are incompatible.')

    paths = []
    for p in znps:
        pp = p.rsplit('/', 1)
        paths.append(pp[0])
    writeNodesFromPaths(fileName, paths, zones, format, mode=0)

# Write all except flowfield containers
def writeZonesWoVars(a, fileName, znp, format=None):
    """Write zones without variables."""
    if isinstance(znp, list): znps = znp
    else: znps = [znp]
    if len(a) < 4: zones = a # suppose single zone in a list
    else:
        if a[3] == 'CGNSTree_t':
            zones = []
            for p in znps: zones.append(Internal.getNodeFromPath(a, p))
        elif a[3] == 'Zone_t': zones = [a]
        else: zones = a

    if len(zones) != len(znps): raise ValueError('writeZones: data and znp are incompatible.')

    ppaths = []
    for p in znps:
        pp = p.rsplit('/', 1)
        ppaths.append(pp[0])
    writeNodesFromPaths(fileName, ppaths, zones, format, maxDepth=0, mode=0)

    nodes = []; paths = []
    # skipType=FlowSolution_t
    for c, z in enumerate(zones):
        children = z[2]
        for i in children:
            if i[3] != 'FlowSolution_t':
                nodes.append(i)
                paths.append(znps[c])
    writeNodesFromPaths(fileName, paths, nodes, format, mode=0)

# Write zone variables only
def writeVariables(a, fileName, var, znp, format=None):
    """Write variables in file."""
    if isinstance(znp, list): znps = znp
    else: znps = [znp]
    if isinstance(var, list): vars = var
    else: vars = [var]
    if len(a) < 4: zones = a # suppose single zone in a list
    else:
        if a[3] == 'CGNSTree_t':
            zones = []
            for p in znps: zones.append(Internal.getNodeFromPath(a, p))
        elif a[3] == 'Zone_t': zones = [a]
        else: zones = a

    if len(zones) != len(znps): raise ValueError('writeZones: data and znp are incompatible.')

    loc = []
    for v in vars:
        vs = v.split(':',1)
        if len(vs) == 2:
            if vs[0] == 'centers': loc.append(Internal.__FlowSolutionCenters__)
            else: loc.append(Internal.__FlowSolutionNodes__)
        else: loc.append(Internal.__FlowSolutionNodes__)

    conts = []; cpaths = []; nodes = []; npaths = []
    for c, z in enumerate(zones):
        for i in z[2]:
            if i[3] == 'FlowSolution_t':
                conts.append(i); cpaths.append(znps[c])
            for j in i[2]:
                if j[3] == 'GridLocation_t':
                    conts.append(j); cpaths.append(znps[c]+'/'+i[0])
        for d, v in enumerate(vars):
            ns = PyTree.getStdNodesFromName(z, v)
            if ns != []:
                npaths.append(znps[c]+'/'+loc[d])
                nodes.append(ns[0])

    writeNodesFromPaths(fileName, cpaths, conts, format, maxDepth=0, mode=0)
    writeNodesFromPaths(fileName, npaths, nodes, format, mode=0)

#==========================================================
class Handle:
    """Handle for partial reading."""
    def __init__(self, fileName):
        self.fileName = fileName
        self.fileVars = None # vars existing in file
        self.fileBCVars = None # BCDataSet vars existing in file
        self.varsN = [] # Variable in nodes in file
        self.varsC = [] # Variable in centers in file
        self.znp = None # zone paths existing in file
        self.size = None # size des zones du fichier
        self.bary = None # Barycentres zones du fichier
        self.bbox = None # BBox zones du fichier
        self.hmoy = None # pas moyen des zones du fichier
        self.format = Converter.checkFileType(fileName) # Real format of file

    # Retourne les chemins des zones de a
    def getZonePaths(self, a):
        """Return the path of zones in a."""
        out = []
        bases = Internal.getBases(a)
        if len(bases) > 0:
            for b in bases:
                zones = Internal.getZones(b)
                for z in zones: out.append('/'+b[0]+'/'+z[0])
        else:
            zones = Internal.getZones(a)
            for z in zones:
                for p in self.znp:
                    r = p.rsplit('/',1)[1]
                    if r == z[0]: out.append(p); break
        return out

    def setZnp(self, a):
        self.znp = []
        bases = Internal.getBases(a)
        for b in bases:
            zones = Internal.getZones(b)
            for z in zones: self.znp.append('/'+b[0]+'/'+z[0])

    # Retourne les variables du fichier
    def getVariables(self, a=None, cont=None):
        """Return the variable names contained in file."""
        if a is not None: p = self.getZonePaths(a)
        else: p = self.znp # all zone vars
        vars = getVariables(self.fileName, p, cont, self.format)
        self.fileVars = vars
        return vars

    # Retourne les variables de BCDatSet du fichier
    def getBCVariables(self, a):
        """Return the variables contained in BCDataSet."""
        if a is not None: p = self.getZonePaths(a)
        else: p = self.znp
        vars = getBCVariables(a, self.fileName, p, None, self.format)
        self.fileBCVars = vars
        return vars

    # Charge un squelette, stocke les chemins des zones du fichier (znp)
    # Stocke le nombre de pts de chaque zone
    def loadSkeleton(self, maxDepth=3, readProcNode=False):
        """Load a skeleton tree."""
        if Cmpi.rank == 0:
            a, self.znp = readZoneHeaders(self.fileName, self.format, maxDepth=maxDepth, readProcNode=readProcNode)
            # evaluation taille des zones
            self.size = {}
            for zn in self.znp:
                z = Internal.getNodeFromPath(a, zn)
                if z[1] is not None:
                    pt = z[1].ravel('k')
                    if len(pt) == 6: s = pt[0]*pt[1]*pt[2]
                    else: s = pt[0]
                else: s = 0
                self.size[zn] = s
        else: a = None
        a = Cmpi.bcast(a)
        return a

    # Charge le squelette, le split et conserve les infos de split
    def loadAndSplitSkeleton(self, NParts=None, NProc=Cmpi.size, splitByBase=False, algorithm='graph'):
        """Load and split skeleton."""
        if Cmpi.rank == 0:
            a, self.znp = readZoneHeaders(self.fileName, self.format)
            # force loc2glob to himself from himself
            for b in Internal.getBases(a):
                for z in Internal.getZones(b):
                    dim = Internal.getZoneDim(z)
                    if dim[0] == 'Structured':
                        Internal._setLoc2Glob(z, z[0], win=[1,dim[1],1,dim[2],1,dim[3]], sourceDim=dim[1:])

            # Lecture ZoneBC + ZoneGC necessaire pour le split
            for b in Internal.getBases(a):
                for z in Internal.getZones(b):
                    paths = []
                    if Internal.getNodeFromName1(z,"ZoneBC") is not None:
                        paths.append(b[0]+'/'+z[0]+'/ZoneBC')
                    if Internal.getNodeFromName1(z,"ZoneGridConnectivity") is not None:
                        paths.append(b[0]+'/'+z[0]+'/ZoneGridConnectivity')
                    if paths != []:
                        _readPyTreeFromPaths(a, self.fileName, paths)

            # Lecture des noms de variables
            varsN = ['%s/CoordinateX'%Internal.__GridCoordinates__,
                     '%s/CoordinateY'%Internal.__GridCoordinates__,
                     '%s/CoordinateZ'%Internal.__GridCoordinates__]
            varsC = []
            for b in Internal.getBases(a):
                for z in Internal.getZones(b):
                    if Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__) is not None:
                        node = readNodesFromPaths(self.fileName, b[0]+'/'+z[0]+'/'+Internal.__FlowSolutionNodes__, maxDepth=1, maxFloatSize=0)
                        for n in node[2]:
                            if n[3] == 'DataArray_t': varsN.append(Internal.__FlowSolutionNodes__+'/'+n[0])
                        break
            for b in Internal.getBases(a):
                for z in Internal.getZones(b):
                    if Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__) is not None:
                        node = readNodesFromPaths(self.fileName, b[0]+'/'+z[0]+'/'+Internal.__FlowSolutionCenters__, maxDepth=1, maxFloatSize=0)
                        for n in node[2]:
                            if n[3] == 'DataArray_t': varsC.append(Internal.__FlowSolutionCenters__+'/'+n[0])
                        break

            import Transform.PyTree as T

            if splitByBase: # split on skeleton par base
                bases = Internal.getBases(a)
                for b in bases:
                    if NParts is not None: T._splitNParts(b, N=NParts)
                    else: T._splitNParts(b, N=NProc)
                    if NProc is not None:
                        import Distributor2.PyTree as D2
                        D2._distribute(b, NProc, algorithm=algorithm)
            else: # split on full skeleton
                if NParts is not None: T._splitNParts(a, N=NParts)
                else: T._splitNParts(a, N=NProc)
                if NProc is not None:
                    import Distributor2.PyTree as D2
                    D2._distribute(a, NProc, algorithm=algorithm)
        else:
            a = None; varsN = None; varsC = None
        a = Cmpi.bcast(a)
        self.varsN = Cmpi.bcast(varsN); self.varsC = Cmpi.bcast(varsC)

        # Force GridLocation in FlowSolutionCenters
        for z in Internal.getZones(a):
            conts = Internal.getNodesFromName1(z, Internal.__FlowSolutionCenters__)
            for c in conts:
                Internal._createUniqueChild(c, 'GridLocation', 'GridLocation_t', 'CellCenter')

        return a

    def loadAndSplit(self, NParts=None, NProc=Cmpi.size, splitByBase=False, algorithm='graph'):
        """Load and split a file."""
        a = self.loadAndSplitSkeleton(NParts, NProc, splitByBase, algorithm)
        _convert2PartialTree(a, rank=Cmpi.rank)
        self._loadContainerPartial(a, variablesN=self.varsN, variablesC=self.varsC)
        return a

    def loadFromProc(self, loadVariables=True, rank=Cmpi.rank):
        """Load and distribute zones from proc node."""
        if Cmpi.rank == 0:
            # Load le squelette niveau2 + les noeuds proc
            t = convertFile2SkeletonTree(self.fileName, self.format, maxDepth=2, maxFloatSize=6)
            zones = Internal.getZones(t)
            paths = []
            for z in zones:
                p = Internal.getPath(t, z)+'/.Solver#Param/proc'
                paths.append(p)
            for z in zones:
                Internal._createUniqueChild(z, '.Solver#Param', 'UserDefinedData_t')
            _readPyTreeFromPaths(t, self.fileName, paths)
            self._loadTreeExtras(t)
        else: t = None
        t = Cmpi.bcast(t)
        # Lit les zones correspondant a proc
        paths = []
        zones = Internal.getZones(t)
        for z in zones:
            proc = Internal.getNodeFromName2(z, 'proc')

            if Internal.getValue(proc) == rank:
                paths.append(Internal.getPath(t, z))
            else: Internal._rmNode(t, z)

        #if loadVariables: skipTypes=None
        #else: skipTypes=['FlowSolution_t']
        if paths != []: _readPyTreeFromPaths(t, self.fileName, paths)
        for z in Internal.getZones(t): PyTree._upgradeZone(z, uncompress=True, upgradeNGon=True)
        return t

    # strategy=strategie pour la distribution (match)
    # algorithm=type d'algorithme pour la distribution
    # loadVariables=True, charge toutes les variables sinon ne charge que les coords
    def loadAndDistribute(self, strategy=None, algorithm='graph', loadVariables=True):
        """Load and distribute zones."""
        if Cmpi.rank == 0:
            if strategy == 'match':
                # Lit le squelette niveau 3 + les zoneGridConnectivity
                t = convertFile2SkeletonTree(self.fileName, self.format, maxDepth=3, maxFloatSize=6)
                paths = []
                bases = Internal.getBases(t)
                for b in bases:
                    zones = Internal.getZones(b)
                    for z in zones:
                        pz = '%s/%s'%(b[0],z[0])
                        gcs = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
                        for gc in gcs:
                            p = pz+'/%s'%gc[0]
                            paths.append(p)
                if paths != []: _readPyTreeFromPaths(t, self.fileName, paths, self.format)
            else:
                # Lit le squelette niveau 2 + les noeuds procs
                t = convertFile2SkeletonTree(self.fileName, self.format, maxDepth=2, maxFloatSize=6)
                paths = []
                bases = Internal.getBases(t)
                for b in bases:
                    zones = Internal.getZones(b)
                    for z in zones:
                        p = '%s/%s/ZoneType'%(b[0],z[0])
                        paths.append(p)
                _readPyTreeFromPaths(t, self.fileName, paths, self.format)
            # Load les extras de l'arbre (autre que base)
            self._loadTreeExtras(t)
            # Distribue
            import Distributor2.PyTree as D2
            D2._distribute(t, Cmpi.size, useCom=strategy, algorithm=algorithm)
        else: t = None
        t = Cmpi.bcast(t)
        # Lit les zones affectees
        paths = []
        bases = Internal.getBases(t)
        for b in bases:
            zones = Internal.getZones(b)
            for z in zones:
                proc = Internal.getNodeFromName2(z, 'proc')
                if Internal.getValue(proc) == Cmpi.rank: paths.append("%s/%s"%(b[0],z[0]))
                else: Internal._rmNodeByPath(t, "%s/%s"%(b[0],z[0]))
        if loadVariables: skipTypes=None
        else: skipTypes=['FlowSolution_t']
        if paths != []: _readPyTreeFromPaths(t, self.fileName, paths, self.format, skipTypes=skipTypes)
        _enforceProcNode(t)
        # Decompression eventuelle
        for z in Internal.getZones(t): PyTree._upgradeZone(z, uncompress=True, upgradeNGon=True)
        return t

    def distributedLoadAndSplitSkeleton(self, NParts=None, NProc=Cmpi.size):
        """
        Load mesh topology in a distributed way: each processus loads a partial data of the
        connectivity, to have a partitionned data among the processes.

        :param      NParts:  The number of parts
        :type       NParts:  Integer
        :param      NProc:   The number of processes
        :type       NProc:   Integer. Value by default : number  of processes runned with MPI
        """
        from XCore import xcore
        a = None
        if Cmpi.rank == 0:
            a, self.znp = readZoneHeaders(self.fileName, self.format, readGridElementRange=True)
        a = Cmpi.bcast(a)
        # force loc2glob to himself from himself
        bases = Internal.getBases(a)
        nb_tot_verts = 0
        for z in Internal.getZones(a):
            dim = Internal.getZoneDim(z)
            nb_tot_verts += dim[1]
        f = {}
        for b in bases:
            for z in Internal.getZones(b):
                connectivities = Internal.getElementNodes(z)
                for c in connectivities:
                    eltName, nb_nodes_per_elt = Internal.eltNo2EltName(c[1][0])
                    rg = Internal.getNodeFromName1(c, 'ElementRange')
                    nb_elts = rg[1][1]
                    nb_loc_elts = nb_elts // Cmpi.size
                    reste_elts = nb_elts % Cmpi.size
                    if Cmpi.rank < reste_elts:
                        nb_loc_elts += 1
                    beg_elt = Cmpi.rank * nb_loc_elts
                    if Cmpi.rank >= reste_elts:
                        beg_elt += reste_elts
                    beg_elt *= nb_nodes_per_elt
                    nb_loc_elts *= nb_nodes_per_elt
                    f[b[0]+'/'+z[0]+'/' + c[0] + '/ElementConnectivity'] = [[0], [1], [nb_loc_elts], [1],
                                                                            [beg_elt], [1], [nb_loc_elts], [1],
                                                                            [nb_elts*nb_nodes_per_elt]]
                dim = Internal.getZoneDim(z)
                nb_verts = dim[1]
                nb_loc_verts = nb_verts//Cmpi.size
                reste_verts = nb_verts % Cmpi.size
                if Cmpi.rank < reste_verts:
                    nb_loc_verts += 1
                beg_verts = nb_loc_verts * Cmpi.rank
                if Cmpi.rank >= reste_verts:
                    beg_verts += reste_verts
                f[b[0]+'/'+z[0]+'/GridCoordinates/CoordinateX'] = [[0], [1], [nb_loc_verts], [1],
                                                                   [beg_verts],[1],[nb_loc_verts], [1], [nb_verts]]
                f[b[0]+'/'+z[0]+'/GridCoordinates/CoordinateY'] = [[0], [1], [nb_loc_verts], [1],
                                                                   [beg_verts],[1],[nb_loc_verts], [1], [nb_verts]]
                f[b[0]+'/'+z[0]+'/GridCoordinates/CoordinateZ'] = [[0], [1], [nb_loc_verts], [1],
                                                                   [beg_verts],[1],[nb_loc_verts], [1], [nb_verts]]
        dvars = readNodesFromFilter(self.fileName, f)
        # Mis en donne par zone pour le decoupeur:
        zones = []
        for b in bases:
            for z in Internal.getZones(b):
                dim = Internal.getZoneDim(z)
                nb_verts     = dim[1]
                nb_glob_elts = dim[2]
                coordinates = (dvars[b[0]+'/'+z[0]+'/GridCoordinates/CoordinateX'],
                               dvars[b[0]+'/'+z[0]+'/GridCoordinates/CoordinateY'],
                               dvars[b[0]+'/'+z[0]+'/GridCoordinates/CoordinateZ'])
                zone = ([], coordinates, nb_verts, nb_glob_elts)
                for c in connectivities:
                    eltName, nb_nodes_per_elt = Internal.eltNo2EltName(c[1][0])
                    zone[0].append((eltName,dvars[b[0]+'/'+z[0]+'/' + c[0] + '/ElementConnectivity']))
                zones.append(zone)
        splitted_data = xcore.split_elements(zones)

        splitted_a = Internal.copyTree(a)
        splitted_bases = Internal.getBases(splitted_a)
        i_zone = 0
        for b in splitted_bases:
            for z in Internal.getZones(b):
                vertexTag = splitted_data[0]["vertexTag"]
                nbElts  = numpy.sum(splitted_data[0]["cellTag"] == i_zone)
                nbVerts = numpy.sum(vertexTag == i_zone)
                #nbElts = splitted_data[0]["cellLN2GN"].shape[0]
                #nbVerts = splitted_data[0]["vertex"].shape[0]
                z[1][0][0] = nbVerts
                z[1][0][1] = nbElts
                connectivities = Internal.getElementNodes(z)
                for c in connectivities:
                    for type_elt in splitted_data[0]["cellVertex"]:
                        num_elt, nb_vertices_per_elt = Internal.eltName2EltNo(type_elt)
                        zone_tag, cellvertex = splitted_data[0]["cellVertex"][type_elt]
                        #print(f"zone tag : {zone_tag}, cellvertex : {cellvertex}, num_elt: {num_elt}")
                        if c[1][0] == num_elt:
                            #print(f"Rajout des connectivites")
                            connectivity = numpy.array([cellvertex[i] for i in range(0,nb_vertices_per_elt*zone_tag.shape[0]) if vertexTag[cellvertex[i]-1] == i_zone])
                            eltconnectivity = Internal.createNode('ElementConnectivity', 'DataArray_t', value=connectivity)
                            eltrange = Internal.createNode('ElementRange', 'IndexRange_t', value=numpy.array([1, connectivity.shape[0]]))
                            Internal.addChild(c, eltrange)
                            Internal.addChild(c, eltconnectivity)
                grdcrd = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
                xcrds = numpy.array(splitted_data[0]["vertex"][:,0])
                ycrds = numpy.array(splitted_data[0]["vertex"][:,1])
                zcrds = numpy.array(splitted_data[0]["vertex"][:,2])
                xcoords = Internal.createNode('CoordinateX', 'DataArray_t', value=xcrds)
                ycoords = Internal.createNode('CoordinateY', 'DataArray_t', value=ycrds)
                zcoords = Internal.createNode('CoordinateZ', 'DataArray_t', value=zcrds)
                Internal.addChild(grdcrd, xcoords)
                Internal.addChild(grdcrd, ycoords)
                Internal.addChild(grdcrd, zcoords)

        return splitted_a

    # Calcul et stocke des infos geometriques sur les zones
    def geomProp(self, znp=None):
        """Store some zone properties."""
        if znp is None: znp = self.znp

        if self.bbox is None: self.bbox = {}
        if self.bary is None: self.bary = {}
        if self.hmoy is None: self.hmoy = {}
        for p in znp:
            if p in self.bbox: continue
            xc = 0.; yc = 0.; zc = 0.
            bbox = [0.,0.,0.,0.,0.,0.]
            # Coordinate X
            path = '%s/GridCoordinates/CoordinateX'%p
            n = readNodesFromPaths(self.fileName, path, maxDepth=0)
            vmin = numpy.min(n[1])
            vmax = numpy.max(n[1])
            bbox[0] = vmin; bbox[3] = vmax
            xc = 0.5*(vmin+vmax)
            # Coordinate Y
            path = '%s/GridCoordinates/CoordinateY'%p
            n = readNodesFromPaths(self.fileName, path, maxDepth=0)
            vmin = numpy.min(n[1])
            vmax = numpy.max(n[1])
            bbox[1] = vmin; bbox[4] = vmax
            yc = 0.5*(vmin+vmax)
            # Coordinate Z
            path = '%s/GridCoordinates/CoordinateZ'%p
            n = readNodesFromPaths(self.fileName, path, maxDepth=0)
            vmin = numpy.min(n[1])
            vmax = numpy.max(n[1])
            bbox[2] = vmin; bbox[5] = vmax
            zc = 0.5*(vmin+vmax)
            self.bary[p] = (xc,yc,zc)
            self.bbox[p] = bbox
            dhx = bbox[3]-bbox[0]
            dhy = bbox[4]-bbox[1]
            dhz = bbox[5]-bbox[2]
            npts = self.size[p]**0.33
            self.hmoy[p] = (dhx+dhy+dhz)*0.33*npts

    # Charge les Grid coordinates + grid connectivity + BC
    # pour toutes les zones de a
    def _loadGCBC(self, a):
        self._loadContainer(a, 'GridCoordinates')
        self._loadConnectivity(a)
        self._loadZoneBCs(a)
        Internal._fixNGon(a)
        return None

    # Charge tous les noeuds extra d'un arbre a
    def _loadTreeExtras(self, a):
        """Load extra data in tree."""
        _loadTreeExtras(a, self.fileName, self.format)
        return None

    # Charge tous les noeuds extra des zones
    def _loadZoneExtras(self, a):
        """Load extra data in zones."""
        _loadZoneExtras(a, self.fileName, self.format)
        return None

    # Charge completement toutes les zones de a ou de chemin fourni
    def loadZones(self, a, znp=None):
        """Fully load zones."""
        b = Internal.copyRef(a)
        self._loadZones(b, znp)
        return b

    def _loadZones(self, a, znp=None):
        """Fully load zones."""
        if znp is None: znp = self.getZonePaths(a)
        _loadZones(a, self.fileName, znp, self.format)
        return None

    # Charge les grid coordinates + grid connectivity + BC
    # pour les zones specifiees
    def loadZonesWoVars(self, a, znp=None, bbox=None):
        """Load zones without loading variables."""
        b = Internal.copyRef(a)
        self._loadZonesWoVars(b, znp, bbox)
        return b

    def _loadZonesWoVars(self, a, znp=None, bbox=None):
        if znp is None: znp = self.getZonePaths(a)
        if bbox is None:
            # Read paths as skeletons
            _readPyTreeFromPaths(a, self.fileName, znp, self.format, maxFloatSize=0)
            _loadContainer(a, self.fileName, znp, Internal.__GridCoordinates__, self.format)
            _loadConnectivity(a, self.fileName, znp, self.format)
            _loadZoneBCs(a, self.fileName, znp, self.format)
            _loadZoneExtras(a, self.fileName, znp, self.format)
            for zp in znp:
                _convert2PartialTree(Internal.getNodeFromPath(a, zp))
        else:
            if self.bbox is None: self.geomProp()
            for zp in znp:
                zbb = self.bbox[zp]
                if zbb[3] >= bbox[0] and zbb[0] <= bbox[3] and zbb[4] >= bbox[1] and zbb[1] <= bbox[4] and zbb[5] >= bbox[2] and zbb[2] <= bbox[5]:
                    print('loading: %s'%zp)
                    _readPyTreeFromPaths(a, self.fileName, [zp], self.format, maxFloatSize=0)
                    _loadContainer(a, self.fileName, [zp], Internal.__GridCoordinates__, self.format)
                    _loadConnectivity(a, self.fileName, [zp], self.format)
                    _loadZoneBCs(a, self.fileName, [zp], self.format)
                    _loadZoneExtras(a, self.fileName, znp, self.format)
                    _convert2PartialTree(Internal.getNodeFromPath(a, zp))
        PyTree._upgradeZone(a, uncompress=True, upgradeNGon=True)
        return None

    # Charge toutes les BCs (avec BCDataSet) des zones de a
    def _loadZoneBCs(self, a, znp=None):
        if znp is None: znp = self.getZonePaths(a)
        _loadZoneBCs(a, self.fileName, znp, self.format)
        return None

    # Charge toutes les BCs (sans BCDataSet) des zones de a
    def _loadZoneBCsWoData(self, a, znp=None):
        if znp is None: znp = self.getZonePaths(a)
        _loadZoneBCsWoData(a, self.fileName, znp, self.format)
        return None

    # Charge la connectivite pour toutes les zones de a
    def _loadConnectivity(self, a, znp=None):
        if znp is None: znp = self.getZonePaths(a)
        _loadConnectivity(a, self.fileName, znp, self.format)
        return None

    # Charge le container "cont" pour toutes les zones de a
    def _loadContainer(self, a, cont, znp=None):
        if znp is None: znp = self.getZonePaths(a)
        _loadContainer(a, self.fileName, znp, cont, self.format)
        return None

    # Charge le container "cont" en partiel (si la zone a loc2glob)
    def _loadContainerPartial(self, a, variablesN=[], variablesC=[], znp=None):
        if znp is None: znp = self.getZonePaths(a)
        _loadContainerPartial(a, self.fileName, znp, variablesN, variablesC, self.format)
        return None

    def loadVariables(self, a, var, znp=None):
        """Load specified variables."""
        b = Internal.copyRef(a)
        self._loadVariables(b, var, znp)
        return b

    # Charge la ou les variables "var" pour toutes les zones de a
    def _loadVariables(self, a, var, znp=None):
        """Load specified variables."""
        if znp is None: znp = self.getZonePaths(a)
        _loadVariables(a, self.fileName, znp, var, self.format)
        return None

    def isInBBox(self, a, bbox, znp=None):
        """Return true if a is in bbox."""
        if znp is None: znp = self.getZonePaths(a)
        return isInBBox(a, self.fileName, self.format, bbox, znp)

    # Ecrit des zones
    def writeZones(self, a, fileName=None, znp=None):
        """Write specified zones."""
        if znp is None: znp = self.getZonePaths(a)
        if fileName is None: fileName = self.fileName
        writeZones(a, fileName, znp, self.format)

    # Ecrit des zones sans les FlowSolution_t
    def writeZonesWoVars(self, a, fileName=None, znp=None):
        """Write specified zones without variables."""
        if znp is None: znp = self.getZonePaths(a)
        if fileName is None: fileName = self.fileName
        writeZonesWoVars(a, fileName, znp, self.format)

    # Ecrit des variables
    def writeVariables(self, a, var, fileName=None, znp=None):
        """Write specified variables."""
        if znp is None: znp = self.getZonePaths(a)
        if fileName is None: fileName = self.fileName
        writeVariables(a, fileName, var, znp, self.format)

    # save zones + field
    # mode 0: parallele, 1: ecriture chacun son tour
    def save(self, a, fileName=None, cartesian=False):
        """Dave zone and fields in file."""
        a2 = Internal.copyRef(a)
        if cartesian: Compressor._compressCartesian(a2)
        if fileName is not None:
            Cmpi.convertPyTree2File(a2, fileName, self.format, ignoreProcNodes=True)

    def mergeAndSave(self, a, fileName=None):
        """Merge and save in file."""
        if fileName is not None:
            # write data without zones
            ap = Internal.copyRef(a)
            bases = Internal.getBases(ap)
            for b in bases: b[2] = []
            ap = Cmpi.allgatherTree(ap)
            if Cmpi.rank == 0: PyTree.convertPyTree2File(ap, fileName)
            Cmpi.barrier()
            # Ecrit les zones partiellement
            fr = {}
            writePyTreeFromFilter(a, fileName, fr, skelData=[])
