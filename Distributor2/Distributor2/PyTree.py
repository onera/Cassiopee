import Distributor2
import Converter.Internal as Internal
import Converter.PyTree as C
import Converter
import Generator as G
import numpy
__version__ = Distributor2.__version__

try: range = xrange
except: pass

#==============================================================================
# Calcul la liste des bbox
# IN: arrays: liste des zones sous forme array
# IN: zoneNames: nom des zones correspondant
#==============================================================================
def computeBBoxes__(arrays, zoneNames):
    bboxes = []; c = 0
    for a in arrays:
        try: bb = G.bbox(a); bb = bb+[zoneNames[c],True] 
        except: bb = [0,0,0,1,1,1,zoneNames[c],False]
        bboxes.append(bb)
        c += 1

    # Parallel eventuel
    try: 
        import Converter.Mpi as Cmpi
        allboxes = Cmpi.allgather(bboxes)
        c = 0
        for bb in bboxes:
            if bb[7] == False:
                for j in allboxes:
                    for k in j:
                        if k[6] == bb[6] and k[7] == True:
                            bboxes[c] = k
            c += 1                
    except: pass
    return bboxes

#==============================================================================
# Distribute t (pyTree) over NProc processors
# IN: NProc: number of processors
# IN: prescribed: dict containing the zones names as key, and the 
# prescribed proc as value
# IN: perfo: describes performance of processors
# IN: weight: weight assigned to zones of t as a list of integers. Must be ordered as the zones in the pyTree
# IN: useCom: use intergrid connectivity in distribution
# if useCom=0, only the number of points is taken into account (full/skel/load skel)
# if useCom='all', take match, overlap into account (full/load skel)
# if useCom='match', take only match into account (full/skel/load skel)
# if useCom='overlap', take only overlap into account (full/load skel)
# if useCom='bbox', take bbox intersection into account (full/load skel)
# IN: algorithm: gradient0, gradient1, genetic, fast
# IN: nghost: nbre de couches de ghost cells ajoutees
#==============================================================================
def distribute(t, NProc, prescribed={}, perfo=[], weight={}, useCom='all', 
               algorithm='graph', mode='nodes', nghost=0):
    """Distribute a pyTree over processors.
    Usage: distribute(t, NProc, prescribed={}, perfo=[], weight={}, useCom='all', algorithm='graph', mode='nodes', nghost=0)"""
    tp = Internal.copyRef(t)
    out = _distribute(tp, NProc, prescribed=prescribed, perfo=perfo,
                      weight=weight, useCom=useCom, algorithm=algorithm,
                      mode='nodes', nghost=nghost)
    return tp, out

# in place version
def _distribute(t, NProc, prescribed={}, perfo=[], weight={}, useCom='all', 
                algorithm='graph', mode='nodes', nghost=0):
    """Distribute a pyTree over processors.
    Usage: _distribute(t, NProc, prescribed={}, perfo=[], weight={}, useCom='all', algorithm='graph', mode='nodes', nghost=0)"""
    zones = Internal.getZones(t)
    # Formation des arrays
    arrays = []; zoneNames = []; aset = []; weightlist = [] # weight for all zones
    for z in zones:
        zname = z[0]
        zoneNames.append(zname)
        aset.append(prescribed.get(zname,-1))

        if zname in weight: weightlist.append(weight[zname])
        else: weightlist.append(1)

        a = C.getFields(Internal.__GridCoordinates__, z, api=2)
        if a == [[]]: # no coord present in z
            dim = Internal.getZoneDim(z)
            if dim[0] == 'Structured':
                ar = Converter.array('x', dim[1], dim[2], dim[3], api=2)
            elif dim[3] == 'NGON':
                ar = Converter.array('x', dim[1], dim[2], 'HEXA', api=2) # hack
            else:
                ar = Converter.array('x', dim[1], dim[2], dim[3], api=2)
            arrays.append(ar)
        else: arrays.append(a[0])
    Nb = len(arrays)
    com = numpy.zeros((Nb, Nb), numpy.int32)

    if useCom == 'match' or useCom == 'all':
        # Formation des coms - raccords coincidents
        tpp, typen = Internal.node2PyTree(t) 
        bases = Internal.getBases(tpp)
        zc = 0; c = 0
        for b in bases:
            zones = Internal.getNodesFromType1(b, 'Zone_t') 
            dict = {}
            pc = 0
            for z in zones: dict[z[0]] = pc; pc += 1
            
            for z in zones:
                match = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                for m in match:
                    donorName = Internal.getValue(m)
                    if donorName in dict: d = dict[donorName]+zc
                    else: d = -1
                    node = Internal.getNodeFromName1(m, 'PointRange')
                    win = node[1]
                    w = Internal.range2Window(win)
                    vol = (w[1]-w[0]+1)*(w[3]-w[2]+1)*(w[5]-w[4]+1)
                    if d != -1: com[c, d] += vol
                c += 1
            zc += len(zones)

    if useCom == 'overlap' or useCom == 'all':
        # Formation des coms - raccords recouvrants
        tol = 1.e-12
        bboxes = computeBBoxes__(arrays, zoneNames)
        c = 0
        zones = Internal.getZones(t)
        for z in zones:
            m1 = Internal.getNodesFromType2(z, 'GridConnectivity_t')
            for m in m1:
                m2 = Internal.getNodesFromType1(m, 'GridConnectivityType_t')
                if m2 != []:
                    v = Internal.getValue(m2[0])
                    if v == 'Overset':
                        node = Internal.getNodeFromName1(m, 'PointRange')
                        win = node[1]
                        w = Internal.range2Window(win)
                        vol = (w[1]-w[0]+1)*(w[3]-w[2]+1)*(w[5]-w[4]+1)
                        d = 0
                        for z2 in zones:
                            if id(z) != id(z2):
                                bboxc = bboxes[c]; bboxd = bboxes[d]
                                xmin1 = bboxc[0]; xmax1 = bboxc[3]
                                ymin1 = bboxc[1]; ymax1 = bboxc[4]
                                zmin1 = bboxc[2]; zmax1 = bboxc[5]
                                xmin2 = bboxd[0]; xmax2 = bboxd[3]
                                ymin2 = bboxd[1]; ymax2 = bboxd[4]
                                zmin2 = bboxd[2]; zmax2 = bboxd[5]
                                if (xmax1 > xmin2-tol and xmin1 < xmax2+tol and
                                    ymax1 > ymin2-tol and ymin1 < ymax2+tol and
                                    zmax1 > zmin2-tol and zmin1 < zmax2+tol):
                                    com[c, d] += vol
                            d += 1
            c += 1

    if useCom == 'bbox':
        # Formation des coms - si les blocs se recouvrent
        tol = 1.e-12
        bboxes = computeBBoxes__(arrays, zoneNames)

        c = 0
        zones = Internal.getZones(t)
        for z in zones:
            d = 0
            for z2 in zones:
                if id(z) != id(z2):
                    dim2 = Internal.getZoneDim(z2)
                    if dim2[0] == 'Structured': np = dim2[1]*dim2[2]*dim2[3]
                    else: np = dim2[1]
                    bboxc = bboxes[c]; bboxd = bboxes[d]
                    xmin1 = bboxc[0]; xmax1 = bboxc[3]
                    ymin1 = bboxc[1]; ymax1 = bboxc[4]
                    zmin1 = bboxc[2]; zmax1 = bboxc[5]
                    xmin2 = bboxd[0]; xmax2 = bboxd[3]
                    ymin2 = bboxd[1]; ymax2 = bboxd[4]
                    zmin2 = bboxd[2]; zmax2 = bboxd[5]
                    if (xmax1 > xmin2-tol and xmin1 < xmax2+tol and
                        ymax1 > ymin2-tol and ymin1 < ymax2+tol and
                        zmax1 > zmin2-tol and zmin1 < zmax2+tol):
                        com[c, d] += np
                d += 1
            c += 1

    if useCom == 'ID' or useCom == 'all':
        dict = {}
        pc = 0
        for z in zones: dict[z[0]] = pc; pc += 1
        for z in zones:
            zname = z[0]
            sr = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
            for s in sr:
                oppname = Internal.getValue(s)
                PL = Internal.getNodeFromName1(s, 'PointList')
                if PL is not None and PL[1] is not None:
                    com[dict[zname],dict[oppname]] = PL[1].size
                else:
                    com[dict[zname],dict[oppname]] = 1.

    # Equilibrage
    out = Distributor2.distribute(arrays, NProc, prescribed=aset, com=com,
                                  perfo=perfo, weight=weightlist, 
                                  algorithm=algorithm, mode=mode, nghost=nghost)
    # Sortie
    zones = Internal.getZones(t)
    procs = out['distrib']
    i = 0
    for z in zones:
        node = Internal.getNodeFromName1(z, '.Solver#Param')
        if node is not None: param = node
        else:
            param = ['.Solver#Param', None, [], 'UserDefinedData_t']
            z[2].append(param)
        v = numpy.zeros((1,1), numpy.int32); v[0,0] = procs[i]
        node = Internal.getNodeFromName1(param, 'proc')
        if node is not None:
            a = node; a[1] = v
        else:
            a = ['proc', v, [], 'DataArray_t']
            param[2].append(a)
        i += 1
    return out

#==============================================================================
# Retourne le dictionnaire proc['blocName']
# a partir d'un arbre distribue contenant des noeuds proc
#==============================================================================
def getProcDict(t, prefixByBase=False):
    """Return the proc of a zone in a dictionary proc['zoneName']."""
    proc = {}
    tp = Internal.node2PyTree(t)
    bases = Internal.getBases(t)
    for b in bases:
        zones = Internal.getNodesFromType1(b, 'Zone_t')
        if prefixByBase: prefix = b[0]+'/'
        else: prefix = ''
        for z in zones:
            nproc = Internal.getNodeFromName2(z, 'proc')
            if nproc is not None: nproc = Internal.getValue(nproc)
            else: nproc = -1
            proc[prefix+z[0]] = nproc
    return proc

#==============================================================================
# Retourne la liste des zones pour un processeur donne
#==============================================================================
def getProcList(t, NProc=None, sort=False):
    """Return the list of zones for each proc."""
    zones = Internal.getNodesFromType2(t, 'Zone_t')
    if NProc is None:
        NProc = 0
        for z in zones:
            proc = Internal.getNodeFromName2(z, 'proc')
            if proc is not None:
               proc = Internal.getValue(proc)
               NProc = max(NProc, proc+1)
    procList = []
    for s in range(NProc): procList.append([])

    if not sort: # pas de tri
        for z in zones: 
           proc = Internal.getNodeFromName2(z, 'proc') 
           if proc is not None:
            procList[Internal.getValue(proc)].append(z[0])

    else:
        # On trie les zones par taille decroissant pour chaque base pour etre conforme avec 
        # mode openmp (see reorder dans Fast.Internal)
        bases = Internal.getNodesFromType1(t,'CGNSBase_t')
        for base in bases: 
           zones = Internal.getNodesFromType1(base,'Zone_t')
           # calcul taille de la zone
           size_zone =[]
           for z in zones:
              dim = Internal.getZoneDim(z)
              if dim[0]=='Structured':
                if dim[3] == 1: kfic = 0
                else          : kfic = 2
                ndimdx = (dim[1]-4)*(dim[2]-4)*(dim[3]-kfic) 
              else: ndimdx = dim[2]
              size_zone.append(ndimdx)

           # Tri les zone par taille decroissante
           new_zones =[]
           for z in range(len(size_zone)):
             vmax    = max(size_zone)
             pos_max = size_zone.index(vmax)
             new_zones.append(zones[pos_max])
             size_zone.pop(pos_max) 
             zones.pop(pos_max) 

           for z in new_zones:
             proc = Internal.getNodeFromName2(z, 'proc')
             if proc is not None: 
                procList[Internal.getValue(proc)].append(z[0])

    return procList

#==============================================================================
# Copie la distribution de b dans a
# Match par noms
#==============================================================================
def copyDistribution(a, b):
    """Copy distribution of a in b."""
    o = Internal.copyRef(a)
    _copyDistribution(o,b)
    return o

def _copyDistribution(a, b):
    procs = getProcDict(b)
    zones = Internal.getZones(a)
    for z in zones:
        if z[0] in procs:
            proc = procs[z[0]]
            node = Internal.getNodeFromName1(z, '.Solver#Param')
            if node is not None: param = node
            else:
                param = ['.Solver#Param', None, [], 'UserDefinedData_t']
                z[2].append(param)
            v = numpy.zeros((1,1), numpy.int32); v[0,0] = proc
            node = Internal.getNodeFromName1(param, 'proc')
            if node is not None:
                a = node; a[1] = v
            else:
                a = ['proc', v, [], 'DataArray_t']
                param[2].append(a)
    return None

#==============================================================================
# addProcNode
#==============================================================================
def addProcNode(t, proc):
    """Add a "proc" node to zones with the given proc value."""
    tp = Internal.copyRef(t)
    _addProcNode(tp, proc)
    return tp

def _addProcNode(t, proc):
    zones = Internal.getZones(t)
    for z in zones:
        Internal.createUniqueChild(z, '.Solver#Param', 'UserDefinedData_t', 
                                   value=None)
        n = Internal.getNodeFromName1(z, '.Solver#Param')
        Internal.createUniqueChild(n, 'proc', 'DataArray_t', value=proc)
    return None

#==============================================================================
# getProc
# Si une seule zone: retourne proc
# Si plusieurs zones: retourne [procs]
# Si non trouve, retourne -1
#==============================================================================
def getProc(t):
    """Return the value of proc node."""
    zones = Internal.getZones(t)
    procs = []
    for z in zones:
        proc = Internal.getNodeFromName2(z, 'proc')
        if proc is None: procs.append(-1)
        else: procs.append(Internal.getValue(proc))
    if len(procs) == 1: return procs[0]
    else: return procs


def printProcStats(t, stats=None, NProc=None):
    """Print stats dictionary."""
    if stats is not None:
        dist = stats['distrib']
        if NProc is None: NProc = max(dist)+1
        if len(list(set(dist))) != NProc:
            print ('Warning: some processors are empty!')
            import sys;sys.exit()

        zones = Internal.getZones(t)
        for proc in range(NProc):
            indexes = [i for i,x in enumerate(dist) if x == proc]
            npts = 0
            lzone = []
            for i in indexes:
                z = zones[i]
                lzone.append(Internal.getName(z))
                dim = Internal.getZoneDim(z)
                npts += dim[1]*dim[2]*dim[3]
            print ('Info: proc '+str(proc)+': '+str(npts)+' points for zones ',lzone)
                
    else: # no dist
        d = getProcList(t)
        if NProc is None: NProc = len(d)
        for proc in range(NProc):
            lzone = d[proc]
            npts = 0
            for i in lzone:
                z = Internal.getNodeFromName2(t, i)
                dim = Internal.getZoneDim(z)
                npts += dim[1]*dim[2]*dim[3]
            print ('Info: proc '+str(proc)+': '+str(npts)+' points for zones ',lzone)

    return None

