"""Distribution module for Cassiopee package.
"""
import Distributor2
import Converter.Internal as Internal
import Converter.PyTree as C
import Generator as G
import numpy
__version__ = Distributor2.__version__

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

# Retourne une cle unique pour le dictionnaire de com
def addCom__(comd, c, d, NBlocs, vol):
    key = c+d*NBlocs
    if key in comd: comd[key] += vol
    else: comd[key] = vol

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
def distribute(t, NProc, prescribed=None, perfo=None, weight=None, useCom='match',
               algorithm='graph', mode='nodes', nghost=0, tbb=None):
    """Distribute a pyTree over processors.
    Usage: distribute(t, NProc, prescribed=None, perfo=None, weight=None, useCom='all', algorithm='graph', mode='nodes', nghost=0)"""
    tp = Internal.copyRef(t)
    out = _distribute(tp, NProc, prescribed=prescribed, perfo=perfo,
                      weight=weight, useCom=useCom, algorithm=algorithm,
                      mode=mode, nghost=nghost, tbb=tbb)
    return tp, out

# in place version
def _distribute(t, NProc, prescribed=None, perfo=None, weight=None, useCom='match',
                algorithm='graph', mode='nodes', nghost=0, tbb=None):
    """Distribute a pyTree over processors.
    Usage: _distribute(t, NProc, prescribed=None, perfo=None, weight=None, useCom='all', algorithm='graph', mode='nodes', nghost=0)"""

    (nbPts, aset, com, comd, weightlist) = getData__(t, NProc, prescribed, weight, useCom, mode, tbb)

    # Equilibrage
    out = Distributor2.distribute(nbPts, NProc, prescribed=aset,
                                  com=com, comd=comd,
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
        v = numpy.zeros((1,1), dtype=Internal.E_NpyInt); v[0,0] = procs[i]
        node = Internal.getNodeFromName1(param, 'proc')
        if node is not None:
            a = node; a[1] = v
        else:
            a = ['proc', v, [], 'DataArray_t']
            param[2].append(a)
        i += 1
    return out

# Internal: get data from tree
def getData__(t, NProc, prescribed=None, weight=None, useCom='match', mode='nodes', tbb=None):
    zones = Internal.getZones(t)
    nbPts = []; arrays = []; zoneNames = []; aset = []; weightlist = [] # weight for all zones
    for z in zones:
        zname = z[0]
        zoneNames.append(zname)
        if prescribed is not None: aset.append(prescribed.get(zname,-1))
        else: aset.append(-1)

        if weight is not None: weightlist.append(weight.get(zname, 1))
        else: weightlist.append(1)

        a = C.getFields(Internal.__GridCoordinates__, z, api=3)
        if a == [[]]: arrays.append(None)
        else: arrays.append(a[0])

        if mode == 'cells': nbPts.append(C.getNCells(z))
        else: nbPts.append(C.getNPts(z))

        #if useCom == 'overlap' or useCom == 'bbox':
        #    a = C.getFields(Internal.__GridCoordinates__, z, api=2)
        #    if a == [[]]: # no coord present in z
        #        print('Warning: no coordinates found. You shouldnt use useCom=overlap or bbox with skeleton tree.')
        #        useCom = 'match'
        #        if mode == 'cells': arrays.append(C.getNCell(z))
        #        else: arrays.append(C.getNPts(z))
        #    else: arrays.append(a[0])
        #else:
        #    if mode == 'cells': arrays.append(C.getNCell(z))
        #    else: arrays.append(C.getNPts(z))

    Nb = len(nbPts)
    com = None; comd = {}

    if useCom == 'match' or useCom == 'all':
        # Formation des coms - raccords coincidents
        tpp, typen = Internal.node2PyTree(t)
        bases = Internal.getBases(tpp)
        zc = 0; c = 0
        for b in bases:
            zones = Internal.getNodesFromType1(b, 'Zone_t')
            mdict = {}
            pc = 0
            for z in zones: mdict[z[0]] = pc; pc += 1

            for z in zones:
                match = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t') # forcement structure
                for m in match:
                    donorName = Internal.getValue(m)
                    if donorName in mdict: d = mdict[donorName]+zc
                    else: d = -1
                    node = Internal.getNodeFromName1(m, 'PointRange')
                    if node is not None and node[1] is not None:
                        win = node[1]
                        w = Internal.range2Window(win)
                        vol = (w[1]-w[0]+1)*(w[3]-w[2]+1)*(w[5]-w[4]+1)
                        if d != -1: addCom__(comd, c, d, Nb, vol)
                match = Internal.getNodesFromType2(z, 'GridConnectivity_t') # non structure
                for m in match:
                    donorName = Internal.getValue(m)
                    if donorName in mdict: d = mdict[donorName]+zc
                    else: d = -1
                    node = Internal.getNodeFromName1(m, 'PointList')
                    if node is not None and node[1] is not None:
                        vol = node[1].size
                        if d != -1: addCom__(comd, c, d, Nb, vol)
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
                                    addCom__(comd, c, d, Nb, vol)
                            d += 1
            c += 1

    if useCom == 'bbox':
        # Formation des coms - si les blocs se recouvrent
        tol = 1.e-12
        if tbb is None:
            # Calcul les bbox a partir de l'arbre
            bboxes = computeBBoxes__(arrays, zoneNames)
        else:
            # Calcul les bbox a partir d'un arbre de bbox global
            bboxes = []
            for z in Internal.getZones(tbb):
                minx = C.getMinValue(z, 'CoordinateX')
                miny = C.getMinValue(z, 'CoordinateY')
                minz = C.getMinValue(z, 'CoordinateZ')
                maxx = C.getMaxValue(z, 'CoordinateX')
                maxy = C.getMaxValue(z, 'CoordinateY')
                maxz = C.getMaxValue(z, 'CoordinateZ')
                bboxes.append([minx,miny,minz,maxx,maxy,maxz,z[0],True])

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
                        addCom__(comd, c, d, Nb, np)
                d += 1
            c += 1

    if useCom == 'ID' or useCom == 'all' or useCom == 'match':
        mdict = {}
        pc = 0
        for z in zones: mdict[z[0]] = pc; pc += 1
        for z in zones:
            zname = z[0]
            sr = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
            for s in sr:
                oppname = Internal.getValue(s)
                PL = Internal.getNodeFromName1(s, 'PointList')
                if PL is not None and PL[1] is not None:
                    addCom__(comd, mdict[zname], mdict[oppname], Nb, PL[1].size)
                else:
                    addCom__(comd, mdict[zname], mdict[oppname], Nb, 1)
    return (nbPts, aset, com, comd, weightlist)

#==============================================================================
# Retourne le dictionnaire proc['blocName']
# a partir d'un arbre distribue contenant des noeuds proc
#==============================================================================
def getProcDict(t, prefixByBase=False):
    """Return the proc of a zone in a dictionary proc['zoneName']."""
    proc = {}
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
        # On trie les zones par taille decroissante
        bases = Internal.getNodesFromType1(t, 'CGNSBase_t')
        for base in bases:
            zones = Internal.getNodesFromType1(base, 'Zone_t')

            # calcul taille de la zone
            sizeZone = []
            for z in zones:
                dim = Internal.getZoneDim(z)
                if dim[0] == 'Structured':
                    if dim[3] == 1: kfic = 0
                    else: kfic = 2
                    ndimdx = (dim[1]-4)*(dim[2]-4)*(dim[3]-kfic)
                else: ndimdx = dim[2]
                sizeZone.append( (z,ndimdx) )

            # Tri les zones par taille decroissante
            newZones = sorted(sizeZone, key=lambda x: x[1], reverse=True)

            for z in newZones:
                z = z[0]
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
    """Copy distribution of a in b."""
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
            v = numpy.zeros((1,1), dtype=Internal.E_NpyInt); v[0,0] = proc
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
    """Add a "proc" node to zones with the given proc value."""
    zones = Internal.getZones(t)
    for z in zones:
        Internal.createUniqueChild(z, '.Solver#Param', 'UserDefinedData_t', value=None)
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

#==============================================================================
# print infos from stats dictionary
#==============================================================================
def printProcStats(t, stats=None, NProc=None):
    """Print stats dictionary."""
    if stats is not None:
        dist = stats['distrib']
        if NProc is None: NProc = max(dist)+1
        if len(list(set(dist))) != NProc:
            print ('Warning: some processors are empty!')

        zones = Internal.getZones(t)
        for proc in range(NProc):
            indexes = [i for i,x in enumerate(dist) if x == proc]
            npts = 0
            lzone = []
            for i in indexes:
                z = zones[i]
                lzone.append(Internal.getName(z))
                dim = Internal.getZoneDim(z)
                if dim[0] == 'Structured': npts += dim[1]*dim[2]*dim[3]
                else: npts += dim[1]
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
                if dim[0] == 'Structured': npts += dim[1]*dim[2]*dim[3]
                else: npts += dim[1]
            print ('Info: proc '+str(proc)+': '+str(npts)+' points for zones ',lzone)
    return None

#================================================================================
# Get stats from tree with proc nodes (redone here from stats.cpp)
#================================================================================
def stats(t, useCom='match', mode='nodes'):
    NProc = 0; nbTot = 0
    zones = Internal.getZones(t)
    nzones = len(zones)
    nbPts = numpy.empty( (nzones), dtype=Internal.E_NpyInt )
    out = numpy.empty( (nzones), dtype=Internal.E_NpyInt )

    for c, z in enumerate(zones):
        if mode == 'nodes': np = C.getNPts(z)
        else: np = C.getNCells(z)
        nbTot += np
        p = getProc(z)
        NProc = max(NProc, p)
        out[c] = p
        nbPts[c] = np
    NProc += 1

    # meanPtsPerProc, nbNodePerProc
    meanPtsPerProc = nbTot*1./NProc
    nbNodePerProc = numpy.zeros( (NProc), dtype=numpy.float64 )
    for c, p in enumerate(out):
        nbNodePerProc[p] += nbPts[c]

    # empty
    empty = 0
    for i in range(NProc):
        if nbNodePerProc[i] == 0: empty = 1

    varMin = 1.e6; varMax = 0.; varRMS = 0.
    for i in range(NProc):
        v = abs(nbNodePerProc[i] - meanPtsPerProc)
        varMin = min(varMin, v)
        varMax = max(varMax, v)
        varRMS = varRMS + v*v

    varMin = varMin / meanPtsPerProc
    varMax = varMax / meanPtsPerProc
    varRMS = numpy.sqrt(varRMS) / (NProc*meanPtsPerProc);
    volRatio = 0.

    (nbPts, aset, com, comd, weightlist) = getData__(t, NProc, None, None, useCom, mode)

    if comd is not None:
        allkeys = comd.keys()
        size = len(allkeys)
        volComd = numpy.empty((2*size), dtype=Internal.E_NpyInt)
        for i, k in enumerate(allkeys):
            volComd[2*i] = k
            volComd[2*i+1] = comd[k]

        volTot = 0.; nptsCom = 0.
        for v in range(size):
            v1 = volComd[2*v]; volcom = volComd[2*v+1];
            k = int(v1/nzones)
            i = v1-k*nzones
            proci = out[i]
            prock = out[k]
            volTot += volcom
            if proci != prock: nptsCom += volcom

    volRatio = nptsCom / volTot

    return (varMin, varMax, varRMS, volRatio, empty)

#==================================================================================
# print stats from tree
#==================================================================================
def printStats(t, useCom='match', mode='nodes'):
    """Print stats from tree distribution."""
    (varMin, varMax, varRMS, volRatio, empty) = stats(t, useCom, mode)
    print("Info: varMin=%f%%, varMax=%f%%, varRMS=%f%%"%(100*varMin, 100*varMax, 100*varRMS))
    print("Info: external com ratio=%f%%"%(volRatio*100))
    if empty == 1: print("Warning: at least one processor is empty!")
    return (varMin, varMax, varRMS, volRatio, empty)

#==================================================================================
# distribute pytrees in files
#==================================================================================
def _checkNcellsNptsPerProc(ts, NP, isAtCenter=False):
    # Calculate and prints the number of cells & points for each proc. Provides a human-readable summary of the MPI distribution.
    import Converter.Mpi as Cmpi
    NPTS   = numpy.zeros(NP, dtype=Internal.E_NpyInt)
    NCELLS = numpy.zeros(NP, dtype=Internal.E_NpyInt)
    # Done by zone for flexibility
    for z in Internal.getZones(ts):
        proc_num        = Cmpi.getProc(z)
        NPTS[proc_num] += C.getNPts(z)
        if not isAtCenter: NCELLS[proc_num] += C.getNCells(z)
        else:              NCELLS[proc_num]  = NPTS[proc_num]

    NPTS     = Cmpi.allreduce(NPTS  ,op=Cmpi.SUM)
    NCELLS   = Cmpi.allreduce(NCELLS,op=Cmpi.SUM)
    NptsTot  = numpy.sum(NPTS)
    NcellsTot= numpy.sum(NCELLS)
    ncells_percent=[]
    if Cmpi.rank == 0:
        for i in range(NP):
            ncells_percent.append(NCELLS[i]/NcellsTot*100.)
            if isAtCenter: print('Rank {} has {} cells & {} % of cells'.format(i,int(NCELLS[i]),round(ncells_percent[i],2)))
            else:          print('Rank {} has {} points & {} cells & {} % of cells'.format(i,int(NPTS[i]),int(NCELLS[i]),round(ncells_percent[i],2)))

        if isAtCenter:   print('All points: {} million cells'.format(NcellsTot/1.e6))
        else:            print('All points: {} million points & {} million cells'.format(NptsTot/1.e6,NcellsTot/1.e6))
        print('Range of % of cells: {} - {}'.format(round(min(ncells_percent),2),round(max(ncells_percent),2)))

    return None


def _write2pathLocal__(tsLocal, tLocal):
    # Modifies the .Solver#Param/proc only in the files
    import Converter.Filter as Filter
    paths = []; ns = []
    bases = Internal.getBases(tsLocal)
    for b in bases:
        zones = Internal.getZones(b)
        for z in zones:
            nodes = Internal.getNodesFromName2(z, 'proc')
            for n in nodes:
                p = 'CGNSTree/%s/%s/.Solver#Param/proc'%(b[0],z[0])
                paths.append(p); ns.append(n)
    Filter.writeNodesFromPaths(tLocal, paths, ns, maxDepth=0, mode=1)


def _distributeSkeletonTree(tIn, NP, algorithm='graph', useCom='ID'):
    """Distribute PyTrees over NP processors. t is a list with the file names.
    Usage: _distributeSkeletonTree(t=[], NP, algorithm='graph', useCom='all')"""
    import Converter.Mpi as Cmpi
    fileNameLength = len(tIn)
    for fileName in tIn:
        ts = Cmpi.convertFile2SkeletonTree(fileName, maxDepth=3)
        if fileName == tIn[0]:
            stats = _distribute(ts, NP, algorithm=algorithm, useCom=useCom)
            tcs = Internal.copyTree(ts)
        if fileNameLength>1 and fileName!=tIn[0]: _copyDistribution(ts, tcs)
        _write2pathLocal__(ts, fileName)
    _checkNcellsNptsPerProc(tcs, NP)
    return None
