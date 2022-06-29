# bufferized parallel probe
import Converter.PyTree as C
import Geom.PyTree as D
import Converter.Mpi as Cmpi
import Converter.Internal as Internal
import Generator.PyTree as G
import Converter.Filter as Filter
import Converter.Distributed as Distributed
import Distributor2.PyTree as D2
import numpy, os

# Probe class
class Probe:
    """Probe for extracting data from fields during computations."""
    # defaullt initialization
    def init0(self):
        # probe mode: 
        # mode=0, probe XYZ, single point
        # mode=1, probe ind, single point
        # mode=2, probe zones
        self._mode = 2

        # -- data probe mode=0 --
        # probe Coordinates
        self._posX = None
        self._posY = None
        self._posZ = None
        # distance probe-node
        self._dist = None

        # -- data probe mode=0 ou mode=1 --
        # probe node/center index in block
        self._ind = None
        # block name containing probe
        self._blockName = None
        # the proc blockName is on
        self._proc = 0
        
        # -- data for all modes --
        # list of extracted field names
        self._fields = []
        # if true, write coords in probe else write only fields
        self._coords = True

        # file attached to
        self._fileName = None
        # current container in file to write to
        self._filecur = 0

        # internal buffer size
        self._bsize = 5000
        # current position in buffer
        self._icur = 0
        # append or not
        self._append = True

        # zones storing probe data
        self._probeZones = None

    # init from position X/ind+blockName/None
    def __init__(self, fileName, t=None, 
                 X=None, 
                 ind=None, blockName=None,
                 fields=None, append=True, 
                 bufferSize=100, writeCoords=True):
        """Create a probe."""
        self.init0()
        self._bsize = bufferSize
        self._fileName = fileName
        self._fields = fields
        self._coords = writeCoords
        self._append = append
        if fields is not None: loc = self.getFieldLoc(fields)
        else: loc = None

        # Localisation a partir de X (mode=0)
        if X is not None:
            self._mode = 0
            self._posX = X[0]
            self._posY = X[1]
            self._posZ = X[2]
        
            if t is not None and loc is not None:
                self.locateProbeXYZ(t, X, loc)

        # Localisation a partir de ind,blockName (mode=1)
        elif ind is not None and blockName is not None:
            self._mode = 1
            self.locateProbeInd(t, ind, blockName)

        else: self._mode = 2
        
        # Cree la probe et on relit le fichier uniquement en mode=0 et 1
        if self._mode == 0 or self._mode == 1:
            if self._proc is not None and self._fields is not None:
                self.createProbeZones()
        
            if self._proc is not None and self._probeZones is not None:
                self.checkFile(append=self._append)
            #if t is not None and fields is not None: 
            #    self.checkVariables(t, fields)
        
    # examine la localisation des champs
    # IN: fields
    def getFieldLoc(self, fields):
        loc = None
        for v in fields:
            vs = v.split(':')
            if len(vs) == 2 and vs[0] == 'centers':
                if loc is None: loc = 'centers'
                elif loc != 'centers': raise ValueError("probe: fields must have the same loc.")
            else:
                if loc is None: loc = 'nodes'
                elif loc != 'nodes': raise ValueError("probe: fields must have the same loc.")
        return loc

    # locate probe in t from position X
    # IN: t: pyTree
    # IN: posX, posY, posZ
    # IN: loc: loc of all fields of probe
    # OUT: ind, blockName, dist, proc
    def locateProbeXYZ(self, t, X, loc):
        P = D.point(X)
        zones = Internal.getZones(t)
        dist = 1.e16; blockName = None; ind = 0
        
        for z in zones:
            if loc == 'centers': zc = C.node2Center(z)
            else: zc = z
            hook = C.createHook(zc, function='nodes')
            (i, d) = C.nearestNodes(hook, P)
            if d[0] < dist: dist = d[0]; blockName = z[0]; ind = i[0]
            C.freeHook(hook)
            if loc == 'centers': zc = None

        # parallel
        ret = Cmpi.allgather(dist)
        dist = 1.e16; proc = 0
        for p, i in enumerate(ret):
            if i is not None and i < dist: 
                dist = i; proc = p
        if Cmpi.rank == 0: print('Info: probe found on proc: %d on block: %s.'%(proc,blockName))
        
        [ind,blockName,dist,proc] = Cmpi.bcast([ind,blockName,dist,proc], root=proc)

        # set
        self._ind = ind
        self._blockName = blockName
        self._dist = dist
        self._proc = proc
        return None
        
    # locate probe from ind and blockName
    # IN: t: pyTree
    # IN: ind, blockName
    # IN: loc: loc of all fields of probe
    # OUT: posX, dist, proc
    def locateProbeInd(self, t, ind, blockName):
        p = Internal.getNodeFromName2(t, blockName)
        if p is not None and not Cmpi.isZoneSkeleton__(p):
            proc = Cmpi.rank
            if Cmpi.rank == 0: print('Info: probe found on proc: %d on block %s.'%(proc,blockName))
        else: 
            proc = -1
            #print('Warning: probe not found on block %s.'%blockName)

        if isinstance(ind, tuple):
            b = Internal.getNodeFromName(t, blockName)
            if b is not None and Cmpi.rank == proc:
                loc = self.getFieldLoc(self._fields)
                dim = Internal.getZoneDim(b)
                if dim[0] == 'Structured':
                    if loc == 'nodes': 
                        ind = (ind[0]-1)+(ind[1]-1)*dim[1]+(ind[2]-1)*dim[2]
                    else: 
                        ind = (ind[0]-1)+(ind[1]-1)*(dim[1]-1)+(ind[2]-1)*(dim[2]-1)
                else: ind = ind[0]
            else: ind = -1

        b = Internal.getNodeFromName(t, blockName)
        if b is not None and Cmpi.rank == proc:
            px = Internal.getNodeFromName2(b, 'CoordinateX')[1]
            px = px.ravel('k')
            py = Internal.getNodeFromName2(b, 'CoordinateY')[1]
            py = py.ravel('k')
            pz = Internal.getNodeFromName2(b, 'CoordinateZ')[1]
            pz = pz.ravel('k')
            self._posX = px[ind]
            self._posY = py[ind]
            self._posZ = pz[ind]

        self._ind = ind
        self._blockName = blockName
        self._dist = 0.
        self._proc = proc
        return None

    # print information on probe
    def printInfo(self):
        """Print information on probe."""
        if Cmpi.rank == 0:
            if self._mode == 0 or self._mode == 1: 
                print('Info: probe: position X: ', self._posX, self._posY, self._posZ)
                print('Info: probe: on block:', self._blockName)
                print('Info: probe: on block global index:', self._ind)
                print('Info: probe: distance probe-node:', self._dist)
                print('Info: probe: located on proc:', self._proc)

            print('Info: probe: extracted vars:', self._fields)
            print('Info: probe: filecur:', self._filecur)
            print('Info: probe: icur:', self._icur)
        return None

    # Create the probe zones with buffer size
    # IN: _proc: proc of probe
    # IN: _fields: champ a extraire
    # OUT: _probeZones: zones de stockage
    def createProbeZones(self, source=None):
        if self._mode == 0 or self._mode == 1:
            if Cmpi.rank != self._proc: return
            z = G.cart((0,0,0), (1,1,1), (self._bsize,1,1))
            z[0] = 'probe'
            D2._addProcNode(z, self._proc)
            self._probeZones = [z]
            C._initVars(self._probeZones, '{time}=-1.') # time sentinel
            # create vars in probe
            for v in self._fields:
                v = v.split(':')
                if len(v) == 2: v = v[1]
                else: v = v[0]
                C._initVars(self._probeZones, '{nodes:%s}=0.'%v)
        elif self._mode == 2 and source is not None:
            zones = Internal.getZones(source)
            self._probeZones = []
            for z in zones:
                npts = C.getNPts(z)
                ncells = C.getNCells(z)
                if self.getFieldLoc(self._fields) == 'nodes':
                    zp = G.cart((0,0,0), (1,1,1), (self._bsize,npts,1))
                else: 
                    zp = G.cart((0,0,0), (1,1,1), (self._bsize,ncells,1))
                zp[0] = '%s'%z[0] # name of probe is identical to name of source
                D2._addProcNode(zp, Cmpi.rank)
                self._probeZones.append(zp)
                C._initVars(zp, '{time}=-1.') # time sentinel

                for v in self._fields:
                    v = v.split(':')
                    if len(v) == 2: v = v[1]
                    else: v = v[0]
                    C._initVars(zp, v, 0.)

        return None

    # Check file, if it doesnt exist, write probe zone in it
    # else get the filecur
    # IN: _proc: proc of probe
    def checkFile(self, append=False):
        if append:
            create = False
            if self._mode == 0 or self._mode == 1:
                isFilePresent = os.path.exists(self._fileName)
            else:
                if Cmpi.rank == 0: isFilePresent = os.path.exists(self._fileName)
                else: isFilePresent = None
                isFilePresent = Cmpi.bcast(isFilePresent, root=0)
    
            if isFilePresent:
                if self._mode == 0 or self._mode == 1:
                    tl = Distributed.convertFile2SkeletonTree(self._fileName)
                else: 
                    tl = Cmpi.convertFile2SkeletonTree(self._fileName)
                zones = Internal.getZones(tl)
                if len(zones) > 0:
                    nodes = Internal.getNodesFromName(zones[0], 'GridCoordinates#*')
                    self._filecur = len(nodes)
                if self._mode == 2:
                    self._filecur = Cmpi.allreduce(self._filecur, op=Cmpi.MAX)
            else: create = True
        else: create = True
        if create:
            if self._mode == 0 or self._mode == 1:
                if Cmpi.rank == self._proc:
                    C.convertPyTree2File(self._probeZones, self._fileName)
            else:
                Cmpi.convertPyTree2File(self._probeZones, self._fileName)
            self._filecur = 0
            return None

        
        #  Nettoyage + ne conserve que les zones du proc
        to = C.newPyTree()
        bases = Internal.getBases(tl)
        for b in bases:
            bl = Internal.newCGNSBase(b[0], parent=to)
            zones = Internal.getZones(b)
            for z in zones:
                Internal._rmNodesByName(z, 'GridCoordinates#*')
                Internal._rmNodesByName(z, 'FlowSolution#*')
                Internal._rmNodesByName(z, 'FlowSolution#Centers#*')
                if D2.getProc(z) == Cmpi.rank: bl[2].append(z)
        
        tl = to
        self._probeZones = Internal.getZones(tl)

        # load GC
        paths = []
        zones = Internal.getZones(tl)
        for z in zones:
            paths += ['CGNSTree/Base/%s/GridCoordinates'%z[0]]
        nodes = Distributed.readNodesFromPaths(self._fileName, paths)
        for c, z in enumerate(self._probeZones):
            cont = Internal.getNodeFromName2(z, 'GridCoordinates')
            if cont is not None: cont[2] = nodes[c][2]

        # load FS
        paths = []
        for z in zones:
            paths += ['CGNSTree/Base/%s/FlowSolution'%z[0]]
        nodes = Distributed.readNodesFromPaths(self._fileName, paths)
        for c, z in enumerate(self._probeZones):
            cont = Internal.getNodeFromName2(z, 'FlowSolution')
            if cont is not None: cont[2] = nodes[c][2]

        # load FS center
        #paths = []
        #for z in zones:
        #    if Internal.getNodeFromName1(z, 'FlowSolution#Centers') is not None:
        #        paths += ['CGNSTree/Base/%s/FlowSolution#Centers'%z[0]]
        #if paths != []:    
        #    nodes = Distributed.readNodesFromPaths(self._fileName, paths)
        #    for c, z in enumerate(self._probeZones):
        #        cont = Internal.getNodeFromName2(z, 'FlowSolution#Centers')
        #        if cont is not None: cont[2] = nodes[c][2]

        if len(self._probeZones) > 0:
            cont = Internal.getNodeFromName2(self._probeZones[0], 'FlowSolution')
            if cont is not None: 
                pt = Internal.getNodeFromName2(cont, 'time')[1]
                sh = pt.shape
                if len(sh) == 2: pt = pt[:,0]
                a = pt > -0.5
                self._icur = numpy.count_nonzero(a)
            else: self._icur = 0
        else: self._icur = 0

        if self._mode == 2:
            self._icur = Cmpi.allreduce(self._icur, op=Cmpi.MAX)

        #print('Info: probe: filecur:', self._filecur, flush=True)
        #print('Info: probe: icur:', self._icur, flush=True)

        return None

    # verifie la var list dans t
    # IN: _proc
    # IN: _blockName
    # IN: varList
    # OUT: _fields
    def checkVariables(self, t, varList):
        if self._mode == 2: return
        if Cmpi.rank != self._proc: return
        if self._blockName is None: return
        block = Internal.getNodeFromName2(t, self._blockName)
        if block is None: return
        if varList is None: # get all vars from blockName
            varList = C.getVarNames(block, excludeXYZ=True)[0]
        self._fields = []
        for v in varList:
            vs = v.split(':')
            contName = Internal.__FlowSolutionNodes__
            if len(vs) == 2 and vs[0] == 'centers': 
                contName = Internal.__FlowSolutionCenters__
                v = vs[1]
            cont = Internal.getNodeFromName1(block, contName)
            if cont is None:
                raise ValueError("probe: can not find solution container in t.")
            var = Internal.getNodeFromName1(cont, v)
            if var is None:
                raise ValueError("probe: can not find field %s in t."%v)
            self._fields.append(v)
        return None

    # trigger extraction on current t
    # IN: _proc: proc du bloc contenant la probe
    # IN: _probeZones: zone de stockage
    # IN: _ind: index of probe (static)
    # IN: _fields: nom des champs a extraire
    def extract(self, t, time):
        """Extract XYZ or Ind fields from t."""

        if self._mode == 0 or self._mode == 1: # single point
            if Cmpi.rank != self._proc: return None

            # time is in "time" of probe zone
            pzone = self._probeZones[0]
            pt = Internal.getNodeFromName2(pzone, 'time')[1]
            pt[self._icur] = time
            px = Internal.getNodeFromName2(pzone, 'CoordinateX')[1]
            py = Internal.getNodeFromName2(pzone, 'CoordinateY')[1]
            pz = Internal.getNodeFromName2(pzone, 'CoordinateZ')[1]
            px[self._icur] = self._posX
            py[self._icur] = self._posY
            pz[self._icur] = self._posZ

            b = Internal.getNodeFromName2(t, self._blockName)

            for c, v in enumerate(self._fields):
                v = v.split(':')
                if len(v) == 2: v = v[1]
                else: v = v[0]
                f = Internal.getNodeFromName2(b, v)[1]
                f = f.ravel('k')
                pf = Internal.getNodeFromName2(pzone, v)[1]
                pf = pf.ravel('k')
                pf[self._icur] = f[self._ind]
                #print('value=',f[self._ind])
            self._icur += 1
            if self._icur >= self._bsize: self.flush()

        elif self._mode == 2: # single zone
            if self._probeZones is None: 
                self.createProbeZones(t)
                self.checkFile(append=self._append)
            Cmpi.barrier()

            # time is in "time" of probe zones
            source = Internal.getZones(t)
            for c, pz in enumerate(self._probeZones):
                sourcez = source[c]
                pt = Internal.getNodeFromName2(pz, 'time')[1]
                pt[self._icur,:] = time

                # Set zone coordinates
                if self._coords:
                    ptx = Internal.getNodeFromName2(pz, 'CoordinateX')[1]
                    pty = Internal.getNodeFromName2(pz, 'CoordinateY')[1]
                    ptz = Internal.getNodeFromName2(pz, 'CoordinateZ')[1]
                    ptx2 = Internal.getNodeFromName2(sourcez, 'CoordinateX')[1].ravel('k')
                    pty2 = Internal.getNodeFromName2(sourcez, 'CoordinateY')[1].ravel('k')
                    ptz2 = Internal.getNodeFromName2(sourcez, 'CoordinateZ')[1].ravel('k')
                    if ptx2.size == ptx.shape[1]:
                        ptx[self._icur,:] = ptx2[:]
                        pty[self._icur,:] = pty2[:]
                        ptz[self._icur,:] = ptz2[:]
                    else:
                        sourcezp = C.node2Center(sourcez)
                        ptx2 = Internal.getNodeFromName2(sourcezp, 'CoordinateX')[1].ravel('k')
                        pty2 = Internal.getNodeFromName2(sourcezp, 'CoordinateY')[1].ravel('k')
                        ptz2 = Internal.getNodeFromName2(sourcezp, 'CoordinateZ')[1].ravel('k')
                        ptx[self._icur,:] = ptx2[:]
                        pty[self._icur,:] = pty2[:]
                        ptz[self._icur,:] = ptz2[:]                    

                # Set zone fields
                for c, v in enumerate(self._fields):
                    v = v.split(':')
                    if len(v) == 2: v = v[1]
                    else: v = v[0]
                    f = Internal.getNodeFromName2(sourcez, v)[1]
                    f = f.ravel('k')
                    pf = Internal.getNodeFromName2(pz, v)[1]
                    pf[self._icur,:] = f[:]
                
            self._icur += 1
            if self._icur >= self._bsize: self.flush()

        return None

    # flush containers of probe
    # IN: _proc: proc du bloc contenant la probe
    # IN: _probeZones: zone de la probe
    # IN: _fileName: nom du fichier
    def flush(self):
        """Flush probe to file."""
        if self._mode == 0 or self._mode == 1:
            if Cmpi.rank != self._proc: return None
            print('Info: probe: flush #%d.'%self._filecur)
            self.flush__()
        else:
            if Cmpi.rank == 0: print('Info: probe: flush #%d.'%self._filecur)
            Cmpi.seq(self.flush__)

    def flush__(self):
        if self._icur >= self._bsize: # because buffer is out
            for pzone in self._probeZones:
                gc = Internal.getNodeFromName1(pzone, 'GridCoordinates')
                fc = Internal.getNodeFromName1(pzone, 'FlowSolution')
                gc = Internal.copyNode(gc)
                gc[0] = 'GridCoordinates#%d'%self._filecur
                fc = Internal.copyNode(fc)
                fc[0] = 'FlowSolution#%d'%self._filecur
                paths = ['CGNSTree/Base/%s'%pzone[0],'CGNSTree/Base/%s'%pzone[0]]
                nodes = [gc,fc]
                Distributed.writeNodesFromPaths(self._fileName, paths, nodes, None, -1, 0)
                C._initVars(pzone, '{time}=-1.') # time sentinel
                for v in self._fields:
                    v = v.split(':')
                    if len(v) == 2: v = v[1]
                    else: v = v[0] 
                    C._initVars(pzone, '{%s}=0.'%v)
            self._filecur += 1
            self._icur = 0

        else: # explicit flush
            for pzone in self._probeZones:
                gc = Internal.getNodeFromName1(pzone, 'GridCoordinates')
                fc = Internal.getNodeFromName1(pzone, 'FlowSolution')
                nodes = [gc,fc]
                paths = ['CGNSTree/Base/%s/GridCoordinates'%pzone[0], 'CGNSTree/Base/%s/FlowSolution'%pzone[0]]
                Distributed.writeNodesFromPaths(self._fileName, paths, nodes, None, -1, 1)
        
        return None

    # read all probe times as a single zone
    # IN: _fileName: nom du fichier
    # IN: time: instant a extraire -> all points
    # IN: index: point a extraire -> all times
    def read(self, time=None, ind=None, probeName=None):
        """Reread all data from probe file."""
        if time is not None:
            return self.readTime(time)
        elif ind is not None:
            return self.readInd(ind, probeName)
        else:
            return self.readInd(0, probeName)

    # Retourne tous les indices d'un seul instant
    def readTime(self, time):
        return None

    # Retourne tous les temps d'un seul indice
    def readInd(self, ind=None, probeName=None):
        tl = Cmpi.convertFile2SkeletonTree(self._fileName)

        # Recupere la probe a recuperer
        zones = Internal.getZones(tl)
        if probeName is not None:
            if isinstance(probeName, int): pz = zones[probeName]
            else: pz = Internal.getNodeFromName2(tl, probeName)
        else: pz = zones[0]
        dim = Internal.getZoneDim(pz)
        sizeNPts = dim[2]

        if ind is not None:
            if isinstance(ind, int): ind = [ind]
            sizeNPts = min(len(ind), sizeNPts)

        # recupere le sizeTime (nbre d'instants dans le fichier)
        nodes = Internal.getNodesFromName(pz, 'FlowSolution#*')
        ncont = len(nodes) # nbre de containers pleins
        paths = ['CGNSTree/Base/%s/FlowSolution'%pz[0]]
        nodesTime = Distributed.readNodesFromPaths(self._fileName, paths)
        pt = Internal.getNodeFromName1(nodesTime[0], 'time')[1]
        sh = pt.shape
        if len(sh) == 2: pt = pt[:,0]
        a = pt > -0.5
        csize = numpy.count_nonzero(a)
        sizeTimeCont = pt.shape[0]
        sizeTime = csize + ncont*sizeTimeCont
        #print('sizeTime', sizeTime)
        #print('sizeNPts', sizeNPts)
        #print('sizeTimeCont', sizeTimeCont)
        #print('ncont', ncont)
        #print('csize', csize)

        # Rebuild full 1D zone (time only)

        out = G.cart((0,0,0), (1,1,1), (sizeTime,sizeNPts,1))
        ptx = Internal.getNodeFromName2(out, 'CoordinateX')[1]
        pty = Internal.getNodeFromName2(out, 'CoordinateY')[1]
        ptz = Internal.getNodeFromName2(out, 'CoordinateZ')[1]
        fcont = Internal.createUniqueChild(out, 'FlowSolution', 'FlowSolution_t')

        # load GCs
        nodes = Internal.getNodesFromName(pz, 'GridCoordinates#*')
        cur = 0
        for n in nodes:
            paths = ['CGNSTree/Base/%s/%s'%(pz[0],n[0])]
            nodesX = Distributed.readNodesFromPaths(self._fileName, paths)[0]
            ptx2 = Internal.getNodeFromName2(nodesX, 'CoordinateX')[1]
            pty2 = Internal.getNodeFromName2(nodesX, 'CoordinateY')[1]
            ptz2 = Internal.getNodeFromName2(nodesX, 'CoordinateZ')[1]
            sh = ptx2.shape
            if len(sh) == 1:
                ptx[cur:cur+sizeTimeCont] = ptx2[0:sizeTimeCont]
            elif ind is None:
                ptx[cur:cur+sizeTimeCont,:] = ptx2[0:sizeTimeCont,:]
            elif sizeNPts == 1:
                ptx[cur:cur+sizeTimeCont] = ptx2[0:sizeTimeCont,ind[0]]
            else:
                for i in ind: ptx[cur:cur+sizeTimeCont,i] = ptx2[0:sizeTimeCont,i]
            cur += sizeTimeCont

        paths = ['CGNSTree/Base/%s/GridCoordinates'%pz[0]]
        nodesX = Distributed.readNodesFromPaths(self._fileName, paths)[0]
        ptx2 = Internal.getNodeFromName2(nodesX, 'CoordinateX')[1]
        pty2 = Internal.getNodeFromName2(nodesX, 'CoordinateY')[1]
        ptz2 = Internal.getNodeFromName2(nodesX, 'CoordinateZ')[1]
        sh = ptx2.shape
        if len(sh) == 1:
            ptx[cur:cur+csize] = ptx2[0:csize]
        elif ind is None:
            ptx[cur:cur+csize,:] = ptx2[0:csize,:]
        elif sizeNPts == 1:
            ptx[cur:cur+csize] = ptx2[0:csize,ind[0]]
        else:
            for i in ind: ptx[cur:cur+csize,i] = ptx2[0:csize,i]

        # load FS
        cur = 0
        nodes = Internal.getNodesFromName(pz, 'FlowSolution#*')
        for n in nodes:
            paths = ['CGNSTree/Base/%s/%s'%(pz[0],n[0])]
            nodesF = Distributed.readNodesFromPaths(self._fileName, paths)[0]
            for i in nodesF[2]:
                if i[3] == 'DataArray_t':
                    pf = Internal.getNodeFromName2(fcont, i[0])
                    if pf is None:
                        if sizeNPts == 1:
                            pf = Internal.newDataArray(i[0], value=numpy.zeros( (sizeTime), dtype=numpy.float64, order='F'), parent=fcont)
                        else:
                            pf = Internal.newDataArray(i[0], value=numpy.zeros( (sizeTime,sizeNPts), dtype=numpy.float64, order='F'), parent=fcont)
                    pf = pf[1]
                    sh = i[1].shape
                    if len(sh) == 1:
                        pf[cur:cur+sizeTimeCont] = i[1][0:sizeTimeCont]
                    elif ind is None:
                        pf[cur:cur+sizeTimeCont,:] = i[1][0:sizeTimeCont,:]
                    elif sizeNPts == 1:
                        pf[cur:cur+sizeTimeCont] = i[1][0:sizeTimeCont,ind[0]]
                    else:
                        for i in ind: pf[cur:cur+sizeTimeCont,i] = i[1][0:sizeTimeCont,i]
            cur += sizeTimeCont

        paths = ['CGNSTree/Base/%s/FlowSolution'%pz[0]]
        nodesF = Distributed.readNodesFromPaths(self._fileName, paths)[0]
        for i in nodesF[2]:
            if i[3] == 'DataArray_t':
                pf = Internal.getNodeFromName2(fcont, i[0])[1]
                sh = i[1].shape
                if len(sh) == 1:
                    pf[cur:cur+csize] = i[1][0:csize]
                elif ind is None:
                    pf[cur:cur+csize,:] = i[1][0:csize,:]
                elif sizeNPts == 1:
                    pf[cur:cur+csize] = i[1][0:csize,ind[0]]
                else:
                    for i in ind: pf[cur:cur+csize,i] = i[1][0:csize,i]
        return out
