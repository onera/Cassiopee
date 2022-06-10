# bufferized parallel probe
import Converter.PyTree as C
import Geom.PyTree as D
import Converter.Mpi as Cmpi
import Converter.Internal as Internal
import Generator.PyTree as G
import Converter.Filter as Filter
import Converter.Distributed as Distributed
import numpy

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
        self._probeZone = None

    # init from position X/ind+blockName/None
    def __init__(self, fileName, t=None, 
                 X=None, 
                 ind=None, blockName=None,
                 fields=None, append=True, bufferSize=5000):
        """Create a probe."""
        self.init0()
        self._bsize = bufferSize
        self._fileName = fileName
        self._fields = fields
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
        
        if t is not None and fields is not None: 
            self.checkVariables(t, fields)
        
        if self._proc is not None and self._fields is not None:
            self.createProbeZone()
        
        if self._proc is not None and self._probeZone is not None:
            self.checkFile(append=self._append)
        
    # examine la localisation des champs
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
    # OUT: ind, blockName, dist
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
    def locateProbeInd(self, t, ind, blockName):
        p = Internal.getNodeFromName2(t, blockName)
        if p is not None and not Cmpi.isZoneSkeleton__(p):
            proc = Cmpi.rank
            print('Info: probe found on proc: %d on block %s.'%(proc,blockName))
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

    # Create the probe zone with buffer size
    # IN: _proc: proc of probe
    # IN: _fields: champ a extraire
    # OUT: _probeZone: zone de stockage
    def createProbeZone(self, source=None):
        if self._mode == 0 or self._mode == 1:
            if Cmpi.rank != self._proc: return
            self._probeZone = G.cart((0,0,0), (1,1,1), (self._bsize,1,1))
            self._probeZone[0] = 'probe'
            C._initVars(self._probeZone, '{time}=-1.') # time sentinel
            # create vars in probe      
            for v in self._fields:
                C._initVars(self._probeZone, '{nodes:%s}=0.'%v)
        elif self._mode == 2 and source is not None:
            npts = C.getNPts(source)
            print('creating', self._bsize, npts)
            self._probeZone = G.cart((0,0,0), (1,1,1), (self._bsize,npts,1))
            self._probeZone[0] = 'probe'
            C._initVars(self._probeZone, '{time}=-1.') # time sentinel
            for v in self._fields:
                C._initVars(self._probeZone, v, 0.)
        return None

    # Check file, if it doesnt exist, write probe zone in it
    # else get the filecur
    # IN: _proc: proc of probe
    def checkFile(self, append=True):
        if Cmpi.rank != self._proc: return
        if append:
            create = False
            try:
                tl = Cmpi.convertFile2SkeletonTree(self._fileName)
                nodes = Internal.getNodesFromName(tl, 'GridCoordinates#*')
                self._filecur = len(nodes)
            except: create = True
        else: create = True
        if create:
            C.convertPyTree2File(self._probeZone, self._fileName)
            self._filecur = 0
        # load GC
        nodes = Distributed.readNodesFromPaths(self._fileName, ['CGNSTree/Base/probe/GridCoordinates'])
        cont = Internal.getNodeFromName2(self._probeZone, 'GridCoordinates')
        if cont is not None: cont[2] = nodes[0][2]
        # load FS
        nodes = Distributed.readNodesFromPaths(self._fileName, ['CGNSTree/Base/probe/FlowSolution'])
        cont = Internal.getNodeFromName2(self._probeZone, 'FlowSolution')
        if cont is not None: 
            cont[2] = nodes[0][2]
            pt = Internal.getNodeFromName2(self._probeZone, 'time')[1]
            a = pt > -0.5
            self._icur = numpy.count_nonzero(a)
        else: self._icur = 0

        print('Info: probe: filecur:', self._filecur)
        print('Info: probe: icur:', self._icur)
        return None

    # verifie la var list dans t
    # IN: _proc
    # IN: _blockName
    # IN: varList
    # OUT: _fields
    def checkVariables(self, t, varList):
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
    # IN: _probeZone: zone de stockage
    # IN: _ind: index of probe (static)
    # IN: _fields: nom des champs a extraire
    def extract(self, t, time):
        """Extract XYZ or Ind fields from t."""

        if self._mode == 0 or self._mode == 1: # single point
            if Cmpi.rank != self._proc: return None

            # time is in "time" of probe zone
            pt = Internal.getNodeFromName2(self._probeZone, 'time')[1]
            pt[self._icur] = time
            px = Internal.getNodeFromName2(self._probeZone, 'CoordinateX')[1]
            py = Internal.getNodeFromName2(self._probeZone, 'CoordinateY')[1]
            pz = Internal.getNodeFromName2(self._probeZone, 'CoordinateZ')[1]
            px[self._icur] = self._posX
            py[self._icur] = self._posY
            pz[self._icur] = self._posZ

            b = Internal.getNodeFromName2(t, self._blockName)

            for c, v in enumerate(self._fields):
                f = Internal.getNodeFromName2(b, v)[1]
                f = f.ravel('k')
                pf = Internal.getNodeFromName2(self._probeZone, v)[1]
                pf = pf.ravel('k')
                pf[self._icur] = f[self._ind]
                #print('value=',f[self._ind])
            self._icur += 1
            if self._icur >= self._bsize: self.flush()

        elif self._mode == 2: # single zone
            if self._probeZone is None: 
                self.createProbeZone(t)
                self.checkFile(append=self._append)

            source = Internal.getZones(t)[0] # only one for now
            
            # time is in "time" of probe zone
            pt = Internal.getNodeFromName2(self._probeZone, 'time')[1]
            pt[self._icur,:] = time

            # Set zone coordinates
            px = Internal.getNodeFromName2(self._probeZone, 'CoordinateX')[1]
            py = Internal.getNodeFromName2(self._probeZone, 'CoordinateY')[1]
            pz = Internal.getNodeFromName2(self._probeZone, 'CoordinateZ')[1]
            
            px2 = Internal.getNodeFromName2(source, 'CoordinateX')[1].ravel('k')
            py2 = Internal.getNodeFromName2(source, 'CoordinateY')[1].ravel('k')
            pz2 = Internal.getNodeFromName2(source, 'CoordinateZ')[1].ravel('k')
            
            px[self._icur,:] = px2[:]
            py[self._icur,:] = py2[:]
            pz[self._icur,:] = pz2[:]

            # Set zone fields
            for c, v in enumerate(self._fields):
                f = Internal.getNodeFromName2(source, v)[1]
                f = f.ravel('k')
                pf = Internal.getNodeFromName2(self._probeZone, v)[1]
                pf[self._icur,:] = f[:]
                
            self._icur += 1
            if self._icur >= self._bsize: self.flush()

        return None

    # flush containers of probe
    # IN: _proc: proc du bloc contenant la probe
    # IN: _probeZone: zone de la probe
    # IN: _fileName: nom du fichier
    def flush(self):
        """Flush probe to file."""
        if Cmpi.rank != self._proc: return None
        # flush containers
        gc = Internal.getNodeFromName1(self._probeZone, 'GridCoordinates')
        fc = Internal.getNodeFromName1(self._probeZone, 'FlowSolution')
        if self._icur >= self._bsize: # because buffer is out
            gc = Internal.copyNode(gc)
            gc[0] = 'GridCoordinates#%d'%self._filecur
            fc = Internal.copyNode(fc)
            fc[0] = 'FlowSolution#%d'%self._filecur
            print('Info: probe: flush #%d (full).'%self._filecur)
            paths = ['CGNSTree/Base/probe','CGNSTree/Base/probe']
            nodes = [gc,fc]
            Distributed.writeNodesFromPaths(self._fileName, paths, nodes, mode=0)
            self._filecur += 1
            C._initVars(self._probeZone, '{time}=-1.') # time sentinel
            for v in self._fields: C._initVars(self._probeZone, '{%s}=0.'%v)
            self._icur = 0
        else: # explicit flush
            nodes = [gc,fc]
            paths = ['CGNSTree/Base/probe/GridCoordinates', 'CGNSTree/Base/probe/FlowSolution']
            Distributed.writeNodesFromPaths(self._fileName, paths, nodes, mode=1)
        return None

    # read all probe fields as a single zone
    # IN: _fileName: nom du fichier
    # IN: time: instant a extraire -> all points
    # IN: index: point a extraire -> all times
    # PAS OPERATIONEL
    def read(self, time=None, ind=None):
        """Reread all data from probe file."""
        if time is not None:
            return self.readTime(time)
        elif ind is not None:
            return self.readInd(ind)
        else:
            return self.readInd(0)

    # Retourne tous les indices d'un seul instant
    def readTime(self, time):
        return None

    # Retourne tous les temps d'un seul indice
    def readInd(self, ind):
        tl = Cmpi.convertFile2SkeletonTree(self._fileName)

        # recupere le sizeTime
        nodes = Internal.getNodesFromName(tl, 'FlowSolution#*')
        ncont = len(nodes) # nbre de containers pleins
        paths = ['CGNSTree/Base/probe/FlowSolution']
        nodesTime = Distributed.readNodesFromPaths(self._fileName, paths)
        pt = Internal.getNodeFromName1(nodesTime[0], 'time')[1]
        a = pt > -0.5
        csize = numpy.count_nonzero(a)
        sizeTimeCont = pt.shape[0]
        sizeTime = csize + ncont*sizeTimeCont
        print('sizeTime', sizeTime)

        # Get all vars

        # Rebuild full zone
        out = G.cart((0,0,0), (1,1,1), (sizeTime,1,1))
        px = Internal.getNodeFromName2(out, 'CoordinateX')
        py = Internal.getNodeFromName2(out, 'CoordinateY')
        pz = Internal.getNodeFromName2(out, 'CoordinateZ')
        
        # load GCs
        nodes = Internal.getNodesFromName(tl, 'GridCoordinates#*')
        cur = 0
        for n in nodes:
            paths = ['CGNSTree/Base/probe/%s'%n[0]]
            nodesX = Distributed.readNodesFromPaths(self._fileName, paths)[0]
            px2 = Internal.getNodeFromName2(nodesX, 'CoordinateX')
            py2 = Internal.getNodeFromName2(nodesX, 'CoordinateY')
            pz2 = Internal.getNodeFromName2(nodesX, 'CoordinateZ')
            px[cur:cur+sizeTimeCont] = px2[ind,0:sizeTimeCont]
            cur += sizeTimeCont

        paths = ['CGNSTree/Base/probe/GridCoordinates']
        nodesX = Distributed.readNodesFromPaths(self._fileName, paths)
        px2 = Internal.getNodeFromName2(nodesX, 'CoordinateX')
        py2 = Internal.getNodeFromName2(nodesX, 'CoordinateY')
        pz2 = Internal.getNodeFromName2(nodesX, 'CoordinateZ')
        px[cur:cur+csize] = px2[ind,0:csize]

        # load FS
        nodes = Internal.getNodesFromName(tl, 'FlowSolution#*')
        for n in nodes:
            paths = ['CGNSTree/Base/probe/%s'%n[0]]
            nodesF = Distributed.readNodesFromPaths(self._fileName, paths)[0]
            
        paths = ['CGNSTree/Base/probe/FlowSolution']
        
        # Enleve la partie sentinelle
        pt = Internal.getNodeFromName1(nodesF[0], 'time')[1]
        a = pt > -0.5
        csize = numpy.count_nonzero(a)
        size = csize
        for n in nodesF[1:]:
            ptn = Internal.getNodeFromName1(n, 'time')
            size += ptn[1].size

        # Rebuild full zone
        out = G.cart((0,0,0), (1,1,1), (size,1,1))
        pt2 = Internal.getNodeFromName2(out, 'time')[1]
        c = 0
        for n in nodesF[1:]:
            ptn = Internal.getNodeFromName1(n, 'time')[1]
            pt2[c:c+self._bsize] = ptn[0:self._bsize]
            c += self._bsize
        pt2[c:c+csize] = pt[0:csize]
        
        pf = Internal.getNodesFromType(nodes[0], 'DataArray_t')
        nfields = len(pf)
        for p in pf:
            C._initVars(out, '%s=0.'%p[0])
            px2 = Internal.getNodeFromName2(out, p[0])[1]
            c = 0
            for n in nodes[1:]:
                pxn = Internal.getNodeFromName1(n, p[0])[1]
                px2[c:c+self._bsize] = pxn[0:self._bsize]
                c += self._bsize
            pxn = Internal.getNodeFromName1(nodes[0], p[0])[1]
            px2[c:c+csize] = pxn[0:csize]
        return out
