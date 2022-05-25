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
        # probe Coordinates
        self._posX = None
        self._posY = None
        self._posZ = None
    
        # probe node/center index in block
        self._ind = None
        # block name containing probe
        self._blockName = None
        # the proc block is on
        self._proc = 0
        # distance probe-node
        self._dist = None

        # pointer on pyTree probe is attached to
        self._t = None

        # list of extracted field names
        self._fields = []
        # list of pointers to corresponding numpys in t
        self._pfields = []

        # file attached to
        self._fileName = None
        # current container in file to write to
        self._filecur = 0

        # internal buffer size
        self._bsize = 5000
        # current position in buffer
        self._icur = 0

        # zone storing probe data
        self._pZone = None

    # init from position X or index
    def __init__(self, fileName, t=None, 
                 X=None, 
                 ind=None, blockName=None, proc=None,  
                 fields=None, append=True, bufferSize=5000):
        """Create a probe."""
        self.init0()
        self._bsize = bufferSize
        self._fileName = fileName
        self._t = t
        self._fields = fields
        if fields is not None: loc = self.getFieldLoc(fields)
        else: loc = None
        
        # Localication a partir de X
        if X is not None:
            self._posX = X[0]
            self._posY = X[1]
            self._posZ = X[2]
        
            if t is not None and loc is not None:
                self.locateProbe(t, X, loc)

        # Localisation a partir de ind,blockName,proc
        if ind is not None and blockName is not None:
            if proc is None: proc = 0
            self.locateProbe2(t, ind, blockName, proc)

        if fields is not None: self.checkVariables(fields)
        
        if self._proc is not None and self._fields is not None:
            self.createProbeZone()
        
        if self._proc is not None and self._pZone is not None:
            self.checkFile(append=append)
        
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
    def locateProbe(self, t, X, loc):
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
        print('Info: probe found on proc:', proc)
        
        [ind,blockName,dist,proc] = Cmpi.bcast([ind,blockName,dist,proc], root=proc)

        # set
        self._ind = ind
        self._blockName = blockName
        self._dist = dist
        self._proc = proc
        return None
        
    # locate probe from ind and blockName
    def locateProbe2(self, t, ind, blockName, proc=0):
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
            ind = Cmpi.bcast(ind, root=proc)

        self._ind = ind
        self._blockName = blockName
        self._dist = 0.
        self._proc = proc
        self.checkVariables(self._fields)
        return None

    # print information on probe
    def printInfo(self):
        """Print information on probe."""
        if Cmpi.rank != self._proc: return
        print('Info: probe: Position X: ', self._posX, self._posY, self._posZ)
        print('Info: probe: Block:', self._blockName)
        print('Info: probe: Block global index:', self._ind)
        print('Info: probe: distance probe-node:', self._dist)
        print('Info: probe: extracted vars:', self._fields)
        print('Info: probe: filecur:', self._filecur)
        print('Info: probe: icur:', self._icur)
        return None

    # Create the probe zone with buffer size
    # IN: _proc: proc of probe
    # IN: _fields: champ a extraire
    # OUT: _pZone: zone de stockage
    def createProbeZone(self):
        if Cmpi.rank != self._proc: return
        self._pZone = G.cart((0,0,0), (1,1,1), (self._bsize,1,1))
        self._pZone[0] = 'probe'
        C._initVars(self._pZone, '{CoordinateX}=-1.') # time sentinel
        # create vars in probe      
        for v in self._fields:
            C._initVars(self._pZone, '{nodes:%s}=0.'%v)
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
            C.convertPyTree2File(self._pZone, self._fileName)
            self._filecur = 0
        # load GC
        nodes = Distributed.readNodesFromPaths(self._fileName, ['CGNSTree/Base/probe/GridCoordinates'])
        px = Internal.getNodeFromName2(self._pZone, 'CoordinateX')
        px2 = Internal.getNodeFromName2(nodes[0], 'CoordinateX')
        px[1] = px2[1]
        px = px[1]
        a = px > -0.5
        self._icur = numpy.count_nonzero(a)
        # load FS
        nodes = Distributed.readNodesFromPaths(self._fileName, ['CGNSTree/Base/probe/FlowSolution'])
        cont = Internal.getNodeFromName2(self._pZone, 'FlowSolution')
        if cont is not None: cont[2] = nodes[0][2]
        print('Info: probe: filecur:', self._filecur)
        print('Info: probe: icur:', self._icur)
        return None

    # verifie la var list dans t, conserve les pointeurs d'acces
    # IN: _proc
    # IN: _blockName
    # IN: varList
    # OUT: _fields
    # OUT: _pFields
    def checkVariables(self, varList):
        if Cmpi.rank != self._proc: return
        if self._blockName is None: return
        block = Internal.getNodeFromName2(self._t, self._blockName)
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
            self._pfields.append(var)
        return None

    # trigger extraction on current t
    # IN: _proc: proc du bloc contenant la probe
    # IN: _pZone: zone de stockage
    # IN: _ind: index of probe (static)
    # IN: _pfields: numpys des champs de la zone de probe
    # IN: _fields: nom des champs a extraire
    def extract(self, time):
        """Extract fields from t."""
        if Cmpi.rank != self._proc: return None
        # set time in CoordinateX
        px = Internal.getNodeFromName2(self._pZone, 'CoordinateX')[1]
        px = px.ravel('k')
        px[self._icur] = time
        for c in range(len(self._fields)):
            v = self._fields[c]
            f = self._pfields[c][1].ravel('k')
            pf = Internal.getNodeFromName2(self._pZone, v)[1]
            pf = pf.ravel('k')
            pf[self._icur] = f[self._ind]
            #print(f[self._ind])
        self._icur += 1
        if self._icur >= self._bsize: self.flush()
        return None

    # flush containers of probe
    # IN: _proc: proc du bloc contenant la probe
    # IN: _pZone: zone de la probe
    # IN: _fileName: nom du fichier
    def flush(self):
        """Flush probe to file."""
        if Cmpi.rank != self._proc: return None
        # flush containers
        gc = Internal.getNodeFromName1(self._pZone, 'GridCoordinates')
        fc = Internal.getNodeFromName1(self._pZone, 'FlowSolution')
        if self._icur >= self._bsize: # because buffer is out
            gc = Internal.copyNode(gc)
            gc[0] = 'GridCoordinates#%d'%self._filecur
            fc = Internal.copyNode(fc)
            fc[0] = 'FlowSolution#%d'%self._filecur
            print('Info: probe: flush %d (full).'%self._filecur)
            paths = ['CGNSTree/Base/probe','CGNSTree/Base/probe']
            nodes = [gc,fc]
            Distributed.writeNodesFromPaths(self._fileName, paths, nodes, mode=0)
            self._filecur += 1
            C._initVars(self._pZone, '{CoordinateX}=-1.') # time sentinel
            for v in self._fields: C._initVars(self._pZone, '{%s}=0.'%v)
            self._icur = 0
        else: # explicit flush
            nodes = [gc,fc]
            paths = ['CGNSTree/Base/probe/GridCoordinates', 'CGNSTree/Base/probe/FlowSolution']
            Distributed.writeNodesFromPaths(self._fileName, paths, nodes, mode=1)
        return None

    # read all probe fields as a single zone
    # IN: _fileName: nom du fichier
    def read(self):
        """Reread all data from probe file."""
        tl = Cmpi.convertFile2SkeletonTree(self._fileName)
        # read time
        nodes = Internal.getNodesFromName(tl, 'GridCoordinates#*')
        # load GCs
        paths = ['CGNSTree/Base/probe/GridCoordinates']
        for n in nodes:
            paths.append('CGNSTree/Base/probe/%s'%n[0])
        nodes = Distributed.readNodesFromPaths(self._fileName, paths)
        
        px = Internal.getNodeFromName1(nodes[0], 'CoordinateX')[1]
        a = px > -0.5
        csize = numpy.count_nonzero(a)
        size = csize
        for n in nodes[1:]:
            pxn = Internal.getNodeFromName1(n, 'CoordinateX')
            size += pxn[1].size

        out = G.cart((0,0,0), (1,1,1), (size,1,1))
        px2 = Internal.getNodeFromName2(out, 'CoordinateX')[1]
        c = 0
        for n in nodes[1:]:
            pxn = Internal.getNodeFromName1(n, 'CoordinateX')[1]
            px2[c:c+self._bsize] = pxn[0:self._bsize]
            c += self._bsize
        px2[c:c+csize] = px[0:csize]
        
        # load FS
        nodes = Internal.getNodesFromName(tl, 'FlowSolution#*')
        paths = ['CGNSTree/Base/probe/FlowSolution']
        for n in nodes:
            paths.append('CGNSTree/Base/probe/%s'%n[0])
        nodes = Distributed.readNodesFromPaths(self._fileName, paths)
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
