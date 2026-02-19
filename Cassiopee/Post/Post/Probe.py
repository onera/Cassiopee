# bufferized parallel probe
import Converter.PyTree as C
import Geom.PyTree as D
import Converter.Mpi as Cmpi
import Converter.Internal as Internal
import Generator.PyTree as G
import Converter.Filter as Filter
import Converter.Distributed as Distributed
import Distributor2.PyTree as D2
import Connector.PyTree as X
import Connector.Mpi as Xmpi
import numpy, os

# Probe class
class Probe:
    """Create Probe object to extract and save data during simulations."""
    # defaullt initialization
    def init0(self):
        # probe mode:
        # mode=0, probe XYZ, single point
        # mode=1, probe ind, single point
        # mode=2, zones donnees
        # mode=3, zones interpolees
        # mode=4, global quantities
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

        # -- data probe mode=3 --
        self._ts = None
        self._graph = None
        self._procDicts = None

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
        self._append = False

        # zones storing probe data
        self._probeZones = None

    # init from position X/ind+blockName/None
    def __init__(self, fileName,
                 t=None,
                 X=None,
                 ind=None, blockName=None,
                 tPermeable=None,
                 fields=None, append=False,
                 bufferSize=100, writeCoords=True, modeForce=0):
        """Create a probe."""
        self.init0()
        self._bsize = bufferSize
        self._fileName = fileName
        self._coords = writeCoords
        self._append = append
        if fields is not None:
            self._fields = list(fields)
            loc = self.getFieldLoc(fields)
        else:
            self._fields = None
            loc = None

        # Localisation a partir de X (mode=0)
        if X is not None:
            self._mode = 0
            self._posX = X[0]
            self._posY = X[1]
            self._posZ = X[2]

            if t is not None and loc is not None:
                self.locateProbeXYZ(t, X, loc)
            else: raise ValueError("Probe: missing PyTree argument (t) for mode %s"%0)

        # Localisation a partir de ind,blockName (mode=1)
        elif ind is not None and blockName is not None:
            self._mode = 1
            if t is not None: self.locateProbeInd(t, ind, blockName)
            else: raise ValueError("Probe: missing PyTree argument (t) for mode %s"%1)

        elif tPermeable is not None:
            self._mode = 3
            self._ts = tPermeable

        elif t is None or modeForce == 4:
            self._mode = 4

        # Empilement de zones
        else:
            self._mode = 2

        # Cree la probe et on relit le fichier uniquement en mode=0 et 1
        if self._mode == 0 or self._mode == 1 or self._mode == 4:
            if self._proc is not None and self._fields is not None:
                self.createProbeZones()

            if self._proc is not None and self._probeZones is not None:
                self.checkFile(append=self._append)
            #if t is not None and fields is not None:
            #    self.checkVariables(t, fields)
        elif self._mode == 3:
            for z in Internal.getZones(tPermeable):
                dimz = Internal.getZoneDim(z)
                if dimz[0] != 'Structured': raise TypeError("Probe: tPermeable must be structured (mode %s)."%self._mode)
        else: # mode 2
            for z in Internal.getZones(t):
                dimz = Internal.getZoneDim(z)
                if dimz[0] != 'Structured': raise TypeError("Probe: t must be structured (mode %s)."%self._mode)

    # examine la localisation des champs
    # IN: fields
    def getFieldLoc(self, fields):
        loc = None
        for v in fields:
            vs = v.split(':')
            if len(vs) == 2 and vs[0] == 'centers':
                if loc is None: loc = 'centers'
                elif loc != 'centers': raise ValueError("Probe: fields must have the same loc.")
            else:
                if loc is None: loc = 'nodes'
                elif loc != 'nodes': raise ValueError("Probe: fields must have the same loc.")
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
                        ni = dim[1]; nj = dim[2]
                        ind = (ind[0]-1)+(ind[1]-1)*ni+(ind[2]-1)*ni*nj
                    else:
                        ni1 = max(dim[1]-1,1)
                        nj1 = max(dim[2]-1,1)
                        ind = (ind[0]-1)+(ind[1]-1)*ni1+(ind[2]-1)*ni1*nj1
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
        """Print information about a new probe."""
        if Cmpi.rank == 0:
            print('Info: probe is in mode %s'%self._mode)
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
        if self._mode == 0 or self._mode == 1 or self._mode == 4:
            if Cmpi.rank != self._proc and self._mode!=4: return
            z = G.cart((0,0,0), (1,1,1), (self._bsize,1,1))
            z[0] = 'probe'
            if self._mode!=4:D2._addProcNode(z, self._proc)
            self._probeZones = [z]
            C._initVars(self._probeZones, '{time}=-1.') # time sentinel
            # create vars in probe
            for c, v in enumerate(self._fields):
                v = v.split(':')
                if len(v) == 2: v = v[1]
                else: v = v[0]
                self._fields[c] = v # change location here because in mode=0, stored in nodes
                C._initVars(self._probeZones, '{nodes:%s}=0.'%v)

        elif self._mode == 2  and source is not None:
            zones = Internal.getZones(source)
            self._probeZones = []
            loc = self.getFieldLoc(self._fields)
            for z in zones:
                dimz = Internal.getZoneDim(z)
                ni = dimz[1]; nj = dimz[2]; nk = dimz[3]
                if loc == 'centers': # attention!
                    zp = G.cart((0,0,0), (1,1,1), (self._bsize+1,ni,nj))
                else: zp = G.cart((0,0,0), (1,1,1), (self._bsize,ni,nj))
                zp[0] = '%s'%z[0] # name of probe zone is identical to name of source
                D2._addProcNode(zp, Cmpi.rank)
                self._probeZones.append(zp)
                C._initVars(zp, 'time', -1) # time sentinel
                for v in self._fields: C._initVars(zp, v, 0.)

        elif self._mode == 3 and source is not None: # surface permeable
            zones = Internal.getZones(source)
            self._probeZones = []
            loc = self.getFieldLoc(self._fields)
            for z in zones:
                dimz = Internal.getZoneDim(z)
                if dimz[0] == 'Structured':
                    ni = dimz[1]; nj = dimz[2]; nk = dimz[3]
                    if loc == 'centers': # attention!
                        zp = G.cart((0,0,0), (1,1,1), (self._bsize+1,ni,nj))
                    else: zp = G.cart((0,0,0), (1,1,1), (self._bsize,ni,nj))
                else: # en non structure, c'est surement pas la bonne solution
                    npts = C.getNPts(z)
                    if loc == 'centers': # attention!
                        zp = G.cart((0,0,0), (1,1,1), (self._bsize+1,npts,1))
                    else: zp = G.cart((0,0,0), (1,1,1), (self._bsize,npts,1))
                zp[0] = '%s'%z[0] # name of probe zone is identical to name of source
                D2._addProcNode(zp, Cmpi.rank)
                self._probeZones.append(zp)
                C._initVars(zp, 'time', -1) # time sentinel
                for v in self._fields: C._initVars(zp, v, 0.)
        return None

    # Prepare for mode=3
    def prepare(self, t, loc='nodes', cartesian=False, extrap=1, nature=1, penalty=1, verbose=2):
        """Prepare the interpolation data for probe extraction in mode 3."""
        for z in Internal.getZones(t):
            dimz = Internal.getZoneDim(z)
            if dimz[0] != 'Structured': raise TypeError("Probe: donor tree (t) must be structured (mode %s)."%self._mode)

        tcs = Internal.copyRef(t)
        Xmpi._setInterpData2(self._ts, tcs, loc=loc, cartesian=cartesian, extrap=extrap, nature=nature, penalty=penalty, verbose=verbose)
        return tcs

    # Check file, if it doesnt exist, write probe zone in it
    # else get the filecur
    # IN: _proc: proc of probe
    def checkFile(self, append=False):
        if append:
            create = False
            if self._mode == 0 or self._mode == 1 or self._mode == 4:
                isFilePresent = os.path.exists(self._fileName)
            else:
                if Cmpi.rank == 0: isFilePresent = os.path.exists(self._fileName)
                else: isFilePresent = None
                isFilePresent = Cmpi.bcast(isFilePresent, root=0)

            if isFilePresent:
                if self._mode == 0 or self._mode == 1:
                    tl = Distributed.convertFile2SkeletonTree(self._fileName)
                elif self._mode == 4:
                    tl = C.convertFile2PyTree(self._fileName)
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
            elif self._mode == 4:
                C.convertPyTree2File(self._probeZones, self._fileName)
            else:
                Cmpi.convertPyTree2File(self._probeZones, self._fileName)

            self._filecur = 0
            return None

        #if self._mode != 4:
        #  Nettoyage + ne conserve que les zones du proc
        to = C.newPyTree()
        bases = Internal.getBases(tl)
        for b in bases:
            bl = Internal.newCGNSBase(b[0], parent=to)
            zones = Internal.getZones(b)
            for z in zones:
                Internal._rmNodesByName(z, 'GridCoordinates#*')
                Internal._rmNodesByName(z, 'FlowSolution#[0123456789]*')
                Internal._rmNodesByName(z, 'FlowSolution#Centers#*')
                if self._mode != 4:
                    if D2.getProc(z) == Cmpi.rank: bl[2].append(z)
                else:bl[2].append(z)
        tl = to

        self._probeZones = Internal.getZones(tl)


        # load GC
        zones = Internal.getZones(tl)
        paths = []
        for z in zones:
            if Internal.getNodeFromName1(z, 'GridCoordinates') is not None:
                paths += ['CGNSTree/Base/%s/GridCoordinates'%z[0]]
        if paths != []:
            nodes = Distributed.readNodesFromPaths(self._fileName, paths)
            for c, z in enumerate(self._probeZones):
                cont = Internal.getNodeFromName2(z, 'GridCoordinates')
                if cont is not None: cont[2] = nodes[c][2]

        # load FS
        paths = []
        for z in zones:
            if Internal.getNodeFromName1(z, 'FlowSolution') is not None:
                paths += ['CGNSTree/Base/%s/FlowSolution'%z[0]]
        if paths != []:
            nodes = Distributed.readNodesFromPaths(self._fileName, paths)
            for c, z in enumerate(self._probeZones):
                cont = Internal.getNodeFromName2(z, 'FlowSolution')
                if cont is not None: cont[2] = nodes[c][2]

        # load FS center
        paths = []
        for z in zones:
            if Internal.getNodeFromName1(z, 'FlowSolution#Centers') is not None:
                paths += ['CGNSTree/Base/%s/FlowSolution#Centers'%z[0]]
        if paths != []:
            nodes = Distributed.readNodesFromPaths(self._fileName, paths)
            for c, z in enumerate(self._probeZones):
                cont = Internal.getNodeFromName2(z, 'FlowSolution#Centers')
                if cont is not None: cont[2] = nodes[c][2]

        if len(self._probeZones) > 0:
            cont = Internal.getNodeFromName2(self._probeZones[0], 'FlowSolution')
            if cont is not None:
                pt = Internal.getNodeFromName2(cont, 'time')[1]
                if pt.ndim == 2: pt = pt[:,0]
                elif pt.ndim == 3: pt = pt[:,0,0]
                elif pt.ndim == 1: pt = pt[:]
                else: pt = pt.ravel('k')[:]
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
                raise ValueError("Probe: cannot find solution container in t.")
            var = Internal.getNodeFromName1(cont, v)
            if var is None:
                raise ValueError("Probe: cannot find field %s in t."%v)
            self._fields.append(v)
        return None

    # trigger extraction on current t
    # IN: _proc: proc du bloc contenant la probe
    # IN: _probeZones: zone de stockage
    # IN: _ind: index of probe (static)
    # IN: _fields: nom des champs a extraire
    # si onlyTransfer=True, les champs ne sont pas stockes dans la probe
    def extract(self, t=None, time=0, value=[], onlyTransfer=False):
        """Extract probe information."""

        if t is None and self._mode != 4:
            raise ValueError("Probe: missing PyTree argument (t) for extract with mode %s"%self._mode)

        if self._mode == 0 or self._mode == 1: # single point
            self.extract1(t, time)

        elif self._mode == 2: # single zone
            self.extract2(t, time)

        elif self._mode == 3:
            # attention ici : t is tcs
            self.extract3(t, time, onlyTransfer)

        elif self._mode == 4: # store values that are given
            self.extract4(time, value)
        return None

    def extract1(self, t, time):
        """Extract for mode=0 or 1."""
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
        return None

    def extract2(self, t, time):
        """Extract for mode=2"""
        if self._probeZones is None:
            self.createProbeZones(t)
            self.checkFile(append=self._append)
        Cmpi.barrier()

        # time is in "time" of probe zones
        source = Internal.getZones(t)
        for c, pz in enumerate(self._probeZones):
            sourcez = source[c]
            pt = Internal.getNodeFromName2(pz, 'time')[1]
            if pt.ndim == 2: pt[self._icur,:] = time
            elif pt.ndim == 1: pt[self._icur] = time
            else: pt[self._icur,:,:] = time
            # Set zone coordinates
            if self._coords:
                ptx = Internal.getNodeFromName2(pz, 'CoordinateX')[1]
                pty = Internal.getNodeFromName2(pz, 'CoordinateY')[1]
                ptz = Internal.getNodeFromName2(pz, 'CoordinateZ')[1]
                ptx2 = Internal.getNodeFromName2(sourcez, 'CoordinateX')[1]
                pty2 = Internal.getNodeFromName2(sourcez, 'CoordinateY')[1]
                ptz2 = Internal.getNodeFromName2(sourcez, 'CoordinateZ')[1]
                if ptx.ndim == 2:
                    if ptx2.size == ptx.shape[1]:
                        ptx[self._icur,:] = ptx2[:]
                        pty[self._icur,:] = pty2[:]
                        ptz[self._icur,:] = ptz2[:]
                    else: # champ en centres
                        #sourcezp = C.node2Center(sourcez)
                        ptx2 = Internal.getNodeFromName2(sourcez, 'CoordinateX')[1]
                        pty2 = Internal.getNodeFromName2(sourcez, 'CoordinateY')[1]
                        ptz2 = Internal.getNodeFromName2(sourcez, 'CoordinateZ')[1]
                        ptx[self._icur,:] = ptx2[:]
                        pty[self._icur,:] = pty2[:]
                        ptz[self._icur,:] = ptz2[:]
                else:
                    if ptx2.size == ptx.shape[1]*ptx.shape[2]:
                        ptx[self._icur,:,:] = ptx2[:,:]
                        pty[self._icur,:,:] = pty2[:,:]
                        ptz[self._icur,:,:] = ptz2[:,:]
                    else: # champ en centres
                        #sourcezp = C.node2Center(sourcez)
                        ptx2 = Internal.getNodeFromName2(sourcez, 'CoordinateX')[1]
                        pty2 = Internal.getNodeFromName2(sourcez, 'CoordinateY')[1]
                        ptz2 = Internal.getNodeFromName2(sourcez, 'CoordinateZ')[1]
                        ptx[self._icur,:,:] = ptx2[:,:]
                        pty[self._icur,:,:] = pty2[:,:]
                        ptz[self._icur,:,:] = ptz2[:,:]
            # Set zone fields
            for c, v in enumerate(self._fields):
                v = v.split(':')
                if len(v) == 2: v = v[1]
                else: v = v[0]
                f = Internal.getNodeFromName2(sourcez, v)[1]
                pf = Internal.getNodeFromName2(pz, v)[1]
                if pf.ndim == 2: pf[self._icur,:] = f[:]
                else: pf[self._icur,:,:] = f[:,:]

        self._icur += 1
        if self._icur >= self._bsize: self.flush()
        return None

    def extract3(self, tcs, time, onlyTransfer=False):
        """Extract for mode=3."""
        if self._probeZones is None:
            if onlyTransfer:
                self._probeZones = [] # create empty probeZones if onlyTransfer=True
                # don't write probe file if onlyTransfer=True
            else:
                self.createProbeZones(self._ts)
                self.checkFile(append=self._append)
            for v in self._fields: C._initVars(self._ts, v, 0.)

        Cmpi.barrier()
        # ts: permeable surface
        # tcs: tc for surface interp must be set
        if self._procDicts is None:
            tsBB = Cmpi.createBBoxTree(self._ts)
            self._procDicts = Cmpi.getProcDict(tsBB)
            if Cmpi.size == 1:
                for i in self._procDicts:
                    if self._procDicts[i] < 0: self._procDicts[i] = 0

        Cmpi.barrier()
        if self._graph is None:
            tcsBB = Cmpi.createBBoxTree(tcs)
            procDictc = Cmpi.getProcDict(tcsBB)
            #interDicts = X.getIntersectingDomains(tsBB, tcsBB)
            interDictD2R = X.getIntersectingDomains(tcsBB, tsBB)
            self._graph = Cmpi.computeGraph(tcsBB, type='bbox3', intersectionsDict=interDictD2R,
                                            procDict=procDictc, procDict2=self._procDicts, t2=tsBB, reduction=True)

        variables = [var.split(':')[-1] for var in self._fields] # self._fields can be located at centers
        Xmpi._setInterpTransfers(self._ts, tcs, variables=variables, graph=self._graph, procDict=self._procDicts)
        if onlyTransfer: return None
        #Cmpi.convertPyTree2File(self._ts, 'out.cgns')
        #Internal.printTree(self._probeZones)

        # Set time
        for z in self._probeZones:
            pf = Internal.getNodeFromName2(z, 'time')[1]
            if pf.ndim == 2: pf[self._icur, :] = time
            else: pf[self._icur,:,:] = time

        # Set Coordinates
        for z in Internal.getZones(self._ts):
            for pz in self._probeZones:
                if pz[0] == z[0]: break
            for n in ['CoordinateX', 'CoordinateY', 'CoordinateZ']:
                px = Internal.getNodeFromName2(z, n)[1]
                px2 = Internal.getNodeFromName2(pz, n)[1]
                if px2.ndim == 2: px2[self._icur, :] = px[:]
                else: px2[self._icur,:,:] = px[:,:]

        # Set zone fields
        for c, v in enumerate(self._fields):
            v = v.split(':')
            if len(v) == 2: v = v[1]
            else: v = v[0]
            for z in Internal.getZones(self._ts):
                for pz in self._probeZones:
                    if pz[0] == z[0]: break
                f = Internal.getNodeFromName2(z, v)[1]
                pf = Internal.getNodeFromName2(pz, v)[1]
                if pf.ndim == 2: pf[self._icur,:] = f[:]
                else: pf[self._icur,:,:] = f[:,:]

        self._icur += 1
        if self._icur >= self._bsize: self.flush()

        return None

    def extract4(self, time, value):
        """Extract for mode=4."""
        # time is in "time" of probe zone
        pzone = self._probeZones[0]
        pt = Internal.getNodeFromName2(pzone, 'time')[1]
        pt[self._icur] = time
        for c,v in enumerate(self._fields):
            pf = Internal.getNodeFromName2(pzone, v)[1]
            pf[self._icur] = value[c]
        self._icur += 1
        if self._icur >= self._bsize: self.flush()
        return None

    def share(self, tcs, tc):
        """Share coordinates and fields between tcs and tc."""
        zones = Internal.getZones(tc)
        for z in zones:
            zp = Internal.getNodeFromName2(tcs, z[0])
            if zp is None: continue
            gc = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
            gcp = Internal.getNodeFromName1(zp, Internal.__GridCoordinates__)
            if gc is not None and gcp is not None:
                px = Internal.getNodeFromName1(gc, 'CoordinateX')
                pxp = Internal.getNodeFromName1(gcp, 'CoordinateX')
                pxp[1] = px[1] # pointers
                py = Internal.getNodeFromName1(gc, 'CoordinateY')
                pyp = Internal.getNodeFromName1(gcp, 'CoordinateY')
                pyp[1] = py[1] # pointers
                pz = Internal.getNodeFromName1(gc, 'CoordinateZ')
                pzp = Internal.getNodeFromName1(gcp, 'CoordinateZ')
                pzp[1] = pz[1] # pointers
            fc = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
            fcp = Internal.getNodeFromName1(zp, Internal.__FlowSolutionNodes__)
            if fc is not None:
                if fcp is None: fcp = Internal.newFlowSolution(Internal.__FlowSolutionNodes__, gridLocation='Vertex', parent=zp)
                pfs = Internal.getNodesFromType(fc, 'DataArray_t')
                for pf in pfs:
                    pfp = Internal.getNodeFromName1(fcp, pf[0])
                    if pfp is not None: pfp[1] = pf[1]
                    else: Internal.newDataArray(pf[0], value=pf[1], parent=fcp)
            fc = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
            fcp = Internal.getNodeFromName1(zp, Internal.__FlowSolutionCenters__)
            if fc is not None:
                if fcp is None: fcp = Internal.newFlowSolution(Internal.__FlowSolutionCenters__, gridLocation='Centers', parent=zp)
                pfs = Internal.getNodesFromType(fc, 'DataArray_t')
                for pf in pfs:
                    pfp = Internal.getNodeFromName1(fcp, pf[0])
                    if pfp is not None: pfp[1] = pf[1]
                    else: Internal.newDataArray(pf[0], value=pf[1], parent=fcp)
        return None

    # flush containers of probe
    # IN: _proc: proc du bloc contenant la probe
    # IN: _probeZones: zone de la probe
    # IN: _fileName: nom du fichier
    def flush(self):
        """Flush buffered data to probe file."""
        if self._mode == 0 or self._mode == 1:
            if Cmpi.rank != self._proc: return None
            print('Info: probe: flush #%d [%s].'%(self._filecur, self._fileName))
            self.flush__()
        elif self._mode == 4:
            print('Info: probe: flush #%d [%s].'%(self._filecur, self._fileName))
            self.flush__()
        else:
            if Cmpi.rank == 0: print('Info: probe: flush #%d [%s].'%(self._filecur, self._fileName))
            Cmpi.seq(self.flush__)

    def flush__(self):
        if self._icur >= self._bsize: # because buffer is out
            for pzone in self._probeZones:
                nodes = []; paths = []
                gc = Internal.getNodeFromName1(pzone, Internal.__GridCoordinates__)
                if gc is not None:
                    gc = Internal.copyNode(gc)
                    gc[0] = 'GridCoordinates#%d'%self._filecur
                    nodes += [gc]
                    paths += ['CGNSTree/Base/%s'%pzone[0]]
                fc = Internal.getNodeFromName1(pzone, Internal.__FlowSolutionNodes__)
                if fc is not None:
                    fc = Internal.copyNode(fc)
                    fc[0] = 'FlowSolution#%d'%self._filecur
                    nodes += [fc]
                    paths += ['CGNSTree/Base/%s'%pzone[0]]
                fcc = Internal.getNodeFromName1(pzone, Internal.__FlowSolutionCenters__)
                if fcc is not None:
                    fcc = Internal.copyNode(fcc)
                    fcc[0] = 'FlowSolution#Centers#%d'%self._filecur
                    nodes += [fcc]
                    paths += ['CGNSTree/Base/%s'%pzone[0]]
                Distributed.writeNodesFromPaths(self._fileName, paths, nodes, None, -1, 0)
                C._initVars(pzone, 'time', -1) # time sentinel
                for v in self._fields: C._initVars(pzone, v, 0.)
            #print('self._filecur=',self._filecur)
            self._filecur += 1
            self._icur = 0

        else: # explicit flush
            for pzone in self._probeZones:
                nodes = []; paths = []
                gc = Internal.getNodeFromName1(pzone, 'GridCoordinates')
                if gc is not None:
                    nodes += [gc]
                    paths += ['CGNSTree/Base/%s/GridCoordinates'%pzone[0]]
                fc = Internal.getNodeFromName1(pzone, 'FlowSolution')
                if fc is not None:
                    nodes += [fc]
                    paths += ['CGNSTree/Base/%s/FlowSolution'%pzone[0]]
                fcc = Internal.getNodeFromName1(pzone, 'FlowSolution#Centers')
                if fcc is not None:
                    nodes += [fcc]
                    paths += ['CGNSTree/Base/%s/FlowSolution#Centers'%pzone[0]]
                Distributed.writeNodesFromPaths(self._fileName, paths, nodes, None, -1, 1)

        return None

    # read all probe times as a single zone
    # IN: _fileName: nom du fichier
    # IN: cont: container a extraire -> all points
    # IN: index: point a extraire -> all times
    def read(self, cont=None, ind=None, probeName=None):
        """Read data from existing probe file."""
        if cont is not None:
            if ind is not None:
                print('Warning: Probe (read): both cont and ind arguments were provided. Only cont argument is considered.')
            return self.readCont(cont)
        elif ind is not None:
            if probeName is None:
                print('Warning: Probe (read): no probeName argument was provided with index value(s). Reading the first probe zone by default.')
                probeName = 0
            return self.readInd(ind, probeName)
        else:
            print('Warning: Probe (read): no ind or cont arguments provided. Reading the first point of the first probe zone by default.')
            return self.readInd(0, 0)

    # Retourne tous les indices de tous les blocs d'un seul instant
    def readCont(self, cont):
        # Read the given cont from file
        tl = Cmpi.convertFile2SkeletonTree(self._fileName)
        zones = Internal.getZones(tl)

        out = []

        # Load full given container par zone
        for z in zones:
            paths = []
            isGC = False; isFS = False; isFC = False

            if Internal.getNodeFromName1(z, 'GridCoordinates#0') is not None:
                isGC = True
                paths += ['CGNSTree/Base/%s/GridCoordinates#%d'%(z[0], cont)]
            if Internal.getNodeFromName1(z, 'FlowSolution#0') is not None:
                isFS = True
                paths += ['CGNSTree/Base/%s/FlowSolution#%d'%(z[0], cont)]
            if Internal.getNodeFromName1(z, 'FlowSolution#Centers#0') is not None:
                isFC = True
                paths += ['CGNSTree/Base/%s/FlowSolution#Centers#%d'%(z[0], cont)]
            nodes2 = Distributed.readNodesFromPaths(self._fileName, paths)

            dimz = Internal.getZoneDim(z)
            nrec = dimz[1]; ni = dimz[2]; nj = dimz[3]
            if isFC: nrec = nrec-1 # if data is stored at the centers
            if nj == 1: zsize = [[ni,ni-1,0]]
            else: zsize = [[ni,ni-1,0],[nj,nj-1,0]]
            for nr in range(nrec):
                cz = Internal.newZone(z[0]+'@'+str(nr), zsize=zsize, ztype='Structured')
                if isGC:
                    gc = Internal.newGridCoordinates(parent=cz)
                    px = Internal.getNodeFromName1(nodes2[0], 'CoordinateX')[1]
                    py = Internal.getNodeFromName1(nodes2[0], 'CoordinateY')[1]
                    pz = Internal.getNodeFromName1(nodes2[0], 'CoordinateZ')[1]
                    if nj > 1:
                        px = px[nr,:]; py = py[nr,:]; pz = pz[nr,:]
                    else:
                        px = px[nr]; py = py[nr]; pz = pz[nr]
                    ox = Internal.newDataArray('CoordinateX', value=px, parent=gc)
                    oy = Internal.newDataArray('CoordinateY', value=py, parent=gc)
                    oz = Internal.newDataArray('CoordinateZ', value=pz, parent=gc)

                if isFS:
                    pos = 0
                    if isGC: pos = 1
                    fs = Internal.newFlowSolution('FlowSolution', 'Vertex', parent=cz)
                    for fields in Internal.getNodesFromType(nodes2[pos], 'DataArray_t'):
                        pf = fields[1]
                        if nj > 1: pf = pf[nr,:]
                        else: pf = pf[nr]
                        op = Internal.newDataArray(fields[0], value=pf, parent=fs)

                if isFC:
                    pos = 0
                    if isGC: pos = 1
                    if isFC: pos += 1
                    fs = Internal.newFlowSolution('FlowSolution#Centers', 'CellCenter', parent=cz)
                    for fields in Internal.getNodesFromType(nodes2[pos], 'DataArray_t'):
                        pf = fields[1]
                        if nr < pf.shape[0]:
                            if nj > 1: pf = pf[nr,:]
                            else: pf = pf[nr]
                            op = Internal.newDataArray(fields[0], value=pf, parent=fs)

                out.append(cz)
        return out

    # Retourne tous les temps d'un seul indice d'un bloc
    def readInd(self, ind, probeName):
        tl = Cmpi.convertFile2SkeletonTree(self._fileName)

        # Recupere la probe a recuperer
        zones = Internal.getZones(tl)
        if isinstance(probeName, int): pz = zones[probeName]
        else: pz = Internal.getNodeFromName2(tl, probeName)
        dim = Internal.getZoneDim(pz)
        sizeNPts = dim[2]

        if isinstance(ind, int) or isinstance(ind, tuple): ind = [ind]
        sizeNPts = min(len(ind), sizeNPts)

        # recupere le sizeTime (nbre d'instants dans le fichier)
        isFC = False
        if Internal.getNodeFromName1(pz, 'FlowSolution#Centers#0') is not None: isFC = True
        ncont = 0 # nbre de containers pleins
        while Internal.getNodeFromName1(pz, 'FlowSolution#%s'%ncont) is not None: ncont += 1
        paths = ['CGNSTree/Base/%s/%s'%(pz[0], Internal.__FlowSolutionNodes__)]
        nodesTime = Distributed.readNodesFromPaths(self._fileName, paths)
        pt = Internal.getNodeFromName1(nodesTime[0], 'time')[1]
        if pt.ndim == 2: pt = pt[:,0]
        elif pt.ndim == 3: pt = pt[:,0,0]
        csize = numpy.count_nonzero(pt > -0.5)
        sizeTimeCont = pt.shape[0]
        if isFC: sizeTimeCont -= 1 # if data is stored at the centers
        sizeTime = csize + ncont*sizeTimeCont # +1 to take into account template cont
        sizeArray = (sizeTime) if sizeNPts == 1 else (sizeTime,sizeNPts)
        # print('sizeTime', sizeTime)
        # print('sizeNPts', sizeNPts)
        # print('sizeTimeCont', sizeTimeCont)
        # print('ncont', ncont)
        # print('csize', csize)

        # Rebuild full 1D zone (time only)
        out = G.cart((0,0,0), (1,1,1), (sizeTime,sizeNPts,1))
        out[0] = 'probe'
        ptx = Internal.getNodeFromName2(out, 'CoordinateX')[1]
        pty = Internal.getNodeFromName2(out, 'CoordinateY')[1]
        ptz = Internal.getNodeFromName2(out, 'CoordinateZ')[1]
        fcont = Internal.createUniqueChild(out, 'FlowSolution', 'FlowSolution_t')

        # load GCs
        cur = 0
        nodes = [Internal.getNodeFromName(pz, 'GridCoordinates#%s'%i) for i in range(ncont)]
        for n in nodes:
            paths = ['CGNSTree/Base/%s/%s'%(pz[0],n[0])]
            nodesX = Distributed.readNodesFromPaths(self._fileName, paths)[0]
            ptx2 = Internal.getNodeFromName2(nodesX, 'CoordinateX')[1]
            pty2 = Internal.getNodeFromName2(nodesX, 'CoordinateY')[1]
            ptz2 = Internal.getNodeFromName2(nodesX, 'CoordinateZ')[1]
            self.slice__(ptx, ptx2, cur, sizeTimeCont, ind, dim[2])
            self.slice__(pty, pty2, cur, sizeTimeCont, ind, dim[2])
            self.slice__(ptz, ptz2, cur, sizeTimeCont, ind, dim[2])
            cur += sizeTimeCont

        paths = ['CGNSTree/Base/%s/%s'%(pz[0], Internal.__GridCoordinates__)]
        nodesX = Distributed.readNodesFromPaths(self._fileName, paths)[0]
        ptx2 = Internal.getNodeFromName2(nodesX, 'CoordinateX')[1]
        pty2 = Internal.getNodeFromName2(nodesX, 'CoordinateY')[1]
        ptz2 = Internal.getNodeFromName2(nodesX, 'CoordinateZ')[1]
        self.slice__(ptx, ptx2, cur, csize, ind, dim[2])
        self.slice__(pty, pty2, cur, csize, ind, dim[2])
        self.slice__(ptz, ptz2, cur, csize, ind, dim[2])

        # load FSs
        cur = 0
        nodes = [Internal.getNodeFromName(pz, 'FlowSolution#%s'%i) for i in range(ncont)]
        for n in nodes:
            paths = ['CGNSTree/Base/%s/%s'%(pz[0],n[0])]
            nodesF = Distributed.readNodesFromPaths(self._fileName, paths)[0]
            for i in nodesF[2]:
                if i[3] == 'DataArray_t':
                    varname, ptf2 = i[0], i[1]
                    pf = Internal.getNodeFromName2(fcont, varname)
                    if pf is None: pf = Internal.newDataArray(varname, value=numpy.zeros(sizeArray, dtype=numpy.float64, order='F'), parent=fcont)
                    ptf = pf[1]
                    self.slice__(ptf, ptf2, cur, sizeTimeCont, ind, dim[2])
            cur += sizeTimeCont

        paths = ['CGNSTree/Base/%s/%s'%(pz[0], Internal.__FlowSolutionNodes__)]
        nodesF = Distributed.readNodesFromPaths(self._fileName, paths)[0]
        for i in nodesF[2]:
            if i[3] == 'DataArray_t':
                varname, ptf2 = i[0], i[1]
                ptf = Internal.getNodeFromName2(fcont, varname)[1]
                self.slice__(ptf, ptf2, cur, csize, ind, dim[2])

        # load FCs (if any)
        if isFC:
            cur = 0
            nodes = [Internal.getNodeFromName(pz, 'FlowSolution#Centers#%s'%i) for i in range(ncont)]
            for n in nodes:
                paths = ['CGNSTree/Base/%s/%s'%(pz[0],n[0])]
                nodesF = Distributed.readNodesFromPaths(self._fileName, paths)[0]
                for i in nodesF[2]:
                    if i[3] == 'DataArray_t':
                        varname, ptf2 = i[0], i[1]
                        pf = Internal.getNodeFromName2(fcont, varname)
                        if pf is None: pf = Internal.newDataArray(varname, value=numpy.zeros(sizeArray, dtype=numpy.float64, order='F'), parent=fcont)
                        ptf = pf[1]
                        self.slice__(ptf, ptf2, cur, sizeTimeCont, ind, dim[2])
                cur += sizeTimeCont

            paths = ['CGNSTree/Base/%s/%s'%(pz[0], Internal.__FlowSolutionCenters__)]
            nodesF = Distributed.readNodesFromPaths(self._fileName, paths)[0]
            for i in nodesF[2]:
                if i[3] == 'DataArray_t':
                    varname, ptf2 = i[0], i[1]
                    ptf = Internal.getNodeFromName2(fcont, varname)[1]
                    self.slice__(ptf, ptf2, cur, csize, ind, dim[2])

        return out

    # slice from ptx2 a list of indices
    def slice__(self, ptx, ptx2, start, size, ind, ni):
        if ptx2.ndim == 1:
            ptx[start:start+size] = ptx2[0:size]
        elif len(ind) == 1:
            i = ind[0]
            if ptx2.ndim == 2: ptx[start:start+size] = ptx2[0:size,i]
            else:
                if isinstance(i, tuple) and len(i) == 2: ii = i[0]; jj = i[1]
                else: jj = i//ni; ii = i-jj*ni
                ptx[start:start+size] = ptx2[0:size,ii,jj]
        else:
            for c, i in enumerate(ind):
                if ptx2.ndim == 2: ptx[start:start+size,c] = ptx2[0:size,i]
                else:
                    if isinstance(i, tuple) and len(i) == 2: ii = i[0]; jj = i[1]
                    else: jj = i//ni; ii = i-jj*ni
                    ptx[start:start+size,c] = ptx2[0:size,ii,jj]
        return None

    # mode 2: lit un container, retourne un numpy du champ pour
    # chaque probe zone (ncells, ntime), concatene avec le precedent
    def readCont2(self, cont, field, pin):
        # field name
        spl = field.split(':')
        if len(spl) == 2:
            fieldName = spl[1]
            if spl[0] == 'centers': fieldCont = 'FlowSolution#Centers'
            else: fcont = 'FlowSolution'
        else:
            fieldName = field; fieldCont = 'FlowSolution'

        # read container
        out = self.read(cont=cont)

        # concatenate all times in p
        zones = Internal.getZones(out)

        # get all probe zone names and number of time
        ntime = 0
        zoneNames = {}
        tmin = 1.e6; tmax = -1.e6
        for z in zones:
            name = z[0]
            name = name.split('@')
            if not name[0] in zoneNames: zoneNames[name[0]] = {}
            time = int(name[1])
            zoneNames[name[0]][time] = z
            tmin = min(time, tmin)
            tmax = max(time, tmax)

        ntime = tmax-tmin+1
        nzones = len(zoneNames.keys())
        #print('INFO: read from probe from time=%d to time=%d'%(tmin, tmax))
        #print(zoneNames.keys())

        allps = []
        for name in zoneNames: # all probe zones
            z = zoneNames[name][tmin]
            npts = C.getNPts(z)
            ncells = C.getNCells(z)
            p2 = numpy.zeros( (ncells, ntime), dtype=numpy.float64)
            for it in range(tmin, tmax+1):
                z = zoneNames[name][it]
                fss = Internal.getNodeFromName1(z, fieldCont)
                pressure = Internal.getNodeFromName1(fss, fieldName)[1]
                p2[:,it] = pressure.ravel('k')[:]
            allps.append(p2)

        # concatenate in input
        if pin == []: pin = [None]*nzones
        for c, r in enumerate(allps):
            if pin[c] is None: pin[c] = r
            else:
                ncells = pin[c].shape[0]
                ntime1 = pin[c].shape[1]
                ntime2 = r.shape[1]
                a = numpy.zeros( (ncells, ntime1+ntime2), dtype=numpy.float64)
                a[:,0:ntime1] = pin[c][:,0:ntime1]
                a[:,ntime1:ntime1+ntime2] = r[:,0:ntime2]
                pin[c] = a
        return pin
