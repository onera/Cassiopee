# POD subclass - numpy implementation
from . import POD
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Converter.Filter as Filter
import Compressor.PyTree as Compressor
import Converter.Internal as Internal
import numpy

class nPOD( POD.POD ):
    def __init__(self, name, type):
        super().__init__(name, type)

    # build Phi from matrix A keeping K modes with svd
    def buildPhi(self, K=-1):
        """Build phi modes with SVD."""
        if self.A is None:
            raise ValueError("build: A matrix is not present.")
        # build modes
        Phi, S, Vt = numpy.linalg.svd(self.A, full_matrices=False)
        # energy of each modes
        #energy = S**2 / numpy.sum(S**2)
        if K > 0: self.K = K
        else: self.K = Phi.shape[1]
        self.Phi = Phi[:, 0:self.K]
        self.Phi = numpy.ascontiguousarray(self.Phi)
        return Phi, S, Vt

    def savePhi(self):
        """Save Phi base modes."""
        if Cmpi.rank > 0: return None
        if self.Phi is None:
            raise ValueError("save: Phi modes are not built.")
        t = C.newPyTree(['POD'])
        node = ["phi", self.Phi, [], 'DataArray_t']
        Compressor._packNode(node)
        t[2][1][2].append(node)
        if self.db is not None:
            node = Internal.createNode('db', 'DataArray_t', value=self.db.dirName)
            t[2][1][2].append(node)
        if self.ref is not None:
            node = Internal.createNode('ref', 'DataArray_t', value=self.ref)
            t[2][1][2].append(node)
        if self.variables is not None:
            node = Internal.createNode('variables', 'UserDefined_t')
            for v in self.variables:
                Internal.createChild(node, v, 'DataArray_t', value=v)
            t[2][1][2].append(node)
        if self.parameters is not None:
            node = Internal.createNode('parameters', 'UserDefined_t')
            for v in self.parameters:
                Internal.createChild(node, v, 'DataArray_t', value=v)
            t[2][1][2].append(node)

        C.convertPyTree2File(t, self.fileName, format='bin_hdf')
        return None

    def loadPhi(self):
        """Load phi from mod file."""
        node = Filter.readNodesFromPaths(self.fileName, ['CGNSTree/POD/phi'])[0]
        Compressor._unpackNode(node)
        self.Phi = node[1]
        nodes = Filter.readNodesFromPaths(self.fileName, ['CGNSTree/POD/db', 'CGNSTree/POD/ref', 'CGNSTree/POD/variables', 'CGNSTree/POD/parameters'])
        self.dbDirName = Internal.getValue(nodes[0])
        self.ref = Internal.getValue(nodes[1])
        self.variables = []
        for i in nodes[2][2]:
            self.variables.append(Internal.getValue(i))
        self.parameters = []
        for i in nodes[3][2]:
            self.parameters.append(Internal.getValue(i))

    def buildCoords(self):
        """Build snapshots coords."""
        if self.A is None:
            raise ValueError("build: A matrix is not present.")
        if self.Phi is None:
            raise ValueError("build: Phi modes are not built.")
        ncols = self.A.shape[1]
        self.coords = numpy.empty( (ncols, self.K), dtype=numpy.float64 )
        coords = numpy.empty( (self.K), dtype=numpy.float64 )
        for j in range(ncols):
            m = self.A[:,j].ravel('k')
            for i in range(self.K):
                c0 = numpy.dot(self.Phi[:,i], m)
                coords[i] = c0
            self.coords[j,:] = coords[:]

    def saveCoords(self):
        """Save snapshot coords."""
        if self.coords is None:
            raise ValueError("build: coords is not present.")
        node = ["coords", self.coords, [], 'DataArray_t']
        Compressor._packNode(node)
        Filter.writeNodesFromPaths(self.fileName, 'CGNSTree/POD', node, format='bin_hdf')

    def loadCoords(self):
        """Load all snapshot coords."""
        node = Filter.readNodesFromPaths(self.fileName, ['CGNSTree/POD/coords'])[0]
        Compressor._unpackNode(node)
        self.coords = node[1]

    def savePoints(self):
        """Save snapshot points."""
        if self.P is None:
            raise ValueError("build: points is not present.")
        node = ["points", self.P, [], 'DataArray_t']
        Compressor._packNode(node)
        Filter.writeNodesFromPaths(self.fileName, 'CGNSTree/POD', node, format='bin_hdf')

    def loadPoints(self):
        """Load all snapshot points."""
        node = Filter.readNodesFromPaths(self.fileName, ['CGNSTree/POD/points'])[0]
        Compressor._unpackNode(node)
        self.P = node[1]

    def instantiate(self, coords):
        """Return field vector from coords."""
        m = self.Phi @ coords
        return m

    def buildTree(self, coords):
        """Rebuild tree from field vector."""
        if self.dbDirName is None:
            raise ValueError("instantiate: db is missing.")
        if self.Phi is None:
            raise ValueError("instantiate: Phi modes are missing.")
        m = self.Phi @ coords
        cgnsName = self.dbDirName+'/%s'%self.ref+'.cgns'
        tref = C.convertFile2PyTree(cgnsName)
        for v in self.variables:
            v = v.split(':')
            if len(v) == 2:
                loc = v[0]; v = v[1]
            else: loc = 'nodes'; v = v[0]
            for z in Internal.getZones(tref):
                npts = C.getNPts(z)
                ncells = C.getNCells(z)
                pt0 = 0
                if loc == 'nodes':
                    FS = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
                    if FS is None:
                        FS = Internal.newFlowSolution(Internal.__FlowSolutionNodes__, 'Vertex', parent=z)
                    b = numpy.zeros( (npts, 1), dtype=numpy.float64)
                    Internal.newDataArray(v, value=b, parent=FS)
                    b[:,0] = m[pt0:pt0+npts]
                    pt0 += npts
                else:
                    FC = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
                    if FC is None:
                        FC = Internal.newFlowSolution(Internal.__FlowSolutionCenters__, 'CellCenter', parent=z)
                    b = numpy.zeros( (ncells, 1), dtype=numpy.float64)
                    Internal.newDataArray(v, value=b, parent=FC)
                    b[:,0] = m[pt0:pt0+ncells]
                    pt0 += ncells
        return tref

    def fetchTree(self, point):
        """Rebuild pyTree for given point."""
        from scipy.interpolate import RBFInterpolator
        if self.P is None:
            raise ValueError("instantiate: snapshot P points are missing.")
        if self.coords is None:
            raise ValueError("instantiate: snapshot coords are missing.")
        xo = self.P
        yo = self.coords
        # Create the RBF interpolator
        rbfInterpolator = RBFInterpolator(xo, yo, kernel='linear')
        # Define the input point
        np = xo.shape[1]
        inputPoint = numpy.empty( (1,np), dtype=numpy.float64 )
        for c, i in enumerate(self.parameters):
            inputPoint[0,c] = point[i]
        # Evaluate the interpolator at the input point
        coords = rbfInterpolator(inputPoint)
        coords = coords.ravel('k')
        t = self.buildTree(coords)
        return t
