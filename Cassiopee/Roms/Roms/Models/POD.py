# POD subclass
from . import Model
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Converter.Filter as Filter
import Compressor.PyTree as Compressor
import numpy

class POD( Model.Model ):
    def __init__(self, name, type):
        super().__init__(name, type)
        # POD base modes
        self.Phi = None
        # number of modes
        self.K = 0

    def savePhi(self):
        """Save Phi base modes."""
        if Cmpi.rank > 0: return None
        t = C.newPyTree(['POD'])
        node = ["phi", self.Phi, [], 'DataArray_t']
        Compressor._packNode(node)
        t[2][1][2].append(node)
        C.convertPyTree2File(t, self.fileName)
        return None

    # build Phi from matrix M keeping K modes with svd
    def buildPhi(self, A, W=None, K=-1):
        """Build phi modes with SVD."""
        # build modes
        Phi, S, Vt = numpy.linalg.svd(A, full_matrices=False)
        # energy of each modes
        #energy = S**2 / numpy.sum(S**2)
        if K > 0: self.K = K
        else: self.K = Phi.shape[1]
        self.Phi = Phi[:, 0:self.K]
        return Phi, S, Vt

    def buildAndSaveCoeffs(self, A):
        """Build snapshots coefficients."""
        coords = numpy.empty( (self.K), dtype=numpy.float64 )
        for j in range(A.shape[1]):
            m = A[:,j].ravel('k')
            for i in range(self.K):
                c0 = numpy.dot(self.Phi[:,i], m)
                coords[i] = c0
            node = ["%05d"%j, coords, [], 'DataArray_t']
            #Compressor._packNode(node)
            #print(node, flush=True)
            Filter.writeNodesFromPaths(self.fileName, 'CGNSTree/POD', [node], format='bin_hdf')

    def instantiate(self):
        return None