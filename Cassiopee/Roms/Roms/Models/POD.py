# POD subclass
from . import Model
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Converter.Filter as Filter
import numpy

class POD(Model):
    def __init__(self, name, type):
        super().__init__(name, type)
        # POD base modes
        self.Phi = None
        # number of modes
        self.K = 0

    def save(self):
        """Save Phi base modes."""
        if Cmpi.rank > 0: return None
        t = C.newPyTree(['POD'])
        node = ["phi", self.Phi, [], 'phi_t']
        t[2][1][2].append(node)
        C.convertPyTree2File(t, self.fileName)
        return None

    # build Phi from matrix M keeping K modes with svd
    def buildSvd(self, A, W, K=-1):
        # build modes
        Phi, S, Vt = numpy.linalg.svd(A, full_matrices=False)
        # energy of each modes
        #energy = S**2 / numpy.sum(S**2)
        if K > 0: self.K = K
        else: self.K = Phi.shape[1]
        self.Phi = Phi[:, 0:self.K]

        # build snapshot coefficients
        coords = numpy.empty( (self.K), dtype=numpy.float64 )
        for i in A.shape[0]:
            m = A[i,:].ravel('k')
            for i in range(self.K):
                c0 = numpy.dot(self.Phi[:,i], m)
                coords[i] = c0
            node = ["%05d"%i, coords, [], 'coeffs_t']
            Filter.writeNodesFromPaths(self.fileName, 'CGNSTree/POD', node)

        return Phi, S, Vt
