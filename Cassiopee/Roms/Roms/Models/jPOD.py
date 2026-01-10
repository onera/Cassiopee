# POD subclass - jax implementation
from . import POD
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Converter.Filter as Filter
import Compressor.PyTree as Compressor
import numpy

class jPOD( POD.POD ):
    def __init__(self, name, type):
        super().__init__(name, type)
        # POD base modes
        self.Phi = None
        # number of modes
        self.K = 0

    def savePhi(self):
        """Save Phi base modes."""
        return None

    # build Phi from matrix M keeping K modes with svd
    def buildPhi(self, A, W=None, K=-1):
        """Build phi modes with SVD."""
        return None

    def buildAndSaveCoeffs(self, A):
        """Build snapshots coefficients."""
        return None

    def readCoeffs(self, j):
        return None

    def instantiate(self, coords):
        """Return snapshot from coords."""
        m = self.Phi @ coords
        return m

    def fetchTree(self, point, db, ref, variables):
        """Rebuild pyTree."""
        return None
