# POD subclass
from . import Model

class POD( Model.Model ):
    def __init__(self, name, type):
        super().__init__(name, type)
        # POD base modes
        self.Phi = None
        # number of modes
        self.K = 0
        # snapshot coords on Phi
        self.coords = None