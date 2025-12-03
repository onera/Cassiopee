# General class of Cassiopee App
# Data oriented
# One method: run

__version__ = '4.1'
__author__ = "Stephanie Peron, Christophe Benoit, Thomas Renaud, Ivan Mary"

import Converter.PyTree as C

class App:
    """General applet for Cassiopee Apps."""
    def __init__(self):
        self.__version__ = "0.0"
        self.authors = ["mickey.mouse@whitehouse.com"]
        self.data = {}
        self.data['dataOut'] = 'out.cgns'
        self.required = []

    def set(self, **kwargs):
        """Set keywords in data dict."""
        for key, value in kwargs.items():
            self.data[key] = value

    def get(self, key):
        """Return a value stored in data dict."""
        return self.data[key]

    def readDataIn(self, name=None):
        """Return a pyTree depending on dataIn."""
        if name is None: name = 'dataIn'
        inp = self.data[name]
        if isinstance(inp, str): t = C.convertFile2PyTree(inp)
        else: t = inp
        return t

    def writeDataOut(self, t, name=None):
        """Export a pyTree to dataOut."""
        if name is None: name = 'dataOut'
        outp = self.data[name]
        if isinstance(outp, str): C.convertPyTree2File(t, outp)
        else: self.data[name] = t
        return None

    def run(self):
        """Run function."""
        # Cette fonction doit etre redefinie par l'applicatif
        print(self.data)
        return

    def requires(self, listKwrd):
        """Set required keywords."""
        self.required = listKwrd

    def required(self):
        """Return a list of required keywords."""
        return self.required

    def check(self):
        """Check if all keywords have been set."""
        for i in self.required:
            if i not in self.data:
                print('key: %s is not present.'%i)
                print('all required keys: %s.'%self.required)
                return False
        return True

    def viewCode(self):
        return
