# The Model class
import Converter.PyTree as C
import Converter.Mpi as Cmpi

class Model:
    def __init__(self, name, type):
        # model name
        if name[-4:] == ".mod": name = name[:-4]
        if name[-5:] == ".mod/": name = name[:-5]
        self.name = name
        # model type
        self.type = type
        # model file
        self.fileName = name+'.mod'
        # pytree link
        self.db = None # link to db
        self.dbDirName = None # db dir name
        self.ref = None # link to ref name in db
        self.parameters = None # parameter names
        self.variables = None # loadded variables in A
        self.filter = None # sub indices in global matrix
        # snapshot matrix
        self.A = None # distributed
        self.W = None # distributed
        # snapshot points
        self.P = None

    def set(self, A, W=None, P=None, db=None, ref=None, parameters=[], variables=[]):
        """Connect Model to matrix A, W or db, ref, variables."""
        self.A = A # snapshot matrix or submatrix
        self.W = W # scalar product weights
        self.P = P # snapshot points
        self.db = db # data base
        self.dbDirName = db.dirName # db dir name
        self.ref = ref # reference for A
        self.parameters = parameters # parameter names
        self.variables = variables # loaded variables in A

    def load(self):
        return None

    def save(self):
        return None

    def build(self):
        return None

    def instantiate(self):
        return None