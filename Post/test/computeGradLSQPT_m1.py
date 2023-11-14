import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Mpi as Cmpi
import Converter.Internal as Internal
import Post.Mpi as Pmpi
import XCore.PyTree as X
import KCore.test as test

LOCAL = test.getLocal()

def f(x, y, z): return 3.*x + 2.*y + z
def g(x, y, z): return 4.*x + 3.*y + 2.*z

if Cmpi.rank == 0:
    t = G.cartHexa((0.,0.,0.), (1./10.,1./10.,1./10.), (11,11,2))
    t = C.convertArray2NGon(t)
    Internal._adaptNGon32NGon4(t)
    t = C.initVars(t, 'centers:f', f, ['centers:CoordinateX', 'centers:CoordinateY', 'centers:CoordinateZ'])
    t = C.initVars(t, 'centers:g', g, ['centers:CoordinateX', 'centers:CoordinateY', 'centers:CoordinateZ'])
    C.convertPyTree2File(t, LOCAL+'/LSQcase.cgns')

Cmpi.barrier()

t, RES = X.loadAndSplitNGon(LOCAL+'/LSQcase.cgns')
t = C.makeParentElements(t)
t = Pmpi.computeGradLSQ(t, ['g', 'f'])

if Cmpi.rank == 0:
    test.testT(t, 1)
