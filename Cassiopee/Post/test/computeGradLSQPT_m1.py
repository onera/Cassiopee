import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Mpi as Cmpi
import Converter.Internal as Internal
import Post.Mpi as Pmpi
import XCore.PyTree as X
import KCore.test as test
import os
import numpy

LOCAL = test.getLocal()

def f(x, y, z): return 3.*x + 2.*y + z
def g(x, y, z): return 4.*x + 3.*y + 2.*z

filepath = os.path.join(LOCAL, 'LSQcase.cgns')
if Cmpi.master:
    t = G.cartHexa((0.,0.,0.), (1./10.,1./10.,1./10.), (11,11,2))
    t = C.convertArray2NGon(t, api=3)
    C._initVars(
        t, 'centers:f', f,
        ['centers:CoordinateX', 'centers:CoordinateY', 'centers:CoordinateZ'],
        isVectorized=True
    )
    C._initVars(
        t, 'centers:g', g,
        ['centers:CoordinateX', 'centers:CoordinateY', 'centers:CoordinateZ'],
        isVectorized=True
    )
    C.convertPyTree2File(t, filepath)

Cmpi.barrier()

t, res = X.loadAndSplitNGon(filepath)
t = C.makeParentElements(t)
t = Pmpi.computeGradLSQ(t, ['g', 'f'])

# Partitioning may differ from that of the reference
nx = Internal.getNodesFromName(t, 'CoordinateX')[0]
minx = numpy.min(nx[1])
if minx == 0.5:
    z = Internal.getNodesFromName(t, 'cartHexa_1')
    if z: Internal._renameNode(t, 'cartHexa_1', 'cartHexa_0')
    zGC = Internal.getNodesFromName(t, 'ZoneGridConnectivity')
    z = Internal.getNodesFromName(zGC, 'Match_0')
    if z: Internal._renameNode(t, 'Match_0', 'Match_1')
    zProc = Internal.getNodesFromName(t, 'proc')
    zProc[0][1] = numpy.array([[0]])
    test.testT(t, 1)

if Cmpi.master: os.remove(filepath)
