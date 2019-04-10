# Interface pour MPI

import Converter.Mpi as Cmpi
import PyTree as P
import Converter.Internal as Internal

#==============================================================================
# extractMesh
# IN: t: maillage source distribue
# IN: extractionMesh: maillage de destination distribue
# IN: graph: graph d'intersection si deja calcule
#==============================================================================
def extractMesh(t, extractionMesh, order=2, extrapOrder=1,
                constraint=40., tol=1.e-6, mode='robust', hook=None, graph=None):
    if graph is None:
        tb = Cmpi.createBBoxTree(t)
        tb2 = Cmpi.createBBoxTree(extractionMesh)
        graph = Cmpi.computeGraph(tb, type='bbox3', t2=tb2)
    tl = Cmpi.addXZones(t, graph)
    tl = Cmpi.convert2PartialTree(tl)
    ext = Cmpi.convert2PartialTree(extractionMesh)
    # print info
    zones = Internal.getZones(tl)
    print ('Rank %d has %d source zones.'%(Cmpi.rank, len(zones)))
    ext = P.extractMesh(tl, ext, order=order, extrapOrder=extrapOrder, constraint=constraint, tol=tol, mode=mode,
                        hook=hook)
    return ext
