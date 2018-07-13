# Interface pour MPI

import Converter.Mpi as Cmpi
import Converter.Internal as Internal
import PyTree as D2

#==============================================================================
# redispatch
# IN: graph: graph 'proc'
#==============================================================================
def redispatch(t, graph=None):
    if graph is None:
        graph = Cmpi.computeGraph(t, type='proc')
    procs = D2.getProcDict(t)
    t = Cmpi.addXZones(t, graph)
    # Enleve les zones envoyees
    zones = Internal.getZones(t)
    for z in zones:
        tag = Internal.getNodesFromName1(z, 'XZone')
        if len(tag) == 0:
            if procs[z[0]] != Cmpi.rank:
                (p, c) = Internal.getParentOfNode(t, z)
                del p[2][c]
        else: # enleve le noeud tag XZone
            (p, c) = Internal.getParentOfNode(z, tag[0])
            del p[2][c]
    return t
