# Interface pour MPI

import Converter.Mpi as Cmpi
import Converter.Internal as Internal
from . import PyTree as D2

#==============================================================================
# redispatch
# IN: graph: graph 'proc'
#==============================================================================
def redispatch(t, graph=None):
    """Redistribute tree from graph."""
    if graph is None:
        graph = Cmpi.computeGraph(t, type='proc')
    procs = D2.getProcDict(t)
    t = Cmpi.addXZones(t, graph)
    # Enleve les zones envoyees
    zones = Internal.getZones(t)
    for z in zones:
        tag = Internal.getNodeFromName1(z, 'XZone')
        if tag is None:
            if procs[z[0]] != Cmpi.rank:
                (p, c) = Internal.getParentOfNode(t, z)
                del p[2][c]
        else: # enleve le noeud tag XZone
            (p, c) = Internal.getParentOfNode(z, tag)
            del p[2][c]
    return t
