# Interface pour MPI

import Converter.Mpi as Cmpi
import Converter.Internal as Internal
from . import PyTree as D2

#==============================================================================
# redispatch
# IN: graph: graph 'proc'
#==============================================================================
def redispatch(t, graph=None, verbose=0):
    """Redistribute tree from graph."""
    tp = Internal.copyRef(t)
    _redispatch(tp, graph, verbose)
    return tp

def _redispatch(t, graph=None, verbose=0):
    """Redistribute tree from graph."""
    import copy

    if graph is None: graph = Cmpi.computeGraph(t, type='proc')

    nzonesMax = len(Internal.getZones(t))

    for i in range(1, Cmpi.size):
        graph_local = copy.deepcopy(graph)

        for rank in range(Cmpi.size): # create a local 2 by 2 connectivity graph to reduce memory consumption
            if rank not in graph_local: graph_local[rank] = {}
            else:
                opp_rank = (rank+i)%Cmpi.size
                if opp_rank in graph_local[rank]: graph_local[rank] = {opp_rank: graph_local[rank][opp_rank]}
                else: graph_local[rank] = {}

        Cmpi._addXZones(t, graph_local)

        nzonesMax = max(nzonesMax, len(Internal.getZones(t)))

        for z in Internal.getZones(t):
            tag = Internal.getNodeFromName1(z, 'XZone')
            if tag is None:
                proc = Cmpi.getProc(z)
                if proc == (Cmpi.rank+i)%Cmpi.size:
                    (p, c) = Internal.getParentOfNode(t, z)
                    del p[2][c]
            else: # remove node tag xzone
                (p, c) = Internal.getParentOfNode(z, tag)
                del p[2][c]

    if verbose > 0:
        nzones = len(Internal.getZones(t))
        Cmpi.barrier()
        print("Info: _redispatch: proc {} has {} zones (max #zones during operation {} (+{:2.2f}%))".format(Cmpi.rank, nzones, nzonesMax, (nzonesMax-nzones)/float(nzones)*100.))
        Cmpi.barrier()

    return None

def redispatch_old(t, graph=None):
    """Redistribute tree from graph."""
    tp = Internal.copyRef(t)
    _redispatch_old(tp, graph)
    return tp

def _redispatch_old(t, graph=None):
    """Redistribute tree from graph."""
    if graph is None:
        graph = Cmpi.computeGraph(t, type='proc')
    procs = D2.getProcDict(t)
    Cmpi._addXZones(t, graph)
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
    return None
