"""Distribution module for Cassiopee package.
"""
__version__ = '4.2'
__author__ = "Christophe Benoit, Xavier Juvigny, Stephanie Peron, Pascal Raud"

from . import distributor2
from Converter.Internal import E_NpyInt
import Converter
import numpy

#==============================================================================
# - distribute -
# IN: arrays: the arrays to balance or the number of points of each zone
# IN: NProc: the number of processors
# IN: prescribed: the array of blocks whose proc is prescribed
# prescribed[i] = 0 means that block i must be placed on proc 0
# IN: perfo: performance of each processor
# perfo[0] = (alpha, beta, gamma) with alpha the solver ratio
# beta the ratio of communications by connection (latency) and gamma the ratio
# of communications by volume (comSpeed)
# IN: weight: relative weight for each block. Useful if the solver is not
# the same on all blocks
# IN: com: the communication volume matrix
# com[i,j] NblocxNbloc matrix indicating the com volume between block i
# and block j
# IN: algorithm: 'gradient0', 'gradient1', 'genetic', 'fast'
# IN: nghost: number of ghost cell layers
#==============================================================================
def distribute(arrays, NProc, prescribed=None, perfo=None, weight=None, com=None, comd=None,
               algorithm='graph', mode='nodes', nghost=0):
    """Distribute zones over NProc processors.
    Usage: distribute(A, NProc, prescribed, perfo, weight, com, algorithm)"""
    if NProc <= 0:
        raise ValueError("distribute: can not distribute on %d (<=0) processors."%NProc)

    # List of the number of points for each array
    # mode: balances the number of pts or cells according to 'nodes','cells'
    nbPts = []
    if isinstance(arrays[0], int): # the number of pts is already in arrays
        nbPts = arrays
    else: # otherwise, we calculate nbPts from the arrays
        for a in arrays:
            if mode == 'cells': c = Converter.getNCells(a)
            else: c = Converter.getNPts(a)
            nbPts.append(c)

    # List of already distributed arrays
    if prescribed is None: # nothing set
        setArrays = [-1]*len(arrays)
    else: setArrays = prescribed

    # List of alpha, beta, gamma for each processor
    if perfo is None:
        # Solver weight (by default)
        alpha = 1.
        # Weight of latency (time for each com)
        beta = 1.e-2
        # Weight of com speed for one unit of com volume
        gamma = 0.1
        perfProcs = [(alpha,beta,gamma)]*NProc
    elif isinstance(perfo, tuple):
        perfProcs = [perfo]*NProc
    else:
        perfProcs = perfo

    # List of solver weights for each block
    Nb = len(arrays)
    if weight is None: weight = [1.]*Nb

    # Matrix of com volumes (volCom or volComd)
    volCom = None; volComd = None
    if com is None and comd is None: volComd = numpy.empty((0), dtype=E_NpyInt)
    elif com is not None and comd is None:
        if isinstance(com, list): volCom = numpy.array(com)
        else: volCom = com
    else: # comd is not None
        if isinstance(comd, dict):
            allkeys = comd.keys()
            size = len(allkeys)
            volComd = numpy.empty((2*size), dtype=E_NpyInt)
            for i, k in enumerate(allkeys):
                volComd[2*i] = k
                volComd[2*i+1] = comd[k]
        else: volComd = comd

    # If algo=graph and no com, force algo=fast
    if volCom is not None:
        if algorithm == 'graph' and numpy.amax(volCom) <= 0: algorithm = 'fast'
    if volComd is not None:
        if algorithm == 'graph' and volComd.size == 0: algorithm = 'fast'
    if volCom is None and volComd is None:
        if algorithm == 'graph': algorithm = 'fast'
    #print('algorithm', algorithm)
    #print('com:\n', com)
    #print('comd', comd)
    #print('volCom:\n', volCom)
    #print('volComd', volComd)

    # Distribution
    out = distributor2.distribute(nbPts, setArrays, perfProcs, weight,
                                  volCom, volComd, NProc, algorithm)
    return out
