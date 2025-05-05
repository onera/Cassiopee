"""Distribution module for Cassiopee package.
"""
__version__ = '4.1'
__author__ = "Christophe Benoit, Xavier Juvigny, Stephanie Peron, Pascal Raud"

from . import distributor2
from Converter.Internal import E_NpyInt
import Converter
import numpy

#==============================================================================
# - distribute -
# IN: arrays: les arrays a equilibrer ou le nbre de pts de chaque zone
# IN: NProc: le nombre de processeurs
# IN: prescribed: le tableaux des blocs dont le proc est impose
# prescribed[i] = 0 veut dire que le bloc i doit etre place sur le proc 0
# IN: perfo: performance de chaque processeur
# perfo[0] = (alpha, beta, gamma) avec alpha le ratio du solveur
# beta le ratio des communications par connection (latence) et gamma le ratio
# des communications par volume (comSpeed)
# IN: weight: poids relatif pour chaque bloc. Utile si le solveur n'est
# pas le meme sur tous les blocs
# IN: com: la matrice du volume de communication
# com[i,j] matrice NblocxNbloc indiquant le volume de com entre le bloc i
# et le bloc j
# IN: algorithm: 'gradient0', 'gradient1', 'genetic', 'fast'
# IN: nghost: nbre de couches de ghost cells
#==============================================================================
def distribute(arrays, NProc, prescribed=None, perfo=None, weight=None, com=None, comd=None,
               algorithm='graph', mode='nodes', nghost=0):
    """Distribute zones over NProc processors.
    Usage: distribute(A, NProc, prescribed, perfo, weight, com, algorithm)"""
    if NProc <= 0:
        raise ValueError("distribute: can not distribute on %d (<=0) processors."%NProc)

    # Liste du nombre de points pour chaque arrays
    # mode: equilibre le nbre de pts ou de cellules suivant 'nodes','cells'
    nbPts = []
    if isinstance(arrays[0], int): # le nbre de pts est deja dans arrays
        nbPts = arrays
    else: # sinon, on calcule nbPts a partir des arrays
        for a in arrays:
            if mode == 'cells': c = Converter.getNCells(a)
            else: c = Converter.getNPts(a)
            nbPts.append(c)

    # Liste des arrays deja distribues
    if prescribed is None: # nothing set
        setArrays = [-1]*len(arrays)
    else: setArrays = prescribed

    # Liste des alpha, beta, gamma pour chaque processeur
    if perfo is None:
        # Poids du solveur (par defaut)
        alpha = 1.
        # Poids de la latence (temps pour chaque com)
        beta = 1.e-2
        # Poids de la vitesse de com pour une unite de volume de com
        gamma = 0.1
        perfProcs = [(alpha,beta,gamma)]*NProc
    elif isinstance(perfo, tuple):
        perfProcs = [perfo]*NProc
    else:
        perfProcs = perfo

    # Liste des poids du solveur pour chaque bloc
    Nb = len(arrays)
    if weight is None: weight = [1.]*Nb

    # Matrice du volume des coms (volCom ou volComd)
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

    # Si algo=graph et pas de com, force algo=fast
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
