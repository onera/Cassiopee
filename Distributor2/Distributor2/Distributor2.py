"""Distribution module for Cassiopee package.
"""
__version__ = '2.9'
__author__ = "Christophe Benoit, Xavier Juvigny, Stephanie Peron, Pascal Raud"
# 
# Python Interface for block distribution over processors
#
from . import distributor2
import numpy

#==============================================================================
# - distribute -
# IN: arrays: les arrays a equilibrer
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
def distribute(arrays, NProc, prescribed=[], perfo=[], weight=[], com=[],
               algorithm='graph', mode='nodes', nghost=0):
    """Distribute zones over NProc processors.
    Usage: distribute(A, NProc, prescribed, perfo, weight, com, algorithm)"""
    if NProc <= 0: 
        raise ValueError("distribute: can not distribute on %d (<=0) processors."%NProc)

    # Liste du nombre de points pour chaque arrays
    # mode: equilibre le nbre de pts ou de cellules suivant 'nodes','cells'
    nbPts = []
    for a in arrays:
        c = 0
        if len(a) == 5:
            if mode == 'cells':
                c = max(a[2]-1-nghost,1)*max(a[3]-1-nghost,1)*max(a[4]-1-nghost,1)
            else: # nodes
                c = max(a[2]-nghost,1)*max(a[3]-nghost,1)*max(a[4]-nghost,1)
        elif len(a) == 4:
            if mode == 'cells':
                if isinstance(a[1], list):
                    if a[3] == 'NGON' or a[3] == 'NGON*':
                        c = a[2][3].size
                else: c = a[2][0].shape[0]
            else: # nodes
                if a[3][-1] != '*': 
                    if isinstance(a[1], list): c = a[1][0].size
                    else: c = a[1].shape[1]
                else: c = 0 # don't know
        nbPts.append(c)

    # Liste des arrays deja distribues
    if prescribed == []: # nothing set
        setArrays = [-1]*len(arrays)
    else: setArrays = prescribed
    
    # Liste des alpha, beta, gamma pour chaque processeur
    if perfo == []:
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
    if weight == []: weight = [1]*Nb
        
    # Matrice du volume des coms
    if com == []: volCom = numpy.zeros((Nb, Nb), numpy.int32)
    else: volCom = com
    
    # Distribution
    out = distributor2.distribute(nbPts, setArrays, perfProcs, weight,
                                  volCom, NProc, algorithm)
    return out
