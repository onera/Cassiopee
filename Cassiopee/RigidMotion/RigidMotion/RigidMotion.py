"""Module performing rigid motion.
"""
__version__ = '4.1'
__author__ = "Stephanie Peron, Christophe Benoit, Pascal Raud"

from . import rigidMotion
import numpy
try: import Converter as C
except ImportError:
    raise ImportError("RigidMotion: requires Converter module.")

#==============================================================================
# IN: time: instant d'evaluation
# IN: F: une fonction du temps de la forme (d, c, R), deplacement, centre,
# de rotation, matrice de rotation.
# si time est different de None, sinon des scalaires.
#==============================================================================
def evalPosition(array, time, F):
    """Move the mesh with defined motion function to time t."""
    b = C.copy(array)
    _evalPosition(b, time, F)
    return b

def _evalPosition(array, time, F):
    """Move the mesh with defined motion function to time t.
    Return an array with moved mesh coordinates.
    Usage: evalPosition(array, time, F)"""
    if time is not None:
        if time < 0:
            f = F(-time)
            if len(f) != 3:
                raise ValueError("evalPosition: f must be a 3D function.")
            c = f[0]; d = f[1]; r0 = f[2]
            r = numpy.transpose(r0)
        else:
            f = F(time)
            if len(f) != 3:
                raise ValueError("evalPosition: f must be a 3D function.")
            d = f[0]; c = f[1]; r = f[2]
    else: d = F[0]; c = F[1]; r = F[2]

    if len(d) != 3: raise ValueError("evalPosition: d must be a 3D vector.")
    if len(c) != 3: raise ValueError("evalPosition: c must be a 3D vector.")
    if len(r) != 3: raise ValueError("evalPosition: rotation matrix must be 3x3.")
    else:
        if len(r[0]) != 3: raise ValueError("evalPosition: rotation matrix must be 3x3.")
    if isinstance(array[0], list):
        for i in array:
            rigidMotion.move(i, d[0], d[1], d[2], \
                             c[0], c[1], c[2], \
                             r[0][0], r[0][1], r[0][2], \
                             r[1][0], r[1][1], r[1][2], \
                             r[2][0], r[2][1], r[2][2])
        return None
    else:
        return rigidMotion.move(array, d[0], d[1], d[2], \
                                c[0], c[1], c[2], \
                                r[0][0], r[0][1], r[0][2], \
                                r[1][0], r[1][1], r[1][2], \
                                r[2][0], r[2][1], r[2][2])

# Applique la formule XP=d+r*(XN-c) sur des numpys de coordonnees
# in place
def _moveN(coordsN, d, c, r):
    return rigidMotion.moveN(coordsN,
                             (d[0], d[1], d[2]), \
                             (c[0], c[1], c[2]), \
                             r[0,0], r[0,1], r[0,2], \
                             r[1,0], r[1,1], r[1,2], \
                             r[2,0], r[2,1], r[2,2])

def moveN(coordsN, d, c, r):
    coordsNP = []
    for i in range(3):
        XP = numpy.copy(coordsN[i])
        coordsNP.append(XP)

    _moveN(coordsNP, d, c, r)
    return coordsNP
