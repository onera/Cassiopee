"""Bridge to Cassiopee Solver.
"""
__version__ = '2.8'
__author__ = "Stephanie Peron, Christophe Benoit"
try:
    import kbridge
except:
    raise ImportError("KBridge module is unvailable.")

def evalKDesFunction(Func1, time):
    return kbridge.evalKDesFunction(Func1, time)

def array2KDesMesh(array, desMesh):
    return kbridge.array2KDesMesh(array, desMesh.this)
