# - nearestElements (array) -
import Converter as C
import Transform as T
import Post as P
import KCore.test as test

def nearest(a):
    f = P.exteriorElts(a)
    f = T.translate(f,(0.15,0.,0.))
    if isinstance(a[0], list): a=a[0]
    hook = C.createHook(a, function='elementCenters')
    elts = C.nearestElements(hook, f)
    C.freeHook(hook)
    return elts

test.stdTestA(nearest)
