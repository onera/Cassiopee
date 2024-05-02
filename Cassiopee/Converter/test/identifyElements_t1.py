# - identifyElements (array) -
import Converter as C
import Post as P
import KCore.test as test

def identify(a):
    f = P.exteriorElts(a)
    if isinstance(a[0], list): a=a[0]
    hook = C.createHook(a, function='elementCenters')
    elts = C.identifyElements(hook, f)
    C.freeHook(hook)
    return elts

test.stdTestA(identify)
