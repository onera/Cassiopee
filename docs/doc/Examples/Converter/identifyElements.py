# - identifyElements (array) -
import Converter as C
import Generator as G
import Post as P

a = G.cartNGon( (0,0,0), (1,1,1), (10,10,10) )
f = P.exteriorElts(a)

# Enregistre les centres des elements dans le hook
hook = C.createHook(a, function='elementCenters')
# Indices des elements de a correspondant aux centres des elts de f
elts = C.identifyElements(hook, f); print(elts)
#>> [  1   2   3   4   5 ... 726 727 728 729]
