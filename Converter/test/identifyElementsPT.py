# - identifyElements (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P

a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
f = P.exteriorElts(a)

# Enregistre les centres des faces dans le hook
hook = C.createHook(a, function='elementCenters')
# Indices des faces de a correspondant aux centres des elts de f
elts = C.identifyElements(hook, f); print(elts)
#>> [  1   2   3  ... 726 727 728 729]
