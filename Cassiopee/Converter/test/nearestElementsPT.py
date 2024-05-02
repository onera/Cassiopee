# - nearestElements (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import Post.PyTree as P

a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
b = T.translate(a,(0.15,0.,0.))
f = P.exteriorElts(b)

# Enregistre les centres des faces dans le hook
hook = C.createHook(a, function='elementCenters')
# Indices des faces de a les plus proches des centres des elts de f
# et distance correspondante
elts,dist = C.nearestElements(hook, f)
print(elts,dist)
