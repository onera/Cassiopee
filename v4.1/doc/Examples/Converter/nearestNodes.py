# - nearestNodes (array) -
import Converter as C
import Generator as G
import Transform as T
import Post as P

a = G.cart((0,0,0), (1,1,1), (10,10,10))
# Enregistre les noeuds de a dans le hook
hook = C.createHook(a, function='nodes')

# Indices des noeuds de a les plus proches des noeuds de f
# et distance correspondante
b = T.translate(a,(0.15,0.,0.))
f = P.exteriorFaces(b)
nodes,dist = C.nearestNodes(hook, f); print(nodes, dist)
#>> [   1    2    3   ...] [0.15  0.15  0.15  ...]
