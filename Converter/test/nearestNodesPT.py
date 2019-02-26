# - nearestNodes (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import Post.PyTree as P

a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = T.translate(a,(0.15,0.,0.))
f = P.exteriorFaces(b)

hook = C.createHook(a, function='nodes')
# Indices des noeuds de a les plus proches des noeuds de f
# et distance correspondante
nodes,dist = C.nearestNodes(hook, f)
print(nodes, dist)
