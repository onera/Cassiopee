# - identifyNodes (array) -
import Converter as C
import Generator as G
import Post as P

a = G.cart((0,0,0), (1,1,1), (10,10,10))
# Enregistre les noeuds de a dans le hook
hook = C.createHook(a, function='nodes')

# Indices des noeuds de a correspondant aux noeuds de f
f = P.exteriorFaces(a)
nodes = C.identifyNodes(hook, f); print(nodes)
#>> [   1    2    3    4    5    6    7  ...]
