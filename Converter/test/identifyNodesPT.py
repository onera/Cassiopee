# - identifyNodes (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P

a = G.cart((0,0,0), (1,1,1), (10,10,10))
hook = C.createHook(a, function='nodes')

# Indices des noeuds de a correspondant aux noeuds de f
f = P.exteriorFaces(a)
nodes = C.identifyNodes(hook, f); print(nodes)
#>> [   1    2    3    4    5    6    7  ...]
