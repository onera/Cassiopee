# - identifyFaces (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
b = G.cartNGon((9,0,0), (1,1,1), (10,10,10))

# Enregistre les centres des faces de a dans le hook
hook = C.createHook(a, function='faceCenters')
# Indices des faces de a correspondant aux faces de b
faces = C.identifyFaces(hook, b); print(faces)
#>> [10 -1 -1 ..., -1 -1 -1]
