# - nearestFaces (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
b = G.cartNGon((9.1,0,0), (1,1,1), (10,10,10))

# Enregistre les centres des faces de a dans le hook
hook = C.createHook(a, function='faceCenters')
# Indices des faces de a les plus proches des faces de b
# et distance correspondante
faces,dist = C.nearestFaces(hook, b)
print(faces,dist)
