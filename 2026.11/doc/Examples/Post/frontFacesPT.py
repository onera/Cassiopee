# - frontFaces (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P

a = G.cart( (0,0,0), (1,1,1), (11,11,11) )

def F(x, y, z):
    return (x + 2*y + z > 20.).astype(float)

a = C.initVars(a, 'tag', F, ['CoordinateX', 'CoordinateY', 'CoordinateZ'], isVectorized=True)
f = P.frontFaces(a, 'tag'); f[0] = 'front'
t = C.newPyTree(['Base']); t[2][1][2] += [a, f]
C.convertPyTree2File(t, 'out.cgns')
