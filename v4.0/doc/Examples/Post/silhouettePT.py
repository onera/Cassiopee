# - silhouette (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Post.PyTree as P

a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,1,30))
t = C.newPyTree(['Base']); t[2][1][2].append(a)

vector=[1.,0.,0.]
res = P.silhouette(a, vector)
t[2][1][2] += res
C.convertPyTree2File(t, 'out.cgns')
