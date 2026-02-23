# - silhouette ( array) -
import Generator as G
import Converter as C
import Post as P

a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,1,30))

vector=[1.,0.,0.]
res = P.silhouette([a], vector)

l = [a]+res
C.convertArrays2File(l, 'out.plt')
