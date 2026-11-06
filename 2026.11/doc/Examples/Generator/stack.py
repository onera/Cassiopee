# - stack (array) -
import Generator as G
import Converter as C
import Transform as T
import Geom as D

# Concatenate 2 structured grids
a = G.cylinder((0,0,0), 1, 1.3, 360, 0, 1., (50,10,1))
b = T.rotate(a, (0,0,0), (1,0,0), 5.)
b = T.translate(b, (0,0,0.5))
c = G.stack(a, b)

# Concatenate a list of structured grids
a = []
for i in range(10):
    a.append(D.circle((0,0,i), 1.))
c = G.stack(a)

C.convertArrays2File(c, 'out.plt')
