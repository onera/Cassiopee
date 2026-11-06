# - patch (array) -
import Transform as T
import Generator as G
import Converter as C
import numpy

c1 = G.cart((0,0,0), (0.01,0.01,1), (201,101,1))
c2 = G.cart((0,0,0), (0.01,0.01,1), (51,81,1))
c2 = T.rotate(c2, (0,0,0),(0,0,1),0.2)
c3 = G.cart((0.0,1.,0), (0.01,0.01,1), (101,1,1))
c3 = T.rotate(c3, (0,0,0),(0,0,1),0.3)
# patch a region at given position
a = T.patch(c1, c2, position=(1,1,1))
# patch some nodes
nodes = numpy.arange(20100, 20201, dtype=numpy.int32)
b = T.patch(c1, c3, nodes=nodes)
C.convertArrays2File([a,b], 'out.plt')
