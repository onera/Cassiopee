# - send (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T

a = G.cart((0,0,0), (1,1,1), (100,100,30))
C.send(a, 'localhost')

for i in range(30):
    a = T.rotate(a, (0,0,0), (0,0,1), 10.)
    C.send(a, 'localhost')
