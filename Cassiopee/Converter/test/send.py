# - send (array) -
import Converter as C
import Generator as G
import Transform as T

a = G.cart((0,0,0), (1,1,1), (100,100,300))
C.send(a, 'localhost')

for i in range(30):
    a = T.rotate(a, (0,0,0), (0,0,1), 10.)
    C.send(a, 'localhost')
