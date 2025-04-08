# - send (array) -
import Converter as C
import Generator as G
import Transform as T
import Converter.Mpi as Cmpi

if Cmpi.rank == 0:
    a = G.cart((0,0,0), (1,1,1), (50,50,30))
elif Cmpi.rank == 1:
    a = G.cart((0,100,0), (1,1,1), (50,50,30))
else:
    a = G.cart((100,0,0), (1,1,1), (50,50,30))
C.send(a, 'localhost', rank=Cmpi.rank)
Cmpi.barrier()

for i in range(10):
    #import time; time.sleep(2)
    a = T.rotate(a, (0,0,0), (0,0,1), 10.)
    C.send(a, 'localhost', rank=Cmpi.rank)
    Cmpi.barrier()
