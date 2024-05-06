# - add (array) -
import Generator as G
import CPlot
import time

a = G.cartTetra( (0,0,0), (1,1,1), (10,10,10) )
A = [a]
CPlot.display(A); time.sleep(1)

for i in range(10):
    b = G.cartTetra( (i*10,0,0), (1,1,1), (10,10,10) )
    CPlot.add(A, 0, b); CPlot.render(); time.sleep(0.5)

for i in range(10):
    b = G.cart( (i*10,10,0), (1,1,1), (10,10,10) )
    CPlot.add(A, 0, b); CPlot.render(); time.sleep(0.5)
