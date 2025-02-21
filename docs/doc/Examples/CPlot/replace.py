# - replace (array) -
import Generator as G
import CPlot
import time

a = G.cart( (0,0,0), (1,1,1), (30,30,30) )
A = [a]
CPlot.display(A); time.sleep(1)

for i in range(10):
    b = G.cart( (i,0,0), (1,1,1), (30,30,30) )
    CPlot.replace(A, 0, b); CPlot.render(); time.sleep(0.1)
