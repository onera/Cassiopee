# - add (array) -
import Generator as G
import CPlot
import time

a = G.cart((0,0,0), (1,1,1), (30,30,30))
A = [a]
CPlot.display(A); time.sleep(1)

b = G.cart((30,0,0), (1,1,1), (30,30,30))
CPlot.add(A, 0, b); CPlot.render(); time.sleep(0.5)
