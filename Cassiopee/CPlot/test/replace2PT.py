# - replace (pyTree) -
import Generator.PyTree as G
import CPlot.PyTree as CPlot
import Converter.PyTree as C
import time
dt = 0.2

a = G.cart( (0,0,0), (1,1,1), (30,30,30) )
b = G.cartHexa( (0,0,40), (1,1,1), (30,30,30) )
t = C.newPyTree(['Base', a, b])
CPlot.display(t); time.sleep(dt)

# Replace struct
for i in range(10):
    a = G.cart( (i,0,0), (1,1,1), (30,30,30) )
    CPlot.replace(t, 1, 0, a); CPlot.render(); time.sleep(dt)

# Replace non struct
for i in range(10):
    a = G.cartHexa( (i,0,40), (1,1,1), (30,30,30) )
    CPlot.replace(t, 1, 1, a); CPlot.render(); time.sleep(dt)
