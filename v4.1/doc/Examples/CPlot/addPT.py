# - add (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import time

a = G.cart((0,0,0), (1,1,1), (30,30,30))

t = C.newPyTree(['Base', a])
CPlot.display(t); time.sleep(1)

b = G.cart((30,0,0), (1,1,1), (30,30,30))
CPlot.add(t, 1, -1, b); CPlot.render(); time.sleep(1)
