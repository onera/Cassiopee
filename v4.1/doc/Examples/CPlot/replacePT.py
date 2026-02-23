# - replace (pyTree) -
import Generator.PyTree as G
import CPlot.PyTree as CPlot
import Converter.PyTree as C
import time

a = G.cartTetra( (0,0,0), (1,1,1), (10,10,10) )
b = G.cartTetra( (11,0,0), (1,1,1), (10,10,10) )
c = G.cart( (0,11,0), (1,1,1), (10,10,10) )

t = C.newPyTree(['Base', a, b, c])
CPlot.display(t); time.sleep(2.)

d = G.cart( (11,11,0), (1,1,1), (10,10,10) )
CPlot.replace(t, 1, 0, d); CPlot.render(); time.sleep(1.)
