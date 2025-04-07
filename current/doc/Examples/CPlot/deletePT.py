# - delete (pyTree) -
import Generator.PyTree as G
import CPlot.PyTree as CPlot
import Converter.PyTree as C
import time

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
b = G.cartTetra( (11,0,0), (1,1,1), (10,10,10) )
c = G.cart( (0,11,0), (1,1,1), (10,10,10) )
t = C.newPyTree(['Base', a, b, c])

CPlot.display(t); time.sleep(1)
CPlot.delete(['Base\cartTetra']); CPlot.render(); time.sleep(1)
