# - changeBlanking (array) -
import Generator as G
import Converter as C
import CPlot
import time

a = G.cart((0,0,0), (1,1,1), (5,5,1))
a = C.initVars(a, 'cellN', 0)
CPlot.display(a, dim=2, mode=0)

CPlot.changeBlanking(); time.sleep(2)
CPlot.changeBlanking(); time.sleep(2)
