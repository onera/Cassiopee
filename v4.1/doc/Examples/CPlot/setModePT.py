# - setMode (pyTree) -
import Generator.PyTree as G
import CPlot.PyTree as CPlot
import time

a = G.cart((0,0,0), (1,1,1), (5,5,1))
CPlot.display(a, dim=2)

CPlot.setMode(1); time.sleep(2)
CPlot.setMode(0); time.sleep(2)
CPlot.setMode(1); time.sleep(2)
