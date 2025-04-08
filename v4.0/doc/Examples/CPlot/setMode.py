# - setMode (array) -
import Generator as G
import CPlot
import time

a = G.cart((0,0,0), (1,1,1), (5,5,1))
CPlot.display(a, mode=0, dim=2); time.sleep(2)
CPlot.setMode(1); time.sleep(2) # solid
CPlot.setMode('mesh'); time.sleep(2)
CPlot.setMode('solid'); time.sleep(2)
