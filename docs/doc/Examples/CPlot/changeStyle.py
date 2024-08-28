# - changeStyle (array) -
import Generator as G
import CPlot
import time

a = G.cart((0,0,0), (1,1,1), (5,5,1))
CPlot.display(a, dim=2, mode=1)

CPlot.changeStyle(); time.sleep(2)
CPlot.changeStyle(); time.sleep(2)
