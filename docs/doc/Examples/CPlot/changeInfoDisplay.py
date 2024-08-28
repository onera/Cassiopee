# - changeInfoDisplay (array) -
import Generator as G
import CPlot
import time

a = G.cart( (0,0,0), (1,1,1), (5,5,1) ); t = [a]
CPlot.display(t, dim=2, mode=0)

CPlot.changeInfoDisplay(); time.sleep(1)
CPlot.changeInfoDisplay(); time.sleep(1)
CPlot.changeInfoDisplay(); time.sleep(1)
