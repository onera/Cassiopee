# - setDim (array) -
import Generator as G
import CPlot
import time

a = G.cart( (0,0,0), (1,1,1), (5,5,1) ); t = [a]
CPlot.display(t); time.sleep(2)
CPlot.setDim(2); time.sleep(2)
