# - setState (array) -
import Generator as G
import CPlot
import time

a = G.cart((0,0,0), (1,1,1), (5,5,5))
CPlot.display(a, mode='solid')
time.sleep(1.)
CPlot.setState(posCam=(8,8,8), posEye=(5,5,5))
