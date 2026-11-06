# - setState (pyTree) -
import Generator.PyTree as G
import CPlot.PyTree as CPlot
import time

a = G.cart((0,0,0), (1,1,1), (5,5,5))
CPlot.display([a], mode='solid')
time.sleep(0.2)
CPlot.setState(posCam=(8,8,8), posEye=(5,5,5))
