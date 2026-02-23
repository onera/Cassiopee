# - changeVariable (array) -
import Generator as G
import Converter as C
import CPlot
import time

def F(x,y): return x*x + y*y

a = G.cart((0,0,0), (1,1,1), (5,5,1))
a = C.addVars(a, 'Density')
a = C.initVars(a, 'F', F, ['x','y'])
CPlot.display(a, dim=2, mode=3)

CPlot.changeVariable(); time.sleep(2)
CPlot.changeVariable(); time.sleep(2)
