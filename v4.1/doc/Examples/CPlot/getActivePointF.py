# - getActivePointF (array) -
import Generator as G
import Converter as C
import CPlot
import time

a = G.cart((0,0,0), (1,1,1), (5,5,5))
a = C.initVars(a, '{F}={x}')
CPlot.display([a])

l = []
while l == []:
    l = CPlot.getActivePointF(); time.sleep(0.1)
print('ActivePoint: ', l)
