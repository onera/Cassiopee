# - display (array) -
# test scalar mode -
import Converter as C
import Generator as G
import CPlot
import time

a = G.cartTetra( (0,0,0), (1,1,1), (1,10,10) )
a = C.initVars(a, '{f}={y}')

CPlot.display([a], mode=3, scalarStyle=0); time.sleep(1.)
CPlot.display([a], mode=3, scalarStyle=1); time.sleep(1.)
CPlot.display([a], mode=3, scalarStyle=0, colormap=24); time.sleep(1.)
