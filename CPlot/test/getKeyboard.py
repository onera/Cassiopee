# - getKeyboard (array) -
import Generator as G
import CPlot
import Geom as D
import time

a = G.cart((0,0,0), (1,1,1), (8,8,1))
CPlot.display(a, dim=2)
CPlot.setState(activateShortCuts=0)

for i in range(50):
    l = ''
    while l == '':
        l = CPlot.getKeyboard(); time.sleep(0.1)
    try:
        a = D.text2D(l)
        CPlot.display(a)
    except:
        v = ord(l[0])
        if v == 1: print('up')
        elif v == 2: print('down')
        elif v == 3: print('left')
        elif v == 4: print('right')
        time.sleep(0.1)
        l = ''
