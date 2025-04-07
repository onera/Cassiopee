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
    l = CPlot.getKeyboard(); CPlot.resetKeyboard()
    try:
        a = D.text2D(l)
        CPlot.display(a)
    except:
        for i in range(len(l)):
            v = ord(l[i]); print(v)
            if v == 1: print('up')
            elif v == 2: print('down')
            elif v == 3: print('left')
            elif v == 4: print('right')
            elif v == 5: print('release up')
            elif v == 6: print('release down')
            elif v == 7: print('release left')
            elif v == 8: print('release right')
        time.sleep(0.1)
        l = ''
