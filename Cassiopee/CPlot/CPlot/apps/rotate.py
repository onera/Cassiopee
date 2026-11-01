# - rotate -
# Rotation autour d'un objet
import Converter as C
import Transform as T
import CPlot

def rotate(a):
    for i in range(90*2):
        a = T.rotate(a, (0.5, 0, 0), (0,0,1), 0.5)
        CPlot.display(a)
    return a

a = C.convertFile2Arrays('dauphin2.plt', 'bin_tp')

CPlot.display(a, win=(700,500), displayBB=0, mode=0, meshStyle=0)
a = rotate(a)

CPlot.display(a, mode=1, solidStyle=0)
a = rotate(a)

CPlot.display(a, mode=2, colormap=1)
a = rotate(a)

CPlot.display(a, mode=2, colormap=0)
a = rotate(a)

CPlot.display(a, mode=2, colormap=1)
a = rotate(a)
