# - display (pyTree) -
# Affichage du shader sphere
import Geom.PyTree as D
import CPlot.PyTree as CPlot
import Converter.PyTree as C
import time

a = D.sphere((0,0,0),1.,N=5)
a = C.convertArray2Hexa(a)

for i in range(10):
    print(i*2./5.)
    CPlot._addRender2Zone(a, material='Sphere', color='White',
                          shaderParameters=[1.5,i*2./5.])

    CPlot.display(a, mode='render', bgColor=1)
    time.sleep(1.)
