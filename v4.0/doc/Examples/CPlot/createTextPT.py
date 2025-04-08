# - Decorator.createArraow (pyTree) -
import CPlot.PyTree as CPlot
import CPlot.Decorator as Decorator
import Generator.PyTree as G
import Converter.PyTree as C

Decorator.setBatch(True)

a = G.cart((0,0,0), (1,1,1), (10,10,10))
CPlot.display(a, mode='mesh',
              export=CPlot.decorator,
              offscreen=2,
              bgColor=1,
              posCam=(16.493430100289256, 23.422964886666588, -13.69258647223435),
              posEye=(4.5, 4.5, 4.5),
              dirCam=(0.7444990345982689, 0.1529242342484992, 0.6498733461696639),)
CPlot.finalizeExport()

fig, ax = Decorator.createSubPlot()

Decorator.createText(ax, 0.5, 0.8, text="Cube", size=20, color='black')
Decorator.createText(ax, 0.3, 0.4, text="Fluctuations", size=20, color='blue', rotation=90.)

Decorator.savefig('out.png')
import os; os._exit(0)
