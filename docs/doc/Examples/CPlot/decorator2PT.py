# - Decorator (pyTree) -
import CPlot.PyTree as CPlot
import CPlot.Decorator as Decorator
import Generator.PyTree as G
import Converter.PyTree as C

Decorator.setBatch(True)

a = G.cart((0,0,0), (1,1,1), (11,10,1))
C._initVars(a, '{F} = {CoordinateX}')

CPlot.display(a, mode='scalar',
              scalarField='F', isoScales=['F',12,0.,10.],
              export=CPlot.decorator,
              offscreen=2,
              isoEdges=1., colormap=16, bgColor=1)
CPlot.finalizeExport()

fig, ax = Decorator.createSubPlot()
Decorator.createText(ax, 0.4, 0.95, "Fast LES", size=40, box=True)
cbar = Decorator.createColorBar(fig, ax, title=r'$\mu_t / \mu$', location="right", color="black", fontSize=10, pad=-2.)

Decorator.savefig('out.png')
import os; os._exit(0)
