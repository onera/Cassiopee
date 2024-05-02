# - Decorator (pyTree) -
import CPlot.PyTree as CPlot
import CPlot.Decorator as Decorator
import Generator.PyTree as G
import Converter.PyTree as C

Decorator.setBatch(True)

a = G.cart((0,0,0), (1,1,1), (10,10,1))
C._initVars(a, '{F} = {CoordinateX}')

CPlot.display(a, mode='Scalar',
              scalarField='F', isoScales=['F',12,0.,10.],
              export=CPlot.decorator, exportResolution='3840x2160', offscreen=1,
              isoEdges=1., colormap=0, bgColor=1)

fig, ax = Decorator.createSubPlot()
ax.set_title('Computation of the year', size=40)
Decorator.createText(ax, 0.02, 0.9, "Fast LES", size=40, box=True)
cbar = Decorator.createColorBar(fig, ax, title=r'$\mu_t / \mu$')
cbar.ax.tick_params(labelcolor='tab:red')

Decorator.savefig('out.png')
