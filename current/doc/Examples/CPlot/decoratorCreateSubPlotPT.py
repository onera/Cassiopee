import Converter.PyTree as C
import Generator.PyTree as G
import CPlot.PyTree as CPlot
import CPlot.Decorator as Decorator

import numpy

offscreen = 2
Decorator.setBatch(True)

a = G.cart((0,0,0), (0.01,0.01,1), (201,101,1))
C._initVars(a, '{F} = {CoordinateX}*{CoordinateY}')

ppw = 1000 # pixels per height
zplane = 0.0 # constant-z value
xlim = [1.0, 2.0] # xmin, xmax
ylim = [0.5, 1.0] # ymin, ymax
posCam, posEye, dirCam, viewAngle, exportResolution = Decorator.getInfo2DMode(xlim, ylim, zplane, ppw)

CPlot.display(a, mode='Scalar', scalarField='F',
    isoScales=['F',25,0.5,1.9], 
    colormap=24, # Jet colorbar
    export=CPlot.decorator, # Plot.decorator = '.decorator.png'
    offscreen=offscreen, 
    viewAngle=viewAngle,
    posCam=posCam, posEye=posEye, dirCam=dirCam,
    exportResolution=exportResolution)

CPlot.finalizeExport(offscreen)

# create raw image
fig, ax = Decorator.createSubPlot()
Decorator.savefig('out0.png', pad=0.0); Decorator.closeAll()

# create image with title and borders
fig, ax = Decorator.createSubPlot(box=True, title='F(x,y) = x*y')
Decorator.savefig('out1.png', pad=0.1); Decorator.closeAll()

# create image with title and graduated axes
fig, ax = Decorator.createSubPlot(box=True, title='F(x,y) = x*y', xlim=xlim, ylim=ylim)
Decorator.savefig('out2.png', pad=0.1)

# add modifications with Matplotlib commands
cx = numpy.linspace(xlim[0], xlim[1], 100)
ax.plot(cx, 1.20/cx, color='k', linestyle='dashed')
bbox = dict(boxstyle="round, pad=0.4, rounding_size=0.5", ec='k', fc="white", linewidth=1.5)
ax.text(1.5, 0.75, 'iso 1.2', size=12, ha='center', va='center', color='k', bbox=bbox, rotation=-27.)
Decorator.savefig('out3.png', pad=0.1)

if offscreen == 2: import os; os._exit(0)