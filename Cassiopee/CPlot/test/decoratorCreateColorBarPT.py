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
              colormap=0, # Default Blue2Red colorbar
              export=CPlot.decorator, # Plot.decorator = '.decorator.png'
              offscreen=offscreen,
              viewAngle=viewAngle,
              posCam=posCam, posEye=posEye, dirCam=dirCam,
              exportResolution=exportResolution)

CPlot.finalizeExport(offscreen)

# create image with title and graduated axes
fig, ax = Decorator.createSubPlot(box=True, title='F(x,y) = x*y', xlim=xlim, ylim=ylim)

# add multicolor colorbar using CPlot information
cbar = Decorator.createColorBar(fig, ax, title='F')
Decorator.savefig('out4.png', pad=0.1); cbar.remove()

# add multicolor colorbar using cmap name
cmap = 'Blue2Red'
cbar = Decorator.createColorBar(fig, ax, title='F', levels=[25,0.5,1.9], cmap=cmap)
Decorator.savefig('out5.png', pad=0.1); cbar.remove()

# add multicolor colorbar using cmap list
cmap = [(0.00,[0,0,1]), (0.25,[0,1,1]), (0.50,[0,1,0]), (0.75,[1,1,0]), (1.00,[1,0,0])]
cbar = Decorator.createColorBar(fig, ax, title='F', levels=[25,0.5,1.9], cmap=cmap)
Decorator.savefig('out6.png', pad=0.1)

if offscreen == 2: import os; os._exit(0)