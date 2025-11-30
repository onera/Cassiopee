import Converter.PyTree as C
import Generator.PyTree as G
import CPlot.PyTree as CPlot
import CPlot.Decorator as Decorator

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

if offscreen == 2: import os; os._exit(0)
