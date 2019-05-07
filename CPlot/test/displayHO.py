# - display HO -
import Geom.PyTree as GP
import Converter.PyTree as C
import Transform.PyTree as T
import Converter.Internal as Ci
import CPlot.PyTree as CPlot

from math import sqrt
import sys

qsph = GP.sphere6((0.,0.,0.), 1.,N=3, ntype='QUAD')
qsph = C.convertLO2HO(qsph, 1)
for i in range(C.getNPts(qsph)):
    x = C.getValue(qsph,'CoordinateX',i)
    y = C.getValue(qsph,'CoordinateY',i)
    z = C.getValue(qsph,'CoordinateZ',i)
    nrm = sqrt(x*x+y*y+z*z);
    x /= nrm;
    y /= nrm;
    z /= nrm;
    C.setValue(qsph,'CoordinateX',i,x)
    C.setValue(qsph,'CoordinateY',i,y)
    C.setValue(qsph,'CoordinateZ',i,z)

tsph = GP.sphere6((0,0,0), 1., N=3, ntype='TRI')
tsph = C.convertLO2HO(tsph, 0)
for i in range(C.getNPts(tsph)):
    x = C.getValue(tsph,'CoordinateX',i)
    y = C.getValue(tsph,'CoordinateY',i)
    z = C.getValue(tsph,'CoordinateZ',i)
    nrm = sqrt(x*x+y*y+z*z);
    x /= nrm;
    y /= nrm;
    z /= nrm;
    C.setValue(tsph,'CoordinateX',i,x)
    C.setValue(tsph,'CoordinateY',i,y)
    C.setValue(tsph,'CoordinateZ',i,z)

tsph = T.translate(tsph, (3.,0,0))
t = C.newPyTree(['Base1', 'Base2'])
t[2][1][2] += Ci.getZones(qsph)
t[2][2][2] += Ci.getZones(tsph)

for shader in ["Chrome","Glass", "Wood", 'Marble', 'XRay', 'Metal', 'Granite', 'Brick', 'Cloud', 'Gooch']:
    s = CPlot.addRender2Zone(t, material=shader, color='White')    
    CPlot.display(s, mode="Render")
    l = ''
    while(l != 'n'):
        l = CPlot.getKeyboard()
        CPlot.resetKeyboard()
