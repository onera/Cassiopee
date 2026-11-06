# - display HO -
import Geom.PyTree as D
import Converter.PyTree as C
import Transform.PyTree as T
import Converter.Internal as Internal
import CPlot.PyTree as CPlot

from math import sqrt

#qsph = D.quadrangle((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.))
qsph = D.sphere6((0.,0.,0.), 1.,N=3, ntype='QUAD')
#qsph = C.convertLO2HO(qsph, 1)
qsph02 = C.convertLO2HO(qsph, 0, 2)
for i in range(C.getNPts(qsph02)):
    x = C.getValue(qsph02,'CoordinateX',i)
    y = C.getValue(qsph02,'CoordinateY',i)
    z = C.getValue(qsph02,'CoordinateZ',i)
    nrm = sqrt(x*x+y*y+z*z)
    x /= nrm; y /= nrm; z /= nrm
    C.setValue(qsph02,'CoordinateX',i,x)
    C.setValue(qsph02,'CoordinateY',i,y)
    C.setValue(qsph02,'CoordinateZ',i,z)

qsph12 = C.convertLO2HO(qsph, 1, 2)
for i in range(C.getNPts(qsph12)):
    x = C.getValue(qsph12,'CoordinateX',i)
    y = C.getValue(qsph12,'CoordinateY',i)
    z = C.getValue(qsph12,'CoordinateZ',i)
    nrm = sqrt(x*x+y*y+z*z)
    x /= nrm; y /= nrm; z /= nrm
    C.setValue(qsph12,'CoordinateX',i,x)
    C.setValue(qsph12,'CoordinateY',i,y)
    C.setValue(qsph12,'CoordinateZ',i,z)

qsph03 = C.convertLO2HO(qsph, 0, 3)
for i in range(C.getNPts(qsph03)):
    x = C.getValue(qsph03,'CoordinateX',i)
    y = C.getValue(qsph03,'CoordinateY',i)
    z = C.getValue(qsph03,'CoordinateZ',i)
    nrm = sqrt(x*x+y*y+z*z)
    x /= nrm; y /= nrm; z /= nrm
    C.setValue(qsph03,'CoordinateX',i,x)
    C.setValue(qsph03,'CoordinateY',i,y)
    C.setValue(qsph03,'CoordinateZ',i,z)

qsph13 = C.convertLO2HO(qsph, 1, 3)
for i in range(C.getNPts(qsph13)):
    x = C.getValue(qsph13,'CoordinateX',i)
    y = C.getValue(qsph13,'CoordinateY',i)
    z = C.getValue(qsph13,'CoordinateZ',i)
    nrm = sqrt(x*x+y*y+z*z)
    x /= nrm; y /= nrm; z /= nrm
    C.setValue(qsph13,'CoordinateX',i,x)
    C.setValue(qsph13,'CoordinateY',i,y)
    C.setValue(qsph13,'CoordinateZ',i,z)

#tsph = D.triangle((0.,0.,0.),(1.,0.,0.),(0.,1.,0.))#D.sphere6((0,0,0), 1., N=3, ntype='TRI')
tsph = D.sphere6((0,0,0), 1., N=3, ntype='TRI')
tsph02 = C.convertLO2HO(tsph, 0, 2)
for i in range(C.getNPts(tsph02)):
    x = C.getValue(tsph02,'CoordinateX',i)
    y = C.getValue(tsph02,'CoordinateY',i)
    z = C.getValue(tsph02,'CoordinateZ',i)
    nrm = sqrt(x*x+y*y+z*z)
    x /= nrm; y /= nrm; z /= nrm
    C.setValue(tsph02,'CoordinateX',i,x)
    C.setValue(tsph02,'CoordinateY',i,y)
    C.setValue(tsph02,'CoordinateZ',i,z)

tsph03 = C.convertLO2HO(tsph, 0, 3)
for i in range(C.getNPts(tsph03)):
    x = C.getValue(tsph03,'CoordinateX',i)
    y = C.getValue(tsph03,'CoordinateY',i)
    z = C.getValue(tsph03,'CoordinateZ',i)
    nrm = sqrt(x*x+y*y+z*z)
    x /= nrm; y /= nrm; z /= nrm
    C.setValue(tsph03,'CoordinateX',i,x)
    C.setValue(tsph03,'CoordinateY',i,y)
    C.setValue(tsph03,'CoordinateZ',i,z)

tsph13 = C.convertLO2HO(tsph, 1, 3)
for i in range(C.getNPts(tsph13)):
    x = C.getValue(tsph13,'CoordinateX',i)
    y = C.getValue(tsph13,'CoordinateY',i)
    z = C.getValue(tsph13,'CoordinateZ',i)
    nrm = sqrt(x*x+y*y+z*z)
    x /= nrm; y /= nrm; z /= nrm
    C.setValue(tsph13,'CoordinateX',i,x)
    C.setValue(tsph13,'CoordinateY',i,y)
    C.setValue(tsph13,'CoordinateZ',i,z)

qsph12 = T.translate(qsph12, (3.,0,0))
qsph03 = T.translate(qsph03, (6.,0,0))
qsph13 = T.translate(qsph13, (9.,0,0))

tsph02 = T.translate(tsph02, (0.,3,0))
tsph03 = T.translate(tsph03, (6.,3,0))
tsph13 = T.translate(tsph13, (9.,3,0))
t = C.newPyTree(['Quadrangle 0 2', 'Quadrangle 1 2', 'Quadrangle 0 3', 'Quadrangle 1 3', 'Triangle 0 2', 'Triangle 0 3', 'Triangle 1 3'])
t[2][1][2] += Internal.getZones(qsph02)
t[2][2][2] += Internal.getZones(qsph12)
t[2][3][2] += Internal.getZones(qsph03)
t[2][4][2] += Internal.getZones(qsph13)

t[2][5][2] += Internal.getZones(tsph02)
t[2][6][2] += Internal.getZones(tsph03)
t[2][7][2] += Internal.getZones(tsph13)

for shader in ["Chrome","Glass", "Wood", 'Marble', 'XRay', 'Metal', 'Granite', 'Brick', 'Cloud', 'Gooch']:
    s = CPlot.addRender2Zone(t, material=shader, color='White')
    CPlot.display(s, mode="Render")
    l = ''
    while l != 'n':
        l = CPlot.getKeyboard()
        CPlot.resetKeyboard()