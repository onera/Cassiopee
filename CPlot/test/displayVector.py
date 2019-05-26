# - display (pyTree) -
# Affichage du shader vector (mode solid)
import Geom.PyTree as D
import CPlot.PyTree as CPlot
import Converter.PyTree as C

a = D.sphere((0,0,0),1.,N=50)
#a = C.convertArray2Hexa(a)
C._initVars(a, '{fx}=-{CoordinateY}')
C._initVars(a, '{fy}={CoordinateX}')
C._initVars(a, '{fz}={CoordinateZ}')

C.convertPyTree2File(a, 'out.cgns')

CPlot.display(a, mode='Vector', vectorField1=0, vectorField2=1, vectorField3=2, vectorStyle=1, vectorScale=80.)
#CPlot.display(a, mode='Vector', vectorField1=0, vectorField2=1, vectorField3=2, vectorStyle=1, vectorScale=80.)
#CPlot.display(a, mode='Vector', vectorField1='fx', vectorField2='fy', vectorField3='fz', vectorStyle=1, vectorScale=80., vectorNormalize=1)
#CPlot.display(a, mode='Vector', vectorField1=0, vectorField2=1, vectorField3=2, vectorStyle=0)
#CPlot.display(a, displayBB=0, mode='Vector', vectorField1=0, vectorField2=1, vectorField3=2, vectorStyle=1, vectorScale=0.25)
#CPlot.display(a, displayBB=0, mode='Vector', vectorField1=0, vectorField2=1, vectorField3=2, vectorStyle=2, vectorScale=0.25)
#CPlot.display(a, displayBB=0, mode='Vector', vectorField1='fx', vectorField2='fy', vectorField3='fz', vectorStyle=1, vectorScale=80., vectorNormalize=1)
#CPlot.display(a, displayBB=0, mode='Vector', vectorField1=0, vectorField2=1, vectorField3=2, vectorStyle=0)


