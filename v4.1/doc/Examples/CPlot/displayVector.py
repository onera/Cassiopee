# - display (pyTree) -
# Affichage du shader vector (mode solid)
import Geom.PyTree as D
import CPlot.PyTree as CPlot
import Converter.PyTree as C

a = D.sphere((0,0,0),1.,N=50)
C._initVars(a, '{fx}={CoordinateY}')
C._initVars(a, '{fy}=0')
C._initVars(a, '{fz}=1')

#a = G.cart((0,0,0), (1,1,1), (10,10,1))
#a = C.convertArray2Tetra(a)
#C._initVars(a, '{fx}=0')
#C._initVars(a, '{fy}=0')
#C._initVars(a, '{fz}=1')

CPlot.display(a, mode='Vector', vectorField1=0, vectorField2=1, vectorField3=2, vectorStyle=1, vectorScale=50., vectorDensity=0)
#CPlot.display(a, mode='Vector', vectorField1=0, vectorField2=1, vectorField3=2, vectorStyle=1, vectorScale=80.)
#CPlot.display(a, mode='Vector', vectorField1='fx', vectorField2='fy', vectorField3='fz', vectorStyle=1, vectorScale=80., vectorNormalize=1)
#CPlot.display(a, mode='Vector', vectorField1=0, vectorField2=1, vectorField3=2, vectorStyle=0)
#CPlot.display(a, displayBB=0, mode='Vector', vectorField1=0, vectorField2=1, vectorField3=2, vectorStyle=1, vectorScale=0.25)
#CPlot.display(a, displayBB=0, mode='Vector', vectorField1=0, vectorField2=1, vectorField3=2, vectorStyle=2, vectorScale=0.25)
#CPlot.display(a, displayBB=0, mode='Vector', vectorField1='fx', vectorField2='fy', vectorField3='fz', vectorStyle=1, vectorScale=80., vectorNormalize=1)
#CPlot.display(a, displayBB=0, mode='Vector', vectorField1=0, vectorField2=1, vectorField3=2, vectorStyle=0)
