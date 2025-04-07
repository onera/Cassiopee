# - computePressureGradients (pyTree) -
import Converter.Internal as Internal
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D
import KCore.test as test
import Post.IBM as P_IBM
import Geom.IBM as D_IBM
import Connector.IBM as X_IBM
import copy
import numpy

tb = D.circle((0.,0.,0.), 1., N=100)
tb = C.newPyTree(['Base', tb])
D_IBM._setSnear(tb, 0.05)
D_IBM._setDfar(tb, 10.)
D_IBM._setIBCType(tb, 'Musker')
C._addState(tb, adim='adim1', MInf=0.2, alphaZ=0., alphaY=0., ReInf=5.e6, EquationDimension=2, GoverningEquations='NSTurbulent')

a, ac = X_IBM.prepareIBMData(tb, t_out=None, tc_out=None, vmin=21, frontType=1, check=False)

ac = P_IBM.extractYplusAtImagePoints(ac)

C.convertPyTree2File(ac, 'out.cgns')
