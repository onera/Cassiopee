# - tkPlotXY (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import tkPlotXY

# Cas test
a = G.cart((0,0,0), (1,1,1), (100,1,1))
C._initVars(a, '{F}={CoordinateX}*{CoordinateX}')
C._initVars(a, '{centers:G}={centers:CoordinateX}')
b = G.cart((100,0,0), (1,1,1), (100,1,1))
C._initVars(b, '{F}={CoordinateX}*{CoordinateX}')
C._initVars(b, '{centers:G}={centers:CoordinateX}')

tkPlotXY.plot([a,b], varx='CoordinateX', vary='centers:G',
              xlabel='x',
              xformat='03.1f', yformat='03.1f',
              legends=['courbe1', 'courbe2'],
              export='fig.png')