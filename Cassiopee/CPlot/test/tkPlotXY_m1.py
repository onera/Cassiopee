# - tkPlotXY (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Mpi as Cmpi
import tkPlotXY
import KCore.test as test

LOCAL = test.getLocal()

# Cas test
if Cmpi.rank == 0:
    a = G.cart((0,0,0), (1,1,1), (101,1,1)); a[0] = 'cart0'
    C._initVars(a, '{F}={CoordinateX}*{CoordinateX}')
    C._initVars(a, '{centers:G}={centers:CoordinateX}')
elif Cmpi.rank == 1:
    a = G.cart((100,0,0), (1,1,1), (100,1,1)); a[0] = 'cart1'
    C._initVars(a, '{F}={CoordinateX}*{CoordinateX}')
    C._initVars(a, '{centers:G}={centers:CoordinateX}')
else:
    a = None

tkPlotXY.plot([a], varx='CoordinateX', vary='centers:G',
              xlabel='x',
              xformat='03.1f', yformat='03.1f',
              legends=['no1', 'no2'],
              xFontSize=11,
              legendFontSize=30,
              export=LOCAL+'/fig.png')
