# Display 1D data from CGNS Tree
import tkPlotXY
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (20,1,1))
C._initVars(a, '{F}={CoordinateX}*{CoordinateX}')
t = C.newPyTree(['Base',a])

desktop, win = tkPlotXY.createTkDesktop()
desktop.setData(t)

graph0 = desktop.createGraph('MyFirstGraph', '1:1')
curve0 = tkPlotXY.Curve(zone=['Base/cart'], varx='CoordinateX', vary='F@FlowSolution', line_color='#7f00ff', marker_face_color='#7f00ff', marker_edge_color='#7f00ff')
graph0.addCurve('1:1', curve0)

desktop.cmd_editCurves()
win.mainloop()
