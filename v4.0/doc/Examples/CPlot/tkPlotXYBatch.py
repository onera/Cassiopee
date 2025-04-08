# - tkPlotXY (batch) -
import Converter.PyTree as C
import Generator.PyTree as G
import tkPlotXY
tkPlotXY.setBatch()

a = G.cart((0,0,0), (1,1,1), (100,1,1))
C._initVars(a, '{F}={CoordinateX}*{CoordinateX}')
C._initVars(a, '{centers:G}={centers:CoordinateX}*{centers:CoordinateX}')

desktop = tkPlotXY.Desktop()
desktop.setData(a)

graph0 = desktop.createGraph('graph', '1:1')
curve0 = tkPlotXY.Curve(zone=['Base/cart'], varx='CoordinateX',
                        vary='F@FlowSolution', line_color='#7f00ff',
                        marker_face_color='#7f00ff',
                        marker_edge_color='#7f00ff')
graph0.addCurve('1:1', curve0)

# Modify axis range
axis = graph0.getAxis('1:1')
axis.x.setValue('axis_autoscale', False)
axis.x.setValue('axis_min', 30.)
axis.x.setValue('axis_max', 50.)
axis.y.setValue('axis_autoscale', False)
axis.y.setValue('axis_min', 1000.)
axis.y.setValue('axis_max', 2000.)
graph0.updateGraph('1:1')

graph0.save('fig.png')
graph0.close()
