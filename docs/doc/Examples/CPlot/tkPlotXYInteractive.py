# - tkPlotXY (interactive) -
import Converter.PyTree as C
import Generator.PyTree as G
import tkPlotXY

a = G.cart((0,0,0), (1.e-3,1,1), (100,1,1))
C._initVars(a, '{F}={CoordinateX}*{CoordinateX}')
C._initVars(a, '{centers:G}={centers:CoordinateX}*{centers:CoordinateX}')

desktop, win = tkPlotXY.createTkDesktop()
desktop.setData(a)

graph0 = desktop.createGraph('graph', '1:1')
# La localisation des champs dans Curve doit etre homogene
curve0 = tkPlotXY.Curve(zone=['Base/cart'], varx='CoordinateX',
                        vary='F@FlowSolution', line_color='#7f00ff',
                        marker_face_color='#7f00ff',
                        marker_edge_color='#7f00ff',
                        legend_label='expe0')
graph0.addCurve('1:1', curve0)

win.mainloop()
