import tkPlotXY
import numpy
import Converter.PyTree as C
import Converter.Internal as Internal

# Change some default parameters of tkPlotXY by using the default_values dictionnary
# AXIS
tkPlotXY.default_values['Axis']['axis_x_min'] = 0.
tkPlotXY.default_values['Axis']['axis_x_max'] = 1.
tkPlotXY.default_values['Axis']['axis_x_autoscale'] = False
# GRID
tkPlotXY.default_values['Grid']['Mx_grid_color'] = 'black'
tkPlotXY.default_values['Grid']['My_grid_color'] = 'black'
# LEGEND
tkPlotXY.default_values['Legend']['legend_border_color'] = 'black'
# CURVE
tkPlotXY.default_values['Curve']['line_width'] = 1.5

# Load cgns files
walls = {}
walls['wall 1'] = C.convertFile2PyTree('wall1.cgns')
walls['wall 2'] = C.convertFile2PyTree('wall2.cgns')
# Create a data zone for each cgns file
# Keys can be later used as legend label
data = {}
for w in walls:
    data[w] = {
    'x/c':Internal.getNodeFromName(walls[w], 'CoordinateX')[1],
    'Cp' :Internal.getNodeFromName(walls[w], 'Cp')[1],
    'Cf' :Internal.getNodeFromName(walls[w], 'Cf')[1],
    }

# Plot
desktop, win = tkPlotXY.createTkDesktop()
desktop.setData(data)

# Create the first graph
graph0 = desktop.createGraph('Graph_Cf', '1:1', figsize=(6,4))
graph0.updateSubPlotParams({'isActive':True, 'right':0.97, 'top':0.97, 'left':0.2, 'bottom':0.13, 'wspace':0.3})
for w in data:
    # Each new curve must be created first and then added to the graph with addCurve()
    curve = tkPlotXY.Curve(zone=[w], varx='x/c', vary='Cf', legend_label=w)
    graph0.addCurve('1:1', curve)


# Create the second graph
graph1 = desktop.createGraph('Graph_Cp', '1:1', figsize=(6,4))
graph1.updateSubPlotParams({'isActive':True, 'right':0.97, 'top':0.97, 'left':0.2, 'bottom':0.13, 'wspace':0.3})
for w in data:
    curve = tkPlotXY.Curve(zone=[w], varx='x/c', vary='Cp', legend_label=w)
    graph1.addCurve('1:1', curve)

# Open window for curve editing when starting tkPlotXY
desktop.cmd_editCurves()
win.mainloop()
