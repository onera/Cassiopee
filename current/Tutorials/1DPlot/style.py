# coding: utf-8
from __future__ import unicode_literals
from tkPlotXY import *
def loadVisu(obj):
    ###################################
    # Graph 0
    #########################
    graph_0=GraphTK(obj,'Graph_Cf','1:1',dpi=None,figsize=(6, 4))

    graph_0.updateSubPlotParams({'isActive':True,'left':0.2,'right':0.97,'top':0.97,'bottom':0.13,'wspace':0.3})
    #########################
    # SubGraph '1:1'
    ###############
    ######
    # Axis
    ######
    axis_0_0_0 = graph_0.getAxis('1:1',0)
    ######
    # Curves
    ######
    curve_0 = Curve(zone=['wall1'],varx='x/c',vary='Cf' , line_color='#e41a1c' , marker_style='pentagon' , marker_face_color='#e41a1c' ,marker_edge_color='#e41a1c' , marker_sampling_step=60, legend_label='wall1' )
    graph_0.addCurve('1:1',curve_0)
    ######
    # Grids
    ######
    grid_0_0_0 = graph_0.getGrid('1:1',axis=axis_0_0_0)
    ######
    # Axis settings
    ######
    axis_0_0_0.setValue('x','axis_label','bonjour')
    axis_0_0_0.setValue('y','axis_label','hello')
    ######
    # Legends
    ######
    legend_0_0 = graph_0.getLegend('1:1')
    ###################################
    # Graph 1
    #########################
    graph_1=GraphTK(obj,'Graph #0','1:1',dpi=None,figsize=None)

    graph_1.updateSubPlotParams({'isActive':True,'left':0.16428571428571437,'right':0.8357142857142857,'top':0.8761904761904763,'hspace':0.05714285714285716,'wspace':0.038095238095238126})
    #########################
    # SubGraph '1:1'
    ###############
    ######
    # Axis
    ######
    axis_1_0_0 = graph_1.getAxis('1:1',0)
    ######
    # Curves
    ######
    curve_1 = Curve(zone=odict_keys(['wall1']),varx='x/c',vary='Cp' , line_color='#e41a1c' , marker_face_color='#e41a1c' ,marker_edge_color='#e41a1c' , legend_label='...' )
    graph_1.addCurve('1:1',curve_1)
    ######
    # Grids
    ######
    grid_1_0_0 = graph_1.getGrid('1:1',axis=axis_1_0_0)
    ######
    # Axis settings
    ######
    ######
    # Legends
    ######
    legend_1_0 = graph_1.getLegend('1:1')
