# coding: utf-8
#from __future__ import unicode_literals

# --- Import Section
# Import python
import copy as cpy
import numpy as np
import os
import re
import subprocess
import shlex
from collections import OrderedDict
import imp

try: xrange= range
except: pass

# Import Tkinter
IMPORTOK = True
try: import Tkinter as TK
except: 
    try: import tkinter as TK
    except: IMPORTOK = False
    
try: import ttk as cttk
except: IMPORTOK = False

try:
    # from tkColorChooser import askcolor
    import tkFileDialog
    import tkMessageBox
except:
    try:
        import tkinter.filedialog as tkFileDialog
        import tkinter.messageBox as tkMessageBox
    except: IMPORTOK = False

# Import matplotlib
try:
    from matplotlib.widgets import SubplotTool
    import matplotlib
    # Fit matplotlib usage
    matplotlib.use('TkAgg')
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
#
#    from matplotlib.backends.backend_tkagg import ToolTip
#    from matplotlib.backends.backend_tkagg import ToolTip,cursord
#    from matplotlib.backend_bases import NavigationToolbar2
#    import matplotlib.backends.windowing as windowing
#    from matplotlib.figure import Figure
#
    import matplotlib.pyplot as plt
    # Will be imported in the right movie class:
    # import matplotlib.animation as animation
    pltLeft = plt.rcParams.get('figure.subplot.left')
    pltRight = plt.rcParams.get('figure.subplot.right')
    pltTop = plt.rcParams.get('figure.subplot.top')
    pltBottom = plt.rcParams.get('figure.subplot.bottom')
    pltWSpace = plt.rcParams.get('figure.subplot.wspace')
    pltHSpace = plt.rcParams.get('figure.subplot.hspace')
except:
    IMPORTOK = False
    pltLeft = 0
    pltRight = 0
    pltTop = 0
    pltBottom = 0
    pltWSpace = 0
    pltHSpace = 0
    class NavigationToolbar2TkAgg:
        def __init__(self): return

# Cassiopee import
try:
    import Converter.PyTree as C
    import Converter.Internal as Internal
    import Transform.PyTree as T
    import CPlot.PyTree as CPlot
    import CPlot.Tk as CTK
    import CPlot.Ttk as TTK
except ImportError:
    CTK = None
    TTK = TK

# import CPlot.iconics as iconics

# local widgets list
WIDGETS = {}
VARS = []
GRAPHS = []
UI_addCurve = []
DESKTOP = None
PREVTPZONES = []
STYLEFILE = "style.py"
EXPORTFILE = "fig.png"

#### TODO : link the following variables with preference in GUI of Cassiopee
# local NUM_COLRS for colormap
NUM_COLORS = 8
# local color map set (cf: http://matplotlib.org/examples/color/colormaps_reference.html)
#COLOR_MAP = 'Set2'
# COLOR_MAP = 'rainbow'
COLOR_MAP = 'Set1'
COLOR_LIST_HISTORY = []
default_base = 'Base'
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# TrueFalseDic = {1: True, 0: False}

marker_dic = {
    'none': 'None',
    'plus': '+',
    'star': '*',
    'pixel': ',',
    'point': '.',
    'star3_down': '1',
    'star3_up': '2',
    'star3_left': '3',
    'star3_right': '4',
    'triangle_left': '<',
    'triangle_right': '>',
    'diamond': 'D',
    'hexagon2': 'H',
    'triangle_up': '^',
    'hline': '_',
    'thin_diamond': 'd',
    'hexagon1': 'h',
    'circle': 'o',
    'pentagon': 'p',
    'square': 's',
    'triangle_down': 'v',
    'x': 'x'
}
markername = sorted(marker_dic.keys())

markerlist = [marker_dic[k] for k in markername]

linestylelist = ['solid', 'dashed', 'dashdot', 'dotted', 'None']

t = np.arange(0., 5., 0.002)

data = {
        'Iteration':t,
        'Residual':np.sin(t),
        'Cf':np.sin(t/2),
        'Debit':t*t
        }

# data = {}

#data = { 'tatatiti':{
#                    'Iteration':t,
#                    'Residual':np.sin(t),
#                    'Cf':np.sin(t/2),
##                    'Debit':t*t
#                },
#        'tatatutu':{
#                    'Iteration':-t,
##                    'Residual':np.sin(-t),
#                    'Cf':np.sin(-t/2),
##                    'Debit':t*t
#                },
#        'tititutu':{
##                    'Iteration':t,
##                    'Residual':np.sin(t),
##                    'Cf':np.sin(t/2),
#                    'Debit':t*t
#                },
#        'tatitatu':{
#                    'Iteration':-t,
#                    'Residual':np.sin(-t),
#                    'Cf':np.sin(-t/2),
#                    'Debit':t*t
#                }
#                    }

# ==============================================================================
# ==============================================================================
default_values = {
'Graph':
        {
            'image_background_color':'#FFFFFF', #White
            'image_background_alpha':1.0,
            'subgraph_background_color':'#FFFFFF', #White
            'subgraph_background_alpha':1.0
        },
'Curve':{
            'varx':'',
            'vary':'',
            'line_color':None,
            'line_style':'solid',
            'line_width':2.0,
            'marker_style':'none',
            'marker_size':6.5,
            'marker_edge_width':0.5,
            'marker_face_color':None,
            'marker_edge_color':None,
            'marker_sampling_start':'',
            'marker_sampling_end':'',
            'marker_sampling_step':1,
            'legend_label':'',
            'legend_display': True,
            'ind_axis' : 0,
            'visible':True
        },
'Grid':{
            'Mx_display' : True,
            'Mx_grid_color' : '#000000',
            'Mx_grid_style' : 'dashed',
            'Mx_grid_width' : 1.,
            'Mx_grid_tick_number':5,
            'Mx_grid_tick_size':10.,
            'My_display' : True,
            'My_grid_color' : '#000000',
            'My_grid_style' : 'dashed',
            'My_grid_width' : 1.,
            'My_grid_tick_number':5,
            'My_grid_tick_size':10.,
            'mx_display' : False,
            'mx_grid_color' : '#000000',
            'mx_grid_style' : 'dashed',
            'mx_grid_width' : 1.,
            'mx_grid_tick_number':5,
            'mx_grid_tick_size':10.,
            'my_display' : False,
            'my_grid_color' : '#000000',
            'my_grid_style' : 'dashed',
            'my_grid_width' : 1.,
            'my_grid_tick_number':5,
            'my_grid_tick_size':10.
        },
'Legend':{
            'legend_display' : True,
            'legend_title' : 'Legend',
            'legend_border_width' : 0.5,
            'legend_border_color' : '#000000',
            'vary':'',
            'legend_background_color' :  '#ffffff',
            'legend_background_color_active' : True,
            'legend_position' : 'best',
            'legend_ncol' : 1,
            'legend_label_weight' : 'normal',
            'legend_label_style' : 'normal',
            'legend_label_size' : 12,
            'legend_label_color' : '#000000',
            'legend_title_weight' : 'normal',
            'legend_title_style' : 'normal',
            'legend_title_size' : 8,
            'legend_title_color' : '#000000'
        },
'Axis':{
            'axis_x_logscale' : False,
            'axis_y_logscale' : False,
            'axis_x_autoscale' : True,
            'axis_y_autoscale' : True,
            'axis_x_min' : 0.,
            'axis_x_max' : 1.,
            'axis_y_min' : 0.,
            'axis_y_max' : 1.,
            'axis_x_label' : '',
            'axis_y_label' : '',
            'axis_x_inverted' : False,
            'axis_y_inverted' : False,
            'axis_x_visible' : 1,
            'axis_y_visible' : 1,
            'axis_x_position' : 'both',
            'axis_y_position' : 'both',
            'axis_x_offset' : 0.,
            'axis_y_offset' : 0.,
            'axis_x_label_fontsize' : plt.rcParams.get('font.size'),
            'axis_y_label_fontsize' : plt.rcParams.get('font.size')
        },
'SubPlotParams':{
                    'left'    : pltLeft,
                    'right'   : pltRight,
                    'top'     : pltTop,
                    'bottom'  : pltBottom,
                    'wspace'  : pltWSpace,
                    'hspace'  : pltHSpace,
                    'isActive':True
                },
'TightLayout':{
                    'pad'   : 1.08, # pad can not be None as default value !!!
                    'hpad'  : None,
                    'wpad'  : None
                }
}
################################################################################
## BEGIN COLOR INTEGRATION
import colorsys
################################################################################
# Note :
# ------
# RGB : convention entre 0 et 255 pour chaque composant
# HSV : convention ?
# HTML : #RGB avec R,G,B en hexadecimal
#
# colorsys renvoie entre 0 et 1
#
################################################################################

# class Example(TTK.Frame):
#     def __init__(self, parent,rgb_l):
#         TTK.Frame.__init__(self, parent)
#         f = HFrame(parent,rgb_l)
#         f.pack(fill="both", expand=True)
################################################################################

################################################################################
class ColorControler(object):
    def __init__(self,parent,color="#3a65c6",colormapName='Set1',nbColorMap=8,colorHistoryList=None ):
        self.parent = parent
        self.color = color
        #
        self.Nstep = 1530
        self.turnOnFlag2ChangeColor()
        # determine rgb value of self.color
        hexa = self.color.lstrip('#')
        (r,g,b) = tuple(int(hexa[i:i+2], 16) for i in (0, 2 ,4))
        # print(r,g,b)
        # (r,g,b) = self.winfo_rgb(self.color)
        r = float(r)/255.
        g = float(g)/255.
        b = float(b)/255.
        (h,s,v) = colorsys.rgb_to_hsv(r, g, b)
        self.hvalue = h
        self.svalue = s
        self.vvalue = v
        # h goes from 0->360
        hmin = 0.
        hmax = 360.
        # Step Value
        step = (hmax - hmin)/self.Nstep
        # Creation of the table of h
        h_l = [i*step for i in range(self.Nstep+1)]
        # Getting the rgb values according to h (rgb 0->255)
        self.rgb_l = [self.discretizedColor(h) for h in h_l]
        # Getting the colormap HTML code liste
        self.colormapList = self.createColorMap(colormapName,nbColorMap)
        # Controler frame positionning
        self.controlerFrame = TTK.Frame(self.parent.frame)
        self.controlerFrame.rowconfigure(0,weight=1)
        self.controlerFrame.rowconfigure(1,weight=0)
        self.controlerFrame.rowconfigure(2,weight=0)
        self.controlerFrame.rowconfigure(3,weight=0)
        self.controlerFrame.rowconfigure(4,weight=0)
        self.controlerFrame.columnconfigure(0,weight=0)
        self.controlerFrame.columnconfigure(1,weight=1)
        self.controlerFrame.grid(row=0,column=0,sticky="NSEW")
        # H-frame
        self.hframe = HFrame(self.controlerFrame,self,self.rgb_l,self.hvalue)
        self.hframe.grid(row=0,column=0,sticky="NSEW")
        # SV-frame
        self.svframe = SVFrame(self.controlerFrame,self,self.hvalue,self.svalue,self.vvalue)
        self.svframe.grid(row=0,column=1,sticky="NSEW")
        #################
        # 2ND & 3RD LINES
        # Old Color Frame
        self.oldcolorFrame = ColorFrame(self.controlerFrame,self,self.hvalue,self.svalue,self.vvalue)
        self.oldcolorFrame.grid(row=1,column=0,sticky="EW")
        # New Color Frame
        self.newcolorFrame = ColorFrame(self.controlerFrame,self,self.hvalue,self.svalue,self.vvalue)
        self.newcolorFrame.grid(row=2,column=0,sticky="EW")
        # subframeOld
        subframeOld = TTK.Frame(self.controlerFrame)
        subframeOld.rowconfigure(0,weight=0)
        subframeOld.rowconfigure(1,weight=0)
        subframeOld.columnconfigure(0,weight=0)
        subframeOld.columnconfigure(1,weight=1)
        subframeOld.columnconfigure(2,weight=1)
        subframeOld.columnconfigure(3,weight=1)
        subframeOld.columnconfigure(4,weight=1)
        subframeOld.grid(row=1,column=1,rowspan=2,sticky="EW")
        # Old Label
        oldLabel = TTK.Label(subframeOld,text='Previous color')
        oldLabel.grid(row=0,column=0,sticky='W')
        # # subframeNew
        # subframeNew = TTK.Frame(self.controlerFrame)
        # subframeNew.rowconfigure(0,weight=0)
        # subframeNew.rowconfigure(1,weight=0)
        # subframeNew.columnconfigure(0,weight=0)
        # subframeNew.columnconfigure(1,weight=1)
        # subframeNew.columnconfigure(2,weight=1)
        # subframeNew.columnconfigure(3,weight=1)
        # subframeNew.columnconfigure(4,weight=1)
        # subframeNew.grid(row=2,column=1,sticky="EW")
        # New Label
        newLabel = TTK.Label(subframeOld,text='New color')
        newLabel.grid(row=1,column=0,sticky='W')
        # RGB OLD
        (r,g,b)=colorsys.hsv_to_rgb(self.hvalue,self.svalue,self.vvalue)
        ### R
        rlblframeOld = TTK.LabelFrame(subframeOld,text='R')
        rlblframeOld.rowconfigure(0,weight=1)
        rlblframeOld.columnconfigure(0,weight=1)
        rlblframeOld.grid(row=0,column=1,sticky='EW')
        r_old = TK.StringVar()
        r_old.set('%d'%(int(round(255.*r))))
        rEntry = TK.Entry(rlblframeOld,state=TK.DISABLED,textvariable=r_old)
        rEntry.grid(row=0,column=0,sticky="NSEW")
        ### G
        glblframeOld = TTK.LabelFrame(subframeOld,text='G')
        glblframeOld.rowconfigure(0,weight=1)
        glblframeOld.columnconfigure(0,weight=1)
        glblframeOld.grid(row=0,column=2,sticky='EW')
        g_old = TK.StringVar()
        g_old.set('%d'%(int(round(255.*g))))
        gEntry = TK.Entry(glblframeOld,state=TK.DISABLED,textvariable=g_old)
        gEntry.grid(row=0,column=0,sticky="NSEW")
        ### B
        blblframeOld = TTK.LabelFrame(subframeOld,text='B')
        blblframeOld.rowconfigure(0,weight=1)
        blblframeOld.columnconfigure(0,weight=1)
        blblframeOld.grid(row=0,column=3,sticky='EW')
        b_old = TK.StringVar()
        b_old.set('%d'%(int(round(255.*b))))
        bEntry = TK.Entry(blblframeOld,state=TK.DISABLED,textvariable=b_old)
        bEntry.grid(row=0,column=0,sticky="NSEW")
        # HTML old
        nr = int(r*255)
        ng = int(g*255)
        nb = int(b*255)
        #
        html_color = '#%02x%02x%02x' % (nr,ng,nb)
        #
        htmllblframeOld = TTK.LabelFrame(subframeOld,text='HTML')
        htmllblframeOld.rowconfigure(0,weight=1)
        htmllblframeOld.columnconfigure(0,weight=1)
        htmllblframeOld.grid(row=0,column=4,sticky='EW')
        html_old = TK.StringVar()
        html_old.set('%s'%html_color)
        htmlEntry = TK.Entry(htmllblframeOld,state=TK.DISABLED,textvariable=html_old)
        htmlEntry.grid(row=0,column=0,sticky="NSEW")

        # RGB NEW
        (r,g,b)=colorsys.hsv_to_rgb(self.svframe.hvalue,self.svframe.svalue,self.svframe.vvalue)
        ### R
        rlblframeNew = TTK.LabelFrame(subframeOld,text='R')
        rlblframeNew.rowconfigure(0,weight=1)
        rlblframeNew.columnconfigure(0,weight=1)
        rlblframeNew.grid(row=1,column=1,sticky='EW')
        self.r_new = TK.StringVar()
        self.r_new.set('%d'%(int(round(255.*r))))
        self.r_new.trace('w',lambda nm, idx, mode,var='rgb': self.event_editNewColorEntry(var))
        rEntry = TK.Entry(rlblframeNew,textvariable=self.r_new)
        rEntry.grid(row=0,column=0,sticky="NSEW")
        ### G
        glblframeNew = TTK.LabelFrame(subframeOld,text='G')
        glblframeNew.rowconfigure(0,weight=1)
        glblframeNew.columnconfigure(0,weight=1)
        glblframeNew.grid(row=1,column=2,sticky='EW')
        self.g_new = TK.StringVar()
        self.g_new.set('%d'%(int(round(255.*g))))
        self.g_new.trace('w',lambda nm, idx, mode,var='rgb': self.event_editNewColorEntry(var))
        gEntry = TK.Entry(glblframeNew,textvariable=self.g_new)
        gEntry.grid(row=0,column=0,sticky="NSEW")
        ### B
        blblframeNew = TTK.LabelFrame(subframeOld,text='B')
        blblframeNew.rowconfigure(0,weight=1)
        blblframeNew.columnconfigure(0,weight=1)
        blblframeNew.grid(row=1,column=3,sticky='EW')
        self.b_new = TK.StringVar()
        self.b_new.set('%d'%(int(round(255.*b))))
        self.b_new.trace('w',lambda nm, idx, mode,var='rgb': self.event_editNewColorEntry(var))
        bEntry = TK.Entry(blblframeNew,textvariable=self.b_new)
        bEntry.grid(row=0,column=0,sticky="NSEW")
        # HTML new
        nr = int(r*255)
        ng = int(g*255)
        nb = int(b*255)
        #
        html_color = '#%02x%02x%02x' % (nr,ng,nb)
        #
        htmllblframeNew = TTK.LabelFrame(subframeOld,text='HTML')
        htmllblframeNew.rowconfigure(0,weight=1)
        htmllblframeNew.columnconfigure(0,weight=1)
        htmllblframeNew.grid(row=1,column=4,sticky='EW')
        self.html_new = TK.StringVar()
        self.html_new.set('%s'%html_color)
        self.html_new.trace('w',lambda nm, idx, mode,var='html': self.event_editNewColorEntry(var))
        htmlEntry = TK.Entry(htmllblframeNew,textvariable=self.html_new)
        htmlEntry.grid(row=0,column=0,sticky="NSEW")

        ##########
        # 4TH Line
        pickUpColorFrame = TTK.Frame(self.controlerFrame)
        pickUpColorFrame.rowconfigure(0,weight=1)
        pickUpColorFrame.columnconfigure(0,weight=1)
        pickUpColorFrame.columnconfigure(1,weight=1)
        pickUpColorFrame.grid(row=3,column=0,columnspan=2,sticky="NSEW")
        #
        classicalColorLblFrame = TTK.LabelFrame(pickUpColorFrame,text="Classical colors")
        classicalColorLblFrame.rowconfigure(0,weight=1)
        for ind in range(nbColorMap):
            classicalColorLblFrame.columnconfigure(ind,weight=1)
        classicalColorLblFrame.grid(row=0,column=0,sticky="NSEW")
        #
        for ind in range(len(self.colormapList)):
            htmlcolor = self.colormapList[ind]
            hexa = htmlcolor.lstrip('#')
            (r,g,b) = tuple(int(hexa[i:i+2], 16) for i in (0, 2 ,4))
            (h,s,v) = colorsys.rgb_to_hsv(float(r)/255.,float(g)/255.,float(b)/255.)
            squareColor = SquareColor(classicalColorLblFrame,self,h,s,v)
            squareColor.grid(row=0,column=ind,sticky='NS')
        #
        lastUsedColorLblFrame = TTK.LabelFrame(pickUpColorFrame,text="Last colors")
        lastUsedColorLblFrame.rowconfigure(0,weight=1)
        for ind in range(len(colorHistoryList)):
            lastUsedColorLblFrame.columnconfigure(ind,weight=1)
        lastUsedColorLblFrame.grid(row=0,column=1,sticky="NSEW")
        #
        for ind in range(len(colorHistoryList)):
            htmlcolor = colorHistoryList[ind]
            hexa = htmlcolor.lstrip('#')
            (r,g,b) = tuple(int(hexa[i:i+2], 16) for i in (0, 2 ,4))
            (h,s,v) = colorsys.rgb_to_hsv(float(r)/255.,float(g)/255.,float(b)/255.)
            squareColor = SquareColor(lastUsedColorLblFrame,self,h,s,v)
            squareColor.grid(row=0,column=ind,sticky='NS')
        ##########
        # 5TH LINE
        buttonFrame = TTK.Frame(self.controlerFrame)
        buttonFrame.rowconfigure(0,weight=0)
        buttonFrame.columnconfigure(0,weight=1)
        buttonFrame.columnconfigure(1,weight=1)
        buttonFrame.grid(row=4,column=0,columnspan=2,sticky='EW')
        #
        OKbutton = TTK.Button(buttonFrame,text="Select",command=self.OKbutton)
        OKbutton.grid(row=0,column=0)
        cancelbutton = TTK.Button(buttonFrame,text="Cancel",command=self.cancelButton)
        cancelbutton.grid(row=0,column=1)

    def discretizedColor(self,h):
        # h goes from 0->360
        x = float(h)*1530./360.
        r = self.getRvalue(x)
        g = self.getGvalue(x)
        b = self.getBvalue(x)
        return (r,g,b)

    def getRvalue(self,h):
        if h<255:
            res = 255
        elif h<510:
            res = 510 - h
        elif h<1020:
            res = 0
        elif h<1275:
            res = h-1020
        else:
            res = 255
        return res

    def getGvalue(self,h):
        if h<255:
            res = h
        elif h<765:
            res = 255
        elif h<1020:
            res = 1020 - h
        else:
            res = 0
        return res

    def getBvalue(self,h):
        if h<510:
            res = 0
        elif h<765:
            res = h - 510
        elif h<1275:
            res = 255
        else:
            res = 1530 - h
        return res

    def cancelButton(self):
        self.parent.cmd_close()

    def OKbutton(self):
        self.parent.selectColor(self.html_new.get())

    def changeColor(self,h,s,v):
        ### Turn off the flag => this will deactivate the event while the entry of rgb & html is changed
        self.turnOffFlag2ChangeColor()
        # Update hsv values
        self.updateHSV4Canvas(h=h,s=s,v=v)
        ### Turn on the flag => this will reactivate the event while the entry of rgb & html is changed
        self.turnOnFlag2ChangeColor()
        ### Update the color of the new color box
        self.updateNewColor()
        ### Update the hline position
        self.redrawHLine()
        ### Redraw the SV space
        self.redrawSVFrame()

    def createColorMap(self,colormapName,nbColorMap):
        colormapList = []
        cm = plt.get_cmap(colormapName)
        for ind in range(nbColorMap):
            rgb = cm(1.*(ind%nbColorMap)/nbColorMap)
            html = '#%02x%02x%02x' % (int(255*rgb[0]),int(255*rgb[1]),int(255*rgb[2]))
            colormapList.append(html)
        return colormapList

    def turnOnFlag2ChangeColor(self):
        self.flag2ChangeColor = True

    def turnOffFlag2ChangeColor(self):
        self.flag2ChangeColor = False

    def event_editNewColorEntry(self,var,*args):
        if not self.flag2ChangeColor:
            return
        if var == 'rgb':
            # Get R val
            r = self.r_new.get()
            g = self.g_new.get()
            b = self.b_new.get()

            re_pattern = re.compile('''^[0-9]*$''')

            # Check the input value (should be integer between 0 and 255)
            for val in [r,g,b]:
                val = val.strip()
                if not re.match(re_pattern,val):
                    return
                if val =='':
                    # val = 0
                    return
                val = int(val)
                if val<0 or val>255:
                    return
        elif var == 'html':
            html = self.html_new.get()

            re_pattern = re_pattern = re.compile('''^\#[0-9A-Fa-f]{6}$''')
            # Check the input to be an html color code #AABBCC where AA,BB,CC are hexadecimal number between 0 and 255
            if not re.match(re_pattern,html):
                return
            # Convert it to rgb
            hexa = html.lstrip('#')
            (r,g,b) = tuple(int(hexa[i:i+2], 16) for i in (0, 2 ,4))
        else:
            return

        # Convert to hsv
        (h,s,v)=colorsys.rgb_to_hsv(float(r)/255.,float(g)/255.,float(b)/255.)
        # Update hsv without editing the text box !
        ### Turn off the flag => this will deactivate the event while the entry of rgb & html is changed
        self.turnOffFlag2ChangeColor()
        ### Update HSV and then rgb & html in display
        self.updateHSV4Canvas(h=h,s=s,v=v)
        ### Turn on the flag => this will reactivate the event while the entry of rgb & html is changed
        self.turnOnFlag2ChangeColor()
        ### Update the color of the new color box
        self.updateNewColor()
        ### Update the hline position
        self.redrawHLine()
        ### Redraw the SV space
        self.redrawSVFrame()

    def updateHSV4Canvas(self,h=None,s=None,v=None):
        if h is not None:
            self.hframe.updateHSV(h=h)
            self.svframe.updateHSV(h=h)
            self.newcolorFrame.updateHSV(h=h)
        if s is not None:
            self.hframe.updateHSV(s=s)
            self.svframe.updateHSV(s=s)
            self.newcolorFrame.updateHSV(s=s)
        if v is not None:
            self.hframe.updateHSV(v=v)
            self.svframe.updateHSV(v=v)
            self.newcolorFrame.updateHSV(v=v)
        #
        self.updateNewColorInfo()
    def updateNewColorInfo(self):
        # if not self.flag2ChangeColor:
        #     return

        (r,g,b) = colorsys.hsv_to_rgb(self.svframe.hvalue,self.svframe.svalue,self.svframe.vvalue)

        self.r_new.set('%d'%(int(round(255.*r))))
        self.g_new.set('%d'%(int(round(255.*g))))
        self.b_new.set('%d'%(int(round(255.*b))))
        nr = int(round(r*255))
        ng = int(round(g*255))
        nb = int(round(b*255))
        #
        self.html_new.set('#%02x%02x%02x' % (nr,ng,nb))
    def updateNewColor(self):
        self.newcolorFrame.update()

    def redrawSVFrame(self):
        self.svframe._draw_all()

    def redrawHLine(self):
        self.hframe._draw_line()


class SVFrame(TK.Canvas):
    '''A gradient frame which uses a canvas to draw the background'''
    def __init__(self, parent,controler,hvalue,svalue,vvalue, borderwidth=1, relief="sunken"):
        self.controler = controler
        self.parent = parent
        self.hvalue = hvalue
        self.svalue = svalue
        self.vvalue = vvalue
        TK.Canvas.__init__(self, parent, borderwidth=borderwidth, relief=relief)
        self.bind("<Configure>", self._draw_all)
        self.bind("<Button-1>",self.left_bt_click)

    def updateHSV(self,h=None,s=None,v=None):
        if h is not None:
            self.hvalue = h
        if s is not None:
            self.svalue = s
        if v is not None:
            self.vvalue = v

    def left_bt_click(self,event):
        # Get the size of the window
        width = self.winfo_width()
        height = self.winfo_height()
        # Update the values of s and v according to the click event
        self.controler.turnOffFlag2ChangeColor()
        self.controler.updateHSV4Canvas(v=float(event.x)/float(width),s=float(event.y)/float(height))
        self.controler.turnOnFlag2ChangeColor()
        #
        self._draw_lines()
        #
        self.controler.updateNewColor()

    def _draw_all(self,event=None):
        self._draw_gradient(event)
        self._draw_lines(event)

    def _draw_lines(self,event=None):
        nValue = 50
        # Delete previous lines before creating new ones
        if self.find_withtag('horizontal2'):
            self.delete(self.find_withtag('horizontal2'))
        if self.find_withtag('vertical2'):
            self.delete(self.find_withtag('vertical2'))
        # Get the size of the wdw info
        width = self.winfo_width()
        height = self.winfo_height()
        # Compute the position values
        s_line = int(self.svalue * height)
        v_line = int(self.vvalue * width)
        # Create new lines
        self.create_line(v_line,0,v_line,height, tags=("horizontal2",), fill="#000000")
        self.create_line(0,s_line,width,s_line, tags=("vertical2",), fill="#000000")




    def _draw_gradient(self, event=None):
        '''Draw the gradient'''
        if self.find_withtag('gradient2'):
            self.delete(self.find_withtag('gradient2'))
        #
        nValue = 50
        #
        width = self.winfo_width()
        height = self.winfo_height()
        #
        limit = height
        i = 0.
        j = 0.
        height_step = float(height)/float(nValue)
        width_step = float(width)/float(nValue)
        while i<height:
            j = 0.
            while j<width:
                s = float(i)/float(height)
                v = float(j)/float(width)
                h = self.hvalue
                (r,g,b) = colorsys.hsv_to_rgb(h, s, v)
                nr = int(round(r*255.))
                ng = int(round(g*255.))
                nb = int(round(b*255.))
                #
                color = "#%02x%02x%02x" % (nr,ng,nb)
                self.create_rectangle(j,i,j+width_step,i+height_step, tags=("gradient2",), fill=color,width=0)
                #
                j = j + width_step
            i = i + height_step
        self.lower("gradient2")

class HFrame(TK.Canvas):
    '''A gradient frame which uses a canvas to draw the background'''
    def __init__(self, parent,controler,rgb_l,hvalue, borderwidth=1, relief="sunken"):
        TK.Canvas.__init__(self, parent, width=40,borderwidth=borderwidth, relief=relief)
        self.controler = controler
        self.parent = parent
        self.hvalue = hvalue
        self.bind("<Configure>", self._draw_all)
        self.rgb_l = rgb_l
        self.bind("<Button-1>",self.left_bt_click)

    def updateHSV(self,h=None,s=None,v=None):
        if h is not None:
            self.hvalue = h
        if s is not None:
            self.svalue = s
        if v is not None:
            self.vvalue = v


    def left_bt_click(self,event):
        # Get the size of the window
        width = 40
        height = self.winfo_height()
        # Update the values of s and v according to the click event
        self.controler.turnOffFlag2ChangeColor()
        self.controler.updateHSV4Canvas(h=float(event.y)/float(height))
        self.controler.turnOnFlag2ChangeColor()
        #
        self._draw_line()
        #
        self.controler.redrawSVFrame()
        #
        self.controler.updateNewColor()

    def _draw_all(self,event=None):
        self._draw_gradient(event)
        self._draw_line(event)

    def _draw_line(self,event=None):
        nValue = 50
        # Delete previous lines before creating new ones
        if self.find_withtag('horizontal'):
            self.delete(self.find_withtag('horizontal'))
        # Get the size of the wdw info
        width = 40
        height = self.winfo_height()
        # Compute the position values
        h_line = int(self.hvalue * height)
        # Create new lines
        self.create_line(0,h_line,width,h_line, tags=("horizontal",), fill="#000000")


    def _draw_gradient(self, event=None):
        '''Draw the gradient'''
        if self.find_withtag('gradient'):
            self.delete(self.find_withtag('gradient'))
        width = 40
        height = self.winfo_height()
        limit = height

        for i in range(limit):
            nbColor = self.controler.Nstep
            ii = nbColor*i/limit

            nr = int(self.rgb_l[ii][0])
            ng = int(self.rgb_l[ii][1])
            nb = int(self.rgb_l[ii][2])
            #
            color = "#%02x%02x%02x" % (nr,ng,nb)
            self.create_line(0,i,width,i, tags=("gradient",), fill=color)
        self.lower("gradient")

class ColorFrame(TK.Canvas):
    '''A gradient frame which uses a canvas to draw the background'''
    def __init__(self, parent,controler,hvalue,svalue,vvalue, borderwidth=1, relief="sunken"):
        TK.Canvas.__init__(self, parent, width=40,height=40,borderwidth=borderwidth, relief=relief)
        self.controler = controler
        self.parent = parent
        self.hvalue = hvalue
        self.svalue = svalue
        self.vvalue = vvalue
        self.bind("<Configure>", self._draw_all)

    def update(self):
        self._draw_all()

    def updateHSV(self,h=None,s=None,v=None):
        if h is not None:
            self.hvalue = h
        if s is not None:
            self.svalue = s
        if v is not None:
            self.vvalue = v
        self._draw_all()

    def _draw_all(self,event=None):
        self._draw_color(event)

    def _draw_color(self, event=None):
        '''Draw the color'''
        if self.find_withtag('color'):
            self.delete(self.find_withtag('color'))
        # width = self.winfo_width()
        width = 80
        height = 40

        (r,g,b) = colorsys.hsv_to_rgb(self.hvalue,self.svalue,self.vvalue)

        nr = int(round(r*255))
        ng = int(round(g*255))
        nb = int(round(b*255))
        color = "#%02x%02x%02x" % (nr,ng,nb)

        self.create_rectangle(0,0,width,height, tags=("color",), fill=color, width=0)

        self.lower("color")

class SquareColor(TK.Canvas):
    '''A gradient frame which uses a canvas to draw the background'''
    def __init__(self, parent,controler,hvalue,svalue,vvalue, borderwidth=1, relief="sunken"):
        self.controler = controler
        self.parent = parent
        self.hvalue = hvalue
        self.svalue = svalue
        self.vvalue = vvalue
        TK.Canvas.__init__(self, parent, borderwidth=borderwidth, relief=relief,width=15,height=15,)
        self.bind("<Configure>", self._draw_all)
        self.bind("<Button-1>",self.left_bt_click)

    def left_bt_click(self,event):
        self.controler.changeColor(self.hvalue,self.svalue,self.vvalue)

    def _draw_all(self,event=None):
        self._draw_gradient(event)

    def _draw_gradient(self, event=None):
        '''Draw the gradient'''
        if self.find_withtag('colorsquare'):
            self.delete(self.find_withtag('colorsquare'))
        #
        width = self.winfo_width()
        height = self.winfo_height()

        (r,g,b) = colorsys.hsv_to_rgb(self.hvalue, self.svalue, self.vvalue)
        nr = int(round(r*255.))
        ng = int(round(g*255.))
        nb = int(round(b*255.))
        #
        color = "#%02x%02x%02x" % (nr,ng,nb)
        self.create_rectangle(0,0,width,height, tags=("colorsquare",), fill=color,width=0)
        self.lower("colorsquare")






class color_dialogWindow(TK.Toplevel):
    # --------------------------------------------------------------------- init
    def __init__(self,parent,B,color,colorHistoryList=COLOR_LIST_HISTORY,colormapName=COLOR_MAP,nbColorMap=NUM_COLORS,extra_data=None):
        TK.Toplevel.__init__(self)
        (x,y) = self.winfo_pointerxy()
        # self.geometry('+%s+%s'%(x,y))
        self.geometry('%dx%d%+d%+d'%(375,375,x,y))
        self.parent = parent
        self.title('Select color')
        #self.protocol("WM_DELETE_WINDOW",self.cmd_close)
        self.extra_data = extra_data
        #
        self.parent = parent
        self.B = B
        #
        self.grid_columnconfigure(0,weight=1)
        self.grid_rowconfigure(0,weight=1)
        #
        self.frame = TTK.Frame(self)
        self.frame.grid(row=0,column=0,sticky='NESW')
        self.frame.rowconfigure(0,weight=1)
        self.frame.columnconfigure(0,weight=1)
        #
        history = list(COLOR_LIST_HISTORY)
        history.reverse()
        colorControler = ColorControler(self,color=color,colorHistoryList=history[:NUM_COLORS],colormapName=colormapName,nbColorMap=nbColorMap)

    def cmd_close(self,event=None):
        self.destroy()

    def selectColor(self,color):
        # Propagate color to B
        self.parent.updateColor(color,self.B,self.extra_data)
        COLOR_LIST_HISTORY.append(color)
        # Close
        self.cmd_close()




## END COLOR INTEGRATION
################################################################################

# ==============================================================================
# ==============================================================================
class editAxisWindow(TK.Toplevel):
    # --------------------------------------------------------------------- init
    def __init__(self):
        TK.Toplevel.__init__(self)
        self.title('Set axis')
        self.protocol("WM_DELETE_WINDOW",self.cmd_close)
        #
        self.list_dialog = None
        self.input_dialog = None
    # --------------------------------------------------------------- initialize
    def initialize(self,parent):
        self.parent = parent
        #
        self.graph = self.parent.activeGraph.val
        self.zone  = self.parent.position.val
        self.ind_axis = 0
        # self.axis_to_twin = self.ind_axis
        try:
            self.subGraph = self.parent.graphWdwL[self.graph].fig.subGraph[self.zone]
        except TypeError: # Exit if no graph & zone available
            self.cmd_close()
            return
#        except IndexError: # Error while closing window ...
#            return
        #
        self.grid_columnconfigure(0,weight=1)
        self.grid_rowconfigure(0,weight=1)
        #
        self.createFrame()
    # -------------------------------------------------------------- createFrame
    def createFrame(self):
        # Close all opened dialog windows
        self.closeAllDialog()
        #
        self.frame = TTK.Frame(self)
        self.frame.grid(row=0,column=0,sticky='NESW')
        #
        self.frame.grid_columnconfigure(0,weight=1)
        self.frame.grid_rowconfigure(0,weight=1)
        self.frame.grid_rowconfigure(1,weight=1)
        self.frame.grid_rowconfigure(2,weight=1)
        #

        ##########################################
        ###### Select axis line
        selectAxisFrame = TTK.Frame(self.frame)
        selectAxisFrame.grid(row=0, column=0, sticky='NSEW')
        selectAxisFrame.grid_columnconfigure(0, weight=1)
        selectAxisFrame.grid_columnconfigure(1, weight=1)
        selectAxisFrame.grid_columnconfigure(2, weight=1)
        # selectAxisFrame.grid_columnconfigure(3, weight=1)
        selectAxisFrame.grid_rowconfigure(0, weight=1)
        #
        # Select method to create axis
        newAxesMethodLblFrame = TTK.LabelFrame(selectAxisFrame,text='Type of new axis')
        newAxesMethodLblFrame.grid_columnconfigure(0, weight=1)
        newAxesMethodLblFrame.grid_rowconfigure(0, weight=1)
        newAxesMethodLblFrame.grid(row=0, column=0, sticky='NSEW')
        #
        self.newAxesMethod = []
        B = TTK.Button(newAxesMethodLblFrame,text='New',command=lambda n=(0,self.newAxesMethod): self.bt_click(n))
        self.newAxesMethod.append(B)
        B.list = ['New','Twin X','Twin Y']
        B.val = 'New'
        B.var = ['']
        B.treatmentId = 0 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        # Create axis button
        B = TTK.Button(selectAxisFrame,text='Add new axis',command=self.cmd_addAxis)
        B.grid(row=0, column=1, sticky='NSEW')
        #
        selectAxisLblFrame = TTK.LabelFrame(selectAxisFrame, text='Select Axis')
        selectAxisLblFrame.grid_columnconfigure(0, weight=1)
        selectAxisLblFrame.grid_rowconfigure(0, weight=1)
        selectAxisLblFrame.grid(row=0, column=2, sticky='NSEW')
        #
        self.addAxisItem = []
        B = TTK.Button(selectAxisLblFrame,text=self.ind_axis,command=lambda n=(0,self.addAxisItem): self.bt_click(n))
        self.addAxisItem.append(B)
        B.list = [i for i in range(len(self.subGraph.axis))]
        B.val = self.ind_axis
        B.var = ['ind_axis']
        B.treatmentId = 0 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        ##########################################
        ###### X axis
        xlblframe = TTK.LabelFrame(self.frame,text="X Axis")
        xlblframe.grid(row=1,column=0,sticky='NSEW')
        xlblframe.grid_columnconfigure(0,weight=1)
        xlblframe.grid_columnconfigure(1,weight=1)
        xlblframe.grid_columnconfigure(2,weight=1)
        xlblframe.grid_columnconfigure(3,weight=1)
        xlblframe.grid_columnconfigure(4,weight=1)
        xlblframe.grid_columnconfigure(5,weight=1)
        xlblframe.grid_columnconfigure(6,weight=1)
        xlblframe.grid_columnconfigure(7,weight=1)
        xlblframe.grid_columnconfigure(8,weight=1)
        xlblframe.grid_columnconfigure(9,weight=1)
        xlblframe.grid_rowconfigure(0,weight=1)
        ## ## --> Visible
        visibleframe = TTK.LabelFrame(xlblframe,text="Visible")
        visibleframe.grid(row=0,column=0,sticky='NSEW')
        visibleframe.grid_columnconfigure(0,weight=1)
        visibleframe.grid_rowconfigure(0,weight=1)
        #
        self.visibleItem = []
        var = TK.IntVar()
        var.set(self.subGraph.axis_property[self.ind_axis].x.axis_visible)
        CB = TTK.Checkbutton(visibleframe,variable=var,command=lambda n=0: self.cb_visible(n))#, variable=var)
        CB.val = var
        CB.var = ['x','axis_visible']
        CB.grid(row=0,column=0,sticky='NSEW')
        self.visibleItem.append(CB)
        ## ## --> Use logscale
        logframe = TTK.LabelFrame(xlblframe,text="Use Log Scale")
        logframe.grid(row=0,column=1,sticky='NSEW')
        logframe.grid_columnconfigure(0,weight=1)
        logframe.grid_rowconfigure(0,weight=1)
        #
        self.logscaleItem = []
        var = TK.IntVar()
        var.set(self.subGraph.axis_property[self.ind_axis].x.axis_logscale)
        CB = TTK.Checkbutton(logframe,variable=var,command=lambda n=0: self.cb_logscale(n))#, variable=var)
        CB.val = var
        CB.var = ['x','axis_logscale']
        CB.grid(row=0,column=0,sticky='NSEW')
        self.logscaleItem.append(CB)
        ## ## --> Use autoscale
        autoscaleframe = TTK.LabelFrame(xlblframe,text="Use Auto Scale")
        autoscaleframe.grid(row=0,column=2,sticky='NSEW')
        autoscaleframe.grid_columnconfigure(0,weight=1)
        autoscaleframe.grid_rowconfigure(0,weight=1)
        #
        self.autoscaleItem = []
        var = TK.IntVar()
        var.set(self.subGraph.axis_property[self.ind_axis].x.axis_autoscale)
        CB = TTK.Checkbutton(autoscaleframe,variable=var,command=lambda n=0: self.cb_autoscale(n))#, variable=var)
        CB.val = var
        CB.var = ['x','axis_autoscale']
        CB.grid(row=0,column=0,sticky='NSEW')
        self.autoscaleItem.append(CB)
        ## ## --> Xmin
        xmin = self.subGraph.axis[self.ind_axis].get_xlim()[0]
        xminframe = TTK.LabelFrame(xlblframe, text="Min")
        xminframe.grid(row=0,column=3,sticky='NESW')
        #
        xminframe.grid_columnconfigure(0,weight=1)
        xminframe.grid_rowconfigure(0,weight=1)
        #
        self.positionItem = []
        B = TTK.Button(xminframe,text=xmin,command=lambda n=(0,self.positionItem): self.bt_click(n))
        B.list = []
        B.val = xmin
        B.var = ['x','axis_min']
        B.treatmentId = 4 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.positionItem.append(B)
        ## ## --> Xmax
        xmax = self.subGraph.axis[self.ind_axis].get_xlim()[1]
        xmaxframe = TTK.LabelFrame(xlblframe, text="Max")
        xmaxframe.grid(row=0,column=4,sticky='NESW')
        #
        xmaxframe.grid_columnconfigure(0,weight=1)
        xmaxframe.grid_rowconfigure(0,weight=1)
        #
        B = TTK.Button(xmaxframe,text=xmax,command=lambda n=(1,self.positionItem): self.bt_click(n))
        B.list = []
        B.val = xmax
        B.var = ['x','axis_max']
        B.treatmentId = 4 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.positionItem.append(B)
        ## ## --> Invert axis
        invertedframe = TTK.LabelFrame(xlblframe,text="Invert Axis")
        invertedframe.grid(row=0,column=5,sticky='NSEW')
        invertedframe.grid_columnconfigure(0,weight=1)
        invertedframe.grid_rowconfigure(0,weight=1)
        #
        self.invertItem = []
        var = TK.IntVar()
        var.set(self.subGraph.axis_property[self.ind_axis].x.axis_inverted)
        CB = TTK.Checkbutton(invertedframe,variable=var,command=lambda n=0: self.cb_invertedAxis(n))#, variable=var)
        CB.val = var
        CB.var = ['x','axis_inverted']
        CB.grid(row=0,column=0,sticky='NSEW')
        self.invertItem.append(CB)
        ## ## --> Label name
        labelframe = TTK.LabelFrame(xlblframe, text="Label name")
        labelframe.grid(row=0,column=6,sticky='NESW')
        #
        labelframe.grid_columnconfigure(0,weight=1)
        labelframe.grid_rowconfigure(0,weight=1)
        #
        self.labelItem = []
        B = TTK.Button(labelframe,text=self.subGraph.axis_property[self.ind_axis].x.axis_label,command=lambda n=(0,self.labelItem): self.bt_click(n))
        B.list = []
        B.val = self.subGraph.axis_property[self.ind_axis].x.axis_label
        B.var = ['x','axis_label']
        B.treatmentId = 3
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.labelItem.append(B)
        ## ## --> Label size
        labelframe = TTK.LabelFrame(xlblframe, text="Label size")
        labelframe.grid(row=0,column=7,sticky='NESW')
        #
        labelframe.grid_columnconfigure(0,weight=1)
        labelframe.grid_rowconfigure(0,weight=1)
        #
        self.labelSizeItem=[]
        B = TTK.Button(labelframe,text=self.subGraph.axis_property[self.ind_axis].x.axis_label_fontsize,command=lambda n=(0,self.labelSizeItem): self.bt_click(n))
        B.list = []
        B.val = self.subGraph.axis_property[self.ind_axis].x.axis_label_fontsize
        B.var = ['x','axis_label_fontsize']
        B.treatmentId = 4 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.labelSizeItem.append(B)
        ## ## --> axis Position
        position = self.subGraph.axis_property[self.ind_axis].x.axis_position
        positionframe = TTK.LabelFrame(xlblframe, text="Axis position")
        positionframe.grid(row=0,column=8,sticky='NESW')
        #
        positionframe.grid_columnconfigure(0,weight=1)
        positionframe.grid_rowconfigure(0,weight=1)
        #
        self.axisPositionItem=[]
        B = TTK.Button(positionframe,text=position,command=lambda n=(0,self.axisPositionItem): self.bt_click(n))
        B.list = ['top','bottom','both']
        B.val = position
        B.var = ['x','axis_position']
        B.treatmentId = 0 #
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.axisPositionItem.append(B)
        ## ## --> Axis offset
        self.offsetItem=[]
        offset = self.subGraph.axis_property[self.ind_axis].x.axis_offset
        offsetframe = TTK.LabelFrame(xlblframe, text="Axis offset")
        offsetframe.grid(row=0,column=9,sticky='NESW')
        #
        offsetframe.grid_columnconfigure(0,weight=1)
        offsetframe.grid_rowconfigure(0,weight=1)
        #
        B = TTK.Button(offsetframe,text=offset,command=lambda n=(0,self.offsetItem): self.bt_click(n))
        B.list = []
        B.val = offset
        B.var = ['x','axis_offset']
        B.treatmentId = 4 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.offsetItem.append(B)

        ##########################################
        ###### Y axis
        ylblframe = TTK.LabelFrame(self.frame,text="Y Axis")
        ylblframe.grid(row=2,column=0,sticky='NSEW')
        ylblframe.grid_columnconfigure(0,weight=1)
        ylblframe.grid_columnconfigure(1,weight=1)
        ylblframe.grid_columnconfigure(2,weight=1)
        ylblframe.grid_columnconfigure(3,weight=1)
        ylblframe.grid_columnconfigure(4,weight=1)
        ylblframe.grid_columnconfigure(5,weight=1)
        ylblframe.grid_columnconfigure(6,weight=1)
        ylblframe.grid_columnconfigure(7,weight=1)
        ylblframe.grid_columnconfigure(8,weight=1)
        ylblframe.grid_columnconfigure(9,weight=1)
        ylblframe.grid_rowconfigure(0,weight=1)
        ## ## --> Visible
        visibleframe = TTK.LabelFrame(ylblframe,text="Visible")
        visibleframe.grid(row=0,column=0,sticky='NSEW')
        visibleframe.grid_columnconfigure(0,weight=1)
        visibleframe.grid_rowconfigure(0,weight=1)
        #
        var = TK.IntVar()
        var.set(self.subGraph.axis_property[self.ind_axis].y.axis_visible)
        CB = TTK.Checkbutton(visibleframe,variable=var,command=lambda n=1: self.cb_visible(n))#, variable=var)
        CB.val = var
        CB.var = ['y','axis_visible']
        CB.grid(row=0,column=0,sticky='NSEW')
        self.visibleItem.append(CB)
        ## ## --> Use logscale
        logframe = TTK.LabelFrame(ylblframe,text="Use Log Scale")
        logframe.grid(row=0,column=1,sticky='NSEW')
        logframe.grid_columnconfigure(0,weight=1)
        logframe.grid_rowconfigure(0,weight=1)
        #
        var = TK.IntVar()
        var.set(self.subGraph.axis_property[self.ind_axis].y.axis_logscale)
        CB = TTK.Checkbutton(logframe,variable=var,command=lambda n=1: self.cb_logscale(n))#, variable=var)
        CB.val = var
        CB.var = ['y','axis_logscale']
        CB.grid(row=0,column=0,sticky='NSEW')
        self.logscaleItem.append(CB)
        ## ## --> Use autoscale
        autoscaleframe = TTK.LabelFrame(ylblframe,text="Use Auto Scale")
        autoscaleframe.grid(row=0,column=2,sticky='NSEW')
        autoscaleframe.grid_columnconfigure(0,weight=1)
        autoscaleframe.grid_rowconfigure(0,weight=1)
        #
        var = TK.IntVar()
        var.set(self.subGraph.axis_property[self.ind_axis].y.axis_autoscale)
        CB = TTK.Checkbutton(autoscaleframe,variable=var,command=lambda n=1: self.cb_autoscale(n))#, variable=var)
        CB.val = var
        CB.var = ['y','axis_autoscale']
        CB.grid(row=0,column=0,sticky='NSEW')
        self.autoscaleItem.append(CB)
        ## ## --> Ymin
        ymin = self.subGraph.axis[self.ind_axis].get_ylim()[0]
        yminframe = TTK.LabelFrame(ylblframe, text="Min")
        yminframe.grid(row=0,column=3,sticky='NESW')
        #
        yminframe.grid_columnconfigure(0,weight=1)
        yminframe.grid_rowconfigure(0,weight=1)
        #
        B = TTK.Button(yminframe,text=ymin,command=lambda n=(2,self.positionItem): self.bt_click(n))
        B.list = []
        B.val = ymin
        B.var = ['y','axis_min']
        B.treatmentId = 4 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.positionItem.append(B)
        ## ## --> Ymax
        ymax = self.subGraph.axis[self.ind_axis].get_ylim()[1]
        ymaxframe = TTK.LabelFrame(ylblframe, text="Max")
        ymaxframe.grid(row=0,column=4,sticky='NESW')
        #
        ymaxframe.grid_columnconfigure(0,weight=1)
        ymaxframe.grid_rowconfigure(0,weight=1)
        #
        B = TTK.Button(ymaxframe,text=ymax,command=lambda n=(3,self.positionItem): self.bt_click(n))
        B.list = []
        B.val = ymax
        B.var = ['y','axis_max']
        B.treatmentId = 4 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.positionItem.append(B)
        ## ## --> Invert axis
        invertedframe = TTK.LabelFrame(ylblframe,text="Invert Axis")
        invertedframe.grid(row=0,column=5,sticky='NSEW')
        invertedframe.grid_columnconfigure(0,weight=1)
        invertedframe.grid_rowconfigure(0,weight=1)
        #
        var = TK.IntVar()
        var.set(self.subGraph.axis_property[self.ind_axis].y.axis_inverted)
        CB = TTK.Checkbutton(invertedframe,variable=var,command=lambda n=1: self.cb_invertedAxis(n))#, variable=var)
        CB.val = var
        CB.var = ['y','axis_inverted']
        CB.grid(row=0,column=0,sticky='NSEW')
        self.invertItem.append(CB)
        ## ## --> Label name
        labelframe = TTK.LabelFrame(ylblframe, text="Label name")
        labelframe.grid(row=0,column=6,sticky='NESW')
        #
        labelframe.grid_columnconfigure(0,weight=1)
        labelframe.grid_rowconfigure(0,weight=1)
        #
        B = TTK.Button(labelframe,text=self.subGraph.axis_property[self.ind_axis].y.axis_label,command=lambda n=(1,self.labelItem): self.bt_click(n))
        B.list = []
        B.val = self.subGraph.axis_property[self.ind_axis].y.axis_label
        B.var = ['y','axis_label']
        B.treatmentId = 3
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.labelItem.append(B)
        ## ## --> Label size
        labelframe = TTK.LabelFrame(ylblframe, text="Label size")
        labelframe.grid(row=0,column=7,sticky='NESW')
        #
        labelframe.grid_columnconfigure(0,weight=1)
        labelframe.grid_rowconfigure(0,weight=1)
        #
        B = TTK.Button(labelframe,text=self.subGraph.axis_property[self.ind_axis].y.axis_label_fontsize,command=lambda n=(1,self.labelSizeItem): self.bt_click(n))
        B.list = []
        B.val = self.subGraph.axis_property[self.ind_axis].y.axis_label_fontsize
        B.var = ['y','axis_label_fontsize']
        B.treatmentId = 4 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.labelSizeItem.append(B)
        ## ## --> axis Position
        position = self.subGraph.axis_property[self.ind_axis].y.axis_position
        positionframe = TTK.LabelFrame(ylblframe, text="Axis position")
        positionframe.grid(row=0,column=8,sticky='NESW')
        #
        positionframe.grid_columnconfigure(0,weight=1)
        positionframe.grid_rowconfigure(0,weight=1)
        #
        B = TTK.Button(positionframe,text=position,command=lambda n=(1,self.axisPositionItem): self.bt_click(n))
        B.list = ['left','right','both']
        B.val = position
        B.var = ['y','axis_position']
        B.treatmentId = 0 #
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.axisPositionItem.append(B)
        ## ## --> Axis offset
        offset = self.subGraph.axis_property[self.ind_axis].y.axis_offset
        offsetframe = TTK.LabelFrame(ylblframe, text="Axis offset")
        offsetframe.grid(row=0,column=9,sticky='NESW')
        #
        offsetframe.grid_columnconfigure(0,weight=1)
        offsetframe.grid_rowconfigure(0,weight=1)
        #
        B = TTK.Button(offsetframe,text=offset,command=lambda n=(1,self.offsetItem): self.bt_click(n))
        B.list = []
        B.val = offset
        B.var = ['y','axis_offset']
        B.treatmentId = 4 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.offsetItem.append(B)
    # -------------------------------------------------------------- cmd_addAxis
    def cmd_addAxis(self,event=None):
        method_to_add = self.newAxesMethod[0].val
        if method_to_add == 'New':
            self.subGraph.addAxis()
        elif method_to_add == 'Twin X':
            self.subGraph.addAxisTwinX(self.ind_axis)
        elif method_to_add == 'Twin Y':
            self.subGraph.addAxisTwinY(self.ind_axis)
        self.createFrame()
        # for others related window
        if self.parent.editCurveWdw is not None:
            self.parent.editCurveWdw.updateAxisList()
        for w in [self.parent.editGridWdw]:
            if w is not None:
                w.createFrame()
    # ----------------------------------------------------------- closeAllDialog
    def closeAllDialog(self):
        if self.list_dialog:
            self.list_dialog.cmd_close()
            self.list_dialog = None
        if self.input_dialog:
            self.input_dialog.cmd_close()
            self.input_dialog = None
    # --------------------------------------------------------- updateButonWidth
    def updateButonWidth(self,B,val):
        B.val = val
        B.config(text=B.val)
        try:
            self.subGraph.axis_property[self.ind_axis].setValue(B.var[0],B.var[1],B.val)
            # Update Graph
            self.parent.graphWdwL[self.graph].updateGraph(self.zone)
        except IndexError:
            return
    # ----------------------------------------------------------------- bt_click
    def bt_click(self,ind):
        bt_list = ind[1]
        B = bt_list[ind[0]]
        self.closeAllDialog()
        if B.treatmentId==0:
            self.list_dialog = list_dialogWindow()
            self.list_dialog.initialize(self,B)
        elif B.treatmentId==1:
            self.list_dialog = color_dialogWindow(self,B,B.val)
            # try:
            #     color = askcolor(B.val,parent=self)
            #     if color[1] is not None:
            #         B.val = color[1]
            #         B.config(bg=B.val)
            #         B.config(activebackground=B.val)
            #         try:
            #             self.subGraph.legend_property.setValue(B.var,B.val)
            #             # Update Graph
            #             self.parent.graphWdwL[self.graph].updateGraph(self.zone)
            #         except IndexError:
            #             return
            # except ValueError:
            #     return
        elif B.treatmentId==2:
            self.input_dialog = input_dialogWindow()
            self.input_dialog.initialize(self,B)
        elif B.treatmentId == 3:
            self.input_dialog = input_dialogStringWindow()
            self.input_dialog.initialize(self,B)
        elif B.treatmentId==4:
            self.input_dialog = inputFloat_dialogWindow()
            self.input_dialog.initialize(self,B)
        else:
            return
    # -------------------------------------------------------------- updateButon
    def updateButon(self,B,val):
        B.val = val
        B.config(text=B.val)
        if B.var[0]=='ind_axis':
            self.ind_axis = int(B.val)
            self.createFrame()
        # elif B.var[0]=='axis_to_twin':
        #     self.axis_to_twin = int(B.val)
        else:
            try:
                self.subGraph.axis_property[self.ind_axis].setValue(B.var[0],B.var[1],B.val)
                # Update Graph
                self.parent.graphWdwL[self.graph].updateGraph(self.zone)
            except IndexError:
                return
    # ------------------------------------------------------------- reloadWindow
    def reloadWindow(self):
        self.frame.destroy()
        self.createFrame()
    # ---------------------------------------------------------------- cmd_close
    def cmd_close(self):
        self.closeAllDialog()
        self.parent.editAxisWdw = None
        self.destroy()
    # -------------------------------------------------------------- cb_logscale
    def cb_logscale(self,ind):
        CB = self.logscaleItem[ind]
        self.subGraph.axis_property[self.ind_axis].setValue(CB.var[0],CB.var[1],CB.val.get())
        # Update Graph
        self.parent.graphWdwL[self.graph].updateGraph(self.zone)
    # -------------------------------------------------------------- cb_logscale
    def cb_autoscale(self,ind):
        CB = self.autoscaleItem[ind]
        self.subGraph.axis_property[self.ind_axis].setValue(CB.var[0],CB.var[1],CB.val.get())
        # Update Graph
        self.parent.graphWdwL[self.graph].updateGraph(self.zone)
    # ---------------------------------------------------------- cb_invertedAxis
    def cb_invertedAxis(self,ind):
        CB = self.invertItem[ind]
        self.subGraph.axis_property[self.ind_axis].setValue(CB.var[0],CB.var[1],CB.val.get())
        # Update Graph
        self.parent.graphWdwL[self.graph].updateGraph(self.zone)
    # --------------------------------------------------------------- cb_visible
    def cb_visible(self,ind):
        CB = self.visibleItem[ind]
        self.subGraph.axis_property[self.ind_axis].setValue(CB.var[0],CB.var[1],CB.val.get())
        # Update Graph
        self.parent.graphWdwL[self.graph].updateGraph(self.zone)

# ==============================================================================
# ==============================================================================
class editLegendWindow(TK.Toplevel):
    # --------------------------------------------------------------------- init
    def __init__(self):
        TK.Toplevel.__init__(self)
        self.title('Set legend')
        self.protocol("WM_DELETE_WINDOW",self.cmd_close)
        #
        self.input_dialog = None
        self.list_dialog = None
    # --------------------------------------------------------------- initialize
    def initialize(self,parent):
        self.parent = parent
        #
        self.graph = self.parent.activeGraph.val
        self.zone  = self.parent.position.val
        try:
            self.subGraph = self.parent.graphWdwL[self.graph].fig.subGraph[self.zone]
        except TypeError: # Exit if no graph & zone available
            self.cmd_close()
            return
        #
        self.grid_columnconfigure(0,weight=1)
        self.grid_rowconfigure(0,weight=1)

        #
        self.createFrame()
    # -------------------------------------------------------------- createFrame
    def createFrame(self):
        # Close all opened dialog windows
        self.closeAllDialog()
        #
        self.frame = TTK.Frame(self)
        self.frame.grid(row=0,column=0,sticky='NESW')
        #
        self.frame.grid_columnconfigure(0,weight=1)
        self.frame.grid_rowconfigure(0,weight=1)
        self.frame.grid_rowconfigure(1,weight=1)
        self.frame.grid_rowconfigure(2,weight=1)
        self.frame.grid_rowconfigure(3,weight=1)
        self.frame.grid_rowconfigure(4,weight=1)
        ########################################################
        ## Legend Title
        lblframe = TTK.LabelFrame(self.frame,text='Legend Title')
        lblframe.grid(row=0,column=0,sticky="NSEW")
        #
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_columnconfigure(1,weight=1)
        lblframe.grid_columnconfigure(2,weight=1)
        lblframe.grid_columnconfigure(3,weight=1)
        lblframe.grid_columnconfigure(4,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        ## --> Legend title text
        titleframe = TTK.LabelFrame(lblframe, text="Text")
        titleframe.grid(row=0,column=0,sticky='NESW')
        #
        titleframe.grid_columnconfigure(0,weight=1)
        titleframe.grid_rowconfigure(0,weight=1)
        #
        self.legend_title_textItem = []
        B = TTK.Button(titleframe,text=self.subGraph.legend_property.legend_title,command=lambda n=(0,self.legend_title_textItem): self.bt_click(n))
        B.list = []
        B.val = self.subGraph.legend_property.legend_title
        B.var = 'legend_title'
        B.treatmentId = 3
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.legend_title_textItem.append(B)

        ## --> Bold Button
        boldframe = TTK.LabelFrame(lblframe,text="Bold")
        boldframe.grid(row=0,column=1,sticky='NSEW')
        #
        boldframe.grid_columnconfigure(0,weight=1)
        boldframe.grid_rowconfigure(0,weight=1)
        #
        self.boldItem=[]
        B = TK.Button(boldframe,text='B',command=lambda n=(0,self.boldItem): self.bt_vaEtVient(n),font='bold') # Can not be TTK because of font
        B.grid(row=0,column=0,sticky='NSEW')
        B.nextVal = [1,0]
        B.list = ['normal','bold']
        B.val = self.subGraph.legend_property.legend_title_weight
        B.indVal = B.list.index(B.val)
        B.color = ['#E6E6E6','#848484']
        B.config(bg=B.color[B.indVal])
        B.var = 'legend_title_weight'
        self.boldItem.append(B)

        ## --> Italic Button
        italicframe = TTK.LabelFrame(lblframe,text="Italic")
        italicframe.grid(row=0,column=2,sticky='NSEW')
        #
        italicframe.grid_columnconfigure(0,weight=1)
        italicframe.grid_rowconfigure(0,weight=1)
        #
        self.italicItem=[]
        B = TK.Button(italicframe,text='I',command=lambda n=(0,self.italicItem): self.bt_vaEtVient(n),font='italic') # Can not be TTK because of font
        B.grid(row=0,column=0,sticky='NSEW')
        B.nextVal = [1,0]
        B.list = ['normal','italic']
        B.val = self.subGraph.legend_property.legend_title_style
        B.indVal = B.list.index(B.val)
        B.color = ['#E6E6E6','#848484']
        B.config(bg=B.color[B.indVal])
        B.var = 'legend_title_style'
        self.italicItem.append(B)

        ## --> Color Button
        self.colorItem = []
        colorFrame = TTK.LabelFrame(lblframe, text="Color")
        colorFrame.grid(row=0,column=3,sticky='NESW')
        #
        colorFrame.grid_columnconfigure(0,weight=1)
        colorFrame.grid_rowconfigure(0,weight=1)
        #
        B = TK.Button(colorFrame,command=lambda n=(0,self.colorItem): self.bt_click(n))
        B.list = []
        B.val = self.subGraph.legend_property.legend_title_color
        B.var = 'legend_title_color'
        B.config(bg=B.val)
        B.config(activebackground=B.val)
        B.ind = 0
        B.treatmentId = 1
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.colorItem.append(B)

        ## --> Size
        sizeframe = TTK.LabelFrame(lblframe, text="Text")
        sizeframe.grid(row=0,column=4,sticky='NESW')
        #
        sizeframe.grid_columnconfigure(0,weight=1)
        sizeframe.grid_rowconfigure(0,weight=1)
        #
        self.legend_title_sizeItem = []
        B = TTK.Button(sizeframe,text=self.subGraph.legend_property.legend_title_size,command=lambda n=(0,self.legend_title_sizeItem): self.bt_click(n))
        B.list = []
        B.val = self.subGraph.legend_property.legend_title_size
        B.var = 'legend_title_size'
        B.treatmentId = 4 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.legend_title_sizeItem.append(B)





        ########################################################
        ## Legend label(s)
        lblframe = TTK.LabelFrame(self.frame,text='Legend Label(s)')
        lblframe.grid(row=1,column=0,sticky="NSEW")
        #
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_columnconfigure(1,weight=1)
        lblframe.grid_columnconfigure(2,weight=1)
        lblframe.grid_columnconfigure(3,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        #
        ## --> Bold Button
        boldframe = TTK.LabelFrame(lblframe,text="Bold")
        boldframe.grid(row=0,column=0,sticky='NSEW')
        #
        boldframe.grid_columnconfigure(0,weight=1)
        boldframe.grid_rowconfigure(0,weight=1)
        #
        B = TK.Button(boldframe,text='B',command=lambda n=(1,self.boldItem): self.bt_vaEtVient(n),font='bold') # Can not be TTK because of font
        B.grid(row=0,column=0,sticky='NSEW')
        B.nextVal = [1,0]
        B.list = ['normal','bold']
        B.val = self.subGraph.legend_property.legend_label_weight
        B.indVal = B.list.index(B.val)
        B.color = ['#E6E6E6','#848484']
        B.config(bg=B.color[B.indVal])
        B.var = 'legend_label_weight'
        self.boldItem.append(B)

        ## --> Italic Button
        italicframe = TTK.LabelFrame(lblframe,text="Italic")
        italicframe.grid(row=0,column=1,sticky='NSEW')
        #
        italicframe.grid_columnconfigure(0,weight=1)
        italicframe.grid_rowconfigure(0,weight=1)
        #
        B = TK.Button(italicframe,text='I',command=lambda n=(1,self.italicItem): self.bt_vaEtVient(n),font='italic') # Can not be TTK because of font
        B.grid(row=0,column=0,sticky='NSEW')
        B.nextVal = [1,0]
        B.list = ['normal','italic']
        B.val = self.subGraph.legend_property.legend_label_style
        B.indVal = B.list.index(B.val)
        B.color = ['#E6E6E6','#848484']
        B.config(bg=B.color[B.indVal])
        B.var = 'legend_label_style'
        self.italicItem.append(B)

        ## --> Color Button
        colorFrame = TTK.LabelFrame(lblframe, text="Color")
        colorFrame.grid(row=0,column=2,sticky='NESW')
        #
        colorFrame.grid_columnconfigure(0,weight=1)
        colorFrame.grid_rowconfigure(0,weight=1)
        #
        B = TK.Button(colorFrame,command=lambda n=(1,self.colorItem): self.bt_click(n))
        B.list = []
        B.val = self.subGraph.legend_property.legend_label_color
        B.var = 'legend_label_color'
        B.config(bg=B.val)
        B.config(activebackground=B.val)
        B.ind = 1
        B.treatmentId = 1
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.colorItem.append(B)

        ## --> Size
        sizeframe = TTK.LabelFrame(lblframe, text="Text")
        sizeframe.grid(row=0,column=3,sticky='NESW')
        #
        sizeframe.grid_columnconfigure(0,weight=1)
        sizeframe.grid_rowconfigure(0,weight=1)
        #
        B = TTK.Button(sizeframe,text=self.subGraph.legend_property.legend_label_size,command=lambda n=(1,self.legend_title_sizeItem): self.bt_click(n))
        B.list = []
        B.val = self.subGraph.legend_property.legend_label_size
        B.var = 'legend_label_size'
        B.treatmentId = 4 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.legend_title_sizeItem.append(B)


        ########################################################
        ## Legend Border
        lblframe = TTK.LabelFrame(self.frame,text='Legend Border')
        lblframe.grid(row=2,column=0,sticky="NSEW")
        #
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_columnconfigure(1,weight=1)
        lblframe.grid_columnconfigure(2,weight=1)
        lblframe.grid_columnconfigure(3,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        #
        ## --> Width
        sizeframe = TTK.LabelFrame(lblframe, text="Width")
        sizeframe.grid(row=0,column=0,sticky='NESW')
        #
        sizeframe.grid_columnconfigure(0,weight=1)
        sizeframe.grid_rowconfigure(0,weight=1)
        #
        self.legend_border_widthItem=[]
        B = TTK.Button(sizeframe,text=self.subGraph.legend_property.legend_border_width,command=lambda n=(0,self.legend_border_widthItem): self.bt_click(n))
        B.list = []
        B.val = self.subGraph.legend_property.legend_border_width
        B.var = 'legend_border_width'
        B.treatmentId = 4 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.legend_border_widthItem.append(B)

        ## --> Background Color Button
        backgroundFrame = TTK.LabelFrame(lblframe, text="Background")
        backgroundFrame.grid(row=0,column=1,sticky='NESW')
        #
        backgroundFrame.grid_columnconfigure(0,weight=1)
        backgroundFrame.grid_columnconfigure(1,weight=1)
        backgroundFrame.grid_rowconfigure(0,weight=1)
        #
        colorFrame = TTK.LabelFrame(backgroundFrame,text="Color")
        colorFrame.grid(row=0,column=0,sticky='NSEW')
        #
        colorFrame.grid_columnconfigure(0,weight=1)
        colorFrame.grid_columnconfigure(1,weight=1)
        colorFrame.grid_rowconfigure(0,weight=1)
        #
        B = TK.Button(colorFrame,command=lambda n=(2,self.colorItem): self.bt_click(n))
        B.list = []
        B.val = self.subGraph.legend_property.legend_background_color
        B.var = 'legend_background_color'
        B.config(bg=B.val)
        B.config(activebackground=B.val)
        B.ind = 2
        B.treatmentId = 1
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.colorItem.append(B)
        #
        displayFrame = TTK.LabelFrame(backgroundFrame, text="Active")
        displayFrame.grid(row=0,column=1,sticky='NESW')
        #
        self.checkboxItem = []
        displayFrame.grid_columnconfigure(0,weight=1)
        displayFrame.grid_rowconfigure(0,weight=1)
        var = TK.IntVar()
        var.set(self.subGraph.legend_property.legend_background_color_active)
        CB = TTK.Checkbutton(displayFrame,variable=var,command=lambda n=0: self.cb_visibility(n))#, variable=var)
        CB.val = var
        CB.var = 'legend_background_color_active'
        CB.grid(row=0,column=0,sticky='NSEW')
        self.checkboxItem.append(CB)

        ## --> Color Border
        colorFrame = TTK.LabelFrame(lblframe, text="Color Border")
        colorFrame.grid(row=0,column=2,sticky='NESW')
        #
        colorFrame.grid_columnconfigure(0,weight=1)
        colorFrame.grid_rowconfigure(0,weight=1)
        #
        B = TK.Button(colorFrame,command=lambda n=(3,self.colorItem): self.bt_click(n))
        B.list = []
        B.val = self.subGraph.legend_property.legend_border_color
        B.var = 'legend_border_color'
        B.config(bg=B.val)
        B.config(activebackground=B.val)
        B.ind = 1
        B.treatmentId = 1
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.colorItem.append(B)

        ## --> Number of columns
        sizeframe = TTK.LabelFrame(lblframe, text="# Columns")
        sizeframe.grid(row=0,column=3,sticky='NESW')
        #
        sizeframe.grid_columnconfigure(0,weight=1)
        sizeframe.grid_rowconfigure(0,weight=1)
        #
        self.legend_ncolItem=[]
        B = TTK.Button(sizeframe,text=self.subGraph.legend_property.legend_ncol,command=lambda n=(0,self.legend_ncolItem): self.bt_click(n))
        B.list = []
        B.val = self.subGraph.legend_property.legend_ncol
        B.var = 'legend_ncol'
        B.treatmentId = 2 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.legend_ncolItem.append(B)

        ########################################################
        ## Legend View
        lblframe = TTK.LabelFrame(self.frame,text='Legend Display')
        lblframe.grid(row=3,column=0,sticky="NSEW")
        #
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_columnconfigure(1,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        #
        posFrame = TTK.LabelFrame(lblframe,text='Activate')
        posFrame.grid(row=0,column=0,sticky="NSEW")
        posFrame.grid_columnconfigure(0,weight=1)
        posFrame.grid_rowconfigure(0,weight=1)
        #
#
        self.legend_position = ['best','upper left','upper center','upper right','center left','center','center right','lower left','lower center','lower right']



        cbox = cttk.Combobox(posFrame,values = self.legend_position,state='readonly')
        cbox.val = self.subGraph.legend_property.legend_position
        cbox.set(cbox.val)
        cbox.bind("<<ComboboxSelected>>",self.cmd_positionChange)
        cbox.grid(row=0,column=0,sticky="NSEW")
        cbox.var = 'legend_position'



        #
        activeFrame = TTK.LabelFrame(lblframe,text='Activate')
        activeFrame.grid(row=0,column=1,sticky="NSEW")
        activeFrame.grid_columnconfigure(0,weight=1)
        activeFrame.grid_rowconfigure(0,weight=1)
        #
        var = TK.IntVar()
        var.set(self.subGraph.legend_property.legend_display)
        CB = TTK.Checkbutton(activeFrame,variable=var,command=lambda n=1: self.cb_visibility(n))#, variable=var)
        CB.val = var
        CB.var = 'legend_display'
        CB.grid(row=0,column=0,sticky='NSEW')
        self.checkboxItem.append(CB)
        ########################################################
        ## Close
        B = TTK.Button(self.frame,text='Close',command=self.cmd_close)
        B.grid(row=4,column=0,sticky='NSEW')

    # ----------------------------------------------------------- closeAllDialog
    def closeAllDialog(self):
        if self.list_dialog:
            self.list_dialog.cmd_close()
            self.list_dialog = None
        if self.input_dialog:
            self.input_dialog.cmd_close()
            self.input_dialog = None
    # ------------------------------------------------------- cmd_positionChange
    def cmd_positionChange(self,event):
        widget = event.widget
        val = widget.get()
        # The two following lines seem to be not necessary
        pos = self.legend_position.index(val)
        widget.val = self.legend_position[pos]
        widget.set(widget.val)
        #
        self.subGraph.legend_property.setValue(widget.var,widget.val)
        # Update Graph
        self.parent.graphWdwL[self.graph].updateGraph(self.zone)
    # ------------------------------------------------------------ cb_visibility
    def cb_visibility(self,ind):
        CB = self.checkboxItem[ind]
        self.subGraph.legend_property.setValue(CB.var,CB.val.get())
        # Update Graph
        self.parent.graphWdwL[self.graph].updateGraph(self.zone)
    # ------------------------------------------------------------- bt_vaEtVient
    def bt_vaEtVient(self,var):
        bt_list = var[1]
        B = bt_list[var[0]]
        #
        B.indVal = B.nextVal[B.indVal]
        B.val = B.list[B.indVal]
        B.config(bg=B.color[B.indVal])
        B.config(activebackground=B.color[B.indVal])
        self.subGraph.legend_property.setValue(B.var,B.val)
        # Update Graph
        self.parent.graphWdwL[self.graph].updateGraph(self.zone)
    # --------------------------------------------------------- updateButonWidth
    def updateButonWidth(self,B,val):
        B.val = val
        B.config(text=B.val)
        try:
            self.subGraph.legend_property.setValue(B.var,B.val)
            # Update Graph
            self.parent.graphWdwL[self.graph].updateGraph(self.zone)
        except IndexError:
            return
    # -------------------------------------------------------------- updateButon
    def updateButon(self,B,val):
        B.val = val
        B.config(text=B.val)
        try:
            self.subGraph.legend_property.setValue(B.var,B.val)
            # Update Graph
            self.parent.graphWdwL[self.graph].updateGraph(self.zone)
        except IndexError:
            return
    # ------------------------------------------------------------- reloadWindow
    def reloadWindow(self):
        self.frame.destroy()
        self.createFrame()
    # ---------------------------------------------------------------- cmd_close
    def cmd_close(self):
        self.closeAllDialog()
        self.parent.editLegendWdw = None
        self.destroy()
    # ----------------------------------------------------------------- bt_click
    def updateColor(self,color,B,extra_data):
        if color is not None:
            B.val = color
            B.config(bg=B.val)
            B.config(activebackground=B.val)
            try:
                self.subGraph.legend_property.setValue(B.var,B.val)
                # Update Graph
                self.parent.graphWdwL[self.graph].updateGraph(self.zone)
            except IndexError:
                return
    # ----------------------------------------------------------------- bt_click
    def bt_click(self,ind):
        bt_list = ind[1]
        B = bt_list[ind[0]]
        self.closeAllDialog()
        if B.treatmentId==0:
            self.list_dialog = list_dialogWindow()
            self.list_dialog.initialize(self,B)
        elif B.treatmentId==1:
            self.list_dialog = color_dialogWindow(self,B,B.val)
            # try:
            #     color = askcolor(B.val,parent=self)
            #     if color[1] is not None:
            #         B.val = color[1]
            #         B.config(bg=B.val)
            #         B.config(activebackground=B.val)
            #         try:
            #             self.subGraph.legend_property.setValue(B.var,B.val)
            #             # Update Graph
            #             self.parent.graphWdwL[self.graph].updateGraph(self.zone)
            #         except IndexError:
            #             return
            # except ValueError:
            #     return
        elif B.treatmentId==2:
            self.input_dialog = input_dialogWindow()
            self.input_dialog.initialize(self,B)
        elif B.treatmentId == 3:
            self.input_dialog = input_dialogStringWindow()
            self.input_dialog.initialize(self,B)
        elif B.treatmentId==4:
            self.input_dialog = inputFloat_dialogWindow()
            self.input_dialog.initialize(self,B)
        else:
            return
            
# ==============================================================================
# ==============================================================================
class editGraphWindow(TK.Toplevel):
    # --------------------------------------------------------------------- init
    def __init__(self):
        TK.Toplevel.__init__(self)
        self.title('Graph settings')
        self.protocol("WM_DELETE_WINDOW",self.cmd_close)
        #
        self.input_dialog = None
        self.list_dialog = None
    # --------------------------------------------------------------- initialize
    def initialize(self,parent):
        self.parent = parent
        #
        self.graph = self.parent.activeGraph.val
        self.zone  = self.parent.position.val
        try:
            self.subGraph = self.parent.graphWdwL[self.graph].fig.subGraph[self.zone]
        except TypeError: # Exit if no graph & zone available
            self.cmd_close()
            return
        #
        self.grid_columnconfigure(0,weight=1)
        self.grid_rowconfigure(0,weight=1)
        #
        self.graphInstance = self.parent.graphWdwL[self.graph]
        self.fig = self.graphInstance.fig.instance
        #
        self.createFrame()
    # -------------------------------------------------------------- createFrame
    def createFrame(self):
        # Close all opened dialog windows
        self.closeAllDialog()
        #
        self.checkboxItem = []
        self.useSubPlotParams = True
        #
        self.frame = TTK.Frame(self)
        self.frame.grid(row=0,column=0,sticky='NESW')
        #
        self.frame.grid_columnconfigure(0,weight=6)
        self.frame.grid_columnconfigure(1,weight=1)
        self.frame.grid_rowconfigure(0,weight=1)
        self.frame.grid_rowconfigure(1,weight=1)
        self.frame.grid_rowconfigure(2,weight=1)
        self.frame.grid_rowconfigure(3,weight=1)
        self.frame.grid_rowconfigure(4,weight=0)
        #
        ### Row 0 -> Colors
        #
        self.colorItem = []
        colorLblFrame = TTK.LabelFrame(self.frame,text="Background colors")
        colorLblFrame.grid_columnconfigure(0,weight=2)
        colorLblFrame.grid_columnconfigure(1,weight=2)
        colorLblFrame.grid_rowconfigure(0,weight=1)
        #
        colorLblFrame.grid(row=0,column=0,columnspan=2,sticky='NSEW')
        # ---> Graphic zone
        graph = self.parent.graphWdwL[self.graph]
        graphicFrame = TTK.LabelFrame(colorLblFrame,text="Graphic Zone")
        graphicFrame.grid_columnconfigure(0,weight=1)
        graphicFrame.grid_columnconfigure(1,weight=1)
        graphicFrame.grid_rowconfigure(0,weight=1)
        #
        graphicFrame.grid(row=0,column=0,sticky='NSEW')
        # ---> ---> Select color
        colorFrame = TTK.LabelFrame(graphicFrame, text="Color")
        colorFrame.grid(row=0,column=0,sticky='NESW')
        #
        colorFrame.grid_columnconfigure(0,weight=1)
        colorFrame.grid_rowconfigure(0,weight=1)
        #
        ind = 0
        B = TK.Button(colorFrame,command=lambda n=(ind,self.colorItem): self.bt_click(n))
        B.list = []
        # B.val = "#FFFFFF" # white
        B.val = graph.subgraph_background_color
        B.var = 'grid_color'
        B.config(bg=B.val)
        B.config(activebackground=B.val)
        B.ind = ind
        B.treatmentId = 1
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        B.which = ('major','x')
        self.colorItem.append(B)
        # ---> ---> alpha Slider
        alphaFrame = TTK.LabelFrame(graphicFrame, text="Alpha")
        alphaFrame.grid(row=0,column=1,sticky='NESW')
        #
        alphaFrame.grid_columnconfigure(0,weight=1)
        alphaFrame.grid_rowconfigure(0,weight=1)
        alphaFrame.grid_rowconfigure(1,weight=1)
        # Create a string var
        self.graphicAlpha = TK.StringVar()
        # its initial Value is the current value of the graph attributes
        val = graph.subgraph_background_alpha
        self.graphicAlpha.set('%3d %%'%(int(val)*100))
        # Link this string var to the function alphaChange
        self.graphicAlpha.trace("w",lambda a,b,c,loc="graphic",var=self.graphicAlpha: self.alphaChange(loc,var))
        # Create the slider linked to the string var
        alphaSlider = TTK.Scale(alphaFrame,from_=0,to=100,orient=TK.HORIZONTAL,showvalue=0,command=lambda s:self.graphicAlpha.set('%3d %%'%float(s)))
        alphaSlider.set(val*100)
        alphaSlider.grid(row=1,column=0,sticky="nsew")
        graphicLabel = TTK.Label(alphaFrame, textvariable=self.graphicAlpha)
        graphicLabel.grid(row=0,column=0,sticky='nsew')   # labels that will update
        # ---> Border zone
        borderFrame = TTK.LabelFrame(colorLblFrame,text="Border Zone")
        borderFrame.grid_columnconfigure(0,weight=1)
        borderFrame.grid_columnconfigure(1,weight=1)
        borderFrame.grid_rowconfigure(0,weight=1)
        #
        borderFrame.grid(row=0,column=1,sticky='NSEW')
        # ---> ---> Select color
        colorFrame = TTK.LabelFrame(borderFrame, text="Color")
        colorFrame.grid(row=0,column=0,sticky='NESW')
        #
        colorFrame.grid_columnconfigure(0,weight=1)
        colorFrame.grid_rowconfigure(0,weight=1)
        #
        ind = 1
        B = TK.Button(colorFrame,command=lambda n=(ind,self.colorItem): self.bt_click(n))
        B.list = []
        # B.val = "#FFFFFF" # White
        B.val = graph.image_background_color
        B.var = 'grid_color'
        B.config(bg=B.val)
        B.config(activebackground=B.val)
        B.ind = ind
        B.treatmentId = 1
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        B.which = ('major','x')
        self.colorItem.append(B)
        # ---> ---> alpha Slider
        alphaFrame = TTK.LabelFrame(borderFrame, text="Alpha")
        alphaFrame.grid(row=0,column=1,sticky='NESW')
        #
        alphaFrame.grid_columnconfigure(0,weight=1)
        alphaFrame.grid_rowconfigure(0,weight=1)
        alphaFrame.grid_rowconfigure(1,weight=1)
        # Create a string var
        self.borderAlpha = TK.StringVar()
        # Set its initial value to the current value
        val = graph.image_background_alpha
        self.borderAlpha.set('%3d %%'%(val*100))
        # Link the modification of its vlue to the function alphaChange
        self.borderAlpha.trace("w",lambda a,b,c,loc="border",var=self.borderAlpha: self.alphaChange(loc,var))
        # Create the slider linked to the string var
        alphaSlider = TTK.Scale(alphaFrame,from_=0,to=100,orient=TK.HORIZONTAL,showvalue=0,command=lambda s:self.borderAlpha.set('%3d %%'%float(s)))
        alphaSlider.set(int(val*100))
        alphaSlider.grid(row=1,column=0,sticky="nsew")
        borderLabel = TTK.Label(alphaFrame, textvariable=self.borderAlpha)   # labels that will update
        borderLabel.grid(row=0,column=0,sticky='nsew')



        #
        ### Row 1 -> SubPlotParams
        #
        SPPlblframe = TTK.LabelFrame(self.frame,text="SubPlotParams")
        SPPlblframe.grid_columnconfigure(0,weight=1)
        SPPlblframe.grid_columnconfigure(1,weight=1)
        SPPlblframe.grid_columnconfigure(2,weight=1)
        SPPlblframe.grid_columnconfigure(3,weight=1)
        SPPlblframe.grid_columnconfigure(4,weight=1)
        SPPlblframe.grid_columnconfigure(5,weight=1)
        SPPlblframe.grid_rowconfigure(0,weight=1)
        #
        SPPlblframe.grid(row=1,column=0,sticky='NSEW')
        self.subPlotParamsBtnItem = {}
        # --->Left
        lblframe = TTK.LabelFrame(SPPlblframe,text="Left")
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        lblframe.grid(row=0,column=0,sticky='NSEW')
        #
        B = TTK.Button(lblframe,text=self.graphInstance.subPlotParams.left,command=lambda n=('left',self.subPlotParamsBtnItem): self.bt_click(n))
        B.list = []
        B.val = self.graphInstance.subPlotParams.left
        B.var = 'left'
        B.treatmentId = 4 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.subPlotParamsBtnItem['left']=B
        # --->Right
        lblframe = TTK.LabelFrame(SPPlblframe,text="Right")
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        lblframe.grid(row=0,column=1,sticky='NSEW')
        #
        B = TTK.Button(lblframe,text=self.graphInstance.subPlotParams.right,command=lambda n=('right',self.subPlotParamsBtnItem): self.bt_click(n))
        B.list = []
        B.val = self.graphInstance.subPlotParams.right
        B.var = 'right'
        B.treatmentId = 4 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.subPlotParamsBtnItem['right']=B
        # --->Top
        lblframe = TTK.LabelFrame(SPPlblframe,text="Top")
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        lblframe.grid(row=0,column=2,sticky='NSEW')
        #
        B = TTK.Button(lblframe,text=self.graphInstance.subPlotParams.top,command=lambda n=('top',self.subPlotParamsBtnItem): self.bt_click(n))
        B.list = []
        B.val = self.graphInstance.subPlotParams.top
        B.var = 'top'
        B.treatmentId = 4 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.subPlotParamsBtnItem['top']=B
        # --->Bottom
        lblframe = TTK.LabelFrame(SPPlblframe,text="Bottom")
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        lblframe.grid(row=0,column=3,sticky='NSEW')
        #
        B = TTK.Button(lblframe,text=self.graphInstance.subPlotParams.bottom,command=lambda n=('bottom',self.subPlotParamsBtnItem): self.bt_click(n))
        B.list = []
        B.val = self.graphInstance.subPlotParams.bottom
        B.var = 'bottom'
        B.treatmentId = 4 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.subPlotParamsBtnItem['bottom']=B
        # --->H.Space
        lblframe = TTK.LabelFrame(SPPlblframe,text="H.Space")
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        lblframe.grid(row=0,column=4,sticky='NSEW')
        #
        B = TTK.Button(lblframe,text=self.graphInstance.subPlotParams.hspace,command=lambda n=('hspace',self.subPlotParamsBtnItem): self.bt_click(n))
        B.list = []
        B.val = self.graphInstance.subPlotParams.hspace
        B.var = 'hspace'
        B.treatmentId = 4 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.subPlotParamsBtnItem['hspace']=B
        # --->W.Space
        lblframe = TTK.LabelFrame(SPPlblframe,text="W.Space")
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        lblframe.grid(row=0,column=5,sticky='NSEW')
        #
        B = TTK.Button(lblframe,text=self.graphInstance.subPlotParams.wspace,command=lambda n=('wspace',self.subPlotParamsBtnItem): self.bt_click(n))
        B.list = []
        B.val = self.graphInstance.subPlotParams.wspace
        B.var = 'wspace'
        B.treatmentId = 4 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.subPlotParamsBtnItem['wspace']=B


        # -> Active/inactive checkbox
        CBframe = TTK.LabelFrame(self.frame,text="Active")
        CBframe.grid_rowconfigure(0,weight=1)
        CBframe.grid_columnconfigure(0,weight=1)
        CBframe.grid(row=1,column=1,sticky='NSEW')
        #
        #
        var = TK.IntVar()
        var.set(1)
        CB = TTK.Checkbutton(CBframe,variable=var,command=lambda n=0: self.cb_active(n))#, variable=var)
        CB.val = var
        CB.var = 'active'
        CB.grid(row=0,column=0,sticky='NSEW')
        self.checkboxItem.append(CB)

        #
        ## Row 2 -> Tight Layout
        #
        TLlblframe = TTK.LabelFrame(self.frame,text="Tight Layout")
        TLlblframe.grid_columnconfigure(0,weight=1)
        TLlblframe.grid_columnconfigure(1,weight=1)
        TLlblframe.grid_columnconfigure(2,weight=1)
        TLlblframe.grid_columnconfigure(3,weight=1)
        TLlblframe.grid_columnconfigure(4,weight=1)
        TLlblframe.grid_columnconfigure(5,weight=1)
        TLlblframe.grid_rowconfigure(0,weight=1)
        #
        TLlblframe.grid(row=2,column=0,sticky='NSEW')
        #
        self.tightLayoutBtnItem = {}
        # --->Pad
        lblframe = TTK.LabelFrame(TLlblframe,text="Pad")
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        lblframe.grid(row=0,column=0,columnspan=2,sticky='NSEW')
        #
        B = TTK.Button(lblframe,text=self.graphInstance.tightLayout.pad,command=lambda n=('pad',self.tightLayoutBtnItem): self.bt_click(n))
        B.list = []
        B.val = self.graphInstance.tightLayout.pad
        B.var = 'pad'
        B.treatmentId = 4 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.tightLayoutBtnItem['pad']=B
        # --->H.Pad
        lblframe = TTK.LabelFrame(TLlblframe,text="H.Pad")
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        lblframe.grid(row=0,column=2,columnspan=2,sticky='NSEW')
        #
        B = TTK.Button(lblframe,text=self.graphInstance.tightLayout.hpad,command=lambda n=('hpad',self.tightLayoutBtnItem): self.bt_click(n))
        B.list = []
        B.val = self.graphInstance.tightLayout.hpad
        B.var = 'hpad'
        B.treatmentId = 4 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.tightLayoutBtnItem['hpad']=B
        # --->W.Pad
        lblframe = TTK.LabelFrame(TLlblframe,text="W.Pad")
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        lblframe.grid(row=0,column=4,columnspan=2,sticky='NSEW')
        #
        B = TTK.Button(lblframe,text=self.graphInstance.tightLayout.wpad,command=lambda n=('wpad',self.tightLayoutBtnItem): self.bt_click(n))
        B.list = []
        B.val = self.graphInstance.tightLayout.wpad
        B.var = 'wpad'
        B.treatmentId = 4 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        self.tightLayoutBtnItem['wpad']=B


        # -> Active/inactive checkbox
        CBframe = TTK.LabelFrame(self.frame,text="Active")
        CBframe.grid_rowconfigure(0,weight=1)
        CBframe.grid_columnconfigure(0,weight=1)
        CBframe.grid(row=2,column=1,sticky='NSEW')
        #
        #
        var = TK.IntVar()
        var.set(0)
        CB = TTK.Checkbutton(CBframe,variable=var,command=lambda n=1: self.cb_active(n))#, variable=var)
        CB.val = var
        CB.var = 'active'
        CB.grid(row=0,column=0,sticky='NSEW')
        self.checkboxItem.append(CB)
        #
        ### Row 3 -> Size and Resolution
        #
        SaRlblframe = TTK.LabelFrame(self.frame,text="Size & resolution")
        SaRlblframe.grid_columnconfigure(0,weight=1)
        SaRlblframe.grid_columnconfigure(1,weight=1)
        SaRlblframe.grid_columnconfigure(2,weight=1)
        SaRlblframe.grid_columnconfigure(3,weight=1)
        SaRlblframe.grid_columnconfigure(4,weight=1)
        SaRlblframe.grid_columnconfigure(5,weight=1)
        SaRlblframe.grid_columnconfigure(6,weight=1)
        SaRlblframe.grid_rowconfigure(0,weight=1)
        #
        SaRlblframe.grid(row=3,column=0,columnspan=2,sticky='NSEW')

        self.sizeItem={}
        # -> Fig.Size
        lblframe = TTK.LabelFrame(SaRlblframe,text="Fig.Size (inches)")
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        lblframe.grid(row=0,column=2,columnspan=1,sticky='NSEW')
        figsize = self.fig.get_size_inches()
        text = (int(figsize[0]),int(figsize[1]))
        B = TTK.Button(lblframe,text=str(text),command=lambda n=('figsize',self.sizeItem): self.bt_click_pattern(n))
        B.val = str(text)
        B.var = 'figsize'
        B.pattern = '^\(\s*[0123456789]*\s*,\s*[0123456789]*\s*\)$'
        B.grid(row=0,column=0,sticky="nsew")
        self.sizeItem['figsize']=B
        # -> DPI
        lblframe = TTK.LabelFrame(SaRlblframe,text="DPI (Dot Per Inch)")
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        lblframe.grid(row=0,column=4,columnspan=1,sticky='NSEW')
        B = TTK.Button(lblframe,text=self.fig.get_dpi(),command=lambda n=('dpi',self.sizeItem): self.bt_click_pattern(n))
        B.val = self.graphInstance.fig.instance.get_dpi()
        B.var = 'dpi'
        B.pattern = '^[0123456789]*$'
        B.grid(row=0,column=0,sticky="nsew")
        self.sizeItem['dpi']=B

        ### TODO : recuperer les valeurs de dpi et de figsize comme dans le desktopframetk

        #
        ### Row 4 -> Buttons
        #
        globalframe = TTK.Frame(self.frame)
        globalframe.grid_rowconfigure(0,weight=1)
        globalframe.grid_columnconfigure(0,weight=1)
        globalframe.grid_columnconfigure(1,weight=1)
        globalframe.grid_columnconfigure(2,weight=1)
        globalframe.grid_columnconfigure(3,weight=1)
        globalframe.grid_columnconfigure(4,weight=1)
        globalframe.grid_columnconfigure(5,weight=1)
        globalframe.grid_columnconfigure(6,weight=1)
        globalframe.grid(row=4,column=0,columnspan=2,sticky='EW')
        # -> Set button
        B = TTK.Button(globalframe,text="Configure",command=self.cmd_configure)
        B.grid(row=0,column=2,columnspan=1,sticky="nsew")
        # -> Close button
        B = TTK.Button(globalframe,text="Close",command=self.cmd_close)
        B.grid(row=0,column=4,columnspan=1,sticky="nsew")
    # -------------------------------------------------------------- alphaChange
    def alphaChange(self,loc,var,*args):
        graph = self.parent.graphWdwL[self.graph]
        val = var.get()
        alpha = float(val.split()[0])/100.
        if loc =="border":
            graph.image_background_alpha = alpha
            fig = self.parent.graphWdwL[self.graph].getFig()
            fig.patch.set_alpha(alpha)
        elif loc == "graphic":
            graph.subgraph_background_alpha = alpha
        self.parent.graphWdwL[self.graph].updateGraph(self.zone)
    # ------------------------------------------------------------- reloadWindow
    def reloadWindow(self):
        self.frame.destroy()
        self.createFrame()
    # --------------------------------------------------------- updateButonWidth
    def updateButonWidth(self,B,val):
        B.val = val
        B.config(text=B.val)

    # ---------------------------------------------------------------- cb_active
    def cb_active(self,ind):
        CB = self.checkboxItem[ind]
        CB_other = self.checkboxItem[abs(ind-1)] # Trick to access the proper checkbox item, not clean but still working
        val = CB.val.get()
        CB_other.val.set(abs(val-1)) ## Trick to set the proper value, not clean but still working
        # Update self.useSubPlotParams
        if self.checkboxItem[0].val.get():
            self.useSubPlotParams = True
        else:
            self.useSubPlotParams = False
#        self.subGraph.grid_property[self.ind_axis].setValue(CB.which[0],CB.which[1],'display',CB.val.get())
#        # Update Graph
#        self.parent.graphWdwL[self.graph].updateGraph(self.zone)
    # ------------------------------------------------------------ cmd_configure
    def cmd_configure(self):
        # print('configure was hit')

        # 1/- Set size and resolution

        # Can not succeed to adapt dinamicaly the size of the canvas !!!
        # Nevertheless, the user, can change the size of the toplevel window (graphtk) with the mouse,
        # it will change the size of the figure in the mean time.
        # TODO : find a real a solution to this issue !!!
#        tkwidget = self.graphInstance.canvas.get_tk_widget()
#        figsizeStr = self.sizeItem['figsize'].val[1:-1]
#        figsize = (int(figsizeStr.split(',')[0]),int(figsizeStr.split(',')[1]))
#        self.fig.set_size_inches(figsize,forward=True)
#        tkwidget = self.graphInstance.canvas.get_tk_widget()
        #
        dpi = int(self.sizeItem['dpi'].val)
        self.fig.set_dpi(dpi)
        # Update
        self.graphInstance.canvas.draw()

        # 2/- Nice view
        if self.useSubPlotParams :
            params = {'left':None,'right':None,'top':None,'bottom':None,'hspace':None,'wspace':None,'isActive':True}
            for k in params:
                if k != 'isActive':
                    if self.subPlotParamsBtnItem[k].val is not None:
                        params[k] = float(self.subPlotParamsBtnItem[k].val)
                    else:
                        params[k] = None
            isParamsOk = self.graphInstance.subPlotParams.checkParams(params)

            #
            if isParamsOk:
                self.graphInstance.updateSubPlotParams(params)
            else:
                #
#                self.lift()
#                self.focus()
                return
        else:
            params = {'pad':None,'hpad':None,'wpad':None,'isActive':True}
            for k in params:
                if k != 'isActive':
                    if self.tightLayoutBtnItem[k].val is not None:
                        params[k] = float(self.tightLayoutBtnItem[k].val)
                    else:
                        params[k] = None
            self.graphInstance.updateTightLayout(params)
        self.cmd_close()

    # ----------------------------------------------------------- closeAllDialog
    def closeAllDialog(self):
        if self.list_dialog:
            self.list_dialog.cmd_close()
            self.list_dialog = None
        if self.input_dialog:
            self.input_dialog.cmd_close()
            self.input_dialog = None
    # --------------------------------------------------------- bt_click_pattern
    def bt_click_pattern(self,ind):
        bt_dict = ind[1]
        B = bt_dict[ind[0]]
        self.closeAllDialog()
        self.input_dialog = inputPattern_dialogWindow(B.pattern)
        self.input_dialog.initialize(self,B)

    # ----------------------------------------------------------------- updateColor
    def updateColor(self,color,B,extra_data):
        if color is not None:
            B.val = color
            B.config(bg=B.val)
            B.config(activebackground=B.val)
            try:
                graph = self.parent.graphWdwL[self.graph]
                if B.ind == 1 :
                    graph.image_background_color = B.val
                    fig = self.parent.graphWdwL[self.graph].getFig()
                    fig.patch.set_facecolor(B.val)
                elif B.ind == 0 :
                    graph.subgraph_background_color = B.val
                    # fig = self.parent.graphWdwL[self.graph].getFig()
                    # for axe in fig.get_axes():
                    #     axe.patch.set_facecolor(B.val)
                        # axe.patch.set_alpha(1.0)
                    # for axes in self.subGraph.axis :
                    #     # axes.patch.set_facecolor(B.val)
                    #     axes.set_axis_bgcolor(B.val)
                # Update Graph
                self.parent.graphWdwL[self.graph].updateGraph(self.zone)
            except IndexError:
                return
    # ----------------------------------------------------------------- bt_click
    def bt_click(self,ind):
        bt_dict = ind[1]
        B = bt_dict[ind[0]]
        self.closeAllDialog()

        if B.treatmentId==1: # Color selector
            self.list_dialog = color_dialogWindow(self,B,B.val)
            # try:
            #     color = askcolor(B.val,parent=self)
            #     if color[1] is not None:
            #         B.val = color[1]
            #         B.config(bg=B.val)
            #         B.config(activebackground=B.val)
            #         try:
            #             graph = self.parent.graphWdwL[self.graph]
            #             if B.ind == 1 :
            #                 graph.image_background_color = B.val
            #                 fig = self.parent.graphWdwL[self.graph].getFig()
            #                 fig.patch.set_facecolor(B.val)
            #             elif B.ind == 0 :
            #                 graph.subgraph_background_color = B.val
            #                 # fig = self.parent.graphWdwL[self.graph].getFig()
            #                 # for axe in fig.get_axes():
            #                 #     axe.patch.set_facecolor(B.val)
            #                     # axe.patch.set_alpha(1.0)
            #                 # for axes in self.subGraph.axis :
            #                 #     # axes.patch.set_facecolor(B.val)
            #                 #     axes.set_axis_bgcolor(B.val)
            #             # Update Graph
            #             self.parent.graphWdwL[self.graph].updateGraph(self.zone)
            #         except IndexError:
            #             return
            # except ValueError:
            #     return
        else:
            self.input_dialog = inputFloat_dialogWindow()
            self.input_dialog.initialize(self,B)

    # ------------------------------------------------------------- reloadWindow
    def reloadWindow(self):
        self.frame.destroy()
        self.createFrame()
    # ---------------------------------------------------------------- cmd_close
    def cmd_close(self):
        self.closeAllDialog()
        self.parent.editGraphWdw = None
        self.destroy()
    # --------------------------------------------------------- updateButon
    def updateButon(self,B,val):
        B.val = val
        B.config(text=B.val)
        if B.var == 'axis':
            self.ind_axis = B.val
            self.createFrame()

        else:
            try:
                self.subGraph.grid_property[self.ind_axis].setValue(B.which[0],B.which[1],'grid_style',B.val)
                # Update Graph
                self.parent.graphWdwL[self.graph].updateGraph(self.zone)
            except IndexError:
                return

# ==============================================================================
# ==============================================================================
class editGridWindow(TK.Toplevel):
    # --------------------------------------------------------------------- init
    def __init__(self):
        TK.Toplevel.__init__(self)
        self.title('Set grid')
        self.protocol("WM_DELETE_WINDOW",self.cmd_close)
        #
        self.input_dialog = None
        self.list_dialog = None
    # --------------------------------------------------------------- initialize
    def initialize(self,parent):
        self.parent = parent
        #
        self.graph = self.parent.activeGraph.val
        self.zone  = self.parent.position.val
        try:
            self.subGraph = self.parent.graphWdwL[self.graph].fig.subGraph[self.zone]
        except TypeError: # Exit if no graph & zone available
            self.cmd_close()
            return
        #
        self.grid_columnconfigure(0,weight=1)
        self.grid_rowconfigure(0,weight=1)

        #
        self.ind_axis = 0
        self.createFrame()
    # -------------------------------------------------------------- createFrame
    def createFrame(self):
        # Close all opened dialog windows
        self.closeAllDialog()
        #
        self.frame = TTK.Frame(self)
        self.frame.grid(row=0,column=0,sticky='NESW')
        #
        self.frame.grid_columnconfigure(0,weight=1)
        self.frame.grid_rowconfigure(0,weight=1)
        self.frame.grid_rowconfigure(1,weight=1)
        self.frame.grid_rowconfigure(2,weight=1)
        ########################################################
        ## Select axis
        #
        self.addAxisItem = []
        selectLblFrame = TTK.LabelFrame(self.frame,text='Select Axis')
        selectLblFrame.grid_columnconfigure(0, weight=1)
        selectLblFrame.grid_rowconfigure(0, weight=1)
        selectLblFrame.grid(row=0, column=0, sticky='NSEW')
        B = TTK.Button(selectLblFrame,text=self.ind_axis,command=lambda n=(0,self.addAxisItem): self.bt_click(n))
        self.addAxisItem.append(B)
        B.list = [i for i in range(len(self.subGraph.axis))]
        B.val = self.ind_axis
        B.var = 'axis'
        B.treatmentId = 0 # 4 is float here
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")

        ########################################################
        ## Major grid
        lblframe = TTK.LabelFrame(self.frame,text='Major grid')
        lblframe.grid(row=1,column=0,sticky="NSEW")
        #
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        lblframe.grid_rowconfigure(1,weight=1)
        ### --> X axis
        ind = 0
        xlblframe = TTK.LabelFrame(lblframe,text='X')
        xlblframe.grid(row=0,column=0,sticky="NSEW")
        #
        xlblframe.grid_columnconfigure(0,weight=1)
        xlblframe.grid_columnconfigure(1,weight=1)
        xlblframe.grid_columnconfigure(2,weight=1)
        xlblframe.grid_columnconfigure(3,weight=1)
        xlblframe.grid_columnconfigure(4,weight=1)
        xlblframe.grid_columnconfigure(5,weight=1)
        xlblframe.grid_rowconfigure(0,weight=1)
        ### ### ---> Display
        self.visibilityItem=[]
        displayFrame = TTK.LabelFrame(xlblframe, text="Display")
        displayFrame.grid(row=0,column=0,sticky='NESW')
        #
        displayFrame.grid_columnconfigure(0,weight=1)
        displayFrame.grid_rowconfigure(0,weight=1)
        var = TK.IntVar()
        var.set(self.subGraph.grid_property[self.ind_axis].major.x.display)
        CB = TTK.Checkbutton(displayFrame,variable=var,command=lambda n=ind: self.cb_visibility(n))#, variable=var)
        CB.val = var
        CB.var = 'display'
        CB.which=('major','x')
        CB.grid(row=0,column=0,sticky='NSEW')
        self.visibilityItem.append(CB)
#        ### ### ---> Grid color
        self.colorItem = []
        colorFrame = TTK.LabelFrame(xlblframe, text="Color")
        colorFrame.grid(row=0,column=1,sticky='NESW')
        #
        colorFrame.grid_columnconfigure(0,weight=1)
        colorFrame.grid_rowconfigure(0,weight=1)
        #
        B = TK.Button(colorFrame,command=lambda n=(ind,self.colorItem): self.bt_click(n))
        B.list = []
        B.val = self.subGraph.grid_property[self.ind_axis].major.x.grid_color
        B.var = 'grid_color'
        B.config(bg=B.val)
        B.config(activebackground=B.val)
        B.ind = ind
        B.treatmentId = 1
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        B.which = ('major','x')
        self.colorItem.append(B)
        ### ### ---> Grid width
        self.widthItem = []
        widthFrame = TTK.LabelFrame(xlblframe, text="Width")
        widthFrame.grid(row=0,column=2,sticky='NESW')
        #
        widthFrame.grid_columnconfigure(0,weight=1)
        widthFrame.grid_rowconfigure(0,weight=1)
        #
        B = TTK.Button(widthFrame,text=self.subGraph.grid_property[self.ind_axis].major.x.grid_width,command=lambda n=(ind,self.widthItem): self.bt_click(n))
        B.list = []
        B.val = self.subGraph.grid_property[self.ind_axis].major.x.grid_width
        B.var = 'grid_width'
        B.ind = ind
        B.treatmentId = 2
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        B.which = ('major','x')
        self.widthItem.append(B)
        ### ### ---> Grid style
        self.styleItem = []
        styleFrame = TTK.LabelFrame(xlblframe, text="Style")
        styleFrame.grid(row=0,column=3,sticky='NESW')
        #
        styleFrame.grid_columnconfigure(0,weight=1)
        styleFrame.grid_rowconfigure(0,weight=1)
        #
        B = TTK.Button(styleFrame,text=self.subGraph.grid_property[self.ind_axis].major.x.grid_style,command=lambda n=(ind,self.styleItem): self.bt_click(n))
        B.list = linestylelist
        B.val = self.subGraph.grid_property[self.ind_axis].major.x.grid_style
        B.var = 'grid_style'
        B.ind = ind
        B.treatmentId = 0
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        B.which = ('major','x')
        self.styleItem.append(B)
        ### ### ---> Tick number
        self.tickNumberItem = []
        ticknumberFrame = TTK.LabelFrame(xlblframe, text="# Ticks")
        ticknumberFrame.grid(row=0,column=4,sticky='NESW')
        #
        ticknumberFrame.grid_columnconfigure(0,weight=1)
        ticknumberFrame.grid_rowconfigure(0,weight=1)
        #
        B = TTK.Button(ticknumberFrame,text=self.subGraph.grid_property[self.ind_axis].major.x.grid_tick_number,command=lambda n=(ind,self.tickNumberItem): self.bt_click(n))
        B.pattern = '''^[0-9]*$'''
        B.val = self.subGraph.grid_property[self.ind_axis].major.x.grid_tick_number
        B.var = 'grid_tick_number'
        B.ind = ind
        B.treatmentId = 3
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        B.which = ('major','x')
        self.tickNumberItem.append(B)

        ### ### ---> Tick size
        self.tickSizeItem = []
        ticksizeFrame = TTK.LabelFrame(xlblframe, text="Ticks Size")
        ticksizeFrame.grid(row=0,column=5,sticky='NESW')
        #
        ticksizeFrame.grid_columnconfigure(0,weight=1)
        ticksizeFrame.grid_rowconfigure(0,weight=1)
        #
        B = TTK.Button(ticksizeFrame,text=self.subGraph.grid_property[self.ind_axis].major.x.grid_tick_size,command=lambda n=(ind,self.tickSizeItem): self.bt_click(n))
        B.list = []
        B.val = self.subGraph.grid_property[self.ind_axis].major.x.grid_tick_size
        B.var = 'grid_tick_size'
        B.ind = ind
        B.treatmentId = 2
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        B.which = ('major','x')
        self.tickSizeItem.append(B)

        ### --> Y axis
        ind = 1
        ylblframe = TTK.LabelFrame(lblframe,text='Y')
        ylblframe.grid(row=1,column=0,sticky="NSEW")
        #
        ylblframe.grid_columnconfigure(0,weight=1)
        ylblframe.grid_columnconfigure(1,weight=1)
        ylblframe.grid_columnconfigure(2,weight=1)
        ylblframe.grid_columnconfigure(3,weight=1)
        ylblframe.grid_columnconfigure(4,weight=1)
        ylblframe.grid_columnconfigure(5,weight=1)
        ylblframe.grid_rowconfigure(0,weight=1)
        ### ### ---> Display
        displayFrame = TTK.LabelFrame(ylblframe, text="Display")
        displayFrame.grid(row=0,column=0,sticky='NESW')
        #
        displayFrame.grid_columnconfigure(0,weight=1)
        displayFrame.grid_rowconfigure(0,weight=1)
        var = TK.IntVar()
        var.set(self.subGraph.grid_property[self.ind_axis].major.y.display)
        CB = TTK.Checkbutton(displayFrame,variable=var,command=lambda n=ind: self.cb_visibility(n))#, variable=var)
        CB.val = var
        CB.var = 'display'
        CB.which = ('major','y')
        CB.grid(row=0,column=0,sticky='NSEW')
        self.visibilityItem.append(CB)
        ### ### ---> Line color
        colorFrame = TTK.LabelFrame(ylblframe, text="Color")
        colorFrame.grid(row=0,column=1,sticky='NESW')
        #
        colorFrame.grid_columnconfigure(0,weight=1)
        colorFrame.grid_rowconfigure(0,weight=1)
        #
        B = TK.Button(colorFrame,command=lambda n=(ind,self.colorItem): self.bt_click(n))
        B.list = []
        B.val = self.subGraph.grid_property[self.ind_axis].major.y.grid_color
        B.var = 'grid_color'
        B.config(bg=B.val)
        B.config(activebackground=B.val)
        B.ind = ind
        B.treatmentId = 1
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        B.which = ('major','y')
        self.colorItem.append(B)
        ### ### ---> Grid width
        widthFrame = TTK.LabelFrame(ylblframe, text="Width")
        widthFrame.grid(row=0,column=2,sticky='NESW')
        #
        widthFrame.grid_columnconfigure(0,weight=1)
        widthFrame.grid_rowconfigure(0,weight=1)
        #
        B = TTK.Button(widthFrame,text=self.subGraph.grid_property[self.ind_axis].major.y.grid_width,command=lambda n=(ind,self.widthItem): self.bt_click(n))
        B.list = []
        B.val = self.subGraph.grid_property[self.ind_axis].major.y.grid_width
        B.var = 'grid_width'
        B.ind = ind
        B.treatmentId = 2
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        B.which = ('major','y')
        self.widthItem.append(B)
        ### ### ---> Grid style
        styleFrame = TTK.LabelFrame(ylblframe, text="Style")
        styleFrame.grid(row=0,column=3,sticky='NESW')
        #
        styleFrame.grid_columnconfigure(0,weight=1)
        styleFrame.grid_rowconfigure(0,weight=1)
        #
        B = TTK.Button(styleFrame,text=self.subGraph.grid_property[self.ind_axis].major.y.grid_style,command=lambda n=(ind,self.styleItem): self.bt_click(n))
        B.list = linestylelist
        B.val = self.subGraph.grid_property[self.ind_axis].major.y.grid_style
        B.var = 'grid_style'
        B.ind = ind
        B.treatmentId = 0
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        B.which = ('major','y')
        self.styleItem.append(B)
        ### ### ---> Tick number
        ticknumberFrame = TTK.LabelFrame(ylblframe, text="# Ticks")
        ticknumberFrame.grid(row=0,column=4,sticky='NESW')
        #
        ticknumberFrame.grid_columnconfigure(0,weight=1)
        ticknumberFrame.grid_rowconfigure(0,weight=1)
        #
        B = TTK.Button(ticknumberFrame,text=self.subGraph.grid_property[self.ind_axis].major.y.grid_tick_number,command=lambda n=(ind,self.tickNumberItem): self.bt_click(n))
        B.pattern = '''^[0-9]*$'''
        B.val = self.subGraph.grid_property[self.ind_axis].major.y.grid_tick_number
        B.var = 'grid_tick_number'
        B.ind = ind
        B.treatmentId = 3
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        B.which = ('major','y')
        self.tickNumberItem.append(B)

        ### ### ---> Tick size
        ticksizeFrame = TTK.LabelFrame(ylblframe, text="Ticks Size")
        ticksizeFrame.grid(row=0,column=5,sticky='NESW')
        #
        ticksizeFrame.grid_columnconfigure(0,weight=1)
        ticksizeFrame.grid_rowconfigure(0,weight=1)
        #
        B = TTK.Button(ticksizeFrame,text=self.subGraph.grid_property[self.ind_axis].major.y.grid_tick_size,command=lambda n=(ind,self.tickSizeItem): self.bt_click(n))
        B.list = []
        B.val = self.subGraph.grid_property[self.ind_axis].major.y.grid_tick_size
        B.var = 'grid_tick_size'
        B.ind = ind
        B.treatmentId = 2
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        B.which = ('major','y')
        self.tickSizeItem.append(B)

        ########################################################
        ## Minor grid
        lblframe = TTK.LabelFrame(self.frame,text='Minor grid')
        lblframe.grid(row=2,column=0,sticky="NSEW")
        #
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        lblframe.grid_rowconfigure(1,weight=1)
        ### --> X axis
        ind = 2
        xlblframe = TTK.LabelFrame(lblframe,text='X')
        xlblframe.grid(row=0,column=0,sticky="NSEW")
        #
        xlblframe.grid_columnconfigure(0,weight=1)
        xlblframe.grid_columnconfigure(1,weight=1)
        xlblframe.grid_columnconfigure(2,weight=1)
        xlblframe.grid_columnconfigure(3,weight=1)
        xlblframe.grid_columnconfigure(4,weight=1)
        xlblframe.grid_columnconfigure(5,weight=1)
        xlblframe.grid_rowconfigure(0,weight=1)
        ### ### ---> Display
        displayFrame = TTK.LabelFrame(xlblframe, text="Display")
        displayFrame.grid(row=0,column=0,sticky='NESW')
        #
        displayFrame.grid_columnconfigure(0,weight=1)
        displayFrame.grid_rowconfigure(0,weight=1)
        var = TK.IntVar()
        var.set(self.subGraph.grid_property[self.ind_axis].minor.x.display)
        CB = TTK.Checkbutton(displayFrame,variable=var,command=lambda n=ind: self.cb_visibility(n))#, variable=var)
        CB.val = var
        CB.var = 'display'
        CB.which=('minor','x')
        CB.grid(row=0,column=0,sticky='NSEW')
        self.visibilityItem.append(CB)
#        ### ### ---> Grid color
        colorFrame = TTK.LabelFrame(xlblframe, text="Color")
        colorFrame.grid(row=0,column=1,sticky='NESW')
        #
        colorFrame.grid_columnconfigure(0,weight=1)
        colorFrame.grid_rowconfigure(0,weight=1)
        #
        B = TK.Button(colorFrame,command=lambda n=(ind,self.colorItem): self.bt_click(n))
        B.list = []
        B.val = self.subGraph.grid_property[self.ind_axis].minor.x.grid_color
        B.var = 'grid_color'
        B.config(bg=B.val)
        B.config(activebackground=B.val)
        B.ind = ind
        B.treatmentId = 1
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        B.which = ('minor','x')
        self.colorItem.append(B)
        ### ### ---> Grid width
        widthFrame = TTK.LabelFrame(xlblframe, text="Width")
        widthFrame.grid(row=0,column=2,sticky='NESW')
        #
        widthFrame.grid_columnconfigure(0,weight=1)
        widthFrame.grid_rowconfigure(0,weight=1)
        #
        B = TTK.Button(widthFrame,text=self.subGraph.grid_property[self.ind_axis].minor.x.grid_width,command=lambda n=(ind,self.widthItem): self.bt_click(n))
        B.list = []
        B.val = self.subGraph.grid_property[self.ind_axis].minor.x.grid_width
        B.var = 'grid_width'
        B.ind = ind
        B.treatmentId = 2
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        B.which = ('minor','x')
        self.widthItem.append(B)
        ### ### ---> Grid style
        styleFrame = TTK.LabelFrame(xlblframe, text="Style")
        styleFrame.grid(row=0,column=3,sticky='NESW')
        #
        styleFrame.grid_columnconfigure(0,weight=1)
        styleFrame.grid_rowconfigure(0,weight=1)
        #
        B = TTK.Button(styleFrame,text=self.subGraph.grid_property[self.ind_axis].minor.x.grid_style,command=lambda n=(ind,self.styleItem): self.bt_click(n))
        B.list = linestylelist
        B.val = self.subGraph.grid_property[self.ind_axis].minor.x.grid_style
        B.var = 'grid_style'
        B.ind = ind
        B.treatmentId = 0
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        B.which = ('minor','x')
        self.styleItem.append(B)
        ### ### ---> Tick number
        ticknumberFrame = TTK.LabelFrame(xlblframe, text="# Ticks")
        ticknumberFrame.grid(row=0,column=4,sticky='NESW')
        #
        ticknumberFrame.grid_columnconfigure(0,weight=1)
        ticknumberFrame.grid_rowconfigure(0,weight=1)
        #
        B = TTK.Button(ticknumberFrame,text=self.subGraph.grid_property[self.ind_axis].minor.x.grid_tick_number,command=lambda n=(ind,self.tickNumberItem): self.bt_click(n))
        B.pattern = '''^[0-9]*$'''
        B.val = self.subGraph.grid_property[self.ind_axis].minor.x.grid_tick_number
        B.var = 'grid_tick_number'
        B.ind = ind
        B.treatmentId = 3
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        B.which = ('minor','x')
        self.tickNumberItem.append(B)

# ! WARNING : tick size has no meaning for the minor grid
#             this is why the following is commented

        # ### ### ---> Tick size
        # ticksizeFrame = TTK.LabelFrame(xlblframe, text="Ticks Size")
        # ticksizeFrame.grid(row=0,column=5,sticky='NESW')
        # #
        # ticksizeFrame.grid_columnconfigure(0,weight=1)
        # ticksizeFrame.grid_rowconfigure(0,weight=1)
        # #
        # B = TTK.Button(ticksizeFrame,text=self.subGraph.grid_property[self.ind_axis].minor.x.grid_tick_size,command=lambda n=(ind,self.tickSizeItem): self.bt_click(n))
        # B.list = []
        # B.val = self.subGraph.grid_property[self.ind_axis].minor.x.grid_tick_size
        # B.var = 'grid_tick_size'
        # B.ind = ind
        # B.treatmentId = 2
        # B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        # B.which = ('minor','x')
        # self.tickSizeItem.append(B)


        ### --> Y axis
        ind = 3
        ylblframe = TTK.LabelFrame(lblframe,text='Y')
        ylblframe.grid(row=1,column=0,sticky="NSEW")
        #
        ylblframe.grid_columnconfigure(0,weight=1)
        ylblframe.grid_columnconfigure(1,weight=1)
        ylblframe.grid_columnconfigure(2,weight=1)
        ylblframe.grid_columnconfigure(3,weight=1)
        ylblframe.grid_columnconfigure(4,weight=1)
        ylblframe.grid_columnconfigure(5,weight=1)
        ylblframe.grid_rowconfigure(0,weight=1)
        ### ### ---> Display
        displayFrame = TTK.LabelFrame(ylblframe, text="Display")
        displayFrame.grid(row=0,column=0,sticky='NESW')
        #
        displayFrame.grid_columnconfigure(0,weight=1)
        displayFrame.grid_rowconfigure(0,weight=1)
        var = TK.IntVar()
        var.set(self.subGraph.grid_property[self.ind_axis].minor.y.display)
        CB = TTK.Checkbutton(displayFrame,variable=var,command=lambda n=ind: self.cb_visibility(n))#, variable=var)
        CB.val = var
        CB.var = 'display'
        CB.which = ('minor','y')
        CB.grid(row=0,column=0,sticky='NSEW')
        self.visibilityItem.append(CB)
        ### ### ---> Line color
        colorFrame = TTK.LabelFrame(ylblframe, text="Color")
        colorFrame.grid(row=0,column=1,sticky='NESW')
        #
        colorFrame.grid_columnconfigure(0,weight=1)
        colorFrame.grid_rowconfigure(0,weight=1)
        #
        B = TK.Button(colorFrame,command=lambda n=(ind,self.colorItem): self.bt_click(n))
        B.list = []
        B.val = self.subGraph.grid_property[self.ind_axis].minor.y.grid_color
        B.var = 'grid_color'
        B.config(bg=B.val)
        B.config(activebackground=B.val)
        B.ind = ind
        B.treatmentId = 1
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        B.which = ('minor','y')
        self.colorItem.append(B)
        ### ### ---> Grid width
        widthFrame = TTK.LabelFrame(ylblframe, text="Width")
        widthFrame.grid(row=0,column=2,sticky='NESW')
        #
        widthFrame.grid_columnconfigure(0,weight=1)
        widthFrame.grid_rowconfigure(0,weight=1)
        #
        B = TTK.Button(widthFrame,text=self.subGraph.grid_property[self.ind_axis].minor.y.grid_width,command=lambda n=(ind,self.widthItem): self.bt_click(n))
        B.list = []
        B.val = self.subGraph.grid_property[self.ind_axis].minor.y.grid_width
        B.var = 'grid_width'
        B.ind = ind
        B.treatmentId = 2
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        B.which = ('minor','y')
        self.widthItem.append(B)
        ### ### ---> Grid style
        styleFrame = TTK.LabelFrame(ylblframe, text="Style")
        styleFrame.grid(row=0,column=3,sticky='NESW')
        #
        styleFrame.grid_columnconfigure(0,weight=1)
        styleFrame.grid_rowconfigure(0,weight=1)
        #
        B = TTK.Button(styleFrame,text=self.subGraph.grid_property[self.ind_axis].minor.y.grid_style,command=lambda n=(ind,self.styleItem): self.bt_click(n))
        B.list = linestylelist
        B.val = self.subGraph.grid_property[self.ind_axis].minor.y.grid_style
        B.var = 'grid_style'
        B.ind = ind
        B.treatmentId = 0
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        B.which = ('minor','y')
        self.styleItem.append(B)
        ### ### ---> Tick number
        ticknumberFrame = TTK.LabelFrame(ylblframe, text="# Ticks")
        ticknumberFrame.grid(row=0,column=4,sticky='NESW')
        #
        ticknumberFrame.grid_columnconfigure(0,weight=1)
        ticknumberFrame.grid_rowconfigure(0,weight=1)
        #
        B = TTK.Button(ticknumberFrame,text=self.subGraph.grid_property[self.ind_axis].minor.y.grid_tick_number,command=lambda n=(ind,self.tickNumberItem): self.bt_click(n))
        B.pattern = '''^[0-9]*$'''
        B.val = self.subGraph.grid_property[self.ind_axis].minor.y.grid_tick_number
        B.var = 'grid_tick_number'
        B.ind = ind
        B.treatmentId = 3
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        B.which = ('minor','y')
        self.tickNumberItem.append(B)

# ! WARNING : tick size has no meaning for the minor grid
#             this is why the following is commented

        # ### ### ---> Tick size
        # ticksizeFrame = TTK.LabelFrame(ylblframe, text="Ticks Size")
        # ticksizeFrame.grid(row=0,column=5,sticky='NESW')
        # #
        # ticksizeFrame.grid_columnconfigure(0,weight=1)
        # ticksizeFrame.grid_rowconfigure(0,weight=1)
        # #
        # B = TTK.Button(ticksizeFrame,text=self.subGraph.grid_property[self.ind_axis].minor.y.grid_tick_size,command=lambda n=(ind,self.tickSizeItem): self.bt_click(n))
        # B.list = []
        # B.val = self.subGraph.grid_property[self.ind_axis].minor.y.grid_tick_size
        # B.var = 'grid_tick_size'
        # B.ind = ind
        # B.treatmentId = 2
        # B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        # B.which = ('minor','y')
        # self.tickSizeItem.append(B)

    # ----------------------------------------------------------- closeAllDialog
    def closeAllDialog(self):
        if self.list_dialog:
            self.list_dialog.cmd_close()
            self.list_dialog = None
        if self.input_dialog:
            self.input_dialog.cmd_close()
            self.input_dialog = None
    # ----------------------------------------------------------------- updateColor
    def updateColor(self,color,B,extra_data):
        if color is not None:
            B.val = color
            B.config(bg=B.val)
            B.config(activebackground=B.val)
            try:
                self.subGraph.grid_property[self.ind_axis].setValue(B.which[0],B.which[1],'grid_color',B.val)
                # Update Graph
                self.parent.graphWdwL[self.graph].updateGraph(self.zone)
            except IndexError:
                return
    # ----------------------------------------------------------------- bt_click
    def bt_click(self,ind):
        bt_list = ind[1]
        B = bt_list[ind[0]]
        self.closeAllDialog()
        if B.treatmentId==0: # List
            self.list_dialog = list_dialogWindow()
            self.list_dialog.initialize(self,B)
        elif B.treatmentId==1: # Color selector
            self.list_dialog = color_dialogWindow(self,B,B.val)
        elif B.treatmentId==2: # Float
            self.input_dialog = inputFloat_dialogWindow()
            self.input_dialog.initialize(self,B)
        elif B.treatmentId==3: # Pattern
            self.input_dialog = inputPattern_dialogWindow(B.pattern)
            self.input_dialog.initialize(self,B)
        else:
            return
    # ------------------------------------------------------------ cb_visibility
    def cb_visibility(self,ind):
        CB = self.visibilityItem[ind]
        self.subGraph.grid_property[self.ind_axis].setValue(CB.which[0],CB.which[1],'display',bool(CB.val.get()))
        # Update Graph
        self.parent.graphWdwL[self.graph].updateGraph(self.zone)
    # ------------------------------------------------------------- reloadWindow
    def reloadWindow(self):
        self.frame.destroy()
        self.createFrame()
    # ---------------------------------------------------------------- cmd_close
    def cmd_close(self):
        self.closeAllDialog()
        self.parent.editGridWdw = None
        self.destroy()
    # --------------------------------------------------------- updateButonWidth
    def updateButonWidth(self,B,val):
        B.val = val
        B.config(text=B.val)
        try:
            self.subGraph.grid_property[self.ind_axis].setValue(B.which[0],B.which[1],B.var,B.val)
            # Update Graph
            self.parent.graphWdwL[self.graph].updateGraph(self.zone)
        except IndexError:
            return
    # # --------------------------------------------------------- updateButonWidth
    # def updateButonWidth(self,B,val):
    #     B.val = val
    #     B.config(text=B.val)
    #     try:
    #         self.subGraph.grid_property[self.ind_axis].setValue(B.which[0],B.which[1],'grid_width',B.val)
    #         # Update Graph
    #         self.parent.graphWdwL[self.graph].updateGraph(self.zone)
    #     except IndexError:
    #         return
    # --------------------------------------------------------- updateButon
    def updateButon(self,B,val):
        B.val = val
        B.config(text=B.val)
        if B.var == 'axis':
            self.ind_axis = B.val
            self.createFrame()

        else:
            try:
                self.subGraph.grid_property[self.ind_axis].setValue(B.which[0],B.which[1],B.var,B.val)
                # Update Graph
                self.parent.graphWdwL[self.graph].updateGraph(self.zone)
            except IndexError:
                return
    # # --------------------------------------------------------- updateButon
    # def updateButon(self,B,val):
    #     B.val = val
    #     B.config(text=B.val)
    #     if B.var == 'axis':
    #         self.ind_axis = B.val
    #         self.createFrame()
    #
    #     else:
    #         try:
    #             self.subGraph.grid_property[self.ind_axis].setValue(B.which[0],B.which[1],'grid_style',B.val)
    #             # Update Graph
    #             self.parent.graphWdwL[self.graph].updateGraph(self.zone)
    #         except IndexError:
    #             return

# ==============================================================================
# ==============================================================================
class input_dialogSelectZoneWindow(TK.Toplevel):
    # --------------------------------------------------------------------- init
    def __init__(self):
        TK.Toplevel.__init__(self)
        (x,y) = self.winfo_pointerxy()
        self.geometry('+%s+%s'%(x,y))
        #self.protocol("WM_DELETE_WINDOW",self.cmd_close)
    # --------------------------------------------------------------- initialize
    def initialize(self,parent,B):
        self.B = B
        self.parent=parent
        self.grid_columnconfigure(0,weight=1)
        self.grid_rowconfigure(0,weight=0)
        self.grid_rowconfigure(1,weight=4)
        self.grid_rowconfigure(2,weight=1)
        ########################################################################
        # # # First line
        label = TTK.LabelFrame(self,text='Filter zone(s)')
        label.grid(row=0,column=0,sticky='EW')
        label.grid_columnconfigure(0,weight=1)
        label.grid_rowconfigure(0,weight=1)
        # # # First line -> Second column
        self.entry_filterSV = TK.StringVar()
        self.entry_filter = TK.Entry(label,textvariable=self.entry_filterSV)
        self.entry_filter.val = ''
        self.entry_filter.insert(TK.END, self.entry_filter.val)
#        self.entry_filterSV.trace('w',lambda nm, idx, mode,var=self.marker_sampling_startSV: self.cmd_editStep(var,self.marker_sampling_start))
        self.entry_filterSV.trace('w',lambda nm, idx, mode,var=self.entry_filterSV: self.cmd_editFilter(var,self.entry_filter))
        self.entry_filter.grid(row=0,column=0,sticky='NSEW')
        re_filter = re.compile(self.entry_filterSV.get())
        ########################################################################
        # # # Second line
        frame = TTK.Frame(self)
        frame.grid(row=1,column=0,sticky='NSEW')
        frame.grid_columnconfigure(0,weight=1)
        frame.grid_columnconfigure(1,weight=1)
        frame.grid_columnconfigure(2,weight=1)
        frame.grid_rowconfigure(0,weight=1)
        frame.grid_rowconfigure(1,weight=1)
        frame.grid_rowconfigure(2,weight=1)
        frame.grid_rowconfigure(3,weight=1)
        # # # Second line -> first column
        lblframe = TTK.LabelFrame(frame,text='Unused zone(s)')
        lblframe.grid(row=0,column=0,rowspan=4,sticky='NSEW')
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        #
        # Create list
        self.unused_zones = TK.Listbox(lblframe,selectmode=TK.EXTENDED)
        self.unused_zones.list = []
        for zone in self.getUnusedZoneList(self.B.val,self.parent.parent.data.keys()):
            self.unused_zones.list.append(zone)
        # Add elements of list to the listbox according to the filter
        self.displayUnusedZones(self.entry_filter.val)
        # Postionning in grid
        self.unused_zones.grid(row=0,column=0,sticky="NSEW")

        # # # Second line -> second column
        # # # Second line -> second column : First button
        Button = TTK.Button(frame,text='>>',command=self.cmd_addAll)
        Button.grid(row=0,column=1,sticky='NSEW')
        # # # Second line -> second column : Second button
        Button = TTK.Button(frame,text='>',command=self.cmd_addSelection)
        Button.grid(row=1,column=1,sticky='NSEW')
        # # # Second line -> second column : Third button
        Button = TTK.Button(frame,text='<',command=self.cmd_removeSelection)
        Button.grid(row=2,column=1,sticky='NSEW')
        # # # Second line -> second column : First button
        Button = TTK.Button(frame,text='<<',command=self.cmd_removeAll)
        Button.grid(row=3,column=1,sticky='NSEW')

        # # # Second line -> third column
        lblframe = TTK.LabelFrame(frame,text='Zone(s) to plot')
        lblframe.grid(row=0,column=2,rowspan=4,sticky='NSEW')
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        #
        self.used_zones = TK.Listbox(lblframe,selectmode=TK.EXTENDED)
        # Create list
        self.used_zones.list = []
        for zone in self.B.val:
            self.used_zones.list.append(zone)
        # Add elements of list to the listbox according to the filter
        self.displayUsedZones(self.entry_filter.val)
        # Postionning in grid
        self.used_zones.grid(row=0,column=0,sticky="NSEW")

        ########################################################################
        # # # Third line
        frame = TTK.Frame(self)
        frame.grid(row=2,column=0,sticky='NSEW')
        frame.grid_columnconfigure(0,weight=1)
        frame.grid_columnconfigure(1,weight=1)
        frame.grid_rowconfigure(0,weight=1)
        #
        Button = TTK.Button(frame,text='OK',command=self.cmd_selectZones)
        Button.grid(row=0,column=0,sticky='NSEW')
        #
        Button = TTK.Button(frame,text='Close',command=self.cmd_close)
        Button.grid(row=0,column=1,sticky='NSEW')

    # -------------------------------------------------------- getUnusedZoneList
    def getUnusedZoneList(self,used,zones):
        res = []
        for zone in zones:
            if not zone in used:
                res.append(zone)
        return res
    # --------------------------------------------------------- displayUsedZones
    def displayUsedZones(self,filter_value):
        self.used_zones.delete(0,TK.END)
        if filter_value.strip() == '':
            for zone in self.used_zones.list:
                self.used_zones.insert(TK.END,zone)
        else:
            re_filter = re.compile(filter_value)
            for zone in self.used_zones.list:
                if re.search(re_filter,zone):
                    self.used_zones.insert(TK.END,zone)
    # ------------------------------------------------------- displayUnusedZones
    def displayUnusedZones(self,filter_value):
        self.unused_zones.delete(0,TK.END)
        if filter_value.strip() == '':
            for zone in self.unused_zones.list:
                self.unused_zones.insert(TK.END,zone)
        else:
            re_filter = re.compile(filter_value)
            for zone in self.unused_zones.list:
                if re.search(re_filter,zone):
                    self.unused_zones.insert(TK.END,zone)
    # ---------------------------------------------------------- cmd_selectZones
    def cmd_selectZones(self):
        # Select zones for curve
        self.B.val = self.used_zones.list
        # Update edit curves window
        self.parent.updateButon(self.B,self.used_zones.list)
        # Close window
        self.cmd_close()

    # --------------------------------------------------------------- cmd_addAll
    def cmd_addAll(self):
        # Add all unused to used list
        for zone in self.unused_zones.get(0,TK.END):
            # Find index in unused listbox
            try:
                index = self.unused_zones.list.index(zone)
            except ValueError:
                continue
            self.used_zones.list.append(zone)
            del self.unused_zones.list[index]
        # Display listbox
        self.displayUsedZones(self.entry_filter.val)
        self.displayUnusedZones(self.entry_filter.val)

    # --------------------------------------------------------- cmd_addSelection
    def cmd_addSelection(self):
        items = map(int, self.unused_zones.curselection())
        zones = []
        for index in items:
            zones.append(self.unused_zones.get(index))
        #
        for zone in zones:
        # Find index in unused listbox
            try:
                index = self.unused_zones.list.index(zone)
            except ValueError:
                continue
            self.used_zones.list.append(zone)
            del self.unused_zones.list[index]
        # Display listbox
        self.displayUsedZones(self.entry_filter.val)
        self.displayUnusedZones(self.entry_filter.val)

    # ------------------------------------------------------------ cmd_removeAll
    def cmd_removeAll(self):
        # Add all used to used list
        for zone in self.used_zones.get(0,TK.END):
        # Find index in used listbox
            try:
                index = self.used_zones.list.index(zone)
            except ValueError:
                continue
            self.unused_zones.list.append(zone)
            del self.used_zones.list[index]
        # Display listbox
        self.displayUsedZones(self.entry_filter.val)
        self.displayUnusedZones(self.entry_filter.val)
    # ------------------------------------------------------ cmd_removeSelection
    def cmd_removeSelection(self):
        items = map(int, self.used_zones.curselection())
        zones = []
        for index in items:
            zones.append(self.used_zones.get(index))
        #
        for zone in zones:
        # Find index in used listbox
            try:
                index = self.used_zones.list.index(zone)
            except ValueError:
                continue
            self.unused_zones.list.append(zone)
            del self.used_zones.list[index]
        # Display listbox
        self.displayUsedZones(self.entry_filter.val)
        self.displayUnusedZones(self.entry_filter.val)
    # ------------------------------------------------------------- cmd_editStep
    # The following function ensures that only integer is used even if defined by user ...
    def cmd_editFilter(self,var,widget):
        widget.val = var.get()
        if widget.val.strip() =='':
            # self.unused_zones
            self.displayUnusedZones(widget.val)
            # self.used_zones
            self.displayUsedZones(widget.val)
        else:
            re_filter = re.compile(widget.val)
            # self.unused_zones
            self.unused_zones.delete(0,TK.END)
            for zone in self.unused_zones.list:
                if re.search(re_filter,zone):
                    self.unused_zones.insert(TK.END,zone)
            # self.used_zones
            self.used_zones.delete(0,TK.END)
            for zone in self.used_zones.list:
                if re.search(re_filter,zone):
                    self.used_zones.insert(TK.END,zone)
            # Display listbox
            self.displayUsedZones(widget.val)
            self.displayUnusedZones(widget.val)
    # ------------------------------------------------------------- cmd_okButton
    def cmd_okButton(self,event=None):
        valString = self.entry.get()
        try:
            valFloat = float(valString)
            self.parent.updateButonWidth(self.B,valFloat)
            self.cmd_close()
        except ValueError:
            self.entry.delete(0,TK.END)
            self.entry.insert(TK.END, self.B.val)
            return
    # ---------------------------------------------------------------- cmd_close
    def cmd_close(self,event=None):
        self.destroy()
# ==============================================================================
# ==============================================================================
class inputPattern_dialogWindow(TK.Toplevel):
    # --------------------------------------------------------------------- init
    def __init__(self,pattern):
        TK.Toplevel.__init__(self)
        (x,y) = self.winfo_pointerxy()
        self.pattern = pattern
        self.geometry('+%s+%s'%(x,y))
        #self.protocol("WM_DELETE_WINDOW",self.cmd_close)
    # --------------------------------------------------------------- initialize
    def initialize(self,parent,B):
        self.B = B
        self.parent=parent
        self.grid_columnconfigure(0,weight=1)
        self.grid_rowconfigure(0,weight=1)
        self.grid_rowconfigure(1,weight=1)
        #
        lblframe = TTK.LabelFrame(self,text='Enter value')
        lblframe.grid(row=0,column=0,sticky='NSEW')
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        #
        self.entry = TK.Entry(lblframe)
        if B.val is not None:
            self.entry.insert(TK.END, B.val)
        self.entry.grid(row=0,column=0,columnspan=1,sticky="NSEW")
        #
        frame = TTK.Frame(self)
        frame.grid(row=1,column=0,sticky="NSEW")
        frame.grid_columnconfigure(0,weight=1)
        frame.grid_columnconfigure(1,weight=1)
        frame.grid_rowconfigure(0,weight=1)

        Button = TTK.Button(frame,text='OK',command=self.cmd_okButton)
        Button.grid(row=0,column=0,sticky='NSEW')

        Button = TTK.Button(frame,text='Cancel',command=self.cmd_close)
        Button.grid(row=0,column=1,sticky='NSEW')
        #
        self.entry.focus_set()
        self.entry.select_from(0)
        self.entry.select_to(TK.END)
        self.bind('<Return>', self.cmd_okButton)
    # ------------------------------------------------------------- cmd_okButton
    def cmd_okButton(self,event=None):
        valString = self.entry.get()
        if re.match(self.pattern,valString):
            self.parent.updateButonWidth(self.B,valString.replace(" ",""))
            self.cmd_close()
        else:
            self.entry.delete(0,TK.END)
            if self.B.val is not None:
                self.entry.insert(TK.END, self.B.val)
            return
    # ---------------------------------------------------------------- cmd_close
    def cmd_close(self):
        self.destroy()
# ==============================================================================
# ==============================================================================
class inputFloat_dialogWindow(TK.Toplevel):
    # --------------------------------------------------------------------- init
    def __init__(self):
        TK.Toplevel.__init__(self)
        (x,y) = self.winfo_pointerxy()
        self.geometry('+%s+%s'%(x,y))
        #self.protocol("WM_DELETE_WINDOW",self.cmd_close)
    # --------------------------------------------------------------- initialize
    def initialize(self,parent,B):
        self.B = B
        self.parent=parent
        self.grid_columnconfigure(0,weight=1)
        self.grid_rowconfigure(0,weight=1)
        self.grid_rowconfigure(1,weight=1)
        #
        lblframe = TTK.LabelFrame(self,text='Enter value')
        lblframe.grid(row=0,column=0,sticky='NSEW')
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        #
        self.entry = TK.Entry(lblframe)
        if B.val is not None:
            self.entry.insert(TK.END, B.val)
        self.entry.grid(row=0,column=0,columnspan=1,sticky="NSEW")
        #
        frame = TTK.Frame(self)
        frame.grid(row=1,column=0,sticky="NSEW")
        frame.grid_columnconfigure(0,weight=1)
        frame.grid_columnconfigure(1,weight=1)
        frame.grid_rowconfigure(0,weight=1)

        Button = TTK.Button(frame,text='OK',command=self.cmd_okButton)
        Button.grid(row=0,column=0,sticky='NSEW')

        Button = TTK.Button(frame,text='Cancel',command=self.cmd_close)
        Button.grid(row=0,column=1,sticky='NSEW')
        #
        self.entry.focus_set()
        self.entry.select_from(0)
        self.entry.select_to(TK.END)
        self.bind('<Return>', self.cmd_okButton)
    # ------------------------------------------------------------- cmd_okButton
    def cmd_okButton(self,event=None):
        valString = self.entry.get()
        try:
            valFloat = float(valString)
            self.parent.updateButonWidth(self.B,valFloat)
            self.cmd_close()
        except ValueError:
            self.entry.delete(0,TK.END)
            if self.B.val is not None:
                self.entry.insert(TK.END, self.B.val)
            return
    # ---------------------------------------------------------------- cmd_close
    def cmd_close(self):
        self.destroy()
# ==============================================================================
# ==============================================================================
class input_dialogStringWindow(TK.Toplevel):
    # --------------------------------------------------------------------- init
    def __init__(self):
        TK.Toplevel.__init__(self)
        (x,y) = self.winfo_pointerxy()
        self.geometry('+%s+%s'%(x,y))
        #self.protocol("WM_DELETE_WINDOW",self.cmd_close)
    # --------------------------------------------------------------- initialize
    def initialize(self,parent,B):
        self.B = B
        self.parent=parent
        self.grid_columnconfigure(0,weight=1)
        self.grid_rowconfigure(0,weight=1)
        self.grid_rowconfigure(1,weight=1)
        #
        lblframe = TTK.LabelFrame(self,text='Enter value')
        lblframe.grid(row=0,column=0,sticky='NSEW')
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        #
        self.entry = TK.Entry(lblframe)
        self.entry.insert(TK.END, B.val)
        self.entry.grid(row=0,column=0,columnspan=1,sticky="NSEW")
        #
        frame = TTK.Frame(self)
        frame.grid(row=1,column=0,sticky="NSEW")
        frame.grid_columnconfigure(0,weight=1)
        frame.grid_columnconfigure(1,weight=1)
        frame.grid_rowconfigure(0,weight=1)

        Button = TTK.Button(frame,text='OK',command=self.cmd_okButton)
        Button.grid(row=0,column=0,sticky='NSEW')

        Button = TTK.Button(frame,text='Cancel',command=self.cmd_close)
        Button.grid(row=0,column=1,sticky='NSEW')
        #
        self.entry.focus_set()
        self.entry.select_from(0)
        self.entry.select_to(TK.END)
        self.bind('<Return>', self.cmd_okButton)
    # ------------------------------------------------------------- cmd_okButton
    def cmd_okButton(self,event=None):
        valString = self.entry.get()
        try:
            valStr = valString
            self.parent.updateButon(self.B,valString)
            self.cmd_close()
        except ValueError:
            self.entry.delete(0,TK.END)
            self.entry.insert(TK.END, self.B.val)
            return
    # ---------------------------------------------------------------- cmd_close
    def cmd_close(self):
        self.destroy()
# ==============================================================================
# ==============================================================================
class input_dialogWindow(TK.Toplevel):
    # --------------------------------------------------------------------- init
    def __init__(self):
        TK.Toplevel.__init__(self)
        (x,y) = self.winfo_pointerxy()
        self.geometry('+%s+%s'%(x,y))
        #self.protocol("WM_DELETE_WINDOW",self.cmd_close)
    # --------------------------------------------------------------- initialize
    def initialize(self,parent,B):
        self.B = B
        self.parent=parent
        self.grid_columnconfigure(0,weight=1)
        self.grid_rowconfigure(0,weight=1)
        self.grid_rowconfigure(1,weight=1)
        #
        lblframe = TTK.LabelFrame(self,text='Enter value')
        lblframe.grid(row=0,column=0,sticky='NSEW')
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        #
        self.entry = TK.Entry(lblframe)
        self.entry.insert(TK.END, B.val)
        self.entry.grid(row=0,column=0,columnspan=1,sticky="NSEW")
        #
        frame = TTK.Frame(self)
        frame.grid(row=1,column=0,sticky="NSEW")
        frame.grid_columnconfigure(0,weight=1)
        frame.grid_columnconfigure(1,weight=1)
        frame.grid_rowconfigure(0,weight=1)

        Button = TTK.Button(frame,text='OK',command=self.cmd_okButton)
        Button.grid(row=0,column=0,sticky='NSEW')

        Button = TTK.Button(frame,text='Cancel',command=self.cmd_close)
        Button.grid(row=0,column=1,sticky='NSEW')
        #
        self.entry.focus_set()
        self.entry.select_from(0)
        self.entry.select_to(TK.END)
        self.bind('<Return>', self.cmd_okButton)
    # ------------------------------------------------------------- cmd_okButton
    def cmd_okButton(self,event=None):
        valString = self.entry.get()
        try:
            valInt = int(valString)
            self.parent.updateButon(self.B,valInt)
            self.cmd_close()
        except ValueError:
            self.entry.delete(0,TK.END)
            self.entry.insert(TK.END, self.B.val)
            return
    # ---------------------------------------------------------------- cmd_close
    def cmd_close(self):
        self.destroy()
# ==============================================================================
# ==============================================================================
class list_dialogWindow(TK.Toplevel):
    # --------------------------------------------------------------------- init
    def __init__(self):
        TK.Toplevel.__init__(self)
        (x,y) = self.winfo_pointerxy()
        self.geometry('+%s+%s'%(x,y))
        #self.protocol("WM_DELETE_WINDOW", self.cmd_close)
    # --------------------------------------------------------------- initialize
    def initialize(self,parent,B):
        self.optionList = B.list
        curValue = B.val
        try:
            indCurValue = self.optionList.index(curValue)
        except ValueError:
            indCurValue = 0
        self.parent = parent
        self.B = B
        #
        self.grid_columnconfigure(0,weight=1)
        self.grid_rowconfigure(0,weight=0)
        self.grid_rowconfigure(1,weight=1)
        #
        filterframe = TTK.LabelFrame(self,text='Filter values')
        filterframe.grid(row=0,column=0,sticky="EW")
        filterframe.grid_columnconfigure(0,weight=1)
        filterframe.grid_rowconfigure(0,weight=1)
        #
        self.entry_filterSV = TK.StringVar()
        self.entry_filter = TK.Entry(filterframe,textvariable=self.entry_filterSV)
        self.entry_filter.val = ''
        self.entry_filter.insert(TK.END, self.entry_filter.val)
        self.entry_filterSV.trace('w',lambda nm, idx, mode,var=self.entry_filterSV: self.cmd_editFilter(var,self.entry_filter))
        self.entry_filter.grid(row=0,column=0,sticky='NSEW')
        re_filter = re.compile(self.entry_filterSV.get())
        #
        self.lblframe = TTK.LabelFrame(self,text='Select value')
        self.lblframe.grid(row=1,column=0,sticky="NSEW")
        self.lblframe.grid_columnconfigure(0,weight=1)
        self.lblframe.grid_columnconfigure(1,weight=0)
        self.lblframe.grid_rowconfigure(0,weight=1)
        #
        self.scrollbar = TK.Scrollbar(self.lblframe)
        self.scrollbar.grid(row=0,column=1,sticky='NS')
        #
        self.lb = TK.Listbox(self.lblframe,selectmode=TK.SINGLE,yscrollcommand=self.scrollbar.set)
        self.scrollbar['command'] = self.lb.yview
        for opt in self.optionList:
            self.lb.insert(TK.END,opt)
        self.lb.select_set(indCurValue)
        self.lb.see(indCurValue)
        self.lb.bind("<<ListboxSelect>>",self.cmd_onSelectChange)
        self.lb.grid(row=0,column=0,sticky="NSEW")
        # #
        # self.lblframe = TTK.LabelFrame(self,text='Select value')
        # self.lblframe.grid(row=1,column=0,sticky="NSEW")
        # self.lblframe.grid_columnconfigure(0,weight=1)
        # self.lblframe.grid_rowconfigure(0,weight=1)
        # #
        # self.lb = TK.Listbox(self.lblframe,selectmode=TK.SINGLE)
        # for opt in self.optionList:
        #     self.lb.insert(TK.END,opt)
        # self.lb.select_set(indCurValue)
        # self.lb.see(indCurValue)
        # self.lb.bind("<<ListboxSelect>>",self.cmd_onSelectChange)
        # self.lb.grid(row=0,column=0,sticky="NSEW")
    # ------------------------------------------------------------- cmd_editStep
    # The following function ensures that only integer is used even if defined by user ...
    def cmd_editFilter(self,var,widget):
        widget.val = var.get()
        self.lb.delete(0,TK.END)
        if widget.val.strip() =='':
            # display all list
            for opt in self.optionList:
                self.lb.insert(TK.END,opt)
        else:
            # display list matching re
            re_filter = re.compile(widget.val)
            # self.unused_zones
            for opt in self.optionList:
                if re.search(re_filter,opt):
                    self.lb.insert(TK.END,opt)
    # ------------------------------------------------------- cmd_onSelectChange
    def cmd_onSelectChange(self,event):
        widget = event.widget
        index = self.lb.curselection()[0]
        newval = self.lb.get(index)
        self.parent.updateButon(self.B,newval)
        self.cmd_close()
    # ---------------------------------------------------------------- cmd_close
    def cmd_close(self):
        self.destroy()
# ==============================================================================
# ==============================================================================

class GraphTK(TK.Toplevel):
    # --------------------------------------------------------------------- init
    def __init__(self,parent,name,conf,dpi=None,figsize=None):
        TK.Toplevel.__init__(self)
        self.dpi = dpi
        self.figsize = figsize
        self.minsize(width=500, height=500)
        self.protocol("WM_DELETE_WINDOW",self.cmd_close)
        self.bind("<FocusIn>",self.clickOnWindow)
        self.grid_columnconfigure(0,weight=1)
        self.grid_rowconfigure(0,weight=1)
        self.parent = parent
        self.name=self.soleName(name)
        self.conf=conf
        self.initialize()
        self.subPlotParams = SubPlotParams()
        self.tightLayout = TightLayout()
        self.image_background_color = default_values['Graph']['image_background_color']
        self.image_background_alpha = default_values['Graph']['image_background_alpha']
        self.subgraph_background_color = default_values['Graph']['subgraph_background_color']
        self.subgraph_background_alpha = default_values['Graph']['subgraph_background_alpha']
        # update la data a l'ouverture du graph
        if IMPORTOK and CTK is not None and CTK.t != [] and len(self.parent.graphWdwL) == 1:
            updateFromTree()

    # --------------------------------------------------------------- initialize
    def initialize(self):
        self.title(string=self.name)
        self.parent.updateGraphName2Id()
        self.parent.graphName2Id[self.name]=len(self.parent.graphWdwL)
        self.index = len(self.parent.graphWdwL)
        self.parent.graphWdwL.append(self)
        # Create a frame
        self.frame = TK.Frame(self) # Can not use TTK here because of line : self.frame.configure(bg="white")
        self.frame.configure(bg="white")
        self.frame.grid(row=0,column=0,sticky="nsew")
        self.frame.grid_columnconfigure(0,weight=1)
        self.frame.grid_rowconfigure(0,weight=1)
        self.frame.grid_rowconfigure(1,weight=0)
        # Figure
        self.fig = MatplotlibFigure(self.parent,self,self.conf,self.dpi,self.figsize)


        # Canvas
        self.canvas = FigureCanvasTkAgg(self.fig.instance,self.frame)
        self.canvas.get_tk_widget().grid_columnconfigure(0,weight=1)
        self.canvas.get_tk_widget().grid_rowconfigure(0,weight=1)
        self.canvas.get_tk_widget().grid(row=0,column=0,sticky="NSEW")
        self.canvas.draw()
        # ax click binding
        self.canvas.mpl_connect('button_press_event',self.clickOnCanvas)


        # Toolbar
        toolbar_frame = TTK.Frame(self)
        toolbar_frame.grid(row=1,column=0,sticky='W')
#        toolbar = NavigationToolbar2TkAgg(self.canvas,toolbar_frame)
        toolbar = CustomToolbar(self.canvas,toolbar_frame,self)
        toolbar.update()

        # Update Main window : graph name list
        self.parent.updateactiveGraph()



        # Configure une zone graphique
    # ------------------------------------------------------------------ getFig
    def getFig(self):
        return self.fig.getFig()
    # ------------------------------------------------------ deleteZoneInCurve
    def deleteZoneInCurve(self,iCurSubGraph,zoneName):
        self.fig.deleteZoneInCurve(iCurSubGraph,zoneName)
    # ------------------------------------------------------ updateGroupCurvesZoneName
    def updateGroupCurvesZoneName(self,iCurSubGraph,oldZoneList,newZoneList):
        self.fig.updateGroupCurvesZoneName(iCurSubGraph,oldZoneList,newZoneList)

    # ------------------------------------------------------ updateSubPlotParams
    def updateSubPlotParams(self,params):
        for var in params:
            self.subPlotParams.setValue(var,params[var])
            if var == 'isActive':
                self.tightLayout.setValue(var,not params[var])
        # Update Figure
        self.applyViewSettings()
    # -------------------------------------------------------- updateTightLayout
    def updateTightLayout(self,params):
        for var in params:
            self.tightLayout.setValue(var,params[var])
            if var == 'isActive':
                self.subPlotParams.setValue(var,not params[var])
        # Update Figure
        self.applyViewSettings()
    # --------------------------------------------------------------- showFigure
    def showFigure(self):
        # We do nothing, but the function needs to exist for interface compatibility with the without TK interface
        return
    # --------------------------------------------------------------- drawFigure
    def drawFigure(self):
        # We do nothing, but the function needs to exist for interface compatibility with the without TK interface
        return
    # --------------------------------------------------------------------- save
    def save(self,path,format=None):
        self.fig.saveFigure(path,format=format)
    # ------------------------------------------------------------------ addGrid
    def getGrid(self,iCurSubGraph,ind=0,axis=None):
        if axis:
            ind = axis.getInd()
        #
        return self.fig.getGrid(iCurSubGraph,ind)
    # ------------------------------------------------------------------ addAxis
    def addAxis(self,iCurSubGraph,shared=None,ind=0,axis=None):
        if axis :
            ind = axis.getInd()
        axis_to_twin = ind
        if shared is None:
            self.fig.subGraph[iCurSubGraph].addAxis()
            newind = len(self.fig.subGraph[iCurSubGraph].axis)-1
            newaxis = self.fig.getAxis(iCurSubGraph,ind=newind)
            newaxis.ind = newind
            return newaxis
        elif ((shared == 'x') or (shared == 'X')):
            self.fig.subGraph[iCurSubGraph].addAxisTwinX(axis_to_twin)
            newind = len(self.fig.subGraph[iCurSubGraph].axis)-1
            newaxis = self.fig.getAxis(iCurSubGraph,ind=newind)
            newaxis.ind = newind
            return newaxis
        elif ((shared == 'y') or (shared == 'Y')):
            self.fig.subGraph[iCurSubGraph].addAxisTwinY(axis_to_twin)
            newind = len(self.fig.subGraph[iCurSubGraph].axis)-1
            newaxis = self.fig.getAxis(iCurSubGraph,ind=newind)
            newaxis.ind = newind
            return newaxis
        else:
            print('''### Error : value used for 'shared' is unknown, must be in [None, 'x', 'X', 'y', 'Y'].''')
    # ---------------------------------------------------------------- getLegend
    def getLegend(self,iCurSubGraph):
        return self.fig.getLegend(iCurSubGraph)
    # ------------------------------------------------------------------ getAxis
    def getAxis(self,iCurSubGraph,ind=0):
        return self.fig.getAxis(iCurSubGraph,ind)
    # -------------------------------------------------------------------- write
    def write(self,ind):
        line = '''graph_%s=GraphTK(obj,'%s','%s',dpi=%s,figsize=%s)\n'''%(ind,self.name,self.conf,self.dpi,self.figsize)
        if self.image_background_color != default_values['Graph']['image_background_color']:
            line += '''graph_%s.setValue('image_background_color',%s)\n'''%(ind,self.image_background_color)
        if self.image_background_alpha != default_values['Graph']['image_background_alpha']:
            line += '''graph_%s.setValue('image_background_alpha',%s)\n'''%(ind,self.image_background_alpha)
        if self.subgraph_background_color != default_values['Graph']['subgraph_background_color']:
            line += '''graph_%s.setValue('subgraph_background_color',%s)\n'''%(ind,self.subgraph_background_color)
        if self.subgraph_background_alpha != default_values['Graph']['subgraph_background_alpha']:
            line += '''graph_%s.setValue('subgraph_background_alpha',%s)\n'''%(ind,self.subgraph_background_alpha)
        return line
    # ------------------------------------------------------------ clickOnCanvas
    def clickOnCanvas(self,event):
        try:
            ax = event.inaxes
            self.parent.selectPositionByName(ax.name)
        except AttributeError: # Click was not on an axe
            return
    # ------------------------------------------------------------ clickOnWindow
    def clickOnWindow(self,event):
        self.parent.selectGraphByName(self.name)
    # ----------------------------------------------------------------- soleName
    def soleName(self,name):
        graphNameList = [g.name for g in self.parent.graphWdwL]
        if not name in graphNameList:
            return name
        else:
            ind=0
            while name+'.%s'%ind in graphNameList:
                ind += 1
            return name+'.%s'%ind

    # ---------------------------------------------------------------- cmd_close
    def cmd_close(self):
        # If we try to close the last graph, then we have to close all the edit windows to avoid their update (that will fail)
        if len(self.parent.graphWdwL) == 1:
            # Close all edit windows if we are closing the last graph:
            for w in [self.parent.editCurveWdw,self.parent.editGridWdw,self.parent.editLegendWdw,self.parent.editAxisWdw,self.parent.editGraphWdw]:
                if w is not None:
                    w.cmd_close()
        #
        self.fig.closeMovie()
        del self.fig
        del self.parent.graphWdwL[self.index]
        self.parent.renumberGraph()
        self.parent.updateactiveGraph()
        self.destroy()
    # ----------------------------------------------------------------- addCurve
    def addCurve(self,iCurSubGraph,curve):
        # If colors are set by default, we use a color map
        # Rmk : this is useless in the workflow using editCurve class but usefull while reloading a file ...
        curve.correctColor(len(self.fig.subGraph[iCurSubGraph].curves))
        self.fig.subGraph[iCurSubGraph].curves.append(curve)
        self.updateGraph(iCurSubGraph)
    # ----------------------------------------------------- removeCurvesZoneName
    def removeCurvesZoneName(self,ax_name,zonename):
        self.fig.removeCurvesZoneName(ax_name,zonename)
    # ----------------------------------------------------- updateCurvesZoneName
    def updateCurvesZoneName(self,iCurSubGraph,oldZoneName,newZoneName):
        self.fig.updateCurvesZoneName(iCurSubGraph,oldZoneName,newZoneName)
    # -------------------------------------------------------------- updateGraph
    def updateGraph(self,iCurSubGraph):
        self.fig.drawOneFigure(iCurSubGraph)
        self.canvas.draw()
    # -------------------------------------------------------- applyViewSettings
    def applyViewSettings(self):
        if self.subPlotParams.isActive:

# Was aimed to turn off eventually the tight_layout, but seems useless ... and not working !
#            try:
#                self.fig.instance.tight_layout(False)
#                print('STOP TIGHT LAYOUT')
#            except AttributeError:
#                pass

            self.fig.instance.subplots_adjust(left=self.subPlotParams.left,
                                             right=self.subPlotParams.right,
                                             top=self.subPlotParams.top,
                                             bottom=self.subPlotParams.bottom,
                                             hspace=self.subPlotParams.hspace,
                                             wspace=self.subPlotParams.wspace)
        if self.tightLayout.isActive:
            self.fig.instance.tight_layout(pad=self.tightLayout.pad,
                                  h_pad=self.tightLayout.hpad,
                                  w_pad=self.tightLayout.wpad)
        self.canvas.draw()
    # ----------------------------------------------------------------- setValue
    def setValue(self,variable,value):
        if variable == 'image_background_color':
            self.image_background_color = value
        elif variable == 'image_background_alpha':
            self.image_background_alpha = value
        elif variable == 'subgraph_background_color':
            self.subgraph_background_color = value
        elif variable == 'subgraph_background_alpha':
            self.subgraph_background_alpha = value
# ==============================================================================
# ==============================================================================

class Graph():
    """
    An object of class Graph corresponds to a window where plots are drawn. A graph window can manage several plots.
    """
    # --------------------------------------------------------------------- init
    def __init__(self,parent,name,conf,dpi=None,figsize=None):

        self.parent = parent
        self.figsize = figsize
        self.dpi = dpi
        self.name=self.soleName(name)
        self.conf=conf
        self.initialize()
        #
        self.useSubPlotParams = True
        self.subPlotParams = SubPlotParams()
        self.tightLayout = TightLayout()
        self.image_background_color = default_values['Graph']['image_background_color']
        self.image_background_alpha = default_values['Graph']['image_background_alpha']
        self.subgraph_background_color = default_values['Graph']['subgraph_background_color']
        self.subgraph_background_alpha = default_values['Graph']['subgraph_background_alpha']
    # --------------------------------------------------------------- initialize
    def initialize(self):
# TODO        # Il faudra regler : le titre de la fenetre

        self.parent.updateGraphName2Id()
        self.parent.graphName2Id[self.name]=len(self.parent.graphWdwL)
        self.index = len(self.parent.graphWdwL)
        self.parent.graphWdwL.append(self)

        # Figure
        self.fig = MatplotlibFigure(self.parent,self,self.conf,self.dpi,self.figsize)
        self.setName(self.name)




        # Update Main window : graph name list
#        self.parent.updateactiveGraph()

        # Configure une zone graphique
    # ------------------------------------------------------ deleteZoneInCurve
    def deleteZoneInCurve(self,iCurSubGraph,zoneName):
        self.fig.deleteZoneInCurve(iCurSubGraph,zoneName)
    # ------------------------------------------------------ updateGroupCurvesZoneName
    def updateGroupCurvesZoneName(self,iCurSubGraph,oldZoneList,newZoneList):
        self.fig.updateGroupCurvesZoneName(iCurSubGraph,oldZoneList,newZoneList)
    # -------------------------------------------------------- applyViewSettings
    def applyViewSettings(self):
        if self.subPlotParams.isActive:

# Was aimed to turn off eventually the tight_layout, but seems useless ... and not working !
#            try:
#                self.fig.instance.tight_layout(False)
#                print('STOP TIGHT LAYOUT')
#            except AttributeError:
#                pass

            self.fig.instance.subplots_adjust(left=self.subPlotParams.left,
                                             right=self.subPlotParams.right,
                                             top=self.subPlotParams.top,
                                             bottom=self.subPlotParams.bottom,
                                             hspace=self.subPlotParams.hspace,
                                             wspace=self.subPlotParams.wspace)
        if self.tightLayout.isActive:
            self.fig.instance.tight_layout(pad=self.tightLayout.pad,
                                  h_pad=self.tightLayout.hpad,
                                  w_pad=self.tightLayout.wpad)
        self.drawFigure()
    # ------------------------------------------------------ updateSubPlotParams
    def updateSubPlotParams(self,params):
        for var in params:
            self.subPlotParams.setValue(var,params[var])
            if var == 'isActive':
                self.tightLayout.setValue(var,not params[var])
        # Update Figure
        self.applyViewSettings()
    # -------------------------------------------------------- updateTightLayout
    def updateTightLayout(self,params):
        for var in params:
            self.tightLayout.setValue(var,params[var])
            if var == 'isActive':
                self.subPlotParams.setValue(var,not params[var])
        # Update Figure
        self.applyViewSettings()
    # --------------------------------------------------------------------- save
    def save(self,path,format=None):
        self.fig.saveFigure(path,format=format)
    # ------------------------------------------------------------------ setName
    def setName(self,name):
        (self.getFig()).canvas.set_window_title(name)
    # --------------------------------------------------------------- drawFigure
    def drawFigure(self):
        self.fig.drawFigure()
    # --------------------------------------------------------------- showFigure
    def showFigure(self):
        self.fig.drawFigure()
        figure = self.getFig()
        figure.show()
        # Other Way of doing, probably better (TODO : test it !!!)
#        plt.figure(figure.number)
#        plt.show()
    # ------------------------------------------------------------------ addGrid
    def getFig(self):
        return self.fig.getFig()
    # ------------------------------------------------------------------ addGrid
    def getGrid(self,iCurSubGraph,ind=0,axis=None):
        if axis:
            ind = axis.getInd()
        #
        return self.fig.getGrid(iCurSubGraph,ind)
    # ------------------------------------------------------------------ addAxis
    def addAxis(self,iCurSubGraph,shared=None,ind=0,axis=None):
        if axis:
            ind = axis.getInd()
        #
        axis_to_twin = ind
        if shared is None:
            self.fig.subGraph[iCurSubGraph].addAxis()
            newind = len(self.fig.subGraph[iCurSubGraph].axis)-1
            newaxis = self.fig.getAxis(iCurSubGraph,ind=newind)
            newaxis.ind = newind
            return newaxis
        elif ((shared == 'x') or (shared == 'X')):
            self.fig.subGraph[iCurSubGraph].addAxisTwinX(axis_to_twin)
            newind = len(self.fig.subGraph[iCurSubGraph].axis)-1
            newaxis = self.fig.getAxis(iCurSubGraph,ind=newind)
            newaxis.ind = newind
            return newaxis
        elif ((shared == 'y') or (shared == 'Y')):
            self.fig.subGraph[iCurSubGraph].addAxisTwinY(axis_to_twin)
            newind = len(self.fig.subGraph[iCurSubGraph].axis)-1
            newaxis = self.fig.getAxis(iCurSubGraph,ind=newind)
            newaxis.ind = newind
            return newaxis
        else:
            print('''### Error: value used for 'shared' is unknown, must be in [None, 'x', 'X', 'y', 'Y'].''')
    # ---------------------------------------------------------------- getLegend
    def getLegend(self,iCurSubGraph):
        return self.fig.getLegend(iCurSubGraph)
    # ------------------------------------------------------------------ getAxis
    def getAxis(self,iCurSubGraph,ind=0):
        return self.fig.getAxis(iCurSubGraph,ind)
    # -------------------------------------------------------------------- write
    def write(self,ind):
        line = '''graph_%s=Graph(obj,'%s','%s',dpi=%s,figsize=%s)'''%(ind,self.name,self.conf,self.dpi,self.figsize)
        return line
    # ----------------------------------------------------------------- soleName
    def soleName(self,name):
        graphNameList = [g.name for g in self.parent.graphWdwL]
        if not name in graphNameList:
            return name
        else:
            ind=0
            while name+'.%s'%ind in graphNameList:
                ind += 1
            return name+'.%s'%ind

#    # ---------------------------------------------------------------- cmd_close
#    def cmd_close(self):
#        self.fig.closeMovie()
#        del self.fig
#        del self.parent.graphWdwL[self.index]
#        self.parent.renumberGraph()
#        self.parent.updateactiveGraph()
#        self.destroy()
    # ----------------------------------------------------------------- addCurve
    def addCurve(self,iCurSubGraph,curve):
        # Without interface, if no color is specified while creating a curve, a color map is used
        curve.correctColor(len(self.fig.subGraph[iCurSubGraph].curves))
        self.fig.subGraph[iCurSubGraph].curves.append(curve)
        self.updateGraph(iCurSubGraph)
    # -------------------------------------------------------------- updateGraph
    def updateGraph(self,iCurSubGraph):
        self.fig.drawOneFigure(iCurSubGraph)
    # ----------------------------------------------------------------- setValue
    def setValue(self,variable,value):
        if variable == 'image_background_color':
            self.image_background_color = value
        elif variable == 'image_background_alpha':
            self.image_background_alpha = value
        elif variable == 'subgraph_background_color':
            self.subgraph_background_color = value
        elif variable == 'subgraph_background_alpha':
            self.subgraph_background_alpha = value
# ==============================================================================
# ==============================================================================
class editCurvesWindow(TK.Toplevel):
    # --------------------------------------------------------------------- init
    def __init__(self):
        TK.Toplevel.__init__(self)
        self.title('Edit  curves')
        self.protocol("WM_DELETE_WINDOW",self.cmd_close)
        #
        self.list_dialog = None
        self.input_dialog = None
        #
        self.frameList = {}
    # --------------------------------------------------------------- initialize
    def initialize(self,parent):
        self.parent = parent
        #
        self.graph = self.parent.activeGraph.val

        self.zone  = self.parent.position.val
        try:
            self.subGraph = self.parent.graphWdwL[self.graph].fig.subGraph[self.zone]
        except TypeError: # Exit if no graph & zone available
            self.cmd_close()
            return
#        except IndexError: # Error while closing window ...
#            return
        #
        self.grid_columnconfigure(0,weight=1)
        self.grid_rowconfigure(0,weight=1)
        #
        self.labelView = [True,True,True,False,False,False]
        self.frame = None
        self.createFrame()
    # -------------------------------------------------------------- update
    def updateBar(self,parent):
        self.frame.grid_forget()
        try:
            self.subGraph = self.parent.graphWdwL[self.graph].fig.subGraph[self.zone]
        except TypeError: # Exit if no graph & zone available
            self.cmd_close()
            return
        try:
            self.frame = self.frameList[self.graph][self.zone]
            self.frame.grid(row=0,column=0,sticky='NSEW')
        except KeyError:
            self.frame=None
            self.createFrame()
            # self.initialize(self.parent)

    # -------------------------------------------------------------- updateData
    def updateData(self):

        # Obviously, nothing to do ...
        # before : self.reloadWindow() ...
        return

    # -------------------------------------------------------------- expandLblFrame
    def expandLblFrame(self,event):
        widget = event.widget
        self.labelView[widget.id]= not self.labelView[widget.id]
        widget.display = self.labelView[widget.id]
        if widget.display:
            widget.frame.grid(row=0,column=0,sticky='NSEW')
            title = widget.title + "v "
        else:
            widget.frame.grid_forget()
            title = widget.title + "> "
        widget.config(text=title)
        self.geometry("")
    # -------------------------------------------------------------- createFrame
    def createFrame(self):
        self.geometry("")
        self.grid_columnconfigure(0,weight=1)
        self.grid_rowconfigure(0,weight=1)
        # Load colormap
        cm = plt.get_cmap(COLOR_MAP)
        # Close all opened dialog windows
        self.closeAllDialog()
        #
        self.frame = TTK.Frame(self)

        self.frame.grid(row=0,column=0,sticky='NESW')
        #
        if self.labelView[0]:
            self.frame.grid_columnconfigure(0,weight=4)
        else:
            self.frame.grid_columnconfigure(0,weight=1)
        if self.labelView[1]:
            self.frame.grid_columnconfigure(1,weight=3)
        else:
            self.frame.grid_columnconfigure(1,weight=1)
        if self.labelView[2]:
            self.frame.grid_columnconfigure(2,weight=3)
        else:
            self.frame.grid_columnconfigure(2,weight=1)
        if self.labelView[3]:
            self.frame.grid_columnconfigure(3,weight=5)
        else:
            self.frame.grid_columnconfigure(3,weight=1)
        if self.labelView[4]:
            self.frame.grid_columnconfigure(4,weight=3)
        else:
            self.frame.grid_columnconfigure(4,weight=1)
        if self.labelView[5]:
            self.frame.grid_columnconfigure(5,weight=2)
        else:
            self.frame.grid_columnconfigure(5,weight=1)
        # self.frame.grid_columnconfigure(6,weight=2)
        # self.frame.grid_columnconfigure(7,weight=1)
        # self.frame.grid_columnconfigure(8,weight=1)
        # self.frame.grid_columnconfigure(9,weight=1)
        # self.frame.grid_columnconfigure(10,weight=1)
        # self.frame.grid_columnconfigure(11,weight=1)
        # self.frame.grid_columnconfigure(12,weight=1)
        # self.frame.grid_columnconfigure(13,weight=1)
        # self.frame.grid_columnconfigure(14,weight=1)
        # self.frame.grid_columnconfigure(15,weight=1)
        # self.frame.grid_columnconfigure(16,weight=1)
        # self.frame.grid_columnconfigure(17,weight=1)
        # self.frame.grid_columnconfigure(18,weight=1)
        # self.frame.grid_columnconfigure(19,weight=1)
        #
        self.frame.grid_rowconfigure(0,weight=1)
        self.frame.grid_rowconfigure(1,weight=0)
        #
        lblframelvl1=[]
        #
        ########################################################################
        #################### -> Level 1
        lblframeManip = TTK.LabelFrame(self.frame, text = "Manip.")
        lblframeManip.grid(row=0,column=0,sticky='NSEW')
        lblframeManip.grid_columnconfigure(0,weight=1)
        lblframeManip.grid_rowconfigure(0,weight=1)
        lblframeManip.id=0
        lblframeManip.display =self.labelView[lblframeManip.id]
        #
        lblframeManip.title = "Manip. "
        if lblframeManip.display:
            title = lblframeManip.title + "v "
        else:
            title = lblframeManip.title + "> "
        lblframeManip.config(text=title)
        #
        lblframeManip.bind("<Button-1>",self.expandLblFrame)
        #
        frameManip = TTK.Frame(lblframeManip)
        frameManip.grid(row=0,column=0,sticky='NSEW')
        if not lblframeManip.display:
            frameManip.grid_forget()
        frameManip.grid_columnconfigure(0,weight=1)
        frameManip.grid_columnconfigure(1,weight=1)
        frameManip.grid_columnconfigure(2,weight=1)
        frameManip.grid_columnconfigure(3,weight=1)
        frameManip.grid_rowconfigure(0,weight=1)
        #
        lblframeManip.frame = frameManip
        # Create a hidden label to set minimum size of the lblframe respecting to its title
        label = TTK.Label(self.frame,text="Manip.")
        label.grid(row=0,column=0,in_=lblframeManip)
        label.lower(lblframeManip)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Selection
        lblframe = TTK.LabelFrame(frameManip, text="Selected")
        lblframe.grid(row=0,column=0,sticky='NESW')
        #
        lblframe.grid_columnconfigure(0,weight=1)
        for ind in range(len(self.subGraph.curves)+1):
            lblframe.grid_rowconfigure(ind,weight=1)
        #
        lblframelvl1.append(lblframe)
        #
        self.frame.selectionItem=[]
        #
        for ind in range(len(self.subGraph.curves)):
            var = TK.IntVar()
            CB = TTK.Checkbutton(lblframe,variable=var,command=lambda n=ind: self.cb_selection(n))#, variable=var)
            CB.val = var
            CB.ind=ind
            CB.var='selection'
            CB.grid(row=ind,column=0,sticky='NSEW')
            CB.container = self.frame.selectionItem
            self.frame.selectionItem.append(CB)
        # Curve to add
        ind = len(self.subGraph.curves)
        var = TK.IntVar()
        CB = TTK.Checkbutton(lblframe,variable=var,command=lambda n=ind: self.cb_selection(n),state=TK.DISABLED)#, variable=var)
        CB.val = var
        CB.ind=ind
        CB.var='selection'
        CB.grid(row=ind,column=0,sticky='NSEW')
        CB.container = self.frame.selectionItem
        self.frame.selectionItem.append(CB)

        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Id
        lblframe = TTK.LabelFrame(frameManip, text="Curve Id")
        lblframe.grid(row=0,column=1,sticky='NESW')
        #
        lblframe.grid_columnconfigure(0,weight=1)
        for ind in range(len(self.subGraph.curves)+1):
            lblframe.grid_rowconfigure(ind,weight=1)
        #
        lblframelvl1.append(lblframe)
        #
        self.frame.IdItem = []
        #
        for ind in range(len(self.subGraph.curves)):
            LBL = TK.Label(lblframe,text='%s'%ind)
            LBL.ind = ind
            LBL.grid(row=ind,column=0,sticky='NSEW')
            LBL.container = self.frame.IdItem
            self.frame.IdItem.append(LBL)
        # Curve to add
        ind = len(self.subGraph.curves)
        LBL = TK.Label(lblframe,text='to add')
        LBL.ind = ind
        LBL.grid(row=ind,column=0,sticky='NSEW')
        LBL.container = self.frame.IdItem
        self.frame.IdItem.append(LBL)

        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Axis to plot
        lblframe = TTK.LabelFrame(frameManip, text="Axis")
        lblframe.grid(row=0,column=2,sticky='NESW')
        #
        lblframe.grid_columnconfigure(0,weight=1)
        for ind in range(len(self.subGraph.curves)+1):
            lblframe.grid_rowconfigure(ind,weight=1)
        #
        lblframelvl1.append(lblframe)
        #
        self.frame.axisItem = []
        for ind in range(len(self.subGraph.curves)):
            c = self.subGraph.curves[ind]
            B = TTK.Button(lblframe,text=c.axis,command=lambda n=(ind,self.frame.axisItem): self.bt_click(n))
            B.list = [i for i in range(len(self.subGraph.axis))]
            B.val = c.axis
            B.var = 'ind_axis'
            B.ind = ind
            B.treatmentId = 0
            B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
            B.container = self.frame.axisItem
            self.frame.axisItem.append(B)
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=0,command=lambda n=(ind,self.frame.axisItem): self.bt_click(n))
        B.list = [i for i in range(len(self.subGraph.axis))]
        indexNone = 0
        B.val = 0
        B.var = 'axis'
        B.ind = ind
        B.treatmentId = 0
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.axisItem
        self.frame.axisItem.append(B)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Visibility
        lblframe = TTK.LabelFrame(frameManip, text="Visibility")
        lblframe.grid(row=0,column=3,sticky='NESW')
        #
        lblframe.grid_columnconfigure(0,weight=1)
        for ind in range(len(self.subGraph.curves)+1):
            lblframe.grid_rowconfigure(ind,weight=1)
        #
        lblframelvl1.append(lblframe)
        #
        self.frame.visibilityItem=[]
        #
        for ind in range(len(self.subGraph.curves)):
            var = TK.IntVar()
            var.set(1)
            # CB = TTK.Checkbutton(lblframe,variable=var,command=lambda n=ind: self.cb_visibility(n))#, variable=var)
            CB = TTK.Checkbutton(lblframe,variable=var,command=lambda n=ind: self.cb_visibility(n))#, variable=var)
            CB.val = var
            CB.var = 'visible'
            CB.ind=ind
            CB.grid(row=ind,column=0,sticky='NSEW')
            CB.container = self.frame.visibilityItem
            self.frame.visibilityItem.append(CB)
        # Curve to add
        ind = len(self.subGraph.curves)
        var = TK.IntVar()
        var.set(1)
        # CB = TTK.Checkbutton(lblframe,variable=var,state=TK.DISABLED,command=lambda n=ind: self.cb_visibility(n))#, variable=var)
        CB = TTK.Checkbutton(lblframe,variable=var,state=TK.DISABLED,command=lambda n=ind: self.cb_visibility(n))#, variable=var)
        CB.val = var
        CB.var = 'visible'
        CB.ind=ind
        CB.grid(row=ind,column=0,sticky='NSEW')
        CB.container = self.frame.visibilityItem
        self.frame.visibilityItem.append(CB)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        lblframeData = TTK.LabelFrame(self.frame, text = "Data")
        lblframeData.grid(row=0,column=1,sticky='NSEW')
        lblframeData.grid_columnconfigure(0,weight=1)
        lblframeData.grid_rowconfigure(0,weight=1)
        lblframeData.id = 1
        lblframeData.display = self.labelView[lblframeData.id]
        #
        lblframeData.title = "Data "
        if lblframeData.display:
            title = lblframeData.title + "v "
        else:
            title = lblframeData.title + "> "
        lblframeData.config(text=title)
        #
        lblframeData.bind("<Button-1>",self.expandLblFrame)
        #
        frameData = TTK.Frame(lblframeData)
        frameData.grid(row=0,column=0,sticky='NSEW')
        if not lblframeData.display:
            frameData.grid_forget()
        frameData.grid_columnconfigure(0,weight=1)
        frameData.grid_columnconfigure(1,weight=1)
        frameData.grid_columnconfigure(2,weight=1)
        frameData.grid_rowconfigure(0,weight=1)
        #
        lblframeData.frame = frameData
        # Create a hidden label to set minimum size of the lblframe respecting to its title
        label = TTK.Label(self.frame,text="Data")
        label.grid(row=0,column=0,in_=lblframeData)
        label.lower(lblframeData)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Zone
        lblframe = TTK.LabelFrame(frameData, text="Zone(s)")
        lblframe.grid(row=0,column=0,sticky='NESW')
        #
        lblframe.grid_columnconfigure(0,weight=1)
        for ind in range(len(self.subGraph.curves)+1):
            lblframe.grid_rowconfigure(ind,weight=1)
        #
        lblframelvl1.append(lblframe)
        #
        self.frame.zoneItem = []
        for ind in range(len(self.subGraph.curves)):
            c = self.subGraph.curves[ind]
            B = TTK.Button(lblframe,text=len(c.zone),command=lambda n=(ind,self.frame.zoneItem): self.bt_click(n))
            B.val = c.zone
#            B.val = c.zoneList
            B.var = 'zone'
            B.ind = ind
            B.treatmentId = 4 # 4 is for selectzones
            B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
            B.container = self.frame.zoneItem
            self.frame.zoneItem.append(B)
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=len(self.parent.data.keys()),command=lambda n=(ind,self.frame.zoneItem): self.bt_click(n))
        B.val = self.parent.data.keys()
#        B.listUnused = []
#        B.val = B.list
        B.var = 'zone'
        B.ind = ind
        B.treatmentId = 4
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.zoneItem
        self.frame.zoneItem.append(B)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& VarX
        lblframe = TTK.LabelFrame(frameData, text="VarX")
        lblframe.grid(row=0,column=1,sticky='NESW')
        #
        lblframe.grid_columnconfigure(0,weight=1)
        for ind in range(len(self.subGraph.curves)+1):
            lblframe.grid_rowconfigure(ind,weight=1)
        #
        lblframelvl1.append(lblframe)
        #
        self.frame.varxItem = []
        for ind in range(len(self.subGraph.curves)):
            c = self.subGraph.curves[ind]
            B = TTK.Button(lblframe,text=c.varx,command=lambda n=(ind,self.frame.varxItem): self.bt_click(n))
            B.ind = ind
            B.list = self.filterVarWithZone(B)
            B.val = c.varx
            B.var = 'varx'
            B.treatmentId = 0
            B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
            B.container = self.frame.varxItem
            self.frame.varxItem.append(B)
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=None,command=lambda n=(ind,self.frame.varxItem): self.bt_click(n))
        B.ind = ind
        B.list = self.filterVarWithZone(B)
        if len(B.list)!=0:
            B.config(text=B.list[0])
            B.val = B.list[0]
        else:
            B.config(text='')
            B.val = None
        B.var = 'varx'
        B.treatmentId = 0
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.varxItem
        self.frame.varxItem.append(B)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& VarY
        lblframe = TTK.LabelFrame(frameData, text="VarY")
        lblframe.grid(row=0,column=2,sticky='NESW')
        #
        lblframe.grid_columnconfigure(0,weight=1)
        for ind in range(len(self.subGraph.curves)+1):
            lblframe.grid_rowconfigure(ind,weight=1)
        #
        lblframelvl1.append(lblframe)
        #
        self.frame.varyItem = []
        for ind in range(len(self.subGraph.curves)):
            c = self.subGraph.curves[ind]
            B = TTK.Button(lblframe,text=c.vary,command=lambda n=(ind,self.frame.varyItem): self.bt_click(n))
            B.ind = ind
            B.list = self.filterVarWithZone(B)
            B.val = c.vary
            B.var = 'vary'
            B.ind = ind
            B.treatmentId = 0
            B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
            B.container = self.frame.varyItem
            self.frame.varyItem.append(B)
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=None,command=lambda n=(ind,self.frame.varyItem): self.bt_click(n))
        B.ind = ind
        B.list = self.filterVarWithZone(B)
        if len(B.list)!=0:
            B.config(text=B.list[0])
            B.val = B.list[0]
        else:
            B.config(text='')
            B.val = None
        B.var = 'vary'
        B.treatmentId = 0
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.varyItem
        self.frame.varyItem.append(B)


        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        lblframeLine = TTK.LabelFrame(self.frame, text="Line")
        lblframeLine.grid(row=0,column=2,sticky='NSEW')
        lblframeLine.grid_columnconfigure(0,weight=1)
        lblframeLine.grid_rowconfigure(0,weight=1)
        lblframeLine.id = 2
        lblframeLine.display = self.labelView[lblframeLine.id]
        #
        lblframeLine.title = "Line "
        if lblframeLine.display:
            title = lblframeLine.title + "v "
        else:
            title = lblframeLine.title + "> "
        lblframeLine.config(text=title)
        #
        lblframeLine.bind("<Button-1>",self.expandLblFrame)
        #
        frameLine = TTK.Frame(lblframeLine)
        frameLine.grid(row=0,column=0,sticky='NSEW')
        if not lblframeLine.display:
            frameLine.grid_forget()
        frameLine.grid_columnconfigure(0,weight=1)
        frameLine.grid_columnconfigure(1,weight=1)
        frameLine.grid_columnconfigure(2,weight=1)
        frameLine.grid_rowconfigure(0,weight=1)
        #
        lblframeLine.frame = frameLine
        #
        # Create a hidden label to set minimum size of the lblframe respecting to its title
        label = TTK.Label(self.frame,text="Line")
        label.grid(row=0,column=0,in_=lblframeLine)
        label.lower(lblframeLine)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Line color
        lblframe = TTK.LabelFrame(frameLine, text="L.color")
        lblframe.grid(row=0,column=0,sticky='NESW')
        #
        lblframe.grid_columnconfigure(0,weight=1)
        for ind in range(len(self.subGraph.curves)+1):
            lblframe.grid_rowconfigure(ind,weight=1)
        #
        lblframelvl1.append(lblframe)
        #
        self.frame.line_colorItem = []
        for ind in range(len(self.subGraph.curves)):
            c = self.subGraph.curves[ind]
            B = TK.Button(lblframe,command=lambda n=(ind,self.frame.line_colorItem): self.bt_click(n))
            B.list = []
            B.val = c.line_color
            B.var = 'line_color'
            B.config(bg=B.val)
            B.config(activebackground=B.val)
            B.ind = ind
            B.treatmentId = 1
            B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
            B.container = self.frame.line_colorItem
            self.frame.line_colorItem.append(B)
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TK.Button(lblframe,command=lambda n=(ind,self.frame.line_colorItem): self.bt_click(n))
        B.list = []
        color = cm(1.*(ind%NUM_COLORS)/NUM_COLORS)
        B.val = '#%02x%02x%02x' % (int(255*color[0]),int(255*color[1]),int(255*color[2]))
        B.var = 'line_color'
        B.config(bg=B.val)
        B.config(activebackground=B.val)
        B.ind = ind
        B.treatmentId = 1
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.line_colorItem
        self.frame.line_colorItem.append(B)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Line size
        lblframe = TTK.LabelFrame(frameLine, text="L.size")
        lblframe.grid(row=0,column=1,sticky='NESW')
        #
        lblframe.grid_columnconfigure(0,weight=1)
        for ind in range(len(self.subGraph.curves)+1):
            lblframe.grid_rowconfigure(ind,weight=1)
        #
        lblframelvl1.append(lblframe)
        #
        self.frame.line_widthItem = []
        for ind in range(len(self.subGraph.curves)):
            c = self.subGraph.curves[ind]
            B = TTK.Button(lblframe,text=c.line_width,command=lambda n=(ind,self.frame.line_widthItem): self.bt_click(n))
            B.list = (np.arange(0.5,100,0.5)).tolist()
            B.val = c.line_width
            B.var = 'line_width'
            B.ind = ind
            B.treatmentId = 0
            B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
            B.container = self.frame.line_widthItem
            self.frame.line_widthItem.append(B)
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=default_values['Curve']['line_width'],command=lambda n=(ind,self.frame.line_widthItem): self.bt_click(n))
        B.list = (np.arange(0.5,100,0.5)).tolist()
        B.val = default_values['Curve']['line_width']
        B.var = 'line_width'
        B.ind = ind
        B.treatmentId = 0
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.line_widthItem
        self.frame.line_widthItem.append(B)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Line style
        lblframe = TTK.LabelFrame(frameLine, text="L.style")
        lblframe.grid(row=0,column=2,sticky='NESW')
        #
        lblframe.grid_columnconfigure(0,weight=1)
        for ind in range(len(self.subGraph.curves)+1):
            lblframe.grid_rowconfigure(ind,weight=1)
        #
        lblframelvl1.append(lblframe)
        #
        self.frame.line_styleItem = []
        for ind in range(len(self.subGraph.curves)):
            c = self.subGraph.curves[ind]
            B = TTK.Button(lblframe,text=c.line_style,command=lambda n=(ind,self.frame.line_styleItem): self.bt_click(n))
            B.list = linestylelist
            B.val = c.line_style
            B.var = 'line_style'
            B.ind = ind
            B.treatmentId = 0
            B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
            B.container = self.frame.line_styleItem
            self.frame.line_styleItem.append(B)
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=default_values['Curve']['line_style'],command=lambda n=(ind,self.frame.line_styleItem): self.bt_click(n))
        B.list = linestylelist
        B.val = default_values['Curve']['line_style']
        B.var = 'line_style'
        B.ind = ind
        B.treatmentId = 0
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.line_styleItem
        self.frame.line_styleItem.append(B)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        lblframeSymbol = TTK.LabelFrame(self.frame, text = "Symbol")
        lblframeSymbol.grid(row=0,column=3,sticky='NSEW')
        lblframeSymbol.grid_columnconfigure(0,weight=1)
        lblframeSymbol.grid_rowconfigure(0,weight=1)
        lblframeSymbol.id = 3
        lblframeSymbol.display = self.labelView[lblframeSymbol.id]
        #
        lblframeSymbol.title = "Symbol "
        if lblframeSymbol.display:
            title = lblframeSymbol.title + "v "
        else:
            title = lblframeSymbol.title + "> "
        lblframeSymbol.config(text=title)
        #
        lblframeSymbol.bind("<Button-1>",self.expandLblFrame)
        #
        frameSymbol = TTK.Frame(lblframeSymbol)
        frameSymbol.grid(row=0,column=0,sticky='NSEW')
        if not lblframeSymbol.display:
            frameSymbol.grid_forget()
        frameSymbol.grid_columnconfigure(0,weight=1)
        frameSymbol.grid_columnconfigure(1,weight=1)
        frameSymbol.grid_columnconfigure(2,weight=1)
        frameSymbol.grid_columnconfigure(3,weight=1)
        frameSymbol.grid_columnconfigure(4,weight=1)
        frameSymbol.grid_rowconfigure(0,weight=1)
        #
        lblframeSymbol.frame = frameSymbol
        #
        # Create a hidden label to set minimum size of the lblframe respecting to its title
        label = TTK.Label(self.frame,text="Symbol")
        label.grid(row=0,column=0,in_=lblframeSymbol)
        label.lower(lblframeSymbol)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Symbol Face color
        lblframe = TTK.LabelFrame(frameSymbol, text="S.color")
        lblframe.grid(row=0,column=0,sticky='NESW')
        #
        lblframe.grid_columnconfigure(0,weight=1)
        for ind in range(len(self.subGraph.curves)+1):
            lblframe.grid_rowconfigure(ind,weight=1)
        #
        lblframelvl1.append(lblframe)
        #
        self.frame.marker_face_colorItem = []
        for ind in range(len(self.subGraph.curves)):
            c = self.subGraph.curves[ind]
            B = TK.Button(lblframe,command=lambda n=(ind,self.frame.marker_face_colorItem): self.bt_click(n))
            B.list = []
            B.val = c.marker_face_color
            B.var = 'marker_face_color'
            B.config(bg=B.val)
            B.config(activebackground=B.val)
            B.ind = ind
            B.treatmentId = 1
            B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
            B.container = self.frame.marker_face_colorItem
            self.frame.marker_face_colorItem.append(B)
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TK.Button(lblframe,command=lambda n=(ind,self.frame.marker_face_colorItem): self.bt_click(n))
        B.list = []
        color = cm(1.*(ind%NUM_COLORS)/NUM_COLORS)
        B.val = '#%02x%02x%02x' % (int(255*color[0]),int(255*color[1]),int(255*color[2]))
        B.var = 'marker_face_color'
        B.config(bg=B.val)
        B.config(activebackground=B.val)
        B.ind = ind
        B.treatmentId = 1
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.marker_face_colorItem
        self.frame.marker_face_colorItem.append(B)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Symbol Face size
        lblframe = TTK.LabelFrame(frameSymbol, text="S.size")
        lblframe.grid(row=0,column=1,sticky='NESW')
        #
        lblframe.grid_columnconfigure(0,weight=1)
        for ind in range(len(self.subGraph.curves)+1):
            lblframe.grid_rowconfigure(ind,weight=1)
        #
        lblframelvl1.append(lblframe)
        #
        self.frame.marker_sizeItem = []
        for ind in range(len(self.subGraph.curves)):
            c = self.subGraph.curves[ind]
            B = TTK.Button(lblframe,text=c.marker_size,command=lambda n=(ind,self.frame.marker_sizeItem): self.bt_click(n))
            B.list = (np.arange(0.5,100,0.5)).tolist()
            B.val = c.marker_size
            B.var = 'marker_size'
            B.ind = ind
            B.treatmentId = 0
            B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
            B.container = self.frame.marker_sizeItem
            self.frame.marker_sizeItem.append(B)
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=default_values['Curve']['marker_size'],command=lambda n=(ind,self.frame.marker_sizeItem): self.bt_click(n))
        B.list = (np.arange(0.5,100,0.5)).tolist()
        B.val = default_values['Curve']['marker_size']
        B.var = 'marker_size'
        B.ind = ind
        B.treatmentId = 0
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.marker_sizeItem
        self.frame.marker_sizeItem.append(B)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Symbol Type
        lblframe = TTK.LabelFrame(frameSymbol, text="S.type")
        lblframe.grid(row=0,column=2,sticky='NESW')
        #
        lblframe.grid_columnconfigure(0,weight=1)
        for ind in range(len(self.subGraph.curves)+1):
            lblframe.grid_rowconfigure(ind,weight=1)
        #
        lblframelvl1.append(lblframe)
        #
        self.frame.marker_styleItem = []
        for ind in range(len(self.subGraph.curves)):
            c = self.subGraph.curves[ind]
            B = TTK.Button(lblframe,text=c.marker_style,command=lambda n=(ind,self.frame.marker_styleItem): self.bt_click(n))
            B.list = markername
            B.val = c.marker_style
            B.var = 'marker_style'
            B.ind = ind
            B.treatmentId = 0
            B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
            B.container = self.frame.marker_styleItem
            self.frame.marker_styleItem.append(B)
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=default_values['Curve']['marker_style'],command=lambda n=(ind,self.frame.marker_styleItem): self.bt_click(n))
        B.list = markername
        B.val = default_values['Curve']['marker_style']
        indexNone = markername.index(B.val)
        B.var = 'marker_style'
        B.ind = ind
        B.treatmentId = 0
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.marker_styleItem
        self.frame.marker_styleItem.append(B)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Symbol edge color
        lblframe = TTK.LabelFrame(frameSymbol, text="S.edge color")
        lblframe.grid(row=0,column=3,sticky='NESW')
        #
        lblframe.grid_columnconfigure(0,weight=1)
        for ind in range(len(self.subGraph.curves)+1):
            lblframe.grid_rowconfigure(ind,weight=1)
        #
        lblframelvl1.append(lblframe)
        #
        self.frame.marker_edge_colorItem = []
        for ind in range(len(self.subGraph.curves)):
            c = self.subGraph.curves[ind]
            B = TK.Button(lblframe,command=lambda n=(ind,self.frame.marker_edge_colorItem): self.bt_click(n))
            B.list = []
            B.val = c.marker_edge_color
            B.var = 'marker_edge_color'
            B.config(bg=B.val)
            B.config(activebackground=B.val)
            B.ind = ind
            B.treatmentId = 1
            B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
            B.container = self.frame.marker_edge_colorItem
            self.frame.marker_edge_colorItem.append(B)
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TK.Button(lblframe,command=lambda n=(ind,self.frame.marker_edge_colorItem): self.bt_click(n))
        B.list = []
        color = cm(1.*(ind%NUM_COLORS)/NUM_COLORS)
        B.val = '#%02x%02x%02x' % (int(255*color[0]),int(255*color[1]),int(255*color[2]))
        B.var = 'marker_edge_color'
        B.config(bg=B.val)
        B.config(activebackground=B.val)
        B.ind = ind
        B.treatmentId = 1
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.marker_edge_colorItem
        self.frame.marker_edge_colorItem.append(B)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Symbol edge size
        lblframe = TTK.LabelFrame(frameSymbol, text="S.edge size")
        lblframe.grid(row=0,column=4,sticky='NESW')
        #
        lblframe.grid_columnconfigure(0,weight=1)
        for ind in range(len(self.subGraph.curves)+1):
            lblframe.grid_rowconfigure(ind,weight=1)
        #
        lblframelvl1.append(lblframe)
        #
        self.frame.marker_edge_widthItem = []
        for ind in range(len(self.subGraph.curves)):
            c = self.subGraph.curves[ind]
            B = TTK.Button(lblframe,text=c.marker_edge_width,command=lambda n=(ind,self.frame.marker_edge_widthItem): self.bt_click(n))
            B.list = (np.arange(0.5,100,0.5)).tolist()
            B.val = c.marker_edge_width
            B.var = 'marker_edge_width'
            B.ind = ind
            B.treatmentId = 0
            B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
            B.container = self.frame.marker_edge_widthItem
            self.frame.marker_edge_widthItem.append(B)
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=default_values['Curve']['marker_edge_width'],command=lambda n=(ind,self.frame.marker_edge_widthItem): self.bt_click(n))
        B.list = (np.arange(0.5,100,0.5)).tolist()
        B.val = default_values['Curve']['marker_edge_width']
        B.var = 'marker_edge_width'
        B.ind = ind
        B.treatmentId = 0
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.marker_edge_widthItem
        self.frame.marker_edge_widthItem.append(B)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        lblframeSymbolSampling = TTK.LabelFrame(self.frame, text = "Symb. Sampling")
        lblframeSymbolSampling.grid(row=0,column=4,sticky='NSEW')
        lblframeSymbolSampling.grid_columnconfigure(0,weight=1)
        lblframeSymbolSampling.grid_rowconfigure(0,weight=1)
        lblframeSymbolSampling.id = 4
        lblframeSymbolSampling.display = self.labelView[lblframeSymbolSampling.id]
        #
        lblframeSymbolSampling.title = "Symb. Sampling "
        if lblframeSymbolSampling.display:
            title = lblframeSymbolSampling.title + "v "
        else:
            title = lblframeSymbolSampling.title + "> "
        lblframeSymbolSampling.config(text=title)
        #
        lblframeSymbolSampling.bind("<Button-1>",self.expandLblFrame)
        #
        frameSymbolSampling = TTK.Frame(lblframeSymbolSampling)
        frameSymbolSampling.grid(row=0,column=0,sticky='NSEW')
        if not lblframeSymbolSampling.display:
            frameSymbolSampling.grid_forget()
        frameSymbolSampling.grid_columnconfigure(0,weight=1)
        frameSymbolSampling.grid_columnconfigure(1,weight=1)
        frameSymbolSampling.grid_columnconfigure(2,weight=1)
        frameSymbolSampling.grid_rowconfigure(0,weight=1)
        #
        lblframeSymbolSampling.frame = frameSymbolSampling
        #
        # Create a hidden label to set minimum size of the lblframe respecting to its title
        label = TTK.Label(self.frame,text="Symb. Sampling")
        label.grid(row=0,column=0,in_=lblframeSymbolSampling)
        label.lower(lblframeSymbolSampling)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Symbol sampling start
        lblframe = TTK.LabelFrame(frameSymbolSampling, text="Start")
        lblframe.grid(row=0,column=0,sticky='NESW')
        #
        lblframe.grid_columnconfigure(0,weight=1)
        for ind in range(len(self.subGraph.curves)+1):
            lblframe.grid_rowconfigure(ind,weight=1)
        #
        lblframelvl1.append(lblframe)
        #
        self.frame.marker_sampling_startItem = []
        for ind in range(len(self.subGraph.curves)):
            c = self.subGraph.curves[ind]
            B = TTK.Button(lblframe,text=c.marker_sampling_start,command=lambda n=(ind,self.frame.marker_sampling_startItem): self.bt_click(n))
            B.list = []
            B.val = c.marker_sampling_start
            B.var = 'marker_sampling_start'
            B.ind = ind
            B.treatmentId = 2
            B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
            B.container = self.frame.marker_sampling_startItem
            self.frame.marker_sampling_startItem.append(B)
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=default_values['Curve']['marker_sampling_start'],command=lambda n=(ind,self.frame.marker_sampling_startItem): self.bt_click(n))
        B.list = []
        B.val = default_values['Curve']['marker_sampling_start']
        B.var = 'marker_sampling_start'
        B.ind = ind
        B.treatmentId = 2
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.marker_sampling_startItem
        self.frame.marker_sampling_startItem.append(B)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Symbol sampling end
        lblframe = TTK.LabelFrame(frameSymbolSampling, text="End")
        lblframe.grid(row=0,column=1,sticky='NESW')
        #
        lblframe.grid_columnconfigure(0,weight=1)
        for ind in range(len(self.subGraph.curves)+1):
            lblframe.grid_rowconfigure(ind,weight=1)
        #
        lblframelvl1.append(lblframe)
        #
        self.frame.marker_sampling_endItem = []
        for ind in range(len(self.subGraph.curves)):
            c = self.subGraph.curves[ind]
            B = TTK.Button(lblframe,text=c.marker_sampling_end,command=lambda n=(ind,self.frame.marker_sampling_endItem): self.bt_click(n))
            B.list = []
            B.val = c.marker_sampling_end
            B.var = 'marker_sampling_end'
            B.ind = ind
            B.treatmentId = 2
            B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
            B.container = self.frame.marker_sampling_endItem
            self.frame.marker_sampling_endItem.append(B)
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=default_values['Curve']['marker_sampling_end'],command=lambda n=(ind,self.frame.marker_sampling_endItem): self.bt_click(n))
        B.list = []
        B.val = default_values['Curve']['marker_sampling_end']
        B.var = 'marker_sampling_end'
        B.ind = ind
        B.treatmentId = 2
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.marker_sampling_endItem
        self.frame.marker_sampling_endItem.append(B)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Symbol sampling step
        lblframe = TTK.LabelFrame(frameSymbolSampling, text="Step")
        lblframe.grid(row=0,column=2,sticky='NESW')
        #
        lblframe.grid_columnconfigure(0,weight=1)
        for ind in range(len(self.subGraph.curves)+1):
            lblframe.grid_rowconfigure(ind,weight=1)
        #
        lblframelvl1.append(lblframe)
        #
        self.frame.marker_sampling_stepItem = []
        for ind in range(len(self.subGraph.curves)):
            c = self.subGraph.curves[ind]
            B = TTK.Button(lblframe,text=c.marker_sampling_step,command=lambda n=(ind,self.frame.marker_sampling_stepItem): self.bt_click(n))
            B.list = []
            B.val = c.marker_sampling_step
            B.var = 'marker_sampling_step'
            B.ind = ind
            B.treatmentId = 2
            B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
            B.container = self.frame.marker_sampling_stepItem
            self.frame.marker_sampling_stepItem.append(B)
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=default_values['Curve']['marker_sampling_step'],command=lambda n=(ind,self.frame.marker_sampling_stepItem): self.bt_click(n))
        B.list = []
        B.val = default_values['Curve']['marker_sampling_step']
        B.var = 'marker_sampling_step'
        B.ind = ind
        B.treatmentId = 2
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.marker_sampling_stepItem
        self.frame.marker_sampling_stepItem.append(B)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        lblframeLegend = TTK.LabelFrame(self.frame, text="Legend")
        lblframeLegend.grid(row=0,column=5,sticky='NSEW')
        lblframeLegend.grid_columnconfigure(0,weight=1)
        lblframeLegend.grid_rowconfigure(0,weight=1)
        lblframeLegend.id = 5
        lblframeLegend.display = self.labelView[lblframeLegend.id]
        #
        lblframeLegend.title = "Legend. "
        if lblframeLegend.display:
            title = lblframeLegend.title + "v "
        else:
            title = lblframeLegend.title + "> "
        lblframeLegend.config(text=title)
        #
        lblframeLegend.bind("<Button-1>",self.expandLblFrame)
        #
        frameLegend = TTK.Frame(lblframeLegend)
        frameLegend.grid(row=0,column=0,sticky='NSEW')
        if not lblframeLegend.display:
            frameLegend.grid_forget()
        frameLegend.grid_columnconfigure(0,weight=1)
        frameLegend.grid_columnconfigure(1,weight=1)
        frameLegend.grid_rowconfigure(0,weight=1)
        #
        lblframeLegend.frame = frameLegend
        #
        # Create a hidden label to set minimum size of the lblframe respecting to its title
        label = TTK.Label(self.frame,text="Legend")
        label.grid(row=0,column=0,in_=lblframeLegend)
        label.lower(lblframeLegend)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Legend label
        lblframe = TTK.LabelFrame(frameLegend, text="Legend label")
        lblframe.grid(row=0,column=0,sticky='NESW')
        #
        lblframe.grid_columnconfigure(0,weight=1)
        for ind in range(len(self.subGraph.curves)+1):
            lblframe.grid_rowconfigure(ind,weight=1)
        #
        lblframelvl1.append(lblframe)
        #
        self.frame.legend_labelItem = []
        for ind in range(len(self.subGraph.curves)):
            c = self.subGraph.curves[ind]
            B = TTK.Button(lblframe,text=c.legend_label,command=lambda n=(ind,self.frame.legend_labelItem): self.bt_click(n))
            B.list = []
            B.val = c.legend_label
            B.var = 'legend_label'
            B.ind = ind
            B.treatmentId = 3
            B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
            B.container = self.frame.legend_labelItem
            self.frame.legend_labelItem.append(B)
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text="...",command=lambda n=(ind,self.frame.legend_labelItem): self.bt_click(n))
        B.list = []
        B.val = "..."
        B.var = 'legend_label'
        B.ind = ind
        B.treatmentId = 3
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.legend_labelItem
        self.frame.legend_labelItem.append(B)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Legend display
        lblframe = TTK.LabelFrame(frameLegend, text="Legend display")
        lblframe.grid(row=0,column=1,sticky='NESW')
        #
        lblframe.grid_columnconfigure(0,weight=1)
        for ind in range(len(self.subGraph.curves)+1):
            lblframe.grid_rowconfigure(ind,weight=1)
        #
        lblframelvl1.append(lblframe)
        #
        self.frame.legend_displayItem=[]
        #
        for ind in range(len(self.subGraph.curves)):
            var = TK.IntVar()
            var.set(1)
            CB = TTK.Checkbutton(lblframe,variable=var,command=lambda n=ind: self.cb_legend_display(n))#, variable=var)
            CB.val = var
            CB.var = 'legend_display'
            CB.ind=ind
            CB.grid(row=ind,column=0,sticky='NSEW')
            CB.container = self.frame.legend_displayItem
            self.frame.legend_displayItem.append(CB)
        # Curve to add
        ind = len(self.subGraph.curves)
        var = TK.IntVar()
        var.set(1)
        CB = TTK.Checkbutton(lblframe,variable=var,state=TK.DISABLED,command=lambda n=ind: self.cb_legend_display(n))#, variable=var)
        CB.val = var
        CB.var = 'legend_display'
        CB.ind=ind
        CB.grid(row=ind,column=0,sticky='NSEW')
        CB.container = self.frame.legend_displayItem
        self.frame.legend_displayItem.append(CB)

        ########################################################################
        #################### -> Button Line
        bottomFrame = TTK.Frame(self.frame)
        bottomFrame.grid(row=1,column=0,columnspan=18,sticky="NSEW")
        #
        bottomFrame.grid_columnconfigure(0,weight=1)
        bottomFrame.grid_columnconfigure(1,weight=1)
        bottomFrame.grid_columnconfigure(2,weight=1)
        bottomFrame.grid_columnconfigure(3,weight=1)
        bottomFrame.grid_columnconfigure(4,weight=1)
        bottomFrame.grid_columnconfigure(5,weight=1)

        #
        bottomFrame.grid_rowconfigure(0,weight=0)
        #
        B = TTK.Button(bottomFrame,text='Move Up',command=self.cmd_moveUp)
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        #
        B = TTK.Button(bottomFrame,text='Move Down',command=self.cmd_moveDown)
        B.grid(row=0,column=1,columnspan=1,sticky="nsew")
        #
        B = TTK.Button(bottomFrame,text='Delete selected curves',command=self.cmd_rmCurves)
        B.grid(row=0,column=2,columnspan=1,sticky="nsew")
        #
        B = TTK.Button(bottomFrame,text='Duplicate selected curves',command=self.cmd_duplicateCurves)
        B.grid(row=0,column=3,columnspan=1,sticky="nsew")
        #
        B = TTK.Button(bottomFrame,text='Create curve to add',command=self.cmd_createCurve)
        B.grid(row=0,column=4,columnspan=1,sticky="nsew")
        #
        B = TTK.Button(bottomFrame,text='Close',command=self.cmd_close)
        B.grid(row=0,column=5,columnspan=1,sticky="nsew")

        #
        try:
            self.frameList[self.graph][self.zone] = self.frame
        except KeyError:
            self.frameList[self.graph] = {self.zone:self.frame}

    # --------------------------------------------------------------- updateAxisList
    def updateAxisList(self):
        for ind in range(len(self.subGraph.curves)+1):
            self.frame.axisItem[ind].list = [i for i in range(len(self.subGraph.axis))]
    # --------------------------------------------------------------- updatelblFrameSize
    def updatelblFrameSize(self):
        for action in [self.frame.selectionItem,self.frame.IdItem,self.frame.axisItem,
                        self.frame.visibilityItem,self.frame.zoneItem,self.frame.varxItem,
                        self.frame.varyItem,self.frame.line_colorItem,self.frame.line_widthItem,
                        self.frame.line_styleItem,self.frame.marker_face_colorItem,self.frame.marker_sizeItem,
                        self.frame.marker_styleItem,self.frame.marker_edge_colorItem,self.frame.marker_edge_widthItem,
                        self.frame.marker_sampling_startItem,self.frame.marker_sampling_endItem,self.frame.marker_sampling_stepItem,
                        self.frame.legend_labelItem,self.frame.legend_displayItem]:
            ### get lblframe
            lblframe = action[0].winfo_parent() # Returns the name of the parent
            lblframe = self.frame.nametowidget(lblframe) # returns the instance of the parent # Returns the name of the parent
        #     ### Loop on curves
        #     for ind in range(len(self.subGraph.curves)):
        #         lblframe.rowconfigure(ind,weight=0)
        #         print('-> ',ind)
        #     lblframe.rowconfigure(len(self.subGraph.curves),weight=0)
        #     print('-> ',len(self.subGraph.curves))
        #     #
        self.frame.grid_rowconfigure(0,weight=len(self.subGraph.curves)+1)
        self.frame.grid_rowconfigure(1,weight=0)
    # --------------------------------------------------------------- popUpCurveLine
    def popUpCurveLine(self,indRemove):
        self.geometry("")
        cm = plt.get_cmap(COLOR_MAP)


        ### grid forget on indRemove
        for action in [self.frame.selectionItem,self.frame.IdItem,self.frame.axisItem,
                       self.frame.visibilityItem,self.frame.zoneItem,self.frame.varxItem,
                       self.frame.varyItem,self.frame.line_colorItem,self.frame.line_widthItem,
                       self.frame.line_styleItem,self.frame.marker_face_colorItem,self.frame.marker_sizeItem,
                       self.frame.marker_styleItem,self.frame.marker_edge_colorItem,self.frame.marker_edge_widthItem,
                       self.frame.marker_sampling_startItem,self.frame.marker_sampling_endItem,self.frame.marker_sampling_stepItem,
                       self.frame.legend_labelItem,self.frame.legend_displayItem]:
            action[indRemove].grid_forget()

        # action[indRemove].destroy()
        ### Loop on curves
        for ind in range(len(self.subGraph.curves)+1):
            ### ### Check index value
            if ind>=indRemove:
                for action in [self.frame.selectionItem,self.frame.IdItem,self.frame.axisItem,
                               self.frame.visibilityItem,self.frame.zoneItem,self.frame.varxItem,
                               self.frame.varyItem,self.frame.line_colorItem,self.frame.line_widthItem,
                               self.frame.line_styleItem,self.frame.marker_face_colorItem,self.frame.marker_sizeItem,
                               self.frame.marker_styleItem,self.frame.marker_edge_colorItem,self.frame.marker_edge_widthItem,
                               self.frame.marker_sampling_startItem,self.frame.marker_sampling_endItem,self.frame.marker_sampling_stepItem,
                               self.frame.legend_labelItem,self.frame.legend_displayItem]:
                    ### ### ### grid_forget
                    action[ind+1].grid_forget()
                    ### ### ### Pop up
                    action[ind] = action[ind+1]
                    action[ind+1] = None
                    action[ind].grid(row=ind,column=0,sticky='NSEW')
                    action[ind].ind = ind
                # Re Link the lambda function
                self.frame.selectionItem[ind].config(command=lambda n=ind: self.cb_selection(n))
                self.frame.IdItem[ind].config(text=ind)
                self.frame.axisItem[ind].config(command=lambda n=(ind,self.frame.axisItem): self.bt_click(n))
                self.frame.visibilityItem[ind].config(command=lambda n=ind: self.cb_visibility(n))
                self.frame.zoneItem[ind].config(command=lambda n=(ind,self.frame.zoneItem): self.bt_click(n))
                self.frame.varxItem[ind].config(command=lambda n=(ind,self.frame.varxItem): self.bt_click(n))
                self.frame.varyItem[ind].config(command=lambda n=(ind,self.frame.varyItem): self.bt_click(n))
                self.frame.line_colorItem[ind].config(command=lambda n=(ind,self.frame.line_colorItem): self.bt_click(n))
                self.frame.line_widthItem[ind].config(command=lambda n=(ind,self.frame.line_widthItem): self.bt_click(n))
                self.frame.line_styleItem[ind].config(command=lambda n=(ind,self.frame.line_styleItem): self.bt_click(n))
                self.frame.marker_face_colorItem[ind].config(command=lambda n=(ind,self.frame.marker_face_colorItem): self.bt_click(n))
                self.frame.marker_sizeItem[ind].config(command=lambda n=(ind,self.frame.marker_sizeItem): self.bt_click(n))
                self.frame.marker_styleItem[ind].config(command=lambda n=(ind,self.frame.marker_styleItem): self.bt_click(n))
                self.frame.marker_edge_colorItem[ind].config(command=lambda n=(ind,self.frame.marker_edge_colorItem): self.bt_click(n))
                self.frame.marker_edge_widthItem[ind].config(command=lambda n=(ind,self.frame.marker_edge_widthItem): self.bt_click(n))
                self.frame.marker_sampling_startItem[ind].config(command=lambda n=(ind,self.frame.marker_sampling_startItem): self.bt_click(n))
                self.frame.marker_sampling_endItem[ind].config(command=lambda n=(ind,self.frame.marker_sampling_endItem): self.bt_click(n))
                self.frame.marker_sampling_stepItem[ind].config(command=lambda n=(ind,self.frame.marker_sampling_stepItem): self.bt_click(n))
                self.frame.legend_labelItem[ind].config(command=lambda n=(ind,self.frame.legend_labelItem): self.bt_click(n))
                self.frame.legend_displayItem[ind].config(command=lambda n=ind: self.cb_legend_display(n))
        for action in [self.frame.selectionItem,self.frame.IdItem,self.frame.axisItem,
                       self.frame.visibilityItem,self.frame.zoneItem,self.frame.varxItem,
                       self.frame.varyItem,self.frame.line_colorItem,self.frame.line_widthItem,
                       self.frame.line_styleItem,self.frame.marker_face_colorItem,self.frame.marker_sizeItem,
                       self.frame.marker_styleItem,self.frame.marker_edge_colorItem,self.frame.marker_edge_widthItem,
                       self.frame.marker_sampling_startItem,self.frame.marker_sampling_endItem,self.frame.marker_sampling_stepItem,
                       self.frame.legend_labelItem,self.frame.legend_displayItem]:

            del action[-1]
            ### get lblframe
            lblframe = action[0].winfo_parent() # Returns the name of the parent
            lblframe = self.frame.nametowidget(lblframe) # returns the instance of the parent # Returns the name of the parent
            lblframe.grid_rowconfigure(len(self.subGraph.curves)+1,weight=0)
        #
        try:
            self.frameList[self.graph][self.zone] = self.frame
        except KeyError:
            self.frameList[self.graph] = {self.zone:self.frame}

    # --------------------------------------------------------------- addLineToFrame
    def addCurveLineToFrame(self,curve):
        # self.geometry("")
        cm = plt.get_cmap(COLOR_MAP)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Selection
        ### delete curve to add
        ind = len(self.subGraph.curves)-1
        ### Add new curve line
        var = TK.IntVar()
        self.frame.selectionItem[ind].config(state=TK.NORMAL)
        self.frame.selectionItem[ind].val = var
        self.frame.selectionItem[ind].ind=ind
        self.frame.selectionItem[ind].var='selection'
        lblframe = self.frame.selectionItem[ind].winfo_parent() # Returns the name of the parent
        lblframe = self.frame.nametowidget(lblframe) # returns the instance of the parent # Returns the name of the parent
        # Curve to add
        ind = len(self.subGraph.curves)
        var = TK.IntVar()
        CB = TTK.Checkbutton(lblframe,variable=var,command=lambda n=ind: self.cb_selection(n),state=TK.DISABLED)#, variable=var)
        CB.val = var
        CB.ind=ind
        CB.var='selection'
        CB.grid(row=ind,column=0,sticky='NSEW')
        CB.container = self.frame.selectionItem
        self.frame.selectionItem.append(CB)
        lblframe.grid_rowconfigure(ind,weight=1)

        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Id
        ### delete curve to add
        ind = len(self.subGraph.curves)-1
        ### Add new curve line
        self.frame.IdItem[ind].config(text='%s'%ind)
        self.frame.IdItem[ind].ind = ind
        lblframe = self.frame.IdItem[ind].winfo_parent() # Returns the name of the parent
        lblframe = self.frame.nametowidget(lblframe) # returns the instance of the parent # Returns the name of the parent
        # Curve to add
        ind = len(self.subGraph.curves)
        LBL = TK.Label(lblframe,text='to add')
        LBL.ind = ind
        LBL.grid(row=ind,column=0,sticky='NSEW')
        LBL.container = self.frame.IdItem
        self.frame.IdItem.append(LBL)
        lblframe.grid_rowconfigure(ind,weight=1)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Axis to plot
        ### delete curve to add
        ind = len(self.subGraph.curves)-1
        ### Add new curve line
        c = self.subGraph.curves[ind]
        self.frame.axisItem[ind].config(text=c.axis)
        self.frame.axisItem[ind].list = [i for i in range(len(self.subGraph.axis))]
        self.frame.axisItem[ind].val = c.axis
        self.frame.axisItem[ind].var = 'ind_axis'
        self.frame.axisItem[ind].ind = ind
        self.frame.axisItem[ind].treatmentId = 0
        lblframe = self.frame.axisItem[ind].winfo_parent() # Returns the name of the parent
        lblframe = self.frame.nametowidget(lblframe) # returns the instance of the parent # Returns the name of the parent
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=0,command=lambda n=(ind,self.frame.axisItem): self.bt_click(n))
        B.list = [i for i in range(len(self.subGraph.axis))]
        indexNone = 0
        B.val = 0
        B.var = 'axis'
        B.ind = ind
        B.treatmentId = 0
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.axisItem
        self.frame.axisItem.append(B)
        lblframe.grid_rowconfigure(ind,weight=1)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Visibility
        ### delete curve to add
        ind = len(self.subGraph.curves)-1
        ### Add new curve line
        var = TK.IntVar()
        var.set(1)
        self.frame.visibilityItem[ind].config(state=TK.NORMAL)
        self.frame.visibilityItem[ind].val = var
        self.frame.visibilityItem[ind].var = 'visible'
        self.frame.visibilityItem[ind].ind=ind
        lblframe = self.frame.visibilityItem[ind].winfo_parent() # Returns the name of the parent
        lblframe = self.frame.nametowidget(lblframe) # returns the instance of the parent # Returns the name of the parent
        # Curve to add
        ind = len(self.subGraph.curves)
        var = TK.IntVar()
        var.set(1)
        # CB = TTK.Checkbutton(lblframe,variable=var,state=TK.DISABLED,command=lambda n=ind: self.cb_visibility(n))#, variable=var)
        CB = TTK.Checkbutton(lblframe,variable=var,state=TK.DISABLED,command=lambda n=ind: self.cb_visibility(n))#, variable=var)
        CB.val = var
        CB.var = 'visible'
        CB.ind=ind
        CB.grid(row=ind,column=0,sticky='NSEW')
        CB.container = self.frame.visibilityItem
        self.frame.visibilityItem.append(CB)
        lblframe.grid_rowconfigure(ind,weight=1)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Zone
        ### delete curve to add
        ind = len(self.subGraph.curves)-1
        ### Add new curve line
        c = self.subGraph.curves[ind]
        self.frame.zoneItem[ind].config(text=len(c.zone))
        self.frame.zoneItem[ind].val = c.zone
        self.frame.zoneItem[ind].var = 'zone'
        self.frame.zoneItem[ind].ind = ind
        self.frame.zoneItem[ind].treatmentId = 4 # 4 is for selectzones
        lblframe = self.frame.zoneItem[ind].winfo_parent() # Returns the name of the parent
        lblframe = self.frame.nametowidget(lblframe) # returns the instance of the parent # Returns the name of the parent
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=len(self.parent.data.keys()),command=lambda n=(ind,self.frame.zoneItem): self.bt_click(n))
        B.val = self.parent.data.keys()
        B.var = 'zone'
        B.ind = ind
        B.treatmentId = 4
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.zoneItem
        self.frame.zoneItem.append(B)
        lblframe.grid_rowconfigure(ind,weight=1)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& VarX
        ### delete curve to add
        ind = len(self.subGraph.curves)-1
        ### Add new curve line
        c = self.subGraph.curves[ind]
        self.frame.varxItem[ind].config(text=c.varx)
        self.frame.varxItem[ind].ind = ind
        self.frame.varxItem[ind].list = self.filterVarWithZone(self.frame.varxItem[ind])
        self.frame.varxItem[ind].val = c.varx
        self.frame.varxItem[ind].var = 'varx'
        self.frame.varxItem[ind].treatmentId = 0
        lblframe = self.frame.varxItem[ind].winfo_parent() # Returns the name of the parent
        lblframe = self.frame.nametowidget(lblframe) # returns the instance of the parent # Returns the name of the parent
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=None,command=lambda n=(ind,self.frame.varxItem): self.bt_click(n))
        B.ind = ind
        B.list = self.filterVarWithZone(B)
        if len(B.list)!=0:
            B.config(text=B.list[0])
            B.val = B.list[0]
        else:
            B.config(text='')
            B.val = None
        B.var = 'varx'
        B.treatmentId = 0
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.varxItem
        self.frame.varxItem.append(B)
        lblframe.grid_rowconfigure(ind,weight=1)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& VarY
        ### delete curve to add
        ind = len(self.subGraph.curves)-1
        ### Add new curve line
        c = self.subGraph.curves[ind]
        self.frame.varyItem[ind].config(text=c.vary)
        self.frame.varyItem[ind].ind = ind
        self.frame.varyItem[ind].list = self.filterVarWithZone(self.frame.varyItem[ind])
        self.frame.varyItem[ind].val = c.vary
        self.frame.varyItem[ind].var = 'vary'
        self.frame.varyItem[ind].ind = ind
        self.frame.varyItem[ind].treatmentId = 0
        lblframe = self.frame.varyItem[ind].winfo_parent() # Returns the name of the parent
        lblframe = self.frame.nametowidget(lblframe) # returns the instance of the parent # Returns the name of the parent
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=None,command=lambda n=(ind,self.frame.varyItem): self.bt_click(n))
        B.ind = ind
        B.list = self.filterVarWithZone(B)
        if len(B.list)!=0:
            B.config(text=B.list[0])
            B.val = B.list[0]
        else:
            B.config(text='')
            B.val = None
        B.var = 'vary'
        B.treatmentId = 0
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.varyItem
        self.frame.varyItem.append(B)
        lblframe.grid_rowconfigure(ind,weight=1)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Line color
        ### delete curve to add
        ind = len(self.subGraph.curves)-1
        ### Add new curve line
        c = self.subGraph.curves[ind]
        self.frame.line_colorItem[ind].list = []
        self.frame.line_colorItem[ind].val = c.line_color
        self.frame.line_colorItem[ind].var = 'line_color'
        self.frame.line_colorItem[ind].config(bg=self.frame.line_colorItem[ind].val)
        self.frame.line_colorItem[ind].config(activebackground=self.frame.line_colorItem[ind].val)
        self.frame.line_colorItem[ind].ind = ind
        self.frame.line_colorItem[ind].treatmentId = 1
        lblframe = self.frame.line_colorItem[ind].winfo_parent() # Returns the name of the parent
        lblframe = self.frame.nametowidget(lblframe) # returns the instance of the parent # Returns the name of the parent
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TK.Button(lblframe,command=lambda n=(ind,self.frame.line_colorItem): self.bt_click(n))
        B.list = []
        color = cm(1.*(ind%NUM_COLORS)/NUM_COLORS)
        B.val = '#%02x%02x%02x' % (int(255*color[0]),int(255*color[1]),int(255*color[2]))
        B.var = 'line_color'
        B.config(bg=B.val)
        B.config(activebackground=B.val)
        B.ind = ind
        B.treatmentId = 1
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.line_colorItem
        self.frame.line_colorItem.append(B)
        lblframe.grid_rowconfigure(ind,weight=1)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Line size
        ### delete curve to add
        ind = len(self.subGraph.curves)-1
        ### Add new curve line
        c = self.subGraph.curves[ind]
        self.frame.line_widthItem[ind].config(text=c.line_width)
        self.frame.line_widthItem[ind].list = (np.arange(0.5,100,0.5)).tolist()
        self.frame.line_widthItem[ind].val = c.line_width
        self.frame.line_widthItem[ind].var = 'line_width'
        self.frame.line_widthItem[ind].ind = ind
        self.frame.line_widthItem[ind].treatmentId = 0
        lblframe = self.frame.line_widthItem[ind].winfo_parent() # Returns the name of the parent
        lblframe = self.frame.nametowidget(lblframe) # returns the instance of the parent # Returns the name of the parent
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=default_values['Curve']['line_width'],command=lambda n=(ind,self.frame.line_widthItem): self.bt_click(n))
        B.list = (np.arange(0.5,100,0.5)).tolist()
        B.val = default_values['Curve']['line_width']
        B.var = 'line_width'
        B.ind = ind
        B.treatmentId = 0
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.line_widthItem
        self.frame.line_widthItem.append(B)
        lblframe.grid_rowconfigure(ind,weight=1)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Line style
        ### delete curve to add
        ind = len(self.subGraph.curves)-1
        ### Add new curve line
        c = self.subGraph.curves[ind]
        self.frame.line_styleItem[ind].config(text=c.line_style)
        self.frame.line_styleItem[ind].list = linestylelist
        self.frame.line_styleItem[ind].val = c.line_style
        self.frame.line_styleItem[ind].var = 'line_style'
        self.frame.line_styleItem[ind].ind = ind
        self.frame.line_styleItem[ind].treatmentId = 0
        lblframe = self.frame.line_styleItem[ind].winfo_parent() # Returns the name of the parent
        lblframe = self.frame.nametowidget(lblframe) # returns the instance of the parent # Returns the name of the parent
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=default_values['Curve']['line_style'],command=lambda n=(ind,self.frame.line_styleItem): self.bt_click(n))
        B.list = linestylelist
        B.val = default_values['Curve']['line_style']
        B.var = 'line_style'
        B.ind = ind
        B.treatmentId = 0
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.line_styleItem
        self.frame.line_styleItem.append(B)
        lblframe.grid_rowconfigure(ind,weight=1)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Symbol Face color
        ### delete curve to add
        ind = len(self.subGraph.curves)-1
        ### Add new curve line
        c = self.subGraph.curves[ind]
        self.frame.marker_face_colorItem[ind].list = []
        self.frame.marker_face_colorItem[ind].val = c.marker_face_color
        self.frame.marker_face_colorItem[ind].var = 'marker_face_color'
        self.frame.marker_face_colorItem[ind].config(bg=self.frame.marker_face_colorItem[ind].val)
        self.frame.marker_face_colorItem[ind].config(activebackground=self.frame.marker_face_colorItem[ind].val)
        self.frame.marker_face_colorItem[ind].ind = ind
        self.frame.marker_face_colorItem[ind].treatmentId = 1
        lblframe = self.frame.marker_face_colorItem[ind].winfo_parent() # Returns the name of the parent
        lblframe = self.frame.nametowidget(lblframe) # returns the instance of the parent # Returns the name of the parent
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TK.Button(lblframe,command=lambda n=(ind,self.frame.marker_face_colorItem): self.bt_click(n))
        B.list = []
        color = cm(1.*(ind%NUM_COLORS)/NUM_COLORS)
        B.val = '#%02x%02x%02x' % (int(255*color[0]),int(255*color[1]),int(255*color[2]))
        B.var = 'marker_face_color'
        B.config(bg=B.val)
        B.config(activebackground=B.val)
        B.ind = ind
        B.treatmentId = 1
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.marker_face_colorItem
        self.frame.marker_face_colorItem.append(B)
        lblframe.grid_rowconfigure(ind,weight=1)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Symbol Face size
        ### delete curve to add
        ind = len(self.subGraph.curves)-1
        ### Add new curve line
        c = self.subGraph.curves[ind]
        self.frame.marker_sizeItem[ind].config(text=c.marker_size)
        self.frame.marker_sizeItem[ind].list = (np.arange(0.5,100,0.5)).tolist()
        self.frame.marker_sizeItem[ind].val = c.marker_size
        self.frame.marker_sizeItem[ind].var = 'marker_size'
        self.frame.marker_sizeItem[ind].ind = ind
        self.frame.marker_sizeItem[ind].treatmentId = 0
        lblframe = self.frame.marker_sizeItem[ind].winfo_parent() # Returns the name of the parent
        lblframe = self.frame.nametowidget(lblframe) # returns the instance of the parent # Returns the name of the parent
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=default_values['Curve']['marker_size'],command=lambda n=(ind,self.frame.marker_sizeItem): self.bt_click(n))
        B.list = (np.arange(0.5,100,0.5)).tolist()
        B.val = default_values['Curve']['marker_size']
        B.var = 'marker_size'
        B.ind = ind
        B.treatmentId = 0
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.marker_sizeItem
        self.frame.marker_sizeItem.append(B)
        lblframe.grid_rowconfigure(ind,weight=1)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Symbol Type
        ### delete curve to add
        ind = len(self.subGraph.curves)-1
        ### Add new curve line
        c = self.subGraph.curves[ind]
        self.frame.marker_styleItem[ind].config(text=c.marker_style)
        self.frame.marker_styleItem[ind].list = markername
        self.frame.marker_styleItem[ind].val = c.marker_style
        self.frame.marker_styleItem[ind].var = 'marker_style'
        self.frame.marker_styleItem[ind].ind = ind
        self.frame.marker_styleItem[ind].treatmentId = 0
        lblframe = self.frame.marker_styleItem[ind].winfo_parent() # Returns the name of the parent
        lblframe = self.frame.nametowidget(lblframe) # returns the instance of the parent # Returns the name of the parent
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=default_values['Curve']['marker_style'],command=lambda n=(ind,self.frame.marker_styleItem): self.bt_click(n))
        B.list = markername
        B.val = default_values['Curve']['marker_style']
        indexNone = markername.index(B.val)
        B.var = 'marker_style'
        B.ind = ind
        B.treatmentId = 0
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.marker_styleItem
        self.frame.marker_styleItem.append(B)
        lblframe.grid_rowconfigure(ind,weight=1)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Symbol edge color
        ### delete curve to add
        ind = len(self.subGraph.curves)-1
        ### Add new curve line
        c = self.subGraph.curves[ind]
        self.frame.marker_edge_colorItem[ind].list = []
        self.frame.marker_edge_colorItem[ind].val = c.marker_edge_color
        self.frame.marker_edge_colorItem[ind].var = 'marker_edge_color'
        self.frame.marker_edge_colorItem[ind].config(bg=self.frame.marker_edge_colorItem[ind].val)
        self.frame.marker_edge_colorItem[ind].config(activebackground=self.frame.marker_edge_colorItem[ind].val)
        self.frame.marker_edge_colorItem[ind].ind = ind
        self.frame.marker_edge_colorItem[ind].treatmentId = 1
        lblframe = self.frame.marker_edge_colorItem[ind].winfo_parent() # Returns the name of the parent
        lblframe = self.frame.nametowidget(lblframe) # returns the instance of the parent # Returns the name of the parent
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TK.Button(lblframe,command=lambda n=(ind,self.frame.marker_edge_colorItem): self.bt_click(n))
        B.list = []
        color = cm(1.*(ind%NUM_COLORS)/NUM_COLORS)
        B.val = '#%02x%02x%02x' % (int(255*color[0]),int(255*color[1]),int(255*color[2]))
        B.var = 'marker_edge_color'
        B.config(bg=B.val)
        B.config(activebackground=B.val)
        B.ind = ind
        B.treatmentId = 1
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.marker_edge_colorItem
        self.frame.marker_edge_colorItem.append(B)
        lblframe.grid_rowconfigure(ind,weight=1)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Symbol edge size
        ### delete curve to add
        ind = len(self.subGraph.curves)-1
        ### Add new curve line
        c = self.subGraph.curves[ind]
        self.frame.marker_edge_widthItem[ind].config(text=c.marker_edge_width)
        self.frame.marker_edge_widthItem[ind].list = (np.arange(0.5,100,0.5)).tolist()
        self.frame.marker_edge_widthItem[ind].val = c.marker_edge_width
        self.frame.marker_edge_widthItem[ind].var = 'marker_edge_width'
        self.frame.marker_edge_widthItem[ind].ind = ind
        self.frame.marker_edge_widthItem[ind].treatmentId = 0
        lblframe = self.frame.marker_edge_widthItem[ind].winfo_parent() # Returns the name of the parent
        lblframe = self.frame.nametowidget(lblframe) # returns the instance of the parent # Returns the name of the parent
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=default_values['Curve']['marker_edge_width'],command=lambda n=(ind,self.frame.marker_edge_widthItem): self.bt_click(n))
        B.list = (np.arange(0.5,100,0.5)).tolist()
        B.val = default_values['Curve']['marker_edge_width']
        B.var = 'marker_edge_width'
        B.ind = ind
        B.treatmentId = 0
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.marker_edge_widthItem
        self.frame.marker_edge_widthItem.append(B)
        lblframe.grid_rowconfigure(ind,weight=1)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Symbol sampling start
        ### delete curve to add
        ind = len(self.subGraph.curves)-1
        ### Add new curve line
        c = self.subGraph.curves[ind]
        self.frame.marker_sampling_startItem[ind].config(text=c.marker_sampling_start)
        self.frame.marker_sampling_startItem[ind].list = []
        self.frame.marker_sampling_startItem[ind].val = c.marker_sampling_start
        self.frame.marker_sampling_startItem[ind].var = 'marker_sampling_start'
        self.frame.marker_sampling_startItem[ind].ind = ind
        self.frame.marker_sampling_startItem[ind].treatmentId = 2
        lblframe = self.frame.marker_sampling_startItem[ind].winfo_parent() # Returns the name of the parent
        lblframe = self.frame.nametowidget(lblframe) # returns the instance of the parent # Returns the name of the parent
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=default_values['Curve']['marker_sampling_start'],command=lambda n=(ind,self.frame.marker_sampling_startItem): self.bt_click(n))
        B.list = []
        B.val = default_values['Curve']['marker_sampling_start']
        B.var = 'marker_sampling_start'
        B.ind = ind
        B.treatmentId = 2
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.marker_sampling_startItem
        self.frame.marker_sampling_startItem.append(B)
        lblframe.grid_rowconfigure(ind,weight=1)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Symbol sampling end
        ### delete curve to add
        ind = len(self.subGraph.curves)-1
        ### Add new curve line
        c = self.subGraph.curves[ind]
        self.frame.marker_sampling_endItem[ind].config(text=c.marker_sampling_end)
        self.frame.marker_sampling_endItem[ind].list = []
        self.frame.marker_sampling_endItem[ind].val = c.marker_sampling_end
        self.frame.marker_sampling_endItem[ind].var = 'marker_sampling_end'
        self.frame.marker_sampling_endItem[ind].ind = ind
        self.frame.marker_sampling_endItem[ind].treatmentId = 2
        lblframe = self.frame.marker_sampling_endItem[ind].winfo_parent() # Returns the name of the parent
        lblframe = self.frame.nametowidget(lblframe) # returns the instance of the parent # Returns the name of the parent
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=default_values['Curve']['marker_sampling_end'],command=lambda n=(ind,self.frame.marker_sampling_endItem): self.bt_click(n))
        B.list = []
        B.val = default_values['Curve']['marker_sampling_end']
        B.var = 'marker_sampling_end'
        B.ind = ind
        B.treatmentId = 2
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.marker_sampling_endItem
        self.frame.marker_sampling_endItem.append(B)
        lblframe.grid_rowconfigure(ind,weight=1)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Symbol sampling step
        ### delete curve to add
        ind = len(self.subGraph.curves)-1
        ### Add new curve line
        c = self.subGraph.curves[ind]
        self.frame.marker_sampling_stepItem[ind].config(text=c.marker_sampling_step)
        self.frame.marker_sampling_stepItem[ind].list = []
        self.frame.marker_sampling_stepItem[ind].val = c.marker_sampling_step
        self.frame.marker_sampling_stepItem[ind].var = 'marker_sampling_step'
        self.frame.marker_sampling_stepItem[ind].ind = ind
        self.frame.marker_sampling_stepItem[ind].treatmentId = 2
        lblframe = self.frame.marker_sampling_stepItem[ind].winfo_parent() # Returns the name of the parent
        lblframe = self.frame.nametowidget(lblframe) # returns the instance of the parent # Returns the name of the parent
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text=default_values['Curve']['marker_sampling_step'],command=lambda n=(ind,self.frame.marker_sampling_stepItem): self.bt_click(n))
        B.list = []
        B.val = default_values['Curve']['marker_sampling_step']
        B.var = 'marker_sampling_step'
        B.ind = ind
        B.treatmentId = 2
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.marker_sampling_stepItem
        self.frame.marker_sampling_stepItem.append(B)
        lblframe.grid_rowconfigure(ind,weight=1)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Legend label
        ### delete curve to add
        ind = len(self.subGraph.curves)-1
        ### Add new curve line
        c = self.subGraph.curves[ind]
        self.frame.legend_labelItem[ind].config(text=c.legend_label)
        self.frame.legend_labelItem[ind].list = []
        self.frame.legend_labelItem[ind].val = c.legend_label
        self.frame.legend_labelItem[ind].var = 'legend_label'
        self.frame.legend_labelItem[ind].ind = ind
        self.frame.legend_labelItem[ind].treatmentId = 3
        lblframe = self.frame.legend_labelItem[ind].winfo_parent() # Returns the name of the parent
        lblframe = self.frame.nametowidget(lblframe) # returns the instance of the parent # Returns the name of the parent
        # Curve to add
        ind = len(self.subGraph.curves)
        B = TTK.Button(lblframe,text="...",command=lambda n=(ind,self.frame.legend_labelItem): self.bt_click(n))
        B.list = []
        B.val = "..."
        B.var = 'legend_label'
        B.ind = ind
        B.treatmentId = 3
        B.grid(row=ind,column=0,columnspan=1,sticky="nsew")
        B.container = self.frame.legend_labelItem
        self.frame.legend_labelItem.append(B)
        lblframe.grid_rowconfigure(ind,weight=1)
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Legend display
        ### delete curve to add
        ind = len(self.subGraph.curves)-1
        ### Add new curve line
        var = TK.IntVar()
        var.set(1)
        self.frame.legend_displayItem[ind].config(state=TK.NORMAL)
        self.frame.legend_displayItem[ind].val = var
        self.frame.legend_displayItem[ind].var = 'legend_display'
        self.frame.legend_displayItem[ind].ind=ind
        lblframe = self.frame.legend_displayItem[ind].winfo_parent() # Returns the name of the parent
        lblframe = self.frame.nametowidget(lblframe) # returns the instance of the parent # Returns the name of the parent
        # Curve to add
        ind = len(self.subGraph.curves)
        var = TK.IntVar()
        var.set(1)
        CB = TTK.Checkbutton(lblframe,variable=var,state=TK.DISABLED,command=lambda n=ind: self.cb_legend_display(n))#, variable=var)
        CB.val = var
        CB.var = 'legend_display'
        CB.ind=ind
        CB.grid(row=ind,column=0,sticky='NSEW')
        CB.container = self.frame.legend_displayItem
        self.frame.legend_displayItem.append(CB)
        lblframe.grid_rowconfigure(ind,weight=1)

        try:
            self.frameList[self.graph][self.zone] = self.frame
        except KeyError:
            self.frameList[self.graph] = {self.zone:self.frame}

    # --------------------------------------------------------------- cmd_moveUp
    def cmd_moveUp(self,event=None):
        # Works only if a single curve is selected !
        nbSelected = 0
        indSelected = None
        for ind in range(len(self.frame.selectionItem)):
            select = self.frame.selectionItem[ind]
            if select.val.get():
                nbSelected += 1
                indSelected = ind
        # Check if only a single curve has been selected
        if nbSelected != 1:
            tkMessageBox.showwarning('Move Up Failed','Need to select a single curve')
            return
        else:
            # if indSelected == 0, can not move up beacause it is already at the top !
            if indSelected!= 0:
                c = cpy.deepcopy(self.subGraph.curves[indSelected])
                self.subGraph.curves[indSelected]=self.subGraph.curves[indSelected-1]
                self.subGraph.curves[indSelected-1]=c
                self.switchCurve(indSelected,-1)
                self.updatelblFrameSize()
                try:
                    self.frameList[self.graph][self.zone] = self.frame
                except KeyError:
                    self.frameList[self.graph] = {self.zone:self.frame}
                #
                # Update Graph
                self.parent.graphWdwL[self.graph].updateGraph(self.zone)
    # ------------------------------------------------------------- switchUpCurve
    def switchCurve(self,indSelected,incr):
        self.geometry("")
        cm = plt.get_cmap(COLOR_MAP)
        if (indSelected==0 and incr == -1):
            return
        if (indSelected==(len(self.frame.selectionItem)-1) and incr==1):
            return
        ### Visual
        for action in [self.frame.selectionItem,self.frame.IdItem,self.frame.axisItem,
                       self.frame.visibilityItem,self.frame.zoneItem,self.frame.varxItem,
                       self.frame.varyItem,self.frame.line_colorItem,self.frame.line_widthItem,
                       self.frame.line_styleItem,self.frame.marker_face_colorItem,self.frame.marker_sizeItem,
                       self.frame.marker_styleItem,self.frame.marker_edge_colorItem,self.frame.marker_edge_widthItem,
                       self.frame.marker_sampling_startItem,self.frame.marker_sampling_endItem,self.frame.marker_sampling_stepItem,
                       self.frame.legend_labelItem,self.frame.legend_displayItem]:
            action[indSelected].grid_forget()
            action[indSelected+incr].grid_forget()
            action[indSelected].grid(row=indSelected+incr,column=0,sticky='NSEW')
            action[indSelected+incr].grid(row=indSelected,column=0,sticky='NSEW')
            action[indSelected].ind = indSelected+incr
            action[indSelected+incr].ind = indSelected
            tmp = action[indSelected+incr]
            action[indSelected+incr] = action[indSelected]
            action[indSelected] = tmp
        ### Edit lambda functions
        for ind in [indSelected+incr,indSelected]:
            self.frame.selectionItem[ind].config(command=lambda n=ind: self.cb_selection(n))
            self.frame.IdItem[ind].config(text=ind)
            self.frame.axisItem[ind].config(command=lambda n=(ind,self.frame.axisItem): self.bt_click(n))
            self.frame.visibilityItem[ind].config(command=lambda n=ind: self.cb_visibility(n))
            self.frame.zoneItem[ind].config(command=lambda n=(ind,self.frame.zoneItem): self.bt_click(n))
            self.frame.varxItem[ind].config(command=lambda n=(ind,self.frame.varxItem): self.bt_click(n))
            self.frame.varyItem[ind].config(command=lambda n=(ind,self.frame.varyItem): self.bt_click(n))
            self.frame.line_colorItem[ind].config(command=lambda n=(ind,self.frame.line_colorItem): self.bt_click(n))
            self.frame.line_widthItem[ind].config(command=lambda n=(ind,self.frame.line_widthItem): self.bt_click(n))
            self.frame.line_styleItem[ind].config(command=lambda n=(ind,self.frame.line_styleItem): self.bt_click(n))
            self.frame.marker_face_colorItem[ind].config(command=lambda n=(ind,self.frame.marker_face_colorItem): self.bt_click(n))
            self.frame.marker_sizeItem[ind].config(command=lambda n=(ind,self.frame.marker_sizeItem): self.bt_click(n))
            self.frame.marker_styleItem[ind].config(command=lambda n=(ind,self.frame.marker_styleItem): self.bt_click(n))
            self.frame.marker_edge_colorItem[ind].config(command=lambda n=(ind,self.frame.marker_edge_colorItem): self.bt_click(n))
            self.frame.marker_edge_widthItem[ind].config(command=lambda n=(ind,self.frame.marker_edge_widthItem): self.bt_click(n))
            self.frame.marker_sampling_startItem[ind].config(command=lambda n=(ind,self.frame.marker_sampling_startItem): self.bt_click(n))
            self.frame.marker_sampling_endItem[ind].config(command=lambda n=(ind,self.frame.marker_sampling_endItem): self.bt_click(n))
            self.frame.marker_sampling_stepItem[ind].config(command=lambda n=(ind,self.frame.marker_sampling_stepItem): self.bt_click(n))
            self.frame.legend_labelItem[ind].config(command=lambda n=(ind,self.frame.legend_labelItem): self.bt_click(n))
            self.frame.legend_displayItem[ind].config(command=lambda n=ind: self.cb_legend_display(n))

    # ------------------------------------------------------------- cmd_moveDown
    def cmd_moveDown(self,event=None):
        # Works only if a single curve is selected !
        nbSelected = 0
        indSelected = None
        for ind in range(len(self.frame.selectionItem)):
            select = self.frame.selectionItem[ind]
            if select.val.get():
                nbSelected += 1
                indSelected = ind
        # Check if only a single curve has been selected
        if nbSelected != 1:
            tkMessageBox.showwarning('Move Up Failed','Need to select a single curve')
            return
        else:
            # if indSelected == len(self.frame.selectionItem)-2, can not move down beacause it is already at the bottom !
            # Rmk : it is "-2" because there is the line for the "curve to add" that counts for one that is at the end
            if indSelected!= len(self.frame.selectionItem)-2:
                c = cpy.deepcopy(self.subGraph.curves[indSelected])
                self.subGraph.curves[indSelected]=self.subGraph.curves[indSelected+1]
                self.subGraph.curves[indSelected+1]=c
                self.switchCurve(indSelected,1)
                self.updatelblFrameSize()
                try:
                    self.frameList[self.graph][self.zone] = self.frame
                except KeyError:
                    self.frameList[self.graph] = {self.zone:self.frame}
                #
                # Update Graph
                self.parent.graphWdwL[self.graph].updateGraph(self.zone)
    # ---------------------------------------------------------- updateVarXYList
    def updateVarXYList(self,ind):
        self.frame.varxItem[ind].list = self.filterVarWithZone(self.frame.varxItem[ind])
        self.frame.varyItem[ind].list = self.filterVarWithZone(self.frame.varyItem[ind])
    # ----------------------------------------------------------- updateVarXYVal
    def updateVarXYVal(self,ind):
        needsUpdate = False
        if not self.frame.varxItem[ind].val in self.frame.varxItem[ind].list:
            self.frame.varxItem[ind].val = None
            self.frame.varxItem[ind].config(text='')
            needsUpdate = True
        if not self.frame.varyItem[ind].val in self.frame.varyItem[ind].list:
            self.frame.varyItem[ind].val = None
            self.frame.varyItem[ind].config(text='')
            needsUpdate = True
        if needsUpdate:
            try:
                self.subGraph.curves[ind].setValue(self.frame.varxItem[ind].var,self.frame.varxItem[ind].val)
                self.subGraph.curves[ind].setValue(self.frame.varyItem[ind].var,self.frame.varyItem[ind].val)
                # Update Graph
                self.parent.graphWdwL[self.graph].updateGraph(self.zone)
            except IndexError:
                return

    # -------------------------------------------------------- filterVarWithZone
    def filterVarWithZone(self,B):
        tmp = {}
        for zone in self.frame.zoneItem[B.ind].val:
            for var in self.parent.data[zone]:
                if not var in tmp:
                    tmp[var]=1
                else:
                    tmp[var]+=1
        #
        res = []
        nbzones = len(self.frame.zoneItem[B.ind].val)
        for var in tmp:
            if tmp[var]==nbzones:
                res.append(var)
        return sorted(res)
    # ----------------------------------------------------------------- bt_click
    def updateColor(self,color,B,extra_data):
        bt_list = extra_data[1]
        # B = bt_list[ind[0]]
        l_ind = [extra_data[0]]
        # If line is selected, apply the modification to all other selected lines
        if self.frame.selectionItem[extra_data[0]].val.get():
            for ind2 in range(len(self.frame.selectionItem)):
                if extra_data[0]!=ind2 and self.frame.selectionItem[ind2].val.get():
                    l_ind.append(ind2)
        if color is not None:
            for ii in l_ind:
                B = bt_list[ii]
                B.val = color
                B.config(bg=B.val)
                B.config(activebackground=B.val)
                try:
                    self.subGraph.curves[B.ind].setValue(B.var,B.val)
                    # Update Graph
                    self.parent.graphWdwL[self.graph].updateGraph(self.zone)

                except IndexError:
                    return
    # ----------------------------------------------------------------- bt_click
    def bt_click(self,ind):
        bt_list = ind[1]
        B = bt_list[ind[0]]
        self.closeAllDialog()
        if B.treatmentId==0:
            self.list_dialog = list_dialogWindow()
            self.list_dialog.initialize(self,B)
        elif B.treatmentId==1:
            self.list_dialog = color_dialogWindow(self,B,B.val,extra_data=ind)
            # try:
            #     color = askcolor(B.val,parent=self)
            #
            #     l_ind = [ind[0]]
            #     # If line is selected, apply the modification to all other selected lines
            #     if self.frame.selectionItem[ind[0]].val.get():
            #         for ind2 in range(len(self.frame.selectionItem)):
            #             if ind[0]!=ind2 and self.frame.selectionItem[ind2].val.get():
            #                 l_ind.append(ind2)
            #     if color[1] is not None:
            #         for ii in l_ind:
            #             B = bt_list[ii]
            #             B.val = color[1]
            #             B.config(bg=B.val)
            #             B.config(activebackground=B.val)
            #             try:
            #                 self.subGraph.curves[B.ind].setValue(B.var,B.val)
            #                 # Update Graph
            #                 self.parent.graphWdwL[self.graph].updateGraph(self.zone)
            #
            #             except IndexError:
            #                 return
            # except ValueError:
            #     return
        elif B.treatmentId==2:
            self.input_dialog = input_dialogWindow()
            self.input_dialog.initialize(self,B)
        elif B.treatmentId == 3:
            self.input_dialog = input_dialogStringWindow()
            self.input_dialog.initialize(self,B)
        elif B.treatmentId == 4:
            self.input_dialog = input_dialogSelectZoneWindow()
            self.input_dialog.initialize(self,B)
        else:
            return
    # ----------------------------------------------------------- closeAllDialog
    def closeAllDialog(self):
        if self.list_dialog:
            self.list_dialog.cmd_close()
            self.list_dialog = None
        if self.input_dialog:
            self.input_dialog.cmd_close()
            self.input_dialog = None
    # ------------------------------------------------------------ cb_visibility
    def cb_visibility(self,ind):
        CB = self.frame.visibilityItem[ind]
        initialValue = CB.val.get()
        CB.val.set(not initialValue)
        self.subGraph.curves[ind].setValue('visible',CB.val.get())
        # If line was selected, apply this modification to all other selected lines
        if self.frame.selectionItem[ind].val.get():
            for ind2 in range(len(self.frame.selectionItem)):
                if ind!=ind2 and self.frame.selectionItem[ind2].val.get():
                    CB = self.frame.visibilityItem[ind2]
                    CB.val.set(not initialValue)
                    if initialValue:
                        CB.state(['!selected'])
                    else:
                        CB.state(['selected'])
                    self.subGraph.curves[ind2].setValue('visible',CB.val.get())
        # Update Graph
        self.parent.graphWdwL[self.graph].updateGraph(self.zone)
    # ---------------------------------------------------------------- cb_legend
    def cb_legend_display(self,ind):
        CB = self.frame.legend_displayItem[ind]
        initialValue = CB.val.get()
        CB.val.set(not initialValue)
        self.subGraph.curves[ind].setValue('legend_display',CB.val.get())
        # If line was selected, apply this modification to all other selected lines
        if self.frame.selectionItem[ind].val.get():
            for ind2 in range(len(self.frame.selectionItem)):
                if ind!=ind2 and self.frame.selectionItem[ind2].val.get():
                    CB = self.frame.legend_displayItem[ind2]
                    CB.val.set(not initialValue)
                    if initialValue:
                        CB.state(['!selected'])
                    else:
                        CB.state(['selected'])
                    self.subGraph.curves[ind2].setValue('legend_display',CB.val.get())
        # Update Graph
        self.parent.graphWdwL[self.graph].updateGraph(self.zone)
    # ------------------------------------------------------------- cb_selection
    def cb_selection(self,ind):
        CB = self.frame.selectionItem[ind]
        CB.val.set(not CB.val.get())
    # ---------------------------------------------------------------cb_rmCurves
    def cmd_rmCurves(self):
        nbDeletion = 0

        deletionList = []
        for ind in range(len(self.frame.selectionItem)):
            select = self.frame.selectionItem[ind]
            if select.val.get():
                deletionList.append(ind)
        deletionList.sort()
        ind=0
        nbDeletion = len(deletionList)-1
        while ind<=nbDeletion:
            indRemove = deletionList[ind]
            del self.subGraph.curves[indRemove]
            self.popUpCurveLine(indRemove)
            ind += 1
            deletionList = [i-1 for i in deletionList]
        self.updatelblFrameSize()

        # print(deletionList)
        # print('*'*80)
        # size2delete = len(self.frame.selectionItem) - 1
        # ind = 0
        # while ind<size2delete:
        #     print('°°° ',ind)
        #     select = self.frame.selectionItem[ind-nbDeletion]
        #     print(ind,' : ',select.val.get())
        #     if select.val.get():
        #         indRemove = ind -nbDeletion
        #         nbDeletion += 1
        #         del self.subGraph.curves[indRemove]
        #         self.popUpCurveLine(indRemove)
        #     ind += 1
            #     print(' - ',ind)
            #     indRemove = ind -nbDeletion
            #     nbDeletion += 1
            #     del self.subGraph.curves[indRemove]
            #     self.popUpCurveLine(indRemove)
            #     self.updatelblFrameSize()
            # else:
            #     print(' + ',ind)
        # self.reloadWindow()
        # Update Graph
        self.parent.graphWdwL[self.graph].updateGraph(self.zone)

    # ---------------------------------------------------------- cmd_createCurve
    def cmd_createCurve(self):

        ind = len(self.subGraph.curves)

        c= Curve(
        zone                   = self.frame.zoneItem[ind].val,
        varx                   = self.frame.varxItem[ind].val,
        vary                   = self.frame.varyItem[ind].val,
        line_color             = self.frame.line_colorItem[ind].val,
        line_style             = self.frame.line_styleItem[ind].val,
        line_width             = self.frame.line_widthItem[ind].val,
        marker_style           = self.frame.marker_styleItem[ind].val,
        marker_size            = self.frame.marker_sizeItem[ind].val,
        marker_edge_width      = self.frame.marker_edge_widthItem[ind].val,
        marker_face_color      = self.frame.marker_face_colorItem[ind].val,
        marker_edge_color      = self.frame.marker_edge_colorItem[ind].val,
        marker_sampling_start  = self.frame.marker_sampling_startItem[ind].val,
        marker_sampling_end    = self.frame.marker_sampling_endItem[ind].val,
        marker_sampling_step   = self.frame.marker_sampling_stepItem[ind].val,
        legend_label           = self.frame.legend_labelItem[ind].val,
        legend_display         = True,
        visible                = True)
        # Add curves
        self.subGraph.curves.append(c)
        # self.frame.destroy()
        self.addCurveLineToFrame(c)
        self.updatelblFrameSize()
        # self.createFrame()
        # Update Graph
        self.parent.graphWdwL[self.graph].updateGraph(self.zone)
    # ------------------------------------------------------------- reloadWindow
    def reloadWindow(self):
        self.frame.destroy()
        self.createFrame()
    # ------------------------------------------------------ cmd_duplicateCurves
    def cmd_duplicateCurves(self):

        for ind in range(len(self.frame.selectionItem)):
            select = self.frame.selectionItem[ind]
            if select.val.get():

                c= Curve(
                zone                   = self.frame.zoneItem[ind].val,
                varx                   = self.frame.varxItem[ind].val,
                vary                   = self.frame.varyItem[ind].val,
                line_color             = self.frame.line_colorItem[ind].val,
                line_style             = self.frame.line_styleItem[ind].val,
                line_width             = self.frame.line_widthItem[ind].val,
                marker_style           = self.frame.marker_styleItem[ind].val,
                marker_size            = self.frame.marker_sizeItem[ind].val,
                marker_edge_width      = self.frame.marker_edge_widthItem[ind].val,
                marker_face_color      = self.frame.marker_face_colorItem[ind].val,
                marker_edge_color      = self.frame.marker_edge_colorItem[ind].val,
                marker_sampling_start  = self.frame.marker_sampling_startItem[ind].val,
                marker_sampling_end    = self.frame.marker_sampling_endItem[ind].val,
                marker_sampling_step   = self.frame.marker_sampling_stepItem[ind].val,
                legend_label           = self.frame.legend_labelItem[ind].val,
                legend_display         = True,
                visible                = True)
                # Add curves
                self.subGraph.curves.append(c)
                # self.frame.destroy()
                self.addCurveLineToFrame(c)
                self.updatelblFrameSize()
        # self.createFrame()
        # Update Graph
        self.parent.graphWdwL[self.graph].updateGraph(self.zone)
    # ---------------------------------------------------------------- cmd_close
    def cmd_close(self):
        self.closeAllDialog()
        self.parent.editCurveWdw = None
        self.destroy()

    # -------------------------------------------------------------- updateButon
    def updateButon(self,B,val):

        l_ind = [B.ind]
        containerOfB = B.container
        # If line is selected, apply the modification to all other selected lines
        if self.frame.selectionItem[B.ind].val.get():
            for ind2 in range(len(self.frame.selectionItem)):
                if B.ind!=ind2 and self.frame.selectionItem[ind2].val.get():
                    l_ind.append(ind2)
        for ind in l_ind:
            B = containerOfB[ind]
            B.val = val
            if isinstance(B.val,list):
                # In this very special case, the B.val is the list of zones to plot
                B.config(text=len(B.val))
                self.updateVarXYList(B.ind)
                self.updateVarXYVal(B.ind)
            else:
                B.config(text=B.val)
            try:
                self.subGraph.curves[B.ind].setValue(B.var,B.val)
                # Update Graph
                self.parent.graphWdwL[self.graph].updateGraph(self.zone)
            except IndexError:
                return
# ==============================================================================
# ==============================================================================
class DesktopFrameTK(TK.Frame):
    def __init__(self,parent):
        TK.Frame.__init__(self,parent)
        self.initialize()
#        self.data = data # TODO : should be self.data = {} or None at the end of development
        self.data = None
        self.thread = None
        self.parent = parent
    # -------------------------------------------------------------- replaceGroupZones
    def replaceGroupZones(self,data,oldZoneList):
        if isinstance(data,list):
            # replace zone according to a tree
            self.replaceGroupZonesWithTree(data,oldZoneList)
        elif isinstance(data,dict):
            # Conformize the input dictionnary
            tmp = {}
            # Check if d structure is 'zone' oriented
            isZoneOriented = True
            for k in data:
                if not isinstance(data[k],dict): # then it is not zone oriented
                    isZoneOriented = False
                    break
            if not isZoneOriented:
                tmp[default_base]=data
            else:
                tmp = data
            # replace zone according to a dict data
            self.replaceGroupZonesWithDict(data,oldZoneList)
    # ------------------------------------------------------ replaceGroupZonesWithTree
    def replaceGroupZonesWithTree(self,d,oldZoneList):
        # add data and determine the list of new data
        newZoneList = self.addDataWithTree(d)
        # Get the curves that are concerned by a group of old zones and change it to the group of new zones
        self.updateGroupCurves(oldZoneList,newZoneList)
        # Compare old zones to new zones group and remove old zones that do not exist anymore
        for zoneName in oldZoneList:
            if zoneName not in newZoneList:
                # Delete old zone from data
                self.deleteZoneFromData(zoneName)
                # Delete old zone from curve
                self.deleteZoneInCurve(zoneName)
        ##### Redraw all
        self.updateAllGraph()
        #
    # ------------------------------------------------------ replaceGroupZonesWithDict
    def replaceGroupZonesWithDict(self,d,oldZoneList):
        # Add new data and determine the list of new zones
        newZoneList = []
        for zoneName in d:
            newZoneList.append(zoneName)
            self.addZoneWithDict(d,zoneName)
        # Get the curves that are concerned by a group of old zones and change it to the group of new zones
        self.updateGroupCurves(oldZoneList,newZoneList)
        # Compare old zones to new zones group and remove old zones that do not exist anymore
        for zoneName in oldZoneList:
            if zoneName not in newZoneList:
                # Delete old zone from data
                self.deleteZoneFromData(zoneName)
                # Delete old zone from curve
                self.deleteZoneInCurve(zoneName)
        ##### Redraw all
        self.updateAllGraph()
    # -------------------------------------------------------------- deleteZoneInCurve
    def deleteZoneInCurve(self,zoneName):
        for graph in self.graphWdwL:
            for ax_name in graph.fig.subGraph:
                graph.deleteZoneInCurve(ax_name,zoneName)
    # -------------------------------------------------------------- updateGroupCurves
    def updateGroupCurves(self,oldZoneList,newZoneList):
        for graph in self.graphWdwL:
            for ax_name in graph.fig.subGraph:
                graph.updateGroupCurvesZoneName(ax_name,oldZoneList,newZoneList)
    # ---------------------------------------------------------- addDataWithTree
    def addDataWithTree(self,t):
        tmp = self.data
        newZoneList = []
        for base in Internal.getNodesFromType1(t, 'CGNSBase_t'):
            basename = Internal.getName(base)
            ## ## Loop on zones
            for zone in Internal.getNodesFromType1(base,'Zone_t'):
                # Grab GridCoordinates
                zonename = Internal.getName(zone)
                ## ## Get GridCoorinates nodes
                try:
                    gridcoord = Internal.getNodesFromType1(zone,'GridCoordinates_t')[0]
                    for child in Internal.getChildren(gridcoord):
                        if not basename+'/'+zonename in tmp:
                            tmp[basename+'/'+zonename]={}
                        tmp[basename+'/'+zonename][Internal.getName(child)]=Internal.getValue(child)
                        newZoneList.append(basename+'/'+zonename)
                except IndexError: # No GridCoorinates node in this zone
                    # print('''gridcoord = Internal.getNodesFromType(zone,'GridCoordinates_t')[0] -----> Can not be loaded''')
                    pass

                # Grab FlowSolution for ZoneBC
                try:
                    zoneBC = Internal.getNodesFromType1(zone,'ZoneBC_t')[0]
                    for bc in Internal.getChildren(zoneBC):
                        bcname = Internal.getName(bc)
                        try:
                            bcdata = Internal.getNodesFromType(zoneBC,'BCData_t')[0]
                            for var in Internal.getChildren(bcdata):
                                if Internal.getType(var) == 'DataArray_t':
                                    if not basename+'/'+zonename in tmp:
                                        tmp[basename+'/'+zonename]={}
                                    tmp[basename+'/'+zonename][Internal.getName(var)+'@'+bcname]=Internal.getValue(var)
                                    newZoneList.append(basename+'/'+zonename)
                        except IndexError:
                            # print('''bcdata = Internal.getNodesFromType(zoneBC,'BCData_t')[0] ------> Can not be loaded''')
                            pass

                except IndexError: # No zoneBC in this zone
                    # print('''zoneBC = Internal.getNodesFromType(zone,'ZoneBC_t')[0] -----> Can not be loaded''')
                    pass
                # Grab FlowSolution (the rest of them)
                for flowsolution in Internal.getNodesFromType1(zone,'FlowSolution_t'):
                    flowsolutionname = Internal.getName(flowsolution)
                    for var in Internal.getChildren(flowsolution):
                        if Internal.getType(var)=='DataArray_t':
                            if not basename+'/'+zonename in tmp:
                                tmp[basename+'/'+zonename]={}
                            tmp[basename+'/'+zonename][Internal.getName(var)+'@'+flowsolutionname]=Internal.getValue(var)
                            newZoneList.append(basename+'/'+zonename)


        self.data = OrderedDict(sorted(tmp.items(),key=lambda t : t[0]))
        return newZoneList
    # ------------------------------------------------------------------ addZone
    def addZone(self,data,zoneName,basename='.*'):
        if self.data is None:
            self.data = {}
        if isinstance(data,list):
            # add zone according to a tree
            self.addZoneWithTree(data,zoneName,basename)
        elif isinstance(data,dict):
            # add zone according to a dict data
            self.addZoneWithDict(data,zoneName)
        if self.editCurveWdw  is not None:
            # self.editCurveWdw.createFrame()
            self.editCurveWdw.updateData()
    # ---------------------------------------------------------- addZoneWithTree
    def addZoneWithTree(self,t,zoneName,basenamefilter):
        for base in Internal.getNodesFromType1(t, 'CGNSBase_t'):
            basename = Internal.getName(base)
            if re.match(basenamefilter,basename):

                zone = Internal.getNodesFromName1(base, zoneName)[0]
                zonename = zoneName
                #
                ## ## Get GridCoorinates nodes
                try:
                    gridcoord = Internal.getNodesFromType2(zone,'GridCoordinates_t')[0]
                    for child in Internal.getChildren(gridcoord):
                        if not basename+'/'+zonename in self.data:
                            self.data[basename+'/'+zonename]={}
                        self.data[basename+'/'+zonename][Internal.getName(child)]=Internal.getValue(child)
                except IndexError: # No GridCoordinates node in this zone
                    # print('''gridcoord = Internal.getNodesFromType(zone,'GridCoordinates_t')[0] -----> Can not be loaded''')
                    pass

                # Grab FlowSolution for ZoneBC
                try:
                    zoneBC = Internal.getNodesFromType1(zone,'ZoneBC_t')[0]
                    for bc in Internal.getChildren(zoneBC):
                        bcname = Internal.getName(bc)
                        try:
                            bcdata = Internal.getNodesFromType(zoneBC,'BCData_t')[0]
                            for var in Internal.getChildren(bcdata):
                                if Internal.getType(var) == 'DataArray_t':
                                    if not basename+'/'+zonename in self.data:
                                        self.data[basename+'/'+zonename]={}
                                    self.data[basename+'/'+zonename][Internal.getName(var)+'@'+bcname]=Internal.getValue(var)
                        except IndexError:
                            # print('''bcdata = Internal.getNodesFromType(zoneBC,'BCData_t')[0] ------> Can not be loaded''')
                            pass

                except IndexError: # No zoneBC in this zone
                    # print('''zoneBC = Internal.getNodesFromType(zone,'ZoneBC_t')[0] -----> Can not be loaded''')
                    pass
                # Grab FlowSolution (the rest of them)
                for flowsolution in Internal.getNodesFromType1(zone,'FlowSolution_t'):
                    flowsolutionname = Internal.getName(flowsolution)
                    for var in Internal.getChildren(flowsolution):
                        if Internal.getType(var)=='DataArray_t':
                            if not basename+'/'+zonename in self.data:
                                self.data[basename+'/'+zonename]={}
                            self.data[basename+'/'+zonename][Internal.getName(var)+'@'+flowsolutionname]=Internal.getValue(var)
    # ---------------------------------------------------------- addZoneWithDict
    def addZoneWithDict(self,d,zoneName):
        if zoneName in d:
            self.data[zoneName]=d[zoneName]
        else:
            print('''### Can not find zone %s in submitted data.'''%newZoneName)
    # ------------------------------------------------------- deleteZoneFromData
    def deleteZoneFromData(self,zoneName,oldBaseName=""):
        for k in self.data:
            re_str = oldBaseName+zoneName.replace('\\','\\\\') # replace \ by \\ for regular expression conversion
            if re.match(re_str,k):
                del self.data[k]
    # -------------------------------------------------------------- replaceZone
    def replaceZone(self,data,oldZoneName,newZoneName,oldBaseName="",newBaseName=""):
        if isinstance(data,list):
            # replace zone according to a tree
            self.replaceZoneWithTree(data,oldBaseName,oldZoneName,newBaseName,newZoneName)
        elif isinstance(data,dict):
            # replace zone according to a dict data
            if oldBaseName != "": oldZoneName = oldBaseName+'/'+oldZoneName
            if newBaseName != "": newZoneName = newBaseName+'/'+newZoneName
            self.replaceZoneWithDict(data,oldZoneName,newZoneName)
        if self.editCurveWdw  is not None:
            # self.editCurveWdw.createFrame()
            self.editCurveWdw.updateData()
    # ------------------------------------------------------ replaceZoneWithDict
    def replaceZoneWithDict(self,d,oldZoneName,newZoneName):
        # Delete old zone from data
        self.deleteZoneFromData(oldZoneName)
        #
        self.addZoneWithDict(d,newZoneName)
        #
        ##### Edit curves : change the zonename
        self.updateAllCurvesZoneName(oldZoneName,newZoneName)
        ##### Redraw all
        self.updateAllGraph()
    # ------------------------------------------------------ replaceZoneWithTree
    def replaceZoneWithTree(self,t,oldBaseName,oldZoneName,newBaseName,newZoneName):
        # Delete old zone from data
        self.deleteZoneFromData(oldZoneName,oldBaseName+'/')
        # Find the newZone in tree
        self.addZoneWithTree(t,newZoneName,newBaseName)
        ##### Edit curves : change the zonename
        self.updateAllCurvesZoneName(oldBaseName+'/'+oldZoneName,newBaseName+'/'+newZoneName)
        ##### Redraw all
        self.updateAllGraph()
    # -------------------------------------------------- updateAllCurvesZoneName
    def updateAllCurvesZoneName(self,oldZoneName,newZoneName):
        for graph in self.graphWdwL:
            for ax_name in graph.fig.subGraph:
                graph.updateCurvesZoneName(ax_name,oldZoneName,newZoneName)

    # -------------------------------------------------------------- replaceData
    def setData(self,data):
        old_zones = []
        if self.data is not None:
            old_zones = self.data.keys()
        self.data = {}
        if isinstance(data,list):
            # set data according to a tree
            self.setDataWithTree(data)
        elif isinstance(data,dict):
            # set data according to a dict data
            self.setDataWithDict(data)
        # Clean curves
        for graph in self.graphWdwL:
            for ax_name in graph.fig.subGraph:
                for zonename in old_zones:
                    if zonename not in self.data:
                        graph.removeCurvesZoneName(ax_name,zonename)
        #
        if self.editCurveWdw  is not None:
            self.editCurveWdw.updateData()
        ##### Redraw all
        self.updateAllGraph()
#    # ------------------------------------------------------------------ setData
#    def setData(self,data):
#        print("self data...")
#        self.data = {}
#        if isinstance(data,list):
#            # set data according to a tree
#            self.setDataWithTree(data)
#        elif isinstance(data,dict):
#            # set data according to a dict data
#            self.setDataWithDict(data)
#        print('-'*50)
#        for k in self.data.keys():
#            print(k)
#            for kk in self.data[k].keys():
#                print('-> ',kk)
#        print('-'*50)
#        if self.editCurveWdw  is not None:
#            print('/'*50)
#            self.editCurveWdw.reloadWindow()
#            print('='*50)
    # ---------------------------------------------------------- setDataWithDict
    def setDataWithDict(self,d):
        tmp = {}
        # Check if d structure is 'zone' oriented
        isZoneOriented = True
        for k in d:
            if not isinstance(d[k],dict): # then it is not zone oriented
                isZoneOriented = False
                break
        if not isZoneOriented:
            tmp[default_base]=d
            self.setDataWithDict(tmp)
            return
        # Here the dict of data is zone oriented
        for k in self.data:
            tmp[k]=self.data[k]
        for k in d:
            tmp[k]=d[k]
        # Order dict
        self.data = OrderedDict(sorted(tmp.items(),key=lambda t : t[0]))
    # ---------------------------------------------------------- setDataWithTree
    def setDataWithTree(self,t):
        tmp = {}
        for base in Internal.getNodesFromType1(t, 'CGNSBase_t'):
            basename = Internal.getName(base)
            ## ## Loop on zones
            for zone in Internal.getNodesFromType1(base, 'Zone_t'):
                # Grab GridCoordinates
                zonename = Internal.getName(zone)
                ## ## Get GridCoorinates nodes
                try:
                    gridcoord = Internal.getNodesFromType2(zone,'GridCoordinates_t')[0]
                    for child in Internal.getChildren(gridcoord):
                        if not basename+'/'+zonename in tmp:
                            tmp[basename+'/'+zonename]={}
                        tmp[basename+'/'+zonename][Internal.getName(child)]=Internal.getValue(child)
                except IndexError: # No GridCoordinates node in this zone
                    # print('''gridcoord = Internal.getNodesFromType(zone,'GridCoordinates_t')[0] -----> Can not be loaded''')
                    pass

                # Grab FlowSolution for ZoneBC
                try:
                    zoneBC = Internal.getNodesFromType1(zone,'ZoneBC_t')[0]
                    for bc in Internal.getChildren(zoneBC):
                        bcname = Internal.getName(bc)
                        try:
                            bcdata = Internal.getNodesFromType(zoneBC,'BCData_t')[0]
                            for var in Internal.getChildren(bcdata):
                                if Internal.getType(var) == 'DataArray_t':
                                    if not basename+'/'+zonename in tmp:
                                        tmp[basename+'/'+zonename]={}
                                    tmp[basename+'/'+zonename][Internal.getName(var)+'@'+bcname]=Internal.getValue(var)
                        except IndexError:
                            # print('''bcdata = Internal.getNodesFromType(zoneBC,'BCData_t')[0] ------> Can not be loaded''')
                            pass

                except IndexError: # No zoneBC in this zone
                    # print('''zoneBC = Internal.getNodesFromType(zone,'ZoneBC_t')[0] -----> Can not be loaded''')
                    pass
                # Grab FlowSolution (the rest of them)
                for flowsolution in Internal.getNodesFromType1(zone,'FlowSolution_t'):
                    flowsolutionname = Internal.getName(flowsolution)
                    for var in Internal.getChildren(flowsolution):
                        if Internal.getType(var)=='DataArray_t':
                            if not basename+'/'+zonename in tmp:
                                tmp[basename+'/'+zonename]={}
                            tmp[basename+'/'+zonename][Internal.getName(var)+'@'+flowsolutionname]=Internal.getValue(var)


        self.data = OrderedDict(sorted(tmp.items(),key=lambda t : t[0]))
    # ----------------------------------------------------------------- setQueue
    def setQueue(self,queue):
        self.queue=queue
    # ---------------------------------------------------------------- setThread
    def setThread(self,thread):
        self.thread=thread
    # ------------------------------------------------------------- killProgramm
    def killProgramm(self):
        if self.thread:
            self.thread.endApplication()
        # # 1/- Close all graphs
        # self.cmd_closeGraph()
        # # 2/- Close all edit windows
        # for win in [self.addCurveWdw,self.editCurveWdw,self.editGridWdw,self.editLegendWdw,self.editAxisWdw,self.editGraphWdw] :
        #     if win:
        #         win.cmd_close()
        # 2/- Close mainwindow
        self.parent.quit()
    # --------------------------------------------------------------- initialize
    def initialize(self):
        #
        self.graphName2Id = {}
        self.graphWdwL = []
        self.addCurveWdw   = None
        self.editCurveWdw  = None
        self.editGridWdw   = None
        self.editLegendWdw = None
        self.editAxisWdw   = None
        self.editGraphWdw  = None
        #

        # # Configure grid for postionning for the main window
        self.grid_columnconfigure(0,weight=1)
        #
        self.grid_rowconfigure(0,weight=1)
        self.grid_rowconfigure(1,weight=8)
        self.grid_rowconfigure(2,weight=1)

        ################## Row 0 of self : create graph ########################

        # Add a labelFrame
        lblframeCreate = TTK.LabelFrame(self, text="Create Graph")
        lblframeCreate.grid(row=0,column=0,sticky="nsew")
        # # Configure grid for postionning for the inside of the LabelFrame
        lblframeCreate.grid_columnconfigure(0,weight=1)
        lblframeCreate.grid_columnconfigure(1,weight=1)
        #
        lblframeCreate.grid_rowconfigure(0,weight=1)

        # Add butons to the label frame

        # ROW 1
        B = TTK.Button(lblframeCreate,text='Add Graph',command=self.cmd_addGraph)
        B.grid(row=0,column=0,sticky="nsew")
        #
        frame = TTK.Frame(lblframeCreate)
        frame.grid(row=0,column=1,sticky="NSEW")
        frame.grid_columnconfigure(0,weight=1)
        frame.grid_columnconfigure(1,weight=1)
        frame.grid_rowconfigure(0,weight=1)
        frame.grid_rowconfigure(1,weight=1)
        #
        lblframe = TTK.LabelFrame(frame,text="Name")
#        lblframe = TK.LabelFrame(lblframeCreate,text="Name")
        lblframe.grid(row=0,column=0,sticky="NSEW")
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        #
        self.entryAddName = TK.Entry(lblframe)
        self.entryAddName.insert(TK.END,"Graph #%s"%len(self.graphWdwL))
        self.entryAddName.bind("<Button-1>",self.eraseEntry)
        self.entryAddName.bind("<FocusOut>",self.checkEntry)
        self.entryAddName.grid(row=0,column=0,sticky="NSEW")
        self.entryAddName.name = 'name'
        #
        lblframe = TTK.LabelFrame(frame,text="Zones")
#        lblframe = TK.LabelFrame(lblframeCreate,text="Zones")
        lblframe.grid(row=0,column=1,sticky="NSEW")
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        #
        self.entryAddGraph = TK.Entry(lblframe)
        self.entryAddGraph.insert(TK.END, "1:1") #TODO : comprendre pourquoi le END
        self.entryAddGraph.bind("<Button-1>",self.eraseEntry)
        self.entryAddGraph.bind("<FocusOut>",self.checkEntry)
        self.entryAddGraph.grid(row=0,column=0,sticky="nsew")
        self.entryAddGraph.name = 'zones'
        #
        lblframe = TTK.LabelFrame(frame,text="Fig. Size (inches)")
#        lblframe = TK.LabelFrame(lblframeCreate,text="Size")
        lblframe.grid(row=1,column=0,sticky="NSEW")
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        #
        self.entrySize = TK.Entry(lblframe)
        self.entrySize.insert(TK.END, "default") #TODO : comprendre pourquoi le END
        self.entrySize.bind("<Button-1>",self.eraseEntry)
        self.entrySize.bind("<FocusOut>",self.checkEntry)
        self.entrySize.grid(row=0,column=0,sticky="nsew")
        self.entrySize.name = 'size'
        #
        lblframe = TTK.LabelFrame(frame,text="dpi (dots per inch)")
#        lblframe = TK.LabelFrame(lblframeCreate,text="Resolution")
        lblframe.grid(row=1,column=1,sticky="NSEW")
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        #
        self.entryResolution = TK.Entry(lblframe)
        self.entryResolution.insert(TK.END, "default") #TODO : comprendre pourquoi le END
        self.entryResolution.bind("<Button-1>",self.eraseEntry)
        self.entryResolution.bind("<FocusOut>",self.checkEntry)
        self.entryResolution.grid(row=0,column=0,sticky="nsew")
        self.entryResolution.name = 'resolution'


        ################## Row 1 of self : Edit graph ##########################
        # Add a labelFrame
        lblframeEdit = TTK.LabelFrame(self, text="Edit Graph")
        lblframeEdit.grid(row=1,column=0,sticky="nsew")
        # # Configure grid for postionning for the inside of the LabelFrame
        lblframeEdit.grid_columnconfigure(0,weight=1)
        lblframeEdit.grid_columnconfigure(1,weight=1)
        #
        lblframeEdit.grid_rowconfigure(0,weight=1)
        lblframeEdit.grid_rowconfigure(1,weight=1)
        lblframeEdit.grid_rowconfigure(2,weight=1)
        lblframeEdit.grid_rowconfigure(3,weight=1)
        lblframeEdit.grid_rowconfigure(4,weight=1)
        lblframeEdit.grid_rowconfigure(5,weight=1)
        lblframeEdit.grid_rowconfigure(6,weight=1)
        lblframeEdit.grid_rowconfigure(7,weight=1)

        #
        ## Select Graph
        #
        lblframe = TTK.LabelFrame(lblframeEdit, text="Graph Name")
        lblframe.grid(row=0,column=0,sticky='NESW')
        # Grid of lblframe
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        #

        self.graphNameList = [g.name for g in self.graphWdwL]
        if len(self.graphNameList)==0:
            self.graphNameList.append('')
        self.activeGraph = cttk.Combobox(lblframe,values = self.graphNameList,state='readonly')
        try:
            self.updateGraphName2Id()
            self.activeGraph.val = self.graphName2Id[self.graphNameList[0]]
        except KeyError:
            self.activeGraph.val = ''
        self.activeGraph.set(self.activeGraph.val)
        self.activeGraph.bind("<<ComboboxSelected>>",self.cmd_graphNameChange)
        self.activeGraph.grid(row=0,column=0,sticky="NSEW")

        #
        ## Select Zone -> Row 0 of lblframeEdit
        #
        lblframe = TTK.LabelFrame(lblframeEdit, text="Zone")
        lblframe.grid(row=0,column=1,sticky='NESW')
        # Grid of lblframe
        lblframe.grid_columnconfigure(0,weight=1)
        lblframe.grid_rowconfigure(0,weight=1)
        #
        try:
            self.positionList = self.graphWdwL[self.activeGraph.val].subGraph.keys()
        except IndexError:
            self.positionList = ['']
        except TypeError:
            self.positionList = ['']
        self.position = cttk.Combobox(lblframe,values = self.positionList,state='readonly')
        self.position.val = self.positionList[0]
        self.position.set(self.position.val)
        self.position.bind("<<ComboboxSelected>>",self.cmd_positionChange)
        self.position.grid(row=0,column=0,sticky="NSEW")

        #
        ## Edit/Delete curves -> Row 1 of lblframeEdit
        #
        B = TTK.Button(lblframeEdit,text='Add/Edit/Delete Curves',command=self.cmd_editCurves)
        B.grid(row=1,column=0,columnspan=1,sticky="nsew")
        #
        ## Set Axis -> Row 2 of lblframeEdit
        #
        B = TTK.Button(lblframeEdit,text='Set Axis',command=self.cmd_setAxis)
        B.grid(row=2,column=0,columnspan=1,sticky="nsew")
        #
        ## Set Legend -> Row 3 of lblframeEdit
        #
        B = TTK.Button(lblframeEdit,text='Set Legend',command=self.cmd_setLegend)
        B.grid(row=1,column=1,columnspan=1,sticky="nsew")
        #
        ## Set Grids -> Row 4 of lblframeEdit
        #
        B = TTK.Button(lblframeEdit,text='Set Grid',command=self.cmd_setGrid)
        B.grid(row=2,column=1,columnspan=1,sticky="nsew")
        #
        ## Set Graph Settings -> Row 5 of lblframeEdit
        #
        B = TTK.Button(lblframeEdit,text='Graph settings',command=self.cmd_setGraph)
        B.grid(row=3,column=0,columnspan=1,sticky="nsew")
        #
        ## Export -> Row 6 of lblframeEdit
        #
        B = TTK.Button(lblframeEdit,text='Export',command=self.cmd_export)
        B.grid(row=3,column=1,columnspan=1,sticky="nsew")
        #
        ## Close Graph -> Row 7 of lblframeEdit
        #
        B = TTK.Button(lblframeEdit,text='Close all graphs',command=self.cmd_closeAllGraph)
        B.grid(row=7,column=0,columnspan=1,sticky="nsew")
        #
        B = TTK.Button(lblframeEdit,text='Close graph',command=self.cmd_closeGraph)
        B.grid(row=7,column=1,columnspan=1,sticky="nsew")

        ################## Row 2 of self : Graph configuration #################
        # Add a labelFrame
        lblframeEdit = TTK.LabelFrame(self, text="Graph configurations management")
        lblframeEdit.grid(row=2,column=0,sticky="nsew")
        # # Configure grid for postionning for the inside of the LabelFrame
        lblframeEdit.grid_columnconfigure(0,weight=1)
        lblframeEdit.grid_columnconfigure(1,weight=1)
        #
        lblframeEdit.grid_rowconfigure(0,weight=1)
        #
        ## Col 0
        #
        B = TTK.Button(lblframeEdit,text='Save',command=self.cmd_confSave)
        B.grid(row=0,column=0,columnspan=1,sticky="nsew")
        #
        B = TTK.Button(lblframeEdit,text='Load',command=self.cmd_confLoad)
        B.grid(row=0,column=1,columnspan=1,sticky="nsew")

    # ----------------------------------------------------- selectPositionByName
    def selectPositionByName(self,name):
        ind = self.positionList.index(name)
        # Avoid updating if clicking occured on the same ax as previous click
        if (self.positionList[ind]==self.position.get()):
            return
        # Update
        self.position.val = self.positionList[ind]
        self.position.set(self.position.val)
        #
        # Update editCurvesWindow
        if self.editCurveWdw is not None:
            self.editCurveWdw.graph = self.activeGraph.val
            self.editCurveWdw.zone  = self.position.val
            self.editCurveWdw.updateBar(self)
        else:
            self.editCurveWdw.initialize(self)
        # We can not merge the two if statements since initialize may set self.eidtCurveWdw to None
        for w in [self.editGridWdw,self.editLegendWdw,self.editAxisWdw]:
            if w is not None:
                w.initialize(self)
            if w is not None:
                w.reloadWindow()
    # -------------------------------------------------------- selectGraphByName
    def selectGraphByName(self,name):
        ind = self.graphNameList.index(name)
        self.updateGraphName2Id()
        self.activeGraph.val = self.graphName2Id[self.graphNameList[ind]]
        self.activeGraph.set(self.graphNameList[ind])
        self.updatePositionList()
    # --------------------------------------------------------------- checkEntry
    def checkEntry(self,event):
        widget = event.widget
        val = widget.get()
        if widget.name == 'resolution':
            pattern = '^[0123456789]*$'
            default = '^default$'
            if not re.match(pattern,val) and not re.match(default,val):
                val=''
#            for i in val:
#                if i not in '0123456789':
#                    val = ''
#                    break
        elif widget.name == 'size':
            pattern = '^\([0123456789]*,[0123456789]*\)$'
            default = '^default$'
            if not re.match(pattern,val) and not re.match(default,val):
                val = ''
        elif widget.name == 'zones':
            pattern = '^[0123456789]*:[0123456789]*$'
            if not re.match(pattern,val):
                val = ''
        if val=='':
            widget.delete(0,TK.END)
            widget.insert(TK.END,widget.lastValue)
    # --------------------------------------------------------------- eraseEntry
    def eraseEntry(self,event):
        widget = event.widget
        widget.lastValue = widget.get()
        widget.delete(0,TK.END)
    # ------------------------------------------------------------- cmd_confSave
    def cmd_confSave(self):
        global STYLEFILE
        # Works only with python 2, for python 3, it seems that the module name has changed to "filedialog"
        filename = tkFileDialog.asksaveasfilename(parent=self, initialdir=os.getcwd(), initialfile=STYLEFILE, filetypes=[('python', ".py")])
        if filename=='' or filename is None: return
        STYLEFILE = filename
        self.confSave(filename,True)
    # ----------------------------------------------------------------- confSave
    def confSave(self,filename,withUI):
        space = ' '*4
        #
        ## 1 /- Check that the directory where we save exists
        #
        directory = os.path.split(filename)[0]
        if not os.path.isdir(directory):
            if withUI:
                tkMessageBox.showerror("Saving configuration","Failed because %s does not exist"%directory)
                return
            else:
                print('''### Error: Saving configuration failed because %s does not exist'''%directory)
                return
        #
        ## 2/- Create the dictionary to save
        #
        instance = open(filename,'w')
        lines  = '''# coding: utf-8\n'''
        lines += '''from __future__ import unicode_literals\n'''
        lines += '''from tkPlotXY import *\n'''
        lines += '''def loadVisu(obj):\n'''
        #
        indgraph = 0
        indcurve = 0

        for graph in self.graphWdwL:
            # Commented lines
            lines += '''%s'''%space+'#'*35+'\n'
            lines += '''%s# Graph %s\n'''%(space,indgraph)
            lines += '''%s'''%space+'#'*25+'\n'
            # Create the graph
            lines += '''    %s\n'''%graph.write(indgraph)
            # Suplotparams
            lines += '''%s'''%graph.subPlotParams.write(indgraph)
            # TightLayout
            lines += '''%s'''%graph.tightLayout.write(indgraph)

            # Get the figure object
            figure = graph.fig
            # Loop on subgraph
            indsubgraph = 0
            for k in figure.subGraph:
                indgrid = 0
                indaxis = 0
                subgraph = figure.subGraph[k]
                # Commented lines
                lines += '''%s'''%space+'#'*25+'\n'
                lines += '''%s# SubGraph '%s'\n'''%(space,subgraph.name)
                lines += '''%s'''%space+'#'*15+'\n'
                #
                # Create all axes
                lines += '''%s######\n%s# Axis\n%s######\n'''%(space,space,space)
                for ind_axis in range(len(subgraph.axis)):
                    if subgraph.axis[ind_axis].type[0]=='main':
                        lines +='''    axis_%s_%s_%s = graph_%s.getAxis('%s',%s)\n'''%(indgraph,indsubgraph,ind_axis,indgraph,subgraph.name,ind_axis)
                        # lines +='''    indAxis_%s_%s_%s = 0\n'''%(indgraph,indsubgraph,ind_axis)
                    elif subgraph.axis[ind_axis].type[0]=='twinx':
                        lines +='''    axis_%s_%s_%s = graph_%s.addAxis('%s',shared='x',axis=axis_%s_%s_%s)\n'''%(indgraph,indsubgraph,
                                                                                                                    ind_axis,indgraph,subgraph.name,
                                                                                                                    indgraph,indsubgraph,
                                                                                                                    subgraph.axis[ind_axis].type[1])
                    elif subgraph.axis[ind_axis].type[0]=='twiny':
                        lines +='''    axis_%s_%s_%s = graph_%s.addAxis('%s',shared='y',axis=axis_%s_%s_%s)\n'''%(indgraph,indsubgraph,
                                                                                                                    ind_axis,indgraph,subgraph.name,
                                                                                                                    indgraph,indsubgraph,
                                                                                                                    subgraph.axis[ind_axis].type[1])
                    elif subgraph.axis[ind_axis].type[0] == 'new':
                        lines +='''    axis_%s_%s_%s = graph_%s.addAxis('%s')\n'''%(indgraph,indsubgraph,
                                                                                    ind_axis,indgraph,subgraph.name)
                # Loop on curves
                lines += '''%s######\n%s# Curves\n%s######\n'''%(space,space,space)
                for curve in subgraph.curves:
                    lines += '''%s\n'''%(curve.write(indgraph,subgraph.name,indsubgraph,indcurve))
                    indcurve += 1
                # Get Grid
                lines += '''%s######\n%s# Grids\n%s######\n'''%(space,space,space)
                for ind_axis in range(len(subgraph.axis)):
                    gridobj = graph.getGrid(subgraph.name, ind_axis)
                    lines += '''    grid_%s_%s_%s = graph_%s.getGrid('%s',axis=axis_%s_%s_%s)\n'''%(indgraph,indsubgraph,ind_axis,indgraph,subgraph.name,indgraph,indsubgraph,ind_axis)
                    lines += '''%s'''%(gridobj.write('grid_%s_%s_%s'%(indgraph,indsubgraph,ind_axis)))

                # Get Axis
                lines += '''%s######\n%s# Axis settings\n%s######\n'''%(space,space,space)
                for ind_axis in range(len(subgraph.axis)):
                    axisobj = graph.getAxis(subgraph.name, ind_axis)
                    # Following line is useless now that addAxis() returns an axis !
                    # lines += '''    axis_%s_%s_%s = graph_%s.getAxis('%s',%s)\n'''%(indgraph,indsubgraph,ind_axis,indgraph,subgraph.name,ind_axis)
                    lines += '''%s'''%(axisobj.write('axis_%s_%s_%s'%(indgraph,indsubgraph,ind_axis)))

                # Get Legend
                lines += '''%s######\n%s# Legends\n%s######\n'''%(space,space,space)
                legendobj = graph.getLegend(subgraph.name)
                lines += '''    legend_%s_%s = graph_%s.getLegend('%s')\n'''%(indgraph,indsubgraph,indgraph,subgraph.name)
                lines += '''%s'''%(legendobj.write('legend_%s_%s'%(indgraph,indsubgraph)))
                #
                indsubgraph += 1
            indgraph += 1

        #
        instance.write(lines)
        instance.close()
    # ------------------------------------------------------------- cmd_confLoad
    def cmd_confLoad(self):
        global STYLEFILE
        # Works only with python 2, for python 3, it seems that the module name has changed to "filedialog"
        filename = tkFileDialog.askopenfilename(parent=self, initialdir=os.getcwd(), initialfile=STYLEFILE, filetypes=[('python', ".py")])
        if filename=='' or filename is None: return
        STYLEFILE = filename
        self.confLoad(filename,True)
    # ----------------------------------------------------------------- confLoad
    def confLoad(self,filename,withUI):
        #
        ## 1/- Check that the file exists
        #
        if not os.path.isfile(filename):
            if withUI:
                tkMessageBox.showerror("Loading configuration","Failed because %s does not exist"%filename)
                return
            else:
                print('''### Error: Loading configuration failed because %s does not exist'''%filename)
                return
        #
        ## 2/- Load data
        #
#        cwd = os.getcwd()
#        os.chdir(os.path.split(filename)[0])
#        modulename = os.path.splitext(os.path.split(filename)[1])[0]
#        print(modulename)
#        loadedModule = __import__(modulename)
#        os.chdir(cwd)
#        loadedModule.loadData(self)
#        self.updateAllGraph()
#

        # WARNING : THE FOLLOWING IS ONLY WORKING FOR PYTHON 2.X
        cwd = os.getcwd()
        modulename = os.path.splitext(os.path.split(filename)[1])[0]
        loadedModule = imp.load_source(modulename, filename)
        loadedModule.loadVisu(self)
        self.updateAllGraph()

        # FOR PYTHON 3.X<3.4 TRY SOMETHING LIKE :
#        from importlib.machinery import SourceFileLoader
#        foo = SourceFileLoader("module.name", "/path/to/file.py").load_module()

        # FOR PYTHON 3.X>3.4 TRY SOMETHING LIKE :
#        import importlib.util
#        spec = importlib.util.spec_from_file_location("module.name", "/path/to/file.py")
#        foo = importlib.util.module_from_spec(spec)
#        spec.loader.exec_module(foo)

    # ------------------------------------------------------ cmd_graphNameChange
    def cmd_graphNameChange(self,event):
        widget = event.widget
        val = widget.get()
        # The two following lines seem to be not necessary
        pos = self.graphNameList.index(val)
        self.updateGraphName2Id()
        widget.val = self.graphName2Id[val]
        widget.set(val)
        self.updatePositionList()
    # ------------------------------------------------------- cmd_positionChange
    def cmd_positionChange(self,event):
        widget = event.widget
        val = widget.get()
        # The two following lines seem to be not necessary
        pos = self.positionList.index(val)
        widget.val = self.positionList[pos]
        self.position.set(widget.val)
        #
        # Update editCurvesWindow
        if self.editCurveWdw is not None:
            self.editCurveWdw.graph = self.activeGraph.val
            self.editCurveWdw.zone  = self.position.val
            self.editCurveWdw.updateBar(self)
        else:
            self.editCurveWdw.initialize(self)
        # We can not merge the two if statements since initialize may set self.eidtCurveWdw to None
        for w in [self.editGridWdw,self.editLegendWdw,self.editAxisWdw]:
            if w is not None:
                w.initialize(self)
            if w is not None:
                w.reloadWindow()
    # -------------------------------------------------------- updateactiveGraph
    def updateactiveGraph(self):
        self.graphNameList = [g.name for g in self.graphWdwL]
        self.updateGraphName2Id()
        if len(self.graphNameList)==0:
            self.graphNameList.append('')
            self.activeGraph['values']=self.graphNameList
            self.activeGraph.val = ''
            self.activeGraph.set('')
            self.updatePositionList()
        else:
            self.activeGraph['values']=self.graphNameList
            self.activeGraph.val = self.graphName2Id[self.graphNameList[0]]
            self.activeGraph.set(self.graphNameList[0])
            self.updatePositionList()
    # ------------------------------------------------------- updateGraphName2Id
    def updateGraphName2Id(self):
        self.graphName2Id = {}
        for ind in range(len(self.graphWdwL)):
            graph = self.graphWdwL[ind]
            self.graphName2Id[graph.name]=ind

    # ------------------------------------------------------- updatePositionList
    def updatePositionList(self):
        # Changes the values available for position according to the graph id selected
        graphId = self.activeGraph.val
        try:
            self.positionList = self.graphWdwL[graphId].fig.subGraph.keys()
        except IndexError:
            self.positionList = ['']
        except TypeError:
            self.positionList = ['']
        self.position['values']=self.positionList
        self.position.val = self.positionList[0]
        self.position.set(self.position.val)
        # Update editCurvesWindow
        if self.editCurveWdw is not None:
            self.editCurveWdw.graph = self.activeGraph.val
            self.editCurveWdw.zone  = self.position.val
            self.editCurveWdw.updateBar(self)
        # else:
        #     print('''no, shouldn't we do something ?''')
            # self.editCurveWdw.initialize(self)
        #
        # We can not merge the two if statements since initialize may set self.editCurveWdw to None
        # for w in [self.editCurveWdw,self.editGridWdw,self.editLegendWdw,self.editAxisWdw,self.editGraphWdw,]:
        for w in [self.editGridWdw,self.editLegendWdw,self.editAxisWdw,self.editGraphWdw,]:
            if w is not None:
                w.initialize(self)
        for w in [self.editGridWdw,self.editLegendWdw,self.editAxisWdw,self.editGraphWdw,]:
            if w is not None:
                w.reloadWindow()
    # ------------------------------------------------------ cmd_editCurvesGraph
    def cmd_editCurvesGraph(self,var,widget):
        try:
            widget.val = int(var.get())
            if widget.val>=len(self.graphWdwL):
                widget.val=len(self.graphWdwL)-1
                var.set(widget.val)
            if widget.val==-1:
                var.set('')
                widget.val=''
        except ValueError:
            var.set(filterInteger(var.get()))
            widget.val = var.get()
        self.updatePositionList()
    # --------------------------------------------------------------- cmd_export
    def cmd_export(self):
        # Get path to save
        global EXPORTFILE
        # Works only with python 2, for python 3, it seems that the module name has changed to "filedialog"
        filename = tkFileDialog.asksaveasfilename(parent=self, initialdir=os.getcwd(), initialfile=EXPORTFILE, filetypes=[('png', ".png")])
        if filename=='' or filename is None: return
        EXPORTFILE = filename
        self.export(filename)

    def export(self, filename):
        # Get active graph
        try:
            val = self.activeGraph.val
            if val<len(self.graphWdwL):
                # Save active graph
                self.graphWdwL[val].save(filename)
                print('Wrote file %s.'%filename)
                return True
        except ValueError: return False

    def getActiveGraphFigSize(self):
        try:
            val = self.activeGraph.val
            if val<len(self.graphWdwL):
                # Save active graph
                f = self.graphWdwL[val].figsize
                if f is None: return (130,100)
                else: return None
            else: return (130,100)
        except:
            return (130,100)

    # ----------------------------------------------------------- cmd_closeGraph
    def cmd_closeGraph(self):
        try:
            val = self.activeGraph.val
            if val<len(self.graphWdwL):
                self.graphWdwL[val].destroy()
                del self.graphWdwL[val]
            self.renumberGraph()

        except ValueError: pass
        self.updateactiveGraph()
    # ------------------------------------------------------------ renumberGraph
    def renumberGraph(self):
        for ind in range(len(self.graphWdwL)):
            graph = self.graphWdwL[ind]
            graph.index = ind
    # -------------------------------------------------------- cmd_closeAllGraph
    def cmd_closeAllGraph(self):
        for graph in self.graphWdwL:
            graph.destroy()
        self.graphWdwL=[]
        self.editCurvesGraphSV.set('')
    # ------------------------------------------------------------- cmd_addGraph
    def cmd_addGraph(self):
        if CTK is not None:
            if CTK.t == []: return

        if self.entryResolution.get() == 'default': dpi = None
        else: dpi = int(self.entryResolution.get())

        sizeStr = self.entrySize.get()
        if sizeStr == 'default':
            figsize = None
        else:
            sizeStrStrip = sizeStr[1:-1]
            figsize = (int(sizeStrStrip.split(',')[0]),int(sizeStrStrip.split(',')[1]))
        self.createGraph(self.entryAddName.get(),self.entryAddGraph.get(),dpi,figsize)
        #
        if len(self.graphWdwL)==1:
            self.cmd_editCurves()
    # -------------------------------------------------------------- createGraph
    def createGraph(self,name,conf,dpi=None,figsize=None):
        new_graph = GraphTK(self,name,conf,dpi,figsize)
        self.lift()
        self.focus()
    # ------------------------------------------------------------- cmd_setGraph
    def cmd_setGraph(self):
        if not self.editGraphWdw:
            self.editGraphWdw = editGraphWindow()
            self.editGraphWdw.initialize(self)
        else:
            self.editGraphWdw.lift()
            self.editGraphWdw.focus()
    # ----------------------------------------------------------- cmd_editCurves
    def cmd_editCurves(self):
        if not self.editCurveWdw:
            self.editCurveWdw = editCurvesWindow()
            self.editCurveWdw.initialize(self)
        else:
            self.editCurveWdw.lift()
            self.editCurveWdw.focus()
    # -------------------------------------------------------------- cmd_setAxis
    def cmd_setAxis(self):
        if not self.editAxisWdw:
            self.editAxisWdw = editAxisWindow()
            self.editAxisWdw.initialize(self)
        else:
            self.editAxisWdw.lift()
            self.editAxisWdw.focus()
    # ------------------------------------------------------------ cmd_setLegend
    def cmd_setLegend(self):
        if not self.editLegendWdw:
            self.editLegendWdw = editLegendWindow()
            self.editLegendWdw.initialize(self)
        else:
            self.editLegendWdw.lift()
            self.editLegendWdw.focus()
    # -------------------------------------------------------------- cmd_setGrid
    def cmd_setGrid(self):
        if not self.editGridWdw:
            self.editGridWdw = editGridWindow()
            self.editGridWdw.initialize(self)
        else:
            self.editGridWdw.lift()
            self.editGridWdw.focus()
    # ---------------------------------------------------------- processIncoming
    def processIncoming(self):
        """
        Handle all the messages currently in the queue (if any).
        """
        # Asynchrone treatment of the queue ...
        while self.queue.qsize():
            try:
                self.data = self.queue.get(0)
                # Check contents of message and do what it says
                # As a test, we simply print it
                print("From client thread, msg = ")
                self.updateAllGraph()
            except Queue.Empty:
                pass
    # ----------------------------------------------------------- updateAllGraph
    def updateAllGraph(self):
        for graph in self.graphWdwL:
            for ax_name in graph.fig.subGraph:
                graph.updateGraph(ax_name)
        self.addFrameAllMovie()
    # -------------------------------------------------------------- addAllMovie
    def addAllMovie(self):
        for graph in self.graphWdwL:
            graph.fig.addMovie(graph.name)
    # --------------------------------------------------------- addFrameAllMovie
    def addFrameAllMovie(self):
        for graph in self.graphWdwL:
            graph.fig.addFrameMovie()
    # ------------------------------------------------------------ closeAllMovie
    def closeAllMovie(self):
        for graph in self.graphWdwL:
            graph.fig.closeMovie()


# ==============================================================================
# ==============================================================================
class DesktopTK(TK.Tk):
    # --------------------------------------------------------------------- init
    def __init__(self,parent):
        TK.Tk.__init__(self,parent)
        self.protocol(name="WM_DELETE_WINDOW", func=self.killProgramm)
        self.title(string="PlotXY by Cassiopee")
        self.desktopFrameTK = DesktopFrameTK(self)
        self.grid_columnconfigure(0,weight=1)
        self.grid_rowconfigure(0,weight=1)
        self.desktopFrameTK.grid(row=0,column=0,sticky='NSEW')
    # ------------------------------------------------------ replaceGroupZonesWithDict
    def replaceGroupZonesWithDict(self,d,oldZoneList):
        self.desktopFrameTK.replaceGroupZonesWithDict(d,oldZoneList)
    # -------------------------------------------------------------- replaceGroupZones
    def replaceGroupZones(self,data,oldZoneList):
        self.desktopFrameTK.replaceGroupZones(d,oldZoneList)
    # ------------------------------------------------------ replaceGroupZonesWithTree
    def replaceGroupZonesWithTree(self,d,oldZoneList):
        self.desktopFrameTK.replaceGroupZonesWithTree(d,oldZoneList)
    # -------------------------------------------------------------- updateGroupCurves
    def updateGroupCurves(self,oldZoneList,newZoneList):
        self.desktopFrameTK.updateGroupCurves(oldZoneList,newZoneList)
    # ---------------------------------------------------------- addDataWithTree
    def addDataWithTree(self,t):
        self.desktopFrameTK.addDataWithTree(t)
    # -------------------------------------------------------------- deleteZoneInCurve
    def deleteZoneInCurve(self,zoneName):
        self.desktopFrameTK.deleteZoneInCurve(zoneName)
    # ------------------------------------------------------------------ addZone
    def addZone(self,data,zoneName,baseName='.*'):
        self.desktopFrameTK.addZone(data,zoneName,baseName)
    # ---------------------------------------------------------- addZoneWithTree
    def addZoneWithTree(self,t,zoneName,baseName='.*'):
        self.desktopFrameTK.addZoneWithTree(t,zoneName,baseName)
    # ---------------------------------------------------------- addZoneWithDict
    def addZoneWithDict(self,d,zoneName):
        self.desktopFrameTK.addZoneWithDict(d,zoneName)
    # ------------------------------------------------------- deleteZoneFromData
    def deleteZoneFromData(self,zoneName,oldBaseName=""):
        self.desktopFrameTK.deleteZoneFromData(zoneName,oldBaseName)
    # -------------------------------------------------------------- replaceZone
    def replaceZone(self,data,oldZoneName,newZoneName,oldBaseName="",newBaseName=""):
        self.desktopFrameTK.replaceZone(data,oldZoneName,newZoneName,oldBaseName,newBaseName)
    # ------------------------------------------------------ replaceZoneWithDict
    def replaceZoneWithDict(self,d,oldZoneName,newZoneName):
        self.desktopFrameTK.replaceZoneWithDict(d,oldZoneName,newZoneName)
    # ------------------------------------------------------ replaceZoneWithTree
    def replaceZoneWithTree(self,t,oldBaseName,oldZoneName,newBaseName,newZoneName):
        self.desktopFrameTK.replaceZoneWithTree(t,oldBaseName,oldZoneName,newBaseName,newZoneName)
    # ------------------------------------------------------------------ setData
    def setData(self,data):
        self.desktopFrameTK.setData(data)
    # ----------------------------------------------------------------- setQueue
    def setQueue(self,queue):
        self.desktopFrameTK.setQueue(queue)
    # ---------------------------------------------------------------- setThread
    def setThread(self,thread):
        self.desktopFrameTK.setThread(thread)
    # ------------------------------------------------------------- killProgramm
    def killProgramm(self):
        self.desktopFrameTK.killProgramm()
        self.quit()
    # ----------------------------------------------------------------- confSave
    def confSave(self,filename,withUI):
        self.desktopFrameTK.confSave(filename,withUI)
    # ----------------------------------------------------------------- confLoad
    def confLoad(self,filename,withUI):
        self.desktopFrameTK.confLoad(filename,withUI)
    # -------------------------------------------------------------- createGraph
    def createGraph(self,name,conf,dpi=None,figsize=None):
        self.desktopFrameTK.createGraph(name,conf,dpi,figsize)
    # ---------------------------------------------------------- processIncoming
    def processIncoming(self):
        self.desktopFrameTK.processIncoming()
    # ----------------------------------------------------------- updateAllGraph
    def updateAllGraph(self):
        self.desktopFrameTK.updateAllGraph()
    # -------------------------------------------------------------- addAllMovie
    def addAllMovie(self):
        self.desktopFrameTK.addAllMovie()
    # --------------------------------------------------------- addFrameAllMovie
    def addFrameAllMovie(self):
        self.desktopFrameTK.addFrameAllMovie()
    # ------------------------------------------------------------ closeAllMovie
    def closeAllMovie(self):
        self.desktopFrameTK.closeAllMovie()
# ==============================================================================
# ==============================================================================
class Desktop():
    """
    An object of class Desktop allows you to create all your graphs.
    """
    # --------------------------------------------------------------------- init
    def __init__(self):
        self.initialize()
#        self.data = data
        self.data = None
        self.thread = None

    # ------------------------------------------------------------------ addZone
    def addZone(self,data,zoneName,baseName='.*'):
        """
        Add a specific zone to the set of data. If pyTree format is used as input, then 'basename' can be specified to add only the zone from a specific base. If 'basename' is not specified in case of pyTree format, then the specified zone for all the bases from the pyTree will be added.
        """
        if self.data is None:
            self.data = {}
        if isinstance(data,list):
            # add zone according to a tree
            self.addZoneWithTree(data,zoneName,baseName)
        elif isinstance(data,dict):
            # add zone according to a dict data
            self.addZoneWithDict(data,zoneName)
    # ---------------------------------------------------------- addZoneWithTree
    def addZoneWithTree(self,t,zoneName,baseNameFilter):
        for base in Internal.getNodesFromType1(t, 'CGNSBase_t'):
            basename = Internal.getName(base)
            if re.match(baseNameFilter,basename):
                zone = Internal.getNodesFromName1(base, zoneName)[0]
                zonename = zoneName

                ## Get GridCoordinates nodes
                try:
                    gridcoord = Internal.getNodesFromType2(zone,'GridCoordinates_t')[0]
                    for child in Internal.getChildren(gridcoord):
                        if not basename+'/'+zonename in self.data:
                            self.data[basename+'/'+zonename]={}
                        self.data[basename+'/'+zonename][Internal.getName(child)]=Internal.getValue(child)
                except IndexError: # No GridCoordinates node in this zone
                    # print('''gridcoord = Internal.getNodesFromType(zone,'GridCoordinates_t')[0] -----> Can not be loaded''')
                    pass

                # Grab FlowSolution for ZoneBC
                try:
                    zoneBC = Internal.getNodesFromType1(zone,'ZoneBC_t')[0]
                    for bc in Internal.getChildren(zoneBC):
                        bcname = Internal.getName(bc)
                        try:
                            bcdata = Internal.getNodesFromType(zoneBC,'BCData_t')[0]
                            for var in Internal.getChildren(bcdata):
                                if Internal.getType(var) == 'DataArray_t':
                                    if not basename+'/'+zonename in self.data:
                                        self.data[basename+'/'+zonename]={}
                                    self.data[basename+'/'+zonename][Internal.getName(var)+'@'+bcname]=Internal.getValue(var)
                        except IndexError:
                            # print('''bcdata = Internal.getNodesFromType(zoneBC,'BCData_t')[0] ------> Can not be loaded''')
                            pass

                except IndexError: # No zoneBC in this zone
                    # print('''zoneBC = Internal.getNodesFromType(zone,'ZoneBC_t')[0] -----> Can not be loaded''')
                    pass
                # Grab FlowSolution (the rest of them)
                for flowsolution in Internal.getNodesFromType1(zone,'FlowSolution_t'):
                    flowsolutionname = Internal.getName(flowsolution)
                    for var in Internal.getChildren(flowsolution):
                        if Internal.getType(var)=='DataArray_t':
                            if not basename+'/'+zonename in self.data:
                                self.data[basename+'/'+zonename]={}
                            self.data[basename+'/'+zonename][Internal.getName(var)+'@'+flowsolutionname]=Internal.getValue(var)
    # ---------------------------------------------------------- addZoneWithDict
    def addZoneWithDict(self,d,zoneName):
        if zoneName in d:
            self.data[zoneName]=d[zoneName]
        else:
            print('''### Can not find zone %s in submitted data.'''%newZoneName)
    # ------------------------------------------------------- deleteZoneFromData
    def deleteZoneFromData(self,zoneName,oldBaseName=""):
        """
        Simply delete data from a given zone and base to the set of data from the Desktop object.
        """
        for k in self.data:
            re_str = oldBaseName+zoneName.replace('\\','\\\\') # replace \ by \\ for regular expression conversion
            if re.match(re_str,k):
                del self.data[k]
        # Shouldn't we add a deleteZoneFromCurve ???????
        # Answer : no it has to be done after deletezonefromdata only if we don't replace the zone ...

    # -------------------------------------------------------------- replaceGroupZones
    def replaceGroupZones(self,data,oldZoneList):
        if isinstance(data,list):
            # replace zone according to a tree
            self.replaceGroupZonesWithTree(data,oldZoneList)
        elif isinstance(data,dict):
            # Conformize the input dictionnary
            tmp = {}
            # Check if d structure is 'zone' oriented
            isZoneOriented = True
            for k in data:
                if not isinstance(data[k],dict): # then it is not zone oriented
                    isZoneOriented = False
                    break
            if not isZoneOriented:
                tmp[default_base]=data
            else:
                tmp = data
            # replace zone according to a dict data
            self.replaceGroupZonesWithDict(data,oldZoneList)
    # ------------------------------------------------------ replaceGroupZonesWithTree
    def replaceGroupZonesWithTree(self,d,oldZoneList):
        # add data and determine the list of new data
        newZoneList = self.addDataWithTree(d)
        # Get the curves that are concerned by a group of old zones and change it to the group of new zones
        self.updateGroupCurves(oldZoneList,newZoneList)
        # Compare old zones to new zones group and remove old zones that do not exist anymore
        for zoneName in oldZoneList:
            if zoneName not in newZoneList:
                # Delete old zone from data
                self.deleteZoneFromData(zoneName)
                # Delete old zone from curve
                self.deleteZoneInCurve(zoneName)
        ##### Redraw all
        self.updateAllGraph()
        #
    # ---------------------------------------------------------- addDataWithTree
    def addDataWithTree(self,t):
        tmp = self.data
        newZoneList = []
        for base in Internal.getNodesFromType1(t, 'CGNSBase_t'):
            basename = Internal.getName(base)
            ## ## Loop on zones
            for zone in Internal.getNodesFromType1(base,'Zone_t'):
                # Grab GridCoordinates
                zonename = Internal.getName(zone)
                ## ## Get GridCoorinates nodes
                try:
                    gridcoord = Internal.getNodesFromType1(zone,'GridCoordinates_t')[0]
                    for child in Internal.getChildren(gridcoord):
                        if not basename+'/'+zonename in tmp:
                            tmp[basename+'/'+zonename]={}
                        tmp[basename+'/'+zonename][Internal.getName(child)]=Internal.getValue(child)
                        newZoneList.append(basename+'/'+zonename)
                except IndexError: # No GridCoorinates node in this zone
                    # print('''gridcoord = Internal.getNodesFromType(zone,'GridCoordinates_t')[0] -----> Can not be loaded''')
                    pass

                # Grab FlowSolution for ZoneBC
                try:
                    zoneBC = Internal.getNodesFromType1(zone,'ZoneBC_t')[0]
                    for bc in Internal.getChildren(zoneBC):
                        bcname = Internal.getName(bc)
                        try:
                            bcdata = Internal.getNodesFromType(zoneBC,'BCData_t')[0]
                            for var in Internal.getChildren(bcdata):
                                if Internal.getType(var) == 'DataArray_t':
                                    if not basename+'/'+zonename in tmp:
                                        tmp[basename+'/'+zonename]={}
                                    tmp[basename+'/'+zonename][Internal.getName(var)+'@'+bcname]=Internal.getValue(var)
                                    newZoneList.append(basename+'/'+zonename)
                        except IndexError:
                            # print('''bcdata = Internal.getNodesFromType(zoneBC,'BCData_t')[0] ------> Can not be loaded''')
                            pass

                except IndexError: # No zoneBC in this zone
                    # print('''zoneBC = Internal.getNodesFromType(zone,'ZoneBC_t')[0] -----> Can not be loaded''')
                    pass
                # Grab FlowSolution (the rest of them)
                for flowsolution in Internal.getNodesFromType1(zone,'FlowSolution_t'):
                    flowsolutionname = Internal.getName(flowsolution)
                    for var in Internal.getChildren(flowsolution):
                        if Internal.getType(var)=='DataArray_t':
                            if not basename+'/'+zonename in tmp:
                                tmp[basename+'/'+zonename]={}
                            tmp[basename+'/'+zonename][Internal.getName(var)+'@'+flowsolutionname]=Internal.getValue(var)
                            newZoneList.append(basename+'/'+zonename)


        self.data = OrderedDict(sorted(tmp.items(),key=lambda t : t[0]))
        return newZoneList
    # ------------------------------------------------------ replaceGroupZonesWithDict
    def replaceGroupZonesWithDict(self,d,oldZoneList):
        # Add new data and determine the list of new zones
        newZoneList = []
        for zoneName in d:
            newZoneList.append(zoneName)
            self.addZoneWithDict(d,zoneName)
        # Get the curves that are concerned by a group of old zones and change it to the group of new zones
        self.updateGroupCurves(oldZoneList,newZoneList)
        # Compare old zones to new zones group and remove old zones that do not exist anymore
        for zoneName in oldZoneList:
            if zoneName not in newZoneList:
                # Delete old zone from data
                self.deleteZoneFromData(zoneName)
                # Delete old zone from curve
                self.deleteZoneInCurve(zoneName)
        ##### Redraw all
        self.updateAllGraph()
    # -------------------------------------------------------------- deleteZoneInCurve
    def deleteZoneInCurve(self,zoneName):
        for graph in self.graphWdwL:
            for ax_name in graph.fig.subGraph:
                graph.deleteZoneInCurve(ax_name,zoneName)
    # -------------------------------------------------------------- updateGroupCurves
    def updateGroupCurves(self,oldZoneList,newZoneList):
        for graph in self.graphWdwL:
            for ax_name in graph.fig.subGraph:
                graph.updateGroupCurvesZoneName(ax_name,oldZoneList,newZoneList)
    # -------------------------------------------------------------- replaceZone
    def replaceZone(self,data,oldZoneName,newZoneName,oldBaseName="",newBaseName=""):
        """
        Allows you to replace a given zone by a new one. If some data from the old zone are plotted, then the plot remains but the data are updated with the data loaded from the new zone. According to the method addZone, old basename and new basename can be given to specify the targeted base.
        """
        if isinstance(data,list):
            # replace zone according to a tree
            self.replaceZoneWithTree(data,oldBaseName,oldZoneName,newBaseName,newZoneName)
        elif isinstance(data,dict):
            if oldBaseName != "": oldZoneName = oldBaseName+'/'+oldZoneName
            if newBaseName != "": newZoneName = newBaseName+'/'+newZoneName
            # replace zone according to a dict data
            self.replaceZoneWithDict(data,oldZoneName,newZoneName)
    # ------------------------------------------------------ replaceZoneWithDict
    def replaceZoneWithDict(self,d,oldZoneName,newZoneName):
        # Delete old zone from data
        self.deleteZoneFromData(oldZoneName)
        #
        self.addZoneWithDict(d,newZoneName)
        #
        ##### Edit curves : change the zonename
        self.updateAllCurvesZoneName(oldZoneName,newZoneName)
        ##### Redraw all
        self.updateAllGraph()
    # ------------------------------------------------------ replaceZoneWithTree
    def replaceZoneWithTree(self,t,oldBaseName,oldZoneName,newBaseName,newZoneName):
        # Delete old zone from data
        self.deleteZoneFromData(oldZoneName,oldBaseName+'/')
        # Find the newZone in tree
        self.addZoneWithTree(t,newZoneName,newBaseName)
        ##### Edit curves : change the zonename
        self.updateAllCurvesZoneName(oldBaseName+'/'+oldZoneName,newBaseName+'/'+newZoneName)
        ##### Redraw all
        self.updateAllGraph()
    # -------------------------------------------------- updateAllCurvesZoneName
    def updateAllCurvesZoneName(self,oldZoneName,newZoneName):
        for graph in self.graphWdwL:
            for ax_name in graph.fig.subGraph:
                graph.updateCurvesZoneName(ax_name,oldZoneName,newZoneName)

#    # ------------------------------------------------------------------ setData
#    def setData(self,data):
#        self.data = {}
#        if isinstance(data,list):
#            # set data according to a tree
#            self.setDataWithTree(data)
#        elif isinstance(data,dict):
#            # set data according to a dict data
#            self.setDataWithDict(data)
    # ------------------------------------------------------------------ setData
    def setData(self,data):
        """
        Set all the data from the input pyTree or dictionnary to the Desktop data manager.
        """
        old_zones = []
        if self.data is not None:
            old_zones = self.data.keys()
        self.data = {}
        if isinstance(data,list):
            # set data according to a tree
            self.setDataWithTree(data)
        elif isinstance(data,dict):
            # set data according to a dict data
            self.setDataWithDict(data)
        # Clean curves
        for graph in self.graphWdwL:
            for ax_name in graph.fig.subGraph:
                for zonename in old_zones:
                    if zonename not in self.data:
                        graph.removeCurvesZoneName(ax_name,zonename)
        #
        if self.editCurveWdw  is not None:
            self.editCurveWdw.updateData()
        ##### Redraw all
        self.updateAllGraph()
    # ---------------------------------------------------------- setDataWithDict
    def setDataWithDict(self,d):
        tmp = {}
        # Check if d structure is 'zone' oriented
        isZoneOriented = True
        for k in d:
            if not isinstance(d[k],dict): # then it is not zone oriented
                isZoneOriented = False
                break
        if not isZoneOriented:
            tmp[default_base]=d
            self.setDataWithDict(tmp)
            return
        # Here the dict of data is zone oriented
        for k in self.data:
            tmp[k]=self.data[k]
        for k in d:
            tmp[k]=d[k]
        # Order dict
        self.data = OrderedDict(sorted(tmp.items(),key=lambda t : t[0]))

    # ---------------------------------------------------------- setDataWithTree
    def setDataWithTree(self,t):
        tmp = {}
#        base =  Internal.getNodesFromType1(t, 'CGNSBase_t')[0] # Assume there is only one base
        for base in Internal.getNodesFromType1(t, 'CGNSBase_t'):
            basename = Internal.getName(base)
            ## ## Loop on zones
            for zone in Internal.getNodesFromType1(base,'Zone_t'):
                # Grab GridCoordinates
                zonename = Internal.getName(zone)
                ## ## Get GridCoorinates nodes
                try:
                    gridcoord = Internal.getNodesFromType1(zone,'GridCoordinates_t')[0]
                    for child in Internal.getChildren(gridcoord):
                        if not basename+'/'+zonename in tmp:
                            tmp[basename+'/'+zonename]={}
                        tmp[basename+'/'+zonename][Internal.getName(child)]=Internal.getValue(child)
                except IndexError: # No GridCoorinates node in this zone
                    # print('''gridcoord = Internal.getNodesFromType(zone,'GridCoordinates_t')[0] -----> Can not be loaded''')
                    pass

                # Grab FlowSolution for ZoneBC
                try:
                    zoneBC = Internal.getNodesFromType1(zone,'ZoneBC_t')[0]
                    for bc in Internal.getChildren(zoneBC):
                        bcname = Internal.getName(bc)
                        try:
                            bcdata = Internal.getNodesFromType(zoneBC,'BCData_t')[0]
                            for var in Internal.getChildren(bcdata):
                                if Internal.getType(var) == 'DataArray_t':
                                    if not basename+'/'+zonename in tmp:
                                        tmp[basename+'/'+zonename]={}
                                    tmp[basename+'/'+zonename][Internal.getName(var)+'@'+bcname]=Internal.getValue(var)
                        except IndexError:
                            # print('''bcdata = Internal.getNodesFromType(zoneBC,'BCData_t')[0] ------> Can not be loaded''')
                            pass

                except IndexError: # No zoneBC in this zone
                    # print('''zoneBC = Internal.getNodesFromType(zone,'ZoneBC_t')[0] -----> Can not be loaded''')
                    pass
                # Grab FlowSolution (the rest of them)
                for flowsolution in Internal.getNodesFromType1(zone,'FlowSolution_t'):
                    flowsolutionname = Internal.getName(flowsolution)
                    for var in Internal.getChildren(flowsolution):
                        if Internal.getType(var)=='DataArray_t':
                            if not basename+'/'+zonename in tmp:
                                tmp[basename+'/'+zonename]={}
                            tmp[basename+'/'+zonename][Internal.getName(var)+'@'+flowsolutionname]=Internal.getValue(var)


        self.data = OrderedDict(sorted(tmp.items(),key=lambda t : t[0]))
    # ----------------------------------------------------------------- setQueue
    def setQueue(self,queue):
        self.queue=queue
    # ---------------------------------------------------------------- setThread
    def setThread(self,thread):
        self.thread=thread
    # ------------------------------------------------------------- killProgramm
    def killProgramm(self):
        if self.thread:
            self.thread.endApplication()
        self.quit()
    # --------------------------------------------------------------- initialize
    def initialize(self):
        #
        self.graphName2Id = {}
        self.graphWdwL = []
        self.addCurveWdw   = None
        self.editCurveWdw  = None
        self.editGridWdw   = None
        self.editLegendWdw = None
        self.editAxisWdw   = None

    # ----------------------------------------------------------------- confSave
    def confSave(self,filename):
        #
        ## 1 /- Check that the directory where we save exists
        #
        directory = os.path.split(filename)[0]
        if not os.path.isdir(directory):
            print('''### Error: Saving configuration failed because %s does not exist'''%directory)
            return
        #
        ## 2/- Create the dictionary to save
        #
        instance = open(filename,'w')
        lines  = '''# coding: utf-8\n'''
        lines += '''from __future__ import unicode_literals\n'''
        lines += '''from tkPlotXY import *\n'''
        lines += '''def loadVisu(obj):\n'''
        #
        indgraph = 0
        indcurve = 0

        for graph in self.graphWdwL:
            # Create the graph
            lines += '''    %s\n'''%graph.write(indgraph)
            # Suplotparams
            lines += '''%s'''%graph.subPlotParams.write(indgraph)
            # TightLayout
            lines += '''%s'''%graph.tightLayout.write(indgraph)
            # Get the figure object
            figure = graph.fig
            # Loop on subgraph
            indsubgraph = 0
            for k in figure.subGraph:
                indgrid = 0
                indaxis = 0
                subgraph = figure.subGraph[k]
                # Create all axes
                for ind_axis in range(len(subgraph.axis)):
                    if subgraph.axis[ind_axis].type[0]=='main':
                        lines +='''    indAxis_%s_%s_%s = 0\n'''%(indgraph,indsubgraph,ind_axis)
                    elif subgraph.axis[ind_axis].type[0]=='twinx':
                        lines +='''    indAxis_%s_%s_%s = graph_%s.addAxis('%s',shared='x',axis=indAxis_%s_%s_%s)\n'''%(indgraph,indsubgraph,
                                                                                                                    ind_axis,indgraph,subgraph.name,
                                                                                                                    indgraph,indsubgraph,
                                                                                                                    subgraph.axis[ind_axis].type[1])
                    elif subgraph.axis[ind_axis].type[0]=='twiny':
                        lines +='''    indAxis_%s_%s_%s = graph_%s.addAxis('%s',shared='y',axis=indAxis_%s_%s_%s)\n'''%(indgraph,indsubgraph,
                                                                                                                    ind_axis,indgraph,subgraph.name,
                                                                                                                    indgraph,indsubgraph,
                                                                                                                    subgraph.axis[ind_axis].type[1])
                    elif subgraph.axis[ind_axis].type[0] == 'new':
                        lines +='''    indAxis_%s_%s_%s = graph_%s.addAxis('%s')\n'''%(indgraph,indsubgraph,
                                                                                    ind_axis,indgraph,subgraph.name)
                # Loop on curves
                for curve in subgraph.curves:
                    lines += '''%s\n'''%(curve.write(indgraph,subgraph.name,indcurve))
                    indcurve += 1
                # Get Grid
                for ind_axis in range(len(subgraph.axis)):
                    gridobj = graph.getGrid(subgraph.name, ind_axis)
                    lines += '''    grid_%s_%s_%s = graph_%s.getGrid('%s',%s)\n'''%(indgraph,indsubgraph,ind_axis,indgraph,subgraph.name,ind_axis)
                    lines += '''%s'''%(gridobj.write('grid_%s_%s_%s'%(indgraph,indsubgraph,ind_axis)))

                # Get Axis
                for ind_axis in range(len(subgraph.axis)):
                    axisobj = graph.getAxis(subgraph.name, ind_axis)
                    lines += '''    axis_%s_%s_%s = graph_%s.getAxis('%s',%s)\n'''%(indgraph,indsubgraph,ind_axis,indgraph,subgraph.name,ind_axis)
                    lines += '''%s'''%(axisobj.write('axis_%s_%s_%s'%(indgraph,indsubgraph,ind_axis)))

                # Get Legend
                legendobj = graph.getLegend(subgraph.name)
                lines += '''    legend_%s_%s = graph_%s.getLegend('%s')\n'''%(indgraph,indsubgraph,indgraph,subgraph.name)
                lines += '''%s'''%(legendobj.write('legend_%s_%s'%(indgraph,indsubgraph)))
                #
                indsubgraph += 1



            indgraph += 1

        #
        instance.write(lines)
        instance.close()

    # ----------------------------------------------------------------- confLoad
    def confLoad(self,filename):
        #
        ## 1/- Check that the file exists
        #)
        if not os.path.isfile(filename):
            print('''### Error: Loading configuration failed because %s does not exist'''%filename)
            return
        #
        ## 2/- Load data
        #
#        cwd = os.getcwd()
#        os.chdir(os.path.split(filename)[0])
#        modulename = os.path.splitext(os.path.split(filename)[1])[0]
#        print(modulename)
#        loadedModule = __import__(modulename)
#        os.chdir(cwd)
#        loadedModule.loadData(self)
#        self.updateAllGraph()
#

        # WARNING : THE FOLLOWING IS ONLY WORKING FOR PYTHON 2.X
        cwd = os.getcwd()
        modulename = os.path.splitext(os.path.split(filename)[1])[0]
        loadedModule = imp.load_source(modulename, filename)
        loadedModule.loadVisu(self)
        self.updateAllGraph()

        # FOR PYTHON 3.X<3.4 TRY SOMETHING LIKE :
#        from importlib.machinery import SourceFileLoader
#        foo = SourceFileLoader("module.name", "/path/to/file.py").load_module()

        # FOR PYTHON 3.X>3.4 TRY SOMETHING LIKE :
#        import importlib.util
#        spec = importlib.util.spec_from_file_location("module.name", "/path/to/file.py")
#        foo = importlib.util.module_from_spec(spec)
#        spec.loader.exec_module(foo)

#    # -------------------------------------------------------- updateactiveGraph
#    def updateactiveGraph(self):
#        self.graphNameList = [g.name for g in self.graphWdwL]
#        self.updateGraphName2Id()
#        if len(self.graphNameList)==0:
#            self.graphNameList.append('')
#            self.activeGraph['values']=self.graphNameList
#            self.activeGraph.val = ''
#            self.activeGraph.set('')
#            self.updatePositionList()
#        else:
#            self.activeGraph['values']=self.graphNameList
#            self.activeGraph.val = self.graphName2Id[self.graphNameList[0]]
#            self.activeGraph.set(self.graphNameList[0])
#            self.updatePositionList()
    # ------------------------------------------------------- updateGraphName2Id
    def updateGraphName2Id(self):
        self.graphName2Id = {}
        for ind in range(len(self.graphWdwL)):
            graph = self.graphWdwL[ind]
            self.graphName2Id[graph.name]=ind

    # ----------------------------------------------------------- cmd_closeGraph
    def cmd_closeGraph(self):
        try:
            val = self.activeGraph.val
            if val<len(self.graphWdwL):
                self.graphWdwL[val].destroy()
                del self.graphWdwL[val]
            self.renumberGraph()

        except ValueError:
            pass
        self.updateactiveGraph()
    # ------------------------------------------------------------ renumberGraph
    def renumberGraph(self):
        for ind in range(len(self.graphWdwL)):
            graph = self.graphWdwL[ind]
            graph.index = ind
    # -------------------------------------------------------- cmd_closeAllGraph
    def cmd_closeAllGraph(self):
        for graph in self.graphWdwL:
            graph.destroy()
        self.graphWdwL=[]
        self.editCurvesGraphSV.set('')
    # -------------------------------------------------------------- createGraph
    def createGraph(self,name,conf,dpi=None,figsize=None):
        """
        Create a window where the plots will be drawn. A matricial description is used to define this window. For instance, here are described some settings for 'conf' variable:

        * '1:1' : a single plot in this graph window

        * '2:2' : 4 plots (2 rows and 2 columns)

        * '2:1' : 2 plots (2 rows and 1 column)

        * '1:2' : 2 plots (1 row and 2 columns)

        figsize and dpi can configured to adapt the size and the resolution of the graph window if needed. You can for instance use figsize=(12,3) to enlarge your image.
        """
        new_graph = Graph(self,name,conf,dpi,figsize)
        return new_graph

    # ---------------------------------------------------------- processIncoming
    def processIncoming(self):
        """
        Handle all the messages currently in the queue (if any).
        """
        # Asynchrone treatment of the queue ...
        while self.queue.qsize():
            try:
                self.data = self.queue.get(0)
                # Check contents of message and do what it says
                # As a test, we simply print it
                print("From client thread, msg = ")
                self.updateAllGraph()
            except Queue.Empty:
                pass
    # ----------------------------------------------------------- updateAllGraph
    def updateAllGraph(self):
        for graph in self.graphWdwL:
            for ax_name in graph.fig.subGraph:
                graph.updateGraph(ax_name)
        self.addFrameAllMovie()
    # -------------------------------------------------------------- addAllMovie
    def addAllMovie(self):
        for graph in self.graphWdwL:
            graph.fig.addMovie(graph.name)
    # --------------------------------------------------------- addFrameAllMovie
    def addFrameAllMovie(self):
        for graph in self.graphWdwL:
            graph.fig.addFrameMovie()
    # ------------------------------------------------------------ closeAllMovie
    def closeAllMovie(self):
        for graph in self.graphWdwL:
            graph.fig.closeMovie()

# ==============================================================================
# ==============================================================================
class MatplotlibFigure():
    def __init__(self,parent,graph,conf,dpi=None,figsize=None):
        self.graph = graph
        # Root window
        self.parent=parent
        # dpi, figsize
        self.dpi = dpi
        self.figsize = figsize
        # Store configuration
        self.conf = conf
        # init movie
        self.movie = None
        # Get nl and nc
        try:
            val = [int(v) for v in (self.conf).split(':')]
        except ValueError:
            self.cmd_close()
            return
        self.nl = val[0]
        self.nc = val[1]
        # Figure
        self.instance = plt.figure(facecolor="white",figsize=self.figsize,dpi=self.dpi)
        # print('''##### Fig instance = ''',self.instance)
        # Ax
        self.subGraph = {}
        for il in range(self.nl):
            for ic in range(self.nc):
                ind = il*self.nc+ic
                self.subGraph['%s:%s'%(il+1,ic+1)] = SubGraph(self.instance,self.nl,self.nc,ind,il,ic)

    # --------------------------------------------------------------- saveFigure
    def saveFigure(self,path,format=None):
        # self.instance.savefig(path,facecolor=self.instance.get_facecolor(),edgecolor='none',format=format)
        self.instance.savefig(path,format=format)
    # ------------------------------------------------------------------- getFig
    def getFig(self):
        return self.instance
    # ------------------------------------------------------------------ addGrid
    def getGrid(self,iCurSubGraph,ind=0):
        return self.subGraph[iCurSubGraph].grid_property[ind]
    # ---------------------------------------------------------------- getLegend
    def getLegend(self,iCurSubGraph):
        return self.subGraph[iCurSubGraph].legend_property
    # ------------------------------------------------------------------ getAxis
    def getAxis(self,iCurSubGraph,ind=0):
        return self.subGraph[iCurSubGraph].axis_property[ind]
    # ----------------------------------------------------------------- addMovie
    def addMovie(self,name):
        self.movie = Movie(self.instance,os.path.join(os.getcwd(),name+'.mp4'))
        self.movie.activate()
    def addFrameMovie(self):
        if self.movie:
            self.movie.write()
    def closeMovie(self):
        if self.movie:
            self.movie.exit()
            self.movie = None
    # --------------------------------------------------------------- drawFigure
    def drawFigure(self):
        for iCurSubGraph in self.subGraph:
            self.drawOneFigure(iCurSubGraph)
    # ------------------------------------------------------removeCurvesZoneName
    def removeCurvesZoneName(self,ax_name,zonename):
        for ind in range(len(self.subGraph[ax_name].curves)):
            c = self.subGraph[ax_name].curves[ind]
            if zonename in c.zone:
                index = c.zone.index(zonename)
                del c.zone[index]
        # delete empty curves
        for c in self.subGraph[ax_name].curves:
            ind = self.subGraph[ax_name].curves.index(c)
            if len(c.zone)==0:
                # Delete current curve
                del self.subGraph[ax_name].curves[ind]

    # ----------------------------------------------------- deleteZoneInCurve
    '''
    Determine if all the old zones are used for a given curve, if it is so, then
    we add all the new zone and remove all the old zones
    -> Keep all the zones that are not concerned by this old/new change
    '''
    def deleteZoneInCurve(self,iCurSubGraph,zoneName):
        for c in self.subGraph[iCurSubGraph].curves:
            if zoneName in c.zone:
                c.zone.remove(zoneName)
    # ----------------------------------------------------- updateGroupCurvesZoneName
    '''
    Determine if all the old zones are used for a given curve, if it is so, then
    we add all the new zone and remove all the old zones
    -> Keep all the zones that are not concerned by this old/new change
    '''
    def updateGroupCurvesZoneName(self,iCurSubGraph,oldZoneList,newZoneList):
        for c in self.subGraph[iCurSubGraph].curves:
            for zoneName in oldZoneList:
                if zoneName not in c.zone:
                    return
            for zoneName in newZoneList:
                if zoneName not in c.zone:
                    c.zone.append(zoneName)
    # ----------------------------------------------------- updateCurvesZoneName
    def updateCurvesZoneName(self,iCurSubGraph,oldZoneName,newZoneName):
        for c in self.subGraph[iCurSubGraph].curves:
            if oldZoneName in c.zone:
                index = c.zone.index(oldZoneName)
                c.zone[index]=newZoneName
    # --------------------------------------------------------------- drawOneFigure
    def drawOneFigure(self,iCurSubGraph):
        curves = self.subGraph[iCurSubGraph].curves
        iCurrentAxis = 0
        for iCurrentAxis in range(len(self.subGraph[iCurSubGraph].axis)):
            self.subGraph[iCurSubGraph].axis[iCurrentAxis].clear()
        ###
        legend_list = []
        legend_text = []
        ###
        for iCurrentAxis in range(len(self.subGraph[iCurSubGraph].axis)):
            plt.sca(self.subGraph[iCurSubGraph].axis[iCurrentAxis])
    #        self.subGraph[iCurSubGraph].axis[iCurrentAxis].clear()
#            plt.cla()
#            self.subGraph[iCurSubGraph].axis[iCurrentAxis].clear()
            self.subGraph[iCurSubGraph].axis[iCurrentAxis].name=self.subGraph[iCurSubGraph].name
    #        self.subGraph[iCurSubGraph].curves = curves

            xaxis_label = ""
            yaxis_label = ""
            for c in self.subGraph[iCurSubGraph].curves:
                if c.axis == iCurrentAxis:

                    if c.marker_sampling_step == '':
                        markingSampling = None
                    else:
                        if c.marker_sampling_start=='':
                            markingSampling = c.marker_sampling_step
                        else:
                            if c.marker_sampling_end=='':
                                markingSampling=(c.marker_sampling_start,c.marker_sampling_step)
                            else:
                                print('''Could not manage to make the 'slice' runs, therefore, it will be necessary to
                                filter data manually ...\n Note that for the moment, the slice becomes: (start,step)''')
        #                        markingSampling = slice(c.marker_sampling_start,c.marker_sampling_end,c.marker_sampling_step)
                                markingSampling=(c.marker_sampling_start,c.marker_sampling_step)
                    if c.visible:
                        if (c.varx is not None) and (c.vary is not None):
                            for zone in c.zone:
                                try:
                                    if self.parent.data[zone][c.varx].shape != self.parent.data[zone][c.vary].shape:
                #                        print('''### Message : can not plot %s vs %s since their shape are different.'''%(c.varx,c.vary))
                                        continue
                                except AttributeError:
                                    try:
                                        if len(self.parent.data[zone][c.varx])!=len(self.parent.data[zone][c.vary]):
                                            # print('#~'*50)
                                            continue
                                    except TypeError:
                                        continue
            #                    except KeyError:
            #                        continue
                                try:
                                    if self.parent.data[zone][c.varx].shape[0]<2:
                                        continue
                                except AttributeError:
                                    try:
                                        if len(self.parent.data[zone][c.varx])<2:
                                            continue
                                    except TypeError:
                                        continue


                                sortedData = self.sortData(self.parent.data[zone][c.varx],self.parent.data[zone][c.vary])
                                # Update xaxis_label & yaxis_label
                                ### xaxis_label
                                re_var = re.compile(c.varx.split('@')[0]+',')
                                if not re.search(re_var,xaxis_label):
                                    xaxis_label+=c.varx.split('@')[0]+','
                                ### yaxis_label
                                re_var = re.compile(c.vary.split('@')[0]+',')
                                if not re.search(re_var,yaxis_label):
                                    yaxis_label+=c.vary.split('@')[0]+','
                                # Create the plot of the curve
                                thisPlot, = self.subGraph[iCurSubGraph].axis[iCurrentAxis].plot(sortedData[0],sortedData[1],color=c.line_color,linestyle=c.line_style,linewidth=c.line_width,
                                                        marker=marker_dic[c.marker_style],markersize=c.marker_size,markerfacecolor=c.marker_face_color,
                                                        markeredgecolor=c.marker_edge_color,markeredgewidth=c.marker_edge_width,markevery=markingSampling)
        #                        thisPlot, = plt.plot(sortedData[0],sortedData[1],color=c.line_color,linestyle=c.line_style,linewidth=c.line_width,
        #                                                marker=marker_dic[c.marker_style],markersize=c.marker_size,markerfacecolor=c.marker_face_color,
        #                                                markeredgecolor=c.marker_edge_color,markeredgewidth=c.marker_edge_width,markevery=markingSampling)

            #                thisPlot, = plt.plot(self.parent.data[c.varx],self.parent.data[c.vary],color=c.line_color,linestyle=c.line_style,linewidth=c.line_width,
            #                marker=marker_dic[c.marker_style],markersize=c.marker_size,markerfacecolor=c.marker_face_color,markeredgecolor=c.marker_edge_color,markeredgewidth=c.marker_edge_width,markevery=markingSampling)

                            if c.legend_display and self.subGraph[iCurSubGraph].legend_property.legend_display:
                                legend_list.append(thisPlot)
                                legend_text.append(c.legend_label)



            ## ## Set Axis
            self.subGraph[iCurSubGraph].axis[iCurrentAxis].relim()
            ## logscale
            if self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].x.axis_logscale:
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].set_xscale("log")
            if self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].y.axis_logscale:
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].set_yscale("log")
            ## Xmin, xmax, ymin and ymax
            self.subGraph[iCurSubGraph].axis[iCurrentAxis].set_xlim((self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].x.axis_min,self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].x.axis_max))
            self.subGraph[iCurSubGraph].axis[iCurrentAxis].set_ylim((self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].y.axis_min,self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].y.axis_max))
            ## auto scale
            count = 0
            if self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].x.axis_autoscale:
                count = count + 1
            if self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].y.axis_autoscale:
                count = count + 2
            conversion = [None,'x','y','both']
            if count == 0:
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].autoscale(enable=False)
            else:
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].autoscale(enable=True,axis=conversion[count])
#            # set label
#            self.subGraph[iCurSubGraph].axis[iCurrentAxis].set_xlabel(r'%s'%self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].x.axis_label)
#            self.subGraph[iCurSubGraph].axis[iCurrentAxis].set_ylabel(r'%s'%self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].y.axis_label)

            ## Invert Axis
            if self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].x.axis_inverted:
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].invert_xaxis()
            if self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].y.axis_inverted:
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].invert_yaxis()


            ##Set Ticks
            xmin,xmax = self.subGraph[iCurSubGraph].axis[iCurrentAxis].get_xlim()
            ymin,ymax = self.subGraph[iCurSubGraph].axis[iCurrentAxis].get_ylim()
            #
            ntickMx = self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].major.x.grid_tick_number
            ntickMy = self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].major.y.grid_tick_number
            ntickmx = self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].minor.x.grid_tick_number
            ntickmy = self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].minor.y.grid_tick_number

            stepx = (xmax-xmin)/(float(ntickMx))
            majorx = np.arange(xmin,xmax,stepx)
            # minorx = np.arange(xmin,xmax,(xmax-xmin)/(float(ntickmx)))
            minorx = np.arange(xmin,xmax,stepx/(float(ntickmx)))
            stepy = (ymax-ymin)/(float(ntickMy))
            majory = np.arange(ymin,ymax,stepy)
            # minory = np.arange(ymin,ymax,(ymax-ymin)/(float(ntickmy)))
            minory = np.arange(ymin,ymax,stepy/(float(ntickmy)))

            # Get the ticks position for debug
            # locs = self.subGraph[iCurSubGraph].axis[iCurrentAxis].xaxis.get_ticklocs()
            #
            self.subGraph[iCurSubGraph].axis[iCurrentAxis].set_xticks(majorx)
            self.subGraph[iCurSubGraph].axis[iCurrentAxis].set_xticks(minorx,minor=True)
            self.subGraph[iCurSubGraph].axis[iCurrentAxis].set_yticks(majory)
            self.subGraph[iCurSubGraph].axis[iCurrentAxis].set_yticks(minory,minor=True)
            # Tick size
            self.subGraph[iCurSubGraph].axis[iCurrentAxis].tick_params(axis='x',labelsize=self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].major.x.grid_tick_size)
            # self.subGraph[iCurSubGraph].axis[iCurrentAxis].tick_params(axis='x',labelsize=self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].minor.x.grid_tick_size,minor=True)
            self.subGraph[iCurSubGraph].axis[iCurrentAxis].tick_params(axis='y',labelsize=self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].major.y.grid_tick_size)
            # self.subGraph[iCurSubGraph].axis[iCurrentAxis].tick_params(axis='y',labelsize=self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].minor.y.grid_tick_size,minor=True)









            ## Spine position
            spines = self.subGraph[iCurSubGraph].axis[iCurrentAxis].spines
            offset_x = float(self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].x.axis_offset)
            spines['top'].set_position(('outward',offset_x))
            spines['bottom'].set_position(('outward',offset_x))
            #
            offset_y = float(self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].y.axis_offset)
            spines['left'].set_position(('outward',offset_y))
            spines['right'].set_position(('outward',offset_y))

            ## make patch and spine invisible
            self.subGraph[iCurSubGraph].axis[iCurrentAxis].set_frame_on(True)
            self.subGraph[iCurSubGraph].axis[iCurrentAxis].patch.set_visible(True)# BLABLA
            if iCurrentAxis==0:
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].patch.set_facecolor(self.graph.subgraph_background_color)# BLABLA
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].patch.set_alpha(self.graph.subgraph_background_alpha)# BLABLA
            else:
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].patch.set_facecolor("none")# BLABLA

            for sp in self.subGraph[iCurSubGraph].axis[iCurrentAxis].spines.values():
                sp.set_visible(False)

            ##### X AXIS
            if self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].x.axis_position == 'both':
                # Spine
                spines['top'].set_visible(True)
#                spines['top'].set_color('r')
                spines['bottom'].set_visible(True)
#                spines['bottom'].set_color('r')
                # Tick
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].xaxis.set_ticks_position('bottom')
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].xaxis.set_ticks_position('both')
                # Label
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].xaxis.set_label_position('bottom')

            elif self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].x.axis_position == 'bottom':
                # Spine
                spines['bottom'].set_visible(True)
#                spines['bottom'].set_color('r')
                # Tick
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].xaxis.set_ticks_position('bottom')
                # Label
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].xaxis.set_label_position('bottom')

            elif self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].x.axis_position == 'top':
                # Spine
                spines['top'].set_visible(True)
#                spines['top'].set_color('r')
                # Tick
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].xaxis.set_ticks_position('top')
                # Label
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].xaxis.set_label_position('top')

            ##### Y AXIS
            if self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].y.axis_position == 'both':
                # Spine
                spines['left'].set_visible(True)
#                spines['left'].set_color('r')
                spines['right'].set_visible(True)
#                spines['right'].set_color('r')
                # Tick
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].yaxis.set_ticks_position('left')
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].yaxis.set_ticks_position('both')
                # Label
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].yaxis.set_label_position('left')

            elif self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].y.axis_position == 'left':
                # Spine
                spines['left'].set_visible(True)
#                spines['left'].set_color('r')
                # Tick
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].yaxis.set_ticks_position('left')
                # Label
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].yaxis.set_label_position('left')

            elif self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].y.axis_position == 'right':
                # Spine
                spines['right'].set_visible(True)
#                spines['right'].set_color('r')
                # Tick
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].yaxis.set_ticks_position('right')
                # Label
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].yaxis.set_label_position('right')


            ## Hide axis
            if not self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].x.axis_visible:
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].get_xaxis().set_visible(False)
                for loc in spines:
                    if loc in ['top','bottom']:
                        spines[loc].set_visible(False)
            else:
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].get_xaxis().set_visible(True)
            if not self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].y.axis_visible:
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].get_yaxis().set_visible(False)
                for loc in spines:
                    if loc in ['left','right']:
                        spines[loc].set_visible(False)
            else:
                self.subGraph[iCurSubGraph].axis[iCurrentAxis].get_yaxis().set_visible(True)




            ## set label
            if len(xaxis_label)>0:
                xaxis_label = xaxis_label[:-1]
            if len(yaxis_label)>0:
                yaxis_label = yaxis_label[:-1]
            if self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].x.axis_label:
                xaxis_label = self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].x.axis_label
            if self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].y.axis_label:
                yaxis_label = self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].y.axis_label
            self.subGraph[iCurSubGraph].axis[iCurrentAxis].set_xlabel(r'%s'%xaxis_label,
                                                                        fontsize=self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].x.axis_label_fontsize)
            self.subGraph[iCurSubGraph].axis[iCurrentAxis].set_ylabel(r'%s'%yaxis_label,
                                                                        fontsize=self.subGraph[iCurSubGraph].axis_property[iCurrentAxis].y.axis_label_fontsize)

            # self.subGraph[iCurSubGraph].axis[iCurrentAxis].grid(color='r', linestyle='-', linewidth=2)


            ### Set Grid
            ### ### major grid
            ### ### ---> X
            self.subGraph[iCurSubGraph].axis[iCurrentAxis].grid(color=self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].major.x.grid_color,
                                                linestyle=self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].major.x.grid_style,
                                                linewidth=self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].major.x.grid_width, which='major', axis='x')
            # self.subGraph[iCurSubGraph].axis[iCurrentAxis].grid(TrueFalseDic[self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].major.x.display], which='major', axis='x')
            self.subGraph[iCurSubGraph].axis[iCurrentAxis].grid(self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].major.x.display, which='major', axis='x')
            ### ### ---> Y
            self.subGraph[iCurSubGraph].axis[iCurrentAxis].grid(color=self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].major.y.grid_color,
                                 linestyle=self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].major.y.grid_style,
                                 linewidth=self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].major.y.grid_width, which='major', axis='y')
            # self.subGraph[iCurSubGraph].axis[iCurrentAxis].grid(TrueFalseDic[self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].major.y.display], which='major', axis='y')
            self.subGraph[iCurSubGraph].axis[iCurrentAxis].grid(self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].major.y.display, which='major', axis='y')
            ### ### minor grid
            ### ### ---> X
            self.subGraph[iCurSubGraph].axis[iCurrentAxis].grid(color=self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].minor.x.grid_color,
                                                 linestyle=self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].minor.x.grid_style,
                                                 linewidth=self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].minor.x.grid_width, which='minor', axis='x')
            # self.subGraph[iCurSubGraph].axis[iCurrentAxis].grid(TrueFalseDic[self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].minor.x.display], which='minor', axis='x')
            self.subGraph[iCurSubGraph].axis[iCurrentAxis].grid(self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].minor.x.display, which='minor', axis='x')
            ### ### ---> Y
            self.subGraph[iCurSubGraph].axis[iCurrentAxis].grid(color=self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].minor.y.grid_color,
                                                 linestyle=self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].minor.y.grid_style,
                                                 linewidth=self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].minor.y.grid_width, which='minor', axis='y')
            # self.subGraph[iCurSubGraph].axis[iCurrentAxis].grid(TrueFalseDic[self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].minor.y.display], which='minor', axis='y')
            self.subGraph[iCurSubGraph].axis[iCurrentAxis].grid(self.subGraph[iCurSubGraph].grid_property[iCurrentAxis].minor.y.display, which='minor', axis='y')

        #### LEGEND CONFIGURATION
        if self.subGraph[iCurSubGraph].legend_property.legend_display:
            legend_label_prop = {'weight':self.subGraph[iCurSubGraph].legend_property.legend_label_weight,'size':self.subGraph[iCurSubGraph].legend_property.legend_label_size,
                                'style':self.subGraph[iCurSubGraph].legend_property.legend_label_style}
            legend = self.subGraph[iCurSubGraph].axis[iCurrentAxis].legend(legend_list,legend_text,loc=self.subGraph[iCurSubGraph].legend_property.legend_position,
                                        ncol=self.subGraph[iCurSubGraph].legend_property.legend_ncol,
                                        prop = legend_label_prop) # Can add the frameon=True option to activate/desactivate the frame (border + background) for legend
        ## ## LEGEND LABEL
            for text in legend.get_texts():
                text.set_color(self.subGraph[iCurSubGraph].legend_property.legend_label_color)
        ## ## LEGEND TITLE
            legend_title_prop = {'weight':self.subGraph[iCurSubGraph].legend_property.legend_title_weight,'size':self.subGraph[iCurSubGraph].legend_property.legend_title_size,
                             'style':self.subGraph[iCurSubGraph].legend_property.legend_title_style}
            legend.set_title(self.subGraph[iCurSubGraph].legend_property.legend_title, prop=legend_title_prop)
            legend.get_title().set_color(self.subGraph[iCurSubGraph].legend_property.legend_title_color)
        ## ## LEGEND BORDER
            legend.get_frame().set_linewidth(self.subGraph[iCurSubGraph].legend_property.legend_border_width)
            legend.get_frame().set_edgecolor(self.subGraph[iCurSubGraph].legend_property.legend_border_color)
        ## ## LEGEND BACKGROUND
            if self.subGraph[iCurSubGraph].legend_property.legend_background_color_active:
                legend.get_frame().set_facecolor(self.subGraph[iCurSubGraph].legend_property.legend_background_color)
            else:
                legend.get_frame().set_facecolor('none')
        #
#        self.instance.tight_layout()


    # --------------------------------------------------------------------------
    def sortData(self,v1,v2):
        return (v1,v2)
#        return zip(*(sorted(zip(v1,v2), key=lambda d: d[0])))



# ==============================================================================
# ==============================================================================
class DirAxis():
    """
    DirAxis contains the settings of the X or Y axis. This settings can directly be accessed from the class Axis.
    """
    def __init__(self,axis_logscale,axis_autoscale,axis_min,axis_max,axis_label,axis_inverted,axis_visible,axis_position,axis_offset,axis_label_fontsize):
        self.axis_logscale=axis_logscale
        self.axis_autoscale=axis_autoscale
        self.axis_min=axis_min
        self.axis_max=axis_max
        self.axis_label=axis_label
        self.axis_inverted = axis_inverted
        self.axis_visible = axis_visible
        self.axis_position = axis_position
        self.axis_offset = axis_offset
        self.axis_label_fontsize = axis_label_fontsize

    def setValue(self,variable,value):
        if variable == 'axis_logscale':
            self.axis_logscale = value
        elif variable == 'axis_autoscale':
            self.axis_autoscale = value
        elif variable == 'axis_min':
            self.axis_min = value
        elif variable == 'axis_max':
            self.axis_max = value
        elif variable == 'axis_label':
            self.axis_label = value
        elif variable == 'axis_inverted':
            self.axis_inverted = value
        elif variable == 'axis_visible':
            self.axis_visible = value
        elif variable == 'axis_position':
            self.axis_position = value
        elif variable == 'axis_offset':
            self.axis_offset = value
        elif variable == 'axis_label_fontsize':
            self.axis_label_fontsize = value
# ==============================================================================
# ==============================================================================
class Axis():
    """
    An Axis object contains the X-DirAxis and the Y-DirAxis of a given plot inside a Graph object. Multiple axis are available for a single plot.
    """
    def __init__(self,ind=0,
                    axis_x_logscale = default_values['Axis']['axis_x_logscale'],
                    axis_y_logscale = default_values['Axis']['axis_y_logscale'],
                    axis_x_autoscale = default_values['Axis']['axis_x_autoscale'],
                    axis_y_autoscale = default_values['Axis']['axis_y_autoscale'],
                    axis_x_min = default_values['Axis']['axis_x_min'],
                    axis_x_max = default_values['Axis']['axis_x_max'],
                    axis_y_min = default_values['Axis']['axis_y_min'],
                    axis_y_max = default_values['Axis']['axis_y_max'],
                    axis_x_label = default_values['Axis']['axis_x_label'],
                    axis_y_label = default_values['Axis']['axis_y_label'],
                    axis_x_inverted = default_values['Axis']['axis_x_inverted'],
                    axis_y_inverted = default_values['Axis']['axis_y_inverted'],
                    axis_x_visible = default_values['Axis']['axis_x_visible'],
                    axis_y_visible = default_values['Axis']['axis_y_visible'],
                    axis_x_position = default_values['Axis']['axis_x_position'],
                    axis_y_position = default_values['Axis']['axis_y_position'],
                    axis_x_offset = default_values['Axis']['axis_x_offset'],
                    axis_y_offset = default_values['Axis']['axis_y_offset'],
                    axis_x_label_fontsize = default_values['Axis']['axis_x_label_fontsize'],
                    axis_y_label_fontsize = default_values['Axis']['axis_y_label_fontsize']):

        self.ind = ind

        self.x = DirAxis(axis_x_logscale,axis_x_autoscale,axis_x_min,axis_x_max,axis_x_label,
                        axis_x_inverted,axis_x_visible,axis_x_position,axis_x_offset,axis_x_label_fontsize)
        self.y = DirAxis(axis_y_logscale,axis_y_autoscale,axis_y_min,axis_y_max,axis_y_label,
                        axis_y_inverted,axis_y_visible,axis_y_position,axis_y_offset,axis_x_label_fontsize)

    def getInd(self):
        return self.ind

    def setValue(self,axis,variable,value):
        if axis == 'x':
            self.x.setValue(variable,value)
        else:
            self.y.setValue(variable,value)
    def write(self,axisobj):
        lines =""

        if self.x.axis_logscale != default_values['Axis']['axis_x_logscale']    :
            lines +='''    %s.setValue('x','axis_logscale',%s)\n'''%(axisobj,bool(self.x.axis_logscale))
        if self.x.axis_autoscale != default_values['Axis']['axis_x_autoscale']    :
            lines +='''    %s.setValue('x','axis_autoscale',%s)\n'''%(axisobj,bool(self.x.axis_autoscale))
        if float(self.x.axis_min) != default_values['Axis']['axis_x_min']    :
            lines +='''    %s.setValue('x','axis_min',%s)\n'''%(axisobj,self.x.axis_min)
        if float(self.x.axis_max) != default_values['Axis']['axis_x_max']    :
            lines +='''    %s.setValue('x','axis_max',%s)\n'''%(axisobj,self.x.axis_max)
        if self.x.axis_label != default_values['Axis']['axis_x_label']    :
            lines +='''    %s.setValue('x','axis_label','%s')\n'''%(axisobj,self.x.axis_label)
        if self.x.axis_inverted != default_values['Axis']['axis_x_inverted']    :
            lines +='''    %s.setValue('x','axis_inverted',%s)\n'''%(axisobj,bool(self.x.axis_inverted))
        if self.x.axis_visible != default_values['Axis']['axis_x_visible']    :
            lines +='''    %s.setValue('x','axis_visible',%s)\n'''%(axisobj,bool(self.x.axis_visible))
        if self.x.axis_position != default_values['Axis']['axis_x_position']    :
            lines +='''    %s.setValue('x','axis_position','%s')\n'''%(axisobj,self.x.axis_position)
        if float(self.x.axis_offset) != default_values['Axis']['axis_x_offset']    :
            lines +='''    %s.setValue('x','axis_offset','%s')\n'''%(axisobj,self.x.axis_offset)
        if float(self.x.axis_label_fontsize) != default_values['Axis']['axis_x_label_fontsize']    :
            lines +='''    %s.setValue('x','axis_label_fontsize','%s')\n'''%(axisobj,self.x.axis_label_fontsize)
        if self.y.axis_logscale != default_values['Axis']['axis_y_logscale']    :
            lines +='''    %s.setValue('y','axis_logscale',%s)\n'''%(axisobj,bool(self.y.axis_logscale))
        if self.y.axis_autoscale != default_values['Axis']['axis_y_autoscale']    :
            lines +='''    %s.setValue('y','axis_autoscale',%s)\n'''%(axisobj,bool(self.y.axis_autoscale))
        if float(self.y.axis_min) != default_values['Axis']['axis_y_min']    :
            lines +='''    %s.setValue('y','axis_min',%s)\n'''%(axisobj,self.y.axis_min)
        if float(self.y.axis_max) != default_values['Axis']['axis_y_max']    :
            lines +='''    %s.setValue('y','axis_max',%s)\n'''%(axisobj,self.y.axis_max)
        if self.y.axis_label != default_values['Axis']['axis_y_label'] :
            lines +='''    %s.setValue('y','axis_label','%s')\n'''%(axisobj,self.y.axis_label)
        if self.y.axis_inverted != default_values['Axis']['axis_y_inverted'] :
            lines +='''    %s.setValue('y','axis_inverted',%s)\n'''%(axisobj,bool(self.y.axis_inverted))
        if self.y.axis_visible != default_values['Axis']['axis_y_visible'] :
            lines +='''    %s.setValue('y','axis_visible',%s)\n'''%(axisobj,bool(self.y.axis_visible))
        if self.y.axis_position != default_values['Axis']['axis_y_position'] :
            lines +='''    %s.setValue('y','axis_position','%s')\n'''%(axisobj,self.y.axis_position)
        if float(self.y.axis_offset) != default_values['Axis']['axis_y_offset'] :
            lines +='''    %s.setValue('y','axis_offset','%s')\n'''%(axisobj,self.y.axis_offset)
        if float(self.y.axis_label_fontsize) != default_values['Axis']['axis_y_label_fontsize'] :
            lines +='''    %s.setValue('y','axis_label_fontsize','%s')\n'''%(axisobj,self.y.axis_label_fontsize)

        return lines
# ==============================================================================
# ==============================================================================
class Legend():
    """
    An object of class Legend configures the legend for a given plot inside a Graph window.
    """
    def __init__(self,
                    legend_display = default_values['Legend']['legend_display'],
                    legend_title = default_values['Legend']['legend_title'],
                    legend_border_width = default_values['Legend']['legend_border_width'],
                    legend_border_color = default_values['Legend']['legend_border_color'],
                    legend_background_color = default_values['Legend']['legend_background_color'],
                    legend_background_color_active = default_values['Legend']['legend_background_color_active'],
                    legend_position = default_values['Legend']['legend_position'],
                    legend_ncol = default_values['Legend']['legend_ncol'],
                    legend_label_weight = default_values['Legend']['legend_label_weight'],
                    legend_label_style = default_values['Legend']['legend_label_style'],
                    legend_label_size = default_values['Legend']['legend_label_size'],
                    legend_label_color = default_values['Legend']['legend_label_color'],
                    legend_title_weight = default_values['Legend']['legend_title_weight'],
                    legend_title_style = default_values['Legend']['legend_title_style'],
                    legend_title_size = default_values['Legend']['legend_title_size'],
                    legend_title_color = default_values['Legend']['legend_title_color']):

        self.legend_display                 = legend_display
        self.legend_title                   = legend_title
        self.legend_border_width            = legend_border_width
        self.legend_border_color            = legend_border_color
        self.legend_background_color        = legend_background_color
        self.legend_background_color_active = legend_background_color_active
        self.legend_position                = legend_position
        self.legend_ncol                    = legend_ncol
        self.legend_label_weight            = legend_label_weight
        self.legend_label_style             = legend_label_style
        self.legend_label_size              = legend_label_size
        self.legend_label_color             = legend_label_color
        self.legend_title_weight            = legend_title_weight
        self.legend_title_style             = legend_title_style
        self.legend_title_size              = legend_title_size
        self.legend_title_color             = legend_title_color

    def setValue(self,variable,value):
        if variable == 'legend_display':
            self.legend_display = value
        elif variable == 'legend_title':
            self.legend_title = value
        elif variable == 'legend_border_width':
            self.legend_border_width = value
        elif variable == 'legend_border_color':
            self.legend_border_color = value
        elif variable == 'legend_background_color':
            self.legend_background_color = value
        elif variable == 'legend_background_color_active':
            self.legend_background_color_active = value
        elif variable == 'legend_position':
            self.legend_position = value
        elif variable == 'legend_ncol':
            self.legend_ncol = value
        elif variable == 'legend_label_weight':
            self.legend_label_weight = value
        elif variable == 'legend_label_style':
            self.legend_label_style = value
        elif variable == 'legend_label_size':
            self.legend_label_size = value
        elif variable == 'legend_label_color':
            self.legend_label_color = value
        elif variable == 'legend_title_weight':
            self.legend_title_weight = value
        elif variable == 'legend_title_style':
            self.legend_title_style = value
        elif variable == 'legend_title_size':
            self.legend_title_size = value
        elif variable == 'legend_title_color':
            self.legend_title_color = value

    def write(self,legendobj):
        lines = ""
        if self.legend_display                 != default_values['Legend']['legend_display']    :
            lines += '''    %s.setValue('legend_display',%s)\n'''%(legendobj,bool(self.legend_display))
        if self.legend_title                   != default_values['Legend']['legend_title']    :
            lines += '''    %s.setValue('legend_title','%s')\n'''%(legendobj,self.legend_title)
        if float(self.legend_border_width)            != default_values['Legend']['legend_border_width']    :
            lines += '''    %s.setValue('legend_border_width',%s)\n'''%(legendobj,self.legend_border_width)
        if self.legend_border_color            != default_values['Legend']['legend_border_color']    :
            lines += '''    %s.setValue('legend_border_color','%s')\n'''%(legendobj,self.legend_border_color)
        if self.legend_background_color        != default_values['Legend']['legend_background_color']    :
            lines += '''    %s.setValue('legend_background_color','%s')\n'''%(legendobj,self.legend_background_color)
        if self.legend_background_color_active != default_values['Legend']['legend_background_color_active']    :
            lines += '''    %s.setValue('legend_background_color_active',%s)\n'''%(legendobj,bool(self.legend_background_color_active))
        if self.legend_position                != default_values['Legend']['legend_position']    :
            lines += '''    %s.setValue('legend_position','%s')\n'''%(legendobj,self.legend_position)
        if int(self.legend_ncol)                    != default_values['Legend']['legend_ncol']    :
            lines += '''    %s.setValue('legend_ncol',%s)\n'''%(legendobj,self.legend_ncol)
        if self.legend_label_weight            != default_values['Legend']['legend_label_weight']    :
            lines += '''    %s.setValue('legend_label_weight','%s')\n'''%(legendobj,self.legend_label_weight)
        if self.legend_label_style             != default_values['Legend']['legend_label_style']    :
            lines += '''    %s.setValue('legend_label_style','%s')\n'''%(legendobj,self.legend_label_style)
        if float(self.legend_label_size)              != default_values['Legend']['legend_label_size']    :
            lines += '''    %s.setValue('legend_label_size',%s)\n'''%(legendobj,self.legend_label_size)
        if self.legend_label_color             != default_values['Legend']['legend_label_color']    :
            lines += '''    %s.setValue('legend_label_color','%s')\n'''%(legendobj,self.legend_label_color)
        if self.legend_title_weight            != default_values['Legend']['legend_title_weight']    :
            lines += '''    %s.setValue('legend_title_weight','%s')\n'''%(legendobj,self.legend_title_weight)
        if self.legend_title_style             != default_values['Legend']['legend_title_style']    :
            lines += '''    %s.setValue('legend_title_style','%s')\n'''%(legendobj,self.legend_title_style)
        if float(self.legend_title_size)              != default_values['Legend']['legend_title_size']    :
            lines += '''    %s.setValue('legend_title_size',%s)\n'''%(legendobj,self.legend_title_size)
        if self.legend_title_color             != default_values['Legend']['legend_title_color']    :
            lines += '''    %s.setValue('legend_title_color','%s')\n'''%(legendobj,self.legend_title_color)


        return lines
# ==============================================================================
# ==============================================================================
class AxisGrid:
    def __init__(self,display,grid_color,grid_style,grid_width,grid_tick_number,grid_tick_size):
        self.display    = display
        self.grid_color = grid_color
        self.grid_style = grid_style
        self.grid_width = grid_width
        self.grid_tick_number = grid_tick_number
        self.grid_tick_size = grid_tick_size
    def setValue(self,variable,value):
        if variable == 'display':
            self.display = value
        elif variable == 'grid_color':
            self.grid_color = value
        elif variable == 'grid_style':
            self.grid_style = value
        elif variable == 'grid_width':
            self.grid_width = value
        elif variable == 'grid_tick_number':
            self.grid_tick_number = value
        elif variable == 'grid_tick_size':
            self.grid_tick_size = value

# ==============================================================================
# ==============================================================================
class LevelGrid():
    def __init__(self,x_display,x_grid_color,x_grid_style,x_grid_width,x_grid_tick_number,x_grid_tick_size,
                      y_display,y_grid_color,y_grid_style,y_grid_width,y_grid_tick_number,y_grid_tick_size):
        self.x = AxisGrid(x_display,x_grid_color,x_grid_style,x_grid_width,x_grid_tick_number,x_grid_tick_size)
        self.y = AxisGrid(y_display,y_grid_color,y_grid_style,y_grid_width,y_grid_tick_number,y_grid_tick_size)
    def setValue(self,direction,variable,value):
        if direction == 'x':
            self.x.setValue(variable,value)
        else:
            self.y.setValue(variable,value)
# ==============================================================================
# ==============================================================================
class Grid:
    """
    Grid contains the main grid and the second grid. They can be configured directly by accessing to Grid or to the proper GridLevel object.
    In case of a multiple axis usage, then multiple Grid objects can be attached to a given plot.
    """
    def __init__(self,
                    Mx_display = default_values['Grid']['Mx_display'],
                    Mx_grid_color = default_values['Grid']['Mx_grid_color'],
                    Mx_grid_style = default_values['Grid']['Mx_grid_style'],
                    Mx_grid_width = default_values['Grid']['Mx_grid_width'],
                    Mx_grid_tick_number = default_values['Grid']['Mx_grid_tick_number'],
                    Mx_grid_tick_size = default_values['Grid']['Mx_grid_tick_size'],
                    My_display = default_values['Grid']['My_display'],
                    My_grid_color = default_values['Grid']['My_grid_color'],
                    My_grid_style = default_values['Grid']['My_grid_style'],
                    My_grid_width = default_values['Grid']['My_grid_width'],
                    My_grid_tick_number = default_values['Grid']['My_grid_tick_number'],
                    My_grid_tick_size = default_values['Grid']['My_grid_tick_size'],
                    mx_display = default_values['Grid']['mx_display'],
                    mx_grid_color = default_values['Grid']['mx_grid_color'],
                    mx_grid_style = default_values['Grid']['mx_grid_style'],
                    mx_grid_width = default_values['Grid']['mx_grid_width'],
                    mx_grid_tick_number = default_values['Grid']['mx_grid_tick_number'],
                    mx_grid_tick_size = default_values['Grid']['mx_grid_tick_size'],
                    my_display = default_values['Grid']['my_display'],
                    my_grid_color = default_values['Grid']['my_grid_color'],
                    my_grid_style = default_values['Grid']['my_grid_style'],
                    my_grid_width = default_values['Grid']['my_grid_width'],
                    my_grid_tick_number = default_values['Grid']['my_grid_tick_number'],
                    my_grid_tick_size = default_values['Grid']['my_grid_tick_size']):

        self.major = LevelGrid(Mx_display,Mx_grid_color,Mx_grid_style,Mx_grid_width,Mx_grid_tick_number,Mx_grid_tick_size,
                               My_display,My_grid_color,My_grid_style,My_grid_width,My_grid_tick_number,My_grid_tick_size)
        self.minor = LevelGrid(mx_display,mx_grid_color,mx_grid_style,mx_grid_width,mx_grid_tick_number,mx_grid_tick_size,
                               my_display,my_grid_color,my_grid_style,my_grid_width,my_grid_tick_number,my_grid_tick_size)
    def setValue(self,level,direction,variable,value):
        if level == 'major':
            self.major.setValue(direction,variable,value)
        else:
            self.minor.setValue(direction,variable,value)
    def write(self,gridobj):
        lines = ""
        if self.major.x.display != default_values['Grid']['Mx_display']    :
            lines +='''    %s.setValue('major','x','display',%s)\n'''%(gridobj,bool(self.major.x.display))
        if self.major.x.grid_color != default_values['Grid']['Mx_grid_color']    :
            lines +='''    %s.setValue('major','x','grid_color','%s')\n'''%(gridobj,self.major.x.grid_color)
        if self.major.x.grid_style != default_values['Grid']['Mx_grid_style']    :
            lines +='''    %s.setValue('major','x','grid_style','%s')\n'''%(gridobj,self.major.x.grid_style)
        if float(self.major.x.grid_width) != default_values['Grid']['Mx_grid_width']    :
            lines +='''    %s.setValue('major','x','grid_width',%s)\n'''%(gridobj,self.major.x.grid_width)
        if int(self.major.x.grid_tick_number) != default_values['Grid']['Mx_grid_tick_number']    :
            lines +='''    %s.setValue('major','x','grid_tick_number',%s)\n'''%(gridobj,self.major.x.grid_tick_number)
        if int(self.major.x.grid_tick_size) != default_values['Grid']['Mx_grid_tick_size']    :
            lines +='''    %s.setValue('major','x','grid_tick_size',%s)\n'''%(gridobj,self.major.x.grid_tick_size)
        if self.major.y.display != default_values['Grid']['My_display']    :
            lines +='''    %s.setValue('major','y','display',%s)\n'''%(gridobj,bool(self.major.y.display))
        if self.major.y.grid_color != default_values['Grid']['My_grid_color']    :
            lines +='''    %s.setValue('major','y','grid_color','%s')\n'''%(gridobj,self.major.y.grid_color)
        if self.major.y.grid_style != default_values['Grid']['My_grid_style']    :
            lines +='''    %s.setValue('major','y','grid_style','%s')\n'''%(gridobj,self.major.y.grid_style)
        if float(self.major.y.grid_width) != default_values['Grid']['My_grid_width']    :
            lines +='''    %s.setValue('major','y','grid_width',%s)\n'''%(gridobj,self.major.y.grid_width)
        if int(self.major.y.grid_tick_number) != default_values['Grid']['My_grid_tick_number']    :
            lines +='''    %s.setValue('major','y','grid_tick_number',%s)\n'''%(gridobj,self.major.y.grid_tick_number)
        if int(self.major.y.grid_tick_size) != default_values['Grid']['My_grid_tick_size']    :
            lines +='''    %s.setValue('major','y','grid_tick_size',%s)\n'''%(gridobj,self.major.y.grid_tick_size)
        if self.minor.x.display != default_values['Grid']['mx_display']    :
            lines +='''    %s.setValue('minor','x','display',%s)\n'''%(gridobj,self.minor.x.display)
        if self.minor.x.grid_color != default_values['Grid']['mx_grid_color']    :
            lines +='''    %s.setValue('minor','x','grid_color','%s')\n'''%(gridobj,self.minor.x.grid_color)
        if self.minor.x.grid_style != default_values['Grid']['mx_grid_style']    :
            lines +='''    %s.setValue('minor','x','grid_style','%s')\n'''%(gridobj,self.minor.x.grid_style)
        if float(self.minor.x.grid_width) != default_values['Grid']['mx_grid_width']    :
            lines +='''    %s.setValue('minor','x','grid_width',%s)\n'''%(gridobj,self.minor.x.grid_width)
        if int(self.minor.x.grid_tick_number) != default_values['Grid']['mx_grid_tick_number']    :
            lines +='''    %s.setValue('minor','x','grid_tick_number',%s)\n'''%(gridobj,self.minor.x.grid_tick_number)
        if int(self.minor.x.grid_tick_size) != default_values['Grid']['mx_grid_tick_size']    :
            lines +='''    %s.setValue('minor','x','grid_tick_size',%s)\n'''%(gridobj,self.minor.x.grid_tick_size)
        if self.minor.y.display != default_values['Grid']['my_display']    :
            lines +='''    %s.setValue('minor','y','display',%s)\n'''%(gridobj,self.minor.y.display)
        if self.minor.y.grid_color != default_values['Grid']['my_grid_color']    :
            lines +='''    %s.setValue('minor','y','grid_color','%s')\n'''%(gridobj,self.minor.y.grid_color)
        if self.minor.y.grid_style != default_values['Grid']['my_grid_style']    :
            lines +='''    %s.setValue('minor','y','grid_style','%s')\n'''%(gridobj,self.minor.y.grid_style)
        if float(self.minor.y.grid_width) != default_values['Grid']['my_grid_width']:
            lines +='''    %s.setValue('minor','y','grid_width',%s)\n'''%(gridobj,self.minor.y.grid_width)
        if int(self.minor.y.grid_tick_number) != default_values['Grid']['my_grid_tick_number']:
            lines +='''    %s.setValue('minor','y','grid_tick_number',%s)\n'''%(gridobj,self.minor.y.grid_tick_number)
        if int(self.minor.y.grid_tick_size) != default_values['Grid']['my_grid_tick_size']:
            lines +='''    %s.setValue('minor','y','grid_tick_size',%s)\n'''%(gridobj,self.minor.y.grid_tick_size)
        return lines

# ==============================================================================
# ==============================================================================
class Curve:
    """
    Curve class describes all the settings concerning a given curve itself.
    """
    def __init__(self,
                zone                          ,
                varx                          = default_values['Curve']['varx'],
                vary                          = default_values['Curve']['vary'],
                line_color                    = default_values['Curve']['line_color'],
                line_style                    = default_values['Curve']['line_style'],
                line_width                    = default_values['Curve']['line_width'],
                marker_style                  = default_values['Curve']['marker_style'],
                marker_size                   = default_values['Curve']['marker_size'],
                marker_edge_width             = default_values['Curve']['marker_edge_width'],
                marker_face_color             = default_values['Curve']['marker_face_color'],
                marker_edge_color             = default_values['Curve']['marker_edge_color'],
                marker_sampling_start         = default_values['Curve']['marker_sampling_start'],
                marker_sampling_end           = default_values['Curve']['marker_sampling_end'],
                marker_sampling_step          = default_values['Curve']['marker_sampling_step'],
                legend_label                  = default_values['Curve']['legend_label'],
                legend_display                = default_values['Curve']['legend_display'],
                ind_axis                      = default_values['Curve']['ind_axis'],
                visible                       = default_values['Curve']['visible'],
                axis                          = None):

        self.zone                  = zone
        self.varx                  = varx
        self.vary                  = vary
        self.line_color            = line_color
        self.line_style            = line_style
        self.line_width            = float(line_width)
        self.marker_style          = marker_style
        self.marker_size           = float(marker_size)
        self.marker_edge_width     = float(marker_edge_width)
        self.marker_face_color     = marker_face_color
        self.marker_edge_color     = marker_edge_color
        self.marker_sampling_start = marker_sampling_start
        self.marker_sampling_end   = marker_sampling_end
        self.marker_sampling_step  = marker_sampling_step
        self.legend_label          = legend_label
        self.legend_display        = legend_display
        #
        if axis:
            axis = axis.getInd()
        else:
            axis = ind_axis
        self.axis                  = axis
        self.visible               = visible
    # ------------------------------------------------------------- correctColor
    def correctColor(self,ind):
        cm = plt.get_cmap(COLOR_MAP)
        color = cm(1.*(ind%NUM_COLORS)/NUM_COLORS)
        html_color = '#%02x%02x%02x' % (int(255*color[0]),int(255*color[1]),int(255*color[2]))
        if self.line_color is None:
            self.line_color = html_color
        if self.marker_face_color is None:
            self.marker_face_color = html_color
        if self.marker_edge_color is None :
            self.marker_edge_color = html_color

    # ------------------------------------------------------------------- setVal
    def setValue(self,variable,value):
        if   variable == 'zone':
            self.zone = value
        elif variable == 'varx':
            self.varx = value
        elif variable =='vary':
            self.vary = value
        elif variable =='line_color':
            self.line_color = value
        elif variable =='line_style':
            self.line_style = value
        elif variable =='line_width':
            self.line_width = value
        elif variable =='marker_style':
            self.marker_style = value
        elif variable =='marker_size':
            self.marker_size = value
        elif variable =='marker_edge_width':
            self.marker_edge_width = value
        elif variable =='marker_face_color':
            self.marker_face_color = value
        elif variable =='marker_edge_color':
            self.marker_edge_color = value
        elif variable =='marker_sampling_start':
            self.marker_sampling_start = value
        elif variable =='marker_sampling_end':
            self.marker_sampling_end = value
        elif variable =='marker_sampling_step':
            self.marker_sampling_step = value
        elif variable =='legend_label':
            self.legend_label = value
        elif variable =='legend_display':
            self.legend_display = value
        elif variable =='ind_axis':
            self.axis = value
        elif variable =='axis':
            ind = value.getInd()
            self.axis = ind
        elif variable =='visible':
            self.visible = value
    # -------------------------------------------------------------------- write
    def write(self,indgraph,iCurSubGraph,indsubgraph,indcurve):
        lines = '''    curve_%s = Curve(zone=%s,varx='%s',vary='%s' '''%(indcurve,self.zone,self.varx,self.vary)
        if self.line_color            != default_values['Curve']['line_color']    :
            lines += ''', line_color='%s' '''%(self.line_color)
        if self.line_style            != default_values['Curve']['line_style']    :
            lines += ''', line_style='%s' '''%(self.line_style)
        if float(self.line_width)     != default_values['Curve']['line_width']    :
            lines += ''', line_width=%s'''%(self.line_width)
        if self.marker_style          != default_values['Curve']['marker_style']    :
            lines += ''', marker_style='%s' '''%(self.marker_style)
        if float(self.marker_size)    != default_values['Curve']['marker_size']    :
            lines += ''', marker_size=%s'''%(self.marker_size)
        if float(self.marker_edge_width) != default_values['Curve']['marker_edge_width']    :
            lines += ''', marker_edge_width=%s'''%(self.marker_edge_width)
        if self.marker_face_color     != default_values['Curve']['marker_face_color']    :
            lines += ''', marker_face_color='%s' '''%(self.marker_face_color)
        if self.marker_edge_color     != default_values['Curve']['marker_edge_color']    :
            lines += ''',marker_edge_color='%s' '''%(self.marker_edge_color)
        if self.marker_sampling_start != default_values['Curve']['marker_sampling_start']    :
            lines += ''', marker_sampling_start=%s'''%(self.marker_sampling_start)
        if self.marker_sampling_end   != default_values['Curve']['marker_sampling_end']    :
            lines += ''', marker_sampling_end=%s'''%(self.marker_sampling_end)
        if self.marker_sampling_step  != default_values['Curve']['marker_sampling_step']    :
            lines += ''', marker_sampling_step=%s'''%(self.marker_sampling_step)
        if self.legend_label          != default_values['Curve']['legend_label']    :
            lines += ''', legend_label='%s' '''%(self.legend_label)
        if self.legend_display        != default_values['Curve']['legend_display']    :
            lines += ''', legend_display=%s'''%(self.legend_display)
        if self.axis                  != default_values['Curve']['ind_axis']    :
            lines += ''', axis=axis_%s_%s_%s'''%(indgraph,indsubgraph,self.axis)
        if self.visible               != default_values['Curve']['visible']    :
            lines += ''', visible=%s'''%(self.visible)
        #
        lines += ''')\n'''
        #
        lines += '''    graph_%s.addCurve('%s',curve_%s)'''%(indgraph,iCurSubGraph,indcurve)
        return lines


# ==============================================================================
# ==============================================================================
class SubGraph:

    def __init__(self,figure,nl,nc,ind,il,ic):

        self.figure                = figure
        self.nl                    = nl
        self.nc                    = nc
        self.ind                   = ind
        self.il                    = il
        self.ic                    = ic
        self.legend_property       = Legend()
        self.axis                  = [self.figure.add_subplot(self.nl,self.nc,self.ind+1)]
        self.axis_property         = [Axis(ind=0)]
        self.grid_property         = [Grid()]
        self.curves                = []
        self.name                  = '%s:%s'%(self.il+1, self.ic+1)
        #
        self.axis[0].type = ['main',None]
        #
        self.axis[0].name = self.name
        self.axis_property[-1].name = self.name
        self.grid_property[-1].name = self.name

    def addAxisTwinX(self,axis_to_twin):

        # Could use twinx() method ...
        self.axis.append(self.figure.add_axes(self.axis[axis_to_twin].get_position(), frameon=False,sharex=self.axis[axis_to_twin],label='%s'%len(self.axis)))
        self.axis[-1].type = ['twinx',axis_to_twin]
        self.axis[-1].name = self.name
        self.axis_property.append(Axis(ind=len(self.axis_property)-1))
        self.axis_property[-1].name = self.name
        self.grid_property.append(Grid())
        #
        # ### Link Axis property
        #
        self.axis_property[-1].x = self.axis_property[axis_to_twin].x
        #
        # ### Link Grid property
        #
        self.grid_property[-1].major.x = self.grid_property[axis_to_twin].major.x
        self.grid_property[-1].minor.x = self.grid_property[axis_to_twin].minor.x
        ###
#        self.axis[-1].set_frame_on(False)
#        self.axis[-1].patch.set_visible(True)
#        for sp in self.axis[-1].spines.values():
#            sp.set_visible(False)
#        self.axis[-1].patch.set_alpha(0.)



    def addAxisTwinY(self,axis_to_twin):
        # Could use twiny() method ...
        self.axis.append(self.figure.add_axes(self.axis[axis_to_twin].get_position(), frameon=False,sharey=self.axis[axis_to_twin],label='%s'%len(self.axis)))
        self.axis[-1].type = ['twiny',axis_to_twin]
        self.axis[-1].name = self.name
        self.axis_property.append(Axis(ind=len(self.axis_property)-1))
        self.axis_property[-1].name = self.name
        self.grid_property.append(Grid())
        #
        # ### Link Axis property
        #
        self.axis_property[-1].y = self.axis_property[axis_to_twin].y
        #
        # ### Link Grid property
        #
        self.grid_property[-1].major.y = self.grid_property[axis_to_twin].major.y
        self.grid_property[-1].minor.y = self.grid_property[axis_to_twin].minor.y
        ###
#        self.axis[-1].set_frame_on(False)
#        self.axis[-1].patch.set_visible(True)
#        for sp in self.axis[-1].spines.values():
#            sp.set_visible(False)




    def addAxis(self):
        self.axis.append(self.figure.add_axes(self.axis[0].get_position(), frameon=False,label='%s'%len(self.axis)))
        self.axis[-1].type = ['new',None]
        self.axis[-1].name = self.name
        self.axis_property.append(Axis(ind=len(self.axis_property)-1))
        self.axis_property[-1].name = self.name
        self.grid_property.append(Grid())

# ==============================================================================
# ==============================================================================
class SubPlotParams:
    """
    SubPlotParams is one way (TightLayout) to set margin, padding for plots positionning inside the Graph window.
    """
    def __init__(self,isActive=default_values['SubPlotParams']['isActive'],
                      left=default_values['SubPlotParams']['left'],
                      right=default_values['SubPlotParams']['right'],
                      top=default_values['SubPlotParams']['top'],
                      bottom=default_values['SubPlotParams']['bottom'],
                      wspace=default_values['SubPlotParams']['wspace'],
                      hspace=default_values['SubPlotParams']['hspace']):

        self.left        = left
        self.right       = right
        self.top         = top
        self.bottom      = bottom
        self.wspace      = wspace
        self.hspace      = hspace
        self.isActive    = isActive
    # -------------------------------------------------------------- checkParams
    def checkParams(self,params):
        isOk = True
        if params['left'] is not None and params['right'] is not None:
            left = params['left']
            right = params['right']
            if right<left:
                tkMessageBox.showwarning('Configure SubPlotParams failed','Left Value must be smaller than right value')
                isOk = False
        if params['bottom'] is not None and params['top'] is not None:
            bottom = params['bottom']
            top = params['top']
            if  top<bottom:
                tkMessageBox.showwarning('Configure SubPlotParams failed','Bottom value must be smaller than top value')
                isOk = False
        return isOk

    # ----------------------------------------------------------------- setValue
    def setValue(self,var,val):
        if var == 'left':
            self.left = val
        elif var == 'right':
            self.right = val
        elif var == 'top':
            self.top = val
        elif var == 'bottom':
            self.bottom = val
        elif var == 'hspace':
            self.hspace = val
        elif var == 'wspace':
            self.wspace = val
        elif var == 'isActive':
            self.isActive = val
    # -------------------------------------------------------------------- write
    def write(self,indgraph):
        line = ""
        if self.isActive:
            string = ''''isActive':True,'''
            if self.left != default_values['SubPlotParams']['left']:
                string+=''''left':%s,'''%(self.left)
            if self.right != default_values['SubPlotParams']['right']:
                string+=''''right':%s,'''%(self.right)
            if self.top != default_values['SubPlotParams']['top']:
                string+=''''top':%s,'''%(self.top)
            if self.bottom != default_values['SubPlotParams']['bottom']:
                string+=''''bottom':%s,'''%(self.bottom)
            if self.hspace != default_values['SubPlotParams']['hspace']:
                string+=''''hspace':%s,'''%(self.hspace)
            if self.wspace != default_values['SubPlotParams']['wspace']:
                string+=''''wspace':%s,'''%(self.wspace)


            line = '''    graph_%s.updateSubPlotParams({%s})\n'''%(indgraph,string[:-1])
        return line
# ==============================================================================
# ==============================================================================
class EventResize:

    def __init__(self,data):

        self.width        = data[0]
        self.height       = data[1]
# ==============================================================================
# ==============================================================================
class TightLayout:
    """
    TightLayout is one way (SubPlotParams) to set margin, padding for plots positionning inside the Graph window.
    """
    def __init__(self,isActive= not default_values['SubPlotParams']['isActive'],
                      pad=default_values['TightLayout']['pad'],
                      hpad=default_values['TightLayout']['hpad'],
                      wpad=default_values['TightLayout']['wpad']):

        self.pad        = pad
        self.hpad       = hpad
        self.wpad       = wpad
        self.isActive    = isActive
    # ----------------------------------------------------------------- setValue
    def setValue(self,var,val):
        if var == 'pad':
            self.pad = val
        elif var == 'hpad':
            self.hpad = val
        elif var == 'wpad':
            self.wpad = val
        elif var == 'isActive':
            self.isActive = val
    # -------------------------------------------------------------------- write
    def write(self,indgraph):
        line = ""
        if self.isActive:
            string = ''''isActive':True,'''
            if self.pad != default_values['TightLayout']['pad']:
                string+=''''pad':%s,'''%(self.pad)
            if self.hpad != default_values['TightLayout']['hpad']:
                string+=''''hpad':%s,'''%(self.hpad)
            if self.wpad != default_values['TightLayout']['wpad']:
                string+=''''wpad':%s,'''%(self.wpad)
            line = '''    graph_%s.updateTightLayout({%s})\n'''%(indgraph,string[:-1])
        return line
# ==============================================================================
# ==============================================================================
# --------------------------------------------------


class Movie(object):
    """
    Class Movie can be used to generate a movie in case of a dynamic plot (Co-processing for instance)
    """

    def __init__(self, fig, filename, fps=10):
        self.fig      = fig
        self.filename = filename
        self.fps      = fps

    def activate(self):
        width, height = self.fig.canvas.get_width_height()
        print("width = ",width,", height = ",height)

############################################################################################################
#BERTRAND :
#        cmd = ('ffmpeg',
#               '-y', '-r', '%d' % self.fps,                   # overwrite, 1fps
#               '-s', '%dx%d' % (width, height),               # size of image string
#               '-pix_fmt', 'argb',                            # format
#               '-f', 'rawvideo',  '-i', '-',                  # tell ffmpeg to expect raw video from the pipe
#               '-vcodec', 'mpeg4', self.filename)             # output encoding
#        self.p = subprocess.Popen(cmd, stdin=subprocess.PIPE)
############################################################################################################
#        cmd = '''ffmpeg -y -r %s -s %dx%d -pix_fmt argb -f rawvideo -r 20 -i - -c:v mpeg1video -vb 50000K -r 30 %s'''%(self.fps,width,height,self.filename)
        cmd = '''ffmpeg -y -r %s -s %dx%d -pix_fmt argb -f rawvideo -i - -c:v mpeg1video -vb 100000K -r 30 %s'''%(self.fps,width,height,self.filename)
        self.p = subprocess.Popen(shlex.split(cmd), stdin=subprocess.PIPE)


    def write(self):
        # write to pipe
        self.p.stdin.write(self.fig.canvas.tostring_argb())
#        self.fig.savefig(self.p.stdin, format='png')

#    def exit(self, type, value, traceback):
    def exit(self):
        print("Finalize(exit) Movie = ",self.filename)
        self.p.communicate()
# ==============================================================================
# ==============================================================================
# --------------------------------------------------
#class Movie(object):

#    def __init__(self, fig, filename, fps=2):
#        self.fig      = fig
#        self.filename = filename
#        self.fps      = fps

#    def activate(self):
#        import matplotlib.animation as animation
#        FFMpegWriter = animation.writers['ffmpeg']
#        self.moviewritter = FFMpegWriter(fps=self.fps)
#        self.moviewritter.setup(self.fig,self.filename,500)

#    def write(self):
#        # write to pipe
#        self.moviewritter.grab_frame()

##    def exit(self, type, value, traceback):
#    def exit(self):
#        print("Finalize(exit) Movie = ",self.filename)
#        self.moviewritter.finish()

# ==============================================================================
# ==============================================================================
# --------------------------------------------------
class CustomToolbar(NavigationToolbar2TkAgg):
    def __init__(self,canvas,parent,graph):
        self.toolitems = (
        ('Home', 'Reset original view', 'initial', 'home'),
        ('Back', 'Back to  previous view', 'previous', 'back'),
        ('Forward', 'Forward to next view', 'next', 'forward'),
        (None, None, None, None),
        ('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan_tkPlotXY'),
        ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom_tkPlotXY'),
        (None, None, None, None),
        ('Subplots', 'Configure subplots', 'subplots', 'configure_subplots_tkPlotXY'),
        ('Save', 'Save the figure', 'filesave', 'save_figure'),
      )
        self.graph = graph
        #
        self.button_dict = {}
        NavigationToolbar2TkAgg.__init__(self,canvas,parent)

    def _Button(self, text, file, command, extension='.ppm'):
        fileimage=file

        im = TK.PhotoImage(data=IMAGE_DICT[fileimage])
        b = TK.Button(
            master=self, text=text, padx=2, pady=2, image=im, command=command)
        b._ntimage = im
        b._image = fileimage
        self.button_dict[fileimage]=b
        b.pack(side=TK.LEFT)
        return b

    def pan_tkPlotXY(self, *args):
        if self._active == 'ZOOM':
            im = TK.PhotoImage(data=IMAGE_DICT['zoom_to_rect'])
            self.button_dict['zoom_to_rect'].config(image=im)
            self.button_dict['zoom_to_rect']._ntimage = im
            im = TK.PhotoImage(data=IMAGE_DICT['move_active'])
            self.button_dict['move'].config(image=im)
            self.button_dict['move']._ntimage = im
        elif self._active == 'PAN':
            im = TK.PhotoImage(data=IMAGE_DICT['move'])
            self.button_dict['move'].config(image=im)
            self.button_dict['move']._ntimage = im
        else:
            im = TK.PhotoImage(data=IMAGE_DICT['move_active'])
            self.button_dict['move'].config(image=im)
            self.button_dict['move']._ntimage = im
        self.pan(self, *args)

    def zoom_tkPlotXY(self, *args):
        if self._active == 'PAN':
            im = TK.PhotoImage(data=IMAGE_DICT['move'])
            self.button_dict['move'].config(image=im)
            self.button_dict['move']._ntimage = im
            im = TK.PhotoImage(data=IMAGE_DICT['zoom_to_rect_active'])
            self.button_dict['zoom_to_rect'].config(image=im)
            self.button_dict['zoom_to_rect']._ntimage = im
        elif self._active == 'ZOOM':
            im = TK.PhotoImage(data=IMAGE_DICT['zoom_to_rect'])
            self.button_dict['zoom_to_rect'].config(image=im)
            self.button_dict['zoom_to_rect']._ntimage = im
        else:
            im = TK.PhotoImage(data=IMAGE_DICT['zoom_to_rect_active'])
            self.button_dict['zoom_to_rect'].config(image=im)
            self.button_dict['zoom_to_rect']._ntimage = im
        self.zoom(self,*args)

    def configure_subplots_tkPlotXY(self):
        toolfig = matplotlib.figure.Figure(figsize=(6,3))
        window = TK.Tk()
        canvas = FigureCanvasTkAgg(toolfig, master=window)
        toolfig.subplots_adjust(top=0.9)
        tool =  CustomSubplotTool(self.canvas.figure, toolfig)
        tool.graph = self.graph
        canvas.show()
        canvas.get_tk_widget().pack(side=TK.TOP, fill=TK.BOTH, expand=1)


# ==============================================================================
# ==============================================================================
# --------------------------------------------------
class CustomSubplotTool(SubplotTool):
    def funcleft(self, val):
        self.targetfig.subplots_adjust(left=val)
        if self.drawon:
            self.targetfig.canvas.draw()
        #
        self.graph.subPlotParams.isActive=True
        self.graph.subPlotParams.left=val

    def funcright(self, val):
        self.targetfig.subplots_adjust(right=val)
        if self.drawon:
            self.targetfig.canvas.draw()
        #
        self.graph.subPlotParams.isActive=True
        self.graph.subPlotParams.right=val

    def funcbottom(self, val):
        self.targetfig.subplots_adjust(bottom=val)
        if self.drawon:
            self.targetfig.canvas.draw()
        #
        self.graph.subPlotParams.isActive=True
        self.graph.subPlotParams.bottom=val

    def functop(self, val):
        self.targetfig.subplots_adjust(top=val)
        if self.drawon:
            self.targetfig.canvas.draw()
        #
        self.graph.subPlotParams.isActive=True
        self.graph.subPlotParams.top=val

    def funcwspace(self, val):
        self.targetfig.subplots_adjust(wspace=val)
        if self.drawon:
            self.targetfig.canvas.draw()
        #
        self.graph.subPlotParams.isActive=True
        self.graph.subPlotParams.wspace=val

    def funchspace(self, val):
        self.targetfig.subplots_adjust(hspace=val)
        if self.drawon:
            self.targetfig.canvas.draw()
        #
        self.graph.subPlotParams.isActive=True
        self.graph.subPlotParams.hspace=val

# ==============================================================================
# ==============================================================================

# Using GUI

def main(data):
    # Ouverture de la fenetre principale
    desktop = DesktopTK(None)
    desktop.setData(data)
    desktop.mainloop()
    desktop.quit()
    print('end of programm')

# ==============================================================================
# ==============================================================================

# Without GUI

class GraphEditor():
    """
    The class GraphEditor is an encapsulation of class Desktop.
    """
    def __init__(self,display):
        self.desktop = Desktop()
        self.desktop.display = display


    def __enter__(self,display):
        print("Enter")
        return self.desktop


    def __exit__(self, type, value, traceback):
        print("Exit")
        for graph in self.desktop.graphWdwL:
            graph.quit()


    def close(self):
        print("Close")
        for graph in self.desktop.graphWdwL:
            graph.quit()


def openGraphEditor(display):
    """
    Create an object of class GraphEditor and returns its Desktop.
    """
    editor = GraphEditor(display)
    return editor.desktop
# ==============================================================================
# ==============================================================================


def filterInteger(string):
    res = ''
    for i in string:
        if i in '0123456789':
            res += i
    return res


#===============================================================================
# Fonction permettant de pointer vers les zones 1D de l'arbre
#===============================================================================
def updateFromTree(event=None):
    if CTK.__MAINTREE__ == 1: tp = CTK.t
    else: tp = CTK.dt

    to = C.newPyTree()
    bases = Internal.getBases(tp)
    for b in bases: C._addBase2PyTree(to, b[0], 1)

    # Get only 1D zones of tp
    for b in bases:
        bn = Internal.getNodeFromName1(to, b[0]) # same base in to
        for z in Internal.getZones(b):
            dim = Internal.getZoneDim(z)
            zname = z[0]
            if dim[0] == 'Structured' and dim[2] == 1 and dim[3] == 1:
                z = C.center2Node(z, Internal.__FlowSolutionCenters__)
                bn[2].append(z)
            elif dim[0] == 'Unstructured' and dim[3] == 'BAR':
                zps = T.splitConnexity(z)
                zps = C.convertBAR2Struct(zps)
                c = 0
                for i in zps:
                    i[0] = zname+str(c); c += 1
                    i = C.center2Node(i, Internal.__FlowSolutionCenters__)
                bn[2] += zps

            elif dim[0] == 'Unstructured' and dim[3] == 'NGON' and dim[4] == 1:
                zp = C.convertArray2Hexa(z)
                zps = T.splitConnexity(zp)
                zps = C.convertBAR2Struct(zps)
                c = 0
                for i in zps:
                    i[0] = zname+str(c); c += 1
                    i = C.center2Node(i, Internal.__FlowSolutionCenters__)
                bn[2] += zps

    if CTK.__MAINTREE__ == 1:
        DESKTOP.setData(to)
    else:
        global PREVTPZONES
        if PREVTPZONES == []: DESKTOP.setData(to)
        else: DESKTOP.replaceGroupZones(to, PREVTPZONES)
        PREVTPZONES = Internal.getZonePaths(to, 2)

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    if not IMPORTOK: return

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkPlotXY', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='My personal applet.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<Control-u>', updateFromTree)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    FrameMenu.add_command(label='Update', accelerator='Ctrl+u', command=updateFromTree)
    CTK.addPinMenu(FrameMenu, 'tkPlotXY')
    WIDGETS['frameMenu'] = FrameMenu

    desktopFrameTK = DesktopFrameTK(Frame)
    desktopFrameTK.grid(row=0,column=0,sticky='NSEW')
#   NEW CHRISTOPHE
    global DESKTOP
    DESKTOP = desktopFrameTK
#   END NEW CHRISTOPHE

#    desktopFrameTK.setData(CTK.t) # commente par CB pour l'instant

##==============================================================================
## Create app widgets
##==============================================================================
#def createApp(win):
#
#    # - Frame -
#    Frame = TK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
#                          text='tkPlotXY', font=CTK.FRAMEFONT, takefocus=1)
#    #BB = CTK.infoBulle(parent=Frame, text='My personal applet.\nCtrl+c to close applet.', temps=0, btype=1)
#    Frame.bind('<Control-c>', hideApp)
#    Frame.bind('<Button-3>', displayFrameMenu)
#    Frame.columnconfigure(0, weight=1)
#    Frame.columnconfigure(1, weight=0)
#    Frame.rowconfigure(0, weight=1)
#    WIDGETS['frame'] = Frame

#    # - Frame menu -
#    FrameMenu = TK.Menu(Frame, tearoff=0)
#    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
#    FrameMenu.add_command(label='Save', command=saveApp)
#    FrameMenu.add_command(label='Reset', command=resetApp)
#    FrameMenu.add_command(label='Update', accelerator='Ctrl+u', command=updateFromTree)
#    CTK.addPinMenu(FrameMenu, 'tkPlotXY')
#    WIDGETS['frameMenu'] = FrameMenu

#    desktopFrameTK = DesktopFrameTK(Frame)
#    desktopFrameTK.grid(row=0,column=0,sticky='NSEW')
#    #
#    frameL = TK.Frame(Frame)
#    frameL.columnconfigure(0,weight=0)
#    frameL.rowconfigure(0,weight=0)
#    frameL.rowconfigure(1,weight=1)
#    frameL.grid(row=0,column=1,sticky='NSEW')
#    closeCross = TK.Button(frameL,text='x',command=cmd_hideApp,height=1, width=1)
#    closeCross.grid(row=0,column=1,sticky='N')
#    blankFrame = TK.Frame(frameL)
#    blankFrame.grid(row=1,column=0,sticky='NSEW')
#
##   NEW CHRISTOPHE
#    global DESKTOP
#    DESKTOP = desktopFrameTK
##   END NEW CHRISTOPHE

##    desktopFrameTK.setData(CTK.t) # commente par CB pour l'instant

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    WIDGETS['frame'].grid(sticky=TK.EW)

#==============================================================================
# Called to hide widgets
#==============================================================================
def cmd_hideApp(event=None):
    hideApp(event)

def hideApp(event=None):
    WIDGETS['frame'].grid_forget()

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp():
    # Update App seulement si une graph window est ouverte
    if IMPORTOK and CTK.t != [] and len(DESKTOP.graphWdwL) > 0:
        updateFromTree()
    return

#==============================================================================
def saveApp():
# A completer ...
    CTK.savePrefFile()

#==============================================================================
def resetApp():
# A completer ...
    CTK.savePrefFile()

#==============================================================================
def displayFrameMenu(event=None):
    WIDGETS['frameMenu'].tk_popup(event.x_root+50, event.y_root, 0)

#==============================================================================
# IMAGE RESSOURCES
IMAGE_DICT = {
'initial':'''R0lGODlhGAAYAOekAAAAAAsNDQ8RERIUFRQXGCgtLjc9P0lSVNd/MrePX9aJTNmLT8OXVNmRVNyRVN6RUtmSW5GiptyWXt2YVtGacOCYWsicfN2gcLyojdqjbOWhZdumfuipSaO2uuStad6uh5G/ya+5uOuva+uweeizZey1VpHFz961l9+3l5PI0u+7WpfJ3qHI0fS8VLDFyp3K0fO9VZ7K0u+6g6DK05/L06PK0qjJ0KTK0qHL1KTL0qTL06nK0KHM1aXL1KrK0KLM1KfL0/TAXqTM1KXM0qTM1fTAYaLN1anL1fXBYOLApajN1LTM1qTR16fQ2a3Q1PbFbarR2arR2vbGb7XP1qvS26zT2LnQ2KzU27fR27LT163U3bnR2a3V2bLU1/fJea7V27TU1/TIl/fKfLDW2/fLfrHX27HX3LXW3bfW2bfW2rHY2+LOvrrX27rX3LbY3/jOiMDX3L/X4PjPirrZ3uLQw7fa373Z3cHY4bjb4L7a3sLa37vc4fbRocDb38nY47/c4MnZ5cnZ5vnTl8Lc4crZ6fjToPnUmc3b6cbe4s7b6cbe4/LVwc/c6dDc6fnXouPZ0dDd6dHd6frYo9Ld6NTe587g8NXf6eXf29vh5+Xi3+Xj4uLk5+Xl5uTm5+bm5tro+Nvo99zp993p+N7q+P///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////yH+EUNyZWF0ZWQgd2l0aCBHSU1QACwAAAAAGAAYAAAI/gA9CRxIsKDBgwgTKlwocI2DBZkYGsx0QUMhESgkEkzyIIwhMkUQFNzUiJHJSJ0GPpIwwpEcKS1IKBjIKZGfHTlu9PARyJKnDxX4CBKDRMWEDZcEbiJ0RMufPnbYtKFixYIMSW+ewPAAgQ7NQyzcIFJENk+aLlFCZPASpESDE5oIQgJyRhEcAwICkC1TJgUGDgwSGNzk58ogRXoidNBLdg+TKRQWHWxkYw5Zsi4YK9qjhkgohIx0sLmsKDPpMjREga6BhrRpsnXKxBgF+gYY15oVlcFB++AkH2myXH69V0dvg5wAQXFCtsAAAAQO7L2jkNKWJmXwkN4LAtTCSUtMIBj5YmYMlyo/lFSSiGkFFiE8XswYEueTxvv48+vfjzAgADs=''',
'move':'''R0lGODlhGAAYAOe1AA5efg9efw5ffg9ffxBgfxFgfxFggBNhgBRhgRRigRRighVigRdkghhlgxplhBxlhBpmhBtmhBxnhR1ohR5ohiJriCNriCRsiSVsiSRtiShuiylviylwjCtwjCxxjS1xjS1yjS9zjjV3kTl5kzp5kz17lUSAmEiDmkuEnE2GnVaMoVeMoViNolmNolqOo1uPpF+RpWOUp2SUqNd/MmaVqGiXqm2arG6brbePX9aJTHSer9mLT3WfsHmisnqissOXVNmRVNyRVN6RUtmSW9yWXt2YVuCYWsicfN2gcLyojZGwvdqjbOWhZZSyv9umfuipSZ24w+Stad6uh5G/ya+5uOuva+uweeizZaa9x+y1VpHFz6m/yd61l9+3l5PI0u+7WpfJ3rDEzKHI0fS8VPO9Ve+6g6TK0qHL1KTL0qTL06nK0KHM1aXL1KrK0KfL0/TAXqTM1PTAYaLN1anL1fXBYOLApbTM1qfQ2a3Q1PbFbarR2fbGb7fR27LT163U3bnR2ffJea7V27TU18TQ1fTIl/fKfPfLfrHX27HX3LfW2sjT2OLOvrrX27bY3/jOiMDX3PjPiuLQw73Z3cHY4bjb4L7a3s7X2/bRocDb38nY47/c4MnZ5fnTl8rZ6fjToPnUmc3b6cbe4s7b6cbe48/c6dDc6fnXouPZ0dDd6frYo9Ld6NTe5+Xf29vh5+Xi3+Xj4uLk5+Xl5ubm5tvo997q+P///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////yH+EUNyZWF0ZWQgd2l0aCBHSU1QACwAAAAAGAAYAAAI/gBlCRxIsKDBgwdfIFxYsAcHGQcXBdnhaiGWFQA+QCHoCgkTT1W6ILRUYwEAACUUCawjhNAnQ3FmFIRVihSICQUAGHBw4RQRK6Yg7RlzJcfAWKIyqUFjho0IABtkSTFyiVMhOl+KOGElEFanOX40YZLEiATUI2VSOcpDJsqQSEdBiWkUapTdSiagUlkC6E0WIFxeEUTlxlKFR3bt5t3gJcmTHzgyqBwIK5MSDAdsaNZsIQAFGzp8iEYQosnAUikYnFy9OgFr1g9YCIyhosHr1rdPFkABkVSYGyMIwBg+HIKACDZoED+gAUYYgaTMCArTIbHdpxsOnaEl0MMWgqraWCTqY/061ENpuB+MtUkPnvKjsB+axHDVnzuHKFnHPmUWQ1mq2KGFHIEg4oEEBgAwgAIX/CdLK2DwAccalriwmgmDOGgQFidktJGGBvGgAUQgHtRCiSgWFBAAOw==''',
'move_active':'''R0lGODlhGAAYAOevAAtLZQxLZgtMZQxMZg1NZg5NZg9OZhBOZxBOaBFOZxJQaBNRaRVRahZRahVSahZSahdTahhTaxtWbRxWbR1Wbh5Wbh1XbiBYbyFZbyFacCJacCNacSRacSRbcSZccipfdC5hdjFidzZmejppezxqfT5rfkVwgUZwgUZxgkdxgkhygklyg0x0hE92hlB2hlJ3hqxmKFN5iFd7ilh8ipJyTKtuPV1+jK5vP15/jWGCjmKCjpx5Q650Q7B0Q7J0Qq51SbB4S7F6RbN6SKB9Y7GAWpaGca6CVnSNl7eBUXaOma+FZbqHOn6TnLaKVLKLbHSZoYyUk7yMVrqPUYWXn7yNYb2RRXSepoeZobKRebKSeXagqL+WSMOWQ3mhsoGgp42do8KXRL+VaYGiqoOiqIOiqYSiqoGjqoeipoaiqYOjqoiipoeiqsOaS4KkqsOaTsSaTbWahJCjq4amroqmqoinrsWeV8WeWZKnr5Snro6prIqqsYuqr8ahYZ2mqpCqrMOgecaiY8aiZY6sr5Krro6ssJWsr6CprbWlmJKtssalbZqssLWmnMambpeusZOvs5qttJiussWngaWsr5qvsqGttpmws6Gut8epecapgKKuuseqeqSvuqWvup6ytZ6ytqawusesgraup6axusitgqixuqqyubeyr6+0ube1srW2ube2tbe3uLi4uK+6xrK7xv///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////yH+EUNyZWF0ZWQgd2l0aCBHSU1QACwAAAAAGAAYAAAI/gBZCRxIsKDBgwdXIFxYMEcGFwcP9biBauGUEwA4MCGIiggSTFGyIJQUIwEAACEMCYTj44+mQG5gFEz16VMHCAUAFGBAIRQQKqAY2eEipcbAVZwonRkzpswHABhYORES6RKgN1uCKDElMFWmNXoqTWpUCATUIWFGJaoDpsmPRUc3eUHUyZNdSCKgQjHCh00VHlhUERSFRpIERXbt5sWgpciSHTQsqByYitKRCgZkaNY8IUAEGTZ0iD7gIcnATyUUnFy9+gBr1g1QCGxhYsHr1rdPFiAB8dOXGSAIsBg+3IGABzJeEDdwgcUXgZ/G+PmiIbHdpxgEiXElcMMVgqTUVwzKY/06VEFkuB9cZYnOnPKesAt6xLAUHjmCHFnH/qQVQ1akxGFFG3sQssEDOQ2AAAX/sXJKF3ekYYYkKqwmQh8NGjTFCBltlKFBOFwA0YcHpUDiiQUFBAA7''',
'zoom_to_rect':'''R0lGODlhGAAYAOfxADqInUKNoUSNo0aOpUeQo0yRqE6Up1OXqtd/MliarVicrFycsF6dsV2fr7ePX9aJTGKisdmLT2aktGSls2elt2mmtmintcOXVNmRVNyRVN6RUmypt9mSW26ruW+rudyWXt2YVtGacOCYWsicfHSvvHixvniyvHqzv92gcLyojdqjbOWhZYO5w9umfuipSXm+yIW8xIi8x4q9x37ByuStaY2/yt6uh4DDzJG/yYPDzYLEz6+5uOuva5LBzITFz+uweeizZey1VpHFz4rH0ZnExZbEzt61l5fF0IjJ1Y3I0d+3l5PI0ozK1IrL14/K1I3L1e+7WpfJ3qHI0fS8VJ3K0fO9VZ7K0u+6g5TN16DK05PO2J/L06PK0pTO2KjJ0KTK0qHL1KTL0qTL06nK0KHM1aXL1KrK0KLM1JfP2KfL0/TAXqTM1KXM0qTM1fTAYZbQ2qLN1anL1fXBYOLApZfR26jN1JnR25jS3LTM1qTR16fQ2a3Q1J3U3vbFbZ7U3qrR2arR2vbGb7XP1p/V4KvS26zT2LnQ2KzU27fR27LT163U3bnR2a3V2bLU1/fJea7V27TU1/TIl/fKfLDW2/fLfrHX27HX3LXW3bfW2bfW2rjW2rHY2+LOvrrX27rX3LbY3/jOiL/X4PjPirrZ3uLQw7fa373Z3cHY4bjb4L7a3rvc4fbRocDb38nY48Hb37/c4MnZ5cnZ5vnTl8Lc4crZ6fjToL7e4/nUmb7f47/f5M3b6cbe4s7b6cbe4/LVwcHg5c/c6dDc6cjf4/nXouPZ0dDd6dHd6frYo9Ld6MTi59Te58Xj587g8Mfj6NXf6cnl6uXf29vh58vm68zn687n7OXi3+Xj4uLk59Dp7tHp7uXl5tPq7+Tm5+bm5tro+Nvo99Tr8NXr8Nzp993p+N7q+Njt8tnu89vv9d7w9t/x9uL0+eP0+e3z+/D0+/L1+/n5/f78/v///////////////////////////////////////////////////////////yH+EUNyZWF0ZWQgd2l0aCBHSU1QACwAAAAAGAAYAAAI/gC7CRxIsKDBgwgTKlwokFOGCNUYGqyGYkUtHkokEpyjIdItSm4QFLwWDJhJY9wGEvvwY5ioQFOAPBiojVerMWG+lDETy1k3GyJWyZIkBwqIFtAEXqMVR9ErVqY6eSJkaMSVY6D6VKHBgRRNXVI+7RLWa1aqTI0A7VDhSE0QDEasESyW5tITPoP82NHCxMeNJSlcXHBApEsXJEpbHZqFBp26dOWyPfulKo+gEL7gaYZ3R2AwL6MYo1sHWTLlTW3EbdbcuRswMZ1cYTF3jkGBAQICADixZdxqzgKBccFkysm2cOCwSUtmq1QlK+R+twb2BZKmIdOoUViQ4IABAibAhkRf3RqZmUyJcixrdqRIjxoxWMAQM35za22w/uyZgSuXhwoSQNCAAjKc4o50AymziB4voIIKCR1sYMEEJeDwzW/ADYQMHkLA8YglkzBSyBl1MPMOhq0NFE0UiKxBBhVZsBGKNyeiaFAS7IgzDjk8toOhfQbp8MaQRA5Jx5FIHtmERkw2WVBAADs=''',
'zoom_to_rect_active':'''R0lGODlhGAAYAOfnAC5tfjVxgTZxgjhyhDlzgj10hj52hkJ5iKxmKEZ7ikZ9ikp9jUt+jkp/jJJyTKtuPU6Cjq5vP1KDkFCEj1KEklSFklOGkZx5Q650Q7B0Q1aHkrJ0Qq51SViJlFmJlLB4S7F6Rad7WrN6SKB9Y12MlmCOlmCOmGKPmbGAWpaGca6CVreBUWmUnK+FZbqHOmGYoGqWnW2Wn26Xn2WaoraKVGaco7KLbHGZonSZoWmcpIyUk7yMVmidpnWao2qeprqPUbyNYb2RRW6fp3qdnnSepnidpbKReXGgp22hqnmeprKSeXagqG6irHCiqnGiqnKiqr+WSMOWQ3mhsoGgp8KXRH6ip7+VaX6iqHakrH+iqYCiqYGiqnalrYahpoKiqIOiqIOiqYSiqoGjqoKjqoeipoaiqYOjqoiipoSjqHimroeiqnmmrcOaS4KkqsOaTsSaTXmnr3qnr4akqrWahHqosJCjq4OnrIamroqmqoinrsWeV36qssWeWX+qs5Gmq4mor5SmrYqprZKnr5Snro6prIqqroqqr4qqsYuqr46qrMahYZCqrMOgeY2rr8aiY8aiZY6sr5KrrpGrsY6ssJOrro6tr5Wsr5WssLWlmJKtssalbZKuspmss7WmnMambpWuspeusZOvs5qttJiussWngZqvspawtKGttpmws5uwtKGut8epeaGuuMapgKKuuseqepiytpmytqSvuqWvusKqmp6ytZqzt56ytqCytqawusesgraup6axuqexusitgqixup21uaqyuaqyuqWzwJ62uZ+2uqG3u7eyr6+0uaK4vKO5vKW5vbe1srW2uaa6vre2tae6vre3uKm7v7a4ubi4uK66xqq8wK+6xrC6xrG6xrK7xq2+wq6+wq+/xLLAxbLBxbXDx7bDx77CycDDycLEycfHysvKy////////////////////////////////////////////////////////////////////////////////////////////////////yH+EUNyZWF0ZWQgd2l0aCBHSU1QACwAAAAAGAAYAAAI/gCpCRxIsKDBgwgTKlwoEFOGCMwYGmSGYkWrHUokEpyzgdGrR24QFGyWq2SuXtMG7voARJcnPlF+PBgYbdYpMl++hDnDShg1GyJIrXL0BgqIFscENnOl5hCqUqAsXfoDaIQVX5r0UKHBoRNNWVMy1cJ1K9WoSIny6FChiE0QDEaeEeRVRpKTPX32xOHSxEeNJSlcXHAwhAsXJEpPGUq1xhu4b9ugGbNlyo6fELTMaTZHR2CuLp8YewsHWTLlSmawbdbcmVouMJZKYeHWjUGBAQICADiRJdtqzp69RAL1RJo1a86SAYO1CdIVbb9b5/qyiJIQZcsoLEhwwACBElugha9u/etMJEI5iBVLUqTHjRgsYIARv7l1NFV58MyAFctDBQkQNKCADKKQE91AwQxyxwuhhEJCBxpYMIEJOFzzG3AD/VIHEW0gMkkjhQQyhhzDlHNhawMhI4UgZohRhRZocFKNiScadIQ42GSjzY7jXFifQTykIeSQQsJh5JFGMqHRkkwWFBAAOw==''',
'previous':'''R0lGODlhGAAYAOe4AAAAAAEBAQMDAwUGBgoMDA8REREUFBQWFxgbHCMoKDI4OTU8PT5FR1JcXlReYV1pa2x5fHeGiXiGidd/MoSUmLePX9aJTNmLT4ydocOXVNmRVNyRVN6RUtmSW9yWXt2YVtGacOCYWsicfJWnqpiqrt2gcLyojdqjbOWhZdumfuipSeStaaS4vN6uh5G/ya+5uOuva+uweeizZey1VpHFz961l9+3l5PI0u+7WpfJ3qHI0fS8VJ3K0fO9VZ7K0u+6g6DK05/L06PK0qjJ0KTK0rHHy6HL1KTL0qTL06nK0KHM1aXL1KrK0KLM1KfL0/TAXqTM1KXM0qTM1fTAYaLN1anL1fXBYOLApajN1LTM1rbM0aTR16fQ2bfN0q3Q1PbFbarR2arR2vbGb7XP1qvS26zT2LnQ1LnQ2KzU27fR27LT163U3bnR2a3V2bLU1/fJea7V27TU1/TIl/fKfLDW2/fLfrHX27HX3LXW3bfW2bfW2rHY2+LOvrrX27rX3LbY3/jOiL/X28DX3L/X4PjPirrZ3uLQw7fa373Z3cHY4bjb4L7a3sLa37vc4fbRocDb38nY47/c4MnZ5cnZ5vnTl8Lc4crZ6fjToMTc4fnUmcXd4s3b6cbe4s7b6cbe4/LVwc/c6dDc6fnXouPZ0dDd6dHd6frYo9Ld6NTe587g8NXf6eXf29vh5+Xi3+Xj4uLk5+Xl5uTm5+bm5tro+Nvo99zp993p+N7q+P///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////yH+EUNyZWF0ZWQgd2l0aCBHSU1QACwAAAAAGAAYAAAI/gBlCRxIsKDBgwgTKlwokM+GC60YGmxVAsUlGDYkErzCQU6mOlMmFHwVCpTJUrEGjvIQQxQhMTtkWBgIqxOkJEeILGEySZWsFiEcUZpjBceHFKsEvrJUZU2kR4j6+CFzRsQPU4C+9FjRwRDNTTr+cPJEdpEeN2FenHjzZIaGGq4IknKCx5Mgsp4YYbBj54YJFRkqGHwFCU0lsppINCgAwFOjLWNAfDoYakghTyweIABwwAEJx3uk1EIICkmfBAAMNBiBCa8nO0FskRaShwAABRG0uD5kx8ct0kTiBKLAYICABRK6kLVj5PfBU0z0qCFrRoKCAAKWI3FuEJYkMF5cOHsqAmF5IoWo2HCxo0j8axe0Fp7KQoMKnDt02pRpgiWVRFY5pAGFEjwAEcUgs2ik4IIMNuggQgEBADs=''',
'next':'''R0lGODlhGAAYAOe5AAAAAAEBAQMDAwUGBgsNDQ8REREUFBMVFhgbHCMoKDI4OTU8PT5GSFFbXVJcXlReYV5qbGx5fHeGiXiGiXmHitd/MoSUmLePX9aJTNmLT4ydocOXVNmRVNyRVN6RUtmSW9yWXt2YVtGacOCYWsicfJWnqpiqrpirr92gcLyojdqjbOWhZdumfuipSeStaaS4vN6uh5G/ya+5uOuva+uweeizZey1VpHFz961l9+3l5PI0u+7WpfJ3qHI0fS8VJ3K0fO9VZ7K0u+6g6DK05/L06PK0qjJ0KTK0rHHy6HL1KTL0qTL06nK0KHM1aXL1KrK0KLM1KfL0/TAXqTM1KXM0qTM1fTAYaLN1anL1fXBYOLApajN1LTM1rbM0aTR16fQ2bfN0q3Q1PbFbarR2arR2vbGb7XP1qvS26zT2LnQ1LnQ2KzU27fR27LT163U3bnR2a3V2bLU1/fJea7V27TU1/TIl/fKfLDW27zT2PfLfrHX27HX3LXW3bfW2bfW2rHY2+LOvrrX27rX3LbY3/jOiL/X27/X4PjPirrZ3uLQw7fa373Z3cHY4bjb4L7a3rvc4fbRocDb38nY47/c4MnZ5cnZ5vnTl8Lc4crZ6fjToMTc4fnUmc3b6cbe4s7b6cbe4/LVwc/c6dDc6fnXouPZ0dDd6dHd6frYo9Ld6NTe587g8NXf6eXf29vh5+Xi3+Xj4uLk5+Xl5uTm5+bm5tro+Nvo99zp993p+N7q+P///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////yH+EUNyZWF0ZWQgd2l0aCBHSU1QACwAAAAAGAAYAAAI/gBnCRxIsKDBgwgTKlwoEFCHDK4YGnSFYkWmGTkkEtTioc6mPFYqFIQlKpRJU7IGkgJBY9ShMj5qYBgYy5MkJkqOOHlSadUsGCMgWbKTZUcIFqwEwsKExc2kSIsCCTqjhoSQU4TEAHHxIRFNTj0GdfpE1pGfOGRkqJAjxQYHHK8IlorCRwMesnj16NGRosWGCwZhSVpzCUABByfwfnrkxYwIUAdFGUH0ycSDAwAQQHix+E8VWwhDLQmEV1OJBgYAJPikh8it0EX6KO4iQQEAAor0BMEV+gidT2AmLBAwgIGFQqyT8D6I6omfNgICKKCQRrGeJcsNxqI0JkwEJIrzMzJSmOrNFz2NwrOOUWshKi43rszZcwcOGihbVElsxYPNlCY/DEGFIbRoZOCBCCaoIEIBAQA7''',
'subplots':'''R0lGODlhGAAYAMZHAA5efg9efw5ffg9ffxBgfxFgfxFggBNhgBRhgRRigRRighVigRdkghhlgxplhBxlhBpmhBtmhBxnhR1ohR5ohiJriCNriCRsiSVsiSRtiShuiylviylwjCtwjCxxjS1xjS1yjS9zjjV3kTl5kzp5kz17lUSAmEiDmkuEnE2GnVaMoVeMoViNolmNolqOo1uPpF+RpWOUp2SUqGaVqGiXqm2arG6brXSer3WfsHmisnqispGwvZSyv524w6a9x6m/ybDEzMTQ1cjT2M7X29/i4+Lj5Obm5v///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////yH+EUNyZWF0ZWQgd2l0aCBHSU1QACwAAAAAGAAYAAAHsoBGgoOEhYaHhkUvRIiNhDkcMo6OPisAHz2Th0M0CwAAJUKahSATBQAGDhdFo4UiAButhySwsoYmtaNDFYW4sYQZooU7GAc1x8cWARQ1NzrPCCE8g0UpDJ/Y2AnZ2Q8sjDEqDdza5J8FKJJGQDYjBDDw8BACETUz8QcaMECFQB2uuQZ5+NHq1S9bgwwiJKRwoQcJBgAMULDK1hAX2EwEWejjxKVMC43g0KAuZJEWjELKCgQAOw==''',
'filesave':'''R0lGODlhGAAYAOfWABgqQwNAYQVBYgBDaQdCYwBEagBFayU9VgBGbAtEZR5AWABHbSY+VwBIbg1FZh9BWQBJbwBKaw9GZxBGaAFLbBJHaQNMbRNIagBNcgBNcypDXCNFXSRGXh5IZQdObyVHXyZIYABRdwpPcABSeCdJYShJYgxQcQBTeSlKYw5RcipLZBBScwBXdjJKZCtMZQBYdxJTdCxNZh1RbRNUdQJZeDlLYS1OZwRZeTpMYhdVdzpNYwdaehhWeDVQYwlbewxcfA5dfSZYdBBefmhLOzJYditdeXlNN3JPR0dZb0xbbE5dbkxfdVFgcYNWRU5hdz5nf1Nic1BjeVVkdVJle1NmfFRnfVpoeltpe1Jtgjt0kUhxiUJ1jER2jUV3jlJ1iV51haBrPqFsP6JtQGN5ill8kKNuQWB7kKRvQaVwQqZxQ1x/k6dyRGh+j2CDl2uBknKFkGuGm6t9WHWIlKx+WXaJla1/Wm2Mm3eKlq6AW26NnLR/XGuOo6+BXLWAXXmMmGyPpLCCXXCPnraBXnqNmXuOmm6RpnKRoHOSoX2QnHSTon6RnXGUqXWUo3KVqnaVpICTn3eWpYGUoIKVoXyXoXmYqIOXonqZqXCcr4SYo3uaqoWZpHybq4aapYGcpn6drYmcp3+froqdqHaitnuisIueqYWgqnejt4yfqo2gq4eirI6hrYijrYmlr4qmsJGksIunsZKlsYyospSns4uqupWotJCrtZeqtpmsuJSvuZWwupuvupywu56yvpm0v5+zv6C0wKG1waS3w6W4xKa5xae6xqi7x6m8yKq+yau/y6zAzK7Czq/Dz7DE0LLG0rPH07jIzrbJ1rvL0r7P1dTZ29vg4/7//P///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////yH+EUNyZWF0ZWQgd2l0aCBHSU1QACwAAAAAGAAYAAAI/gCrCRxIsKBBgj8SKvTBsKHDhjRo/Kj2w0dEI2vSoNlYpmNHMWHCgBlCgcKLhE/MwGgiqI8eQHxi4qlDc86cOEdkwIFBw0ebaSJESYvmbJmxYsJ+7dJ1y5YsWFmIUCvSM0cQGKaiPWM2q9GiRIcK/dmTx04nLgMEOOg5A8YMUdCcKfMCAsSHuxw4bNighosBAxMijhhxwlRRY8OA7cJlixYsVqlOneIyAXDEFCJWjFp2zJiWB6BDKxhNpgsCyzREZDbVeVisTZYoOWJ06FCgPJNMoxaBYcSlXrlu1XrValUpUJ4yUYKUiMvpwDQsUMAgBEjCiNix83iRI8df6BQixVyosKDy3/Po00O/sKCC+fTw0QfmQf50/PvfafA47eAvgQAABiigAOgtQIN36GnwBRZTUHGFFVVMEYUTDBR4YHotGCMFAAD4QYsrqHxSSQ8WInheC8RAwaEcsKjyCSaP6FAihsIwweEbqoSiSSSI4DAjei74sgSHbpCyoyJ++HiegTSk58IuSnDIBieSKDIIHTWgB0FEMKDXwi1jIJGEHJUoQsgdb9SQwF9bTqRfdzl0YEMMLqiAQgkl2MXBASy88IJE1QQEADs='''
}
#IMAGE_DICT = {
#'initial':'''UDYKIyBDUkVBVE9SOiBHSU1QIFBOTSBGaWx0ZXIgVmVyc2lvbiAxLjEKMjQgMjQKMjU1Cubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5uLOvtyRVNmLT+Xi3+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5uXi392gcOWhZfjToOuva9+3l+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5uLApd6RUvTIl/nUmffLfvTAYdd/Mubm5ubm5ubm5ubm5ubm5uLk59Dc6c/c6c/c6c/c6dHd6eTm5+bm5ubm5ubm5uPZ0dyWXuuwefnXovjPivbGb/S8VOizZdaJTObm5ubm5ubm5uXl5s7b6cnY46nK0KTL0qTK0qXL1KrK0MnZ5tXf6ebm5t6uh+CYWvbRofnTl/fKfPXBYO+7Wt2YVtumfuXf2+bm5ubm5uLk58rZ6anL1a3U3b/c4MDb373Z3brX27rX3KvS27nQ2MicfO+6g/rYo/jOiPbFbfO9VeStadmSW+LQw+bm5ubm5ubm5uXl5s3b6aHI0bbY38be4sbe48be48be477a3rfW2rLU16rR2q+5uNqjbPfJefTAXuy1VtmRVN61l+Xj4ubm5ubm5ubm5ubm5tDd6afL07XW3cbe48DX3Dc9Pw8REQsNDcbe48be47HX27HX25PI0ryojeipScOXVLePX+bm5ubm5ubm5ubm5ubm5ubm5uLk58nY46zU28Lc4cbe48La35GipqO2ugsNDcbe48be48be47vc4aTR17XP1tGacPLVwebm5ubm5ubm5ubm5ubm5ubm5ubm5tDc6ajJ0LrZ3sbe48be48be48be47DFygsNDcbe48be48be47vc4bHY26TM1dzp9+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5s/c6aTL07rX28be48be48be48be47DFygsNDcbe48be48be48be47HX25/L093p+Obm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5s/c6aPK0rfW2cbe48be48be48be47DFygsNDcbe48be48be47fa37HX257K0t7q+Obm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5s/c6aTK0rTU18be48be48be48be47DFygsNDcbe48be48be48be47HX26HL1N7q+Obm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5tLd6KrK0LfW2rLT18be48be48be47DFygsNDcbe48be48be48be47HX26TL097q+Obm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5uXl5snZ5arR2a3Q1Mbe48be4ygtLhIUFQAAABQXGElSVMbe48be47HX28HY4ebm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5tTe57nR2afQ2bHX27jb4Mbe48be48be48be48be48be47HX25G/ydvo9+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5tLd6LTM1pHFz6LN1a7V27HX3LDW263V2azT2KLM1KjN1M7g8Obm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5tvh55fJ3rfR26TM1KHM1Z3K0aDK06XM0r/X4Nro+Obm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5g==''',
#'move':'''UDYKIyBDUkVBVE9SOiBHSU1QIFBOTSBGaWx0ZXIgVmVyc2lvbiAxLjEKMjQgMjQKMjU1Cubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5luPpObm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5nmisilwjGSUqObm5ubm5ubm5ubm5ubm5ubm5ubm5uLOvtyRVNmLT+Xi3+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5qa9x1eMoQ5efi1xjZ24w+bm5ubm5ubm5ubm5uXi392gcOWhZfjToOuva9+3l+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5s7X22iXqhVigQ5efg5efj17lcjT2Obm5ubm5uLApd6RUvTIl/nUmffLfvTAYdd/Mubm5ubm5ubm5ubm5ubm5uLk59Dc6c/c6S1yjR1ohRFgfw5efhFggBplhCRsiePZ0dyWXuuwefnXovjPivbGb/S8VOizZdaJTObm5ubm5ubm5uXl5s7b6cnY46nK0KTL0qTK0qXL1DV3kQ5efilvi+bm5t6uh+CYWvbRofnTl/fKfPXBYO+7Wt2YVtumfuXf2+bm5ubm5uLk58rZ6anL1a3U3b/c4MDb373Z3brX2zp5kw5efilvi8icfO+6g/rYo/jOiPbFbfO9VeStadmSW+LQw+bm5ubm5ubm5uXl5s3b6aHI0bbY38be4sbe48be48be477a3kSAmA5efilvi6+5uNqjbPfJefTAXuy1VtmRVN61l+Xj4ubm5ubm5ubm5ubm5tDd6afL087X2yJriMDX3Mbe48be48be48be40SAmA5efilvi5PI0ryojeipScOXVLePXyRticjT2Obm5ubm5ubm5ubm5uLk58nY45GwvSVsiRNhgG2arG2arG2arG2arG2arCNriA9efx5ohm2arHSer3qisnqisnqishRhgS9zjpSyv+bm5ubm5ubm5tDc6U2GnRdkgg5efg5efg5efg5efg5efg5efg5efg5efhRigQ5efg5efg5efg5efg5efg5efg5efg5efhxlhFiNoubm5ubm5mOUp1aMoRhlgw5efg5efg5efg5efg5efg5efg5efg5efhRigQ5efg5efg5efg5efg5efg5efg5efg5efhFgf0uEnGSUqObm5s/c6bDEzG6brTl5kxBgf1+RpV+RpV+RpV+RpV+RpRpmhA5ffhtmhG2arGaVqF+RpV+RpV+RpRNhgChui1+RpbDEzObm5ubm5s/c6aTK0rTU17DEzCtwjMbe48be48be48be48be4zV3kQ5efilvi7HX26HL1N7q+Obm5ubm5ixxjam/yebm5ubm5ubm5ubm5tLd6KrK0LfW2rLT18be48be48be48be48be48be4zV3kQ5efilvi7HX26TL097q+Obm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5uXl5snZ5arR2a3Q1Mbe48be48be48be48be48be4zV3kQ5efilvi7HX28HY4ebm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5tTe57nR2afQ2bHX27jb4Mbe48be48be48be4zV3kQ5efilvi5G/ydvo9+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5tLd6LTM1pHFz6LN1a7V27HX3CxxjRxnhRFggA5efg9ffxRigiRsiebm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5tvh55fJ3rfR26TM1KHM1c7X21qOow5efg5efg5efkSAmMTQ1ebm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5qa9x0iDmg5efi1xjZ24w+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5nWfsChui2SUqObm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5lmNoubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5g==''',
#'move_active':'''UDYKIyBDUkVBVE9SOiBHSU1QIFBOTSBGaWx0ZXIgVmVyc2lvbiAxLjEKMjQgMjQKMjU1Cri4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uElyg7i4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uGGCjiFacFB2hri4uLi4uLi4uLi4uLi4uLi4uLi4uLWlmLB0Q65vP7e1sri4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uIWXn0ZwgQtLZSRacX6TnLi4uLi4uLi4uLi4uLe1srGAWreBUcapgLyMVrKSebi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uKWsr1N5iBFOZwtLZQtLZTFid6Cprbi4uLi4uLWahLJ0QsOgeceqesaiZcOaTqxmKLi4uLi4uLi4uLi4uLi4uLW2uaawuqawuiRbcRdTag5NZgtLZQ5NZhVRah1Wbraup7B4S7yNYcesgsambsWeWcOWQ7qPUatuPbi4uLi4uLi4uLe3uKWvuqGttoeipoOiqIOiqISiqipfdAtLZSFZb7i4uLKLbLN6SMWngcepecaiY8SaTb+WSLF6Ra+FZbeyr7i4uLi4uLW2uaKuuoeiqoqqsZmws5qvspeusZWsry5hdgtLZSFZb6B9Y7+VacitgsalbcWeV8KXRLaKVK51SbWmnLi4uLi4uLi4uLe3uKSvuoGgp5Ktsp6ytZ6ytp6ytp6ytpiusjZmegtLZSFZb4yUk66CVsahYcOaS72RRa50Q7KRebe2tbi4uLi4uLi4uLi4uKaxuoaiqaWsrxtWbZqssJ6ytp6ytp6ytp6ytjZmegtLZSFZb3agqJaGcbqHOpx5Q5JyTB1XbqCprbi4uLi4uLi4uLi4uLW2uaGttnSNlx5Wbg9OZld7ild7ild7ild7ild7ihxWbQxLZhhTa1d7il1+jGKCjmKCjmKCjhBOZyZccnaOmbi4uLi4uLi4uKawuj5rfhJQaAtLZQtLZQtLZQtLZQtLZQtLZQtLZQtLZRBOZwtLZQtLZQtLZQtLZQtLZQtLZQtLZQtLZRZRakZxgri4uLi4uE92hkVwgRNRaQtLZQtLZQtLZQtLZQtLZQtLZQtLZQtLZRBOZwtLZQtLZQtLZQtLZQtLZQtLZQtLZQtLZQ5NZjxqfVB2hri4uKawuo2do1h8ii5hdg1NZkx0hEx0hEx0hEx0hEx0hBVSagtMZRZSald7ilJ3hkx0hEx0hEx0hA9OZiBYb0x0hI2do7i4uLi4uKawuoOiqJCqrI2doyJacJ6ytp6ytp6ytp6ytp6ytipfdAtLZSFZb46sr4GiqrK7xri4uLi4uCNacYeZobi4uLi4uLi4uLi4uKixuoiippKrro6prJ6ytp6ytp6ytp6ytp6ytp6ytipfdAtLZSFZb46sr4OiqbK7xri4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLe3uKGut4inroqmqp6ytp6ytp6ytp6ytp6ytp6ytipfdAtLZSFZb46sr5qttLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uKqyuZSnroamro6sr5Ovs56ytp6ytp6ytp6ytipfdAtLZSFZb3SZoa+6xri4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uKixupCjq3SepoKkqouqr46ssCNacRZSag5NZgtLZQxMZhBOaB1Wbri4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uK+0uXmhspKnr4OjqoGjqqWsr0hyggtLZQtLZQtLZTZmep2mqri4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uIWXnzppewtLZSRacX6TnLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uF5/jSBYb1B2hri4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uEdxgri4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uA==''',
#'zoom_to_rect':'''UDYKMjQgMjQKMjU1Cubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5uLOvtyRVNmLT+Xi3+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5uXi392gcOWhZfjToOuva9+3l+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5uLApd6RUvTIl/nUmffLfvTAYdd/Mubm5ubm5ubm5ubm5ubm5uLk59Dc6c/c6c/c6c/c6dHd6eTm5+bm5ubm5ubm5uPZ0dyWXuuwefnXovjPivbGb/S8VOizZdaJTObm5ubm5ubm5uXl5s7b6cnY46nK0KTL0qTK0qXL1KrK0MnZ5tXf6ebm5t6uh+CYWvbRofnTl/fKfPXBYO+7Wt2YVtumfuXf2+bm5ubm5uLk58rZ6anL1a3U3b/c4MDb373Z3brX27rX3KvS27nQ2MicfO+6g/rYo/jOiPbFbfO9VeStadmSW+LQw+bm5ubm5ubm5uXl5s3b6aHI0bbY38be4sjf48be48Lc4b7a3rfW2rLU16rR2q+5uNqjbPfJefTAXuy1VtmRVN61l+Xj4ubm5ubm5ubm5ubm5tDd6afL07XW3Y3L1Z3U3p/V4J7U3pnR25PO2IzK1ITFz4DDzJPI0ryojeipScOXVLePX5nExZTO2JTO2IjJ1ebm5ubm5uLk58nY46zU28Lc4ZfP2N7w9uL0+d/x9tjt8tHp7snl6sHg5bvc4aTR17XP1tGacPLVwf78/v78/v78/v78/pjS3Obm5ubm5tDc6ajJ0LrZ3sLc4ZfP2N7w9uP0+d/x9tjt8tHp7snl6sHg5bvc4bHY26TM1dzp9/78/v78/v78/v78/v78/pjS3Obm5ubm5s/c6aTL07rX28Hb35TN19nu89vv9V6dsUyRqEaOpUSNo0KNoTqInXqzv5/L093p+P78/v78/v78/v78/v78/pjS3Obm5ubm5s/c6aPK0rfW2b3Z3Y/K1NPq79Xr8NTr8NDp7svm68Ti577e47fa37HX257K0t7q+P78/v78/v78/v78/v78/pjS3Obm5ubm5s/c6aTK0rTU17jW2orH0czn687n7Gelt1ycsFiarVOXqk6Up0eQo3iyvKHL1N7q+P78/v78/v78/v78/v78/pjS3Obm5ubm5tLd6KrK0LfW2rLT14PDzcXj58fj6JfF0JbEzpLBzI2/yoi8x4O5w4W8xKTL097q+P78/v78/v78/v78/v78/pjS3Obm5ubm5uXl5snZ5arR2a3Q1H7Byr7f47/f5G+ruWmmtmaktGKisV2fr1icrIq9x8HY4fL1+/78/v78/v78/v78/v78/pjS3Obm5ubm5ubm5tTe57nR2afQ2Xm+yLjb4Ljb4HSvvG6ruWypt2intWSls3ixvpG/ydvo9/78/v78/v78/v78/v78/v78/pjS3Obm5ubm5ubm5ubm5tLd6LTM1pHFz6LN1a7V27HX3LDW263V2azT2KLM1KjN1M7g8Pn5/f78/v78/v78/v78/v78/v78/pjS3Obm5ubm5ubm5ubm5ubm5tvh55fJ3rfR26TM1KHM1Z3K0aDK06XM0r/X4Nro+Pn5/f78/v78/v78/v78/v78/v78/v78/pjS3Obm5ubm5ubm5ubm5ubm5ubm5o3I0e3z+9zp993p+N7q+N7q+N7q+PD0+/78/v78/v78/v78/v78/v78/v78/v78/v78/pjS3Obm5ubm5ubm5ubm5ubm5ubm5oLEz5bQ2pbQ2pbQ2pbQ2pbQ2pbQ2pbQ2pbQ2pfR25fR25fR25fR25fR25fR25fR25fR24rL1+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5g==''',
#'zoom_to_rect_active':'''UDYKIyBDUkVBVE9SOiBHSU1QIFBOTSBGaWx0ZXIgVmVyc2lvbiAxLjEKMjQgMjQKMjU1Cri4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLWlmLB0Q65vP7e1sri4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLe1srGAWreBUcapgLyMVrKSebi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLWahLJ0QsOgeceqesaiZcOaTqxmKLi4uLi4uLi4uLi4uLi4uLW2uaawuqawuqawuqawuqexura4ubi4uLi4uLi4uLaup7B4S7yNYcesgsambsWeWcOWQ7qPUatuPbi4uLi4uLi4uLe3uKWvuqGttoeipoOiqIOiqISiqoiipqGuuKqyuri4uLKLbLN6SMWngcepecaiY8SaTb+WSLF6Ra+FZbeyr7i4uLi4uLW2uaKuuoeiqoqqsZmws5qvspeusZWsr5WssImor5SmraB9Y7+VacitgsalbcWeV8KXRLaKVK51SbWmnLi4uLi4uLi4uLe3uKSvuoGgp5Ktsp6ytaCytp6ytpuwtJiuspKrro6qrIinroyUk66CVsahYcOaS72RRa50Q7KRebe2tbi4uLi4uLi4uLi4uKaxuoaiqZGrsXGiqn6qsn+qs36qsnqnr3alrXCiqmqepmaco3agqJaGcbqHOpx5Q5JyTHqdnnalrXalrW2hqri4uLi4uLW2uaGttoqqr5uwtHmmrbLAxbXDx7LBxa2+wqe6vqG3u5qzt5awtIOnrJGmq6d7WsKqmsvKy8vKy8vKy8vKy3qosLi4uLi4uKawuoahppWuspuwtHmmrbLAxbbDx7LBxa2+wqe6vqG3u5qzt5awtI6tr4OjqrC6xsvKy8vKy8vKy8vKy8vKy3qosLi4uLi4uKawuoOiqZWsr5qvsnakrK6+wq+/xEt+jj10hjhyhDZxgjVxgS5tfmKPmX+iqbG6xsvKy8vKy8vKy8vKy8vKy3qosLi4uLi4uKawuoKiqJKrrpeusXKiqqm7v6q8wKq8wKa6vqK4vJ21uZiytpKuso6sr36iqLK7xsvKy8vKy8vKy8vKy8vKy3qosLi4uLi4uKawuoOiqJCqrJOrrm6fp6O5vKW5vVKEkkp9jUZ7ikJ5iD52hjlzgmCOloGiqrK7xsvKy8vKy8vKy8vKy8vKy3qosLi4uLi4uKixuoiippKrro6prGmcpJ62uZ+2unmepnidpXWao3GZom2Wn2mUnGqWnYOiqbK7xsvKy8vKy8vKy8vKy8vKy3qosLi4uLi4uLe3uKGut4inroqmqmWaopiytpmytlmJlFSFklKDkE6Cjkp/jEZ9im6Xn5qttMLEycvKy8vKy8vKy8vKy8vKy3qosLi4uLi4uLi4uKqyuZSnroamrmGYoJOvs5Ovs12MlliJlFaHklOGkVCEj2COmHSZoa+6xsvKy8vKy8vKy8vKy8vKy8vKy3qosLi4uLi4uLi4uLi4uKixupCjq3SepoKkqouqr46ssI2rr4qqroqprYKjqoakqqWzwMfHysvKy8vKy8vKy8vKy8vKy8vKy3qosLi4uLi4uLi4uLi4uLi4uK+0uXmhspKnr4OjqoGjqn6ip4CiqYSjqJmss666xsfHysvKy8vKy8vKy8vKy8vKy8vKy8vKy3qosLi4uLi4uLi4uLi4uLi4uLi4uHGgp77CybC6xrG6xrK7xrK7xrK7xsDDycvKy8vKy8vKy8vKy8vKy8vKy8vKy8vKy8vKy3qosLi4uLi4uLi4uLi4uLi4uLi4uGidpnimrnimrnimrnimrnimrnimrnimrnimrnmnr3mnr3mnr3mnr3mnr3mnr3mnr3mnr26irLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uLi4uA==''',
#'previous':'''UDYKIyBDUkVBVE9SOiBHSU1QIFBOTSBGaWx0ZXIgVmVyc2lvbiAxLjEKMjQgMjQKMjU1Cubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5uLOvtyRVNmLT+Xi3+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5uXi392gcOWhZfjToOuva9+3l+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5uLApd6RUvTIl/nUmffLfvTAYdd/Mubm5ubm5ubm5ubm5ubm5uLk59Dc6c/c6c/c6c/c6dHd6eTm5+bm5ubm5ubm5uPZ0dyWXuuwefnXovjPivbGb/S8VOizZdaJTObm5ubm5ubm5uXl5s7b6cnY46nK0KTL0qTK0qXL1KrK0MnZ5tXf6ebm5t6uh+CYWvbRofnTl/fKfPXBYO+7Wt2YVtumfuXf2+bm5ubm5uLk58rZ6anL1a3U3b/c4MDb373Z3brX27rX3KvS27nQ2MicfO+6g/rYo/jOiPbFbfO9VeStadmSW+LQw+bm5ubm5ubm5uXl5s3b6aHI0bbY38be4sbe48be48be477a3rfW2rLU16rR2q+5uNqjbPfJefTAXuy1VtmRVN61l+Xj4ubm5ubm5ubm5ubm5tDd6afL07XW3cbe48DX3Mbe48be48be48La34ydobHX27HX25PI0ryojeipScOXVLePX+bm5ubm5ubm5ubm5ubm5ubm5uLk58nY46zU28Lc4cbe48be48Xd4piqrlJcXg8REQAAAMbe47vc4aTR17XP1tGacPLVwebm5ubm5ubm5ubm5ubm5ubm5ubm5tDc6ajJ0LrZ3sbe46S4vF1paxgbHAAAABQWF1ReYZiqrsbe47vc4bHY26TM1dzp9+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5s/c6aTL07rX2yMoKAAAABEUFFJcXpWnqsTc4cbe48be48be48be47HX25/L093p+Obm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5s/c6aPK0rfW2QoMDAAAADI4OXeGibbM0cbe48be48be48be47fa37HX257K0t7q+Obm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5s/c6aTK0rTU17/X24SUmD5FRwUGBgMDAzU8PXiGibfN0sbe48be47HX26HL1N7q+Obm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5tLd6KrK0LfW2rLT18be48be47nQ1HiGiTI4OQEBAQMDA8be48be47HX26TL097q+Obm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5uXl5snZ5arR2a3Q1Mbe48be48be48be48be47HHy2x5fMbe48be47HX28HY4ebm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5tTe57nR2afQ2bHX27jb4Mbe48be48be48be48be48be47HX25G/ydvo9+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5tLd6LTM1pHFz6LN1a7V27HX3LDW263V2azT2KLM1KjN1M7g8Obm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5tvh55fJ3rfR26TM1KHM1Z3K0aDK06XM0r/X4Nro+Obm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5g==''',
#'next':'''UDYKIyBDUkVBVE9SOiBHSU1QIFBOTSBGaWx0ZXIgVmVyc2lvbiAxLjEKMjQgMjQKMjU1Cubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5uLOvtyRVNmLT+Xi3+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5uXi392gcOWhZfjToOuva9+3l+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5uLApd6RUvTIl/nUmffLfvTAYdd/Mubm5ubm5ubm5ubm5ubm5uLk59Dc6c/c6c/c6c/c6dHd6eTm5+bm5ubm5ubm5uPZ0dyWXuuwefnXovjPivbGb/S8VOizZdaJTObm5ubm5ubm5uXl5s7b6cnY46nK0KTL0qTK0qXL1KrK0MnZ5tXf6ebm5t6uh+CYWvbRofnTl/fKfPXBYO+7Wt2YVtumfuXf2+bm5ubm5uLk58rZ6anL1a3U3b/c4MDb373Z3brX27rX3KvS27nQ2MicfO+6g/rYo/jOiPbFbfO9VeStadmSW+LQw+bm5ubm5ubm5uXl5s3b6aHI0bbY38be4sbe48be48be477a3rfW2rLU16rR2q+5uNqjbPfJefTAXuy1VtmRVN61l+Xj4ubm5ubm5ubm5ubm5tDd6afL07XW3YydobzT2Mbe48be48be48be48be47HX27HX25PI0ryojeipScOXVLePX+bm5ubm5ubm5ubm5ubm5ubm5uLk58nY46zU28Lc4QAAAA8REVJcXpirr8be48be48be48be47vc4aTR17XP1tGacPLVwebm5ubm5ubm5ubm5ubm5ubm5ubm5tDc6ajJ0LrZ3sbe45iqrlReYRMVFgAAABgbHF5qbKS4vMbe47vc4bHY26TM1dzp9+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5s/c6aTL07rX28be48be48be48Tc4ZWnqlFbXREUFAAAACMoKMbe47HX25/L093p+Obm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5s/c6aPK0rfW2cbe48be48be48be47bM0XeGiTI4OQAAAAsNDbfa37HX257K0t7q+Obm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5s/c6aTK0rTU18be47fN0niGiTU8PQMDAwUGBj5GSISUmL/X28be47HX26HL1N7q+Obm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5tLd6KrK0LfW2rLT1wMDAwEBATI4OXmHirnQ1Mbe48be48be48be47HX26TL097q+Obm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5uXl5snZ5arR2a3Q1Gx5fLHHy8be48be48be48be48be48be48be47HX28HY4ebm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5tTe57nR2afQ2bHX27jb4Mbe48be48be48be48be48be47HX25G/ydvo9+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5tLd6LTM1pHFz6LN1a7V27HX3LDW263V2azT2KLM1KjN1M7g8Obm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5tvh55fJ3rfR26TM1KHM1Z3K0aDK06XM0r/X4Nro+Obm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5g==''',
#'subplots':'''UDYKMjQgMjQKMjU1Cubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5uLj5FuPpN/i4+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5nmisilwjGSUqObm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5qa9x1eMoQ5efi1xjZ24w+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5s7X22iXqhVigQ5efg5efj17lcjT2Obm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5i1yjR1ohRFgfw5efhFggBplhCRsieLj5Obm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5jV3kQ5efilvi+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5jp5kw5efilvi+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5kSAmA5efilvi+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5s7X2yJriObm5ubm5ubm5ubm5ubm5kSAmA5efilvi+bm5ubm5ubm5ubm5ubm5iRticjT2Obm5ubm5ubm5ubm5ubm5ubm5pGwvSVsiRNhgG2arG2arG2arG2arG2arCNriA9efx5ohm2arHSer3qisnqisnqishRhgS9zjpSyv+bm5ubm5ubm5uLj5E2GnRdkgg5efg5efg5efg5efg5efg5efg5efg5efhRigQ5efg5efg5efg5efg5efg5efg5efg5efhxlhFiNot/i4+bm5mOUp1aMoRhlgw5efg5efg5efg5efg5efg5efg5efg5efhRigQ5efg5efg5efg5efg5efg5efg5efg5efhFgf0uEnGSUqObm5ubm5rDEzG6brTl5kxBgf1+RpV+RpV+RpV+RpV+RpRpmhA5ffhtmhG2arGaVqF+RpV+RpV+RpRNhgChui1+RpbDEzObm5ubm5ubm5ubm5ubm5rDEzCtwjObm5ubm5ubm5ubm5ubm5jV3kQ5efilvi+bm5ubm5ubm5ubm5ubm5ixxjam/yebm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5jV3kQ5efilvi+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5jV3kQ5efilvi+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5jV3kQ5efilvi+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ixxjRxnhRFggA5efg9ffxRigiRsieLj5Obm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5s7X21qOow5efg5efg5efkSAmMTQ1ebm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5qa9x0iDmg5efi1xjZ24w+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5nWfsChui2SUqObm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5uLj5FmNot/i4+bm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5g==''',
#'filesave':'''UDYKMjQgMjQKMjU1Cv///////////////////////////////////////////////////////////////////////////////////////////////////wxcfQxcfQtcfAxcfQtcfQtcfQtbfQtbfAtbfAtbfAtbfAtbfAtbfAtbfAtbfAtbfAtbfAtbfAxbfAtaewtaewxcff///wxcfQtbfApYegtYeXpMOKZzQ6VyQ6RxQ6RxQ6RxQ6NwQaJvQqJuQaFuQaFuQaFtQKBsQJ9sP59rP2hKPQNLbgRLbglXeAxcfQxcfTpnf197jwpTdIRWRrSCXbOBXbKAXLGBXbGAXLGAXbCAXK9/XK5/XK5+XK5+W619XK1+XKx8W3BQRx9RbmyGmSBScQtZegxbfGWCldXY3AtQcHejtLzP17nL1LTH0bDDzay/yam7xqS4w6G0v52yvZuvu5ituJaqtpOntJClsDx0jjBZdd3g5DNcdwpZegtbfBhWdChYdAhTdHmjtbrM1bbI0bHEzo2quHOWqHKVp3GUpXGSom+Ro26PoW2OoWyNnmuMnYOcqT51jQJDZwRBYwNFaApZegtbfAhUdgdTdQlUdnqitLbJ07LG0K7Cy092iylHXyhHXyhGXydGXydGXidFXidFXiZFXSZFXVx/kkB2jQFFaQBEaAFHagpZegtaewdSdQdRdAhTdnujtbPG0K/DzKy/yaS5xZ+1wJuwvZauupSrt5KptY+ms42kr4qirYifrImfq0F2jwFGaQBEaAFHagpZegtaewZRcwZQcgdSdHyis7DDzKzAyam9xkpxhyBBWyBBWyBBWx9BWx9BWx9BWh9AWh9AWh9AWlh8jkN3jgFGagBEaAFHagpZegtZewZPcgVOcQdSdHyjta3Ayam9xqa5w42ntHybqXqZqXiYp3eVpXWUo3SSoXORn3GPnXCNnX6WokV3jgFGawBEaAFHagpZegpZegVOcQRNcAZRdHKcrpu0wZewvZWtupGqt42ms4ulsomjsIagrYOeq4Gdqn+ap32YpXuVo3iToEF1jQFGagBEaAFHagpZegpZegRMbwRLbgRMcBlefhldfRdcfBdcexVaehVaeRRZeRRZeRRZeRNZeRNYeBNYeBJXeBJXdxJWdxBWdgFFaQBEaAFHagpZegpYegNLbgNKbQNJbAJIawFHawFGaQBEaABEaABEaABEaABEaABEaABEaABEaABEaABEaABEaABEaABEaABEaABEaAFHagpZegpYeQNJbAJIbAJIawFHagFFaQBEaABEaABEaABEaABEaABEaABEaABEaABEaABEaABEaABEaABEaABEaABEaABEaAFHagpYegpXeQJIawFHawFGagFFaQBEaABEaABEaABEaABEaABEaABEaABEaABEaABEaABEaABEaABEaABEaABEaABEaABEaAFHagpYeglXeQFGagFFaQBFaABEaABEaANCZQRAYwRAYwRAYwRAYwRAYwRAYwRAYwRAYwRAYwRAYwRBYwBEaABEaABEaABEaAFHawpYeglWeABEaABEaABEaABEaABEaChDXV90hldsf1NleVJmeVZqfVVpfFNnelFleE9jdk1hdCk+VgBEaABEaABEaABEaAFHawpYeglWeABEaABEaABEaABEaABEaC9KY6u+yFRjdh0qQh0qQnqMmZWps5Gkro2gqomcpoSXoTdPZABEaABEaABEaABEaAJHawpYeglWeABEaABEaABEaABEaABEaC9KY6e7xVNhdB0qQh0qQnaIlpKlr46hq4mcp4WYooGUnjdOYwBEaABEaABEaABEaAJHawpYeglWeABEaABEaABEaABEaABEaC9KY6S3wVFgcx0qQh0qQnOFk46hrIqdp4aZo4GUn32QmjdNYwBEaABEaABEaABEaAJHawpYeglWeABEaABEaABEaABEaABEaDBLY6CzvVBech0qQh0qQnCBkIueqIaZpIKVn36Rm3mNlzZNYgBEaABEaABEaABEaAJHawpYegpYeQBEaABEaABEaABEaABEaDBLY52wuk5dcB0qQh0qQm1+jYeapIOWoH6SnHqNl3aJkzZMYgBEaABEaABEaABEaANJbQpZegtaewhTdQBEaABEaABEaABEaDBKY5mstmh5iEpZbUtbbXaIlYOXoX+SnHuOmHeKlHKFkDdMYQBEZwBEaABEaANJbQpZegxcff///wtZewpXeQlWeAlWeAlWeCJIYzBOZi5NZS1MZC1LYyxKYitJYSpIYClHXyhGXidFXSE+VwlWdwlXeAlXeApZegxcff///w=='''
#}

#==============================================================================
if __name__ == "__main__":
    if False:
        import sys
        if len(sys.argv) == 2:
            CTK.FILE = sys.argv[1]
            try: CTK.t = C.convertFile2PyTree(CTK.FILE)
            except: pass

        # Main window
        (win, menu, file, tools) = CTK.minimal('tkPlotXY '+C.__version__)

        createApp(win); showApp()

        # - Main loop -
        win.mainloop()
    else:
        main(data)
