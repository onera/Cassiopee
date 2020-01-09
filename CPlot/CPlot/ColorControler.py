# Fenetre de dialogue pour le choix de la couleur
# Attention : dependance aux colormaps de matplotlib
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import colorsys
import matplotlib.pyplot as plt

# local NUM_COLORS for colormap
NUM_COLORS = 8
# local color map set (cf: http://matplotlib.org/examples/color/colormaps_reference.html)
#COLOR_MAP = 'Set2'
# COLOR_MAP = 'rainbow'
COLOR_MAP = 'Set1'
COLOR_LIST_HISTORY = []

################################################################################
# Note :
# ------
# RGB : convention entre 0 et 255 pour chaque composant
# HSV : convention ?
# HTML : #RGB avec R,G,B en hexadecimal
# colorsys renvoie entre 0 et 1
#
################################################################################
class ColorControler(object):
    def __init__(self,parent,color="#3a65c6",colormapName='Set1',nbColorMap=8,colorHistoryList=None ):
        self.parent = parent
        self.color = color
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
        for ind in range(nbColorMap): classicalColorLblFrame.columnconfigure(ind,weight=1)
        classicalColorLblFrame.grid(row=0,column=0,sticky="NSEW")
        #
        for ind in range(len(self.colormapList)):
            htmlcolor = self.colormapList[ind]
            hexa = htmlcolor.lstrip('#')
            (r,g,b) = tuple(int(hexa[i:i+2], 16) for i in (0, 2 ,4))
            (h,s,v) = colorsys.rgb_to_hsv(float(r)/255.,float(g)/255.,float(b)/255.)
            squareColor = SquareColor(classicalColorLblFrame,self,h,s,v)
            squareColor.grid(row=0,column=ind,sticky='NS')

        lastUsedColorLblFrame = TTK.LabelFrame(pickUpColorFrame,text="Last colors")
        lastUsedColorLblFrame.rowconfigure(0,weight=1)
        for ind in range(len(colorHistoryList)):
            lastUsedColorLblFrame.columnconfigure(ind, weight=1)
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
        if h < 255: res = 255
        elif h < 510: res = 510 - h
        elif h < 1020: res = 0
        elif h < 1275: res = h-1020
        else: res = 255
        return res

    def getGvalue(self,h):
        if h < 255: res = h
        elif h < 765: res = 255
        elif h < 1020: res = 1020 - h
        else: res = 0
        return res

    def getBvalue(self, h):
        if h < 510: res = 0
        elif h < 765: res = h - 510
        elif h < 1275: res = 255
        else: res = 1530 - h
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
                if not re.match(re_pattern,val): return
                if val =='': return
                val = int(val)
                if val < 0 or val > 255: return
        elif var == 'html':
            html = self.html_new.get()

            re_pattern = re_pattern = re.compile('''^\#[0-9A-Fa-f]{6}$''')
            # Check the input to be an html color code #AABBCC where AA,BB,CC are hexadecimal number between 0 and 255
            if not re.match(re_pattern,html):
                return
            # Convert it to rgb
            hexa = html.lstrip('#')
            (r,g,b) = tuple(int(hexa[i:i+2], 16) for i in (0, 2 ,4))
        else: return

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
        if h is not None: self.hvalue = h
        if s is not None: self.svalue = s
        if v is not None: self.vvalue = v

    def left_bt_click(self,event):
        # Get the size of the window
        width = self.winfo_width()
        height = self.winfo_height()
        # Update the values of s and v according to the click event
        self.controler.turnOffFlag2ChangeColor()
        self.controler.updateHSV4Canvas(v=float(event.x)/float(width),s=float(event.y)/float(height))
        self.controler.turnOnFlag2ChangeColor()
        self._draw_lines()
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
        if self.find_withtag('gradient2'): self.delete(self.find_withtag('gradient2'))
        #
        nValue = 50
        #
        width = self.winfo_width()
        height = self.winfo_height()
        #
        limit = height
        i = 0.; j = 0.
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
                color = "#%02x%02x%02x" % (nr,ng,nb)
                self.create_rectangle(j,i,j+width_step,i+height_step, tags=("gradient2",), fill=color,width=0)
                j += width_step
            i += height_step
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
        if h is not None: self.hvalue = h
        if s is not None: self.svalue = s
        if v is not None: self.vvalue = v

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
            ii = nbColor*i//limit

            nr = int(self.rgb_l[ii][0])
            ng = int(self.rgb_l[ii][1])
            nb = int(self.rgb_l[ii][2])

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
        if h is not None: self.hvalue = h
        if s is not None: self.svalue = s
        if v is not None: self.vvalue = v
        self._draw_all()

    def _draw_all(self,event=None):
        self._draw_color(event)

    def _draw_color(self, event=None):
        '''Draw the color'''
        if self.find_withtag('color'):
            self.delete(self.find_withtag('color'))
        # width = self.winfo_width()
        width = 80; height = 40

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
