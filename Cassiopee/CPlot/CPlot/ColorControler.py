# Fenetre de dialogue pour le choix de la couleur
# Attention : dependance aux colormaps de matplotlib
try: import Tkinter as TK
except: import tkinter as TK
from . import Ttk as TTK
import colorsys
import tkPlotXY

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
    def __init__(self, parent, color="#3a65c6", colormapName='Set1', nbColorMap=8, colorHistoryList=None):

        if color == 'black' or color == 'Black': color = "#000000"

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

    def discretizedColor(self, h):
        # h goes from 0->360
        x = float(h)*1530./360.
        r = self.getRvalue(x)
        g = self.getGvalue(x)
        b = self.getBvalue(x)
        return (r,g,b)

    def getRvalue(self, h):
        if h < 255: res = 255
        elif h < 510: res = 510 - h
        elif h < 1020: res = 0
        elif h < 1275: res = h-1020
        else: res = 255
        return res

    def getGvalue(self, h):
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

    def createColorMap(self, colormapName, nbColorMap):
        colormapList = []
        try:
            import matplotlib.pyplot as plt
            cm = plt.get_cmap(colormapName)
            for ind in range(nbColorMap):
                rgb = cm(1.*(ind%nbColorMap)/nbColorMap)
                html = '#%02x%02x%02x' % (int(255*rgb[0]),int(255*rgb[1]),int(255*rgb[2]))
                colormapList.append(html)
        except:
            if colormapName == 'Set1':
                colormapList = ['#e41a1c', '#377eb7', '#4dae4a', '#994ea1', '#ff8100', '#fdfb32', '#a7572b', '#f481bd']
            elif colormapName == 'Set2':
                colormapList = ['#66c2a5', '#e9936a', '#a79bb1', '#c692c5', '#c5b289', '#c8d845', '#f7d34a', '#ddc198']
            elif colormapName == 'Set3':
                colormapList = ['#8dd3c7', '#e6e4c1', '#ec8d8a', '#91b1c3', '#d6c965', '#f4ced8', '#d0bfd1', '#c6c6c2']
            elif colormapName == 'Greys':
                colormapList = ['#ffffff', '#efefef', '#d8d8d8', '#bcbcbc', '#959595', '#727272', '#505050', '#232323']
            elif colormapName == 'Purples':
                colormapList = ['#fcfbfd', '#eeecf4', '#d9d9ea', '#bbbcdb', '#9d99c7', '#7f7cb9', '#6950a2', '#53258e']
            elif colormapName == 'Blues':
                colormapList = ['#f7fbff', '#ddeaf6', '#c5daee', '#9dc9e0', '#6aadd5', '#4191c5', '#2070b4', '#08509a']
            elif colormpaName == 'Greens':
                colormapList = ['#f7fcf5', '#e4f4df', '#c6e8bf', '#a0d89a', '#73c375', '#40aa5c', '#228a44', '#006b2b']
            elif colormapName == 'Oranges':
                colormapList = ['#fff5eb', '#fde5cd', '#fdcfa1', '#fdad6a', '#fc8c3b', '#f06812', '#d74701', '#a43503']
            elif colormapName == 'Reds':
                colormapList = ['#fff5f0', '#fddfd1', '#fcbaa0', '#fb9171', '#fa6949', '#ee3a2b', '#ca171c', '#a30e14']
            elif colormapName == 'YlOrBr':
                colormapList = ['#ffffe5', '#fef6bb', '#fee290', '#fec34e', '#fd9828', '#eb6f13', '#ca4b02', '#973304']
            elif colormapName == 'YlOrRd':
                colormapList = ['#ffffcc', '#feec9f', '#fed875', '#fdb14b', '#fc8c3b', '#fb4c29', '#e2191c', '#bb0026']
            elif colormapName == 'PuRd':
                colormapList = ['#f7f4f9', '#e6e0ee', '#d3b8d9', '#c993c6', '#df64af', '#e62888', '#cc1155', '#960042']
            elif colormapName == 'RdPu':
                colormapList = ['#fff7f3', '#fcdfdc', '#fbc4bf', '#f99eb4', '#f667a0', '#dc3396', '#ac017d', '#780076']
            elif colormapName == 'BuPu':
                colormapList = ['#f7fcfd', '#dfebf3', '#bed2e5', '#9dbbd9', '#8c95c5', '#8b6ab0', '#873f9c', '#7f0e7a']
            elif colormapName == 'GnBu':
                colormapList = ['#f7fcf0', '#dff2da', '#cbeac4', '#a7dcb5', '#7acbc4', '#4db2d2', '#2a8bbd', '#0866aa']
            elif colormapName == 'PuBu':
                colormapList = ['#fff7fb', '#ebe6f1', '#cfd0e5', '#a5bcda', '#73a8ce', '#358fbf', '#046faf', '#03598b']
            elif colormapName == 'YlGnBu':
                colormapList = ['#ffffd9', '#ecf7b1', '#c6e8b4', '#7eccbb', '#40b5c3', '#1d90bf', '#225da7', '#243392']
            elif colormapName == 'PuBuGn':
                colormapList = ['#fff7fb', '#ebe1ef', '#cfd0e5', '#a5bcda', '#66a8ce', '#348fbe', '#018088', '#016a58']
            elif colormapName == 'BuGn':
                colormapList = ['#f7fcfd', '#e4f4f8', '#cbebe5', '#98d7c8', '#65c1a3', '#40ad75', '#228a44', '#006b2b']
            elif colormapName == 'YlGn':
                colormapList = ['#ffffe5', '#f6fbb8', '#d8efa2', '#acdc8d', '#77c578', '#40aa5c', '#228342', '#006736']
            elif colormapName == 'binary':
                colormapList = ['#ffffff', '#dfdfdf', '#bfbfbf', '#9f9f9f', '#7f7f7f', '#5f5f5f', '#3f3f3f', '#1f1f1f']
            elif colormapName == 'gist_yarg':
                colormapList = ['#ffffff', '#dfdfdf', '#bfbfbf', '#9f9f9f', '#7f7f7f', '#5f5f5f', '#3f3f3f', '#1f1f1f']
            elif colormapName == 'gist_gray':
                colormapList = ['#000000', '#202020', '#404040', '#606060', '#808080', '#a0a0a0', '#c0c0c0', '#e0e0e0']
            elif colormapName == 'gray':
                colormapList = ['#000000', '#202020', '#404040', '#606060', '#808080', '#a0a0a0', '#c0c0c0', '#e0e0e0']
            elif colormapName == 'bone':
                colormapList = ['#000000', '#1c1b26', '#38374d', '#545473', '#707b8f', '#8ca1ab', '#a8c7c7', '#d4e3e3']
            elif colormapName == 'pink':
                colormapList = ['#1e0000', '#744949', '#a16868', '#c2827f', '#d0ab93', '#ddcda4', '#e9e9b6', '#f4f4de']
            elif colormapName == 'spring':
                colormapList = ['#ff00ff', '#ff20df', '#ff40bf', '#ff609f', '#ff807f', '#ffa05f', '#ffc03f', '#ffe01f']
            elif colormapName == 'summer':
                colormapList = ['#007f66', '#208f66', '#409f66', '#60af66', '#80bf66', '#a0cf66', '#c0df66', '#e0ef66']
            elif colormapName == 'autumn':
                colormapList = ['#ff0000', '#ff2000', '#ff4000', '#ff6000', '#ff8000', '#ffa000', '#ffc000', '#ffe000']
            elif colormapName == 'winter':
                colormapList = ['#0000ff', '#0020ef', '#0040df', '#0060cf', '#0080bf', '#00a0af', '#00c09f', '#00e08f']
            elif colormapName == 'cool':
                colormapList = ['#00ffff', '#20dfff', '#40bfff', '#609fff', '#807fff', '#a05fff', '#c03fff', '#e01fff']
            elif colormapName == 'hot':
                colormapList = ['#0a0000', '#5e0000', '#b20000', '#ff0700', '#ff5b00', '#ffaf00', '#ffff06', '#ffff84']
            elif colormapName == 'afmhot':
                colormapList = ['#000000', '#400000', '#800000', '#c04000', '#ff8000', '#ffc041', '#ffff81', '#ffffc1']
            elif colormapName == 'gist_heat':
                colormapList = ['#000000', '#300000', '#600000', '#900000', '#c00000', '#f04100', '#ff8102', '#ffc183']
            elif colormapName == 'copper':
                colormapList = ['#000000', '#27180f', '#4f311f', '#764a2f', '#9e633f', '#c57c4f', '#ed955f', '#ffae6f']
            elif colormapName == 'PRGn':
                colormapList = ['#40004b', '#7e3b8d', '#ad8bbd', '#dec8e2', '#f6f6f6', '#cbeac5', '#7dc37e', '#288340']
            elif colormapName == 'BrBG':
                colormapList = ['#543005', '#995d12', '#cfa255', '#f0dfb2', '#f4f4f4', '#b3e2db', '#58b0a6', '#0c7068']
            elif colormapName == 'PuOr':
                colormapList = ['#7f3b08', '#be6209', '#ee9d3c', '#fdd6a2', '#f6f6f6', '#cdcde4', '#978dbd', '#5d378f']
            elif colormapName == 'RdGy':
                colormapList = ['#67001f', '#bb2a33', '#e58368', '#faceb6', '#fefefe', '#d5d5d5', '#9f9f9f', '#595959']
            elif colormapName == 'RdBu':
                colormapList = ['#67001f', '#bb2a33', '#e58368', '#faceb6', '#f6f6f6', '#bfdceb', '#68aacf', '#286fb0']
            elif colormpaName == 'RdYlGn':
                colormapList = ['#a50026', '#de3f2e', '#f88e52', '#fdd481', '#fefebd', '#cbe881', '#84ca66', '#2a9f54']
            elif colormapName == 'Spectral':
                colormapList = ['#9e0142', '#dc494b', '#f88e52', '#fdd481', '#fefebe', '#d5ee9b', '#86cea4', '#3d94b7']
            elif colormapName == 'coolwarm':
                colormapList = ['#3a4cc0', '#6182ea', '#8daffd', '#b8cff8', '#dddcdb', '#f4c3ab', '#f39879', '#dc5e4b']
            elif colormapName == 'bwr':
                colormapList = ['#0000ff', '#4040ff', '#8080ff', '#c0c0ff', '#fffefe', '#ffbebe', '#ff7e7e', '#ff3e3e']
            elif colormapName == 'seismic':
                colormapList = ['#00004c', '#0000a6', '#0101ff', '#8181ff', '#fffdfd', '#ff7d7d', '#fd0000', '#bd0000']
            elif colormapName == 'hsv':
                colormapList = ['#ff0000', '#ffbd00', '#83ff00', '#00ff39', '#00fff5', '#004bff', '#7100ff', '#ff00cf']
            elif colormapName == 'Pastel1':
                colormapList = ['#fbb4ae', '#b3cde2', '#cceac5', '#decbe3', '#fed9a6', '#fefecb', '#e5d8be', '#fcdaec']
            elif colormapName == 'Pastel2':
                colormapList = ['#b3e2cd', '#f4cfb0', '#d7d3d9', '#e5cee5', '#ece0d6', '#eff3be', '#fbedb6', '#ebdecc']
            elif colormapName == 'Paired':
                colormapList = ['#a6cee3', '#569fa4', '#51af42', '#f78787', '#f07047', '#fe850a', '#ae90c5', '#ccbd99']
            elif colormapName == 'Accent':
                colormapList = ['#7fc97f', '#b6b1c9', '#edbb98', '#fee892', '#98b3a4', '#80429c', '#e21a62', '#b15c22']
            elif colormapName == 'Dark2':
                colormapList = ['#1b9e77', '#c16610', '#8d6b87', '#bd4298', '#a46952', '#98a713', '#d49c09', '#9c7327']
            elif colormapName == 'flag':
                colormapList = ['#ff0000', '#f10000', '#d90000', '#bf0000', '#a40000', '#880000', '#6c0000', '#500000']
            elif colormapName == 'prism':
                colormapList = ['#ff0000', '#f5ff00', '#0056c3', '#f00070', '#ffd700', '#00a362', '#a200ce', '#ff9400']
            elif colormapName == 'ocean':
                colormapList = ['#007f00', '#004f20', '#001f40', '#001060', '#004080', '#0070a0', '#41a0c0', '#a2d0e0']
            elif colormapName == 'gist_earth':
                colormapList = ['#000000', '#153878', '#2a737e', '#3b8d62', '#5da04b', '#99ae58', '#bcaa62', '#dab69f']
            elif colormapName == 'terrain':
                colormapList = ['#333399', '#0888ee', '#01cc66', '#81e57f', '#fefd98', '#beab75', '#815d56', '#c1afab']
            elif colormapName == 'gist_stern':
                colormapList = ['#000000', '#a52040', '#404080', '#6060c0', '#8080fc', '#a0a074', '#c0c011', '#e0e08a']
            elif colormapName == 'gnuplot':
                colormapList = ['#000000', '#5a00b4', '#7f04fe', '#9c0db2', '#b42000', '#c93e00', '#dd6c00', '#eeac00']
            elif colormapName == 'gnuplot2':
                colormapList = ['#000000', '#000080', '#0000ff', '#6400ff', '#c829d5', '#ff6995', '#ffa955', '#ffe915']
            elif colormapName == 'CMRmap':
                colormapList = ['#000000', '#26267f', '#4d26be', '#9a337e', '#fe4025', '#e58000', '#e5c01b', '#e6e683']
            elif colormapName == 'cubehelix':
                colormapList = ['#000000', '#1a1d3b', '#15534b', '#437730', '#a1794a', '#d383a9', '#c6b4ed', '#cbe7ef']
            elif colormapName == 'brg':
                colormapList = ['#0000ff', '#4000bf', '#80007f', '#c0003f', '#fe0100', '#be4100', '#7e8100', '#3ec100']
            elif colormapName == 'gist_rainbow':
                colormapList = ['#ff0028', '#ff8300', '#cdff00', '#20ff00', '#00ff8b', '#00c5ff', '#0017ff', '#9600ff']
            elif colormapName == 'jet':
                colormapList = ['#00007f', '#0000ff', '#0080ff', '#15ffe1', '#7cff79', '#e4ff12', '#ff9400', '#ff1d00']
            elif colormapName == 'nipy_spectral':
                colormapList = ['#000000', '#4200a1', '#0077dd', '#00aa97', '#00bc00', '#66ff00', '#ffc900', '#eb0000']
            else: colormapList = ['#e41a1c','#377eb7','#4dae4a','#994ea1','#ff8100','#fdfb32','#a7572b','#f481bd']
        return colormapList

    def turnOnFlag2ChangeColor(self):
        self.flag2ChangeColor = True

    def turnOffFlag2ChangeColor(self):
        self.flag2ChangeColor = False

    def event_editNewColorEntry(self,var,*args):
        import re
        if not self.flag2ChangeColor: return
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
    def __init__(self,parent,B,color,colorHistoryList=COLOR_LIST_HISTORY,colormapName=None,nbColorMap=None,extra_data=None):
        if nbColorMap is None: nbColorMap = tkPlotXY.NUM_COLORS
        if colormapName is None: colormapName = tkPlotXY.COLOR_MAP
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
        colorControler = ColorControler(self,color=color,colorHistoryList=history[:tkPlotXY.NUM_COLORS],colormapName=colormapName,nbColorMap=nbColorMap)

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
