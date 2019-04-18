try: import Tkinter as TK
except: import tkinter as TK
try: import Tk as CTK 
except: from . import Tk as CTK

ttk = None
try: import ttk
except: 
    try: import tkinter.ttk as ttk
    except: ttk = None
#Uncomment that for a pure Tk interface
#ttk = None

#=================================================================
# Installe des themes livres avec CPlot
#=================================================================
def installLocalThemes(win):
    try:
        import KCore.installPath
        folder = KCore.installPath.installPath+'/CPlot'
        win.call("lappend", "auto_path", "[%s]"%folder + "/themes")
        win.eval("source %s/themes/pkgIndex.tcl"%folder)
        #win.call("package", "require", "Img")
        win.call("package", "require", "Tk", "8.5")
    except: pass

#=================================================================
# Get available ttk themes
def getAvailableThemes():
    if ttk is not None:
        return ttk.Style().theme_names()
    else: return ['default']

#================================================================
# Set theme if possible, otherwise try default themes
def setTheme(myTheme):
    if ttk is not None:
        available = ttk.Style().theme_names()
        if myTheme in available:
            ttk.Style().theme_use(myTheme); createStyles(); return
        tryThemes = ["xpnative", "clearlooks", "plastik",
                     "black", "aquativa", "classic", "default"]
        for t in tryThemes:
            if t in available:
                ttk.Style().theme_use(t); createStyles(); return

#=============================================================
# Create specific styles
def createStyles():
    if ttk is not None:
        style = ttk.Style()
        # All fonts
        style.configure('.', font=CTK.GENERALFONT)
        # K1 LabelFrameStyle
        style.configure("K1.TLabelframe.Label", font=CTK.FRAMEFONT)
        # K1 Scale Style
        style.configure("Horizontal.K1.TScale", showvalue=0, borderwidth=1)
        # Sunken Button Style
        style.configure("SUNKEN.TButton", relief=TK.SUNKEN)
        # Raised Button Style
        style.configure("RAISED.TButton", relief=TK.RAISED)
        # Iconic button
        style.configure("ICONIC.TButton", borderwidth=0)
        # Radiobutton of left menu
        style.configure("MENU.TRadiobutton", relief=TK.GROOVE)
        # ScrollBar with width
        style.configure("Vertical.TKTREE.TScrollbar", width=10)
        style.configure("Horizontal.TKTREE.TScrollbar", width=10)
        # Red bg button#117864
        style.configure("RED.TButton", foreground='red')
        # Red bg button
        style.configure("GREEN.TButton", foreground='#117864')

#===========================================================
# Set global style and create specific styles
def setStyle():
    if ttk is not None:
        # pour windows : xpnative, vista, winnative
        myStyle = 'xpnative'
        # pour linux : clearlook, aquativa, arc
        #myStyle = 'default'
        available = ttk.Style().theme_names()
        # Custom theme
        if myStyle in available:
            ttk.Style().theme_use(myStyle)
        # - Custom styles -
        createStyles()

def LabelFrame(*args, **kwargs):
    if ttk is None: return TK.LabelFrame(*args, **kwargs)
    else:
        style = 0
        if 'font' in kwargs:
            style = 1; kwargs.pop('font', None)
        l = ttk.LabelFrame(*args, **kwargs)
        if style == 1: l.configure(style="K1.TLabelframe")
        return l

def Frame(*args, **kwargs):
    if ttk is None: return TK.Frame(*args, **kwargs)
    else: return ttk.Frame(*args, **kwargs)

def Button(*args, **kwargs):
    if ttk is None: return TK.Button(*args, **kwargs)
    else:
        style = 0
        if 'padx' in kwargs:
            kwargs.pop('padx', None)
            kwargs.pop('pady', None)
            style = 1
        b = ttk.Button(*args, **kwargs)
        if style == 1: b.configure(style="ICONIC.TButton")
        return b

def configButton(B, color):
    if ttk is None:
        B.config(bg=color)
        B.config(activebackground=color)

def OptionMenu(*args, **kwargs):
    if ttk is None: return TK.OptionMenu(*args, **kwargs)
    else:
        # add default arg
        largs = (args[0],args[1],None)+args[2:] 
        return ttk.OptionMenu(*largs, **kwargs)

def ComboBox(*args, **kwargs):
    if ttk is None: # process args here
        raise ValueError('No comboxbox.')
    else: return ttk.Combobox(*args, **kwargs)

def superOptionMenu(F, var, itemList, command, 
                    updateCommand1, updateCommand2):
    if ttk is None:
        B = TK.OptionMenu(F, var, *itemList,
                          command=com)
        F.bind('<Enter>', updateCommand1)
    else:
        B = ttk.Combobox(F, textvariable=VARS[0], 
                         values=itemList, 
                         state='readonly', width=10)
        B.bind('<<ComboboxSelected>>', command)
        F.bind('<Enter>', updateCommand2)
    return B

def Scale(*args, **kwargs):
    if ttk is None:
        val = -1
        if 'padx' in kwargs:
           val = kwargs['value']
           kwargs.pop('value', None)
        s = TK.Scale(*args, **kwargs)
        if val > 0: s.set(val)
        return s
    else:
        style = 0
        kwargs.pop('showvalue', None) # transmit in style?
        if 'borderwidth' in kwargs:
            style = 1
            kwargs.pop('borderwidth')
        s = ttk.Scale(*args, **kwargs)
        #if style == 1: s.configure(style="K1.TScale")
        return s

def Checkbutton(*args, **kwargs):
    if ttk is None: return TK.Checkbutton(*args, **kwargs)
    else: return ttk.Checkbutton(*args, **kwargs)

def Entry(*args, **kwargs):
    if ttk is None: return TK.Entry(*args, **kwargs)
    else: return ttk.Entry(*args, **kwargs)

def Menu(*args, **kwargs):
    return TK.Menu(*args, **kwargs)
    
def Text(*args, **kwargs):
    if ttk is None: return TK.Text(*args, **kwargs) 
    else: return ttk.Text(*args, **kwargs)

def Radiobutton(*args, **kwargs):
    if ttk is None: TK.Radiobutton(*args, **kwargs) 
    else:
        if 'offrelief' in kwargs: kwargs.pop('offrelief')
        if 'selectcolor' in kwargs: kwargs.pop('selectcolor')
        if 'indicatoron' in kwargs: kwargs.pop('indicatoron')
        b = ttk.Radiobutton(*args, **kwargs)
        b.configure(style='MENU.TRadiobutton')
        return b

def selectRadioButton(B):
    if ttk is None: B.select()
    else: B.invoke()

def deselectRadioButton(B):
    if ttk is None: B.deselect()
    else: B.state(["!selected"])

def Label(*args, **kwargs):
    if ttk is None: return TK.Label(*args, **kwargs)
    else: return ttk.Label(*args, **kwargs)

def Scrollbar(*args, **kwargs):
    if ttk is None: return TK.Scrollbar(*args, **kwargs)
    else:
        width = 0
        if 'width' in kwargs:
            width = kwargs.get('width', 10)         
            kwargs.pop('width')
        b = ttk.Scrollbar(*args, **kwargs)
        #if width != 0: b.configure(style='TKTREE.TScrollbar')
        return b

def Listbox(*args, **kwargs):
    if ttk is None: return TK.Listbox(*args, **kwargs)
    else: return TK.Listbox(*args, **kwargs)

def raiseButton(B):
    if ttk is None: B.config(relief=TK.RAISED)
    else: B.configure(style='RAISED.TButton')

def sunkButton(B):
    if ttk is None: B.config(relief=TK.SUNKEN)
    else: B.configure(style='SUNKEN.TButton')   

def setButtonRed(B):
    if ttk is None: B.config(bg='red')
    else: B.configure(style='RED.TButton')   

def setButtonGreen(B):
    if ttk is None: B.config(bg='green')
    else: B.configure(style='GREEN.TButton')   
