try: import tkinter as TK
except: import Tkinter as TK
try: import Tk as CTK
except: from . import Tk as CTK

ttk = None
try: import tkinter.ttk as ttk
except:
    try: import ttk
    except: ttk = None
#Uncomment that for a pure Tk interface
#ttk = None

# Couleur du fond (theme)
BACKGROUNDCOLOR = '#FFFFFF'
FOREGROUNDCOLOR = '#000000'

# Icon light (0) or dark (1) - deduced from theme colors
ICONMODE = 0

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

# Themes a charger dynamiquement
dynamicThemes = {
    "aquativo": "aquativo.tcl",
    "arc": "arc.tcl",
    "black": "black.tcl",
    #"blue": "blue.tcl",
    "breeze": "breeze.tcl",
    "clearlooks": "clearlooks.tcl",
    #"elegance": "elegance.tcl",
    "equilux": "equilux.tcl",
    "itft1": "itft1.tcl",
    #"keramik": "keramik.tcl",
    #"plastik": "plastik.tcl",
    "radiance": "radiance.tcl",
    #"scid": "scid.tcl",
    "smog": "smog.tcl",
    "ubuntu": "ubuntu.tcl",
    "winxpblue": "winxpblue.tcl",
    #"forest": "forest-dark.tcl"
}

#=================================================================
# Get available ttk themes
#=================================================================
def getAvailableThemes():
    if ttk is not None:
        l = ttk.Style().theme_names()
        l = list(l)
        # Ajout des themes dynamiques
        l += dynamicThemes.keys()
        l = set(l)
        l = list(l)
        l.sort()
        return l
    else: return ['default']

#================================================================
# Set theme if possible, otherwise try default themes
#================================================================
def setTheme(myTheme):
    if ttk is not None:

        if myTheme == "None": # no pref set
            myTheme = "clearlooks"

        # Load des themes dynamiques
        if myTheme in dynamicThemes.keys():
            import KCore.installPath
            folder = KCore.installPath.installPath+'/CPlot'
            try: # if already loaded
                CTK.WIDGETS['masterWin'].call("lappend", "auto_path", "[%s]"%folder + "/themes/%s"%myTheme)
                CTK.WIDGETS['masterWin'].eval("source %s/themes/%s/%s"%(folder,myTheme,dynamicThemes[myTheme]))
            except: pass
            ttk.Style().theme_use(myTheme); createStyles(); return

        # Load des themes availables
        available = ttk.Style().theme_names()
        if myTheme in available:
            ttk.Style().theme_use(myTheme); createStyles(); return

#=============================================================
# Create specific styles
#=============================================================
def createStyles():
    if ttk is not None:

        global BACKGROUNDCOLOR, FOREGROUNDCOLOR, ICONMODE
        style = ttk.Style()
        # Get theme colors
        ret = style.lookup('TFrame', 'background')
        if ret != "": BACKGROUNDCOLOR = str(ret)
        ret = style.lookup('TFrame', 'foreground')
        if ret != "": FOREGROUNDCOLOR = str(ret)
        #ret = style.lookup('TFrame', 'font')
        #if ret != "": CTK.GENERALFONT = ret

        # Deduce icon style (light or dark)
        if BACKGROUNDCOLOR == 'black': BACKGROUNDCOLOR = "#000000"
        elif BACKGROUNDCOLOR == 'white': BACKGROUNDCOLOR = "#FFFFFF"
        if FOREGROUNDCOLOR == 'black': FOREGROUNDCOLOR = "#000000"
        elif FOREGROUNDCOLOR == 'white': FOREGROUNDCOLOR = "#FFFFFF"
        if BACKGROUNDCOLOR[0] == "#" and FOREGROUNDCOLOR[0] == "#":
            R = int(BACKGROUNDCOLOR[1:3], 16)
            G = int(BACKGROUNDCOLOR[3:5], 16)
            B = int(BACKGROUNDCOLOR[5:7], 16)
            intensityBackground = 0.299*R + 0.587*G + 0.114*B
            R = int(FOREGROUNDCOLOR[1:3], 16)
            G = int(FOREGROUNDCOLOR[3:5], 16)
            B = int(FOREGROUNDCOLOR[5:7], 16)
            intensityForeground = 0.299*R + 0.587*G + 0.114*B
            if intensityBackground > intensityForeground: ICONMODE = 0
            else: ICONMODE = 1
        else: ICONMODE = 0

        # Set all fonts
        style.configure('.', font=CTK.GENERALFONT)
        # Set all frame backgrounds
        #style.configure('TFrame', bgColor=BACKGROUNDCOLOR)

        # K1 LabelFrameStyle
        style.configure("K1.TLabelframe.Label", font=CTK.FRAMEFONT)
        # K1 Scale Style
        style.configure("Horizontal.K1.TScale", showvalue=0, borderwidth=1)
        # Sunken Button Style
        style.configure("SUNKEN.TButton", relief=TK.SUNKEN)
        # Raised Button Style
        style.configure("RAISED.TButton", relief=TK.RAISED)
        # Iconic button
        style.configure("ICONIC.TButton", borderwidth=0, relief=TK.FLAT)
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
#===========================================================
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
        if 'borderwidth' in kwargs:
            kwargs.pop('borderwidth', None)
            style = 1
        b = ttk.Button(*args, **kwargs)
        if style == 1: b.configure(style="ICONIC.TButton")
        return b

# TK button avec background du theme
def Button2(*args, **kwargs):
    B = TK.Button(*args, **kwargs)
    B.config(bg=BACKGROUNDCOLOR)
    return B

def configButton(B, color):
    if ttk is None:
        B.config(bg=color)
        B.config(activebackground=color)

def OptionMenu(*args, **kwargs):
    if ttk is None: return TK.OptionMenu(*args, **kwargs)
    else:
        # add default arg
        largs = (args[0],args[1],None)+args[2:]
        o = ttk.OptionMenu(*largs, **kwargs)
        o["menu"].config(bg=BACKGROUNDCOLOR, fg=FOREGROUNDCOLOR)
        return o

# Si ttk, renvoie une combobox avec accelerateur clavier
def Combobox(*args, **kwargs):
    if ttk is None: # process args here
        raise ValueError('No comboxbox.')
    else:
        return ComboboxAuto(*args, **kwargs)

def superOptionMenu(F, var, itemList, command,
                    updateCommand1, updateCommand2):
    if ttk is None:
        B = TK.OptionMenu(F, var, *itemList,
                          command=command)
        F.bind('<Enter>', updateCommand1)
    else:
        B = ttk.Combobox(F, textvariable=var,
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

# Pas de menu specifique en ttk
def Menu(*args, **kwargs):
    if ttk is None: return TK.Menu(*args, **kwargs)
    else:
        M = TK.Menu(*args, **kwargs)
        M.config(bg=BACKGROUNDCOLOR, fg=FOREGROUNDCOLOR)
        return M

def Menubutton(*args, **kwargs):
    if ttk is None: return TK.Menubutton(*args, **kwargs)
    else: return ttk.Menubutton(*args, **kwargs)

# Pas de texte specifique en ttk
def Text(*args, **kwargs):
    if ttk is None: return TK.Text(*args, **kwargs)
    else: return TK.Text(*args, **kwargs)

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

def Toplevel(*args, **kwargs):
    if ttk is None: return TK.Toplevel(*args, **kwargs)
    else:
        t = TK.Toplevel(*args, **kwargs)
        t.config(bg=BACKGROUNDCOLOR)
        return t

# Si ttk, renvoie une listbox avec accelerateur clavier
def Listbox(*args, **kwargs):
    if ttk is None: return TK.Listbox(*args, **kwargs)
    else:
        return ListboxAuto(*args, **kwargs)
        #return TK.Listbox(*args, **kwargs)

# si ttk, renvoie un NoteBook
def NoteBook(F):
    if ttk is None: return TK.Frame(F)
    else: return ttk.Notebook(F)

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

if ttk is not None:
    # Combobox avec accelerateur clavier
    class ComboboxAuto(ttk.Combobox):
        def __init__(self, *args, **kwargs):
            ttk.Combobox.__init__(self, *args, **kwargs)
            self.bind('<KeyRelease>', self.keyRelease)
            self.bind('<Enter>', self.onEnter)
            self.bind('<Leave>', self.onLeave)
            self._searchString = ''

        def onEnter(self, event):
            self.focus_set()
            self._searchString = ''

        def onLeave(self, event):
            self.focus_set()
            self._searchString = ''

        def keyRelease(self, event):
            keysym = event.keysym
            if keysym == 'underscore': keysym = "_"
            values = self['values']

            if len(keysym) == 1: # not code
                key = keysym[0]
                self._searchString += key
                ss = self._searchString
                ss = ss.lower()
                ls = len(ss)

                # Match first chars
                #index = 0
                #for i in values:
                #    if ss <= i[0:ls].lower(): break
                #    index += 1
                #if index >= len(values): index = 0

                # Match any place
                index = 0
                for i in values:
                    if ss in i.lower(): break
                    index += 1
                if index < len(values): self.set(values[index])
            elif keysym == "Right":
                index = self.current()
                index += 1
                if index >= len(values): index = 0
                self.set(values[index])
            elif keysym == "Left":
                index = self.current()
                index -= 1
                if index < 0: index = len(values)-1
                self.set(values[index])
            elif keysym == "Return":
                self._searchString = ''
            elif keysym == "Escape":
                self._searchString = ''
            elif keysym == "Delete":
                self._searchString = ''
            elif keysym == "BackSpace":
                ls = len(self._searchString)
                if ls > 0: self._searchString = self._searchString[0:ls-1]
                else: self._searchString = ''
                #print(self._searchString)

    # Listbox avec accelerateur clavier
    class ListboxAuto(TK.Listbox):
        def __init__(self, *args, **kwargs):
            TK.Listbox.__init__(self, *args, **kwargs)
            TK.Listbox.config(self, bg=BACKGROUNDCOLOR, fg=FOREGROUNDCOLOR)
            self.bind('<KeyRelease>', self.keyRelease)
            self.bind('<Enter>', self.onEnter)
            self.bind('<Leave>', self.onLeave)
            self.bind('<<ListboxSelect>>', self.onSelect)
            self._searchString = ''

        def onEnter(self, event):
            self.focus_set()
            self._searchString = ''

        def onLeave(self, event):
            self.focus_set()
            self._searchString = ''

        def onSelect(self, event):
            self.focus_set()
            self._searchString = ''

        def keyRelease(self, event):
            keysym = event.keysym
            if keysym == 'underscore': keysym = "_"
            values = self.get(0, TK.END)

            if len(keysym) == 1: # not code
                key = keysym[0]
                self._searchString += key
                ss = self._searchString
                ss = ss.lower()
                ls = len(ss)

                # Match first chars
                #index = 0
                #for i in values:
                #    if ss <= i[0:ls].lower(): break
                #    index += 1
                #if index >= len(values): index = 0

                # Match any place
                index = 0
                for i in values:
                    if ss in i.lower(): break
                    index += 1
                if index < len(values): self.see(index)
            elif keysym == "Return":
                self._searchString = ''
            elif keysym == "Delete":
                self._searchString = ''
            elif keysym == "BackSpace":
                ls = len(self._searchString)
                if ls > 0: self._searchString = self._searchString[0:ls-1]
                else: self._searchString = ''
                #print(self._searchString)
