# - The Cassiopee installer -
# Check/Set necessary paths for compilation

WIDTH=40

#==============================================================================
# Doc:
# Cette app permet de mettre facilement les chemins necessaires a
# la compilation de Cassiopee.
# Machine doit correspondre au uname de la machine ou a un nom generique
# matchant le uname ou a un nom ELSAPROD definit dans l'environement.
# Si les chemins sont remplis, faire entree dans le widget de texte permet
# de verifier l'existence de la librairie et des includes.
# Si le texte est vide, l'app essai de deviner le chemin a partir du PATH
# et LD_LIBRARY_PATH de l'environnement.
#
# Si le python installe n'a pas tkinter, il reste possible d'ajouter
# tous les chemins dans config.py (ancien systeme).
#==============================================================================
try: import Tkinter as TK
except: raise ImportError('Cassiopee installer requires tkinter.\nIf your python dont run tkinter, edit the installBase.py file manually and add your machine settings.')

import Dist
import os.path

#==============================================================================
def quit(event=None):
    import os; os._exit(0)

#==============================================================================
# Met les valeurs par defaut dans l'interface
# On met gcc car c'est le plus rependu
#==============================================================================
def setDefaultVars():
    Vdescription.set("")
    Vf77compiler.set("gfortran")
    Vf90compiler.set("gfortran")
    VCppcompiler.set("gcc")

    VuseOMP.set(1)
    VuseStatic.set(0)
    VCPlotOffScreen.set(0)
    
    VadditionalIncludePaths.set("[]")
    VadditionalLibs.set("[]")
    VadditionalLibPaths.set("[]")

#==============================================================================
# Remplace les variables de l'interfaces par celles fournies dans la liste v
#==============================================================================
def setVars(v):
    Vdescription.set(v[0])
    Vf77compiler.set(v[1])
    Vf90compiler.set(v[2])
    VCppcompiler.set(v[3])

    VuseOMP.set(v[6])
    VuseStatic.set(v[7])
    VCPlotOffScreen.set(v[8])
    
    VadditionalIncludePaths.set(v[9])
    VadditionalLibs.set(v[10])
    VadditionalLibPaths.set(v[11])
    
#==============================================================================
def saveConfigFile(event=None):
    machine = Vmachine.get()
    
    try: additionalIncludePaths = eval(VadditionalIncludePaths.get())
    except: additionalIncludePaths = []

    try: additionalLibs = eval(VadditionalLibs.get())
    except: additionalLibs = []

    try: additionalLibPaths = eval(VadditionalLibPaths.get())
    except: additionalLibPaths = []

    dict[machine] = [Vdescription.get(),
                     Vf77compiler.get(), 
                     Vf90compiler.get(),
                     VCppcompiler.get(),
                     VuseOMP.get(), VuseStatic.get(), VCPlotOffScreen.get(),
                     additionalIncludePaths, additionalLibs,
                     additionalLibPaths]
    Dist.writeInstallBase(dict)
    return

#==============================================================================
def readConfigFile(event=None):
    machine = Vmachine.get()
    return

#==============================================================================
def setMachineName(name):
    Vmachine.set(name)
    changeMachineName()

#==============================================================================
def changeMachineName(event=None):
    machine = Vmachine.get()
    
    # Reset color
    if WIDGETS != {}:
        for name in ['Cppcompiler', 'f77compiler',
                     'f90compiler', 
                     'additionalIncludePaths', 'additionalLibs',
                     'additionalLibPaths']:
            entry = WIDGETS[name]
            entry.config(background=entry.master.cget('bg'))
            entry.config(foreground=entry.master.cget('fg'))
    
    key = ''
    for i in dict:
        if re.compile(i).search(machine) is not None:
            key = i; break

    if key != '': # already defined in install base
        v = dict[key]
        setVars(v)
    else: setDefaultVars()

#==============================================================================
# Check lib
# IN: name: name of the entry
#==============================================================================
def check__(name):
    ok = False
    if name == 'Python':
        try:
            (pythonVersion, pythonIncDir, pythonLibDir, pythonLibs) = Dist.checkPython()
            out = ['Python: %s'%pythonVersion]
        except: out = ['Python: include are missing.']
    elif name == 'Numpy':
        try:
            (numpyVersion, numpyIncDir, ti) = Dist.checkNumpy()
            out = ['Numpy: %s'%numpyVersion]
        except: out = ['Numpy: is missing or numpy includes are missing.']
    elif name == 'Cpp': # C++
        additionalLibPaths = eval(VadditionalLibPaths.get())
        (ok, CppLibs, CppLibPaths) = Dist.checkCppLibs(
            [], additionalLibPaths, VCppcompiler.get(), VuseOMP.get())
        if ok: out = ['C++: OK']
        else: out = ['C++: Fail']
    elif name == 'f77': # Fortran
        additionalLibPaths = eval(VadditionalLibPaths.get())
        (ok, FLibs, FLibPaths) = Dist.checkFortranLibs(
            [], additionalLibPaths, Vf77compiler.get(), VuseOMP.get())
        if ok: out = ['f77: OK']
        else: out = ['f77: Fail']
    elif name == 'f90': # Fortran
        additionalLibPaths = eval(VadditionalLibPaths.get())
        (ok, FLibs, FLibPaths) = Dist.checkFortranLibs(
            [], additionalLibPaths, Vf90compiler.get(), VuseOMP.get())
        if ok: out = ['f90: OK']
        else: out = ['f90: Fail']   
    
    elif name == 'png':
        additionalLibPaths = eval(VadditionalLibPaths.get())
        additionalIncludePaths = eval(VadditionalIncludePaths.get())
        (ok, pngIncDir, pngLib) = Dist.checkPng(additionalLibPaths,
                                                additionalIncludePaths)
        if ok: out = ['png: OK']
        else: out = ['png: libpng or png.h is missing']

    elif name == 'mpeg':
        additionalLibPaths = eval(VadditionalLibPaths.get())
        additionalIncludePaths = eval(VadditionalIncludePaths.get())
        (ok, mpgIncDir, mpgLib) = Dist.checkMpeg(additionalLibPaths,
                                                 additionalIncludePaths)
        if ok: out = ['mpeg: OK']
        else: out = ['mpeg: missing']
        
    elif name == 'hdf':
        additionalLibPaths = eval(VadditionalLibPaths.get())
        additionalIncludePaths = eval(VadditionalIncludePaths.get())
        (ok, hdfIncDir, hdfLib) = Dist.checkHdf(additionalLibPaths,
                                                additionalIncludePaths)
        if ok: out = ['hdf: OK']
        else: out = ['hdf: missing']
    return out

#==============================================================================
# Verifie tout!
#==============================================================================
def check(event=None):

    out = []
    # Check python
    out += check__('Python')
    out += check__('Numpy')

    out += check__('Cpp')
    out += check__('f77')
    out += check__('f90')
    out += check__('png')
    out += check__('mpeg')

    out += check__('hdf')
    return out

#==============================================================================
# guess/check Adf
#==============================================================================
def checkAdf(event=None):
    path = VadfPath.get()
    entry = WIDGETS['adfPath']
    (ok, adfIncDir, adfLib) = Dist.checkAdf(path)
    if (path == ''): # guess
        if ok:
            adfLib = os.path.split(adfLib)
            VadfPath.set('[\''+adfLib[0]+'\']')
    if ok:
        entry.config(background='green')
        entry.config(foreground='black')
    else:
        entry.config(background='red')
        entry.config(foreground='white')
    return ok

#==============================================================================
# guess/check Hdf
#==============================================================================
def checkHdf(event=None):
    path = VhdfPath.get()
    entry = WIDGETS['hdfPath']
    (ok, hdfIncDir, hdfLib) = Dist.checkHdf(path)
    if path == '': # guess
        if ok:
            hdfLib = os.path.split(hdfLib)
            VhdfPath.set('[\''+hdfLib[0]+'\']')
    if ok:
        entry.config(background='green')
        entry.config(foreground='black')
    else:
        entry.config(background='red')
        entry.config(foreground='white')
    return ok

#==============================================================================
# guess/check png
#==============================================================================
def checkPng(event=None):
    path = VpngPath.get()
    entry = WIDGETS['pngPath']
    (ok, pngIncDir, pngLib) = Dist.checkPng(path)
    if path == '': # guess
        if ok:
            pngLib = os.path.split(pngLib)
            VpngPath.set('[\''+pngLib[0]+'\']')
    if ok:
        entry.config(background='green')
        entry.config(foreground='black')
    else:
        entry.config(background='red')
        entry.config(foreground='white')
    return ok

#==============================================================================
# guess/check ffmpeg
#==============================================================================
def checkMpeg(event=None):
    path = VmpegPath.get()
    entry = WIDGETS['mpegPath']
    (ok, mpegIncDir, mpegLib) = Dist.checkMpeg(path)
    if path == '': # guess
        if ok:
            mpegLib = os.path.split(mpegLib)
            VmpegPath.set('[\''+mpegLib[0]+'\']')
    if ok:
        entry.config(background='green')
        entry.config(foreground='black')
    else:
        entry.config(background='red')
        entry.config(foreground='white')
    return ok

#==============================================================================
# Load install base
import installBase
dict = installBase.installDict

#==============================================================================
# Get machine
import platform, re, os
a = platform.uname()
system = a[0] # Linux, Windows
host = a[1]   # host name
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

#==============================================================================
master = TK.Tk()

# - VARS -
Vdescription = TK.StringVar(master) ; Vdescription.set('')
Vmachine = TK.StringVar(master) ; Vmachine.set(host)
VCppcompiler = TK.StringVar(master) ; VCppcompiler.set('')
Vf77compiler = TK.StringVar(master) ; Vf77compiler.set('')
Vf90compiler = TK.StringVar(master) ; Vf90compiler.set('')
VuseOMP = TK.IntVar(master) ; VuseOMP.set(0)
VuseStatic = TK.IntVar(master) ; VuseStatic.set(0)
VCPlotOffScreen = TK.IntVar(master) ; VCPlotOffScreen.set(0)
VadditionalIncludePaths = TK.StringVar(master) ; VadditionalIncludePaths.set('')
VadditionalLibs = TK.StringVar(master) ; VadditionalLibs.set('')
VadditionalLibPaths = TK.StringVar(master) ; VadditionalLibPaths.set('')

# Widgets
WIDGETS = {}

# Init
key = ''
for i in dict.keys():
    if (re.compile(i).search(host) is not None or
        re.compile(i).search(prod) is not None):
        key = i ; break

if (re.compile(key).search(host) is not None): setMachineName(host)
if (re.compile(key).search(prod) is not None): setMachineName(prod)

if key != '': # already defined in install base
    v = dict[key]
    setVars(v)

#==============================================================================
def instructions():
    winl = TK.Toplevel(border=0)
    winl.title('Instructions')
    winl.columnconfigure(0, weight=1)
    winl.rowconfigure(0, weight=1)
    # position de la fenetre parent
    xpos = winl.master.winfo_rootx()
    ypos = winl.master.winfo_rooty()
    winl.geometry("%+d%+d" % (xpos+40, ypos)) 
    scrollbar = TK.Scrollbar(winl, orient=TK.VERTICAL, width=10)
    scrollbar.grid(sticky=TK.NSEW, row=0, column=1)
    
    textWidget = TK.Text(winl, yscrollcommand=scrollbar.set,
                         width=45, height=20, background='White')
    textWidget.tag_config('title', justify=TK.CENTER, foreground='blue')
    textWidget.grid(sticky=TK.NSEW, row=0, column=0)
    scrollbar.config(command=textWidget.yview)
    myText = "- Cassiopee - \n"
    textWidget.insert(TK.END, myText, 'title')
    myText = "A CFD pre- and post-processing tool"
    textWidget.insert(TK.END, myText, 'title')
    myText = "\n\nInstallation from sources under unix or mingw.\n"
    textWidget.insert(TK.END, myText)
    myText = "Check requirements:\n\n"
    myText += "Press File/Check.\n\n"
    myText += "Press: File/save installBase\n"
    myText += "Then exit installer and type:\n\n"
    myText += "./install --prefix=yourInstallationPath\n\n"
    myText += "If something fails, then install the dev package of the missing libraries or change the path if libraries are installed to non standard directories.\n"
    textWidget.insert(TK.END, myText)
    return

#==============================================================================
def printCheck():
    out = check()
    winl = TK.Toplevel(border=0)
    winl.title('Check')
    winl.columnconfigure(0, weight=1)
    winl.rowconfigure(0, weight=1)
    # position de la fenetre parent
    xpos = winl.master.winfo_rootx()
    ypos = winl.master.winfo_rooty()
    winl.geometry("%+d%+d" % (xpos+40, ypos)) 
    scrollbar = TK.Scrollbar(winl, orient=TK.VERTICAL, width=10)
    scrollbar.grid(sticky=TK.NSEW, row=0, column=1)
    
    textWidget = TK.Text(winl, yscrollcommand=scrollbar.set,
                         width=45, height=20, background='White')
    textWidget.tag_config('title', justify=TK.CENTER, foreground='blue')
    textWidget.grid(sticky=TK.NSEW, row=0, column=0)
    scrollbar.config(command=textWidget.yview)
    myText = ''
    for i in out: myText += i + '\n'
    textWidget.insert(TK.END, myText, 'title')

#==============================================================================
# Main window
#==============================================================================
master.title('*Cassiopee* installer')
master.columnconfigure(0, weight=1)

menu = TK.Menu(master)
file = TK.Menu(menu, tearoff=0)
menu.add_cascade(label='File', menu=file)
file.add_command(label='Check', command=printCheck)
file.add_command(label='Save installBase', command=saveConfigFile)
file.add_command(label='Quit', command=quit)
machines = TK.Menu(menu, tearoff=0)
menu.add_cascade(label='Machines', menu=machines)
for i in dict:
    machines.add_command(label=i, command=lambda i=i : setMachineName(i))
master.config(menu=menu)
help = TK.Menu(menu, tearoff=0)
menu.add_cascade(label='Help', menu=help)
help.add_command(label='Instructions', command=instructions)

B = TK.Label(master, text="Machine")
B.grid(row=0, column=0)
B = TK.Entry(master, textvariable=Vmachine, background='White', width=WIDTH)
B.bind('<Return>', changeMachineName)
B.grid(row=0, column=1)

#==============================================================================
# Required frame
Frame = TK.LabelFrame(master, borderwidth=2, relief=TK.RIDGE,
                      text='Required')
Frame.columnconfigure(0, weight=1)
Frame.columnconfigure(1, weight=1)
Frame.grid(row=1, column=0, columnspan=2, sticky=TK.EW)

# C compiler
B = TK.Label(Frame, text="C++ compiler")
B.grid(row=0, column=0)
B = TK.Entry(Frame, textvariable=VCppcompiler, background='White', width=WIDTH)
WIDGETS['Cppcompiler'] = B
B.grid(row=0, column=1, sticky=TK.EW)

# f77 compiler
B = TK.Label(Frame, text="f77 compiler")
B.grid(row=2, column=0)
B = TK.Entry(Frame, textvariable=Vf77compiler, background='White', width=WIDTH)
WIDGETS['f77compiler'] = B
B.grid(row=2, column=1, sticky=TK.EW)

# f90 compiler
B = TK.Label(Frame, text="f90 compiler")
B.grid(row=5, column=0)
B = TK.Entry(Frame, textvariable=Vf90compiler, background='White', width=WIDTH)
WIDGETS['f90compiler'] = B
B.grid(row=5, column=1, sticky=TK.EW)

# open MP
F = TK.Frame(Frame)
B = TK.Label(F, text="useOMP")
B.grid(row=0, column=0)
B = TK.Checkbutton(F, text="", variable=VuseOMP)
B.grid(row=0, column=1)

# Static production
B = TK.Label(F, text="staticLibs")
B.grid(row=0, column=2)
B = TK.Checkbutton(F, text="", variable=VuseStatic)
B.grid(row=0, column=3)
F.grid(row=8, column=0, columnspan=2, sticky=TK.EW)

# additionalIncludePaths
B = TK.Label(Frame, text="+includePaths")
B.grid(row=9, column=0)
B = TK.Entry(Frame, textvariable=VadditionalIncludePaths, background='White',
             width=WIDTH)
WIDGETS['additionalIncludePaths'] = B
B.grid(row=9, column=1, sticky=TK.EW)

# additionalLibs
B = TK.Label(Frame, text="+libs")
B.grid(row=10, column=0)
B = TK.Entry(Frame, textvariable=VadditionalLibs, background='White',
             width=WIDTH)
WIDGETS['additionalLibs'] = B
B.grid(row=10, column=1, sticky=TK.EW)

# additionalLibPaths
B = TK.Label(Frame, text="+libPaths")
B.grid(row=11, column=0)
B = TK.Entry(Frame, textvariable=VadditionalLibPaths, background='White',
             width=WIDTH)
WIDGETS['additionalLibPaths'] = B
B.grid(row=11, column=1, sticky=TK.EW)

#==============================================================================
#B = TK.Button(master, text="Check", command=check)
#B.grid(row=5, column=0)
#B = TK.Button(master, text="Save config file", command=saveConfigFile)
#B.grid(row=5, column=1)

master.mainloop()
