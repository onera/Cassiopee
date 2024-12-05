# *Cassiopee* Computation Center
import Tkinter as TK
import tkMessageBox
import os, sys, re, time
import subprocess
import CPlot.Tk as CTK
import KCore.Dist as Dist
from threading  import Thread
try: from Queue import Queue, Empty
except ImportError:
    from queue import Queue, Empty
ON_POSIX = 'posix' in sys.builtin_module_names

ttk = CTK.importTtk()

# CASSIOPEE var
CASSIOPEE = os.getenv('CASSIOPEE')
if CASSIOPEE == '':
    print('Error: CASSIOPEE must be present in your environement.')
    sys.exit()

# Systeme
mySystem = Dist.getSystem()[0]

# Success
regSuccess = re.compile('Writing restart.cgns')
regReal = re.compile('real')

# Separator
separator = ':'
separatorl = separator+' '

# Liste des cas existants dans le rep. Cases
CASES = []
# Info pour chaque cas (si il y en a)
INFOS = []
# Status de chaque cas
STATUS = []
# Widgets VARS
VARS = []
# Widgets
CSWINDOW = None
WIDGETS = {}

# RUNNING PROCESS
RUNPROC = 1

#==============================================================================
# - Case selector -

def _deleteCSWindow():
    try: CSWINDOW.destroy()
    except: pass

def _destroyCSWindow(event):
    global CSWINDOW
    CSWINDOW = None

def okCaseSelector(event=None):
    lb = WIDGETS['selectListBox']
    selection = lb.curselection()
    if (len(selection) == 0): return
    s = selection[0]
    s = int(s)
    t = CASES[s]
    VARS[0].set(t)
    _deleteCSWindow()
    _destroyCSWindow(event)

def cancelCaseSelector(event=None):
    _deleteCSWindow()
    _destroyCSWindow(event)

def caseSelector():
    global CSWINDOW
    if CSWINDOW is None:
        CSWINDOW = TK.Toplevel()
        CSWINDOW.columnconfigure(0, weight=1)
        CSWINDOW.rowconfigure(0, weight=1)
        CSWINDOW.title("Select a case...")
        xpos = CSWINDOW.master.winfo_rootx()
        ypos = CSWINDOW.master.winfo_rooty()
        CSWINDOW.geometry("%+d%+d" % (xpos, ypos))
        CSWINDOW.protocol("WM_DELETE_WINDOW", _deleteCSWindow)
        CSWINDOW.bind("<Destroy>", _destroyCSWindow)
        F = TK.Frame(CSWINDOW, borderwidth=0)
        F.grid(row=0, column=0)
        F.columnconfigure(0, weight=1)
        lb = TK.Listbox(F, selectmode=TK.SINGLE, 
                        width=120, height=22,
                        background='White')
        WIDGETS['selectListBox'] = lb
        lb.grid(row=0, column=0, columnspan=10, sticky=TK.EW)
        sb = TK.Scrollbar(F, orient=TK.VERTICAL)
        sb.grid(row=0, column=10, sticky=TK.NSEW)
        B = TK.Button(F, text='OK', command=okCaseSelector)
        B.grid(row=1, column=0, columnspan=8, sticky=TK.EW)
        B = TK.Button(F, text='Cancel', command=cancelCaseSelector)
        B.grid(row=1, column=9)
        # populate list
        getCases()
        for i in range(len(CASES)):
            txt = CASES[i].ljust(20) + separatorl + INFOS[i].ljust(75) + STATUS[i].ljust(6) 
            lb.insert(TK.END, txt)
    else:
        CSWINDOW.withdraw(); CSWINDOW.deiconify(); CSWINDOW.focus_set()

#==============================================================================
# Quit
#==============================================================================
def Quit(event=None):
    import os; os._exit(0)

#==============================================================================
# Ajoute la sortie de stdout dans une queue
#==============================================================================
def enqueueOutput(out, queue):
    for line in iter(out.readline, ''):
        queue.put(line)
    out.close()

#==============================================================================
# Update le no dans la liste box avec les donnees fournies
# no: numero du cas dans la liste box
# info = [case, nit, scale, gfx, status]
#==============================================================================
def updateListBox(no, info):
    s = buildString(info)
    listbox.delete(no, no)
    listbox.insert(no, s)
    listbox.update()

#==============================================================================
# Analyse une ligne de retour de run contenant le no d'it
# [Retourne -1: unknown, sinon retourne le %: it/nit]
# Retourne None: unknow, sinon retourne la string it/nit
#==============================================================================
def analyse(line, nit):
    s1 = line.split('/')
    if (len(s1) == 2): # it / nit
        it1 = s1[0].replace('-', '')
        it1 = it1.strip()
        s2 = s1[1].split('-')
        if (len(s2) == 1):
            it2 = s1[1].replace('-', '')
            it2 = it2.strip()
        else:
            it2 = s2[0].replace('-', '')
            it2 = it2.strip()
        i1 = int(it1)*1.; i2 = int(it2)*1.
        pourcentage = (nit-i2+i1)/nit*100.
        line = line.replace('-', '')
        line = line.strip()
        return (pourcentage, line)
    else:
        if regSuccess.search(line) is not None: return (-2., None)
        if regReal.search(line) is not None:
            line = line.replace('real', '')
            line = line.strip()
            return (-3., line)
    return (-1., None)

#==============================================================================
# Lance la commande cmd correspondant au cas no dans la liste box
#==============================================================================
def launchAndSurvey(cmd, no, info):
    global RUNPROC
    proc = subprocess.Popen(cmd,
                            stderr=subprocess.STDOUT,
                            stdout=subprocess.PIPE, shell=True,
                            bufsize=1, close_fds=ON_POSIX)
    RUNPROC = proc
    q = Queue()
    t = Thread(target=enqueueOutput, args=(proc.stdout, q))
    t.daemon = True # thread dies with the program
    t.start()
    nit = info[1]
    status = 0
    CPU = ''; iterf = ''
    while proc.poll() is None:
        time.sleep(0.1)
        #line = proc.stdout.readline() # bloquant
        try: line = q.get_nowait() # non bloquant
        except: line = ''

        if line != '':
            ret = analyse(line, int(nit))
            if (ret[0] > -0.5):
                info[4] = 'RUNNING (%2.2f'%ret[0]
                info[4] += ' % - '+ret[1]+')'
                iterf = ret[1]
                updateListBox(no, info)
            elif (ret[0] == -2.): status = 1
            elif (ret[0] == -3.): CPU = ret[1]
        master.update()

    # Vide le buffer
    sys.stdout.flush()
    line = 't'
    while (line != ''):
        try: line = q.get_nowait() # non bloquant
        except: line = ''
        if (line != ''):
            ret = analyse(line, int(nit))
            if (ret[0] > -0.5): iterf = ret[1]
            elif (ret[0] == -2.): status = 1
            elif (ret[0] == -3.): CPU = ret[1]

    if status == 1:
        iterf = iterf.split('/')[1]
        info[4] = 'DONE - ['+iterf+'] - '+CPU
    else:
        if (RUNPROC == -1): info[4] = 'SCHEDULED' # has been canceled
        else: info[4] = 'FAILED'
    updateListBox(no, info)

    return status

#==============================================================================
# Retourne la liste des cas de calcul situes dans CASES
# Un cas de calcul doit avoir un script run.py
#==============================================================================
def getCases():
    global CASES, INFOS, STATUS
    if CASES != []: return CASES
    path = CASSIOPEE+'/Cases'
    reps = os.listdir(path)
    for i in reps:
        a = os.access(path+'/'+i+'/run.py', os.F_OK)
        if a: 
            CASES.append(i)
            a = os.access(path+'/'+i+'/Readme.txt', os.F_OK)
            if a:
                f = open(path+'/'+i+'/Readme.txt', 'r')
                text = f.read()
                f.close()
                text = text[0:75]
                text = text.replace('\n', ' - ')
                INFOS.append(text)
            else: INFOS.append('No information.')
            a = os.access(path+'/'+i+'/restart.cgns', os.F_OK)
            if a: STATUS.append('RESTART')
            else: STATUS.append('NEW')
    return CASES

#==============================================================================
# Update du menu case
#==============================================================================
def updateCases1(event=None):
    getCases()
    m = WIDGETS['case'].children['menu']
    m.delete(0, TK.END)
    for i in CASES:
        m.add_command(label=i, command=lambda v=VARS[0],l=i:v.set(l))

#==============================================================================
# Update du menu case
#==============================================================================
def updateCases2(event=None):
    getCases()
    if 'case' in WIDGETS:
        WIDGETS['case']['values'] = CASES

#==============================================================================
# Run cases
# Lance les cas de la liste si leur status est scheduled
#==============================================================================
def runCases():
    global RUNPROC
    RUNPROC = 1
    selection = listbox.get(0, TK.END)
    launched = []; no = []
    c = 0
    for s in selection:
        info = getInfoFromString(s)
        if (info[4] == 'SCHEDULED'): launched.append(s); no.append(c)
        c += 1

    while (len(launched)>0):
        if (RUNPROC == -1): break
        runSingleCase(launched[0], no[0])
        # Check list box for update
        selection = listbox.get(0, TK.END)
        launched = []; no = []
        c = 0
        for s in selection:
            info = getInfoFromString(s)
            if (info[4] == 'SCHEDULED'): launched.append(s); no.append(c)
            c += 1

#==============================================================================
def stopRun():
    global RUNPROC
    if (RUNPROC != 1 and RUNPROC != -1): RUNPROC.kill()
    RUNPROC = -1
    return

#==============================================================================
# Lance un seul cas (en local)
#==============================================================================
def runSingleCase(string, no):
    # Case info
    info = getInfoFromString(string)
    info[4] = 'RUNNING (0%)'

    # Update status
    updateListBox(no, info)

    # Launch
    path = CASSIOPEE+'/Cases/'+info[0]
    if (mySystem == 'mingw' or mySystem == 'windows'):
        # Commande Dos (sans time)
        path = path.replace('/', '\\')
        cmd = 'cd '+path+' && python run.py '+info[1]+' '+info[2]+' '+info[3]
        cmd2 = 'echo %time%'
    else:
        # Unix - le shell doit avoir l'environnement cassiopee
        cmd = 'cd '+path+'; time python run.py '+info[1]+' '+info[2]+' '+info[3]
    ret = launchAndSurvey(cmd, no, info)

#==============================================================================
# IN: path: file path
# IN: file: file name
#==============================================================================
def rmFile(path, file):
    if (mySystem == 'mingw' or mySystem == 'windows'):
        path = path.replace('/', '\\')
        cmd = 'cd '+path+' && del '+file
    else:
        cmd = 'cd '+path+'; rm -f '+file
    try:
        subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT)
        subprocess.call(cmd2, shell=True, stderr=subprocess.STDOUT)
    except: pass   
    return

#==============================================================================
# Fait un rm de restart.cgns
#==============================================================================
def cleanCases():
    res = tkMessageBox.askquestion("Clean", "Are You Sure to rm restart.cgns\nfor selected cases?", 
                                   icon='warning')
    if (res == 'no'): return
    path = CASSIOPEE+'/Cases'
    selection = listbox.curselection()
    for s in selection:
        t = listbox.get(s)
        info = getInfoFromString(t)
        rmFile(path+'/'+info[0], 'restart.cgns')

#==============================================================================
def addCase2List():
    string = buildString([VARS[0].get(), VARS[1].get(), VARS[2].get(),
                         VARS[3].get(), 'SCHEDULED'])
    WIDGETS['listbox'].insert(TK.END, string)

#==============================================================================
# IN: info
# info = [case, nit, scale, gfx, status]
# case: nom du cas (string)
# nit: nbre d'iterations a effectuer (string)
# scale: scaling (string)
# gfx: 1 ou 0 (string)
# status: statut du cas (string)
#==============================================================================
def buildString(info):
    case = info[0]
    nit = info[1]
    scale = info[2]
    gfx = info[3]
    status = info[4]
    string = case.ljust(20) + separatorl + nit.ljust(10) + separatorl + \
    scale.ljust(5) + separatorl + gfx.ljust(2) + separatorl + \
    status.ljust(10)
    return string

#==============================================================================
def getInfoFromString(string):
    split1 = string.split(separator)
    case = split1[0].strip()
    nit = split1[1].strip()
    scale = split1[2].strip()
    gfx = split1[3].strip()
    status = split1[4].strip()
    return [case, nit, scale, gfx, status]

#==============================================================================
# Lance l'editeur sur run.py
#==============================================================================
def editCase(event=None):
    path = CASSIOPEE+'/Cases/'
    selection = listbox.curselection()
    for s in selection:
        t = listbox.get(s)
        splits = t.split(separator)
        case = splits[0]
        case = case.strip()
        pathl = path+'/'+case
        if (mySystem == 'mingw' or mySystem == 'windows'):
            pathl = pathl.replace('/', '\\')
            cmd = 'cd '+pathl+' && emacs run.py'
        else:
            cmd = 'cd '+pathl+'; emacs run.py'
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

#==============================================================================
# Archive un cas
#==============================================================================
def archiveCases(event=None):
    path = CASSIOPEE+'/Cases/'
    selection = listbox.curselection()
    for s in selection:
        t = listbox.get(s)
        splits = t.split(separator)
        case = splits[0]
        case = case.strip()
        pathl = path+'/'+case
        dirname = time.strftime("%d-%m-%Y_%H-%M")
        try: os.mkdir(pathl+'/'+dirname)
        except: pass
        import shutil
        shutil.copyfile(pathl+'/run.py', pathl+'/'+dirname+'/run.py')
        shutil.copyfile(pathl+'/restart.cgns', pathl+'/'+dirname+'/restart.cgns')

#==============================================================================
# Remonte le premier selectionne d'un cran
#==============================================================================
def Up(event=None):
    selection = listbox.curselection()
    return

#==============================================================================
def Down(event=None):
    return

#==============================================================================
# Delete une ligne de la liste box
#==============================================================================
def Del(event=None):
    selection = listbox.curselection()
    for s in selection: listbox.delete(s)

#==============================================================================
# Lance la visu sur restart.cgns
#==============================================================================
def viewCase(event=None):
    path = CASSIOPEE+'/Cases/'
    selection = listbox.curselection()
    for s in selection:
        t = listbox.get(s)
        splits = t.split(separator)
        case = splits[0]
        case = case.strip()
        pathl = path+'/'+case
        if (mySystem == 'mingw' or mySystem == 'windows'):
            pathl = pathl.replace('/', '\\')
            cmd = 'cd '+pathl+' && cassiopee restart.cgns'
        else:
            cmd = 'cd '+pathl+'; cassiopee restart.cgns'
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

#==============================================================================
# Main
#==============================================================================

# Main window
master = TK.Tk()
master.title('*Cassiopee* Computation Center')
#GENERALFONT = ('Courier', 9)
GENERALFONT = ('Andale Mono', 9)
master.option_add('*Font', GENERALFONT)

# Main menu
menu = TK.Menu(master)
file = TK.Menu(menu, tearoff=0)
menu.add_cascade(label='File', menu=file)
file.add_command(label='Quit', command=Quit, accelerator='Ctrl+Q')
master.config(menu=menu)
master.bind_all("<Control-q>", Quit)

# - VARS -
# -0- Case -
V = TK.StringVar(master); V.set(''); VARS.append(V)
# -1- Nit
V = TK.StringVar(master); V.set('10'); VARS.append(V)
# -2- scale
V = TK.StringVar(master); V.set('1.'); VARS.append(V)
# -3- gfx
V = TK.StringVar(master); V.set('0'); VARS.append(V)

# Main frame
Frame = TK.Frame(master)
Frame.columnconfigure(0, weight=1)

# list box
listbox = TK.Listbox(Frame, selectmode=TK.EXTENDED, width=120, height=22,
                     background='White')
listbox.grid(row=0, column=0, columnspan=10, sticky=TK.EW)
WIDGETS['listbox'] = listbox
scrollbar = TK.Scrollbar(Frame, orient=TK.VERTICAL)
scrollbar.grid(row=0, column=10, sticky=TK.NSEW)

# - premiere ligne -

# Case label display
F = TK.Frame(Frame, borderwidth=0)
F.columnconfigure(0, weight=1)
B = TK.Entry(F, textvariable=VARS[0], background='White', width=10)
B.grid(sticky=TK.EW)
WIDGETS['case'] = B
#F.bind('<Enter>', updateCases1)
F.grid(row=1, column=0, sticky=TK.EW)

#if ttk is None:
#    B = TK.OptionMenu(F, VARS[0], '')
#    B.grid(sticky=TK.EW)
#    WIDGETS['case'] = B
#    F.bind('<Enter>', updateCases1)
#    F.grid(row=1, column=0, sticky=TK.EW)
#else:
#    B = ttk.Combobox(F, textvariable=VARS[0], 
#                     values=[], state='normal')
#    B.grid(sticky=TK.EW)
#    WIDGETS['case'] = B
#    F.bind('<Enter>', updateCases2)
#    F.grid(row=1, column=0, sticky=TK.EW)

# Open case
B = TK.Button(Frame, text='..', padx=0, pady=0, command=caseSelector)
BB = CTK.infoBulle(parent=B, text='Pick a case.')
B.grid(row=1, column=1)

# nit
B = TK.Entry(Frame, textvariable=VARS[1], background='White', width=9)
BB = CTK.infoBulle(parent=B, text='Number of iterations.')
B.grid(row=1, column=2)

# scale
B = TK.Entry(Frame, textvariable=VARS[2], background='White', width=5)
BB = CTK.infoBulle(parent=B, text='Case scale.')
B.grid(row=1, column=3)

# gfx
B = TK.Checkbutton(Frame, variable=VARS[3])
BB = CTK.infoBulle(parent=B, text='Toggle gfx.')
B.grid(row=1, column=4)

# Add to list
B = TK.Button(Frame, text='Add', command=addCase2List)
B.grid(row=1, column=5)

# - deuxieme ligne -
F = TK.Frame(Frame)
F.grid(row=2, column=0, columnspan=10, sticky=TK.EW)

# Up in list
B = TK.Button(F, text='Up', command=Up, fg='blue')
B.grid(row=0, column=0, sticky=TK.EW)

# Down in list
B = TK.Button(F, text='Down', command=Down, fg='blue')
B.grid(row=0, column=1, sticky=TK.EW)

# Delete from list
B = TK.Button(F, text='Del', command=Del, fg='blue')
BB = CTK.infoBulle(parent=B, text='Del from list')
B.grid(row=0, column=2, sticky=TK.EW)

# Separator
label = TK.Label(F, text='-')
label.grid(row=0, column=3, sticky=TK.EW)

# Run
B = TK.Button(F, text='Run', command=runCases, fg='blue')
BB = CTK.infoBulle(parent=B, text='Run cases')
B.grid(row=0, column=4, sticky=TK.EW)

# Stop
B = TK.Button(F, text='Stop', command=stopRun, fg='red')
BB = CTK.infoBulle(parent=B, text='Stop run')
B.grid(row=0, column=5, sticky=TK.EW)

# Edit run.py
B = TK.Button(F, text='Edit', command=editCase)
BB = CTK.infoBulle(parent=B, text='Edit run.py')
B.grid(row=0, column=6, sticky=TK.EW)

# View
B = TK.Button(F, text='View', command=viewCase)
BB = CTK.infoBulle(parent=B, text='Display restart.cgns')
B.grid(row=0, column=7, sticky=TK.EW)

# Separator
label = TK.Label(F, text='-')
label.grid(row=0, column=8, sticky=TK.EW)

# Clean
B = TK.Button(F, text='Clean', command=cleanCases, fg='red')
BB = CTK.infoBulle(parent=B, text='Rm restart.cgns')
B.grid(row=0, column=10, sticky=TK.EW)

# Archive
B = TK.Button(F, text='Archive', command=archiveCases)
BB = CTK.infoBulle(parent=B, text='Archive case.')
B.grid(row=0, column=9, sticky=TK.EW)

Frame.grid()
TK.mainloop()
