# Validation des solveurs de Fast
import Tkinter as TK
import tkMessageBox
import os, sys, re, time
import subprocess
import CPlot.Tk as CTK
import KCore.Dist as Dist
import KCore.kcore
from threading  import Thread
try: from Queue import Queue, Empty
except ImportError:
    from queue import Queue, Empty
ON_POSIX = 'posix' in sys.builtin_module_names

ttk = CTK.importTtk()

# CASSIOPEE var
CASSIOPEE = os.getenv('CASSIOPEE')
if (CASSIOPEE == ''):
    print 'Error: CASSIOPEE must be present in your environement.'
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
# Retourne la liste des cas de calcul situe dans CASSIOPEE/Validation/Cases
# Un cas de calcul doit avoir un shell valid
#==============================================================================
def getCases():
    global CASES, INFOS, STATUS
    if CASES != []: return CASES
    path = CASSIOPEE+'/Validation/Cases'
    reps = os.listdir(path)
    for i in reps:
        a = os.access(path+'/'+i+'/valid', os.F_OK)
        if (a == True): 
            CASES.append(i)
            a = os.access(path+'/'+i+'/Readme.txt', os.F_OK)
            if (a == True):
                f = open(path+'/'+i+'/Readme.txt', 'r')
                text = f.read()
                f.close()
                text = text[0:75]
                text = text.replace('\n', ' - ')
                INFOS.append(text)
            else: INFOS.append('No information.')
            a = os.access(path+'/'+i+'/restart.cgns', os.F_OK)
            if (a == True): STATUS.append('RESTART')
            else: STATUS.append('NEW')
    return CASES

#==============================================================================
def runCases():
    return

#==============================================================================
def stopCases():
    return

#==============================================================================
def updateCases():
    return

#==============================================================================
def editCase():
    return

#==============================================================================
# Modifie le nbre de threads utilises pour la valid
#==============================================================================
def setThreads(event=None):
    nt = Threads.get()
    try:
        nti = int(nt)
        KCore.kcore.setOmpMaxThreads(nti)
        print 'Num threads set to %d.\n'%nti
    except:
        print 'Bad thread number.\n'
    return

#==============================================================================
# Recupere le nbre de threads (OMP_NUM_THREADS)
#==============================================================================
def getThreads():
    nt = KCore.kcore.getOmpMaxThreads()
    Threads.set(str(nt))
    text.update()
    return

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
# Quit
#==============================================================================
def Quit(event=None):
    import os; os._exit(0)

#==============================================================================
# Main
#==============================================================================

# Main window
master = TK.Tk()
master.title('*Fast* validation')
#GENERALFONT = ('Courier', 9)
GENERALFONT = ('Andale Mono', 9)
master.option_add('*Font', GENERALFONT)

# Main menu
menu = TK.Menu(master)
file = TK.Menu(menu, tearoff=0)
#menu.add_cascade(label='File', menu=file)
file.add_command(label='Quit', command=Quit, accelerator='Ctrl+Q')
master.config(menu=menu)
master.bind_all("<Control-q>", Quit)

# Main frame
frame = TK.Frame(master)
frame.columnconfigure(0, weight=1)
frame.rowconfigure(0, weight=1)
frame.columnconfigure(1, weight=1)
frame.grid(row=0, column=0, sticky=TK.EW)

listbox = TK.Listbox(frame, selectmode=TK.EXTENDED, width=120, height=47,
                     background='White')
listbox.grid(row=0, column=0, columnspan=10, sticky=TK.NSEW)
scrollbar = TK.Scrollbar(frame, orient=TK.VERTICAL)
scrollbar.grid(row=0, column=10, sticky=TK.NSEW)

Status = TK.StringVar(master)
label = TK.Label(frame, textvariable=Status)
Status.set('Stopped'); label.config(bg='red')
label.grid(row=1, column=0, sticky=TK.EW)

# Check
button = TK.Button(frame, text='Run check', command=runCases, fg='blue')
BB = CTK.infoBulle(parent=button, text='Run selected case from restart.')
button.grid(row=1, column=3, sticky=TK.EW)
button = TK.Button(frame, text='Update check', command=runCases, fg='blue')
BB = CTK.infoBulle(parent=button, text='Run selected case from restart.')
button.grid(row=1, column=4, sticky=TK.EW)

# Stop
button = TK.Button(frame, text='Stop', command=stopCases, fg='red')
button.grid(row=1, column=5, sticky=TK.EW)

# Full
button = TK.Button(frame, text='Run full', command=runCases, fg='blue')
BB = CTK.infoBulle(parent=button, text='Run selected case from scratch.')
button.grid(row=1, column=6, sticky=TK.EW)
button = TK.Button(frame, text='Update full', command=updateCases, fg='blue')
BB = CTK.infoBulle(parent=button, text='Update case (replace data base files).')
button.grid(row=1, column=7, sticky=TK.EW)

# Edit
button = TK.Button(frame, text='Edit', command=editCase)
button.grid(row=1, column=8, sticky=TK.EW)

Threads = TK.StringVar(master)
text = TK.Entry(frame, textvariable=Threads, background='White', width=3)
text.grid(row=1, column=9, sticky=TK.EW)
text.bind('<Return>', setThreads)
getThreads()
BB = CTK.infoBulle(parent=text, text='Number of threads.')
frame.grid(sticky=TK.NSEW)

# Recupere les cas
getCases()



TK.mainloop()
