# - Gestion des familles -
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# Cree un liste des family zone names (optionMenu)
#==============================================================================
def updateFamilyZoneNameList(event=None):
    if CTK.t == []: return
    varsl = C.getFamilyZoneNames(CTK.t)
    m = WIDGETS['zones'].children['menu']
    m.delete(0, TK.END)
    for i in varsl:
        m.add_command(label=i, command=lambda v=VARS[2],l=i:v.set(l))
        
#==============================================================================
# Cree une list des family zone names (combobox)
#==============================================================================
def updateFamilyZoneNameList2(event=None):
    if CTK.t == []: return
    varsl = C.getFamilyZoneNames(CTK.t)
    varsl = list(set(varsl))
    if 'zones' in WIDGETS:
        WIDGETS['zones']['values'] = varsl

#==============================================================================
def tagWithZoneFamily(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    name = VARS[2].get()
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
        
    CTK.saveTree()
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        CTK.t[2][nob][2][noz] = C.tagWithFamily(z, name)
    CTK.TXT.insert('START', 'Zones are tagged with '+name+'.\n')
    CTK.TKTREE.updateApp()

#==============================================================================
def createZoneFamily(event=None):
    if CTK.t == []: return
    name = VARS[0].get()
    if name == '':
        CTK.TXT.insert('START', 'FamilyZone name is invalid.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    CTK.saveTree()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        bases = CTK.t[2]
        for b in bases:
            if b[3] == 'CGNSBase_t':
                C._addFamily2Base(b, name)
        CTK.TXT.insert('START', 'Zone Family '+name+' added to all bases.\n')
    else:
        nob = CTK.Nb[nzs[0]]+1
        noz = CTK.Nz[nzs[0]]
        z = CTK.t[2][nob][2][noz]
        C._addFamily2Base(CTK.t[2][nob], name)
        CTK.TXT.insert('START', 'Zone Family '+name+' added to base '+CTK.t[2][nob][0]+'.\n')
    CTK.TKTREE.updateApp()

#==============================================================================
def createBCFamily(event=None):
    if CTK.t == []: return
    name = VARS[1].get()
    if name == '':
        CTK.TXT.insert('START', 'FamilyZone name is invalid.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    CTK.saveTree()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        bases = CTK.t[2]
        for b in bases:
            if b[3] == 'CGNSBase_t':
                C._addFamily2Base(b, name, 'UserDefined')
        CTK.TXT.insert('START', 'BC Family '+name+' added to all bases.\n')
    else:
        nob = CTK.Nb[nzs[0]]+1
        noz = CTK.Nz[nzs[0]]
        z = CTK.t[2][nob][2][noz]
        C._addFamily2Base(CTK.t[2][nob], name, 'UserDefined')
        CTK.TXT.insert('START', 'BC Family '+name+' added to base '+CTK.t[2][nob][0]+'.\n')
    CTK.TKTREE.updateApp()

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):

    ttk = CTK.importTtk()

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkFamily', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Create famlies of\nblocks or BCs.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=4)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkFamily')
    WIDGETS['frameMenu'] = FrameMenu
    
    # - VARS -
    # -0- Nom de la FamilyZone (new) -
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -1- Nom de la FamilyBC (new) -
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -2- Nom de la famille zone pour le tag
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    
    # - Create zone family name -
    B = TTK.Button(Frame, text="NewZoneFamily", command=createZoneFamily)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Create a new zone family tag.')
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=15)
    BB = CTK.infoBulle(parent=B, text='Zone family name.')
    B.grid(row=0, column=1, sticky=TK.EW)
    B.bind('<Return>', createZoneFamily)

    # - Tag with zone family -
    B = TTK.Button(Frame, text="Tag zone", command=tagWithZoneFamily)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Tag a zone with a zone family.')
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)

    if ttk is None:
        B = TK.OptionMenu(F, VARS[2], '')
        B.grid(sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Family zone name.')
        F.bind('<Enter>', updateFamilyZoneNameList)
        F.grid(row=1, column=1, sticky=TK.EW)
        WIDGETS['zones'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[2], 
                         values=[], state='readonly')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateFamilyZoneNameList2)
        F.grid(row=1, column=1, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Family zone name.')
        WIDGETS['zones'] = B

    # - Create BC family -
    B = TTK.Button(Frame, text="NewBCFamily", command=createBCFamily)
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Create a new BC family tag.')
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=15)
    BB = CTK.infoBulle(parent=B, text='BC family name.')
    B.grid(row=2, column=1, sticky=TK.EW)
    B.bind('<Return>', createBCFamily)
    
#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    WIDGETS['frame'].grid(sticky=TK.EW)

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    WIDGETS['frame'].grid_forget()

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(): return

#==============================================================================
def displayFrameMenu(event=None):
    WIDGETS['frameMenu'].tk_popup(event.x_root+50, event.y_root, 0)
    
#==============================================================================
if (__name__ == "__main__"):
    import sys
    if (len(sys.argv) == 2):
        CTK.FILE = sys.argv[1]
        try:
            CTK.t = C.convertFile2PyTree(CTK.FILE)
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.display(CTK.t)
        except: pass

    # Main window
    (win, menu, file, tools) = CTK.minimal('tkFamily '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
