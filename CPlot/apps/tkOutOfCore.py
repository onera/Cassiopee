# - out/in core memory management -
try: import Tkinter as TK
except: import tkinter as TK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal
import Generator.PyTree as G
import numpy

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def inCore():
    if CTK.t == []: return
    if (CTK.__MAINTREE__ <= 0):
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        ooc = Internal.getNodesFromName1(z, 'OutOfCore')
        if (len(ooc) != 0):
            # Read zone
            base = CTK.t[2][nob]
            name = '.'+base[0]+'#'+z[0]+'.cgns'
            t = C.convertFile2PyTree(name)
            CTK.replace(CTK.t, nob, noz, t[2][1][2][0])
    
    CTK.saveTree() # il est apres pour forcer le flush
    CTK.TXT.insert('START', 'Selected zones in core.\n')
    CTK.t = C.fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    
#==============================================================================
def outCore():
    if (CTK.t == []): return
    if (CTK.__MAINTREE__ <= 0):
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if (nzs == []):
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        ooc = Internal.getNodesFromName1(z, 'OutOfCore')
        if (len(ooc) == 0):
            # Save zone
            base = CTK.t[2][nob]
            name = '.'+base[0]+'#'+z[0]+'.cgns'
            t = C.newPyTree(['Base'])
            t[2][1][2].append(z)
            C.convertPyTree2File(t, name)
            # Replace zone
            bb = G.BB(z)
            bb[2].append(['OutOfCore', numpy.array([1]), [], \
                          'UserDefinedData_t'])
            bb = CPlot.addRender2Zone(bb, blending=0.2)
            CTK.replace(CTK.t, nob, noz, bb)
    CTK.saveTree() # il est apres pour forcer le flush
    CTK.TXT.insert('START', 'Selected zones out of core.\n')
    CTK.t = C.fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
 
#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                          text='tkOutOfCore', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Push blocks out/in core memory.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    WIDGETS['frame'] = Frame

    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkOutOfCore')
    WIDGETS['frameMenu'] = FrameMenu
    
    # - IN/OUT -
    B = TK.Button(Frame, text="IN Core", command=inCore)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Get the selection in core memory.')
    B = TK.Button(Frame, text="OUT of Core", command=outCore)
    B.grid(row=0, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Get the selection out of core memory.')
    
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
    (win, menu, file, tools) = CTK.minimal('tkOutOfCore '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
