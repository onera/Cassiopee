# - tkStereo -
"""Activate stereo rendering."""
import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def setStereo(event=None):
    if CTK.t == []: return
    mode = VARS[0].get()
    if mode == 'None': CPlot.setState(stereo=0)
    elif mode == 'Anaglyph (b&w)': CPlot.setState(stereo=1)
    elif mode == 'Anaglyph (color)': CPlot.setState(stereo=2)

#==============================================================================
def setDist(event=None):
    if CTK.t == []: return
    dist = WIDGETS['dist'].get() / 100.
    VARS[1].set('Stereo distance [%.2f].'%dist)
    CPlot.setState(stereoDist=dist)

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkStereo  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Set anaglyph mode.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=0)
    Frame.columnconfigure(1, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkStereo')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Stero mode -
    V = TK.StringVar(win); V.set('None'); VARS.append(V)
    if 'tkStereoMode' in CTK.PREFS: V.set(CTK.PREFS['tkStereoMode'])
    # -1- Stereo dist info bulle
    V = TK.StringVar(win); V.set('Stereo distance.'); VARS.append(V)

    # - Stereo mode -
    B = TTK.Label(Frame, text="Stereo")
    BB = CTK.infoBulle(parent=B, text='Choose a stereo mode.')

    B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.OptionMenu(Frame, VARS[0], 'None', 'Anaglyph (b&w)',
                       'Anaglyph (color)', command=setStereo)
    B.grid(row=0, column=1, sticky=TK.EW)

    # - Stereo Dist -
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL, showvalue=0,
                  borderwidth=1, command=setDist, value=50)
    WIDGETS['dist'] = B
    B.grid(row=1, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, textVariable=VARS[1])

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['RenderNoteBook'].add(WIDGETS['frame'], text='tkStereo')
    except: pass
    CTK.WIDGETS['RenderNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['RenderNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(): return

#==============================================================================
def saveApp():
    CTK.PREFS['tkStereoMode'] = VARS[0].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('None')
    CTK.PREFS['tkStereoMode'] = VARS[0].get()
    CTK.savePrefFile()

#==============================================================================
def displayFrameMenu(event=None):
    WIDGETS['frameMenu'].tk_popup(event.x_root+50, event.y_root, 0)

#==============================================================================
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 2:
        CTK.FILE = sys.argv[1]
        try:
            CTK.t = C.convertFile2PyTree(CTK.FILE)
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.display(CTK.t)
        except: pass

    # Main window
    (win, menu, file, tools) = CTK.minimal('tkStereo '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
