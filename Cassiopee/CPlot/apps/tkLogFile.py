# - Log file management -
try: import tkinter as TK
except: import Tkinter as TK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Log as Log
import os

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# Start a new recording session
#==============================================================================
def record():
    dir = os.path.dirname(CTK.FILE)
    fileName = VARS[0].get()
    fileName = os.path.join(dir, fileName)
    Log.openLogFile(fileName)
    CTK.TXT.insert('START', 'Recording session to %s.\n'%fileName)

#==============================================================================
# Pause recording session
#==============================================================================
def pause():
    Log.pause()
    CTK.TXT.insert('START', 'Recording session paused.\n')

#==============================================================================
# Stop recording session
#==============================================================================
def stop():
    Log.closeLogFile()
    CTK.TXT.insert('START', 'Recording session ended.\n')

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                          text='tkLogFile  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Manage python log file.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    WIDGETS['frame'] = Frame

    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkLogFile')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- File name -
    V = TK.StringVar(win); V.set('SessionLog.py'); VARS.append(V)

    # - FileName -
    B = TK.Entry(Frame, textvariable=VARS[0], background='White')
    BB = CTK.infoBulle(parent=B, text='Log file name.')
    B.grid(row=0, column=0, columnspan=3, sticky=TK.EW)

    # - Record -
    B = TK.Button(Frame, text="Record", command=record)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Start a new recording session.')
    B = TK.Button(Frame, text="Pause", command=pause)
    B.grid(row=1, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Pause recording.')
    B = TK.Button(Frame, text="Stop", command=stop)
    B.grid(row=1, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Stop recording session.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    WIDGETS['frame'].grid(sticky=TK.NSEW)

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
    (win, menu, file, tools) = CTK.minimal('tkLogFile '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
