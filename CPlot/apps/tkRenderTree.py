# - tkRenderTree -
# Add global rendering options to tree
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import Converter.Internal as Internal
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def addBillboardFile(event=None):
    if CTK.t == []: return
    v = VARS[0].get()
    v.replace(' ', '')
    v = v.split(';')
    vo = []
    for i in v:
        if i != '': vo.append(i) 
    CPlot._addRender2PyTree(CTK.t, billBoards=vo)
    CTK.TKTREE.updateApp()
    
    renderInfo = Internal.getNodeFromName1(CTK.t, '.RenderInfo')
    pos = Internal.getNodeFromName1(renderInfo, 'billBoards')
    out = []
    for i in pos[2]: out.append(Internal.getValue(i))
    CPlot.setState(billBoards=out, billBoardSize=0.8)

#==============================================================================
def addBumpMapFile(event=None):
    if CTK.t == []: return
    v = VARS[0].get()
    v.replace(' ', '')
    v = v.split(';')
    vo = []
    for i in v:
        if i != '': vo.append(i) 
    CPlot._addRender2PyTree(CTK.t, bumpMaps=vo)
    CTK.TKTREE.updateApp()
    
    renderInfo = Internal.getNodeFromName1(CTK.t, '.RenderInfo')
    pos = Internal.getNodeFromName1(renderInfo, 'bumpMaps')
    out = []
    for i in pos[2]: out.append(Internal.getValue(i))
    CPlot.setState(bumpMaps=out)
    
#==============================================================================
def addTextureFile(event=None):
    if CTK.t == []: return
    v = VARS[0].get()
    v.replace(' ', '')
    v = v.split(';')
    vo = []
    for i in v:
        if i != '': vo.append(i) 
    CPlot._addRender2PyTree(CTK.t, materials=vo)
    CTK.TKTREE.updateApp()
    
    renderInfo = Internal.getNodeFromName1(CTK.t, '.RenderInfo')
    pos = Internal.getNodeFromName1(renderInfo, 'materials')
    out = []
    for i in pos[2]: out.append(Internal.getValue(i))
    CPlot.setState(materials=out)
    
#==============================================================================
def chooseFile():
    try: import tkFileDialog
    except: import tkinter.filedialog as tkFileDialog
    import os.path
    init = VARS[0].get()
    init = init.split(';')[0]
    files = tkFileDialog.askopenfilenames(
        filetypes=CTK.fileTypes, initialfile=init, multiple=1)
    if files == '' or files == None or files == (): # user cancel
        return
    # strangely, initfile is part of the return
    files = CTK.fixFileString__(files, init)
    s = ''
    for f in files:
        bname = os.path.basename(f)
        if os.path.exists(bname): s += bname+';' # short
        else: s += f+';' # full name
    s = s[:-1]
    VARS[0].set(s)

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkRenderTree', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Add global rendering options to tree.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkRenderTree')
    WIDGETS['frameMenu'] = FrameMenu
    
    # - VARS -
    # -0- Name of file to add -
    V = TK.StringVar(win); V.set('file.png'); VARS.append(V)
    
    # - Buttons -
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    F.columnconfigure(1, weight=0)
    B = TTK.Entry(F, textvariable=VARS[0], background='White')
    BB = CTK.infoBulle(parent=B, text='Name of image file to add as billboard or texture.')
    B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.Button(F, text="...", padx=0, command=chooseFile)
    BB = CTK.infoBulle(parent=B, text='Select an image file (png).')
    B.grid(row=0, column=1, sticky=TK.EW)
    F.grid(row=0, column=0, sticky=TK.EW)
    
    B = TTK.Button(Frame, text="Add billboard", command=addBillboardFile)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='Add file name as billboard file to be referenced\n in Sphere material of tkRenderSet.')
    B = TTK.Button(Frame, text="Add Texture", command=addTextureFile)
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='Add file name as texture file to be referenced\n in Texmat material of tkRenderSet..')
    B = TTK.Button(Frame, text="Add BumpMap", command=addBumpMapFile)
    B.grid(row=3, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='Add file name as bumpmap file to be referenced\n in Texmat material of tkRenderSet..')
    
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
    if len(sys.argv) == 2:
        CTK.FILE = sys.argv[1]
        try:
            CTK.t = C.convertFile2PyTree(CTK.FILE)
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.display(CTK.t)
        except: pass

    # Main window
    (win, menu, file, tools) = CTK.minimal('tkRenderTree '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
