# - tkNodeEdit -
"""Edit pyTree nodes."""
import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal
import Converter.Filter as Filter
import Compressor.PyTree as Compressor
import numpy

# local widgets list
WIDGETS = {}; VARS = []

# format string pour ecriture casual (6)
def strFormat(v):
    if isinstance(v, int): return "%ld"%v
    elif isinstance(v, float): return "%.5g"%v
    else: return "%s"%v

# format string pour ecriture precise (12)
def strFormat2(v):
    if isinstance(v, int): return "%ld"%v
    elif isinstance(v, float): return "%.12g"%v
    else: return "%s"%v

#==============================================================================
# Affiche la valeur de noeud et le type du noeud
#==============================================================================
def showNodeValue(event=None):
    if CTK.t == []: return
    node = CTK.TKTREE.getCurrentSelectedNode()
    if node == []:
        index = 0; VARS[1].set(str(index))
        VARS[2].set(''); VARS[3].set(''); VARS[4].set(''); return

    VARS[3].set(node[3]); VARS[4].set(node[0])

    index = VARS[1].get()
    try: index = int(index)
    except: index = 0; VARS[1].set('0')
    if index < 0: index = 0; VARS[1].set(str(index))
    v = node[1]
    flatView2 = [None, None, None, None] # for display info window: [dim, start, index, last]
    if isinstance(v, float): v = strFormat2(v); flatView2[2] = v; index = 0
    elif isinstance(v, int): v = strFormat2(v); flatView2[2] = v; index = 0
    elif isinstance(v, numpy.ndarray):
        if v.dtype == 'c':
            v = Internal.getValue(node); flatView2[2] = v; index = 0
        else:
            sh = v.shape
            if len(sh) == 3:
                dim = '('+str(sh[0])+'x'+str(sh[1])+'x'+str(sh[2])+'): '
            elif len(sh) == 2:
                dim = '('+str(sh[0])+'x'+str(sh[1])+'): '
            elif len(sh) == 1:
                dim = '('+str(sh[0])+'): '
            else:
                dim = str(sh)+': '

            try: vf = v.ravel(order='F')
            except: vf = v.flat
            if index >= v.size: index = v.size-1; VARS[1].set(str(index))
            # Construit une chaine representant le tableau a plat
            if index > 3:
                if v.size <= index+1: flatView = ''
                elif v.size <= index+2:
                    flatView = ' '+strFormat(vf[index+1])+''
                elif v.size <= index+3:
                    flatView = ' '+strFormat(vf[index+1])+' '+strFormat(vf[index+2])+''
                else:
                    flatView = ' '+strFormat(vf[index+1])+' '+strFormat(vf[index+2])+'...'
                flatView += '\n'
                CTK.TXT.insert('START', flatView)
                flatView = strFormat(vf[index])
                CTK.TXT.insert('START', flatView, 'Warning')
                flatView = dim+'...'+strFormat(vf[index-2])+' '+strFormat(vf[index-1])+' '
                flatView2[0] = dim
                CTK.TXT.insert('START', flatView)

            elif index == 2:
                if v.size <= index+1: flatView = ''
                elif v.size <= index+2:
                    flatView = ' '+strFormat(vf[index+1])+''
                elif v.size <= index+3:
                    flatView = ' '+strFormat(vf[index+1])+' '+strFormat(vf[index+2])+''
                else:
                    flatView = ' '+strFormat(vf[index+1])+' '+strFormat(vf[index+2])+'...'
                flatView += '\n'
                CTK.TXT.insert('START', flatView)
                flatView = strFormat(vf[index])
                CTK.TXT.insert('START', flatView, 'Warning')
                flatView = dim+strFormat(vf[index-2])+' '+strFormat(vf[index-1])+' '
                CTK.TXT.insert('START', flatView)

            elif index == 1:
                if v.size <= index+1: flatView = ''
                elif v.size <= index+2:
                    flatView = ' '+strFormat(vf[index+1])+''
                elif v.size <= index+3:
                    flatView = ' '+strFormat(vf[index+1])+' '+strFormat(vf[index+2])+''
                else:
                    flatView = ' '+strFormat(vf[index+1])+' '+strFormat(vf[index+2])+'...'
                flatView += '\n'
                CTK.TXT.insert('START', flatView)
                flatView = strFormat(vf[index])
                CTK.TXT.insert('START', flatView, 'Warning')
                flatView = dim+strFormat(vf[index-1])+' '
                CTK.TXT.insert('START', flatView)

            elif index == 0:
                if v.size <= index+1: flatView = ''
                elif v.size <= index+2:
                    flatView = ' '+strFormat(vf[index+1])+''
                elif v.size <= index+3:
                    flatView = ' '+strFormat(vf[index+1])+' '+strFormat(vf[index+2])+''
                else:
                    flatView = ' '+strFormat(vf[index+1])+' '+strFormat(vf[index+2])+'...'
                flatView += '\n'
                CTK.TXT.insert('START', flatView)
                flatView = strFormat(vf[index])
                CTK.TXT.insert('START', flatView, 'Warning')
                flatView = dim
                CTK.TXT.insert('START', flatView)

            v = strFormat2(vf[index])

            SIZE = 10
            flatView2[0] = dim
            if index-SIZE//2 >= 0 and index+SIZE//2 < vf.size:
                left = index - SIZE//2
                right = index + SIZE//2
            elif index - SIZE//2 < 0:
                left = max(0, index-SIZE//2)
                right = min(SIZE - left, vf.size)
            else:
                right = min(index+SIZE//2, vf.size)
                left = max(right - SIZE, 0)

            rep = ''
            for i in range(left, index):
                rep += strFormat(vf[i])+' '
            if rep != '':
                if left > 0: rep = '...'+rep
                flatView2[1] = rep

            rep = strFormat(vf[index])+' '
            flatView2[2] = rep

            rep = ''
            for i in range(index+1, right):
                rep += strFormat(vf[i])+' '

            if rep != '':
                if right < vf.size: rep += '...'
                flatView2[3] = rep

    else: # fall back
        v = str(node[1]); index = 0

    VARS[1].set(str(index))
    VARS[2].set(v)

    # update info window
    #pos = WIDGETS['infoText'].index(TK.INSERT); print(pos)
    WIDGETS['infoText'].delete(1.0, TK.END)
    if flatView2[3] is not None: # last
        WIDGETS['infoText'].insert('START', flatView2[3])
    if flatView2[2] is not None: # active
        WIDGETS['infoText'].insert('START', flatView2[2], 'Active')
    if flatView2[1] is not None: # start
        WIDGETS['infoText'].insert('START', flatView2[1])
    if flatView2[0] is not None: # dim
        WIDGETS['infoText'].insert('START', flatView2[0]+'\n', 'Dim')
    #WIDGETS['infoText'].mark_set("insert", pos)

#==============================================================================
# increase index in inspector
#==============================================================================
def increaseIndex(event=None):
    if CTK.t == []: return
    index = VARS[1].get()
    index = int(index)
    index += 1
    VARS[1].set(str(index))
    showNodeValue()

#==============================================================================
# decrease index in inspector
#==============================================================================
def decreaseIndex(event=None):
    if CTK.t == []: return
    index = VARS[1].get()
    index = int(index)
    index -= 1
    VARS[1].set(str(index))
    showNodeValue()

#==============================================================================
# Changement de la valeur du noeud
#==============================================================================
def setNodeValue(event=None):
    if CTK.t == []: return
    node = CTK.TKTREE.getCurrentSelectedNode()
    if node == []: return
    index = VARS[1].get()
    try: index = int(index)
    except: index = 0
    val = VARS[2].get()
    v = node[1]
    if isinstance(v, float):
        Internal.setValue(node, val)
        CTK.TXT.insert('START', 'Value "%d" set in node.\n.'%val) # not done by showNodeValue
    elif isinstance(v, int):
        Internal.setValue(node, val)
        CTK.TXT.insert('START', 'Value "%g" set in node.\n.'%val) # not done by showNodeValue
    elif isinstance(v, numpy.ndarray):
        if v.dtype == 'c':
            Internal.setValue(node, val)
            CTK.TXT.insert('START', 'Value "%s" set in node.\n.'%val) # not done by showNodeValue
        elif v.dtype == 'int32':
            try: vf = v.ravel(order='F')
            except: vf = v.flat
            vf[index] = int(val)
            #Internal.setValue(node, node[1])
        else:
            try: vf = v.ravel(order='F')
            except: vf = v.flat
            vf[index] = float(val)
            #Internal.setValue(node, node[1])
            CTK.display(CTK.t) # view may be impacted
    showNodeValue()

#==============================================================================
# Change le type du noeud. Ca peut etre tres dangereux.
#==============================================================================
def setNodeType(event=None):
    if CTK.t == []: return
    node = CTK.TKTREE.getCurrentSelectedNode()
    if node == []: return
    node[3] = VARS[3].get()

#==============================================================================
# Change le nom du noeud.
#==============================================================================
def setNodeName(event=None):
    if CTK.t == []: return
    node = CTK.TKTREE.getCurrentSelectedNode()
    if node == []: return
    node[0] = VARS[4].get()
    CTK.TKTREE.updateApp()

#==============================================================================
# Appele de tkTree, pour informer qu'un noeud a ete selectionne
#==============================================================================
def updateNode(node):
    showNodeValue()

#==============================================================================
def setByLevel(node1, node2, depth, maxDepth):
    #print("setting", node1[0])
    zdata = Internal.getNodeFromName1(node2, 'ZData')
    if zdata is not None: Compressor._unpackNode(node2)
    node1[1] = node2[1]
    if depth < maxDepth or maxDepth == -1:
        for c, n in enumerate(node1[2]):
            if len(node2[2]) > c:
                setByLevel(n, node2[2][c], depth+1, maxDepth)

#==============================================================================
def freeByLevel(node1, depth, maxDepth):
    #print("freeing", node1[0])
    node1[1] = None
    if depth < maxDepth or maxDepth == -1:
        for n in node1[2]:
            freeByLevel(n, depth+1, maxDepth)

#==============================================================================
def loadNode():
    if CTK.t == []: return
    if CTK.HANDLE is None: fileName = CTK.FILE
    else: fileName = CTK.HANDLE.fileName
    node = CTK.TKTREE.getCurrentSelectedNode()
    if node == []: return
    path = Internal.getPath(CTK.t, node)
    if path is None: return
    depth = VARS[7].get()
    depth = int(depth)

    # check if node is compressed
    zdata = Internal.getNodeFromName1(node, 'ZData')
    if zdata is not None:
        depth += 1 # load zdata also

    # read node
    nodes = Filter.readNodesFromPaths(fileName, [path], maxDepth=depth)

    # check if node is compressed
    if zdata is not None:
        depth -= 1

    # depth replace
    setByLevel(node, nodes[0], 0, depth)

    updateNode(node)

#==============================================================================
def freeNode():
    if CTK.t == []: return
    node = CTK.TKTREE.getCurrentSelectedNode()
    if node == []: return
    depth = VARS[7].get()
    depth = int(depth)
    freeByLevel(node, 0, depth)
    updateNode(node)

#==============================================================================
def saveNode():
    if CTK.t == []: return
    if CTK.HANDLE is None: fileName = CTK.FILE
    else: fileName = CTK.HANDLE.fileName
    node = CTK.TKTREE.getCurrentSelectedNode()
    if node == []: return
    path = Internal.getPath(CTK.t, node)
    depth = VARS[7].get()
    depth = int(depth)
    Filter.writeNodesFromPaths(fileName, [path], [node], maxDepth=depth, mode=1)

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):

    ttk = CTK.importTtk()

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkNodeEdit  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)

    #BB = CTK.infoBulle(parent=Frame, text='Edit pyTree nodes.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    Frame.columnconfigure(2, weight=1)
    Frame.columnconfigure(3, weight=0)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkNodeEdit')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Base to move to -
    V = TK.StringVar(win); V.set('New Base'); VARS.append(V)
    # -1- Displayed index
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -2- Displayed node value
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -3- Displayed node type
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -4- Displayed node name
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -5- Name/Type to be deleted
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -6- Type of deletion
    V = TK.StringVar(win); V.set('Name'); VARS.append(V)
    # -7- Maxdepth for load
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -8- display text value
    V = TK.StringVar(win); VARS.append(V)

    # - Small node editor -
    F = TTK.Frame(Frame, borderwidth=0, relief=CTK.FRAMESTYLE, takefocus=1)
    F.columnconfigure(0, weight=1)
    F.columnconfigure(1, weight=4)
    F.grid(row=0, column=0, columnspan=4, sticky=TK.EW)
    B = TTK.Entry(F, textvariable=VARS[4], background='White', width=5)
    B.grid(row=0, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Selected node name.\nYou can change this by pressing <Enter>.')
    B.bind('<Return>', setNodeName)

    B = TTK.Entry(F, textvariable=VARS[3], background='White', width=5)
    B.grid(row=0, column=1, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Selected node type.\nYou can change this by pressing <Enter>.')
    B.bind('<Return>', setNodeType)

    FrameA = TK.Frame(F)
    FrameA.columnconfigure(0, weight=0)
    FrameA.columnconfigure(1, weight=1)
    FrameA.columnconfigure(2, weight=0)
    FrameA.grid(row=1, column=0, sticky=TK.EW)
    B = TTK.Entry(FrameA, textvariable=VARS[1], background='White', width=5)
    B.grid(row=0, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Displayed index (start: 0).')
    B.bind('<Return>', showNodeValue)
    B = TTK.Button(FrameA, text="<", command=decreaseIndex, width=2)
    B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.Button(FrameA, text=">", command=increaseIndex, width=2)
    B.grid(row=0, column=2, sticky=TK.EW)

    B = TTK.Entry(F, textvariable=VARS[2], background='White', width=5, foreground="green")
    B.grid(row=1, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Selected node value at given index.\nYou can change this by pressing <Enter>.')
    B.bind('<Return>', setNodeValue)

    # display info text
    B = TTK.Text(Frame, background='White', width=5, height=5)
    B.grid(row=1, column=0, columnspan=4, sticky=TK.EW)
    B.tag_config("Dim", foreground="blue")
    B.tag_config("Active", foreground="green")
    B.mark_set('START', TK.INSERT)
    B.mark_gravity('START', TK.LEFT)
    WIDGETS['infoText'] = B

    # - Load/free Nodes -
    B = TTK.Button(Frame, text="Load node", command=loadNode)
    BB = CTK.infoBulle(parent=B, text='Load selected node from file.')
    B.grid(row=2, column=0, sticky=TK.EW)
    B = TTK.Button(Frame, text="Free", command=freeNode)
    BB = CTK.infoBulle(parent=B, text='Free selected node.')
    B.grid(row=2, column=1, columnspan=1, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[7], background='White', width=2)
    B.grid(row=2, column=2, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Max depth for load, save and free (-1: full).')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['TreeNoteBook'].add(WIDGETS['frame'], text='tkNodeEdit')
    except: pass
    CTK.WIDGETS['TreeNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['TreeNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(): return

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
    (win, menu, file, tools) = CTK.minimal('tkNodeEdit '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
