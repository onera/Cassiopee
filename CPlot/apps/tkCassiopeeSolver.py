# - Cassiopee solver app -
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal
import os, math, os.path
# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# Write setup file
#==============================================================================
def writeSetupFile():
    if CTK.t == []: return
    
    # EquationDimension
    nodes = Internal.getNodesFromName(CTK.t, 'EquationDimension')
    if nodes != []:
        dim = nodes[0][1][0]
    else:
        CTK.TXT.insert('START', 'EquationDimension is missing (tkState).\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return

    # GoverningEquations
    nodes = Internal.getNodesFromName(CTK.t, 'GoverningEquations')
    if nodes != []:
        equations = Internal.getValue(nodes[0])
        if equations == 'Euler': model = 'euler'
        if equations == 'NSLaminar': model = 'nslam'
        if equations == 'NSTurbulent': model = 'nstur'
    else:
        CTK.TXT.insert('START', 'GoverningEquations is missing (tkState).\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return

    # Turbulence model
    if (model == 'nstur'):
        nodes = Internal.getNodesFromName(CTK.t, 'TurbulenceModel')
        if (nodes != []):
            tm = Internal.getValue(nodes[0])
            if (tm == 'OneEquation_SpalartAllmaras'): model += '_sa'
            elif (tm == 'TwoEquation_Wilcox'): model += '_kw'
            elif (tm == 'TwoEquation_MenterSST'): model += '_kw'
            else:
                CTK.TXT.insert('START', 'This turbulence model is not accepted by Cassiopee solver.\n')
                CTK.TXT.insert('START', 'Error: ', 'Error')
                return
    
    # ReferenceState
    nodes = Internal.getNodesFromName(CTK.t, 'ReferenceState')
    if (nodes == []):
        CTK.TXT.insert('START', 'Reference state is missing (tkState).\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return
    state = nodes[0]

    # Mach
    nodes = Internal.getNodesFromName(state, 'Mach')
    if (nodes != []): Mach = Internal.getValue(nodes[0])
    else:
        CTK.TXT.insert('START', 'Mach is missing (tkState).\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return
    
    # Reynolds
    nodes = Internal.getNodesFromName(state, 'Reynolds')
    if (nodes != []): Reynolds = Internal.getValue(nodes[0])
    elif (equations == 'NSLaminar' or equations == 'NSTurbulent'):
        CTK.TXT.insert('START', 'Reynolds is missing (tkState).\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return
    else: Reynolds = 1.
    if (Reynolds <= 0.): Reynolds = 1.
    
    # Incidences
    nodes = Internal.getNodesFromName(state, 'VelocityX')
    if (nodes != []): UInf = Internal.getValue(nodes[0])
    else: UInf = 0.
    nodes = Internal.getNodesFromName(state, 'VelocityY')
    if (nodes != []): VInf = Internal.getValue(nodes[0])
    else: VInf = 0.
    nodes = Internal.getNodesFromName(state, 'VelocityZ')
    if (nodes != []): WInf = Internal.getValue(nodes[0])
    else: WInf = 0.
    if (UInf != 0.):
        aly = math.atan(WInf/UInf)
        alz = math.atan( math.cos(aly)*VInf/UInf )
    else:
        aly = 0.; alz = 0.
    alphaZ = alz*180./math.pi
    alphaY = aly*180./math.pi
    
    treeFile = os.path.splitext(CTK.FILE)[0]+'.cgns'
    
    f = open('dump.py', 'w')
    st =  'from Cassiopee import dump\n'
    st += 'import Converter.PyTree as C\n'
    st += 'FILE = \''+treeFile+'\'\n'
    st += 't = C.convertFile2PyTree(FILE, skeleton=1)\n'
    st += 'dump(t, FILE)\n'
    f.write(st)
    f.close()
    
    f = open('setup.py', 'w')
    st =  'import Cassiopee as K\n'
    st += 'import elsA_user as E\n'
    st += 'import Converter.PyTree as C\n'
    st += 'import Converter.Mpi as Cmpi\n'
    st += 'rank = Cmpi.rank; size = Cmpi.size\n'
    st += '# -- Input file --\n'
    st += 'FILE = \''+treeFile+'\'\n'
    st += '# -- Define fluid model --\n'
    st += 'MODEL = \''+model+'\'\n'
    st += '# -- Define integration --\n'
    st += 'INTEG = \''+VARS[0].get()+'\'\n'
    st += '# -- Define numerical scheme --\n'
    st += 'SCHEME = \''+VARS[1].get()+'\'\n'
    f.write(st)

    st =      '\n'
    st += '# -- Define if the computation is a restart --\n'
    st += '# If RESTART = \'active\', the output.plt file and the\n'
    st += '# refinementIndicator file\n'
    st += '# must correspond to the result of the previous computation.\n'
    v = int(VARS[2].get())
    if (v == 0 or v == 1):
        st += 'RESTART = \'inactive\'\n'
    else:
        st += 'RESTART = \'active\'\n'
    f.write(st)

    st = '\n'
    st += '# -- Define iterations, time (only usefull for unsteady computations) --\n'
    st += 'INITITER = '+VARS[2].get()+'\n'
    st += 'NITER = '+VARS[3].get()+'\n'
    st += 'INITIME = 0.\n'
    st += 'FINTIME = 100000.\n'
    st += 'TIMESTEP = 0.000001\n'
    f.write(st)

    st = '\n'
    st += '# -- Define flow conditions --\n'
    st += 'Mach = '+str(Mach)+'\n'
    st += 'adim = K.adimStateForStandardGas(MachInf=Mach, alphaZ='+str(alphaZ)+', alphaY='+str(alphaY)+', ReInf='+str(Reynolds)+')\n'
    if (dim == 2):
        st += 'Pb = E.cfdpb(\'CFD_Pb_2D\')\n'
        st += 'Pb.set(\'config\', \'2d\')\n'
    else:
        st += 'Pb = E.cfdpb(\'CFD_Pb_3D\')\n'
        st += 'Pb.set(\'config\', \'3d\')\n'
    f.write(st)

    st = '\n'
    st += '# -- Define fluid properties --\n'
    st += 'mod = K.createModel(MODEL, adim)\n'
    f.write(st)

    st = '\n'
    st += '# -- Define numerics --\n'
    st += 'num = K.createNumerics(INTEG, SCHEME, NITER, RESTART, INITITER, TIMESTEP, INITIME, FINTIME, adim, MODEL)\n'
    if (VARS[8].get() == 'Cartesian mesh generation'):
        st += 'num.set(\'automesh\', \'active\')\n'
        st += 'num.set(\'dfar\', '+VARS[9].get()+')\n'
    else:
        st += 'num.set(\'automesh\', \'inactive\')\n'
    f.write(st)

    st = '\n'
    st += 't = C.convertFile2PyTree(FILE, skeleton=1)\n'
    #st += 'if (size == 1):\n'
    #st += '    t = C.convertFile2PyTree(FILE)\n'
    #st += 'else:\n'
    #st += '    t = C.convertFile2PyTree(FILE, skeleton=1)\n'
    #st += '    t = Cmpi.es(t, FILE, rank=rank)\n'
    st += 'stateInf = K.createAdimState(adim, MODEL)\n'
    if (VARS[5].get() == 'Solver performs hole cutting'): delta = '-1'
    else: delta = VARS[6].get()
    st += 'K.setup(t, mod, num, Pb, stateInf, delta='+delta+')\n'
    st += 't = [] # to free memory\n'
    st += 'Pb.compute()\n'
    st += 'Pb.extract()\n'
    #st += '\n'
    #st += 'import Post.PyTree as P\n'
    #st += 't1 = C.convertFile2PyTree(\'output.plt\')\n'
    #st += 't = P.importVariables(t1, t)\n'
    #st += 'C.convertPyTree2File(t, \'output.cgns\')\n'
    #st += 't = C.center2Node(t, \'FlowSolution#Centers\')\n'
    #st += 't = C.rmVars(t, \'FlowSolution#Centers\')\n'
    #st += 'C.convertPyTree2File(t, \'outputn.cgns\')\n'
    f.write(st)
    CTK.TXT.insert('START', 'File setup.py written.\n')
    f.close()

#==============================================================================
# Run
#==============================================================================
def run():
    if (CTK.t == []): return
    treeFile = os.path.splitext(CTK.FILE)[0]+'.cgns'
    outputFile = 'xterm'
    C.convertPyTree2File(CTK.t, treeFile)
    writeSetupFile()
    if (outputFile == 'xterm'):
        os.system('python dump.py')
        os.system('xterm -e \'python setup.py; sleep 1000000\'&')
    else:
        os.system('python dump.py')
        os.system('python setup.py\ >& '+outputFile+'\'&')
    return

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkCassiopeeSolver', font=CTK.FRAMEFONT, 
                           takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Export to Cassiopee \nCartesian solver.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=0)
    Frame.columnconfigure(1, weight=1)
    WIDGETS['frame'] = Frame
    
    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkCassiopeeSolver')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Integ -
    V = TK.StringVar(win); V.set('steady'); VARS.append(V)
    # -1- Scheme -
    V = TK.StringVar(win); V.set('jameson'); VARS.append(V)
    # -2- Inititer
    V = TK.StringVar(win); V.set('1'); VARS.append(V)
    # -3- Niter
    V = TK.StringVar(win); V.set('10');VARS.append(V)
    # -4- Chimera option
    V = TK.StringVar(win)
    V.set('Solver performs hole cutting'); VARS.append(V)
    # -5- double wall tolerance -
    V = TK.StringVar(win); V.set('10.'); VARS.append(V)
    # -6- delta XRay -
    V = TK.StringVar(win); V.set('1.e-10'); VARS.append(V)    
    # -7- tolerance XRay -
    V = TK.StringVar(win); V.set('1.e-8'); VARS.append(V)
    # -8- Cartesian option -
    V = TK.StringVar(win); V.set('No Cartesian generation'); VARS.append(V)
    # -9- Dfar -
    V = TK.StringVar(win); V.set('10.'); VARS.append(V)
    
    # - Integ -
    B = TTK.Label(Frame, text="Integ")
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Kind of integration.')
    B = TTK.OptionMenu(Frame, VARS[0], 'steady', 'unsteady', 'gear')
    B.grid(row=0, column=1, sticky=TK.EW)

    # - Scheme -
    B = TTK.Label(Frame, text="Scheme")
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Numerical scheme used.')
    B = TTK.OptionMenu(Frame, VARS[1], 'jameson', 'jameson_ho',
                       'jameson_ho5', 'roe', 'ausmppmiles', 'ausmpup')
    B.grid(row=1, column=1, sticky=TK.EW)

    # - Inititer -
    B = TTK.Label(Frame, text="Inititer")
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='First iteration of this computation.')
    B = TTK.Entry(Frame, textvariable=VARS[2])
    B.grid(row=2, column=1, sticky=TK.EW)

    # - Niter -
    B = TTK.Label(Frame, text="Niter")
    B.grid(row=3, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='Number of iterations of this computation.')
    B = TTK.Entry(Frame, textvariable=VARS[3])
    B.grid(row=3, column=1, sticky=TK.EW)

    # - Chimera settings -
    B = TTK.OptionMenu(Frame, VARS[4],
                       'Solver performs hole cutting', 'Use OversetHole nodes')
    B.grid(row=4, column=0, columnspan=2, sticky=TK.EW)

    #- XRay settings  -
    B = TTK.Label(Frame, text="XRay delta")
    B.grid(row=5, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[6])
    B.grid(row=5, column=1, sticky=TK.EW)

    # - Cartesian option -
    B = TTK.OptionMenu(Frame, VARS[8],
                       'Cartesian mesh generation', 'No Cartesian generation')
    B.grid(row=6, column=0, columnspan=2, sticky=TK.EW)

    #- Dfar settings  -
    B = TTK.Label(Frame, text="Dfar")
    B.grid(row=7, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[9])
    B.grid(row=7, column=1, sticky=TK.EW)

    # - Run -
    B = TTK.Button(Frame, text="Run", command=run)
    B.grid(row=8, column=0, sticky=TK.EW)

    # - Write setup file -
    B = TTK.Button(Frame, text="Write setup", command=writeSetupFile)
    B.grid(row=8, column=1, sticky=TK.EW)
    
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
    (win, menu, file, tools) = CTK.minimal('tkCassiopeeSolver '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
