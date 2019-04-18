# - Flow Equation and Reference State manipulator -
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# get state from pyTree, l'affiche dans les widgets.
# Dans la norme CGNS, le state est stocke par Base ou par zone.
# Ici: si la selection est une base complete : on stocke dans la base
# Sinon, dans les zones.
#==============================================================================
def getState():
    if CTK.t == []: return

    import math
    referenceState = []
    flowEquationSet = []

    if CTK.__MAINTREE__ == 1:
        nzs = CPlot.getSelectedZones()
        fullBase = CPlot.isSelAFullBase(CTK.t, CTK.Nb, nzs)
    else: nzs = []; fullBase = 0

    if CTK.__MAINTREE__ <= 0 or nzs == []:
        A = Internal.getBases(CTK.t)
        if A != []: A = A[0]
        else: return
        for n in A[2]:
            if n[0] == 'FlowEquationSet': flowEquationSet = n; break
        for n in A[2]:
            if n[0] == 'ReferenceState': referenceState = n; break
    elif fullBase > 0:
        A = CTK.t[2][fullBase]
        for n in A[2]:
            if n[0] == 'FlowEquationSet': flowEquationSet = n; break
        for n in A[2]:
            if n[0] == 'ReferenceState': referenceState = n; break
    else:
         nob = CTK.Nb[nzs[0]]+1
         noz = CTK.Nz[nzs[0]]
         A = CTK.t[2][nob][2][noz]
         for n in A[2]:
             if n[0] == 'FlowEquationSet': flowEquationSet = n; break
         if flowEquationSet == []:
             B, r = Internal.getParentOfNode(CTK.t, A)
             for n in B[2]:
                 if (n[0] == 'FlowEquationSet'): flowEquationSet = n; break
         for n in A[2]:
             if n[0] == 'ReferenceState': referenceState = n; break
         if referenceState == []:
             B, r = Internal.getParentOfNode(CTK.t, A)
             for n in B[2]:
                 if (n[0] == 'ReferenceState'): referenceState = n; break

    if flowEquationSet != []:
        # EquationDimension
        node = Internal.getNodeFromName1(flowEquationSet, 'EquationDimension')
        if node is not None:
            dim = Internal.getValue(node)
            if (dim == 2): VARS[0].set('2D')
            elif (dim == 3): VARS[0].set('3D')

        # GoverningEquations
        node = Internal.getNodeFromName1(flowEquationSet, 'GoverningEquations')
        if node is not None:
            eq = Internal.getValue(node)
            VARS[1].set(eq)

    # ReferenceState
    if referenceState == []: return
    else: state = referenceState

    # Mach
    node = Internal.getNodeFromName1(state, 'Mach')
    if node is not None:
        mach = Internal.getValue(node)
        VARS[2].set(str(mach))
    else: mach = 0.5

    # Reynolds
    node = Internal.getNodeFromName1(state, 'Reynolds')
    if node is not None:
        reynolds = Internal.getValue(node) 
        VARS[3].set(str(reynolds))
    else: reynolds = 1.e6

    # Incidences
    node = Internal.getNodeFromName1(state, 'VelocityX')
    if node is not None: UInf = Internal.getValue(node)
    else: UInf = 0.
    node = Internal.getNodeFromName1(state, 'VelocityY')
    if node is not None: VInf = Internal.getValue(node)
    else: VInf = 0.
    node = Internal.getNodeFromName1(state, 'VelocityZ')
    if node is not None: WInf = Internal.getValue(node)
    else: WInf = 0.
    if abs(UInf) > 1.e-12:
        aly = math.atan(WInf/UInf)
        alz = math.atan(math.cos(aly)*VInf/UInf)
    else:
        aly = 0.; alz = 0.
    VARS[4].set(str(alz*180./math.pi))
    VARS[5].set(str(aly*180./math.pi))

    # Modele de turbulence
    node = Internal.getNodeFromName(state, 'Density')
    if node is not None: Density = Internal.getValue(node)
    else: Density = 1.
    node = Internal.getNodeFromName(state, 'Rok')
    if node is not None: RokInf = Internal.getValue(node)
    else: RokInf = 1.e-6
    Vit = UInf*UInf+VInf*VInf+WInf*WInf
    Vit = math.sqrt(Vit)
    if Vit > 1.e-10: TurbLevel = math.sqrt(2*RokInf/(3*Vit*Vit*Density))
    else: TurbLevel = 1.e-4
    MuInf = Density*Vit / max(reynolds,1.e-10) # L=1
    node = Internal.getNodeFromName(state, 'TurbulentSANuTildeDensity')
    if node is not None: RoNuTilde = Internal.getValue(node)
    else: RoNuTilde = 1.e-6
    MutInf = RoNuTilde
    MutSMu = MutInf / max(MuInf, 1.e-12)
    VARS[9].set(str(MutSMu))
    VARS[10].set(str(TurbLevel))

    CTK.TXT.insert('START', 'State displayed.\n')

#==============================================================================
# set dim in pyTree (all bases)
#==============================================================================
def setDim(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    # Pb dim
    if VARS[0].get() == '2D': dim = 2
    else: dim = 3
    CTK.saveTree()

    nodes = Internal.getBases(CTK.t)
    for b in nodes:
        C.addState2Node__(b, 'EquationDimension', dim)

    CTK.TKTREE.updateApp()
    CTK.TXT.insert('START', 'Dim set in all bases.\n')

#==============================================================================
# set state in pyTree
#==============================================================================
def setState(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()

    # Mach
    mach = VARS[2].get()
    try: mach = float(mach)
    except: mach = 0.
    
    # Reynolds
    Re = VARS[3].get()
    try: Re = float(Re)
    except: Re = 1.e6
    
    # Incidences
    alphaZ = VARS[4].get()
    try: alphaZ = float(alphaZ)
    except: alphaZ = 0.
    alphaY = VARS[5].get()
    try: alphaY = float(alphaY)
    except: alphaY = 0.

    # Grandeurs turb
    MutSMuInf = VARS[9].get()
    try: MutSMuInf = float(MutSMuInf)
    except: MutSMuInf = 0.2
    TurbLevelInf = VARS[10].get()
    try: TurbLevelInf = float(TurbLevelInf)
    except: TurbLevelInf = 1.e-4
    
    adim = ''; ADIM = VARS[11].get()
    if ADIM == 'adim1(Ro,A,T)': adim = 'adim1'
    elif ADIM == 'adim2(Ro,U,T)': adim = 'adim2'
    elif ADIM == 'dim1(real UInf,TInf,PInf,Rgp=287.053)': adim = 'dim1'
    elif ADIM == 'dim2(real UInf,TInf,RoInf,Rgp=287.053)': adim = 'dim2'
    elif ADIM == 'dim3(real UInf,PInf,RoInf,Rgp=287.053)': adim = 'dim3'
    CTK.saveTree()

    if CTK.__MAINTREE__ <= 0 or nzs == []:
        nodes = Internal.getBases(CTK.t)
        fullBase = False
    else:
        fullBase = CPlot.isSelAFullBase(CTK.t, CTK.Nb, nzs)
        if fullBase > 0:
            nodes = [CTK.t[2][fullBase]]
        else:
            nodes = []
            for nz in nzs:
                nob = CTK.Nb[nz]+1
                noz = CTK.Nz[nz]
                nodes.append(CTK.t[2][nob][2][noz])
    
    for b in nodes:
        p, r = Internal.getParentOfNode(CTK.t, b)
        C.addState2Node__(p[2][r], 'GoverningEquations', VARS[1].get())
        if VARS[1].get() == 'NSTurbulent':
            if VARS[6].get() == 'SpalartAllmaras':
                C.addState2Node__(p[2][r],
                                  'TurbulenceModel',
                                  'OneEquation_SpalartAllmaras')
            elif VARS[6].get() == 'JonesLaunder(k-eps)':
                C.addState2Node__(p[2][r],
                                  'TurbulenceModel',
                                  'TwoEquation_JonesLaunder')
            elif VARS[6].get() == 'Wilcox(k-w)':
                C.addState2Node__(p[2][r],
                                  'TurbulenceModel',
                                  'TwoEquation_Wilcox')
            elif VARS[6].get() == 'MenterSST(k-w)':
                C.addState2Node__(p[2][r],
                                  'TurbulenceModel',
                                  'TwoEquation_MenterSST')
        C._addState(p[2][r], MInf=mach, alphaZ=alphaZ, alphaY=alphaY, 
                    ReInf=Re, MutSMuInf=MutSMuInf, TurbLevelInf=TurbLevelInf,
                    adim=adim)

    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    if nzs == []: CTK.TXT.insert('START', 'State set in all bases.\n')
    elif fullBase > 0: CTK.TXT.insert('START', 'State set in selected base.\n')
    else: CTK.TXT.insert('START', 'State set in selected zones.\n')
    
#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkState', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Set information on case.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    WIDGETS['frame'] = Frame
    
    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkState')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Pb dim -
    V = TK.StringVar(win); V.set('3D'); VARS.append(V)
    # -1- Model -
    V = TK.StringVar(win); V.set('Euler'); VARS.append(V)
    # -2- Mach -
    V = TK.StringVar(win); V.set('0.5'); VARS.append(V)
    # -3- Reinf -
    V = TK.StringVar(win); V.set('1.e8'); VARS.append(V)
    # -4- Incidence/Z -
    V = TK.StringVar(win); V.set('0.'); VARS.append(V)
    # -5- Incidence/Y -
    V = TK.StringVar(win); V.set('0.'); VARS.append(V)
    # -6- Type de modele de turbulence
    V = TK.StringVar(win); V.set('SpalartAllmaras')
    VARS.append(V)
    # -7- Valeurs pour les grandeurs turbulentes
    V = TK.StringVar(win); V.set('1.e-6'); VARS.append(V)
    # -8- Jeux de variables definissant l'etat de reference
    V = TK.StringVar(win); V.set('[adim1] ro,T,a'); VARS.append(V)
    # -9- MutSMuInf
    V = TK.StringVar(win); V.set('0.2'); VARS.append(V)
    # -10- TurbLevelInf
    V = TK.StringVar(win); V.set('1.e-4'); VARS.append(V)
    # -11- Adim
    V = TK.StringVar(win); V.set('adim1(Ro,A,T)'); VARS.append(V)

    # - Pb dim -
    F = TTK.Frame(Frame, borderwidth=2, relief=CTK.FRAMESTYLE)
    F.columnconfigure(0, weight=1)
    F.columnconfigure(1, weight=2)
    F.columnconfigure(1, weight=2)
    B = TTK.Label(F, text="Pb dim")
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Dimension of the problem.')
    B = TTK.OptionMenu(F, VARS[0], '3D', '2D', command=setDim)
    B.grid(row=0, column=1, sticky=TK.EW)
    B = TTK.Button(F, text="Set dim", command=setDim)
    BB = CTK.infoBulle(parent=B, text='Set dim in all bases.')
    B.grid(row=0, column=2, sticky=TK.EW)
    F.grid(row=0, column=0, sticky=TK.EW)

    F = TTK.Frame(Frame, borderwidth=2, relief=CTK.FRAMESTYLE)
    F.columnconfigure(0, weight=1)
    F.columnconfigure(1, weight=4)

    # - Modele de fluide -
    B = TTK.Label(F, text="GovEquations")
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Governing equations to solve.')
    B = TTK.OptionMenu(F, VARS[1], 'Euler', 'NSLaminar', 'NSTurbulent')
    B.grid(row=1, column=1, sticky=TK.EW)

    # - Jeu de grandeurs definissant l'etat de reference -
    #B = TK.Label(F, text="Variables")
    #B.grid(row=2, column=0, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Change the set of variables defining the reference state.')
    #B = TK.OptionMenu(F, VARS[8], '[adim1] ro,a,T', '[adim2] ro,u,T', '[dim1] U,T,P,L', '[dim2] ro,u,T,L', '[dim3] ro,u,P,L')
    #B.grid(row=2, column=1, sticky=TK.EW)

    # - Mach -
    B = TTK.Label(F, text="Mach")
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Reference (infinite) mach number.')
    B = TTK.Entry(F, textvariable=VARS[2], background='White')
    B.grid(row=2, column=1, sticky=TK.EW)

    # - Reynolds -
    B = TTK.Label(F, text="Reynolds")
    B.grid(row=3, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Reference (infinite) Reynolds number.')
    B = TTK.Entry(F, textvariable=VARS[3], background='White')
    B.grid(row=3, column=1, sticky=TK.EW)

    # - Incidence/Z -
    B = TTK.Label(F, text="Incidence/Z (deg)")
    BB = CTK.infoBulle(parent=B, text='Angle of incident flow around z axis.')
    B.grid(row=4, column=0, sticky=TK.EW)
    B = TTK.Entry(F, textvariable=VARS[4], background='White')
    B.grid(row=4, column=1, sticky=TK.EW)

    # - Incidence/Y -
    B = TTK.Label(F, text="Incidence/Y (deg)")
    B.grid(row=5, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Angle of incident flow around y axis.')
    B = TTK.Entry(F, textvariable=VARS[5], background='White')
    B.grid(row=5, column=1, sticky=TK.EW)

    # - Modele de turbulence -
    B = TTK.Label(F, text="TubulenceModel")
    B.grid(row=6, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Type of turbulence model.')
    B = TTK.OptionMenu(F, VARS[6], 'SpalartAllmaras', 'JonesLaunder(k-eps)',
                      'Wilcox(k-w)', 'MenterSST(k-w)')
    B.grid(row=6, column=1, sticky=TK.EW)

    # - Valeurs des grandeurs turbulentes -
    B = TTK.Label(F, text="MutSMu")
    B.grid(row=7, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Ratio between turbulent viscosity and molecular viscosity.')
    B = TTK.Entry(F, textvariable=VARS[9], background='White')
    B.grid(row=7, column=1, sticky=TK.EW)

    B = TTK.Label(F, text="TurbLevel")
    B.grid(row=8, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Level of turbulence.')
    B = TTK.Entry(F, textvariable=VARS[10], background='White')
    B.grid(row=8, column=1, sticky=TK.EW)
    
    # - Adim -
    B = TTK.Label(F, text="Adim")
    B.grid(row=9, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Type of adimensionalization.')
    B = TTK.OptionMenu(F, VARS[11], 'adim1(Ro,A,T)','adim2(Ro,U,T)')#,'dim1(real UInf,TInf,PInf,Rgp=287.053)', 'dim2(real UInf,TInf,RoInf,Rgp=287.053)','dim3(real UInf,PInf,RoInf,Rgp=287.053)')
    B.grid(row=9, column=1, sticky=TK.EW)
    
    # - get state, inutile a mon avis -
    B = TTK.Button(F, text="Get state", command=getState)
    B.grid(row=10, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Read state from pyTree.')

    # - set state -
    B = TTK.Button(F, text="Set state", command=setState)
    B.grid(row=10, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Write state to pyTree.')

    F.grid(row=1, column=0, sticky=TK.EW)

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    WIDGETS['frame'].grid(sticky=TK.EW); getState()

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    WIDGETS['frame'].grid_forget()

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(event=None):
    getState()

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
    (win, menu, file, tools) = CTK.minimal('tkState '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
