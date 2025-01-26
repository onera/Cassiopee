# - tkState -
"""Flow Equation and Reference State manipulator."""
try: import tkinter as TK
except: import Tkinter as TK
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
                if n[0] == 'ReferenceState': referenceState = n; break

    if flowEquationSet != []:
        # EquationDimension
        node = Internal.getNodeFromName1(flowEquationSet, 'EquationDimension')
        if node is not None:
            dim = Internal.getValue(node)
            if dim == 2: VARS[0].set('2D')
            elif dim == 3: VARS[0].set('3D')

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

    # Reynolds
    node = Internal.getNodeFromName1(state, 'Reynolds')
    if node is not None:
        reynolds = Internal.getValue(node)
        VARS[3].set(str(reynolds))
    else: reynolds = None

    # UInf
    node = Internal.getNodeFromName1(state, 'VelocityX')
    if node is not None: UInf = Internal.getValue(node)
    else: UInf = None
    node = Internal.getNodeFromName1(state, 'VelocityY')
    if node is not None: VInf = Internal.getValue(node)
    else: VInf = None
    node = Internal.getNodeFromName1(state, 'VelocityZ')
    if node is not None: WInf = Internal.getValue(node)
    else: WInf = None
    if UInf is not None and VInf is not None and WInf is not None:
        Vit = UInf*UInf+VInf*VInf+WInf*WInf
        Vit = math.sqrt(Vit)
        VARS[12].set(str(Vit))
    else: Vit = None

    # TInf
    node = Internal.getNodeFromName1(state, 'Temperature')
    if node is not None:
        TInf = Internal.getValue(node)
        VARS[13].set(str(TInf))

    # PInf
    node = Internal.getNodeFromName1(state, 'Pressure')
    if node is not None:
        PInf = Internal.getValue(node)
        VARS[14].set(str(PInf))

    # RoInf
    node = Internal.getNodeFromName1(state, 'Density')
    if node is not None:
        RoInf = Internal.getValue(node)
        VARS[16].set(str(RoInf))

    # Incidences
    if UInf is not None and VInf is not None and WInf is not None and abs(UInf) > 1.e-12:
        aly = math.atan(WInf/UInf)
        alz = math.atan(math.cos(aly)*VInf/UInf)
        VARS[4].set(str(alz*180./math.pi))
        VARS[5].set(str(aly*180./math.pi))

    # Modele de turbulence
    node = Internal.getNodeFromName(state, 'Density')
    if node is not None: Density = Internal.getValue(node)
    else: Density = None
    node = Internal.getNodeFromName(state, 'Rok')
    if node is not None: RokInf = Internal.getValue(node)
    else: RokInf = None

    if (reynolds is not None and Density is not None and
            RokInf is not None and Vit is not None and Vit > 1.e-10):
        TurbLevel = math.sqrt(2*RokInf/(3*Vit*Vit*Density))
        VARS[10].set(str(TurbLevel))

        MuInf = Density*Vit / max(reynolds,1.e-10) # L=1

        node = Internal.getNodeFromName(state, 'TurbulentSANuTildeDensity')
        if node is not None:
            RoNuTilde = Internal.getValue(node)
            MutInf = RoNuTilde
            MutSMu = MutInf / max(MuInf, 1.e-12)
            VARS[9].set(str(MutSMu))

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

    # Grandeurs primaires
    RoInf = VARS[16].get()
    try: RoInf = float(RoInf)
    except: RoInf = 1.
    UInf = VARS[12].get()
    try: UInf = float(UInf)
    except: UInf = 1.
    TInf = VARS[13].get()
    try: TInf = float(TInf)
    except: TInf = 298.
    PInf = VARS[14].get()
    try: PInf = float(PInf)
    except: PInf = 101325.
    LInf = VARS[15].get()
    try: LInf = float(LInf)
    except: LInf = 1.

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
    elif ADIM == 'dim1(real UInf,TInf,PInf)': adim = 'dim1'
    elif ADIM == 'dim2(real UInf,TInf,RoInf)': adim = 'dim2'
    elif ADIM == 'dim3(real UInf,PInf,RoInf)': adim = 'dim3'
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
                    UInf=UInf, TInf=TInf, LInf=LInf, RoInf=RoInf, PInf=PInf,
                    adim=adim)

    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    if nzs == []: CTK.TXT.insert('START', 'State set in all bases.\n')
    elif fullBase > 0: CTK.TXT.insert('START', 'State set in selected base.\n')
    else: CTK.TXT.insert('START', 'State set in selected zones.\n')

# Called when adim is switched -> change buttons
def switchAdim(event=None):
    adim = VARS[11].get()
    if adim == 'adim1(Ro,A,T)' or adim == 'adim2(Ro,U,T)':
        # switch UInf and Mach
        for w in WIDGETS['UInf']: w.grid_forget()
        for c, w in enumerate(WIDGETS['Mach']):
            w.grid(row=2, column=c, sticky=TK.EW)
        # switch TInf or PInf to Reynolds
        for w in WIDGETS['TInf']: w.grid_forget()
        for w in WIDGETS['PInf']: w.grid_forget()
        for w in WIDGETS['RoInf']: w.grid_forget()
        for w in WIDGETS['LInf']: w.grid_forget()
        for c, w in enumerate(WIDGETS['Reynolds']):
            w.grid(row=4, column=c, sticky=TK.EW)
    elif adim == 'dim1(real UInf,TInf,PInf)':
        # switch Mach and UInf
        for w in WIDGETS['Mach']: w.grid_forget()
        for c, w in enumerate(WIDGETS['UInf']):
            w.grid(row=2, column=c, sticky=TK.EW)
        for w in WIDGETS['Reynolds']: w.grid_forget()
        for w in WIDGETS['RoInf']: w.grid_forget()
        for c, w in enumerate(WIDGETS['TInf']):
            w.grid(row=4, column=c, sticky=TK.EW)
        for c, w in enumerate(WIDGETS['PInf']):
            w.grid(row=4, column=c+2, sticky=TK.EW)
        for c, w in enumerate(WIDGETS['LInf']):
            w.grid(row=5, column=c, sticky=TK.EW)
    elif adim == 'dim2(real UInf,TInf,RoInf)':
        # switch Mach and UInf
        for w in WIDGETS['Mach']: w.grid_forget()
        for c, w in enumerate(WIDGETS['UInf']):
            w.grid(row=2, column=c, sticky=TK.EW)
        for w in WIDGETS['Reynolds']: w.grid_forget()
        for w in WIDGETS['PInf']: w.grid_forget()
        for c, w in enumerate(WIDGETS['TInf']):
            w.grid(row=4, column=c, sticky=TK.EW)
        for c, w in enumerate(WIDGETS['RoInf']):
            w.grid(row=4, column=c+2, sticky=TK.EW)
        for c, w in enumerate(WIDGETS['LInf']):
            w.grid(row=5, column=c, sticky=TK.EW)
    elif adim == 'dim3(real UInf,PInf,RoInf)':
        # switch Mach and UInf
        for w in WIDGETS['Mach']: w.grid_forget()
        for c, w in enumerate(WIDGETS['UInf']):
            w.grid(row=2, column=c, sticky=TK.EW)
        for w in WIDGETS['Reynolds']: w.grid_forget()
        for w in WIDGETS['TInf']: w.grid_forget()
        for c, w in enumerate(WIDGETS['PInf']):
            w.grid(row=4, column=c, sticky=TK.EW)
        for c, w in enumerate(WIDGETS['RoInf']):
            w.grid(row=4, column=c+2, sticky=TK.EW)
        for c, w in enumerate(WIDGETS['LInf']):
            w.grid(row=5, column=c, sticky=TK.EW)

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkState  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Set information on case.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkState')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Pb dim -
    V = TK.StringVar(win); V.set('3D'); VARS.append(V)
    # -1- Model -
    V = TK.StringVar(win); V.set('NSLaminar'); VARS.append(V)
    # -2- Mach -
    V = TK.StringVar(win); V.set('0.5'); VARS.append(V)
    # -3- Reinf -
    V = TK.StringVar(win); V.set('1.e8'); VARS.append(V)
    # -4- Incidence/Z -
    V = TK.StringVar(win); V.set('0.'); VARS.append(V)
    # -5- Incidence/Y -
    V = TK.StringVar(win); V.set('0.'); VARS.append(V)
    # -6- Type de modele de turbulence
    V = TK.StringVar(win); V.set('SpalartAllmaras'); VARS.append(V)
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
    # -12- UInf
    V = TK.StringVar(win); V.set('2.8'); VARS.append(V)
    # -13- TInf
    V = TK.StringVar(win); V.set('298.'); VARS.append(V)
    # -14- PInf
    V = TK.StringVar(win); V.set('101325.'); VARS.append(V)
    # -15- LInf
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)
    # -16- RoInf
    V = TK.StringVar(win); V.set('1.225'); VARS.append(V)
    # -17- MachTip
    V = TK.StringVar(win); V.set('None'); VARS.append(V)

    # - Pb dim -
    F = TTK.Frame(Frame, borderwidth=2, relief=CTK.FRAMESTYLE)
    F.columnconfigure(0, weight=1)
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
    F.columnconfigure(2, weight=1)
    F.columnconfigure(3, weight=4)

    # - Modele de fluide -
    B = TTK.Label(F, text="Equations")
    B.grid(row=1, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Governing equations to solve.')
    B = TTK.OptionMenu(F, VARS[1], 'Euler', 'NSLaminar', 'NSTurbulent')
    B.grid(row=1, column=1, columnspan=3, sticky=TK.EW)

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
    WIDGETS['Mach'] = [B]
    B = TTK.Entry(F, textvariable=VARS[2], background='White', width=6)
    B.grid(row=2, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Reference (infinite) mach number.')
    WIDGETS['Mach'] += [B]

    # - MachTip -
    B = TTK.Label(F, text="MachTip")
    B.grid(row=2, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='If this mach number is defined, it is used for turbulence definition.\nThis is the case of rotors.')
    WIDGETS['MachTip'] = [B]
    B = TTK.Entry(F, textvariable=VARS[17], background='White', width=6)
    B.grid(row=2, column=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='If this mach number is defined, it is used for turbulence definition.\nThis is the case of rotors.')
    WIDGETS['MachTip'] += [B]

    # - UInf -
    B = TTK.Label(F, text="UInf (m/s)")
    BB = CTK.infoBulle(parent=B, text='Reference (infinite) speed (m/s).')
    WIDGETS['UInf'] = [B]
    B = TTK.Entry(F, textvariable=VARS[12], background='White', width=6)
    BB = CTK.infoBulle(parent=B, text='Reference (infinite) speed (m/s).')
    WIDGETS['UInf'] += [B]

    # - Incidence/Z -
    B = TTK.Label(F, text="Inc./Z")
    BB = CTK.infoBulle(parent=B, text='Angle of incident flow around z axis (deg).')
    B.grid(row=3, column=0, sticky=TK.EW)
    WIDGETS['IncZ'] = [B]
    B = TTK.Entry(F, textvariable=VARS[4], background='White', width=6)
    B.grid(row=3, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Angle of incident flow around z axis (deg).')
    WIDGETS['IncZ'] += [B]

    # - Incidence/Y -
    B = TTK.Label(F, text="Inc./Y")
    B.grid(row=3, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Angle of incident flow around y axis (deg).')
    WIDGETS['IncY'] = [B]
    B = TTK.Entry(F, textvariable=VARS[5], background='White', width=6)
    B.grid(row=3, column=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Angle of incident flow around y axis (deg).')
    WIDGETS['IncY'] += [B]

    # - Reynolds -
    B = TTK.Label(F, text="Reynolds")
    B.grid(row=4, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Reference (infinite) Reynolds number.')
    WIDGETS['Reynolds'] = [B]
    B = TTK.Entry(F, textvariable=VARS[3], background='White', width=6)
    BB = CTK.infoBulle(parent=B, text='Reference (infinite) Reynolds number.')
    B.grid(row=4, column=1, sticky=TK.EW)
    WIDGETS['Reynolds'] += [B]

    # - TInf -
    B = TTK.Label(F, text="TInf (K)")
    BB = CTK.infoBulle(parent=B, text='Reference (infinite) temperature (K).')
    WIDGETS['TInf'] = [B]
    B = TTK.Entry(F, textvariable=VARS[13], background='White', width=6)
    BB = CTK.infoBulle(parent=B, text='Reference (infinite) temperature (K).')
    WIDGETS['TInf'] += [B]

    # - PInf -
    B = TTK.Label(F, text="PInf (Pa)")
    BB = CTK.infoBulle(parent=B, text='Reference (infinite) pressure (Pa).')
    WIDGETS['PInf'] = [B]
    B = TTK.Entry(F, textvariable=VARS[14], background='White', width=6)
    BB = CTK.infoBulle(parent=B, text='Reference (infinite) pressure (Pa).')
    WIDGETS['PInf'] += [B]

    # - LInf -
    B = TTK.Label(F, text="LInf (m)")
    BB = CTK.infoBulle(parent=B, text='Reference length (m).')
    WIDGETS['LInf'] = [B]
    B = TTK.Entry(F, textvariable=VARS[15], background='White', width=6)
    BB = CTK.infoBulle(parent=B, text='Reference length (m).')
    WIDGETS['LInf'] += [B]

    # - RoInf -
    B = TTK.Label(F, text="RoInf (kg/m3)")
    BB = CTK.infoBulle(parent=B, text='Reference (infinite) density (kg/m3).')
    WIDGETS['RoInf'] = [B]
    B = TTK.Entry(F, textvariable=VARS[16], background='White', width=6)
    BB = CTK.infoBulle(parent=B, text='Reference (infinite) density (kg/m3).')
    WIDGETS['RoInf'] += [B]

    # - Modele de turbulence -
    B = TTK.Label(F, text="TurbModel")
    B.grid(row=6, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Type of turbulence model.')
    B = TTK.OptionMenu(F, VARS[6], 'SpalartAllmaras', 'JonesLaunder(k-eps)',
                       'Wilcox(k-w)', 'MenterSST(k-w)')
    B.grid(row=6, column=1, columnspan=3, sticky=TK.EW)

    # - Valeurs des grandeurs turbulentes -
    B = TTK.Label(F, text="MutSMu")
    B.grid(row=7, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Ratio between turbulent viscosity and molecular viscosity.')
    B = TTK.Entry(F, textvariable=VARS[9], background='White', width=5)
    B.grid(row=7, column=1, sticky=TK.EW)

    B = TTK.Label(F, text="TurbLevel")
    B.grid(row=7, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Level of turbulence.')
    B = TTK.Entry(F, textvariable=VARS[10], background='White', width=5)
    B.grid(row=7, column=3, sticky=TK.EW)

    # - Adim -
    B = TTK.Label(F, text="Adim")
    B.grid(row=8, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Type of adimensionalization.')
    B = TTK.OptionMenu(F, VARS[11], 'adim1(Ro,A,T)','adim2(Ro,U,T)',
                       'dim1(real UInf,TInf,PInf)',
                       'dim2(real UInf,TInf,RoInf)',
                       'dim3(real UInf,PInf,RoInf)',
                       command=switchAdim)
    B.grid(row=8, column=1, columnspan=3, sticky=TK.EW)

    # - get state, inutile a mon avis -
    B = TTK.Button(F, text="Get state", command=getState)
    B.grid(row=9, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Read state from pyTree.')

    # - set state -
    B = TTK.Button(F, text="Set state", command=setState)
    B.grid(row=9, column=2, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Write state to pyTree.')

    F.grid(row=1, column=0, sticky=TK.EW)

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['StateNoteBook'].add(WIDGETS['frame'], text='tkState')
    except: pass
    CTK.WIDGETS['StateNoteBook'].select(WIDGETS['frame'])
    getState()

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['StateNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(event=None):
    getState()

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
    (win, menu, file, tools) = CTK.minimal('tkState '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
