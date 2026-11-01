# - rigid motion settings -
import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal
import RigidMotion.PyTree as RM

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# Retourne le no du motion (1,2,..) suivant VARS[11]
#==============================================================================
def getMotionType():
    ntype = VARS[11].get()
    if ntype == '1:MotionStrings': return 1
    elif ntype == '2:KMotionRotor': return 2
    elif ntype == '3:ConstantMotion': return 3

#==============================================================================
def setTransOrigin():
    point = CPlot.getActivePoint()
    if point != []:
        VARS[1].set(str(point[0]))
        VARS[2].set(str(point[1]))
        VARS[3].set(str(point[2]))

#==============================================================================
def setCenterRotation():
    point = CPlot.getActivePoint()
    if point != []:
        VARS[4].set(str(point[0]))
        VARS[5].set(str(point[1]))
        VARS[6].set(str(point[2]))

#==============================================================================
def set2Axis():
    point = CPlot.getActivePoint()
    if point != []:
        VARS[7].set(str(point[0]))
        VARS[8].set(str(point[1]))
        VARS[9].set(str(point[2]))

#==============================================================================
def resetVars():
    type = getMotionType()
    if type == 1: resetVars1()
    elif type == 2: resetVars2()
    elif type == 3: resetVars3()

#==============================================================================
def resetVars1():
    VARS[1].set('0')
    VARS[2].set('0')
    VARS[3].set('0')
    VARS[4].set('0')
    VARS[5].set('0')
    VARS[6].set('0')
    VARS[7].set('0')
    VARS[8].set('0')
    VARS[9].set('0')
    VARS[10].set('0')

#==============================================================================
def resetVars2():
    VARS[17].set('0.; 0.; 0.')
    VARS[12].set(0.)
    VARS[13].set(0.)

#==============================================================================
def resetVars3():
    VARS[43].set('0.; 0.; 0.')
    VARS[44].set('0.; 0.; 0.')
    VARS[45].set('0.; 0.; 1.')
    VARS[46].set(0.)

#==============================================================================
def getVars(event=None):
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nob = CTK.Nb[nzs[0]]+1
    noz = CTK.Nz[nzs[0]]
    z = CTK.t[2][nob][2][noz]
    cont = Internal.getNodesFromName1(z, 'TimeMotion')
    if cont == []:
        CTK.TXT.insert('START', 'No motion in this zone.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    name = VARS[0].get()
    motion = Internal.getNodesFromName1(cont[0], name)
    if motion == []:
        name = cont[0][2][0][0]
        motion = Internal.getNodesFromName1(cont[0], name)
        VARS[0].set(name)

    motion = motion[0]
    t = Internal.getNodesFromName1(motion, 'MotionType')
    ntype = t[0][1][0]

    if ntype == 1:
        VARS[11].set('1:MotionStrings'); changeMotionType(); getVars1()
    elif ntype == 2:
        VARS[11].set('2:KMotionRotor'); changeMotionType(); getVars2()
    elif ntype == 3:
        VARS[11].set('3:ConstantMotion'); changeMotionType(); getVars3()

#==============================================================================
# Si Name existe deja dans la selection, affiche
# Sinon, ne fait rien
#==============================================================================
def getVars1(event=None):
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    nob = CTK.Nb[nzs[0]]+1
    noz = CTK.Nz[nzs[0]]
    z = CTK.t[2][nob][2][noz]

    cont = Internal.getNodesFromName1(z, 'TimeMotion')
    if cont == []:
        CTK.TXT.insert('START', 'No motion in this zone.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    name = VARS[0].get()

    motion = Internal.getNodesFromName1(cont[0], name)
    if motion == []:
        CTK.TXT.insert('START', 'No motion named '+name+' in this zone.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    motion = motion[0]
    t = Internal.getNodesFromName1(motion, 'tx')
    VARS[1].set(t[0][1].tobytes().decode())
    t = Internal.getNodesFromName1(motion, 'ty')
    VARS[2].set(t[0][1].tobytes().decode())
    t = Internal.getNodesFromName1(motion, 'tz')
    VARS[3].set(t[0][1].tobytes().decode())
    t = Internal.getNodesFromName1(motion, 'cx')
    VARS[4].set(t[0][1].tobytes().decode())
    t = Internal.getNodesFromName1(motion, 'cy')
    VARS[5].set(t[0][1].tobytes().decode())
    t = Internal.getNodesFromName1(motion, 'cz')
    VARS[6].set(t[0][1].tobytes().decode())
    t = Internal.getNodesFromName1(motion, 'ex')
    VARS[7].set(t[0][1].tobytes().decode())
    t = Internal.getNodesFromName1(motion, 'ey')
    VARS[8].set(t[0][1].tobytes().decode())
    t = Internal.getNodesFromName1(motion, 'ez')
    VARS[9].set(t[0][1].tobytes().decode())
    t = Internal.getNodesFromName1(motion, 'angle')
    VARS[10].set(t[0][1].tobytes().decode())

#==============================================================================
# Si Name existe deja dans la selection, affiche
# Sinon, ne fait rien
#==============================================================================
def getVars2(event=None):
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    nob = CTK.Nb[nzs[0]]+1
    noz = CTK.Nz[nzs[0]]
    z = CTK.t[2][nob][2][noz]

    cont = Internal.getNodesFromName1(z, 'TimeMotion')
    if cont == []:
        CTK.TXT.insert('START', 'No motion in this zone.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    name = VARS[0].get()
    motion = Internal.getNodesFromName1(cont[0], name)
    if motion == []:
        CTK.TXT.insert('START', 'No motion named '+name+' in this zone.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    motion = motion[0]
    t = Internal.getNodesFromName1(motion, 'transl_speed')
    VARS[17].set(str(t[0][1][0])+'; '+str(t[0][1][1])+'; '+str(t[0][1][2]))
    t = Internal.getNodesFromName1(motion, 'psi0')
    VARS[12].set(t[0][1][0])
    t = Internal.getNodesFromName1(motion, 'psi0_b')
    VARS[13].set(t[0][1][0])
    t = Internal.getNodesFromName1(motion, 'alp_pnt')
    VARS[14].set(str(t[0][1][0])+'; '+str(t[0][1][1])+'; '+str(t[0][1][2]))
    t = Internal.getNodesFromName1(motion, 'alp_vct')
    VARS[15].set(str(t[0][1][0])+'; '+str(t[0][1][1])+'; '+str(t[0][1][2]))
    t = Internal.getNodesFromName1(motion, 'alp0')
    VARS[16].set(t[0][1][0])
    t = Internal.getNodesFromName1(motion, 'rot_pnt')
    VARS[18].set(str(t[0][1][0])+'; '+str(t[0][1][1])+'; '+str(t[0][1][2]))
    t = Internal.getNodesFromName1(motion, 'rot_vct')
    VARS[19].set(str(t[0][1][0])+'; '+str(t[0][1][1])+'; '+str(t[0][1][2]))
    t = Internal.getNodesFromName1(motion, 'rot_omg')
    VARS[20].set(t[0][1][0])

    t = Internal.getNodesFromName1(motion, 'del_pnt')
    VARS[21].set(str(t[0][1][0])+'; '+str(t[0][1][1])+'; '+str(t[0][1][2]))
    t = Internal.getNodesFromName1(motion, 'del_vct')
    VARS[22].set(str(t[0][1][0])+'; '+str(t[0][1][1])+'; '+str(t[0][1][2]))
    t = Internal.getNodesFromName1(motion, 'del0')
    VARS[23].set(t[0][1][0])
    t = Internal.getNodesFromName1(motion, 'delc')
    VARS[24].set(str(t[0][1][0])+'; '+str(t[0][1][1])+'; '+str(t[0][1][2]))
    t = Internal.getNodesFromName1(motion, 'dels')
    VARS[25].set(str(t[0][1][0])+'; '+str(t[0][1][1])+'; '+str(t[0][1][2]))

    t = Internal.getNodesFromName1(motion, 'bet_pnt')
    VARS[26].set(str(t[0][1][0])+'; '+str(t[0][1][1])+'; '+str(t[0][1][2]))
    t = Internal.getNodesFromName1(motion, 'bet_vct')
    VARS[27].set(str(t[0][1][0])+'; '+str(t[0][1][1])+'; '+str(t[0][1][2]))
    t = Internal.getNodesFromName1(motion, 'bet0')
    VARS[28].set(t[0][1][0])
    t = Internal.getNodesFromName1(motion, 'betc')
    VARS[29].set(str(t[0][1][0])+'; '+str(t[0][1][1])+'; '+str(t[0][1][2]))
    t = Internal.getNodesFromName1(motion, 'bets')
    VARS[30].set(str(t[0][1][0])+'; '+str(t[0][1][1])+'; '+str(t[0][1][2]))

    t = Internal.getNodesFromName1(motion, 'tet_pnt')
    VARS[31].set(str(t[0][1][0])+'; '+str(t[0][1][1])+'; '+str(t[0][1][2]))
    t = Internal.getNodesFromName1(motion, 'tet_vct')
    VARS[32].set(str(t[0][1][0])+'; '+str(t[0][1][1])+'; '+str(t[0][1][2]))
    t = Internal.getNodesFromName1(motion, 'tet0')
    VARS[33].set(t[0][1][0])
    t = Internal.getNodesFromName1(motion, 'tetc')
    VARS[34].set(str(t[0][1][0]))
    t = Internal.getNodesFromName1(motion, 'tets')
    VARS[35].set(str(t[0][1][0]))

    t = Internal.getNodesFromName1(motion, 'span_vct')
    VARS[36].set(str(t[0][1][0])+'; '+str(t[0][1][1])+'; '+str(t[0][1][2]))

    t = Internal.getNodesFromName1(motion, 'pre_lag_pnt')
    VARS[37].set(str(t[0][1][0])+'; '+str(t[0][1][1])+'; '+str(t[0][1][2]))
    t = Internal.getNodesFromName1(motion, 'pre_lag_vct')
    VARS[38].set(str(t[0][1][0])+'; '+str(t[0][1][1])+'; '+str(t[0][1][2]))
    t = Internal.getNodesFromName1(motion, 'pre_lag_ang')
    VARS[39].set(t[0][1][0])

    t = Internal.getNodesFromName1(motion, 'pre_con_pnt')
    VARS[40].set(str(t[0][1][0])+'; '+str(t[0][1][1])+'; '+str(t[0][1][2]))
    t = Internal.getNodesFromName1(motion, 'pre_con_vct')
    VARS[41].set(str(t[0][1][0])+'; '+str(t[0][1][1])+'; '+str(t[0][1][2]))
    t = Internal.getNodesFromName1(motion, 'pre_con_ang')
    VARS[42].set(t[0][1][0])

#==============================================================================
# Si Name existe deja dans la selection, affiche
# Sinon, ne fait rien
#==============================================================================
def getVars3(event=None):
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    nob = CTK.Nb[nzs[0]]+1
    noz = CTK.Nz[nzs[0]]
    z = CTK.t[2][nob][2][noz]

    cont = Internal.getNodesFromName1(z, 'TimeMotion')
    if cont == []:
        CTK.TXT.insert('START', 'No motion in this zone.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    name = VARS[0].get()

    motion = Internal.getNodesFromName1(cont[0], name)
    if motion == []:
        CTK.TXT.insert('START', 'No motion named '+name+' in this zone.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    motion = motion[0]
    t = Internal.getNodesFromName1(motion, 'transl_speed')
    VARS[43].set(str(t[0][1][0])+'; '+str(t[0][1][1])+'; '+str(t[0][1][2]))
    t = Internal.getNodesFromName1(motion, 'axis_pnt')
    VARS[44].set(str(t[0][1][0])+'; '+str(t[0][1][1])+'; '+str(t[0][1][2]))
    t = Internal.getNodesFromName1(motion, 'axis_vct')
    VARS[45].set(str(t[0][1][0])+'; '+str(t[0][1][1])+'; '+str(t[0][1][2]))
    t = Internal.getNodesFromName1(motion, 'omega')
    VARS[46].set(t[0][1][0])

#==============================================================================
def setVars(event=None):
    ntype = getMotionType()
    if ntype == 1: setVars1()
    elif ntype == 2: setVars2()
    elif ntype == 3: setVars3()

#==============================================================================
# set vars pour le motion type 1
#==============================================================================
def setVars1(event=None):
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    name = VARS[0].get()
    CTK.saveTree()

    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        CTK.t[2][nob][2][noz] = RM.setPrescribedMotion1(
            z, name,
            VARS[1].get(), VARS[2].get(),
            VARS[3].get(), VARS[4].get(),
            VARS[5].get(), VARS[6].get(),
            VARS[7].get(), VARS[8].get(),
            VARS[9].get(), VARS[10].get())
    CTK.TKTREE.updateApp()
    CTK.TXT.insert('START', 'Motion set in selected zones.\n')

#==============================================================================
# set vars pour le motion type 2
#==============================================================================
def setVars2(event=None):
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    name = VARS[0].get()
    CTK.saveTree()

    transl_speed = CTK.varsFromWidget(VARS[17].get(), 1)
    psi0 = VARS[12].get()
    psi0_b = VARS[13].get()
    alp_pnt = CTK.varsFromWidget(VARS[14].get(), 1)
    alp_vct = CTK.varsFromWidget(VARS[15].get(), 1)
    alp0 = VARS[16].get()
    rot_pnt = CTK.varsFromWidget(VARS[18].get(), 1)
    rot_vct = CTK.varsFromWidget(VARS[19].get(), 1)
    rot_omg = VARS[20].get()

    del_pnt = CTK.varsFromWidget(VARS[21].get(), 1)
    del_vct = CTK.varsFromWidget(VARS[22].get(), 1)
    del0 = VARS[23].get()
    delc = CTK.varsFromWidget(VARS[24].get(), 1)
    dels = CTK.varsFromWidget(VARS[25].get(), 1)

    bet_pnt = CTK.varsFromWidget(VARS[26].get(), 1)
    bet_vct = CTK.varsFromWidget(VARS[27].get(), 1)
    bet0 = VARS[28].get()
    betc = CTK.varsFromWidget(VARS[29].get(), 1)
    bets = CTK.varsFromWidget(VARS[30].get(), 1)

    tet_pnt = CTK.varsFromWidget(VARS[31].get(), 1)
    tet_vct = CTK.varsFromWidget(VARS[32].get(), 1)
    tet0 = VARS[33].get()
    tetc = CTK.varsFromWidget(VARS[34].get(), 1)
    tets = CTK.varsFromWidget(VARS[35].get(), 1)

    span_vct = CTK.varsFromWidget(VARS[36].get(), 1)

    pre_lag_pnt = CTK.varsFromWidget(VARS[37].get(), 1)
    pre_lag_vct = CTK.varsFromWidget(VARS[38].get(), 1)
    pre_lag_ang = VARS[39].get()

    pre_con_pnt = CTK.varsFromWidget(VARS[40].get(), 1)
    pre_con_vct = CTK.varsFromWidget(VARS[41].get(), 1)
    pre_con_ang = VARS[42].get()

    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        CTK.t[2][nob][2][noz] = RM.setPrescribedMotion2(
            z, name,
            transl_speed=transl_speed,
            psi0=psi0, psi0_b=psi0_b,
            alp_pnt=alp_pnt, alp_vct=alp_vct,
            alp0=alp0,
            rot_pnt=rot_pnt, rot_vct=rot_vct,
            rot_omg=rot_omg,
            del_pnt=del_pnt, del_vct=del_vct,
            del0=del0,
            delc=delc, dels=dels,
            bet_pnt=bet_pnt, bet_vct=bet_vct,
            bet0=bet0,
            betc=betc, bets=bets,
            tet_pnt=tet_pnt, tet_vct=tet_vct,
            tet0=tet0,
            tetc=tetc, tets=tets,
            span_vct=span_vct,
            pre_lag_pnt=pre_lag_pnt, pre_lag_vct=pre_lag_vct,
            pre_lag_ang=pre_lag_ang,
            pre_con_pnt=pre_con_pnt, pre_con_vct=pre_con_vct,
            pre_con_ang=pre_con_ang)
    CTK.TKTREE.updateApp()
    CTK.TXT.insert('START', 'Motion set in selected zones.\n')

#==============================================================================
# set vars pour le motion type 3
#==============================================================================
def setVars3(event=None):
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    name = VARS[0].get()
    CTK.saveTree()

    transl_speed = CTK.varsFromWidget(VARS[43].get(), 1)
    axis_pnt = CTK.varsFromWidget(VARS[44].get(), 1)
    axis_vct = CTK.varsFromWidget(VARS[45].get(), 1)
    omega = VARS[46].get()

    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        CTK.t[2][nob][2][noz] = RM.setPrescribedMotion3(
            z, name,
            transl_speed=transl_speed, axis_pnt=axis_pnt,
            axis_vct=axis_vct, omega=omega)
    CTK.TKTREE.updateApp()
    CTK.TXT.insert('START', 'Motion set in selected zones.\n')

#==============================================================================
def changeMotionType(event=None):
    ntype = getMotionType()
    WIDGETS['frame1'].grid_forget()
    WIDGETS['frame2'].grid_forget()
    WIDGETS['frame3'].grid_forget()
    if ntype == 1:
        WIDGETS['frame1'].grid(row=2, column=0, columnspan=3, sticky=TK.EW)
    elif ntype == 2:
        WIDGETS['frame2'].grid(row=2, column=0, columnspan=3, sticky=TK.EW)
    elif ntype == 3:
        WIDGETS['frame3'].grid(row=2, column=0, columnspan=3, sticky=TK.EW)
    return

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkRigidMotion  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Define block rigid motions.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=0)
    Frame.columnconfigure(1, weight=1)
    Frame.columnconfigure(2, weight=0)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkRigidMotion')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- motion name -
    V = TK.StringVar(win); V.set('New'); VARS.append(V)
    #---------------------
    # -- Motion strings --
    #---------------------
    # -1- tx -
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -2- ty -
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -3- tz -
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -4- cx -
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -5- cy -
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -6- cz -
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -7- ex -
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -8- ey -
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -9- ez -
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -10- angle -
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -11- motion type
    V = TK.StringVar(win); V.set('1:MotionStrings'); VARS.append(V)
    #-----------------------------
    # -- Cassiopee Motion rotor --
    #-----------------------------
    # -12- psi0 -
    V = TK.DoubleVar(win); V.set(0.); VARS.append(V)
    # -13- psi0_b -
    V = TK.DoubleVar(win); V.set(0.); VARS.append(V)
    # -14- alp_pnt -
    V = TK.StringVar(win); V.set('0.; 0.; 0.'); VARS.append(V)
    # -15- alp_vct -
    V = TK.StringVar(win); V.set('0.; 1.; 0.'); VARS.append(V)
    # -16- alp0 -
    V = TK.DoubleVar(win); V.set(0.); VARS.append(V)
    # -17- transl_speed -
    V = TK.StringVar(win); V.set('0.; 0.; 0.'); VARS.append(V)
    # -18- rot_pnt -
    V = TK.StringVar(win); V.set('0.; 0.; 0.'); VARS.append(V)
    # -19- rot_vct -
    V = TK.StringVar(win); V.set('0.; 0.; 1.'); VARS.append(V)
    # -20- rot_omg -
    V = TK.DoubleVar(win); V.set(0.); VARS.append(V)
    # -21- del_pnt -
    V = TK.StringVar(win); V.set('0.; 0.; 0.'); VARS.append(V)
    # -22- del_vct -
    V = TK.StringVar(win); V.set('0.; 0.; 1.'); VARS.append(V)
    # -23- del0 -
    V = TK.DoubleVar(win); V.set(0.); VARS.append(V)
    # -24- delc -
    V = TK.StringVar(win); V.set('0.; 0.; 0.'); VARS.append(V)
    # -25- dels -
    V = TK.StringVar(win); V.set('0.; 0.; 0.'); VARS.append(V)
    # -26- bet_pnt -
    V = TK.StringVar(win); V.set('0.; 0.; 0.'); VARS.append(V)
    # -27- bet_vct -
    V = TK.StringVar(win); V.set('0.; 0.; 1.'); VARS.append(V)
    # -28- bet0 -
    V = TK.DoubleVar(win); V.set(0.); VARS.append(V)
    # -29- betc -
    V = TK.StringVar(win); V.set('0.; 0.; 0.'); VARS.append(V)
    # -30- bets -
    V = TK.StringVar(win); V.set('0.; 0.; 0.'); VARS.append(V)
    # -31- tet_pnt -
    V = TK.StringVar(win); V.set('0.; 0.; 0.'); VARS.append(V)
    # -32- tet_vct -
    V = TK.StringVar(win); V.set('1.; 0.; 0.'); VARS.append(V)
    # -33- tet0 -
    V = TK.DoubleVar(win); V.set(0.); VARS.append(V)
    # -34- tetc -
    V = TK.StringVar(win); V.set('0.'); VARS.append(V)
    # -35- tets -
    V = TK.StringVar(win); V.set('0.'); VARS.append(V)
    # -36- span_vct -
    V = TK.StringVar(win); V.set('1.; 0.; 0.'); VARS.append(V)
    # -37- pre_lag_pnt -
    V = TK.StringVar(win); V.set('0.; 0.; 0.'); VARS.append(V)
    # -38- pre_lag_vct -
    V = TK.StringVar(win); V.set('1.; 0.; 0.'); VARS.append(V)
    # -39- pre_lag_ang -
    V = TK.DoubleVar(win); V.set(0.); VARS.append(V)
    # -40- pre_con_pnt -
    V = TK.StringVar(win); V.set('0.; 0.; 0.'); VARS.append(V)
    # -41- pre_con_vct -
    V = TK.StringVar(win); V.set('1.; 0.; 0.'); VARS.append(V)
    # -42- pre_con_ang -
    V = TK.DoubleVar(win); V.set(0.); VARS.append(V)
    #----------------------
    # -- Constant motion --
    #----------------------
    # -43- transl_speed -
    V = TK.StringVar(win); V.set('0.; 0.; 0.'); VARS.append(V)
    # -44- axis_pnt -
    V = TK.StringVar(win); V.set('0.; 0.; 0.'); VARS.append(V)
    # -45- axis_vct -
    V = TK.StringVar(win); V.set('0.; 0.; 1.'); VARS.append(V)
    # -46- omega -
    V = TK.DoubleVar(win); V.set(0.); VARS.append(V)

    # - Settings -
    B = TTK.OptionMenu(Frame, VARS[11], '1:MotionStrings', '2:KMotionRotor',
                       '3:ConstantMotion',
                       command=changeMotionType)
    B.grid(row=0, column=0, columnspan=3, sticky=TK.EW)

    B = TTK.Label(Frame, text="Name:")
    BB = CTK.infoBulle(parent=B, text='Name of motion.')
    B.grid(row=1, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White')
    B.grid(row=1, column=1, columnspan=2, sticky=TK.EW)

    # - Set/Get -
    B = TTK.Button(Frame, text="Get", command=getVars)
    B.grid(row=3, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Get motion from zone (if exists).')
    B = TTK.Button(Frame, text="Set", command=setVars)
    B.grid(row=3, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set motion in zone.')
    B = TTK.Button(Frame, text="Reset", command=resetVars)
    B.grid(row=3, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Reset all to default.')

    # Frame pour chaque type de motion
    Frame1 = TTK.Frame(Frame, borderwidth=0)
    Frame1.columnconfigure(0, weight=0)
    Frame1.columnconfigure(1, weight=1)
    WIDGETS['frame1'] = Frame1
    Frame2 = TTK.Frame(Frame, borderwidth=0)
    Frame2.columnconfigure(0, weight=0)
    Frame2.columnconfigure(1, weight=1)
    WIDGETS['frame2'] = Frame2
    Frame3 = TTK.Frame(Frame, borderwidth=0)
    Frame3.columnconfigure(0, weight=0)
    Frame3.columnconfigure(1, weight=1)
    WIDGETS['frame3'] = Frame3

    # - Motion Type1: motion strings -
    B = TTK.Label(Frame1, text="tx:")
    BB = CTK.infoBulle(parent=B,
                       text='Translation of origin.\nCan depend on {t}.')
    B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame1, textvariable=VARS[1], background='White')
    BB = CTK.infoBulle(parent=B,
                       text='Translation of origin.\nCan depend on {t}.')
    B.grid(row=0, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars1)

    B = TTK.Label(Frame1, text="ty:")
    BB = CTK.infoBulle(parent=B,
                       text='Translation of origin.\nCan depend on {t}.')
    B.grid(row=1, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame1, textvariable=VARS[2], background='White')
    BB = CTK.infoBulle(parent=B,
                       text='Translation of origin.\nCan depend on {t}.')
    B.grid(row=1, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars1)

    B = TTK.Label(Frame1, text="tz:")
    BB = CTK.infoBulle(parent=B,
                       text='Translation of origin.\nCan depend on {t}.')
    B.grid(row=2, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame1, textvariable=VARS[3], background='White')
    BB = CTK.infoBulle(parent=B,
                       text='Translation of origin.\nCan depend on {t}.')
    B.grid(row=2, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars1)

    B = TTK.Button(Frame1, text="Set", command=setTransOrigin)
    B.grid(row=1, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set translation of origin from mouse.')

    B = TTK.Label(Frame1, text="cx:")
    BB = CTK.infoBulle(parent=B,
                       text='Center of rotation.\nCan depend on {t}.')
    B.grid(row=4, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame1, textvariable=VARS[4], background='White')
    BB = CTK.infoBulle(parent=B,
                       text='Center of rotation.\nCan depend on {t}.')
    B.grid(row=4, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars1)

    B = TTK.Label(Frame1, text="cy:")
    BB = CTK.infoBulle(parent=B,
                       text='Center of rotation.\nCan depend on {t}.')
    B.grid(row=5, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame1, textvariable=VARS[5], background='White')
    BB = CTK.infoBulle(parent=B,
                       text='Center of rotation.\nCan depend on {t}.')
    B.grid(row=5, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars1)

    B = TTK.Label(Frame1, text="cz:")
    BB = CTK.infoBulle(parent=B,
                       text='Center of rotation.\nCan depend on {t}.')
    B.grid(row=6, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame1, textvariable=VARS[6], background='White')
    B.grid(row=6, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars1)

    B = TTK.Button(Frame1, text="Set", command=setCenterRotation)
    B.grid(row=5, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set rotation center from mouse.')

    B = TTK.Label(Frame1, text="ex:")
    BB = CTK.infoBulle(parent=B,
                       text='Second rot axis point.\nCan depend on {t}.')
    B.grid(row=7, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame1, textvariable=VARS[7], background='White')
    BB = CTK.infoBulle(parent=B,
                       text='Second rot axis point.\nCan depend on {t}.')
    B.grid(row=7, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars1)

    B = TTK.Label(Frame1, text="ey:")
    BB = CTK.infoBulle(parent=B,
                       text='Second rot axis point.\nCan depend on {t}.')
    B.grid(row=8, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame1, textvariable=VARS[8], background='White')
    BB = CTK.infoBulle(parent=B,
                       text='Second rot axis point.\nCan depend on {t}.')
    B.grid(row=8, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars1)

    B = TTK.Label(Frame1, text="ez:")
    BB = CTK.infoBulle(parent=B,
                       text='Second rot axis point.\nCan depend on {t}.')
    B.grid(row=9, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame1, textvariable=VARS[9], background='White')
    BB = CTK.infoBulle(parent=B,
                       text='Second rot axis point.\nCan depend on {t}.')
    B.grid(row=9, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars1)

    B = TTK.Button(Frame1, text="Set", command=set2Axis)
    B.grid(row=8, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set rot axis from mouse.')

    B = TTK.Label(Frame1, text="angle:")
    BB = CTK.infoBulle(parent=B,
                       text='Rotation angle.\nCan depend on {t}. In degrees.')
    B.grid(row=10, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame1, textvariable=VARS[10], background='White')
    BB = CTK.infoBulle(parent=B,
                       text='Rotation angle.\nCan depend on {t}. In degrees.')
    B.grid(row=10, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars1)

    # - Motion Type2: KMotionRotor -
    B = TTK.Label(Frame2, text="transl_speed:")
    BB = CTK.infoBulle(parent=B, text='Translation speed and span_vect.')
    B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame2, textvariable=VARS[17], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='transl_speed.')
    B.grid(row=0, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars2)

    B = TTK.Entry(Frame2, textvariable=VARS[36], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='span vector.')
    B.grid(row=0, column=2, sticky=TK.EW)
    B.bind('<Return>', setVars2)

    B = TTK.Label(Frame2, text="psi0:")
    BB = CTK.infoBulle(parent=B, text='Azymuth at t=0.')
    B.grid(row=1, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame2, textvariable=VARS[12], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='psi0. In degrees.')
    B.grid(row=1, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars2)
    B = TTK.Entry(Frame2, textvariable=VARS[13], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='psi0_b. In degrees.')
    B.grid(row=1, column=2, sticky=TK.EW)
    B.bind('<Return>', setVars2)

    B = TTK.Label(Frame2, text="alp:")
    BB = CTK.infoBulle(parent=B, text='Shaft angle.')
    B.grid(row=2, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame2, textvariable=VARS[14], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='alp_pnt.')
    B.grid(row=2, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars2)
    B = TTK.Entry(Frame2, textvariable=VARS[15], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='alp_vct.')
    B.grid(row=2, column=2, sticky=TK.EW)
    B.bind('<Return>', setVars2)
    B = TTK.Entry(Frame2, textvariable=VARS[16], background='White', width=4)
    BB = CTK.infoBulle(parent=B, text='alp0. In degrees.')
    B.grid(row=2, column=3, sticky=TK.EW)
    B.bind('<Return>', setVars2)

    B = TTK.Label(Frame2, text="rot:")
    BB = CTK.infoBulle(parent=B, text='Rotation.')
    B.grid(row=3, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame2, textvariable=VARS[18], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='rot_pnt.')
    B.grid(row=3, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars2)
    B = TTK.Entry(Frame2, textvariable=VARS[19], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='rot_vct.')
    B.grid(row=3, column=2, sticky=TK.EW)
    B.bind('<Return>', setVars2)
    B = TTK.Entry(Frame2, textvariable=VARS[20], background='White', width=4)
    BB = CTK.infoBulle(parent=B, text='rot_omg. In rad/time unit.')
    B.grid(row=3, column=3, sticky=TK.EW)
    B.bind('<Return>', setVars2)

    B = TTK.Label(Frame2, text="del:")
    BB = CTK.infoBulle(parent=B, text='..')
    B.grid(row=4, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame2, textvariable=VARS[21], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='del_pnt.')
    B.grid(row=4, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars2)
    B = TTK.Entry(Frame2, textvariable=VARS[22], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='del_vct.')
    B.grid(row=4, column=2, sticky=TK.EW)
    B.bind('<Return>', setVars2)
    B = TTK.Entry(Frame2, textvariable=VARS[23], background='White', width=4)
    B.grid(row=4, column=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='del0. In degrees.')
    B.bind('<Return>', setVars2)
    B = TTK.Label(Frame2, text="delc,s:")
    BB = CTK.infoBulle(parent=B, text='del modes.')
    B.grid(row=5, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame2, textvariable=VARS[24], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='delc.')
    B.grid(row=5, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars2)
    B = TTK.Entry(Frame2, textvariable=VARS[25], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='dels.')
    B.grid(row=5, column=2, sticky=TK.EW)
    B.bind('<Return>', setVars2)

    B = TTK.Label(Frame2, text="bet:")
    BB = CTK.infoBulle(parent=B, text='..')
    B.grid(row=6, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame2, textvariable=VARS[26], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='bet_pnt.')
    B.grid(row=6, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars2)
    B = TTK.Entry(Frame2, textvariable=VARS[27], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='bet_vct.')
    B.grid(row=6, column=2, sticky=TK.EW)
    B.bind('<Return>', setVars2)
    B = TTK.Entry(Frame2, textvariable=VARS[28], background='White', width=4)
    B.grid(row=6, column=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='bet0. In degrees.')
    B.bind('<Return>', setVars2)
    B = TTK.Label(Frame2, text="betc,s:")
    BB = CTK.infoBulle(parent=B, text='bet modes.')
    B.grid(row=7, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame2, textvariable=VARS[29], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='betc.')
    B.grid(row=7, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars2)
    B = TTK.Entry(Frame2, textvariable=VARS[30], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='bets.')
    B.grid(row=7, column=2, sticky=TK.EW)
    B.bind('<Return>', setVars2)

    B = TTK.Label(Frame2, text="tet:")
    BB = CTK.infoBulle(parent=B, text='..')
    B.grid(row=8, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame2, textvariable=VARS[31], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='tet_pnt.')
    B.grid(row=8, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars2)
    B = TTK.Entry(Frame2, textvariable=VARS[32], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='tet_vct.')
    B.grid(row=8, column=2, sticky=TK.EW)
    B.bind('<Return>', setVars2)
    B = TTK.Entry(Frame2, textvariable=VARS[33], background='White', width=4)
    B.grid(row=8, column=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='tet0. In degrees.')
    B.bind('<Return>', setVars2)
    B = TTK.Label(Frame2, text="tetc,s:")
    BB = CTK.infoBulle(parent=B, text='tet modes.')
    B.grid(row=9, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame2, textvariable=VARS[34], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='tetc.')
    B.grid(row=9, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars2)
    B = TTK.Entry(Frame2, textvariable=VARS[35], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='tets.')
    B.grid(row=9, column=2, sticky=TK.EW)
    B.bind('<Return>', setVars2)

    B = TTK.Label(Frame2, text="pre_lag:")
    BB = CTK.infoBulle(parent=B, text='..')
    B.grid(row=10, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame2, textvariable=VARS[37], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='pre_lag_pnt.')
    B.grid(row=10, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars2)
    B = TTK.Entry(Frame2, textvariable=VARS[38], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='pre_lag_vct.')
    B.grid(row=10, column=2, sticky=TK.EW)
    B.bind('<Return>', setVars2)
    B = TTK.Entry(Frame2, textvariable=VARS[39], background='White', width=4)
    B.grid(row=10, column=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='pre_lag_ang. In degrees.')
    B.bind('<Return>', setVars2)

    B = TTK.Label(Frame2, text="pre_con:")
    BB = CTK.infoBulle(parent=B, text='..')
    B.grid(row=11, column=0, sticky=TK.EW)
    B = TK.Entry(Frame2, textvariable=VARS[40], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='pre_con_pnt.')
    B.grid(row=11, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars2)
    B = TK.Entry(Frame2, textvariable=VARS[41], background='White', width=10)
    BB = CTK.infoBulle(parent=B, text='pre_con_vct.')
    B.grid(row=11, column=2, sticky=TK.EW)
    B.bind('<Return>', setVars2)
    B = TK.Entry(Frame2, textvariable=VARS[42], background='White', width=4)
    B.grid(row=11, column=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='pre_con_ang. In degrees.')
    B.bind('<Return>', setVars2)

    # - Motion Type3: ConstantMotion -
    B = TTK.Label(Frame3, text="transl_speed:")
    BB = CTK.infoBulle(parent=B, text='Translation speed.')
    B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame3, textvariable=VARS[43], background='White')
    BB = CTK.infoBulle(parent=B, text='Translation speed.')
    B.grid(row=0, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars3)
    B = TTK.Label(Frame3, text="axis_pnt:")
    BB = CTK.infoBulle(parent=B, text='Rotation center.')
    B.grid(row=1, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame3, textvariable=VARS[44], background='White')
    BB = CTK.infoBulle(parent=B, text='Rotation center.')
    B.grid(row=1, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars3)
    B = TTK.Label(Frame3, text="axis_vct:")
    BB = CTK.infoBulle(parent=B, text='Rotation axis.')
    B.grid(row=2, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame3, textvariable=VARS[45], background='White')
    BB = CTK.infoBulle(parent=B, text='Rotation axis.')
    B.grid(row=2, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars3)
    B = TTK.Label(Frame3, text="omega:")
    BB = CTK.infoBulle(parent=B, text='Rotation speed. In rad/time unit.')
    B.grid(row=3, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame3, textvariable=VARS[46], background='White')
    BB = CTK.infoBulle(parent=B, text='Rotation speed. In rad/time unit.')
    B.grid(row=3, column=1, sticky=TK.EW)
    B.bind('<Return>', setVars3)

    changeMotionType()

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['MotionNoteBook'].add(WIDGETS['frame'], text='tkRigidMotion')
    except: pass
    CTK.WIDGETS['MotionNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['MotionNoteBook'].hide(WIDGETS['frame'])

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
    (win, menu, file, tools) = CTK.minimal('tkRigidMotion '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
