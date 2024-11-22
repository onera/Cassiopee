"""Module performing rigid motion.
"""
from . import RigidMotion
from . import rigidMotion
__version__ = RigidMotion.__version__

try:
    import Converter
    import Converter.PyTree as C
    import Converter.Internal as Internal
    import Transform.PyTree as T
except ImportError:
    raise ImportError("RigidMotion: requires Converter, Transform modules.")

import numpy
from math import cos, sin, sqrt, pi
# Stocke les functions deja definies
DEFINEDMOTIONS = {}

__DEG2RAD__= pi/180.
__RAD2DEG__= 180./pi

#=============================================================================
# Permet de definir un mouvement solide par des chaines (dependant de {t})
# definissant:
# la translation (tx,ty,tz), le centre de rotation (cx,cy,cz),
# et l'axe de rotation (ex,ey,ez)
# l angle correspond toujours a des degres
#=============================================================================
def setPrescribedMotion1(t, name, tx="0", ty="0", tz="0",
                         cx="0", cy="0", cz="0",
                         ex="0", ey="0", ez="1", angle="0"):
    """Define a motion of type 1 (time strings)."""
    tp = Internal.copyRef(t)
    _setPrescribedMotion1(tp, name, tx=tx, ty=ty, tz=tz, cx=cx, cy=cy, cz=cz, 
                          ex=ex, ey=ey, ez=ez, angle=angle)
    return tp

def _setPrescribedMotion1(t, name, tx="0", ty="0", tz="0",
                         cx="0", cy="0", cz="0",
                         ex="0", ey="0", ez="1", angle="0"):
    for z in Internal.getZones(t):
        # Recupere le conteneur TimeMotion
        cont = Internal.getNodeFromName1(z, 'TimeMotion')
        if cont is None:
            cont = ['TimeMotion', None, [], 'UserDefinedData_t']
            z[2].append(cont)

        # Le TimeRigidMotion name existe-t-il?
        motion = Internal.getNodeFromName1(cont, name)
        if motion is None:
            motion = [name, None, [], 'TimeRigidMotion_t']
            cont[2].append(motion)

        # Set it
        motion[2] = []
        motion[2].append(['MotionType', numpy.array([1], dtype=Internal.E_NpyInt), [], 'DataArray_t'])
        Internal._createChild(motion, 'tx', 'DataArray_t', value=tx)
        Internal._createChild(motion, 'ty', 'DataArray_t', value=ty)
        Internal._createChild(motion, 'tz', 'DataArray_t', value=tz)
        Internal._createChild(motion, 'cx', 'DataArray_t', value=cx)
        Internal._createChild(motion, 'cy', 'DataArray_t', value=cy)
        Internal._createChild(motion, 'cz', 'DataArray_t', value=cz)
        Internal._createChild(motion, 'ex', 'DataArray_t', value=ex)
        Internal._createChild(motion, 'ey', 'DataArray_t', value=ey)
        Internal._createChild(motion, 'ez', 'DataArray_t', value=ez)
        Internal._createChild(motion, 'angle', 'DataArray_t', value=angle)

    return None

#==================================================================
# Permet de definir un mouvement de RotorMotion 
#==================================================================
def setPrescribedMotion2(t, name,
                         transl_speed=(0.,0.,0.), psi0=0., psi0_b=0.,
                         alp_pnt=(0.,0.,0.), alp_vct=(0.,1.,0.), alp0=0.,
                         rot_pnt=(0.,0.,0.), rot_vct=(0.,0.,1.), rot_omg=0.,
                         del_pnt=(0.,0.,0.), del_vct=(0.,0.,1.), del0=0.,
                         delc=(0.,0.,0.), dels=(0.,0.,0.),
                         bet_pnt=(0.,0.,0.), bet_vct=(0.,1.,0.), bet0=0.,
                         betc=(0.,0.,0.), bets=(0.,0.,0.),
                         tet_pnt=(0.,0.,0.), tet_vct=(1.,0.,0.), tet0=0.,
                         tetc=(0.,), tets=(0.,),
                         span_vct=(1.,0.,0.),
                         pre_lag_pnt=(0.,0.,0.), pre_lag_vct=(1.,0.,0.), pre_lag_ang=0.,
                         pre_con_pnt=(0.,0.,0.), pre_con_vct=(1.,0.,0.), pre_con_ang=0.):
    """Define a motion of type 2 (rotor)"""
    tp = Internal.copyRef(t)
    _setPrescribedMotion2(tp, name, transl_speed=transl_speed, psi0=psi0, psi0_b=psi0_b,
                          alp_pnt=alp_pnt, alp_vct=alp_vct, alp0=alp0,
                          rot_pnt=rot_pnt, rot_vct=rot_vct, rot_omg=rot_omg,
                          del_pnt=del_pnt, del_vct=del_vct, del0=del0,
                          delc=delc, dels=dels,
                          bet_pnt=bet_pnt, bet_vct=bet_vct, bet0=bet0,
                          betc=betc, bets=bets,
                          tet_pnt=tet_pnt, tet_vct=tet_vct, tet0=tet0,
                          tetc=tetc, tets=tets,
                          span_vct=span_vct,
                          pre_lag_pnt=pre_lag_pnt, pre_lag_vct=pre_lag_vct, pre_lag_ang=pre_lag_ang,
                          pre_con_pnt=pre_con_pnt, pre_con_vct=pre_con_vct, pre_con_ang=pre_con_ang)
    return tp

def _setPrescribedMotion2(t, name,
                          transl_speed=(0.,0.,0.),# forward velocity in the abs frame
                          psi0=0., psi0_b=0.,
                          alp_pnt=(0.,0.,0.), alp_vct=(0.,1.,0.), alp0=0.,
                          rot_pnt=(0.,0.,0.), rot_vct=(0.,0.,1.), rot_omg=0.,
                          del_pnt=(0.,0.,0.), del_vct=(0.,0.,1.), del0=0.,
                          delc=(0.,), dels=(0.,),# harmonics for lead lag
                          bet_pnt=(0.,0.,0.), bet_vct=(0.,1.,0.), bet0=0.,
                          betc=(0.,), bets=(0.,),#harmonics for flapping
                          tet_pnt=(0.,0.,0.), tet_vct=(1.,0.,0.), tet0=0.,
                          tetc=(0.,), tets=(0.,),# harmonics for pitching
                          span_vct=(1.,0.,0.), # blade spanwise vct
                          pre_lag_pnt=(0.,0.,0.), pre_lag_vct=(1.,0.,0.), pre_lag_ang=0.,
                          pre_con_pnt=(0.,0.,0.), pre_con_vct=(1.,0.,0.), pre_con_ang=0.):
    for z in Internal.getZones(t):
        # Recupere le conteneur TimeMotion
        cont = Internal.getNodeFromName1(z, 'TimeMotion')
        if cont is None:
            cont = ['TimeMotion', None, [], 'UserDefinedData_t']
            z[2].append(cont)

        # Le TimeRigidMotion name existe-t-il?
        motion = Internal.getNodeFromName1(cont, name)
        if motion is None:
            motion = [name, None, [], 'TimeRigidMotion_t']
            cont[2].append(motion)

        # Set it
        motion[2] = []
        motion[2].append(['MotionType', numpy.array([2], dtype=Internal.E_NpyInt),
                          [], 'DataArray_t'])
        motion[2].append(['transl_speed', numpy.array([transl_speed[0], transl_speed[1], transl_speed[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['psi0', numpy.array([psi0], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['psi0_b', numpy.array([psi0_b], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['alp_pnt', numpy.array([alp_pnt[0], alp_pnt[1], alp_pnt[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['alp_vct', numpy.array([alp_vct[0], alp_vct[1], alp_vct[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['alp0', numpy.array([alp0], numpy.float64),
                          [], 'DataArray_t'])

        motion[2].append(['rot_pnt', numpy.array([rot_pnt[0], rot_pnt[1], rot_pnt[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['rot_vct', numpy.array([rot_vct[0], rot_vct[1], rot_vct[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['rot_omg', numpy.array([rot_omg], numpy.float64),
                          [], 'DataArray_t'])

        # LEAD LAG
        motion[2].append(['del_pnt', numpy.array([del_pnt[0], del_pnt[1], del_pnt[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['del_vct', numpy.array([del_vct[0], del_vct[1], del_vct[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['del0', numpy.array([del0], numpy.float64),
                          [], 'DataArray_t'])
        # harmonics for lead-lag
        DELHC = []
        for noh in range(len(delc)): DELHC.append(delc[noh])
        motion[2].append(['delc', numpy.array(DELHC, numpy.float64), [], 'DataArray_t'])
        DELHS = []
        for noh in range(len(dels)): DELHS.append(dels[noh])
        motion[2].append(['dels', numpy.array(DELHS, numpy.float64), [], 'DataArray_t'])

        # FLAPPING
        motion[2].append(['bet_pnt', numpy.array([bet_pnt[0], bet_pnt[1], bet_pnt[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['bet_vct', numpy.array([bet_vct[0], bet_vct[1], bet_vct[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['bet0', numpy.array([bet0], numpy.float64),
                          [], 'DataArray_t'])
        # harmonics for flapping
        DELHC = []
        for noh in range(len(betc)): DELHC.append(betc[noh])
        motion[2].append(['betc', numpy.array(DELHC, numpy.float64), [], 'DataArray_t'])
        DELHS = []
        for noh in range(len(bets)): DELHS.append(bets[noh])
        motion[2].append(['bets', numpy.array(DELHS, numpy.float64), [], 'DataArray_t'])

        # PITCHING
        motion[2].append(['tet_pnt', numpy.array([tet_pnt[0], tet_pnt[1], tet_pnt[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['tet_vct', numpy.array([tet_vct[0], tet_vct[1], tet_vct[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['tet0', numpy.array([tet0], numpy.float64),
                          [], 'DataArray_t'])
        DELHC = []
        for noh in range(len(tetc)): DELHC.append(tetc[noh])
        motion[2].append(['tetc', numpy.array(DELHC, numpy.float64), [], 'DataArray_t'])
        DELHS = []
        for noh in range(len(tets)): DELHS.append(tets[noh])
        motion[2].append(['tets', numpy.array(DELHS, numpy.float64), [], 'DataArray_t'])

        # blade spanwise vector 
        motion[2].append(['span_vct', numpy.array([span_vct[0], span_vct[1], span_vct[2]], numpy.float64),
                          [], 'DataArray_t'])

        # pre-lag 
        motion[2].append(['pre_lag_pnt', numpy.array([pre_lag_pnt[0], pre_lag_pnt[1], pre_lag_pnt[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['pre_lag_vct', numpy.array([pre_lag_vct[0], pre_lag_vct[1], pre_lag_vct[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['pre_lag_ang', numpy.array([pre_lag_ang], numpy.float64),
                          [], 'DataArray_t'])
        # pre-conicity
        motion[2].append(['pre_con_pnt', numpy.array([pre_con_pnt[0], pre_con_pnt[1], pre_con_pnt[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['pre_con_vct', numpy.array([pre_con_vct[0], pre_con_vct[1], pre_con_vct[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['pre_con_ang', numpy.array([pre_con_ang], numpy.float64),
                          [], 'DataArray_t'])
    return None

#=============================================================================
# Permet de definir un mouvement de rotation + translation a vitesse constante
#=============================================================================
def setPrescribedMotion3(t, name,
                         transl_speed=(0.,0.,0.),  
                         axis_pnt=(0.,0.,0.), axis_vct=(0.,0.,1.),
                         omega=0.):
    """Define a motion of type 3 (constant rotation+translation speed)."""
    tp = Internal.copyRef(t)
    _setPrescribedMotion3(tp, name, transl_speed=transl_speed, axis_pnt=axis_pnt,
                          axis_vct=axis_vct, omega=omega)
    return tp

def _setPrescribedMotion3(t, name, transl_speed=(0.,0.,0.),
                         axis_pnt=(0.,0.,0.), axis_vct=(0.,0.,1.),
                         omega=0.):
    """Define a motion of type 3 (constant rotation+translation speed)."""
    for z in Internal.getZones(t):
        # Recupere le conteneur TimeMotion
        cont = Internal.getNodeFromName1(z, 'TimeMotion')
        if cont is None:
            cont = ['TimeMotion', None, [], 'UserDefinedData_t']
            z[2].append(cont)

        # Le TimeRigidMotion name existe-t-il?
        motion = Internal.getNodeFromName1(cont, name)
        if motion is None:
            motion = [name, None, [], 'TimeRigidMotion_t']
            cont[2].append(motion)

        # Set it
        motion[2] = []
        motion[2].append(['MotionType', numpy.array([3], dtype=Internal.E_NpyInt),
                          [], 'DataArray_t'])
        motion[2].append(['transl_speed', numpy.array([transl_speed[0], transl_speed[1], transl_speed[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['axis_pnt', numpy.array([axis_pnt[0], axis_pnt[1], axis_pnt[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['axis_vct', numpy.array([axis_vct[0], axis_vct[1], axis_vct[2]], numpy.float64),
                          [], 'DataArray_t'])
        motion[2].append(['omega', numpy.array([omega], numpy.float64),
                          [], 'DataArray_t'])
    return None

#==============================================================================
# IN: m: motion node
# IN: string: chaine dependant de {t} a evaluer
# IN: time: instant
# OUT: valeur
#==============================================================================
def evalTimeString__(m, string, time):
    st = Internal.getNodeFromName1(m, string)
    st = Internal.getValue(st)
    st = st.replace('{t}', str(time))
    tx = eval(st)
    return tx

def getTimeString__(m, string):
    st = Internal.getNodeFromName1(m, string)
    st = Internal.getValue(st)
    return st

def evalTimeDerivativeString__(m, string, time):
    import Converter.expression as expr
    st = Internal.getNodeFromName1(m, string)
    st = Internal.getValue(st)
    a = expr.ast('{%s}'%string+'='+st); da = expr.derivate(a)
    #print(da)
    crds = Converter.array('t',2,1,1); crds[1][0,0] = time; crds[1][0,1] = time
    Converter._addVars(crds, string); a.run(crds)
    Converter._addVars(crds, 'd_'+string); da.run(crds, d_t=numpy.ones(1))
    #print(crds[1][2,0])
    return crds[1][2,0]

#==============================================================================
# IN: m: motion node
# IN: string: nom du noeud dont on retourne la valeur
#==============================================================================
def getNodeValue__(m, string):
    st = Internal.getNodeFromName1(m, string)
    return st[1]

def _moveZone__(z, time):
    cont = Internal.getNodeFromName1(z, 'TimeMotion')
    if cont is not None:
        motions = Internal.getNodesFromType1(cont, 'TimeRigidMotion_t')
        motions.reverse()
        for m in motions:
            mtype = Internal.getNodeFromName1(m, 'MotionType')
            dtype = mtype[1][0]
            if dtype == 1: # type 1: time string
                tx = evalTimeString__(m, 'tx', time)
                ty = evalTimeString__(m, 'ty', time)
                tz = evalTimeString__(m, 'tz', time)
                cx = evalTimeString__(m, 'cx', time)
                cy = evalTimeString__(m, 'cy', time)
                cz = evalTimeString__(m, 'cz', time)
                ex = evalTimeString__(m, 'ex', time)
                ey = evalTimeString__(m, 'ey', time)
                ez = evalTimeString__(m, 'ez', time)
                angle = evalTimeString__(m, 'angle', time)
                T._translate(z, (tx,ty,tz))
                if angle != 0:
                    angle = angle #*__RAD2DEG__
                    T._rotate(z, (cx,cy,cz), (ex-cx,ey-cy,ez-cz), angle, vectors=[])

            elif dtype == 2: # type 2: rotor_motion for helicopters in FF                
                transl_speed = Internal.getValue(Internal.getNodeFromName(m, 'transl_speed'))
                psi0 = Internal.getValue(Internal.getNodeFromName(m, 'psi0'))
                psi0_b = Internal.getValue(Internal.getNodeFromName(m, 'psi0_b'))
                alp_pnt = Internal.getValue(Internal.getNodeFromName(m, 'alp_pnt'))
                alp_vct = Internal.getValue(Internal.getNodeFromName(m, 'alp_vct'))
                alp0 = Internal.getValue(Internal.getNodeFromName(m, 'alp0'))
                rot_pnt = Internal.getValue(Internal.getNodeFromName(m, 'rot_pnt'))
                rot_vct = Internal.getValue(Internal.getNodeFromName(m, 'rot_vct'))
                rot_omg = Internal.getValue(Internal.getNodeFromName(m, 'rot_omg'))
                del_pnt = Internal.getValue(Internal.getNodeFromName(m, 'del_pnt'))
                del_vct = Internal.getValue(Internal.getNodeFromName(m, 'del_vct'))
                del0 = Internal.getValue(Internal.getNodeFromName(m, 'del0'))
                delc = getNodeValue__(m, 'delc')
                dels = getNodeValue__(m, 'dels')
                bet_pnt = Internal.getValue(Internal.getNodeFromName(m, 'bet_pnt'))
                bet_vct = Internal.getValue(Internal.getNodeFromName(m, 'bet_vct'))
                bet0 = Internal.getValue(Internal.getNodeFromName(m, 'bet0'))
                betc = getNodeValue__(m, 'betc')
                bets = getNodeValue__(m, 'bets')                
                betc = getNodeValue__(m, 'betc')
                bets = getNodeValue__(m, 'bets')
                tet_pnt = Internal.getValue(Internal.getNodeFromName(m, 'tet_pnt'))
                tet_vct = Internal.getValue(Internal.getNodeFromName(m, 'tet_vct'))
                tet0 = Internal.getValue(Internal.getNodeFromName(m, 'tet0'))
                tetc = getNodeValue__(m, 'tetc')
                tets = getNodeValue__(m, 'tets')
                span_vct = Internal.getValue(Internal.getNodeFromName(m, 'span_vct'))
                pre_lag_ang = Internal.getValue(Internal.getNodeFromName(m, 'pre_lag_ang'))
                pre_lag_pnt = Internal.getValue(Internal.getNodeFromName(m, 'pre_lag_pnt'))
                pre_lag_vct = Internal.getValue(Internal.getNodeFromName(m, 'pre_lag_vct'))
                pre_con_ang = Internal.getValue(Internal.getNodeFromName(m, 'pre_con_ang'))
                pre_con_pnt = Internal.getValue(Internal.getNodeFromName(m, 'pre_con_pnt'))
                pre_con_vct = Internal.getValue(Internal.getNodeFromName(m, 'pre_con_vct'))
                [r0,x0,rotMat,s0,omega] = rigidMotion._computeRotorMotionInfo(
                    time, transl_speed.tolist(), psi0, psi0_b,
                    alp_pnt.tolist(),alp_vct.tolist(),alp0,
                    rot_pnt.tolist(),rot_vct.tolist(),rot_omg,
                    del_pnt.tolist(),del_vct.tolist(),del0, delc.tolist(), dels.tolist(),
                    bet_pnt.tolist(), bet_vct.tolist(), bet0, betc.tolist(), bets.tolist(),
                    tet_pnt.tolist(), tet_vct.tolist(), tet0, tetc.tolist(), tets.tolist(),
                    span_vct.tolist(),
                    pre_lag_ang, pre_lag_pnt.tolist(), pre_lag_vct.tolist(),
                    pre_con_ang, pre_con_pnt.tolist(), pre_con_vct.tolist())
                coordsC = [x0[0], x0[1], x0[2]]
                coordsD = [r0[0], r0[1], r0[2]]
                GC = Internal.getNodeFromName(z,Internal.__GridCoordinates__)
                XN = Internal.getNodeFromName(GC,'CoordinateX')
                YN = Internal.getNodeFromName(GC,'CoordinateY')
                ZN = Internal.getNodeFromName(GC,'CoordinateZ')
                XI = Internal.getValue(XN)
                YI = Internal.getValue(YN)
                ZI = Internal.getValue(ZN)
                _moveN([XI,YI,ZI], coordsD, coordsC, rotMat)
                Internal.setValue(XN,XI)
                Internal.setValue(YN,YI)
                Internal.setValue(ZN,ZI)

            elif dtype == 3: # type 3: constant rotation+translation speed
                transl_speed = getNodeValue__(m, 'transl_speed')
                axis_pnt = getNodeValue__(m, 'axis_pnt')
                axis_vct = getNodeValue__(m, 'axis_vct')
                omega = getNodeValue__(m, 'omega')               
                cx = axis_pnt[0]
                cy = axis_pnt[1]
                cz = axis_pnt[2]
                ex = axis_vct[0]
                ey = axis_vct[1]
                ez = axis_vct[2]
                angle = omega[0]*time*__RAD2DEG__
                T._rotate(z, (cx,cy,cz), (ex,ey,ez), angle, vectors=[])
                T._translate(z, (transl_speed[0]*time,transl_speed[1]*time,transl_speed[2]*time))
            else:
                print("Warning: Motion type not found. Nothing done.")
    return None

# Recopie GridCoordinates#Init (s'il existe) dans GridCoordinates 
def copyGridInit2Grid(t):
    tp = Internal.copyRef(t)
    _copyGridInit2Grid(tp)
    return tp

def _copyGridInit2Grid(t):
    """Copy GridCoordinates#Init to GridCoordinates."""
    zones = Internal.getZones(t)
    for z in zones:
        gridInit = Internal.getNodeFromName1(z, 'GridCoordinates#Init')
        if gridInit:
            xcoord0 = Internal.getNodeFromName1(gridInit, 'CoordinateX')[1]
            ycoord0 = Internal.getNodeFromName1(gridInit, 'CoordinateY')[1]
            zcoord0 = Internal.getNodeFromName1(gridInit, 'CoordinateZ')[1]
            grid = Internal.getNodeFromName1(z, 'GridCoordinates')
            xcoord = Internal.getNodeFromName1(grid, 'CoordinateX')[1]
            ycoord = Internal.getNodeFromName1(grid, 'CoordinateY')[1]
            zcoord = Internal.getNodeFromName1(grid, 'CoordinateZ')[1]
            rigidMotion.copyCoords(xcoord0, ycoord0, zcoord0, xcoord, ycoord, zcoord)
    return None

# Copy GridCoordinates dans GridCoordinates#Init
# si mode=0, ne recopie que si TimeMotion est present
# si mode=1, recopie toujours
def copyGrid2GridInit(t, mode=0):
    tp = Internal.copyRef(t)
    _copyGrid2GridInit(tp, mode)
    return tp

def _copyGrid2GridInit(t, mode=0):
    """Copy GridCoordinates to GridCoordinates#Init."""
    zones = Internal.getZones(t)
    for z in zones:
        if mode == 0:
            timeMotion = Internal.getNodeFromName1(z, 'TimeMotion')
            if timeMotion is None: continue
        grid = Internal.getNodeFromName1(z, 'GridCoordinates')
        xcoord = Internal.getNodeFromName1(grid, 'CoordinateX')
        ycoord = Internal.getNodeFromName1(grid, 'CoordinateY')
        zcoord = Internal.getNodeFromName1(grid, 'CoordinateZ')
        if grid:
            gridInit = Internal.getNodeFromName1(z, 'GridCoordinates#Init')
            if not gridInit:
                gridInit = Internal.createNode('GridCoordinates#Init', 'GridCoordinates_t', parent=z)
            xcoord0 = Internal.getNodeFromName1(gridInit, 'CoordinateX')
            if xcoord0 is None:
                xcoord0 = Internal.copyNode(xcoord); xcoord[2] = []; gridInit[2].append(xcoord0)
            ycoord0 = Internal.getNodeFromName1(gridInit, 'CoordinateY')
            if ycoord0 is None:
                ycoord0 = Internal.copyNode(ycoord); ycoord[2] = []; gridInit[2].append(ycoord0)
            zcoord0 = Internal.getNodeFromName1(gridInit, 'CoordinateZ')
            if zcoord0 is None:
                zcoord0 = Internal.copyNode(zcoord); zcoord[2] = []; gridInit[2].append(zcoord0)
            rigidMotion.copyCoords(xcoord[1], ycoord[1], zcoord[1], xcoord0[1], ycoord0[1], zcoord0[1])
            # copie CartesianData si il y en a
            cd = Internal.getNodeFromName1(grid, 'CartesianData')
            if cd is not None:
                Internal._createUniqueChild(gridInit, 'CartesianData', 'DataArray_t', value=Internal.getValue(cd))
    return None

# Switch GridCoordinates et GridCoordinates#Init si ils existent
def _switchGridAndGridInit(t):
    """Switch GridCoordinates and GridCoordinates#Init."""
    zones = Internal.getZones(t)
    for z in zones:
        gc1 = Internal.getNodeFromName1(z, 'GridCoordinates#Init')
        gc2 = Internal.getNodeFromName1(z, 'GridCoordinates')
        if gc1 is not None and gc2 is not None:
            gc1[0] = 'GridCoordinates'
            gc2[0] = 'GridCoordinates#Init'
    return None

# eval position by deformation: 
# ajoute Displacement#0 dans Coordinates#Init
def _evalDeformationPosition(t):
    zones = Internal.getZones(t)
    for z in zones:
        defcont = Internal.getNodeFromName1(z, "Displacement#0")
        if defcont is not None:
            cont = Internal.getNodeFromName1(z, "GridCoordinates#Init")
            XA1 = Internal.getNodeFromName1(cont, "CoordinateX")[1]
            XA2 = Internal.getNodeFromName1(cont, "CoordinateY")[1]
            XA3 = Internal.getNodeFromName1(cont, "CoordinateZ")[1]

            cont = Internal.getNodeFromName1(z, "GridCoordinates")
            XB1 = Internal.getNodeFromName1(cont, "CoordinateX")[1]
            XB2 = Internal.getNodeFromName1(cont, "CoordinateY")[1]
            XB3 = Internal.getNodeFromName1(cont, "CoordinateZ")[1]

            DA1 = Internal.getNodeFromName1(defcont, "DisplacementX")[1]
            DA2 = Internal.getNodeFromName1(defcont, "DisplacementY")[1]
            DA3 = Internal.getNodeFromName1(defcont, "DisplacementZ")[1]

            XB1[:] = XA1 + DA1
            XB2[:] = XA2 + DA2
            XB3[:] = XA3 + DA3
    return None

# copie Displacement#0/GridVelocityX dans Motion/VelocityX
def _evalDeformationSpeed(t):
    zones = Internal.getZones(t)
    for z in zones:
        defcont = Internal.getNodeFromName1(z, "Displacement#0")
        if defcont is not None:
            DA1 = Internal.getNodeFromName1(defcont, "VelocityX")[1]
            DA2 = Internal.getNodeFromName1(defcont, "VelocityY")[1]
            DA3 = Internal.getNodeFromName1(defcont, "VelocityZ")[1]
            motion = Internal.getNodeFromName1(z, "Displacement#0")
            if motion is None:
                motion = Internal.createNode('Motion', 'UserDefined_t', parent=z)
            sx = Internal.getNodeFromName1(motion, 'VelocityX')
            if sx is None: sx = Internal.copyNode(DA1); sx[0] = 'VelocityX'; motion[2].append(sx); sx[1] = sx[1].reshape((sx[1].size))
            sy = Internal.getNodeFromName1(motion, 'VelocityY')
            if sy is None: sy = Internal.copyNode(DA1); sy[0] = 'VelocityY'; motion[2].append(sy); sy[1] = sy[1].reshape((sy[1].size))
            sz = Internal.getNodeFromName1(motion, 'VelocityZ')
            if sz is None: sz = Internal.copyNode(DA1); sz[0] = 'VelocityZ'; motion[2].append(sz); sz[1] = sz[1].reshape((sz[1].size))
            sx[1][:] = DA1[1][:]
            sy[1][:] = DA2[1][:]
            sz[1][:] = DA3[1][:]
    return None

#==============================================================================
# Evalue la position reelle de la zone a l'instant t
# Le mouvement est stocke dans chaque zone.
# Les coordonnees de la zone sont modifiees
#==============================================================================
def _evalPosition__(a, time):
    _copyGridInit2Grid(a)
    zones = Internal.getZones(a)
    for z in zones: _moveZone__(z, time)
    return None

def evalPosition__(t, time):
    a = Internal.copyRef(t)
    zones = Internal.getZones(a)
    # copie GridCoordinates
    for z in zones:
        n = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
        np = Internal.copyTree(n)
        #for c, p in enumerate(n[2]): n[2][c] = np[2][c]
        n[2][:] = np[2][:]
    for z in zones: 
        _copyGridInit2Grid(z)
        _moveZone__(z, time)
    return a

#==============================================================================
# Evalue la position reelle de la zone a l'instant t
# Le mouvement est defini dans la fonction F.
# Les coordonnees de la zone sont modifiees
#==============================================================================
def evalPosition___(a, time, F):
    """Move the mesh with defined motion to time t. Return an array with
    moved mesh coordinates.
    Usage: evalPosition(a, time, F)"""
    return C.TZGC2(a, 'nodes', False, RigidMotion.evalPosition, time, F)

def _evalPosition___(a, time, F):
    """Move the mesh with defined motion to time t. Return an array with
    moved mesh coordinates.
    Usage: evalPosition(a, time, F)"""
    return C.__TZGC2(a, 'nodes', False, RigidMotion._evalPosition, time, F)

#==============================================================================
# Evalue la position reelle de la zone a l'instant t
# Si GridCoordinates#Init est present, il est utilise pour calculer la position
# sinon les coords dans a doivent correspondre a time=0
#==============================================================================
def evalPosition(a, time, F=None):
    """Move the mesh with defined motion to time t. Return an array with moved mesh coordinates."""
    if F is None: return evalPosition__(a, time)
    else: return evalPosition___(a, time, F)

def _evalPosition(a, time, F=None):
    """Move the mesh with defined motion to time t. Return an array with moved mesh coordinates."""
    if F is None: return _evalPosition__(a, time)
    else: return _evalPosition___(a, time, F)

#=========================================================
# Matrice de rotation a partir des donnees de l'arbre
#=========================================================
def getDictOfMotionMatrix(a, time, F=None):
    dictRM = {}
    for z in Internal.getZones(a):
        zname = Internal.getName(z)
        dictRM[zname] = getMotionMatrixForZone(z, time, F)
    return dictRM

def getMotionMatrixForZone(z, time, F=None):    
    if F is None:
        cont = Internal.getNodeFromName1(z, 'TimeMotion')
        if cont is not None:
            motions = Internal.getNodesFromType1(cont, 'TimeRigidMotion_t')
            for m in motions:
                mtype = Internal.getNodeFromName1(m, 'MotionType')
                dtype = mtype[1][0]
                if dtype == 1: # type 1: time string
                    cx = evalTimeString__(m, 'cx', time)
                    cy = evalTimeString__(m, 'cy', time)
                    cz = evalTimeString__(m, 'cz', time)
                    ex = evalTimeString__(m, 'ex', time)
                    ey = evalTimeString__(m, 'ey', time)
                    ez = evalTimeString__(m, 'ez', time)
                    theta = evalTimeString__(m, 'angle', time)
                    # theta doit etre en radians
                    theta = theta*__DEG2RAD__
                    return getRotationMatrix__(cx,cy,cz,ex,ey,ez,theta)

                elif dtype == 3: # type 3: constant transl + rotation speed
                    axis_pnt = getNodeValue__(m, 'axis_pnt')
                    axis_vct = getNodeValue__(m, 'axis_vct')
                    omega = getNodeValue__(m, 'omega')
                    cx = axis_pnt[0]
                    cy = axis_pnt[1]
                    cz = axis_pnt[2]
                    ex = axis_vct[0]
                    ey = axis_vct[1]
                    ez = axis_vct[2]
                    theta = omega[0]*time
                    # getRotationMatrix: theta en radians
                    return getRotationMatrix__(cx,cy,cz,ex,ey,ez,theta)

                elif dtype == 4: # ?
                    Rot = numpy.zeros((3,3), numpy.float64)
                    Rot[0,0]=1.; Rot[1,1]=1.; Rot[2,2]=1.
                    return Rot
                else:
                    raise ValueError("getMotionMatrixForZone: MotionType invalid.")
    else:
        raise ValueError("getMotionMatrixForZone: not yet implemented with a function.")

    Rot = numpy.zeros((3,3), numpy.float64)
    Rot[0,0]=1.; Rot[1,1]=1.; Rot[2,2]=1.
    return Rot

# getRotationMatrix: l'angle theta doit etre en radians
def getRotationMatrix__(cx,cy,cz,ex,ey,ez,theta):
    Rot = numpy.zeros((3,3), numpy.float64)
    vnorm = sqrt(ex*ex+ey*ey+ez*ez)
    if vnorm < 1.e-12: return Rot

    vnorm = 1./vnorm
    v1 = ex*vnorm; v2 = ey*vnorm; v3 = ez*vnorm

    t1 =  cos(theta)
    t2 =  1.-t1
    t3 =  v1*v1
    t6 =  t2*v1
    t7 =  t6*v2
    t8 =  sin(theta)
    t9 =  t8*v3
    t11 = t6*v3
    t12 = t8*v2
    t15 = v2*v2
    t19 = t2*v2*v3
    t20 = t8*v1
    t24 = v3*v3
    Rot[0,0] = t1 + t2*t3
    Rot[0,1] = t7 - t9
    Rot[0,2] = t11 + t12
    Rot[1,0] = t7 + t9
    Rot[1,1] = t1 + t2*t15
    Rot[1,2] = t19 - t20
    Rot[2,0] = t11 - t12
    Rot[2,1] = t19 + t20
    Rot[2,2] = t1 + t2*t24
    #print('Matrice de rotation de RigidMotion en x: ', Rot[0,:])
    #print('Matrice de rotation de RigidMotion en y: ', Rot[1,:])
    #print('Matrice de rotation de RigidMotion en z: ', Rot[2,:])
    return Rot

# Applique la formule XP=d+r*(XN-c) sur des numpys de coordonnees
# in place
def _moveN(coordsN, d, c, r):
    return RigidMotion._moveN(coordsN, d, c, r)

def moveN(coordsN, d, c, r):
    return RigidMotion.moveN(coordsN, d, c, r)

# evaluation de la vitesse de grille par des formules analytiques
def evalGridSpeed(a, time):
    """Eval grid speed at given time. Position must be already at time."""
    ap = Internal.copyRef(a)
    _evalGridSpeed(ap, time)
    return ap

# Evalue la vitesse a un instant t
# Les coords dans a doivent correspondre a time=time
# si out=1, sort aussi dans le FlowSolution (pour visu)
def _evalGridSpeed(a, time, out=0):
    """Eval grid speed at given time. Position must be already at time."""
    zones = Internal.getZones(a)
    for z in zones:
    # Find Coordinates pointers (must already be updated)
        grid = Internal.getNodeFromName1(z, 'GridCoordinates')
        xcoord = Internal.getNodeFromName1(grid, 'CoordinateX')
        ycoord = Internal.getNodeFromName1(grid, 'CoordinateY')
        zcoord = Internal.getNodeFromName1(grid, 'CoordinateZ')

        # Get translation vector
        cont = Internal.getNodeFromName1(z, 'TimeMotion')
        if cont is not None:
        # Get speed pointers
            mmo = Internal.getNodeFromName1(z, 'Motion')
            if mmo is None: mmo = Internal.createNode('Motion', 'UserDefined_t', parent=z)  
            sx = Internal.getNodeFromName1(mmo, 'VelocityX')
            if sx is None: sx = Internal.copyNode(xcoord); sx[0] = 'VelocityX'; mmo[2].append(sx); sx[1] = sx[1].reshape((sx[1].size))
            sy = Internal.getNodeFromName1(mmo, 'VelocityY')
            if sy is None: sy = Internal.copyNode(xcoord); sy[0] = 'VelocityY'; mmo[2].append(sy); sy[1] = sy[1].reshape((sy[1].size))
            sz = Internal.getNodeFromName1(mmo, 'VelocityZ')
            if sz is None: sz = Internal.copyNode(xcoord); sz[0] = 'VelocityZ'; mmo[2].append(sz); sz[1] = sz[1].reshape((sz[1].size))

            motions = Internal.getNodesFromType1(cont, 'TimeRigidMotion_t')
            for m in motions:
                mtype = Internal.getNodeFromName1(m, 'MotionType')
                dtype = mtype[1][0]

                if dtype == 1: # type 1: time string
                    tx = getTimeString__(m, 'tx')
                    ty = getTimeString__(m, 'ty')
                    tz = getTimeString__(m, 'tz')
                    cx = getTimeString__(m, 'cx')
                    cy = getTimeString__(m, 'cy')
                    cz = getTimeString__(m, 'cz')
                    ex = getTimeString__(m, 'ex')
                    ey = getTimeString__(m, 'ey')
                    ez = getTimeString__(m, 'ez')
                    theta = getTimeString__(m, 'angle')

                    # dtx, dty, dtz
                    dtx = evalTimeDerivativeString__(m, 'tx', time)
                    dty = evalTimeDerivativeString__(m, 'ty', time)
                    dtz = evalTimeDerivativeString__(m, 'tz', time)

                    # dcx, dcy, dcz
                    dcx = evalTimeDerivativeString__(m, 'cx', time)
                    dcy = evalTimeDerivativeString__(m, 'cy', time)
                    dcz = evalTimeDerivativeString__(m, 'cz', time)

                    sx[1][:] = 0.; sy[1][:] = 0.; sz[1][:] = 0.

                elif dtype == 2: # rotor motion
                    transl_speed=Internal.getValue(Internal.getNodeFromName(m, 'transl_speed'))
                    psi0 = Internal.getValue(Internal.getNodeFromName(m, 'psi0'))
                    psi0_b = Internal.getValue(Internal.getNodeFromName(m, 'psi0_b'))
                    alp_pnt = Internal.getValue(Internal.getNodeFromName(m, 'alp_pnt'))
                    alp_vct = Internal.getValue(Internal.getNodeFromName(m, 'alp_vct'))
                    alp0 = Internal.getValue(Internal.getNodeFromName(m, 'alp0'))
                    rot_pnt = Internal.getValue(Internal.getNodeFromName(m, 'rot_pnt'))
                    rot_vct = Internal.getValue(Internal.getNodeFromName(m, 'rot_vct'))
                    rot_omg = Internal.getValue(Internal.getNodeFromName(m, 'rot_omg'))
                    del_pnt = Internal.getValue(Internal.getNodeFromName(m, 'del_pnt'))
                    del_vct = Internal.getValue(Internal.getNodeFromName(m, 'del_vct'))
                    del0 = Internal.getValue(Internal.getNodeFromName(m, 'del0'))
                    delc = getNodeValue__(m, 'delc')
                    dels = getNodeValue__(m, 'dels')
                    bet_pnt = Internal.getValue(Internal.getNodeFromName(m, 'bet_pnt'))
                    bet_vct = Internal.getValue(Internal.getNodeFromName(m, 'bet_vct'))
                    bet0 = Internal.getValue(Internal.getNodeFromName(m, 'bet0'))
                    betc = getNodeValue__(m, 'betc')
                    bets = getNodeValue__(m, 'bets')                
                    betc = getNodeValue__(m, 'betc')
                    bets = getNodeValue__(m, 'bets')
                    tet_pnt = Internal.getValue(Internal.getNodeFromName(m, 'tet_pnt'))
                    tet_vct = Internal.getValue(Internal.getNodeFromName(m, 'tet_vct'))
                    tet0 = Internal.getValue(Internal.getNodeFromName(m, 'tet0'))
                    tetc = getNodeValue__(m, 'tetc')
                    tets = getNodeValue__(m, 'tets')
                    span_vct = Internal.getValue(Internal.getNodeFromName(m, 'span_vct'))
                    pre_lag_ang = Internal.getValue(Internal.getNodeFromName(m, 'pre_lag_ang'))
                    pre_lag_pnt = Internal.getValue(Internal.getNodeFromName(m, 'pre_lag_pnt'))
                    pre_lag_vct = Internal.getValue(Internal.getNodeFromName(m, 'pre_lag_vct'))
                    pre_con_ang = Internal.getValue(Internal.getNodeFromName(m, 'pre_con_ang'))
                    pre_con_pnt = Internal.getValue(Internal.getNodeFromName(m, 'pre_con_pnt'))
                    pre_con_vct = Internal.getValue(Internal.getNodeFromName(m, 'pre_con_vct'))
                    [r0,x0,rotMat,s0,omega]=rigidMotion._computeRotorMotionInfo(
                        time, transl_speed.tolist(), psi0, psi0_b,
                        alp_pnt.tolist(),alp_vct.tolist(),alp0,
                        rot_pnt.tolist(),rot_vct.tolist(),rot_omg,
                        del_pnt.tolist(),del_vct.tolist(),del0, delc.tolist(), dels.tolist(),
                        bet_pnt.tolist(), bet_vct.tolist(), bet0, betc.tolist(), bets.tolist(),
                        tet_pnt.tolist(), tet_vct.tolist(), tet0, tetc.tolist(), tets.tolist(),
                        span_vct.tolist(),
                        pre_lag_ang, pre_lag_pnt.tolist(), pre_lag_vct.tolist(),
                        pre_con_ang, pre_con_pnt.tolist(), pre_con_vct.tolist())
                    #print('omega', omega)
                    #print('s0', s0)
                    # sx = s0 - omega ^ x0 + omega ^ x
                    rigidMotion.evalGridMotionN([xcoord[1],ycoord[1],zcoord[1]],[sx[1],sy[1],sz[1]],(s0[0],s0[1],s0[2]),(x0[0],x0[1],x0[2]),(omega[0],omega[1],omega[2]))

                elif dtype == 3: # type 3: constant rotation / translation
                    transl_speed = Internal.getNodeFromName1(m, 'transl_speed')
                    if transl_speed is None: transl_speed = [0.,0.,0.]
                    else: transl_speed = transl_speed[1]
                    axis_pnt = Internal.getNodeFromName1(m, 'axis_pnt')
                    if axis_pnt is None: axis_pnt = [0.,0.,0.]
                    else: axis_pnt = axis_pnt[1]
                    axis_vct = Internal.getNodeFromName1(m, 'axis_vct')
                    if axis_vct is None: axis_vct = [0.,0.,1.]
                    else: axis_vct = axis_vct[1]
                    omega = Internal.getNodeFromName1(m, 'omega')
                    if omega is None: omega = 0.
                    else: omega = Internal.getValue(omega)

                    rigidMotion.evalSpeed3(xcoord[1], ycoord[1], zcoord[1],
                      sx[1], sy[1], sz[1], omega, omega*time,
                      transl_speed[0], transl_speed[1], transl_speed[2],
                      axis_pnt[0], axis_pnt[1], axis_pnt[2],
                      axis_vct[0], axis_vct[1], axis_vct[2])

            # Recopie le champ dans un FlowSolution pour visu
            if out == 1:
                C._initVars(z, 'GridVelocityX', 0.)
                C._initVars(z, 'GridVelocityY', 0.)
                C._initVars(z, 'GridVelocityZ', 0.)
                px = Internal.getNodeFromName2(z, 'GridVelocityX')[1]
                py = Internal.getNodeFromName2(z, 'GridVelocityY')[1]
                pz = Internal.getNodeFromName2(z, 'GridVelocityZ')[1]
                px.ravel('k')[:] = sx[1][:]
                py.ravel('k')[:] = sy[1][:]
                pz.ravel('k')[:] = sz[1][:]

    return None

# Evaluation de la position inverse
# coords: liste de numpys a inverser attaches a la zone correspondant a time=time
# z: zone contenant le motion
# retourne les coords a time=0
def evalPositionM1(coords, z, time):
    cont = Internal.getNodeFromName1(z, 'TimeMotion')
    coordsO = [numpy.copy(coords[0]),numpy.copy(coords[1]),numpy.copy(coords[2])]
    if cont is not None:
        motions = Internal.getNodesFromType1(cont, 'TimeRigidMotion_t')
        for m in motions:
            mtype = Internal.getNodeFromName1(m, 'MotionType')
            dtype = mtype[1][0]

            if dtype == 3: # constant rotation + translation speeds
                axis_pnt = getNodeValue__(m, 'axis_pnt')
                axis_vct = getNodeValue__(m, 'axis_vct')
                omega = getNodeValue__(m, 'omega')
                speed = getNodeValue__(m, 'transl_speed')
                coordsD = [-speed[0]*time+axis_pnt[0], -speed[1]*time+axis_pnt[1], -speed[2]*time+axis_pnt[2]]
                coordsC = [axis_pnt[0], axis_pnt[1], axis_pnt[2]]
                mat = getRotationMatrix__(axis_pnt[0],axis_pnt[1],axis_pnt[2],
                                          axis_vct[0],axis_vct[1],axis_vct[2],omega*time)
                mat = numpy.transpose(mat)
                _moveN(coordsO, coordsD, coordsC, mat)

            elif dtype == 2: # type 2: rotor_motion for helicopters in FF
                transl_speed=Internal.getValue(Internal.getNodeFromName(m, 'transl_speed'))
                psi0 = Internal.getValue(Internal.getNodeFromName(m, 'psi0'))
                psi0_b = Internal.getValue(Internal.getNodeFromName(m, 'psi0_b'))
                alp_pnt = Internal.getValue(Internal.getNodeFromName(m, 'alp_pnt'))
                alp_vct = Internal.getValue(Internal.getNodeFromName(m, 'alp_vct'))
                alp0 = Internal.getValue(Internal.getNodeFromName(m, 'alp0'))
                rot_pnt = Internal.getValue(Internal.getNodeFromName(m, 'rot_pnt'))
                rot_vct = Internal.getValue(Internal.getNodeFromName(m, 'rot_vct'))
                rot_omg = Internal.getValue(Internal.getNodeFromName(m, 'rot_omg'))
                del_pnt = Internal.getValue(Internal.getNodeFromName(m, 'del_pnt'))
                del_vct = Internal.getValue(Internal.getNodeFromName(m, 'del_vct'))
                del0 = Internal.getValue(Internal.getNodeFromName(m, 'del0'))
                delc = getNodeValue__(m, 'delc')
                dels = getNodeValue__(m, 'dels')
                bet_pnt = Internal.getValue(Internal.getNodeFromName(m, 'bet_pnt'))
                bet_vct = Internal.getValue(Internal.getNodeFromName(m, 'bet_vct'))
                bet0 = Internal.getValue(Internal.getNodeFromName(m, 'bet0'))
                betc = getNodeValue__(m, 'betc')
                bets = getNodeValue__(m, 'bets')                
                betc = getNodeValue__(m, 'betc')
                bets = getNodeValue__(m, 'bets')
                tet_pnt = Internal.getValue(Internal.getNodeFromName(m, 'tet_pnt'))
                tet_vct = Internal.getValue(Internal.getNodeFromName(m, 'tet_vct'))
                tet0 = Internal.getValue(Internal.getNodeFromName(m, 'tet0'))
                tetc = getNodeValue__(m, 'tetc')
                tets = getNodeValue__(m, 'tets')
                span_vct = Internal.getValue(Internal.getNodeFromName(m, 'span_vct'))
                pre_lag_ang = Internal.getValue(Internal.getNodeFromName(m, 'pre_lag_ang'))
                pre_lag_pnt = Internal.getValue(Internal.getNodeFromName(m, 'pre_lag_pnt'))
                pre_lag_vct = Internal.getValue(Internal.getNodeFromName(m, 'pre_lag_vct'))
                pre_con_ang = Internal.getValue(Internal.getNodeFromName(m, 'pre_con_ang'))
                pre_con_pnt = Internal.getValue(Internal.getNodeFromName(m, 'pre_con_pnt'))
                pre_con_vct = Internal.getValue(Internal.getNodeFromName(m, 'pre_con_vct'))
                [r0,x0,rotMat,s0,omega]=rigidMotion._computeRotorMotionInfo(
                    time, transl_speed.tolist(), psi0, psi0_b,
                    alp_pnt.tolist(), alp_vct.tolist(), alp0,
                    rot_pnt.tolist(), rot_vct.tolist(), rot_omg,
                    del_pnt.tolist(), del_vct.tolist(), del0, delc.tolist(), dels.tolist(),
                    bet_pnt.tolist(), bet_vct.tolist(), bet0, betc.tolist(), bets.tolist(),
                    tet_pnt.tolist(), tet_vct.tolist(), tet0, tetc.tolist(), tets.tolist(),
                    span_vct.tolist(),
                    pre_lag_ang, pre_lag_pnt.tolist(), pre_lag_vct.tolist(),
                    pre_con_ang, pre_con_pnt.tolist(), pre_con_vct.tolist())
                coordsD = [x0[0], x0[1], x0[2]]
                coordsC = [r0[0], r0[1], r0[2]]
                rotMat = numpy.transpose(rotMat)

                # XP=d+r*(XN-c)
                _moveN(coordsO, coordsD, coordsC, rotMat)
    return coordsO

# Computes the new coordinates and grid velocity for rotor motion
def _setRotorMotionCoordsAndGridVel(a, time):
    _copyGridInit2Grid(a)
    zones = Internal.getZones(a)
    for z in zones:
        cont = Internal.getNodeFromName1(z, 'TimeMotion')
        if cont is not None:
            motions = Internal.getNodesFromType1(cont, 'TimeRigidMotion_t')
            motions.reverse()
            for m in motions:
                mtype = Internal.getNodeFromName1(m, 'MotionType')
                dtype = mtype[1][0]
                if dtype == 2:
                    # Find Coordinates pointers (must already be updated)
                    grid = Internal.getNodeFromName1(z, 'GridCoordinates#Init')
                    if grid is None: grid = Internal.getNodeFromName1(z, 'GridCoordinates')
                    xcoord = Internal.getNodeFromName1(grid, 'CoordinateX')

                    # Get grid velocity pointers
                    name = 'Motion'
                    mmo = Internal.getNodeFromName1(z,name)
                    if mmo is None: mmo = Internal.createNode(name, 'UserDefined_t', parent=z)  
                    sx = Internal.getNodeFromName1(mmo, 'VelocityX')
                    if sx is None: sx = Internal.copyNode(xcoord); sx[0] = 'VelocityX'; mmo[2].append(sx); sx[1] = sx[1].reshape((sx[1].size))
                    sy = Internal.getNodeFromName1(mmo, 'VelocityY')
                    if sy is None: sy = Internal.copyNode(xcoord); sy[0] = 'VelocityY'; mmo[2].append(sy); sy[1] = sy[1].reshape((sy[1].size))
                    sz = Internal.getNodeFromName1(mmo, 'VelocityZ')
                    if sz is None: sz = Internal.copyNode(xcoord); sz[0] = 'VelocityZ'; mmo[2].append(sz); sz[1] = sz[1].reshape((sz[1].size))
                    transl_speed = Internal.getValue(Internal.getNodeFromName(m, 'transl_speed'))
                    psi0 = Internal.getValue(Internal.getNodeFromName(m, 'psi0'))
                    psi0_b = Internal.getValue(Internal.getNodeFromName(m, 'psi0_b'))
                    alp_pnt = Internal.getValue(Internal.getNodeFromName(m, 'alp_pnt'))
                    alp_vct = Internal.getValue(Internal.getNodeFromName(m, 'alp_vct'))
                    alp0 = Internal.getValue(Internal.getNodeFromName(m, 'alp0'))
                    rot_pnt = Internal.getValue(Internal.getNodeFromName(m, 'rot_pnt'))
                    rot_vct = Internal.getValue(Internal.getNodeFromName(m, 'rot_vct'))
                    rot_omg = Internal.getValue(Internal.getNodeFromName(m, 'rot_omg'))
                    del_pnt = Internal.getValue(Internal.getNodeFromName(m, 'del_pnt'))
                    del_vct = Internal.getValue(Internal.getNodeFromName(m, 'del_vct'))
                    del0 = Internal.getValue(Internal.getNodeFromName(m, 'del0'))
                    delc = getNodeValue__(m, 'delc')
                    dels = getNodeValue__(m, 'dels')
                    bet_pnt = Internal.getValue(Internal.getNodeFromName(m, 'bet_pnt'))
                    bet_vct = Internal.getValue(Internal.getNodeFromName(m, 'bet_vct'))
                    bet0 = Internal.getValue(Internal.getNodeFromName(m, 'bet0'))
                    betc = getNodeValue__(m, 'betc')
                    bets = getNodeValue__(m, 'bets')                
                    betc = getNodeValue__(m, 'betc')
                    bets = getNodeValue__(m, 'bets')
                    tet_pnt = Internal.getValue(Internal.getNodeFromName(m, 'tet_pnt'))
                    tet_vct = Internal.getValue(Internal.getNodeFromName(m, 'tet_vct'))
                    tet0 = Internal.getValue(Internal.getNodeFromName(m, 'tet0'))
                    tetc = getNodeValue__(m, 'tetc')
                    tets = getNodeValue__(m, 'tets')
                    span_vct = Internal.getValue(Internal.getNodeFromName(m, 'span_vct'))
                    pre_lag_ang = Internal.getValue(Internal.getNodeFromName(m, 'pre_lag_ang'))
                    pre_lag_pnt = Internal.getValue(Internal.getNodeFromName(m, 'pre_lag_pnt'))
                    pre_lag_vct = Internal.getValue(Internal.getNodeFromName(m, 'pre_lag_vct'))
                    pre_con_ang = Internal.getValue(Internal.getNodeFromName(m, 'pre_con_ang'))
                    pre_con_pnt = Internal.getValue(Internal.getNodeFromName(m, 'pre_con_pnt'))
                    pre_con_vct = Internal.getValue(Internal.getNodeFromName(m, 'pre_con_vct'))
                    rigidMotion._computeRotorMotionZ(
                        z,  sx[1], sy[1], sz[1], time, transl_speed.tolist(), psi0, psi0_b,
                        alp_pnt.tolist(), alp_vct.tolist(), alp0,
                        rot_pnt.tolist(), rot_vct.tolist(), rot_omg,
                        del_pnt.tolist(), del_vct.tolist(), del0, delc.tolist(), dels.tolist(),
                        bet_pnt.tolist(), bet_vct.tolist(), bet0, betc.tolist(), bets.tolist(),
                        tet_pnt.tolist(), tet_vct.tolist(), tet0, tetc.tolist(), tets.tolist(),
                        span_vct.tolist(),
                        pre_lag_ang, pre_lag_pnt.tolist(), pre_lag_vct.tolist(),
                        pre_con_ang, pre_con_pnt.tolist(), pre_con_vct.tolist(),
                        Internal.__GridCoordinates__,
                        Internal.__FlowSolutionNodes__,
                        Internal.__FlowSolutionCenters__)
    return None

# Evaluation de la vitesse par differences finies
# Il faut que les coordonnees dans a correspondent a time=0
# si il n'y a pas de GridCoordinates#Init
# si out=1, sort aussi dans le FlowSolution
def _evalGridSpeed2(a, time, dtime=1.e-6, out=0):
    zones = Internal.getZones(a)
    for z in zones:
        # evalue la position a l'instant time
        b = evalPosition(z, time)
        # evalue la position a l'instant time+dtime
        c = evalPosition(z, time+dtime)
        # evalue la difference
        coord1X = Internal.getNodeFromName(b, 'CoordinateX')[1] 
        coord1Y = Internal.getNodeFromName(b, 'CoordinateY')[1]
        coord1Z = Internal.getNodeFromName(b, 'CoordinateZ')[1]
        coord2X = Internal.getNodeFromName(c, 'CoordinateX')[1] 
        coord2Y = Internal.getNodeFromName(c, 'CoordinateY')[1] 
        coord2Z = Internal.getNodeFromName(c, 'CoordinateZ')[1]
        speedX = coord2X[:]-coord1X[:]
        speedY = coord2Y[:]-coord1Y[:]
        speedZ = coord2Z[:]-coord1Z[:]
        speedX[:] = speedX[:]/dtime
        speedY[:] = speedY[:]/dtime
        speedZ[:] = speedZ[:]/dtime
        # ajoute dans mmo
        cont = Internal.getNodeFromName1(z, 'TimeMotion')
        if cont is not None:
            # Get speed pointers
            xcoord = Internal.getNodeFromName2(z, 'CoordinateX')
            name  = 'Motion'
            mmo = Internal.getNodeFromName1(z, name)
            if mmo is None: mmo = Internal.createNode(name, 'UserDefined_t', parent=z)  
            sx = Internal.getNodeFromName1(mmo, 'VelocityX')
            if sx is None: sx = Internal.copyNode(xcoord); sx[0] = 'VelocityX'; mmo[2].append(sx); sx[1] = sx[1].reshape((sx[1].size))
            sy = Internal.getNodeFromName1(mmo, 'VelocityY')
            if sy is None: sy = Internal.copyNode(xcoord); sy[0] = 'VelocityY'; mmo[2].append(sy); sy[1] = sy[1].reshape((sy[1].size))
            sz = Internal.getNodeFromName1(mmo, 'VelocityZ')
            if sz is None: sz = Internal.copyNode(xcoord); sz[0] = 'VelocityZ'; mmo[2].append(sz); sz[1] = sz[1].reshape((sz[1].size))
            sx[1][:] = speedX.ravel('k')[:]
            sy[1][:] = speedY.ravel('k')[:]
            sz[1][:] = speedZ.ravel('k')[:]

        # ajoute dans FlowSolution
        if out == 1:
            C._initVars(z, 'GridVelocityX', 0.)
            C._initVars(z, 'GridVelocityY', 0.)
            C._initVars(z, 'GridVelocityZ', 0.)
            Internal.getNodeFromName(z, 'GridVelocityX')[1][:] = speedX[:]
            Internal.getNodeFromName(z, 'GridVelocityY')[1][:] = speedY[:]
            Internal.getNodeFromName(z, 'GridVelocityZ')[1][:] = speedZ[:]
    return None


#Updates the Coordinates of the immersed boundary points as the object rotates.
#Changes the values of the Target, Image, and Walls points based on the motion of the object.
#Only supports motion type 3 at the moment
#Prelim. test show that using vectorial Python has a negligible impact on the total run time.
#FastS.compute and the chimera blanking at each iteration require significantly more time.
#This routine is on the order of the time needed for the standard position routins for chimera; e.g.
#    R._evalPosition(tb, timeit)
#    R._evalPosition(t, timeit)
#    R._evalPosition(tc, timeit)
#    R._evalGridSpeed(t, timeit)

def _evalPositionIBC(tc,time):
    list_vars= ['C','I','W']
    for z in Internal.getZones(tc):
        cont = Internal.getNodeFromName1(z, 'TimeMotion')
        if cont is not None:
            motions = Internal.getNodesFromType1(cont, 'TimeRigidMotion_t')
            motions.reverse()
            transl_speed = Internal.getNodeFromName(motions[0], 'transl_speed')[1]
            axis_pnt     = Internal.getNodeFromName(motions[0], 'axis_pnt')[1]
            axis_vct     = Internal.getNodeFromName(motions[0], 'axis_vct')[1]
            omega        = Internal.getNodeFromName(motions[0], 'omega')[1]            
            cx           = axis_pnt[0]
            cy           = axis_pnt[1]
            cz           = axis_pnt[2]
            ex           = axis_vct[0]
            ey           = axis_vct[1]
            ez           = axis_vct[2]
            angle        = omega[0]*time

            a11         = numpy.cos(angle)            + ex**2*(1.-numpy.cos(angle))
            a12         = ex*ey*(1.-numpy.cos(angle)) - ez*numpy.sin(angle)
            a13         = ex*ez*(1.-numpy.cos(angle)) + ey*numpy.sin(angle)

            a21         = ex*ey*(1.-numpy.cos(angle)) + ez*numpy.sin(angle)
            a22         = numpy.cos(angle)            + ey**2*(1.-numpy.cos(angle))
            a23         = ey*ez*(1.-numpy.cos(angle)) - ex*numpy.sin(angle)

            a31         = ex*ez*(1.-numpy.cos(angle)) - ey*numpy.sin(angle)
            a32         = ey*ez*(1.-numpy.cos(angle)) + ex*numpy.sin(angle)
            a33         = numpy.cos(angle)            + ez**2*(1.-numpy.cos(angle))

            for subz in Internal.getNodesFromType(z,'ZoneSubRegion_t'):
                if subz[0][0:3]=="IBC":
                    for var in list_vars:
                        coordX_PC     = Internal.getNodeFromName(subz,'CoordinateX_P'+var)[1]
                        coordX_PCInit = Internal.getNodeFromName(subz,'CoordinateX_P'+var+'#Init')[1]
                        coordY_PC     = Internal.getNodeFromName(subz,'CoordinateY_P'+var)[1]
                        coordY_PCInit = Internal.getNodeFromName(subz,'CoordinateY_P'+var+'#Init')[1]
                        coordZ_PC     = Internal.getNodeFromName(subz,'CoordinateZ_P'+var)[1]
                        coordZ_PCInit = Internal.getNodeFromName(subz,'CoordinateZ_P'+var+'#Init')[1]

                        coordX_PC[:] = transl_speed[0]*time + a11*(coordX_PCInit[:]-cx)+a12*(coordY_PCInit[:]-cy)+a13*(coordZ_PCInit[:]-cz)+cx
                        coordY_PC[:] = transl_speed[1]*time + a21*(coordX_PCInit[:]-cx)+a22*(coordY_PCInit[:]-cy)+a23*(coordZ_PCInit[:]-cz)+cy
                        coordZ_PC[:] = transl_speed[2]*time + a31*(coordX_PCInit[:]-cx)+a32*(coordY_PCInit[:]-cy)+a33*(coordZ_PCInit[:]-cz)+cz
    return None

