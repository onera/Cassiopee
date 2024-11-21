# - tkFastSolver -
"""Interface to use Fast solvers."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal
try: import Fast.PyTree as Fast
except: pass
import CPlot.iconics as iconics
import math

# local widgets list
WIDGETS = {}; VARS = []

# Store main tree
MAINTREE = None
# Store body tree
BODY = None
# tkSlice module
TKSLICE = None
# WALL extraction zones
WALL = None
# tkPlotXY Desktop keeping wall data
DESKTOP1 = None
# The current number of run
NITRUN = 0
# LOAD extraction zones
LOAD = None
# tkPlotXY Desktop keeping load data
DESKTOP2 = None
# CFL
CFL = 0.7
# TIMESTEP
TIMESTEP = 0.002

#==============================================================================
def changeMode(event=None):
    global BODY
    global MAINTREE
    mode = VARS[10].get()
    if mode == 'Body':
        MAINTREE = CTK.t
        if BODY is not None:
            CTK.t = BODY
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.TKTREE.updateApp()
            CTK.display(CTK.t)
        CTK.TXT.insert('START', 'Revert to body tree.\n')
    elif mode == 'PrevStep':
        # Reload from restart.cgns
        try:
            CTK.t = C.convertFile2PyTree('restart.cgns')
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.TKTREE.updateApp()
            CTK.display(CTK.t)
            CTK.TXT.insert('START', 'Revert to previous solution step.\n')
        except: pass
        VARS[10].set('Main')
    else:
        if MAINTREE is not None:
            CTK.t = MAINTREE
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.TKTREE.updateApp()
            CTK.display(CTK.t)
        CTK.TXT.insert('START', 'Display main computation tree.\n')

def bodyMode(event=None):
    global BODY
    global MAINTREE

    MAINTREE = CTK.t
    if BODY is not None:
        CTK.t = BODY
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CTK.TKTREE.updateApp()
        CTK.display(CTK.t)
    CTK.TXT.insert('START', 'Revert to body tree.\n')
    VARS[10].set('Body')

def mainTreeMode(event=None):
    global BODY
    global MAINTREE

    if MAINTREE is not None:
        CTK.t = MAINTREE
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CTK.TKTREE.updateApp()
        CTK.display(CTK.t)
        CPlot.setState(mode=3, scalarField='Density')
        CTK.TXT.insert('START', 'Display main computation tree.\n')
        VARS[10].set('Main')

def reloadPrevStep(event=None):
    global BODY
    global MAINTREE
    global NITRUN

    try:
        CTK.t = C.convertFile2PyTree('restart.cgns')
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CTK.TKTREE.updateApp()
        CTK.display(CTK.t)
        CTK.TXT.insert('START', 'Revert to previous solution step.\n')
        VARS[10].set('Main')
        NITRUN -= 1
    except: 
        CTK.TXT.insert('START', 'restart file not found.\n')
        VARS[10].set('Body')

#==============================================================================
# Set data in selected zones
#==============================================================================
def setData():
    global BODY
    if CTK.t == []: return
    mode = VARS[10].get()
    if mode == 'Body': BODY = Internal.copyRef(CTK.t)

    temporal_scheme = VARS[0].get()
    ss_iteration = int(VARS[1].get())
    scheme = VARS[4].get()
    time_step = VARS[5].get()
    timeVal = VARS[11].get() # "cfl" or "time_step"

    numb = {'temporal_scheme':temporal_scheme,
            'ss_iteration':ss_iteration}
    numz = {'scheme':scheme, 'senseurType':0}

    if timeVal == "time_step": numz['time_step'] = time_step
    if timeVal == "cfl": numz['cfl'] = time_step

    nzs = CPlot.getSelectedZones()
    CTK.saveTree()

    if nzs == []:
        Fast._setNum2Base(CTK.t, numb)
        Fast._setNum2Zones(CTK.t, numz)
        CTK.TXT.insert('START', 'Solver data set in all bases.\n')

    else:
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            z = CTK.t[2][nob][2][noz]
            b, c = Internal.getParentOfNode(CTK.t, z)
            Fast._setNum2Base(b, numb)
            Fast._setNum2Zones(z, numz)
        CTK.TXT.insert('START', 'Solver data set in selection.\n')

#=============================================================================
# Modifie le body, met le body modifie dans CTK.t
#=============================================================================
def updateBodyAndPrepare():
    import RigidMotion.PyTree as RigidMotion

    bodyTime = VARS[13].get()
    steps = VARS[14].get()
    delta = (1.-0.)/steps
    bodyTime += delta
    VARS[13].set(bodyTime)
    print(bodyTime)

    global BODY
    if BODY is None: return
    # Creation de l'arbre de reprise
    tinit = CTK.t
    zvars = ['centers:Density_M1', 'centers:VelocityX_M1', 'centers:VelocityY_M1', 'centers:VelocityZ_M1', 'centers:Temperature_M1']
    zvars += ['centers:Density_P1', 'centers:VelocityX_P1', 'centers:VelocityY_P1', 'centers:VelocityZ_P1', 'centers:Temperature_P1']
    zvars += ['centers:TurbulentDistance', 'centers:cellN', 'centers:ViscosityEddy']
    C._rmVars(tinit, zvars)
    CTK.t = None

    # Apply motion or something to body
    CTK.t = RigidMotion.evalPosition(BODY, bodyTime)

    # Regenere le prep de CTK.t
    prepare(tinit)

    dim = getDim(CTK.t)
    if dim == 2: CPlot.display(CTK.t)
    else: displaySlices()

#==============================================================================
# Get data from selected zone
#==============================================================================
def getData():
    global TIMESTEP, CFL
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        zone = Internal.getNodeFromType2(CTK.t, 'Zone_t')
    else: # get first of selection
        nz = nzs[0]
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        zone = CTK.t[2][nob][2][noz]
    if zone is not None:
        n = Internal.getNodeFromPath(zone, '.Solver#define/time_step')
        if n is not None:
            val = Internal.getValue(n)
            VARS[11].set("time_step")
            VARS[5].set(val)
            TIMESTEP = val
        n = Internal.getNodeFromPath(zone, '.Solver#define/cfl')
        if n is not None:
            val = Internal.getValue(n)
            VARS[11].set("cfl")
            VARS[5].set(val)
            CFL = val
        d, c = Internal.getParentOfNode(CTK.t, zone)
        n = Internal.getNodeFromPath(d, '.Solver#define/temporal_scheme')
        if n is not None:
            val = Internal.getValue(n)
            VARS[0].set(val)
        n = Internal.getNodeFromPath(zone, '.Solver#define/scheme')
        if n is not None:
            val = Internal.getValue(n)
            VARS[4].set(val)


#==============================================================================
# get dim from equation set
#==============================================================================
def getDim(t):
    dim = Internal.getNodeFromName2(t, 'FlowEquationSet')
    if dim is not None:
        dim = Internal.getNodeFromName1(dim, 'EquationDimension')
        dim = Internal.getValue(dim)
    else: dim = 3
    return dim

#==============================================================================
def run(event=None):
    CTK.setCursor(2, WIDGETS['compute'])
    dim = getDim(CTK.t)
    mode = VARS[10].get()
    if mode == 'Body':
        global BODY
        BODY = Internal.copyRef(CTK.t)
        try:
            prepare() # save t, tc
            CTK.TXT.insert('START', 'Prepare OK.\n')
        except:
            CTK.setCursor(0, WIDGETS['compute'])    
            CTK.TXT.insert('START', 'Prepare failed.\n')
            return

        VARS[10].set('Main')
        CTK.t = CTK.upgradeTree(CTK.t)
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CTK.TKTREE.updateApp()
        if dim == 2: CTK.display(CTK.t)
        else: displaySlices()

    # Set CPlot to scalar mode to monitor solution
    if CPlot.getState('mode') == 0: # mesh
        CPlot.setState(mode=3, scalarField='Density')

    # get number of runs
    nruns = VARS[15].get()

    # Launch computations of nit iterations
    for n in range(nruns):
        compute(n)
        # stop if NAN to enable previous
        isfinite = C.isFinite(CTK.t)
        if not isfinite:
            CTK.TXT.insert('START', 'NAN detected. Use PrevStep in mode.\n')
            break

        if dim == 2: displayByReplace(CTK.t)
        else: displaySlices()

    CTK.TKTREE.updateApp()
    CTK.setCursor(0, WIDGETS['compute'])
    return None

# Display replacing all zones in place
# Prend moins de memoire
def displayByReplace(t):
    for nob, b in enumerate(t[2]):
        if b[3] == 'CGNSBase_t':
            for noz, z in enumerate(b[2]):
                if z[3] == 'Zone_t':
                    CPlot.replace(t, nob, noz, z)
    CPlot.render()

#==============================================================================
# A partir de CTK.t considere comme les bodies
# tinit est un arbre de reprise eventuel
def prepare(tinit=None):
    if CTK.t == []: return

    # Save preventif
    C.convertPyTree2File(CTK.t, 'body.cgns')

    # check state
    state = Internal.getNodeFromName2(CTK.t, 'ReferenceState')
    if state is None:
        CTK.TXT.insert('START', 'No state in case (tkState)')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # check dim
    dim = Internal.getNodeFromName2(CTK.t, 'FlowEquationSet')
    if dim is not None:
        dim = Internal.getNodeFromName1(dim, 'EquationDimension')
        dim = Internal.getValue(dim)
        if dim == 2: 
            # in 2D, the case must be in XY plane
            C._initVars(CTK.t, 'CoordinateZ', 0.)

    # Recupere la base REFINE
    b = Internal.getNodeFromName1(CTK.t, 'REFINE')
    if b is not None:
        tbox = C.newPyTree()
        tbox[2].append(b)
        Internal._rmNodesFromName1(CTK.t, 'REFINE')
    else: tbox = None

    import Apps.Fast.IBM as App
    myApp = App.IBM(format='single')
    myApp.input_var.vmin = 21
    myApp.input_var.tbox = tbox
    myApp.input_var.check = False
    myApp.input_var.tinit = tinit

    CTK.t, tc = myApp.prepare(CTK.t, t_out='t.cgns', tc_out='tc.cgns')

    # Preparation pour le front 42
    #CTK.t, tc = myApp.prepare(CTK.t, t_out='t.cgns', tc_out='tc.cgns', vmin=21, 
    #                          tbox=tbox, check=False, tinit=tinit, frontType=42, yplus=150.)

    return None

#==============================================================================
# Lance des iterations
#==============================================================================
def compute(nrun):
    if CTK.t == []: return

    import Apps.Fast.IBM as App
    global NITRUN # numero courant du run

    # Save preventif avec compression cartesienne
    Fast.saveFile(CTK.t, 'restart.cgns', compress=2)

    # check state
    state = Internal.getNodeFromName2(CTK.t, 'ReferenceState')
    if state is None:
        CTK.TXT.insert('START', 'No state in case (tkState)\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # set numerics
    temporal_scheme = VARS[0].get()
    scheme = VARS[4].get()
    a = VARS[11].get()
    if a == 'cfl': time_step_nature = 'local'
    else: time_step_nature = 'global'
    val = float(VARS[5].get())
    if time_step_nature == 'local': ss_iteration = 1
    else: ss_iteration = 20

    myApp = App.IBM(format='single')
    myApp.set(numb={
    "temporal_scheme": temporal_scheme,
    "ss_iteration": ss_iteration,
    "omp_mode": 1,
    "modulo_verif": 50
    })

    myApp.set(numz={
    "time_step": val,
    "scheme": scheme,
    "time_step_nature": time_step_nature,
    "cfl": val
    })

    nit = VARS[9].get() # nbre d'iterations a faire
    nruns = VARS[15].get() # nbre de runs a faire
    moduloVerif = 50

    # Force recalculation of MX_OMP_SIZE_INT if temporal scheme changed
    import FastC.PyTree as FastC
    FastC.MX_OMP_SIZE_INT = -1

    # open compute
    CTK.t, tc, ts, metrics, graph = myApp.setup('restart.cgns', 'tc.cgns')

    import FastS.PyTree as FastS
    it0 = 0; time0 = 0.
    first = Internal.getNodeFromName1(CTK.t, 'Iteration')
    if first is not None: it0 = Internal.getValue(first)
    first = Internal.getNodeFromName1(CTK.t, 'Time')
    if first is not None: time0 = Internal.getValue(first)
    time_step = Internal.getNodeFromName(CTK.t, 'time_step')
    time_step = Internal.getValue(time_step)
    time_step_nature = Internal.getNodeFromName(CTK.t, 'time_step_nature')
    time_step_nature = Internal.getValue(time_step_nature)
    if time_step_nature == 'local': time_step = 0.

    for it in range(1, nit+1):
        FastS._compute(CTK.t, metrics, it, tc, graph)
        time0 += time_step
        if it%50 == 0:
            CTK.TXT.insert('START', '%d / %d - %f\n'%(it0+it,it0+nit*(nruns-nrun),time0))
            CTK.TXT.update()
        if it%moduloVerif == 0:
            FastS.display_temporal_criteria(CTK.t, metrics, it, format='single', stopAtNan=False)
            #CTK.display(CTK.t)
    NITRUN += 1
    Internal.createUniqueChild(CTK.t, 'Iteration', 'DataArray_t', value=it0+nit)
    Internal.createUniqueChild(CTK.t, 'Time', 'DataArray_t', value=time0)

    # Wall extraction
    import Connector.ToolboxIBM as TIBM
    import Post.IBM as P_IBM
    global WALL
    # extract one BAR wall for each base (body)
    #WALL = TIBM.extractIBMWallFields(tc, tb=BODY) # avec surface
    #WALL = TIBM.extractIBMWallFields(tc) # seulement en node
    #WALL = Internal.getZones(WALL)

    (WALL, CL, CD) = P_IBM.loads(tb_in=BODY, tc_in=tc, alpha=getAlphaAngle(BODY), beta=0., Sref=1.)
    C.convertPyTree2File(WALL, 'walls.cgns')
    C.convertPyTree2File(tc, 'tc_restart.cgns')

    # optional plots
    if CTK.TKPLOTXY is not None: 
        updateWallPlot(WALL)
        updateLoadPlot(CL, CD, NITRUN, nit)
    return None

#===================================================================
# Retourne l'angle d'attaque a partir du refstate (en degres)
#===================================================================
def getAlphaAngle(t):
    refState = Internal.getNodeFromName2(t, 'ReferenceState')
    alpha = 0.
    if refState is not None:
        vx = Internal.getNodeFromName1(refState, 'VelocityX')
        vx = Internal.getValue(vx)
        vy = Internal.getNodeFromName1(refState, 'VelocityY')
        vy = Internal.getValue(vy)
        if abs(vx) < 1.e-12 and abs(vy) < 1.e-12: alpha = 0.
        elif abs(vx) < 1.e-12 and vy > 0: alpha = 90.
        elif abs(vx) < 1.e-12 and vy < 0: alpha = -90.
        else: alpha = math.atan2(vy, vx)*180./math.pi
    return alpha

#========================================================================
# update 1D plots graphs from walls tree
#========================================================================
def updateWallPlot(walls):
    import tkPlotXY
    if not tkPlotXY.IMPORTOK: return
    global DESKTOP1
    # rename zones to keep the same names through computation
    for c, z in enumerate(Internal.getZones(walls)): z[0] = 'wall'+str(c)

    # filter walls following extractWalls tag
    outwalls = []
    wallsz = Internal.getZones(walls)
    for c, b in enumerate(Internal.getBases(BODY)):
        zones = Internal.getZones(b)
        if zones: 
            n = Internal.getNodeFromPath(zones[0], '.Solver#define/extractWalls')
            if n is not None:
                v = Internal.getValue(n)
                if v == 1:
                    wallsz[c][0] = b[0] # set body base name to zone name
                    # remove useless fields
                    Internal._rmNodesFromName(wallsz[c], 'VelocityX')
                    Internal._rmNodesFromName(wallsz[c], 'VelocityY')
                    Internal._rmNodesFromName(wallsz[c], 'VelocityZ')
                    Internal._rmNodesFromName(wallsz[c], 'yplus')
                    Internal._rmNodesFromName(wallsz[c], 'friction*')
                    outwalls.append(wallsz[c])

    if outwalls == []: return

    # create desktop if needed
    if DESKTOP1 is None:
        DESKTOP1 = tkPlotXY.DesktopFrameTK(CTK.WIDGETS['masterWin'])
        DESKTOP1.setData(outwalls)
        graph = DESKTOP1.createGraph('Wall fields', '1:1')
        for z in DESKTOP1.data:
            # print("z", z)
            curve = tkPlotXY.Curve(zone=[z], varx='CoordinateX', vary='Pressure@FlowSolution',
                                   legend_label=z)
            graph.addCurve('1:1', curve)
    else:
        DESKTOP1.setData(outwalls)

#========================================================================
# update 1D plots graphs from CL and CD values
#========================================================================
def updateLoadPlot(CL, CD, nitrun, nit):
    import tkPlotXY
    if not tkPlotXY.IMPORTOK: return
    import Generator.PyTree as G
    global DESKTOP2, LOAD

    # filter CL, CD following extractLoads tag
    outCL = []; outCD = []
    for c, b in enumerate(Internal.getBases(BODY)):
        zones = Internal.getZones(b)
        if zones: 
            n = Internal.getNodeFromPath(zones[0], '.Solver#define/extractLoads')
            if n is not None:
                v = Internal.getValue(n)
                if v == 1:
                    outCL.append(CL[c])
                    outCD.append(CD[c])

    if outCL == []: return

    if DESKTOP2 is None:
        LOAD = []
        for c, v in enumerate(outCL): # par base = component
            z = G.cart((0,0,0), (1,1,1), (2,1,1))
            # X is iteration, Y is CL, Z is CD
            Internal._renameNode(z, 'CoordinateX', 'it')
            Internal._renameNode(z, 'CoordinateY', 'CL')
            Internal._renameNode(z, 'CoordinateZ', 'CD')
            C.setValue(z, 'CL', nitrun-1, outCL[c])
            C.setValue(z, 'CD', nitrun-1, outCD[c])
            C.setValue(z, 'it', nitrun-1, 0)
            C.setValue(z, 'CL', nitrun, outCL[c])
            C.setValue(z, 'CD', nitrun, outCD[c])
            C.setValue(z, 'it', nitrun, nit)
            z[0] = Internal.getBases(BODY)[c][0]
            LOAD.append(z)
        DESKTOP2 = tkPlotXY.DesktopFrameTK(CTK.WIDGETS['masterWin'])
        DESKTOP2.setData(LOAD)
        graph = DESKTOP2.createGraph('Lift/Drag', '1:1')
        for z in DESKTOP2.data:
            curve = tkPlotXY.Curve(zone=[z], varx='it', vary='CL',
                                   legend_label=z)
            graph.addCurve('1:1', curve)
    else:
        z0 = LOAD[0]
        npts0 = C.getNPts(z0)
        if npts0 <= NITRUN:
            for c, v in enumerate(outCL): # par base = component
                z = G.cart((0,0,0), (1,1,1), (npts0+1,1,1))
                # X is iteration, Y is CL, Z is CD
                Internal._renameNode(z, 'CoordinateX', 'it')
                Internal._renameNode(z, 'CoordinateY', 'CL')
                Internal._renameNode(z, 'CoordinateZ', 'CD')
                z0 = LOAD[c]
                z0p = Internal.getNodeFromName2(z0, 'it')[1]
                zp = Internal.getNodeFromName2(z, 'it')[1]
                zp[0:-1] = z0p[:]
                z0p = Internal.getNodeFromName2(z0, 'CL')[1]
                zp = Internal.getNodeFromName2(z, 'CL')[1]
                zp[0:-1] = z0p[:]
                z0p = Internal.getNodeFromName2(z0, 'CD')[1]
                zp = Internal.getNodeFromName2(z, 'CD')[1]
                zp[0:-1] = z0p[:]
                C.setValue(z, 'CL', nitrun, outCL[c])
                C.setValue(z, 'CD', nitrun, outCD[c])
                itprev = C.getValue(z, 'it', nitrun-1)
                C.setValue(z, 'it', nitrun, itprev+nit)
                z[0] = Internal.getBases(BODY)[c][0]
                LOAD[c] = z
        else:
            for c, z in enumerate(LOAD): # par base = component
                C.setValue(z, 'CL', nitrun, outCL[c])
                C.setValue(z, 'CD', nitrun, outCD[c])
                itprev = C.getValue(z, 'it', nitrun-1)
                C.setValue(z, 'it', nitrun, itprev+nit)
        DESKTOP2.setData(LOAD)

#===============================================================
def writeFiles():
    writePrepFile()
    writeComputeFile()
    CTK.TXT.insert('START', 'Write prep.py and compute.py.\n')
    return None

#===============================================================
def changeCflTimeStep(event=None):
    global TIMESTEP, CFL
    ntype = VARS[11].get() # nouveau type?
    v = VARS[5].get()
    # Save previous (not sure)
    if ntype == 'time_step': CFL = v
    else: TIMESTEP = v
    if ntype == 'time_step': VARS[5].set(TIMESTEP)
    else: VARS[5].set(CFL)

#==============================================================================
# Write prep file
#==============================================================================
def writePrepFile():
    if CTK.t == []: return

    mode = VARS[10].get()
    if mode == 'Body': tbody = CTK.t
    elif BODY is not None: tbody = BODY

    # Recupere la base REFINE si presente
    b = Internal.getNodeFromName1(tbody, 'REFINE')
    if b is not None:
        tbox = C.newPyTree()
        tbox[2].append(b)
        tbody = Internal.rmNodesFromName1(tbody, 'REFINE')
    else: tbox = None

    # Save preventif
    C.convertPyTree2File(tbody, 'body.cgns')
    if tbox is not None:
        C.convertPyTree2File(tbox, 'tbox.cgns')

    f = open('prep.py', 'w')

    text= """
import Apps.Fast.IBM as App
myApp = App.IBM(format='single')
myApp.input_var.vmin = 21
myApp.input_var.check = False    
"""

    if tbox is None: text +="myApp.prepare('body.cgns', t_out='t.cgns', tc_out='tc.cgns')"
    else: text += "myApp.prepare('body.cgns', t_out='t.cgns', tbox='tbox.cgns', tc_out='tc.cgns')"

    f.write(text)
    CTK.TXT.insert('START', 'File prep.py written.\n')
    f.close()

#==============================================================================
# Write compute file
#==============================================================================
def writeComputeFile():
    if CTK.t == []: return

    temporal_scheme = VARS[0].get()
    scheme = VARS[4].get()
    a = VARS[11].get()
    if a == 'cfl': time_step_nature = 'local'
    else: time_step_nature = 'global'
    val = float(VARS[5].get())
    if time_step_nature == 'global': time_step = val; cfl = 4.
    else: time_step = 0.1; cfl = val
    if time_step_nature == 'local': ss_iteration = 3
    else: ss_iteration = 30
    nit = VARS[9].get()

    f = open('compute.py', 'w')

    text= """
import Apps.Fast.IBM as App

myApp = App.IBM(format='single')
myApp.set(numb={
    "temporal_scheme": "%s",
    "ss_iteration": %d,
    "omp_mode": 1,
    "modulo_verif": 50
})
myApp.set(numz={
    "time_step": %g,
    "scheme": "%s",
    "time_step_nature": "%s",
    "cfl": %g,
})

# Compute
myApp.compute('t.cgns', 'tc.cgns', t_out='restart.cgns', tc_out='tc_restart.cgns', nit=%d)
"""%(temporal_scheme, ss_iteration, time_step, scheme, time_step_nature, cfl, nit)

    f.write(text)
    CTK.TXT.insert('START', 'File compute.py written.\n')
    f.close()

#==============================================================================
def displaySlices():
    global TKSLICE
    if TKSLICE is None:
        frame = CTK.TKMODULEFRAMES['tkSlice']
        try: TKSLICE = __import__('tkSlice'); TKSLICE.createApp(frame)
        except: TKSLICE = None
    if TKSLICE is not None:
        if WALL is not None: TKSLICE.WALL = WALL
        zones = Internal.getZones(CTK.t)
        CTK.__MAINACTIVEZONES__ = range(len(zones))
        TKSLICE.NODE2CENTER = True
        TKSLICE.VARS[0].set('X')
        TKSLICE.VARS[1].set(TKSLICE.XVALUE)
        TKSLICE.view()
        TKSLICE.VARS[0].set('Y')
        TKSLICE.VARS[1].set(TKSLICE.YVALUE)
        TKSLICE.view()

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkFastSolver  [ + ]  ', font=CTK.FRAMEFONT, 
                           takefocus=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    Frame.columnconfigure(2, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkFastSolver')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- temporal_scheme -
    V = TK.StringVar(win); V.set('explicit'); VARS.append(V)
    # -1- ss_iteration -
    V = TK.StringVar(win); V.set('20'); VARS.append(V)
    # -2- modulo_verif -
    V = TK.StringVar(win); V.set('200'); VARS.append(V)
    # -3- restart_fields -
    V = TK.StringVar(win); V.set('1');VARS.append(V)
    # -4- scheme -
    V = TK.StringVar(win); V.set('roe_min'); VARS.append(V)
    # -5- Time step or cfl -
    V = TK.DoubleVar(win); V.set(TIMESTEP); VARS.append(V)
    # -6- Snear -
    V = TK.DoubleVar(win); V.set(0.01); VARS.append(V)
    # -7- IBC type -
    V = TK.StringVar(win); V.set('Musker'); VARS.append(V)
    # -8- dfar local -
    V = TK.DoubleVar(win); V.set(20.); VARS.append(V)
    # -9- nbre d'iterations -
    V = TK.IntVar(win); V.set(100); VARS.append(V)
    # -10 - body or main mode
    V = TK.StringVar(win); V.set('Body'); VARS.append(V)
    # -11- time_step or cfl
    V = TK.StringVar(win); V.set('time_step'); VARS.append(V)
    # -12- mask inv or not -
    V = TK.StringVar(win); V.set('out'); VARS.append(V)
    # -13- body time
    V = TK.DoubleVar(win); V.set(0.); VARS.append(V)
    # -14- body time steps
    V = TK.DoubleVar(win); V.set(10); VARS.append(V)
    # -15- nbre de runs -
    V = TK.IntVar(win); V.set(1); VARS.append(V)

    #- Mode -
    B = TTK.Label(Frame, text="Mode")
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Currently viewed tree.\nIn body mode:\nOne Base per body\nBase REFINE for refinement zones.')
    # B = TTK.OptionMenu(Frame, VARS[10], 'Body', 'Main', 'PrevStep', command=changeMode)
    # B.grid(row=0, column=1, columnspan=2, sticky=TK.EW)

    B = TTK.Button(Frame, text="Body", command=bodyMode)
    B.grid(row=0, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Display IBC bodies.')
    B = TTK.Button(Frame, text="Main", command=mainTreeMode)
    B.grid(row=0, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Display main tree.')

    # - temporal scheme -
    B = TTK.Label(Frame, text="time_scheme")
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Time integration.')
    B = TTK.OptionMenu(Frame, VARS[0], 'explicit', 'implicit', 'implicit_local')
    B.grid(row=1, column=1, columnspan=2, sticky=TK.EW)

    # - ss_iteration -
    #B = TTK.Label(Frame, text="ss_iteration")
    #B.grid(row=1, column=0, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Nbre de sous iterations max.')
    #B = TTK.Entry(Frame, textvariable=VARS[1])
    #B.grid(row=1, column=1, columnspan=2, sticky=TK.EW)

    # - scheme -
    B = TTK.Label(Frame, text="scheme")
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Numerical scheme.')
    B = TTK.OptionMenu(Frame, VARS[4], 'roe_min', 'ausmpred', 'senseur')
    B.grid(row=2, column=1, columnspan=2, sticky=TK.EW)

    # - time_step -
    B = TTK.OptionMenu(Frame, VARS[11], "time_step", "cfl", command=changeCflTimeStep)
    B.grid(row=3, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Time step.')
    B = TTK.Entry(Frame, textvariable=VARS[5], background='White')
    B.grid(row=3, column=1, columnspan=2, sticky=TK.EW)

    # - Set data -
    B = TTK.Button(Frame, text="Set data", command=setData)
    BB = CTK.infoBulle(parent=B, text='Set data into selected zone.')
    B.grid(row=4, column=0, columnspan=1, sticky=TK.EW)
    B = TTK.Button(Frame, text="Get data", command=getData,
                   image=iconics.PHOTO[8], padx=0, pady=0, compound=TK.RIGHT)
    BB = CTK.infoBulle(parent=B, text='Get data from selected zone.')
    B.grid(row=4, column=1, columnspan=2, sticky=TK.EW)

    # - save & reload -
    B = TTK.Button(Frame, text="Reload", command=reloadPrevStep)
    BB = CTK.infoBulle(parent=B, text='Reload the main tree from the previous step.')
    B.grid(row=5, column=0, sticky=TK.EW)
    B = TTK.Button(Frame, text="Save files", command=writeFiles)
    BB = CTK.infoBulle(parent=B, text='Write python and cgns files to run elsewhere.')
    B.grid(row=5, column=1, columnspan=2, sticky=TK.EW)

    # - compute -
    B = TTK.Button(Frame, text="Compute", command=run)
    WIDGETS['compute'] = B
    BB = CTK.infoBulle(parent=B, text='Launch computation.')
    B.grid(row=6, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[9], width=5, background='White')
    BB = CTK.infoBulle(parent=B, text='Number of iterations for each run.')
    B.grid(row=6, column=1, columnspan=1, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[15], width=5, background='White')
    BB = CTK.infoBulle(parent=B, text='Number of runs.')
    B.grid(row=6, column=2, columnspan=1, sticky=TK.EW)

    # - Body time -
    #B = TTK.Button(Frame, text="Step", command=updateBodyAndPrepare)
    #BB = CTK.infoBulle(parent=B, text='Apply one motion step to body.')
    #B.grid(row=6, column=0, columnspan=1, sticky=TK.EW)
    #B = TTK.Entry(Frame, textvariable=VARS[13], width=4, background="White")
    #B.grid(row=6, column=1, columnspan=1, sticky=TK.EW)
    #B.bind('<Return>', updateBodyAndPrepare)
    #BB = CTK.infoBulle(parent=B, text='Current body time.')
    #B = TTK.Entry(Frame, textvariable=VARS[14], width=4, background="White")
    #B.grid(row=6, column=2, columnspan=1, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Number of steps to reach body time 1.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['SolverNoteBook'].add(WIDGETS['frame'], text='tkFastSolver')
    except: pass
    CTK.WIDGETS['SolverNoteBook'].select(WIDGETS['frame'])
    getData()

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['SolverNoteBook'].hide(WIDGETS['frame'])

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
    (win, menu, file, tools) = CTK.minimal('tkFastSolver '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
