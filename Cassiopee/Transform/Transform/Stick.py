# Stick
import Converter.PyTree as C
import Converter.Internal as Internal
import Transform.PyTree as T
import numpy

#=======================================================================================
# Stick t on tp
# IN: t: volume tree with "FamilySpecified::stick" BCs
# IN: tp: surface projection tree
# Done to stick blade meshes on hubs
# Projection can introduce artifacts
# modifies t
#=======================================================================================
def _stick(t, tp, stickBCName='FamilySpecified:stick', nitSmooth=0):

    import Ael.Quantum as KDG

    # merge tp in a single projection surface
    zones = Internal.getZones(tp)
    surf = []
    for z in zones:
        p = C.convertArray2Hexa(z)
        surf.append(p)
    surf = T.join(surf)

    # Deform "stickBC"
    defo = {}
    zones = Internal.getZones(t)
    for z in zones:
        walls = C.extractBCOfName(t, stickBCName, reorder=False)
        C._initVars(walls, '{hx} = 0.')
        C._initVars(walls, '{hy} = 0.')
        C._initVars(walls, '{hz} = 0.')

        for w in walls:
            if nitSmooth == 0: wp = T.projectOrtho(w, surf)
            else: wp = T.projectOrthoSmooth(w, surf, niter=nitSmooth)
            #wp = T.projectDir(w, surf, (0,1,0))
            hp = Internal.getNodeFromName2(w, 'hx')
            cp = Internal.getNodeFromName2(w, 'CoordinateX')
            cpp = Internal.getNodeFromName2(wp, 'CoordinateX')
            hp[1][:] = cpp[1][:]-cp[1][:]
            hp = Internal.getNodeFromName2(w, 'hy')
            cp = Internal.getNodeFromName2(w, 'CoordinateY')
            cpp = Internal.getNodeFromName2(wp, 'CoordinateY')
            hp[1][:] = cpp[1][:]-cp[1][:]
            hp = Internal.getNodeFromName2(w, 'hz')
            cp = Internal.getNodeFromName2(w, 'CoordinateZ')
            cpp = Internal.getNodeFromName2(wp, 'CoordinateZ')
            hp[1][:] = cpp[1][:]-cp[1][:]

            C._initVars(w, '{xn} = {hx} + {CoordinateX}')
            C._initVars(w, '{yn} = {hy} + {CoordinateY}')
            C._initVars(w, '{zn} = {hz} + {CoordinateZ}')

            name = w[0].replace('/', '#')
            hx = Internal.getNodeFromName2(w, 'hx')[1]
            hy = Internal.getNodeFromName2(w, 'hy')[1]
            hz = Internal.getNodeFromName2(w, 'hz')[1]
            defo[name] = [hx,hy,hz]

    DeformationArgs={"Approach"          :  "Quaternions",
                     "Epsilon"           :  0.15,
                     "Leafsize"          :  4,
                     "OmpAllInOne"       :  True,
                     "Ndivision"         :  100,
                     "NullDisplacements" :  "Weighted",
                     "Smoothing"         :  False }

    defTree = KDG.KeDefGrid(t, **DeformationArgs)
    defTree.set_Amplitude(1.)

    # Imposed boundaries
    for k in defo:
        defTree.setBndSurfTo(k, "imposed", defo[k])

    # Free boundaries
    free = C.extractBCOfName(t, 'FamilySpecified:free', reorder=False)
    for f in free:
        name = f[0].replace('/', '#')
        defTree.setBndSurfTo(name, "free")

    # Fixed boundaries
    fixed = C.extractBCOfName(t, 'FamilySpecified:fixed', reorder=False)
    for f in fixed:
        name = f[0].replace('/', '#')
        defTree.setBndSurfTo(name, "null")

    # Sliding boundaries
    sliding = C.extractBCOfName(t, 'FamilySpecified:sliding', reorder=False)
    for f in sliding:
        name = f[0].replace('/', '#')
        defTree.setBndSurfTo(name, "slidingonsurface")

    defTree.makeSources()
    defTree.computeMeshDisplacement()

    Internal.__FlowSolutionNodes__ = 'Displacement#0'
    C._initVars(t, '{CoordinateX}={CoordinateX}+{DisplacementX}')
    C._initVars(t, '{CoordinateY}={CoordinateY}+{DisplacementY}')
    C._initVars(t, '{CoordinateZ}={CoordinateZ}+{DisplacementZ}')
    Internal.__FlowSolutionNodes__ = 'FlowSolution'

    return None