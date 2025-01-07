"""CFD Post-processing module."""
from . import Post
from . import post
__version__ = Post.__version__

try: range = xrange
except: pass

import math

try:
    import Converter
    import Converter.PyTree as C
    import Converter.Internal as Internal
    import numpy
    import Generator.PyTree as G
except ImportError:
    raise ImportError("Post.PyTree: requires Converter.PyTree and Generator.PyTree modules.")

# Extraction de la liste des pts
def extractPoint(t, Pts, order=2, extrapOrder=1,
                 constraint=40., tol=1.e-6, hook=None, mode='robust'):
    """Extract the solution in one point.
    Usage: extractPoint(t, (x,y,z), order, tol, hook, mode)"""

    if not isinstance(Pts, list):
        a = Converter.array('CoordinateX,CoordinateY,CoordinateZ',1,1,1)
        a[1][0,0] = Pts[0]; a[1][1,0] = Pts[1]; a[1][2,0] = Pts[2]
    else:
        npts = len(Pts)
        a = Converter.array('CoordinateX,CoordinateY,CoordinateZ',npts,1,1)
        for i in range(npts):
            a[1][0,i] = Pts[i][0]
            a[1][1,i] = Pts[i][1]
            a[1][2,i] = Pts[i][2]

    z = C.convertArrays2ZoneNode('extractPt', [a])
    z = C.convertArray2Node(z)
    _extractMesh(t, z, order, extrapOrder, constraint, tol, hook, mode)

    soln = C.getFields(Internal.__FlowSolutionNodes__, z)[0]
    solc = C.getFields(Internal.__FlowSolutionCenters__, z)[0]
    if not isinstance(Pts, list):
        resn = []; resc = []
        solnN = soln[1]
        for nov in range(solnN.shape[0]):
            resn.append(solnN[nov][0])

        if solc != []: # que les noeuds
            solcN = solc[1]
            for nov in range(solcN.shape[0]):
                resc.append(solcN[nov][0])
            return [resn,resc]
        else: return resn

    else:
        allresN = []; allresC = []
        solnN = soln[1]
        npts = solnN.shape[1]
        nvar = solnN.shape[0]
        for nop in range(npts):
            resn = []
            for nov in range(nvar):
                resn.append(solnN[nov][nop])
            allresN.append(resn)

        if solc != []: # que les noeuds
            solcN = solc[1]
            npts = solcN.shape[1]
            nvar = solcN.shape[0]
            for nop in range(npts):
                resc = []
                for nov in range(nvar):
                    resc.append(solcN[nov][nop])
                allresC.append(resc)

            return [allresN,allresC]
        else: return allresN

def extractPlane(t, T, order=2, tol=1.e-6):
    """Slice solution with a plane.
    Usage: extractPlane(t, (coefa, coefb, coefc, coefd), order)"""
    #t = C.deleteFlowSolutions__(t, 'centers')
    t = C.center2Node(t, Internal.__FlowSolutionCenters__)
    A = C.getAllFields(t, 'nodes')
    a = Post.extractPlane(A, T, order, tol)
    return C.convertArrays2ZoneNode('extractedPlane', [a])

def projectCloudSolution(cloud, surf, dim=3, loc='nodes', ibm=False, isPreProjectOrtho=False, old=False):
    """Project the solution defined on a set of points to a TRI surface."""
    surf2 = Internal.copyRef(surf)
    _projectCloudSolution(cloud, surf2, dim=dim, loc=loc, ibm=ibm, isPreProjectOrtho=isPreProjectOrtho, old=old)
    return surf2

def _projectCloudSolution(cloud, surf, dim=3, loc='nodes', ibm=False, isPreProjectOrtho=False, old=False):
    """Project the solution defined on a set of points to a TRI surface."""

    ##old=True for wind tunnel test cases. Has shown to yield more physically accurate results.
    ##old means we are reverting back to predominant extrapolations for the projectCloudSolution.
    ##When a more stable & robust solution is obtained for these test cases this argument will be removed.
    ##See Antoine J. @ DAAA/DEFI for more questions. - error appears at 90 edges of the wind tunnels.

    ##This orthogonal projection does a projection of the cloud solution onto the surface.
    ##Provides better results and is needed when the surface mesh for post processing is different than
    ##that used for the preprocessing (automatic grid generation).
    ##default=False to pass validation database. In due time: default=True
    ##Suggestion proposed by Stephanie P. @ DAAA/DEFI.
    ##See Antoine J. @ DAAA/DEFI for more questions.
    if isPreProjectOrtho:
        print("Project Cloud Solution::performing pre orthogonal projection")
        import Transform.PyTree as T
        cloud = T.projectOrtho(cloud, surf);


    fc = C.getAllFields(cloud, 'nodes')[0]
    zones = Internal.getZones(surf)
    for noz in range(len(zones)):
        interpData = Internal.getNodeFromName(zones[noz], 'POST_MLS')
        if interpData is not None:
            offset = Internal.getNodeFromName(interpData, 'offset')[1]
            interpDonor = Internal.getNodeFromName(interpData, 'interpDonor')[1]
            interpCoef = Internal.getNodeFromName(interpData, 'interpCoef')[1]
            fs = C.getAllFields(zones[noz], 'nodes')[0]
            res = Post.projectCloudSolutionWithInterpData(fc, fs, offset, interpDonor, interpCoef, dim=dim)
            C.setFields([res], zones[noz], 'nodes')
        else:
            fs = C.getAllFields(zones[noz], 'nodes')[0]
            res = Post.projectCloudSolution(fc, fs, dim=dim, ibm=ibm, old=old)
            C.setFields([res], zones[noz], 'nodes')
    return None

def prepareProjectCloudSolution(cloud, surf, dim=3, loc='nodes', ibm=False):
    """Compute the MLS interpolation data for projectCloudSolutionWithInterpData."""
    surf2 = Internal.copyRef(surf)
    _prepareProjectCloudSolution(cloud, surf2, dim=dim, loc=loc, ibm=ibm)
    return surf2

def _prepareProjectCloudSolution(cloud, surf, dim=3, loc='nodes', ibm=False):
    """Compute the MLS interpolation data for projectCloudSolutionWithInterpData."""
    fc = C.getAllFields(cloud, 'nodes')[0]
    zones = Internal.getZones(surf)
    for noz in range(len(zones)):
        fs = C.getAllFields(zones[noz], 'nodes')[0]
        res = Post.prepareProjectCloudSolution(fc, fs, dim=dim, ibm=ibm)

        interpDonor, interpCoef = res
        offset = [0]; cpt = 0
        for listOfDonorPoints in interpDonor:
            cpt += len(listOfDonorPoints)
            offset.append(cpt)
        offset = numpy.array(offset, dtype=Internal.E_NpyInt)
        interpDonor = numpy.concatenate(interpDonor)
        interpCoef = numpy.concatenate(interpCoef)

        children = []
        children.append(Internal.createNode('offset' , 'DataArray_t', offset))
        children.append(Internal.createNode('interpDonor' , 'DataArray_t', interpDonor))
        children.append(Internal.createNode('interpCoef' , 'DataArray_t', interpCoef))
        Internal._createChild(zones[noz], 'POST_MLS', 'UserDefinedData_t', value=None, children=children)

    return None

# hook is a list of pointers on ADT for donor zones of t - created by C.createHook(a,'extractMesh')
def extractMesh(t, extractionMesh, order=2, extrapOrder=1,
                constraint=40., tol=1.e-6, hook=None, mode='robust'):
    """Extract the solution on a given mesh.
    Usage: extractMesh(t, extractMesh, order, extrapOrder, constraint, tol, hook)"""
    te = Internal.copyRef(extractionMesh)
    _extractMesh(t, te, order, extrapOrder, constraint, tol, hook, mode)
    return te

def _extractMesh(t, extractionMesh, order=2, extrapOrder=1,
                 constraint=40., tol=1.e-6, hook=None, mode='robust'):
    """Extract the solution on a given mesh.
    Usage: extractMesh(t, extractMesh, order, extrapOrder, constraint, tol, hook)"""

    # we sort structured then unstructured
    orderedZones=[]
    for i,z in enumerate(Internal.getZones(extractionMesh)):
        if Internal.getZoneType(z)==1: orderedZones.append(i)
    for i,z in enumerate(Internal.getZones(extractionMesh)):
        if Internal.getZoneType(z)==2: orderedZones.append(i)

    if mode == 'robust':
        tc = C.center2Node(t, Internal.__FlowSolutionCenters__)
        C._orderVariables(tc, varsn=[], varsc=[])

        if hook is not None:
            if not isinstance(hook,list): raise TypeError("_extractMesh: hook must be a list of hooks on ADTs.")
        fa = C.getAllFields(tc, 'nodes')
        del tc
        an = C.getFields(Internal.__GridCoordinates__, extractionMesh)
        res = Post.extractMesh(fa, an, order, extrapOrder, constraint, tol, hook)
        if len(res) != len(orderedZones):
            raise ValueError("_extractMesh: invalid number of zones.")
        nor = 0
        for r in res:
            nozorig = orderedZones[nor]
            z = Internal.getZones(extractionMesh)[nozorig]
            C.setFields([r], z, 'nodes', writeDim=False)
            nor += 1

    else: # accurate: extract les centres sur le maillage en centres
        C._orderVariables(t, varsn=[], varsc=[])
        varsC = C.getVarNames(t, excludeXYZ=True, loc='centers')
        varsN = C.getVarNames(t, excludeXYZ=True, loc='nodes')
        if len(varsC) == 0 or len(varsN) == 0: return extractionMesh
        varsC = varsC[0]; varsN = varsN[0]
        zones = Internal.getZones(extractionMesh)
        if len(varsN) != 0:
            an = C.getFields(Internal.__GridCoordinates__, zones, api=1)
            allf = C.getAllFields(t, 'nodes', api=1)
            if allf != []:
                res = Post.extractMesh(allf, an, order, extrapOrder,
                                       constraint, tol)
                if len(res) != len(orderedZones):
                    raise ValueError("_extractMesh: invalid number of zones.")
                nor = 0
                for r in res:
                    nozorig = orderedZones[nor]
                    z = Internal.getZones(extractionMesh)[nozorig]
                    C.setFields([r], z, 'nodes', writeDim=False)
                    nor += 1


        if len(varsC) != 0:
            if hook is not None:
                print("Warning: _extractMesh: hook is not used in 'accurate' mode.")
            tp = Internal.addGhostCells(t, t, 1, adaptBCs=0)
            an = C.getFields(Internal.__GridCoordinates__, zones)
            ac = Converter.node2Center(an)
            fc = C.getFields(Internal.__GridCoordinates__, tp)
            nzones = len(fc)
            for i in range(len(fc)):
                if len(fc[i]) == 4:
                    try: import Transform
                    except: raise ImportError('extractMesh: Transform module is required.')
                    fc[i] = Transform.dual(fc[i], extraPoints=0)
                else: fc[i] = Converter.node2Center(fc[i])
            fa = C.getFields(Internal.__FlowSolutionCenters__, tp)
            del tp
            allf = []
            nzones = len(fc)
            for i in range(nzones):
                if fc[i] != [] and fa[i] != 0:
                    allf.append(Converter.addVars([fc[i], fa[i]]))
                elif fa[i] != []: allf.append(fa[i])
                elif fc[i] != []: allf.append(fc[i])
            if allf != []:
                res = Post.extractMesh(allf, ac, order, extrapOrder,
                                       constraint, tol)
                res = Converter.rmVars(res, ['CoordinateX','CoordinateY','CoordinateZ'])
                if len(res) != len(orderedZones):
                    raise ValueError("_extractMesh: invalid number of zones.")
                nor = 0
                for r in res:
                    nozorig = orderedZones[nor]
                    z = Internal.getZones(extractionMesh)[nozorig]
                    C.setFields([r], z, 'centers', writeDim=False)
                    nor += 1
    return None

def coarsen(t, indicName='indic', argqual=0.1, tol=1.e6):
    """Coarsen a surface TRI-type mesh given a coarsening indicator for each
    element.
    Usage: coarsen(t, indicName, argqual, tol)"""
    tp = Internal.copyRef(t)
    C._deleteZoneBC__(tp)
    C._deleteGridConnectivity__(tp)
    nodes = Internal.getZones(tp)
    for z in nodes:
        taga = C.getFields(Internal.__FlowSolutionCenters__, z)
        taga = Converter.extractVars(taga, [indicName])[0]
        fc = C.getFields(Internal.__GridCoordinates__, z)[0]
        fa = C.getFields(Internal.__FlowSolutionNodes__, z)[0]
        if fc != [] and fa != []:
            f = Converter.addVars([fc, fa])
            fp = Post.coarsen(f, taga, argqual, tol)
            C.setFields([fp], z, 'nodes')
        elif fa != []:
            fp = Post.coarsen(fa, taga, argqual, tol)
            C.setFields([fp], z, 'nodes')
        elif fc != []:
            fp = Post.coarsen(fc, taga, argqual, tol)
            C.setFields([fp], z, 'nodes')
    C._deleteFlowSolutions__(tp, 'centers')
    return tp

def refine(t, indicName='indic', w=-1):
    """Refine a surface TRI-type mesh given a refinement indicator for each
    element or refine with butterfly algorithm.
    Usage: refine(t, indicName, w)"""
    tp = Internal.copyRef(t)
    _refine(tp, indicName, w)
    return tp

def _refine(t, indicName='indic', w=-1):
    """Refine a surface TRI-type mesh given a refinement indicator for each element."""
    C._deleteZoneBC__(t)
    C._deleteGridConnectivity__(t)
    if w < 0: # linear with indicator field
        zones = Internal.getZones(t)
        indicName = indicName.split(':')
        if len(indicName) != 1: indicName = indicName[1]
        else: indicName = indicName[0]
        for z in zones:
            taga = C.getFields(Internal.__FlowSolutionCenters__, z)
            taga = Converter.extractVars(taga, [indicName])[0]
            fc = C.getFields(Internal.__GridCoordinates__, z)[0]
            fa = C.getFields(Internal.__FlowSolutionNodes__, z)[0]
            if fc != [] and fa != []:
                f = Converter.addVars([fc, fa])
                fp = Post.refine(f, taga)
                C.setFields([fp], z, 'nodes')
            elif fa != []:
                fp = Post.refine(fa, taga)
                C.setFields([fp], z, 'nodes')
            elif fc != []:
                fp = Post.refine(fc, taga)
                C.setFields([fp], z, 'nodes')
        C._deleteFlowSolutions__(t, 'centers')
    else: # butterfly everywhere
        C._deleteFlowSolutions__(t, 'centers')
        C._TZA(t, 'nodes', 'nodes', Post.refine, None, None, w)
    return None

#==============================================================================
# Verifie une chaine de variable
# Retourne 0: toutes les variables sont aux noeuds
# Retourne 1: toutes les variables sont aux centres
# Retourne -1: mixte
#==============================================================================
def checkVariables__(vars):
    if isinstance(vars, list):
        return checkVariables___(vars)
    elif isinstance(vars, str):
        variables = vars.split(',')
        return checkVariables___(variables)
    else: return -1

def checkVariables___(vars):
    loc = -1 # unknown
    for i in vars:
        s = i.find('centers:')
        if s != -1: # found
            if loc == -1: loc = 1
            elif loc == 0: return -1
        else:
            if loc == -1: loc = 0
            elif loc == 1: return -1
    return loc

# SelectCells preserving center flow field solutions
# only for tag defined as center variable
def selectCells(t, F, varStrings=[], strict=0, cleanConnectivity=True):
    """Select cells in a given pyTree.
    Usage: selectCells(t, F, varStrings)"""
    tp = Internal.copyRef(t)
    if varStrings == []: # formula
        ret = checkVariables__(F)
    else: ret = checkVariables__(varStrings)

    if ret == 0: # nodes
        C._deleteZoneBC__(tp)
        C._deleteGridConnectivity__(tp)

        zones = Internal.getZones(tp)
        for z in zones:

            # Get ParentElement node if any
            GEl = Internal.getElementNodes(z)
            NGON = 0; found = False; PE = None
            for c in GEl:
                if c[1][0] == 22: found = True; break
                NGON += 1
            if found:
                node = GEl[NGON]
                PE = Internal.getNodeFromName1(node, 'ParentElements')

            fc = C.getFields(Internal.__GridCoordinates__, z)[0]
            fa = C.getFields(Internal.__FlowSolutionNodes__, z)[0]
            fb = C.getFields(Internal.__FlowSolutionCenters__, z)[0]

            if fc != [] and fa != []:
                ft = Converter.addVars([fc, fa])

                if fb != []:
                    if PE is not None:
                        (PE2,fp,fq) = Post.selectCells(ft, F, fb, varStrings, strict, PE[1], cleanConnectivity)
                    else:
                        (fp,fq) = Post.selectCells(ft, F, fb, varStrings, strict, None, cleanConnectivity)

                    C._deleteFlowSolutions__(z, 'centers')
                    C.setFields([fq], z, 'centers')

                else:
                    if PE is not None:
                        (PE2,fp) = Post.selectCells(ft, F, fb, varStrings, strict, PE[1], cleanConnectivity)
                    else:
                        fp = Post.selectCells(ft, F, fb, varStrings, strict, None, cleanConnectivity)

                    C._deleteFlowSolutions__(z, 'centers')

                C.setFields([fp], z, 'nodes')

            elif fa != []:
                if fb != []:
                    if PE is not None:
                        (PE2,fp,fq) = Post.selectCells(fa, F, fb, varStrings, strict, PE[1], cleanConnectivity)
                    else:
                        (fp,fq) = Post.selectCells(fa, F, fb, varStrings, strict, None, cleanConnectivity)

                    C.setFields([fq], z, 'centers')
                else:
                    if PE is not None:
                        (PE2,fp) = Post.selectCells(fa, F, fb, varStrings, strict, PE[1], cleanConnectivity)
                    else:
                        fp = Post.selectCells(fa, F, fb, varStrings, strict, None, cleanConnectivity)

                C.setFields([fp], z, 'nodes')

            elif fc != []:
                if fb != []:
                    if PE is not None:
                        (PE2,fp,fq) = Post.selectCells(fc, F, fb, varStrings, strict, PE[1], cleanConnectivity)
                    else:
                        (fp,fq) = Post.selectCells(fc, F, fb, varStrings, strict, None, cleanConnectivity)

                    C.setFields([fq], z, 'centers')
                else:
                    if PE is not None:
                        (PE2, fp) = Post.selectCells(fc, F, fb, varStrings, strict, PE[1], cleanConnectivity)
                    else:
                        fp = Post.selectCells(fc, F, fb, varStrings, strict, None, cleanConnectivity)

                C.setFields([fp], z, 'nodes')

            # Set ParentElement
            if PE is not None:
                GEl = Internal.getElementNodes(z)
                NGON = 0; found = False
                for c in GEl:
                    if c[1][0] == 22: found = True; break
                    NGON += 1
                if found: Internal.createUniqueChild(GEl[NGON], 'ParentElements', 'DataArray_t', value=PE2)

        if Internal.isTopTree(tp): C._deleteEmptyZones(tp)
        return tp

    elif ret == 1: # centers
        C._deleteZoneBC__(tp)
        C._deleteGridConnectivity__(tp)
        C._addVars(tp, 'centers:__tag__')
        fields = C.getFields(Internal.__FlowSolutionCenters__, tp)

        tags = []
        if varStrings == []:
            for f in fields:
                tag = Post.buildTag2__(f, F)
                tags.append(tag)
        else:
            for f in fields:
                tag = Post.buildTag1__(f, F, varStrings)
                tags.append(tag)
        C.setFields(tags, tp, 'centers', False)
        tp = selectCells2(tp, 'centers:__tag__')

        C._rmVars(tp, 'centers:__tag__')
        if Internal.isTopTree(tp): C._deleteEmptyZones(tp)
        return tp
    else:
        raise ValueError("selectCells: variables are not co-localized.")


def selectCells2(t, tagName, strict=0, cleanConnectivity=True):
    """Select cells in a given array.
    Usage: selectCells2(t, tagName)"""
    import KCore
    tp = Internal.copyRef(t)
    C._deleteZoneBC__(tp)
    C._deleteGridConnectivity__(tp)
    res = tagName.split(':')
    if len(res) == 2 and res[0] == 'centers': loc = 1
    else: loc = 0

    zones = Internal.getZones(tp)
    for z in zones:

        # Get ParentElement node if any
        # =============================
        GEl = Internal.getElementNodes(z)
        NGON = 0; found = False; PE = None
        for c in GEl:
            if c[1][0] == 22: found = True; break
            NGON += 1
        if found:
            node = GEl[NGON]
            PE = Internal.getNodeFromName1(node, 'ParentElements')

        if loc == 0: # noeuds
            fb   = C.getFields(Internal.__FlowSolutionCenters__, z)[0]
            taga = C.getFields(Internal.__FlowSolutionNodes__, z)
            taga = Converter.extractVars(taga, [tagName])[0]
        else:
            fb   = C.getFields(Internal.__FlowSolutionCenters__, z)[0]
            taga = Converter.extractVars(fb, [res[1]])

        fc = C.getFields(Internal.__GridCoordinates__, z)[0]
        fa = C.getFields(Internal.__FlowSolutionNodes__, z)[0]

        if loc != 0: # centres
            if KCore.isNamePresent(fb,res[1]) > -1:
                if fb[1].shape[0] == 1: fb = None
                else: fb = Converter.rmVars(fb, res[1])

        if fc != [] and fa != []:
            f = Converter.addVars([fc, fa])

            if fb != [] and fb is not None: # il y a des champs en centres
                if PE is not None:
                    (PE2, fp,fq) = Post.selectCells2(f, taga, fb, strict, loc, PE[1], cleanConnectivity)
                else:
                    (fp,fq) = Post.selectCells2(f, taga, fb, strict, loc, None, cleanConnectivity)

                C._deleteFlowSolutions__(z, 'centers')
                C.setFields([fq], z, 'centers')

            else:  # pas de champ en centres
                if PE is not None:
                    (PE2, fp) = Post.selectCells2(f, taga, [], strict, loc, PE[1], cleanConnectivity)
                else:
                    fp = Post.selectCells2(f, taga, [], strict, loc, None, cleanConnectivity)

                Internal._rmNodesFromName(z,Internal.__FlowSolutionCenters__)

            C.setFields([fp], z, 'nodes')

        elif fa != []:
            if fb != [] and fb is not None: # il y a des champs en centres
                if PE is not None:
                    (PE2, fp,fq) = Post.selectCells2(fa, taga, fb, strict, loc, PE[1], cleanConnectivity)
                else:
                    (fp,fq) = Post.selectCells2(fa, taga, fb, strict, loc, None, cleanConnectivity)
                C.setFields([fq], z, 'centers')
            else:        # pas de champ en centres
                if PE is not None:
                    (PE2, fp) = Post.selectCells2(fa, taga, [], strict, loc, PE[1], cleanConnectivity)
                else:
                    fp = Post.selectCells2(fa, taga, [], strict, loc, None, cleanConnectivity)
                Internal._rmNodesFromName(z,Internal.__FlowSolutionCenters__)
            C.setFields([fp], z, 'nodes')

        elif fc != []:
            if fb != [] and fb is not None: # il y a des champs en centres
                if PE is not None:
                    (PE2, fp,fq) = Post.selectCells2(fc, taga, fb, strict, loc, PE[1], cleanConnectivity)
                else:
                    (fp,fq) = Post.selectCells2(fc, taga, fb, strict, loc, None, cleanConnectivity)

                C.setFields([fq], z, 'centers')
            else:        # pas de champ en centres

                if PE is not None:
                    (PE2, fp) = Post.selectCells2(fc, taga, [], strict, loc, PE[1], cleanConnectivity)
                else:
                    fp = Post.selectCells2(fc, taga, [], strict, loc, None, cleanConnectivity)
                Internal._rmNodesFromName(z,Internal.__FlowSolutionCenters__)
            C.setFields([fp], z, 'nodes')

        # Set ParentElement
        if PE is not None:
            GEl = Internal.getElementNodes(z)
            NGON = 0; found = False
            for c in GEl:
                if c[1][0] == 22: found = True; break
                NGON += 1
            if found: Internal.createUniqueChild(GEl[NGON], 'ParentElements', 'DataArray_t', value=PE2)

    if Internal.isTopTree(tp): C._deleteEmptyZones(tp)

    return tp

def selectCells3(t, tagName):
    try: import Transform.PyTree as T
    except: raise ImportError('selectCells: Transform module is required.')
    tp = Internal.copyRef(t)
    C._deleteZoneBC__(tp)
    C._deleteGridConnectivity__(tp)
    tpp, type = Internal.node2PyTree(tp)
    res = tagName.split(':')
    if len(res) == 2: name = res[1]
    else: name = tagName
    bases = Internal.getBases(tpp)
    for b in bases:
        c = 0
        for z in b[2]:
            if z[3] == 'Zone_t':
                centers = Internal.getNodeFromName(z, Internal.__FlowSolutionCenters__)
                if centers is not None:
                    taga = Internal.getNodeFromName(centers, name)
                    if taga is not None:
                        ret = post.selectCells3(taga[1], 1)
                        ret = T.subzone(z, ret, type='elements')
                        b[2][c] = ret
            c += 1
    tp = Internal.pyTree2Node(tpp, type)
    return tp

def frontFaces(t, tagName):
    """Select faces located in the front of tag=0 and tag=1.
    Usage: frontFaces(t, tagName)"""
    tp = Internal.copyRef(t)
    C._deleteZoneBC__(tp)
    C._deleteGridConnectivity__(tp)
    C._deleteFlowSolutions__(tp, 'centers')
    nodes = Internal.getZones(tp)
    for z in nodes:
        taga = C.getFields(Internal.__FlowSolutionNodes__, z)
        taga = Converter.extractVars(taga, [tagName])[0]
        fc = C.getFields(Internal.__GridCoordinates__, z)[0]
        fa = C.getFields(Internal.__FlowSolutionNodes__, z)[0]
        if fc != [] and fa != []:
            f = Converter.addVars([fc, fa])
            fp = Post.frontFaces(f, taga)
            C.setFields([fp], z, 'nodes')
        elif fa != []:
            fp = Post.frontFaces(fa, taga)
            C.setFields([fp], z, 'nodes')
        elif fc != []:
            fp = Post.frontFaces(fc, taga)
            C.setFields([fp], z, 'nodes')
    return tp

def interiorFaces(t, strict=0):
    """Interior faces of an array. The argument strict equal to 1 means
    that interior faces with only interior nodes are taken into account.
    Usage: interiorFaces( t, strict)"""
    t = C.deleteFlowSolutions__(t, 'centers')
    t = C.TZA1(t, 'nodes', 'nodes', True, Post.interiorFaces, strict)
    bases = Internal.getBases(t)
    for b in bases: b[1][0] = 1
    return t

def exteriorFacesStructured(t):
    """Return the list of exterior faces for a structured mesh
    Usage: exteriorFacesStructured(a)"""
    t = C.deleteFlowSolutions__(t, 'centers')
    zones = Internal.getZones(t)
    listzones = []
    for z in zones:
        field = C.getAllFields(z, 'nodes')[0]
        A = Post.exteriorFacesStructured(field)
        c = 1
        for a in A:
            name = z[0]+str(c)
            listzones.append(C.convertArrays2ZoneNode(name, [a]))
            c += 1
    return listzones

def exteriorFaces(t, indices=None):
    """Exterior faces of a mesh.
    Usage: exteriorFaces(t,indices)"""
    tp = Internal.copyRef(t)
    _exteriorFaces(tp, indices)
    return tp

def _exteriorFaces(t, indices=None):
    C._deleteZoneBC__(t)
    C._deleteFlowSolutions__(t, 'centers')
    return C._TZA1(t, 'nodes', 'nodes', True, Post.exteriorFaces, indices)

def exteriorElts(t):
    """Exterior (border) elts of a mesh.
    Usage: exteriorElts(t)"""
    tp = Internal.copyRef(t)
    _exteriorElts(tp)
    return tp

def _exteriorElts(t):
    """Exterior (border) elts of a mesh."""
    C._deleteFlowSolutions__(t, 'centers')
    return C._TZA1(t, 'nodes', 'nodes', True, Post.exteriorElts)

def exteriorEltsStructured(t, depth=1):
    """Exterior (border) elts of a mesh as a structured grid.
    Usage: exteriorEltsStructured(t, depth)"""
    tp = Internal.copyRef(t)
    _exteriorEltsStructured(tp, depth)
    return tp

def _exteriorEltsStructured(t, depth=1):
    """Exterior (border) elts of a mesh as a structured grid."""
    C._TZA2(t, 'nodes', 'nodes', True, Post.exteriorEltsStructured, depth)
    C._TZA2(t, 'centers', 'centers', False, Post.exteriorEltsStructured, depth)
    return None

def computeVariables(t, varList,
                     gamma=-1., rgp=-1., s0=0., betas=-1.,
                     Cs=-1., mus=-1., Ts=-1.):
    """Compute variables defined in varList.
    Usage: computeVariables(array, varList, gamma=1.4, rgp=287.053, s0=0., betas=1.458e-6, Cs=110.4, mus=0., Ts=0.)"""
    tp = Internal.copyRef(t)
    try:
        state = C.getState(t)
        if gamma < 0: gamma = state[11]
        cvInf = state[7]
        if rgp < 0: rgp = (gamma-1.)*cvInf
        if Cs < 0: Cs = state[16]
        if mus < 0: mus = state[15]
        if Ts < 0: Ts = state[17]
        if betas < 0: betas = mus*(Ts+Cs)/(Ts**(3/2))
    except: pass
    if gamma < 0: gamma = 1.4
    if rgp < 0: rgp=287.053
    if betas < 0: betas=1.458e-6
    if Cs < 0: Cs=110.4
    if mus < 0: mus=1.76e-5
    if Ts < 0: Ts=273.15

    if isinstance(varList, str): varList = [varList]
    varnamesn = [] # variables aux noeuds
    varnamesc = [] # variables aux centres
    for var in varList:
        v = var.split(':')
        if len(v) > 1:
            if v[0] == 'centers': name = v[1]; varnamesc.append(name)
            else: name = v[1]; varnamesn.append(name)
        else: varnamesn.append(var)
    presn = C.isNamePresent(t, 'Density')
    presc = C.isNamePresent(t, 'centers:Density')
    if presc == -1 and varnamesc != []:
        tp = C.node2Center(tp, ['Density', 'MomentumX', 'MomentumY',
                                'MomentumZ', 'EnergyStagnationDensity'])
    if presn == -1 and varnamesn != []:
        tp = C.center2Node(tp, ['centers:Density', 'centers:MomentumX',
                                'centers:MomentumY', 'centers:MomentumZ',
                                'centers:EnergyStagnationDensity'])

    tp = C.TZAGC(tp, 'both', 'both', False,
                 Post.computeVariables, Post.computeVariables,
                 varnamesn, gamma, rgp, s0, betas, Cs, mus, Ts,
                 varnamesc, gamma, rgp, s0, betas, Cs, mus, Ts)
    if presc == -1 and varnamesc != []:
        tp = C.rmVars(tp, ['centers:Density', 'centers:MomentumX',
                           'centers:MomentumY', 'centers:MomentumZ',
                           'centers:EnergyStagnationDensity'])
    if presn == -1 and varnamesn != []:
        tp = C.rmVars(tp, ['Density', 'MomentumX', 'MomentumY',
                           'MomentumZ', 'EnergyStagnationDensity'])
    return tp

def _computeVariables(t, varList,
                      gamma=-1., rgp=-1., s0=0., betas=-1.,
                      Cs=-1., mus=-1., Ts=-1.):
    """Compute variables defined in varList.
    Usage: computeVariables(t, varList, gamma=1.4, rgp=287.053, s0=0., betas=1.458e-6, Cs=110.4, mus=0., Ts=0.)"""
    try:
        state = C.getState(t)
        if gamma < 0: gamma = state[11]
        cvInf = state[7]
        if rgp < 0: rgp = (gamma-1.)*cvInf
        if Cs < 0: Cs = state[16]
        if mus < 0: mus = state[15]
        if Ts < 0: Ts = state[17]
        if betas < 0: betas = mus*(Ts+Cs)/(Ts**(3/2))
    except: pass
    if gamma < 0: gamma = 1.4
    if rgp < 0: rgp=287.053
    if betas < 0: betas=1.458e-6
    if Cs < 0: Cs=110.4
    if mus < 0: mus=1.76e-5
    if Ts < 0: Ts=273.15
    if isinstance(varList, str): varList = [varList]

    varnamesn = [] # variables aux noeuds
    varnamesc = [] # variables aux centres
    for var in varList:
        v = var.split(':')
        if len(v) > 1:
            if v[0] == 'centers': name = v[1]; varnamesc.append(name)
            else: name = v[1]; varnamesn.append(name)
        else: varnamesn.append(var)
    presn = C.isNamePresent(t, 'Density')
    presc = C.isNamePresent(t, 'centers:Density')
    if presc == -1 and varnamesc != []:
        raise ValueError("computeVariables: conservative variables missing.")
    if presn == -1 and varnamesn != []:
        raise ValueError("computeVariables: conservative variables missing.")

    C._TZAGC(t, 'both', 'both', False,
             Post.computeVariables, Post.computeVariables,
             varnamesn, gamma, rgp, s0, betas, Cs, mus, Ts,
             varnamesc, gamma, rgp, s0, betas, Cs, mus, Ts)
    if presc == -1 and varnamesc != []:
        C._rmVars(t, ['centers:Density', 'centers:MomentumX',
                      'centers:MomentumY', 'centers:MomentumZ',
                      'centers:EnergyStagnationDensity'])
    if presn == -1 and varnamesn != []:
        C._rmVars(t, ['Density', 'MomentumX', 'MomentumY',
                      'MomentumZ', 'EnergyStagnationDensity'])
    return None

def computeVariables2(t, varList, gamma=-1., rgp=-1., s0=0., betas=-1.,
                      Cs=-1., mus=-1., Ts=-1.):
    """Compute variables (in place version) defined in varList.
    Usage: computeVariables2(array, varList, gamma=1.4, rgp=287.053, s0=0., betas=1.458e-6, Cs=110.4, mus=0., Ts=0.)"""
    tp = Internal.copyRef(t)
    _computeVariables2(tp, varList, gamma, rgp, s0, betas, Cs, mus, Ts)
    return tp

def _computeVariables2(t, varList, gamma=-1., rgp=-1., s0=0., betas=-1.,
                       Cs=-1., mus=-1., Ts=-1.):
    if gamma < 0:
        try: gamma = C.getState(t, 'Gamma')
        except: pass
    if gamma < 0: gamma = 1.4
    if rgp < 0:
        try:
            cv  = C.getState(t, 'Cv')
            rgp = cv*(gamma-1.)
        except: pass
    if rgp   < 0: rgp = 287.053
    if s0 < 0:
        try: s0 = C.getState(t, 's0')
        except: pass
    if betas < 0:
        try: betas = C.getState(t, 'betas')
        except: pass
    if betas < 0: betas = 1.458e-6
    if Cs < 0:
        try: Cs = C.getState(t, 'Cs')
        except: pass
    if mus < 0:
        try: mus = C.getState(t, 'Mus')
        except: pass
    if Ts < 0:
        try: Ts = C.getState(t, 'Ts')
        except: pass
    if Cs >= 0. and mus >= 0. and Ts >= 0.:
        betas = mus*(1.+Cs/Ts)/math.sqrt(Ts)

    if Cs  < 0: Cs  = 110.4
    if mus < 0: mus = 1.76e-5
    if Ts  < 0: Ts  = 273.15

    varnamesn = [] # variables aux noeuds
    varnamesc = [] # variables aux centres
    for var in varList:
        v = var.split(':')
        if len(v) > 1:
            if v[0] == 'centers':
                name = v[1]
                varnamesc.append(name)
            else:
                name = v[1]
                varnamesn.append(name)
        else: varnamesn.append(var)

    if varnamesn != []:
        C.__TZC2(t, 'nodes', False,
                 Post._computeVariables2, varnamesn, gamma, rgp, s0, betas, Cs, mus, Ts)

    if varnamesc != []:
        C.__TZC2(t, 'centers', False,
                 Post._computeVariables2, varnamesc, gamma, rgp, s0, betas, Cs, mus, Ts)

    return None

def _computeVariablesBC(t, varList, gamma=1.4, rgp=287.053, s0=0.,
                        betas=1.458e-6, Cs=110.4, mus=1.76e-5, Ts=273.15):
    zones = Internal.getZones(t)
    if zones == []: zones = [t] # must be a BC node
    for z in zones:
        bcs = Internal.getNodesFromType2(z, 'BC_t')
        for b in bcs:
            datas = Internal.getBCDataSet(z, b)
            fields = []; connects = []
            for d in datas: # build array
                np = d[1].size; ne = np-1
                dim = ['Unstructured', np, ne, 'NODE', 1]
                f = Internal.convertDataNode2Array(d, dim, connects)
                fields.append(f[1])
            if fields != []:
                fields = Converter.addVars(fields)
                fn = Post.computeVariables(fields, varList, gamma, rgp, s0, betas, Cs, mus, Ts)
                nofld = fn[1].shape[0]-1
                for varName in varList:
                    varName = varName.split('=')[0]
                    varName = varName.replace('{', '')
                    varName = varName.replace('}', '')
                    varName = varName.replace('centers:', '')
                    varName = varName.replace('nodes:', '')
                    varName = varName.replace(' ', '')
                    f = Converter.extractVars(fn, [varName])
                    fieldFaceNode = Internal.createDataNode(varName, f, 0, cellDim=1)
                    cont, c = Internal.getParentOfNode(z, datas[0])
                    Internal._createUniqueChild(cont, varName, 'DataArray_t', value=fieldFaceNode[1])
    return None

#===============================================================================
# INPUT: t: tree of skin/wall borders (velocity gradients must be defined yet)
#===============================================================================
def computeWallShearStress(t):
    """Compute wall shear stress."""
    from . import extraVariablesPT
    tp = Internal.copyRef(t)
    extraVariablesPT._computeWallShearStress(tp)
    return tp

def _computeWallShearStress(t):
    """Compute wall shear stress."""
    from . import extraVariablesPT
    return extraVariablesPT._computeWallShearStress(t)

def computeExtraVariable(t, varname, gamma=-1, rgp=-1.,
                         Cs=-1., mus=-1., Ts=-1.):
    """Compute variables that requires a change of location."""
    from . import extraVariablesPT
    try:
        state = C.getState(t)
        if gamma < 0: gamma = state[11]
        cvInf = state[7]
        if rgp < 0: rgp = (gamma-1.)*cvInf
        if Cs < 0: Cs = state[16]
        if mus < 0: mus = state[15]
        if Ts < 0: Ts = state[17]
    except: pass
    if gamma < 0: gamma = 1.4
    if rgp < 0: rgp=287.053
    if Cs < 0: Cs=110.4
    if mus < 0: mus=1.76e-5
    if Ts < 0: Ts=273.15

    if varname == 'Vorticity' or varname == 'nodes:Vorticity':
        t2 = extraVariablesPT.computeVorticity(t)
        zvars = ['centers:VorticityX', 'centers:VorticityY', 'centers:VorticityZ']
        t2 = C.center2Node(t2, zvars)
        C._rmVars(t2, zvars)
        return t2
    elif varname == 'centers:Vorticity':
        return extraVariablesPT.computeVorticity(t)
    elif varname == 'VorticityMagnitude' or varname == 'nodes:VorticityMagnitude':
        t2 = extraVariablesPT.computeVorticityMagnitude(t)
        t2 = C.center2Node(t2, 'centers:VorticityMagnitude')
        C._rmVars(t2, 'centers:VorticityMagnitude')
        return t2
    elif varname == 'centers:VorticityMagnitude':
        return extraVariablesPT.computeVorticityMagnitude(t)
    elif varname == 'QCriterion' or varname == 'nodes:QCriterion':
        t2 = extraVariablesPT.computeQCriterion(t)
        t2 = C.center2Node(t2, 'centers:QCriterion')
        C._rmVars(t2, 'centers:QCriterion')
        return t2
    elif varname == 'centers:QCriterion':
        return extraVariablesPT.computeQCriterion(t)
    elif varname == 'ShearStress' or varname == 'nodes:ShearStress':
        t2 = extraVariablesPT.computeShearStress(t, gamma, rgp, Cs, mus, Ts)
        zvars = ['centers:ShearStressXX',
                 'centers:ShearStressXY',
                 'centers:ShearStressXZ',
                 'centers:ShearStressYX',
                 'centers:ShearStressYY',
                 'centers:ShearStressYZ',
                 'centers:ShearStressZX',
                 'centers:ShearStressZY',
                 'centers:ShearStressZZ']
        t2 = C.center2Node(t2, zvars)
        C._rmVars(t2, zvars)
        return t2
    elif varname == 'centers:ShearStress':
        return extraVariablesPT.computeShearStress(t, gamma, rgp, Cs, mus, Ts)
    elif varname == 'SkinFriction' or varname == 'nodes:SkinFriction':
        return extraVariablesPT.computeSkinFriction(t, centers=0, tangent=0)
    elif varname == 'centers:SkinFriction':
        return extraVariablesPT.computeSkinFriction(t, centers=1, tangent=0)
    elif varname == 'SkinFrictionTangential' or varname == 'nodes:SkinFrictionTangential':
        return extraVariablesPT.computeSkinFriction(t, centers=0, tangent=1)
    elif varname == 'centers:SkinFrictionTangential':
        return extraVariablesPT.computeSkinFriction(t, centers=1, tangent=1)
    else:
        print('Warning: computeExtraVariable: unknown variable: %s.'%varname)

#==============================================================================
# Importe les variables de t1 dans t2, retourne t2 modifie
# si method=0, utilise les noms de zones pour les reconnaitre
# si method=1, utilise les coordonnees (a eps pres)
# si method=2, utilise le meme ordre des zones dans les deux arbres
# addExtra=1, ajoute les zones non reconnues dans une base EXTRA
#==============================================================================
def importVariables(t1, t2, method=0, eps=1.e-6, addExtra=1):
    """Import the variables of tree t1 to tree t2."""
    a2 = Internal.copyRef(t2)
    zones1 = Internal.getZones(t1)
    zones2 = Internal.getZones(a2)
    nzones1 = len(zones1); nzones2 = len(zones2)
    dejaVu = numpy.zeros(nzones2, dtype=Internal.E_NpyInt)

    locDict={}
    if method == 2:
        for noz1 in range(nzones1):
            z1 = zones1[noz1]; found = 0
            noz2 = noz1
            if dejaVu[noz2] == 0 and found == 0:
                z2 = zones2[noz2]
                (parent2, d) = Internal.getParentOfNode(a2, z2)
                loc = identifyZonesLoc__(z1, z2, method, eps)
                locDict[z1[0]]=loc
                if abs(loc) == 1: # nodes
                    dejaVu[noz2] = noz1+1; found = 1
                    sol1 = Internal.getNodesFromType1(z1, 'FlowSolution_t')
                    dim1 = Internal.getZoneDim(z1)
                    cn1 = Internal.getElementNodes(z1)
                    for s1 in sol1:
                        loc1 = 'nodes'
                        loci = Internal.getNodesFromType1(s1, 'GridLocation_t')
                        if len(loci) > 0:
                            v = Internal.getValue(loci[0])
                            if v == 'CellCenter': loc1 = 'centers'
                        sol11 = Internal.getNodesFromType(s1, 'DataArray_t')
                        for s11 in sol11:
                            ar1 = Internal.convertDataNode2Array(s11, dim1, cn1)[1]
                            z2 = C.setFields([ar1], z2, loc1)
                    parent2[2][d] = z2

                elif abs(loc) == 2: # centers
                    dejaVu[noz2] = noz1+1; found = 1
                    sol1 = Internal.getNodesFromType1(z1, 'FlowSolution_t')
                    dim1 = Internal.getZoneDim(z1)
                    cn1 = Internal.getNodeFromName1(z1, 'GridElements')
                    for s1 in sol1:
                        loc1 = 0
                        loci = Internal.getNodesFromType1(s1, 'GridLocation_t')
                        if len(loci) > 0:
                            v = loci[0]
                            if v[0]!='V': loc1 = 1# 'Vertex'

                        if loc1==0:
                            sol11 = Internal.getNodesFromType(s1, 'DataArray_t')
                            for s11 in sol11:
                                ar1 = Internal.convertDataNode2Array(s11, dim1, cn1)[1]
                                z2 = C.setFields([ar1], z2, 'centers')
                    parent2[2][d] = z2

    elif method in (0,1):
        for noz1 in range(nzones1):
            z1 = zones1[noz1]; found = 0
            for noz2 in range(nzones2):
                if dejaVu[noz2] == 0 and found == 0:
                    z2 = zones2[noz2]
                    (parent2, d) = Internal.getParentOfNode(a2, z2)
                    loc = identifyZonesLoc__(z1, z2, method, eps)
                    locDict[z1[0]]=loc
                    if loc == 1: # nodes
                        dejaVu[noz2] = noz1+1; found = 1
                        sol1 = Internal.getNodesFromType1(z1, 'FlowSolution_t')
                        dim1 = Internal.getZoneDim(z1)
                        cn1 = Internal.getElementNodes(z1)

                        for s1 in sol1:
                            loc1 = 'nodes'
                            loci = Internal.getNodesFromType1(s1, 'GridLocation_t')
                            if len(loci) > 0:
                                v = Internal.getValue(loci[0])
                                if v == 'CellCenter': loc1 = 'centers'

                            sol11 = Internal.getNodesFromType(s1, 'DataArray_t')
                            for s11 in sol11:
                                ar1 = Internal.convertDataNode2Array(s11, dim1, cn1)[1]
                                z2 = C.setFields([ar1], z2, loc1)
                        parent2[2][d] = z2

                    elif loc == 2: # centers
                        dejaVu[noz2] = noz1+1; found = 1
                        sol1 = Internal.getNodesFromType1(z1, 'FlowSolution_t')
                        dim1 = Internal.getZoneDim(z1)
                        cn1 = Internal.getNodeFromName1(z1, 'GridElements')
                        for s1 in sol1:
                            loc1 = 0
                            loci = Internal.getNodesFromType1(s1, 'GridLocation_t')
                            if len(loci) > 0:
                                v = loci[0]
                                if v[0] != 'V': loc1 = 1 # 'Vertex'

                            if loc1 == 0:
                                sol11 = Internal.getNodesFromType(s1, 'DataArray_t')
                                for s11 in sol11:
                                    ar1 = Internal.convertDataNode2Array(s11, dim1, cn1)[1]
                                    z2 = C.setFields([ar1], z2, 'centers')

                        parent2[2][d] = z2
    else: raise NotImplementedError("Method {0!r} is not implemented. Please refer to the documentation".format(method))

    if not Internal.isTopTree(a2): return a2

    tag = numpy.zeros(nzones1, dtype=Internal.E_NpyInt)
    for noz2 in range(nzones2):
        if dejaVu[noz2] > 0: tag[dejaVu[noz2]-1] = 1

    extra = False
    for noz1 in range(nzones1):
        if tag[noz1] == 0: extra = True; break

    if extra and addExtra == 1:
        print('Warning: importVariables: extra grid(s) in t2 detected, added to EXTRA base.')
        C._addBase2PyTree(a2, 'EXTRA', 3)
        base = Internal.getNodesFromName1(a2, 'EXTRA')
        for noz1 in range(nzones1):
            if tag[noz1] == 0:
                z = zones1[noz1]
                loc = locDict[z[0]]
                if abs(loc)==1: # ajout direct
                    base[0][2].append(z)

                elif abs(loc) == 2: # ajout coord en noeud et champ direct (qui correspond aux centres du  nouveau)
                    C._rmVars(z,Internal.__FlowSolutionCenters__)
                    zc = C.center2Node(z)
                    coords = Internal.getNodesFromType1(zc, 'GridCoordinates_t')
                    dim = Internal.getZoneDim(zc)
                    cn = Internal.getElementNodes(zc)
                    for x in coords:
                        ax = Internal.getNodesFromType(x, 'DataArray_t')
                        for sx in ax:
                            ar = Internal.convertDataNode2Array(sx, dim, cn)[1]
                            z = C.setFields([ar], z, 'nodes')
                    fields = Internal.getNodesFromType1(z, 'FlowSolution_t')
                    for x in fields:
                        gloc= Internal.getNodeFromType(x,'GridLocation_t')
                        vloc = 1
                        if gloc is None: vloc = 0
                        else:
                            if gloc[1][0]=='V': vloc=0# =='Vertex'
                        if vloc == 0:
                            ax = Internal.getNodesFromType(x, 'DataArray_t')
                            for sx in ax:
                                ar = Internal.convertDataNode2Array(sx, dim, cn)[1]
                                z = C.setFields([ar], z, 'centers')
                    C._rmVars(z, Internal.__FlowSolutionNodes__)
                    base[0][2].append(z)
    return a2

#==============================================================================
# Retourne 0 si zones differentes,
#          1 si z1 et z2 correspondent a la meme zone en noeuds
#          2 si z1 correspond aux centres de z2
#         -1 si z1 et z2 correspondent par les dimensions en noeuds
#         -2 si z1 et z2 correspondent par les dimensions en centres
# method=0: identification par les noms des zones
# method=1: identification par les sommets des bounding boxes
#==============================================================================
def identifyZonesLoc__(z1, z2, method, eps):
    dim2 = Internal.getZoneDim(z2)
    dim1 = Internal.getZoneDim(z1)

    if method == 0:
        if dim1 == dim2:
            if z1[0] == z2[0]: return 1
            else: return -1
        else:
            if (dim1[0] == 'Structured' and dim2[0] == 'Structured'):
                # Pour les grilles non structurees, il n'est pas possible
                # de stocker les coord en centres actuellement
                ni1 = dim1[1]; nj1 = dim1[2]; nk1 = dim1[3]
                ni2 = dim2[1]; nj2 = dim2[2]; nk2 = dim2[3]
                ni21 = max(ni2-1,1); nj21 = max(nj2-1,1); nk21 = max(nk2-1,1)
                if (ni1 == ni21 and nj1 == nj21 and nk1 == nk21):
                    if z1[0] == z2[0]: return 2
                    else: return -2
            return 0

    # method=1, loc = noeuds ?
    if dim1 == dim2:
        if C.isNamePresent(z1, 'CoordinateX') == -1: return -1
        if C.isNamePresent(z2, 'CoordinateX') == -1: return -1
        (xmin1,ymin1,zmin1) = C.getMinValue(z1, 'GridCoordinates')
        (xmax1,ymax1,zmax1) = C.getMaxValue(z1, 'GridCoordinates')
        (xmin2,ymin2,zmin2) = C.getMinValue(z2, 'GridCoordinates')
        (xmax2,ymax2,zmax2) = C.getMaxValue(z2, 'GridCoordinates')
        dx1 = abs(xmin1-xmin2); dx2 = abs(xmax1-xmax2)
        dy1 = abs(ymin1-ymin2); dy2 = abs(ymax1-ymax2)
        dz1 = abs(zmin1-zmin2); dz2 = abs(zmax1-zmax2)

        # Identification des zones a partir des 8 sommets
        if (dx1 < eps and dx2 < eps and dy1 < eps and
                dy2 < eps and dz1 < eps and dz2 < eps):
            return 1
    # loc = centres ?
    else:
        if dim1[0] == 'Structured' and dim2[0] == 'Structured':
            # Pour les grilles non structurees, il n'est pas possible
            # de stocker les coord en centres actuellement
            ni1 = dim1[1]; nj1 = dim1[2]; nk1 = dim1[3]
            ni2 = dim2[1]; nj2 = dim2[2]; nk2 = dim2[3]
            ni21 = max(ni2-1,1); nj21 = max(nj2-1,1); nk21 = max(nk2-1,1)
            if (ni1 == ni21 and nj1 == nj21 and nk1 == nk21):
                if C.isNamePresent(z1, 'CoordinateX') == -1: return -2
                if C.isNamePresent(z2, 'CoordinateX') == -1: return -2
                z2c = C.node2Center(z2)
                (xmin1,ymin1,zmin1) = C.getMinValue(z1, 'GridCoordinates')
                (xmax1,ymax1,zmax1) = C.getMaxValue(z1, 'GridCoordinates')
                (xmin2,ymin2,zmin2) = C.getMinValue(z2c, 'GridCoordinates')
                (xmax2,ymax2,zmax2) = C.getMaxValue(z2c, 'GridCoordinates')
                dx1 = abs(xmin1-xmin2); dx2 = abs(xmax1-xmax2)
                dy1 = abs(ymin1-ymin2); dy2 = abs(ymax1-ymax2)
                dz1 = abs(zmin1-zmin2); dz2 = abs(zmax1-zmax2)

                # Identification des zones a partir des 8 sommets
                if (dx1 < eps and dx2 < eps and dy1 < eps and
                        dy2 < eps and dz1 < eps and dz2 < eps):
                    return 2
    return 0

#==============================================================================
def zipper(t, options=[]):
    """Extract Chimera surface as an unique unstructured surface.
    Usage: zipper(t, options)"""
    t = C.deleteFlowSolutions__(t,'centers')
    arrays = C.getAllFields(t, 'nodes')
    a = Post.zipper(arrays, options)
    return C.convertArrays2ZoneNode('zipper', [a])

def extractArraysForScalarInteg__(t, var=''):
    zones = Internal.getZones(t)
    zvars = var.split(':')
    if len(zvars) == 2: loc = 'centers'
    else: loc = 'nodes'

    coords = C.getFields(Internal.__GridCoordinates__, zones)
    fields = []; fieldsc = []
    if var == '':
        fields = C.getFields(Internal.__FlowSolutionNodes__, zones)
        fields = Converter.addVars([coords,fields])
        fieldsc = C.getFields(Internal.__FlowSolutionCenters__, zones)
    elif var == Internal.__GridCoordinates__:
        fields = coords
    elif var == Internal.__FlowSolutionNodes__:
        fields = C.getFields(Internal.__FlowSolutionNodes__, zones)
    elif var == Internal.__FlowSolutionCenters__:
        fieldsc = C.getFields(Internal.__FlowSolutionCenters__, zones)
    else: # un champ specifique
        if loc == 'nodes': fields = C.getField(var, zones)
        else: fieldsc = C.getField(var,zones)

    # mise a jour de fields [[],[],[],..,[]] -> []
    foundn = 0; foundc = 0
    for nof in range(len(fields)):
        if fields[nof] != []: foundn = 1
    if foundn == 0: fields = []
    for nof in range(len(fieldsc)):
        if fieldsc[nof] != []: foundc = 1
    if foundc == 0: fieldsc = []
    return [coords, fields, fieldsc]

def extractArraysForVectorInteg__(t, vector):
    zones = Internal.getZones(t)
    loc = 'unknown'
    if len(vector) != 3: raise ValueError("extractArraysForVectorInteg: vector must be of size 3.")
    zvars = []
    for var in vector:
        v = var.split(':')
        if len(v) == 1 and loc != 'centers': loc = 'nodes'; zvars.append(v[0])
        elif len(v)  > 1 and loc != 'nodes': loc = 'centers'; zvars.append(v[1])
        else: raise ValueError("extractArraysForVectorInteg: all the components of the vector must have the same location.")

    coords = C.getFields(Internal.__GridCoordinates__, zones)
    fields = []; fieldsc = []
    if vector == ['CoordinateX','CoordinateY','CoordinateZ']: fields = coords
    else:
        if loc == 'nodes':
            fields = C.getFields(Internal.__FlowSolutionNodes__, zones)
            if fields != []: fields = Converter.extractVars(fields, zvars)
        else:
            fieldsc = C.getFields(Internal.__FlowSolutionCenters__, zones)
            if fieldsc != []: fieldsc = Converter.extractVars(fieldsc, zvars)

    # mise a jour de fields [[],[],[],..,[]] -> []
    foundn = 0; foundc = 0
    for nof in range(len(fields)):
        if fields[nof] != []: foundn = 1
    if foundn == 0: fields = []
    for nof in range(len(fieldsc)):
        if fieldsc[nof] != []: foundc = 1
    if foundc == 0: fieldsc = []
    return [coords, fields, fieldsc]

def extractRatioForInteg__(t):
    zones = Internal.getZones(t)
    # extraction des ratios
    ration = C.getField('ratio', zones)
    ratioc = C.getField('centers:ratio', zones)
    foundn = 0; foundc = 0
    for nof in range(len(ration)):
        if ration[nof] != []: foundn = 1
    if foundn == 0: ration = []
    for nof in range(len(ratioc)):
        if ratioc[nof] != []: foundc = 1
    if foundn == 0: ration = []
    if foundc == 0: ratioc = []
    return [ration, ratioc]

def integ2(t, var=''):
    """Integral of fields defined in t.
    Usage: integ(t, var)"""
    info = extractArraysForScalarInteg__(t, var)
    infor = extractRatioForInteg__(t)
    coords = info[0]; fieldsn = info[1]; fieldsc = info[2]
    ration = infor[0]; ratioc = infor[1]

    resn = []; resc = []
    if fieldsn != []: resn = Post.integ(coords, fieldsn, ration)
    if fieldsc != []: resc = Post.integ(coords, fieldsc, ratioc)
    if resn != [] and resc != []: return resn+resc
    elif resn != []: return resn
    elif resc != []: return resc
    else: return []

def integ(t, var=''):
    """Integral of fields defined in t."""
    try: import Generator.PyTree as G
    except ImportError: raise ImportError("integ: requires Generator module.")
    zvars = var.split(':')
    if len(zvars) == 2:
        if zvars[0] == 'centers': loc = 1; varName = zvars[1]
        else: loc = 0; varName = zvars[1]
    else: loc = 0; varName = var
    ret = 0.
    for z in Internal.getZones(t):
        vol = None; removeVol = False; removeVar = False
        cont = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
        if cont is not None:
            vol = Internal.getNodeFromName1(cont, 'vol')
            ratio = Internal.getNodeFromName1(cont, 'ratio')
        if vol is None:
            removeVol = True
            G._getVolumeMap(z)
        if loc == 0:
            try:
                z2 = C.node2Center(z, varName); removeVar = True
                ret += post.integ2(z2, varName, Internal.__GridCoordinates__,
                                   Internal.__FlowSolutionNodes__, Internal.__FlowSolutionCenters__)
            except: pass
        else:
            ret += post.integ2(z, varName, Internal.__GridCoordinates__,
                               Internal.__FlowSolutionNodes__, Internal.__FlowSolutionCenters__)
        if removeVol: C._rmVars(z, 'centers:vol')
        if removeVar: C._rmVars(z, 'centers:'+varName)
    return [ret] # for compatibility

def integNorm(t, var=''):
    """Integral of fields times normal.
    Usage: integNorm(t, var)"""
    info = extractArraysForScalarInteg__(t, var)
    infor = extractRatioForInteg__(t)
    coords = info[0]; fieldsn = info[1]; fieldsc = info[2]
    ration = infor[0]; ratioc = infor[1]

    resn = []; resc = []
    if fieldsn != []: resn = Post.integNorm(coords, fieldsn, ration)
    if fieldsc != []: resc = Post.integNorm(coords, fieldsc, ratioc)
    if resn != [] and resc != []: return resn+resc
    elif resn != []: return resn
    elif resc != []: return resc
    return [[0.,0.,0.]]

def integNormProduct(t, vector=[]):
    """Integral of scalar product fields times normal.
    Usage: integNormProduct(t, vector)"""
    info = extractArraysForVectorInteg__(t, vector)
    infor = extractRatioForInteg__(t)
    coords = info[0]; fieldsn = info[1]; fieldsc = info[2]
    ration = infor[0]; ratioc = infor[1]
    resn = []; resc = []
    if fieldsn != []: resn = Post.integNormProduct(coords, fieldsn, ration)
    if fieldsc != []: resc = Post.integNormProduct(coords, fieldsc, ratioc)
    if resn != [] and resc != []: return resn+resc
    elif resn != []: return resn
    elif resc != []: return resc
    return []

def integMoment(t, center=(0.,0.,0.), vector=[]):
    """Integral of moments.
    Usage: integMoment(t, (xc, yc, zc), vector)"""
    info = extractArraysForVectorInteg__(t, vector)
    infor = extractRatioForInteg__(t)
    coords = info[0]; fieldsn = info[1]; fieldsc = info[2]
    ration = infor[0]; ratioc = infor[1]

    resn = []; resc = []
    if fieldsn != []: resn = Post.integMoment(coords, fieldsn, ration, center)
    if fieldsc != []: resc = Post.integMoment(coords, fieldsc, ratioc, center)
    if resn != [] and resc != []: return resn+resc
    elif resn != []: return resn
    elif resc != []: return resc
    return [0.,0.,0.]

def integMomentNorm(t, center=(0.,0.,0.), var=''):
    """Integral of moments (OM^f.vect(n)).
    Usage: integMomentNorm( t, (xc, yc, zc), var )"""
    info = extractArraysForScalarInteg__(t, var)
    infor = extractRatioForInteg__(t)
    coords = info[0]; fieldsn = info[1]; fieldsc = info[2]
    ration = infor[0]; ratioc = infor[1]

    resn = []; resc = []
    if fieldsn != []: resn = Post.integMomentNorm(coords, fieldsn, ration, center)
    if fieldsc != []: resc = Post.integMomentNorm(coords, fieldsc, ratioc, center)
    if resn != [] and resc != []: return resn+resc
    elif resn != []: return resn
    elif resc != []: return resc
    return [0.,0.,0.]

def usurp(t):
    """Extract unique surfaces using ranked polygons.
    Usage: usurp(t)"""
    a = Internal.copyRef(t)
    _usurp(a)
    return a

def _usurp(t):
    """Extract unique surfaces using ranked polygons.
    Usage: usurp(t)"""
    coords = C.getFields(Internal.__GridCoordinates__, t)
    fields = C.getFields(Internal.__FlowSolutionCenters__, t)
    r = Post.usurp(coords, fields)
    C.setFields(r, t, 'centers')
    return None

def computeGrad(t, var):
    """Compute the gradient of a variable defined in array.
    Usage: computeGrad(t,var) """
    tp = Internal.copyRef(t)
    gradVar = var # variable on which gradient is computed
    v = var.split(':')
    posv = 0
    if len(v) > 1:
        if v[0] == 'centers':
            posv = C.isNamePresent(tp, v[1])
            tp = C.center2Node(tp, var)
        gradVar = v[1]

    nodes = C.getAllFields(tp, 'nodes')
    if posv == -1: C._rmVars(tp, gradVar)
    centers = Post.computeGrad(nodes, gradVar)
    C.setFields(centers, tp, 'centers')
    return tp

def computeGradLSQ(t, fldNames, parRun=0, fcenters=None, ptlists=None,
                   rfields=None):
    tp = Internal.copyRef(t)
    _computeGradLSQ(tp, fldNames, parRun, fcenters, ptlists, rfields)
    return tp

def _computeGradLSQ(t, fldNames, parRun=0, fcenters=None, ptlists=None,
                    rfields=None):
    fc = None
    if parRun == 0:
        fc, fa = G.getFaceCentersAndAreas(t)
        centers = G.getCellCenters(t, fc, fa)

    zones = Internal.getZones(t)

    for i in range(len(zones)):
        zone = zones[i]

        arr = C.getFields(Internal.__GridCoordinates__, zone, api=3)[0]
        if arr == None: continue

        fsolc = Internal.getNodeFromName1(zone,
                                          Internal.__FlowSolutionCenters__)
        if fsolc == None: raise ValueError("FlowSolutionCenters not found.")

        flds = []
        for fldName in fldNames:
            fsol = Internal.getNodeFromName1(fsolc, fldName)
            if fsol == None:
                raise ValueError('Field ' + fldName + ' not found.')
            flds.append(fsol[1])

        pe = Internal.getNodeFromName(zone, 'ParentElements')[1]

        rflds = None
        if parRun == 0:
            cc = centers[i]
            Grads = post.computeGradLSQ(arr, flds, pe, cc[0], cc[1], cc[2], fc[i],
                                        None, rflds)
        else:
            rflds = rfields[i]
            cx = Internal.getNodeFromName1(fsolc, 'CCx')[1]
            cy = Internal.getNodeFromName1(fsolc, 'CCy')[1]
            cz = Internal.getNodeFromName1(fsolc, 'CCz')[1]
            Grads = post.computeGradLSQ(arr, flds, pe, cx, cy, cz, fcenters[0],
                                        ptlists, rflds)

        for j in range(len(fldNames)-3):
            Grad = Grads[j]
            Internal.createNode('grad' + fldNames[j] + 'x', 'DataArray_t',
                                Grad[0], None, fsolc)
            Internal.createNode('grad' + fldNames[j] + 'y', 'DataArray_t',
                                Grad[1], None, fsolc)
            Internal.createNode('grad' + fldNames[j] + 'z', 'DataArray_t',
                                Grad[2], None, fsolc)

    return None

def computeGrad2(t, var, ghostCells=False, withCellN=True, withTNC=False):
    """Compute the gradient of a variable defined in array.
    Usage: computeGrad2(t, var)"""
    tp = Internal.copyRef(t)
    _computeGrad2(tp, var, ghostCells, withCellN, withTNC)
    return tp

def _computeGrad2(t, var, ghostCells=False, withCellN=True, withTNC=False):
    """Compute the gradient of a variable defined in array.
    Usage: computeGrad2(t, var)"""

    if type(var) == list:
        raise ValueError("computeGrad2: not available for lists of variables.")
    vare = var.split(':')
    if len(vare) > 1: vare = vare[1]

    # Test if field exist
    solc = C.getFields(Internal.__FlowSolutionCenters__, t)[0]
    if solc == []:
        raise ValueError("_computeGrad2: no field detected (check container).")

    # Compute fields on BCMatch (for all match connectivities)
    if not ghostCells:
        allMatch = C.extractAllBCMatch(t, vare)
        if withTNC:
            allMatchTNC = C.extractAllBCMatchTNC(t,vare)
        else:
            allMatchTNC = {}
    else:
        allMatch    = {}
        allMatchTNC = {}

    zones = Internal.getZones(t)
    for z in zones:

        # Test if vol present
        cont = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
        if cont is not None:
            vol  = Internal.getNodeFromName1(cont, 'vol')
            if vol is not None: vol = vol[1]

        # Test if cellN present
        if withCellN:
            if cont is not None:
                cellN  = Internal.getNodeFromName1(cont, 'cellN')
                if cellN is not None: cellN = cellN[1]
        else: cellN = None

        f = C.getField(var, z)[0]
        x = C.getFields(Internal.__GridCoordinates__, z)[0]
        # Get BCDataSet if any
        indices=None; BCField=None

        isghost = Internal.getNodeFromType1(z, 'Rind_t')
        if isghost is None or not ghostCells: # not a ghost cells zone : add BCDataSet
            zoneBC = Internal.getNodesFromType1(z, 'ZoneBC_t')
            if zoneBC is not None:
                BCs = Internal.getNodesFromType1(zoneBC, 'BC_t')
                for b in BCs:
                    datas = Internal.getBCDataSet(z, b)
                    inds = Internal.getBCFaceNode(z, b)
                    if datas != [] and inds != []:
                        bcf = None
                        for i in datas:
                            if i[0] == vare: bcf = i; break
                        if bcf is not None:
                            indsp = inds[1].ravel(order='K')
                            bcfp = bcf[1].ravel(order='K')
                            if indices is None: indices = indsp
                            else: indices = numpy.concatenate((indices, indsp))
                            if BCField is None: BCField = bcfp
                            else: BCField = numpy.concatenate((BCField, bcfp))

        # compute field on BCMatch for current zone
        if allMatch != {}:
            indFace, fldFace = C.computeBCMatchField(z, allMatch, vare)
            # --------------------------------------------------------------
            if fldFace is not None:
                fldp = None
                for fgc in fldFace:
                    fgc = fgc[1][0]
                    if fldp is None: fldp = fgc
                    else: fldp = numpy.concatenate((fldp,fgc))

                indp    = indFace.ravel(order='K')
                fldp    = fldp.ravel(order='K')

                if indices is None: indices = indp
                else: indices = numpy.concatenate((indices, indp))

                if BCField is None: BCField = fldp
                else: BCField = numpy.concatenate((BCField, fldp))
            # --------------------------------------------------------------

        # compute field on TNC match for current zone
        if allMatchTNC != {}:
            indFaceTNC, fldFaceTNC = C.computeBCMatchField(z, allMatchTNC, vare)
            # --------------------------------------------------------------
            if fldFaceTNC is not None:
                fldp = None
                for fgc in fldFaceTNC:
                    fgc = fgc[1][0]
                    if fldp is None: fldp = fgc
                    else: fldp = numpy.concatenate((fldp, fgc))

                indp = indFaceTNC.ravel(order='K')
                fldp = fldp.ravel(order='K')

                if indices is None: indices = indp
                else: indices = numpy.concatenate((indices, indp))

                if BCField is None: BCField = fldp
                else: BCField = numpy.concatenate((BCField, fldp))
            # --------------------------------------------------------------

        if f != []:
            centers = Post.computeGrad2(x, f, vol, cellN, indices=indices, BCField=BCField)
            C.setFields([centers], z, 'centers')

    return None

def computeNormGrad(t, var):
    """Compute the norm of gradient of a variable defined in array.
    Usage: computeNormGrad(t) """
    tp = Internal.copyRef(t)
    gradVar = var             # variable on which gradient is computed
    v = var.split(':')
    posv = 0
    if len(v) > 1:
        if v[0] == 'centers':
            posv = C.isNamePresent(tp, v[1])
            tp = C.center2Node(tp, var)
            gradVar = v[1]
    nodes = C.getAllFields(tp, 'nodes')
    if posv == -1: C._rmVars(tp, gradVar)
    centers = Post.computeNormGrad(nodes, gradVar)
    C.setFields(centers, tp, 'centers')
    return tp

def computeDiv(t,var):
    """Compute the divergence of a variable defined in array.
    Usage: computeDiv(t, var) """
    if isinstance(var, list):
        raise ValueError("computeDiv: not available for lists of variables.")
    tp = Internal.copyRef(t)
    sdirlist = ['X', 'Y', 'Z'] # The order is important!
    dim = 3
    posv = [0]*dim
    divVector = list()
    for i in range(dim):
        name = var+sdirlist[i]
        v = name.split(':')
        if len(v) > 1:
            if v[0] == 'centers':
                posv[i] = C.isNamePresent(tp, v[1])
                tp = C.center2Node(tp, name)
            name = v[1]
        divVector.append(name)

    tn = C.getAllFields(tp, 'nodes') # < get all fields (?) We only need a few...
    for i in range(dim):
        if posv[i] == -1: C._rmVars(tp, divVector[i])
    tc = Post.computeDiv(tn, divVector)
    C.setFields(tc, tp, 'centers')
    return tp

def isArrayEmpty(arr):
    """Return whether an input python list or numpy.ndarray is empty"""
    if isinstance(arr, list):
        return all(map(isArrayEmpty, arr))
    elif isinstance(arr, numpy.ndarray):
        return arr.size == 0
    return False

def computeDiv2(t, var, ghostCells=False, withTNC=False, rmVar=False):
    """Compute the divergence of a variable defined in array."""
    tp = Internal.copyRef(t)
    _computeDiv2(tp, var, ghostCells, withTNC, rmVar=rmVar)
    return tp

def _computeDiv2(t, var, ghostCells=False, withTNC=False, rmVar=False):
    """Compute the divergence at the cell centers of a series of variables."""

    if not isinstance(var, list): var = [var]
    nvars = len(var)
    vare = [v.split(':')[-1] for v in var]
    # Compute fields on BCMatch (for all match connectivities)
    zones = Internal.getZones(t)
    ndims = Internal.getZoneDim(zones[0])[-1]
    if ndims == 1:
        raise ValueError("computeDiv2: not available for 1-dimensional elements.")

    sdirlist = ['X', 'Y', 'Z'] # The order is important!
    varList  = ' '.join('{0}X {0}Y {0}Z'.format(v) for v in vare).split()

    if not ghostCells:
        allMatch = C.extractAllBCMatch(t,varList)
        if withTNC: allMatchTNC = C.extractAllBCMatchTNC(t,varList)
        else: allMatchTNC = {}
    else:
        allMatch    = {}; allMatchTNC = {}

    for z in zones:

        # Test if vol present
        cont = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
        vol  = Internal.getNodeFromName1(cont, 'vol')
        if vol is not None: vol = vol[1]

        cellN  = Internal.getNodeFromName1(cont, 'cellN')
        if cellN is not None: cellN = cellN[1]

        flist, f = [], []
        for v in var:
            for sdir in sdirlist:
                flist.append(C.getField('%s%s' % (v,sdir), z)[0])
        flist = [x for x in flist if x]

        if flist and len(flist)%ndims == 0:
            f = [Converter.addVars(flist[i:i+ndims]) for i in range(0, len(flist), ndims)]

        x = C.getFields(Internal.__GridCoordinates__, z)[0]
        # Get BCDataSet if any
        indices  = None
        BCFieldX = None; BCFieldY = None; BCFieldZ = None

        isghost = Internal.getNodeFromType1(z, 'Rind_t')
        if isghost is None: # not a ghost cells zone : add BCDataSet
            zoneBC = Internal.getNodesFromType1(z, 'ZoneBC_t')
            if zoneBC is not None:
                BCs = Internal.getNodesFromType1(zoneBC, 'BC_t')
                for b in BCs:
                    datas = Internal.getBCDataSet(z, b)
                    inds = Internal.getBCFaceNode(z, b)

                    if datas and inds:
                        bcfpx = [[] for _ in range(nvars)]
                        bcfpy = [[] for _ in range(nvars)]
                        bcfpz = [[] for _ in range(nvars)]
                        for v in range(3*nvars):
                            j = v//3
                            k = v%3
                            for i in datas:
                                if varList[v] == i[0]:
                                    if k == 0: bcfpx[j] = i[1].ravel(order='K')
                                    elif k == 1: bcfpy[j] = i[1].ravel(order='K')
                                    else: bcfpz[j] = i[1].ravel(order='K')
                                    break

                        if not (isArrayEmpty(bcfpx) and isArrayEmpty(bcfpy)):
                            indsp = inds[1].ravel(order='K')
                            if indices is None: indices = indsp
                            else: indices = numpy.concatenate((indices, indsp))

                        # BCFieldX
                        if BCFieldX is None: BCFieldX = bcfpx
                        else: BCFieldX = numpy.concatenate((BCFieldX, bcfpx), axis=1)
                        # BCFieldY
                        if BCFieldY is None: BCFieldY = bcfpy
                        else: BCFieldY = numpy.concatenate((BCFieldY, bcfpy), axis=1)
                        # BCFieldZ
                        if BCFieldZ is None: BCFieldZ = bcfpz
                        else: BCFieldZ = numpy.concatenate((BCFieldZ, bcfpz), axis=1)

        # Compute field on BCMatch for current zone
        if allMatch:
            indFace, fldFace = C.computeBCMatchField(z, allMatch, varList)

            if fldFace is not None:

                fldX = None; fldY = None; fldZ = None
                foundVar = fldFace[0][0].split(',')

                if foundVar and (len(foundVar)/nvars)%3 == 0: # 3D
                    for fgc in fldFace:
                        fgcX = numpy.asarray(fgc[1][0::3], order='K')
                        fgcY = numpy.asarray(fgc[1][1::3], order='K')
                        fgcZ = numpy.asarray(fgc[1][2::3], order='K')

                        if fldX is None: fldX = fgcX
                        else: fldX = numpy.concatenate((fldX,fgcX), axis=1)

                        if fldY is None: fldY = fgcY
                        else: fldY = numpy.concatenate((fldY,fgcY), axis=1)

                        if fldZ is None: fldZ = fgcZ
                        else: fldZ = numpy.concatenate((fldZ,fgcZ), axis=1)

                    indp    = indFace.ravel(order='K')
                    fldX    = fldX.ravel(order='K').reshape(nvars,-1)
                    fldY    = fldY.ravel(order='K').reshape(nvars,-1)
                    fldZ    = fldZ.ravel(order='K').reshape(nvars,-1)

                    if indices is None: indices = indp
                    else: indices = numpy.concatenate((indices, indp))

                    if BCFieldX is None: BCFieldX = fldX
                    else: BCFieldX = numpy.concatenate((BCFieldX, fldX), axis=1)
                    if BCFieldY is None: BCFieldY = fldY
                    else: BCFieldY = numpy.concatenate((BCFieldY, fldY), axis=1)
                    if BCFieldZ is None: BCFieldZ = fldZ
                    else: BCFieldZ = numpy.concatenate((BCFieldZ, fldZ), axis=1)

                elif foundVar and (len(foundVar)/nvars)%2 == 0: # 2D
                    # Config (XY)
                    if foundVar[0][-1] == 'X' and foundVar[1][-1] == 'Y':
                        for fgc in fldFace:
                            fgcX = numpy.asarray(fgc[1][0::2], order='K')
                            fgcY = numpy.asarray(fgc[1][1::2], order='K')

                            if fldX is None: fldX = fgcX
                            else: fldX = numpy.concatenate((fldX,fgcX), axis=1)

                            if fldY is None: fldY = fgcY
                            else: fldY = numpy.concatenate((fldY,fgcY), axis=1)

                        indp    = indFace.ravel(order='K')
                        fldX    = fldX.ravel(order='K').reshape(nvars,-1)
                        fldY    = fldY.ravel(order='K').reshape(nvars,-1)

                        if indices is None: indices = indp
                        else: indices = numpy.concatenate((indices, indp))

                        if BCFieldX is None: BCFieldX = fldX
                        else: BCFieldX = numpy.concatenate((BCFieldX, fldX), axis=1)
                        if BCFieldY is None: BCFieldY = fldY
                        else: BCFieldY = numpy.concatenate((BCFieldY, fldY), axis=1)

                    # Config (XZ)
                    if foundVar[0][-1] == 'X' and foundVar[1][-1] == 'Z':
                        for fgc in fldFace:
                            fgcX = numpy.asarray(fgc[1][0::2], order='K')
                            fgcZ = numpy.asarray(fgc[1][1::2], order='K')

                            if fldX is None: fldX = fgcX
                            else: fldX = numpy.concatenate((fldX,fgcX), axis=1)

                            if fldZ is None: fldZ = fgcZ
                            else: fldZ = numpy.concatenate((fldZ,fgcZ), axis=1)

                        indp    = indFace.ravel(order='K')
                        fldX    = fldX.ravel(order='K').reshape(nvars,-1)
                        fldZ    = fldZ.ravel(order='K').reshape(nvars,-1)

                        if indices is None: indices = indp
                        else: indices = numpy.concatenate((indices, indp))

                        if BCFieldX is None: BCFieldX = fldX
                        else: BCFieldX = numpy.concatenate((BCFieldX, fldX), axis=1)
                        if BCFieldZ is None: BCFieldZ = fldZ
                        else: BCFieldZ = numpy.concatenate((BCFieldZ, fldZ), axis=1)

                    # Config (YZ)
                    if foundVar[0][-1] == 'Y' and foundVar[1][-1] == 'Z':
                        for fgc in fldFace:
                            fgcY = numpy.asarray(fgc[1][0::2], order='K')
                            fgcZ = numpy.asarray(fgc[1][1::2], order='K')

                            if fldY is None: fldY = fgcY
                            else: fldY = numpy.concatenate((fldY,fgcY), axis=1)

                            if fldZ is None: fldZ = fgcZ
                            else: fldZ = numpy.concatenate((fldZ,fgcZ), axis=1)

                        indp    = indFace.ravel(order='K')
                        fldY    = fldY.ravel(order='K').reshape(nvars,-1)
                        fldZ    = fldZ.ravel(order='K').reshape(nvars,-1)

                        if indices is None: indices = indp
                        else: indices = numpy.concatenate((indices, indp))

                        if BCFieldY is None: BCFieldY = fldY
                        else: BCFieldY = numpy.concatenate((BCFieldY, fldY), axis=1)
                        if BCFieldZ is None: BCFieldZ = fldZ
                        else: BCFieldZ = numpy.concatenate((BCFieldZ, fldZ), axis=1)

        # =================================================================================
        # Compute field on TNC match for current zone
        if allMatchTNC:
            indFaceTNC, fldFaceTNC = C.computeBCMatchField(z, allMatchTNC, varList)

            if fldFaceTNC is not None:

                fldX = None; fldY = None; fldZ = None
                foundVar = fldFaceTNC[0][0].split(',')

                if foundVar and (len(foundVar)/nvars)%3 == 0: # 3D
                    for fgc in fldFaceTNC:
                        fgcX = numpy.asarray(fgc[1][0::3], order='K')
                        fgcY = numpy.asarray(fgc[1][1::3], order='K')
                        fgcZ = numpy.asarray(fgc[1][2::3], order='K')

                        if fldX is None: fldX = fgcX
                        else: fldX = numpy.concatenate((fldX,fgcX), axis=1)

                        if fldY is None: fldY = fgcY
                        else: fldY = numpy.concatenate((fldY,fgcY), axis=1)

                        if fldZ is None: fldZ = fgcZ
                        else: fldZ = numpy.concatenate((fldZ,fgcZ), axis=1)

                    indp    = indFaceTNC.ravel(order='K')
                    fldX    = fldX.ravel(order='K').reshape(nvars,-1)
                    fldY    = fldY.ravel(order='K').reshape(nvars,-1)
                    fldZ    = fldZ.ravel(order='K').reshape(nvars,-1)

                    if indices is None: indices = indp
                    else: indices = numpy.concatenate((indices, indp))

                    if BCFieldX is None: BCFieldX = fldX
                    else: BCFieldX = numpy.concatenate((BCFieldX, fldX), axis=1)

                    if BCFieldY is None: BCFieldY = fldY
                    else: BCFieldY = numpy.concatenate((BCFieldY, fldY), axis=1)

                    if BCFieldZ is None: BCFieldZ = fldZ
                    else: BCFieldZ = numpy.concatenate((BCFieldZ, fldZ), axis=1)

                elif foundVar and (len(foundVar)/nvars)%2 == 0: # 2D
                    # Config (XY)
                    if foundVar[0][-1] == 'X' and foundVar[1][-1] == 'Y':
                        for fgc in fldFaceTNC:
                            fgcX = numpy.asarray(fgc[1][0::2], order='K')
                            fgcY = numpy.asarray(fgc[1][1::2], order='K')

                            if fldX is None: fldX = fgcX
                            else: fldX = numpy.concatenate((fldX,fgcX), axis=1)

                            if fldY is None: fldY = fgcY
                            else: fldY = numpy.concatenate((fldY,fgcY), axis=1)

                        indp    = indFaceTNC.ravel(order='K')
                        fldX    = fldX.ravel(order='K').reshape(nvars,-1)
                        fldY    = fldY.ravel(order='K').reshape(nvars,-1)

                        if indices is None: indices = indp
                        else: indices = numpy.concatenate((indices, indp))

                        if BCFieldX is None: BCFieldX = fldX
                        else: BCFieldX = numpy.concatenate((BCFieldX, fldX), axis=1)

                        if BCFieldY is None: BCFieldY = fldY
                        else: BCFieldY = numpy.concatenate((BCFieldY, fldY), axis=1)

                    # Config (XZ)
                    if foundVar[0][-1] == 'X' and foundVar[1][-1] == 'Z':
                        for fgc in fldFaceTNC:
                            fgcX = numpy.asarray(fgc[1][0::2], order='K')
                            fgcZ = numpy.asarray(fgc[1][1::2], order='K')

                            if fldX is None: fldX = fgcX
                            else: fldX = numpy.concatenate((fldX,fgcX), axis=1)

                            if fldZ is None: fldZ = fgcZ
                            else: fldZ = numpy.concatenate((fldZ,fgcZ), axis=1)

                        indp    = indFaceTNC.ravel(order='K')
                        fldX    = fldX.ravel(order='K').reshape(nvars,-1)
                        fldZ    = fldZ.ravel(order='K').reshape(nvars,-1)

                        if indices is None: indices = indp
                        else: indices = numpy.concatenate((indices, indp))

                        if BCFieldX is None: BCFieldX = fldX
                        else: BCFieldX = numpy.concatenate((BCFieldX, fldX), axis=1)

                        if BCFieldZ is None: BCFieldZ = fldZ
                        else: BCFieldZ = numpy.concatenate((BCFieldZ, fldZ), axis=1)

                    # Config (YZ)
                    if foundVar[0][-1] == 'Y' and foundVar[1][-1] == 'Z':
                        for fgc in fldFaceTNC:
                            fgcY = numpy.asarray(fgc[1][0::2], order='K')
                            fgcZ = numpy.asarray(fgc[1][1::2], order='K')

                            if fldY is None: fldY = fgcY
                            else: fldY = numpy.concatenate((fldY,fgcY), axis=1)

                            if fldZ is None: fldZ = fgcZ
                            else: fldZ = numpy.concatenate((fldZ,fgcZ), axis=1)

                        indp    = indFaceTNC.ravel(order='K')
                        fldY    = fldY.ravel(order='K').reshape(nvars,-1)
                        fldZ    = fldZ.ravel(order='K').reshape(nvars,-1)

                        if indices is None: indices = indp
                        else: indices = numpy.concatenate((indices, indp))

                        if BCFieldY is None: BCFieldY = fldY
                        else: BCFieldY = numpy.concatenate((BCFieldY, fldY), axis=1)

                        if BCFieldZ is None: BCFieldZ = fldZ
                        else: BCFieldZ = numpy.concatenate((BCFieldZ, fldZ), axis=1)

        # =================================================================================

        if BCFieldX is None: BCFieldX = [None]*nvars
        if BCFieldY is None: BCFieldY = [None]*nvars
        if BCFieldZ is None: BCFieldZ = [None]*nvars

        if not isArrayEmpty(f):
            for i in range(nvars):
                if isArrayEmpty(f[i]): continue
                centers = Post.computeDiv2(x, f[i], vol, cellN, indices=indices,
                                           BCFieldX=BCFieldX[i], BCFieldY=BCFieldY[i],
                                           BCFieldZ=BCFieldZ[i])
                C.setFields([centers], z, 'centers')

    # Conditional clean up of partial derivatives
    if rmVar:
        for i in range(0, 3*nvars, 3):
            C._rmVars(t, ['{}:{}'.format('centers', v) for v in varList[i:i+ndims]])
            C._rmBCDataVars(t, varList[i:i+ndims])

    return None

def computeCurl(t, vector):
    """Compute the curl of a vector defined in array.
    Usage: computeCurl(t, vector) """
    tp = Internal.copyRef(t)
    curlVector = []
    n = len(vector)
    posv = [0]*n
    for i in range(n):
        var = vector[i]
        curlVar = var
        v = var.split(':')
        if len(v) > 1:
            if v[0] == 'centers':
                posv[i] = C.isNamePresent(tp, v[1])
                tp = C.center2Node(tp, var)
                curlVar = v[1]
        curlVector.append(curlVar)

    nodes = C.getAllFields(tp, 'nodes')
    for i in range(n):
        if posv[i] == -1: tp = C.rmVars(tp, curlVector[i])

    centers = Post.computeCurl(nodes, curlVector)
    C.setFields(centers, tp, 'centers')
    return tp

def computeNormCurl(t, vector):
    """Compute the norm of the curl of a vector defined in array.
    Usage: computeNormCurl(t, vector) """
    tp = Internal.copyRef(t)
    curlVector = []
    n = len(vector)
    posv = [0]*n
    for i in range(n):
        var = vector[i]
        curlVar = var
        v = var.split(':')
        if len(v) > 1:
            if v[0] == 'centers':
                posv[i] = C.isNamePresent(tp, v[1])
                tp = C.center2Node(tp, var)
                curlVar = v[1]
        curlVector.append(curlVar)

    nodes = C.getAllFields(tp, 'nodes')
    for i in range(n):
        if posv[i] == -1: tp = C.rmVars(tp, curlVector[i])
    centers = Post.computeNormCurl(nodes, curlVector)
    C.setFields(centers, tp, 'centers')
    return tp

def computeDiff(t, var):
    """Compute the difference of a variable defined in array.
    Usage: computeDiff(t, var) """
    v = var.split(':'); loc = 'nodes'
    if len(v) > 1:
        if v[0] == 'centers':
            loc = 'centers'; var = v[1]
    tp = Internal.copyRef(t)
    Internal._addGhostCells(tp, tp, 1, adaptBCs=0)
    nodes = C.getAllFields(tp, loc)
    res = Post.computeDiff(nodes, var)
    C.setFields(res, tp, loc)
    tp = Internal.rmGhostCells(tp, tp, 1, adaptBCs=0)
    return tp

def perlinNoise(t, alpha=2., beta=2., n=8):
    """Generate a perlin noise."""
    return C.TZGC1(t, 'nodes', True, Post.perlinNoise, alpha, beta, n)

def _perlinNoise(t, alpha=2., beta=2., n=8):
    return C._TZGC1(t, 'nodes', True, Post.perlinNoise, alpha, beta, n)

def streamLine(t, X0, vector, N=2000, dir=2):
    """Compute a streamline starting from (x0,y0,z0) given
    a list of arrays containing 'vector' information.
    Usage: streamLine(t, (x0,y0,z0), (vx,vy,vz), N, dir)"""
    t = C.center2Node(t, Internal.__FlowSolutionCenters__)
    arrays = C.getAllFields(t, 'nodes')
    a = Post.streamLine(arrays, X0, vector, N, dir)
    return C.convertArrays2ZoneNode('streamLine', [a])

# Retourne une liste de zones
def streamLine2(t, X0, vector, N=2000, dir=2, eps=1.e-2):
    """Compute a streamline starting from (x0,y0,z0) given
    a list of arrays containing 'vector' information.
    Usage: streamLine2(t, (x0,y0,z0), (vx,vy,vz), N, dir)"""
    t = C.center2Node(t, Internal.__FlowSolutionCenters__)
    arrays = C.getAllFields(t, 'nodes')
    for c, v in enumerate(vector): vector[c] = v.replace('centers:', '')
    a = Post.streamLine2(arrays, X0, vector, N, dir, eps)
    out = []
    for i in a:
        if i is not None and i[2] > 0: # if not empty
            out.append(C.convertArrays2ZoneNode('streamLine', [i]))
    return out

def streamSurf(t, b, vector, N=2000, dir=1):
    """Compute a streamsurf starting from BAR b given
    a list of arrays containing 'vector' information.
    Usage: streamSurf(t, b, (vx,vy,vz), N, dir)"""
    t = C.center2Node(t, Internal.__FlowSolutionCenters__)
    arrays = C.getAllFields(t, 'nodes')
    bar = C.getFields(Internal.__GridCoordinates__, b)[0]
    a = Post.streamSurf(arrays, bar, vector, N, dir)
    return C.convertArrays2ZoneNode('streamSurf', [a])

def streamRibbon(t, X0, N0, vector, N=2000, dir=2):
    """Compute a streamribbon starting from (x0,y0,z0) given
    a list of arrays containing 'vector' information. The width and orientation
    of the ribbon is given by N0 the normal vector to the streamline starting
    from X0.
    Usage: streamRibbon(arrays, (x0,y0,z0), (n0x,n0y,n0z), (vx,vy,vz), N, dir)"""
    t = C.center2Node(t, Internal.__FlowSolutionCenters__)
    arrays = C.getAllFields(t, 'nodes')
    a = Post.streamRibbon(arrays, X0, N0, vector, N, dir)
    return C.convertArrays2ZoneNode('streamRibbon', [a])

def isoLine(t, var, value):
    """Compute an isoline on surface zones.
    Usage: isoLine(t, var, value)"""
    t = C.center2Node(t, Internal.__FlowSolutionCenters__)
    v, loc = Internal.fixVarName(var)
    arrays = C.getAllFields(t, 'nodes')
    a = Post.isoLine(arrays, v, value)
    return C.convertArrays2ZoneNode('isoLine', [a])

def isoSurf(t, var, value, vars=None, split='simple'):
    """Compute an iso surface in volume zones using marching tetra.
    Usage: isoSurf(t, var, value, vars, split)"""
    if vars is None:
        t = C.center2Node(t, Internal.__FlowSolutionCenters__)
    else:
        for v in vars:
            vs = v.split(':')
            vc = []; vn = ['CoordinateX','CoordinateY','CoordinateY']
            if len(vs) == 2 and vs[0] == 'centers': vc.append(vs[1])
            elif len(vs) == 2 and vs[0] == 'nodes': vn.append(vs[1])
            else: vn.append(v)
        if vc != []: t = C.center2Node(t, vc)

    zones = Internal.getZones(t)
    var, loc = Internal.fixVarName(var)
    ret = []
    for z in zones:
        if vars is None: array = C.getAllFields(z, 'nodes')[0]
        else: array = C.getFields([Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__], z, vn+vc)[0]
        try:
            a = Post.isoSurf(array, var, value, split)
            if a != []:
                zp = C.convertArrays2ZoneNode(z[0], a)
                zp[0] = z[0] # pour identification
                ret.append(zp)
        except: raise
    return ret

def isoSurfMC(t, var, value, vars=None, split='simple'):
    """Compute an iso surface in volume zones using marching cubes.
    Usage: isoSurfMC(t, var, value, vars, split)"""
    rmvar = False
    if vars is  None:
        target_var = Internal.__FlowSolutionCenters__
    else:
        target_var=[]; vc = []; vn = ['CoordinateX','CoordinateY','CoordinateZ'];
        if var not in vars:
            vars.append(var)
            rmvar = True
        for v in vars:
            vs = v.split(':')
            if len(vs) == 2 and vs[0] == 'centers':
                vc.append(vs[1])
                target_var.append(v)
            elif len(vs) == 2 and vs[0] == 'nodes': vn.append(vs[1])
            else: vn.append(v)

    zones = Internal.getZones(t)
    var, loc = Internal.fixVarName(var)
    ret = []
    for z in zones:
        purge = False
        if Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__) is None: purge=True
        if vars is  None or vc != []:
            z = C.center2Node(z, target_var )
        if vars is None:
            array = C.getAllFields(z, 'nodes')[0]
        else:
            array = C.getFields([Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__], z, vars=vn+vc)[0]
        try:
            a = Post.isoSurfMC(array, var, value, split)
            if purge: Internal._rmNodesByName1(z, Internal.__FlowSolutionNodes__)
            if a != []:
                zp = C.convertArrays2ZoneNode(z[0], a)
                zp[0] = z[0] # pour identification

                if rmvar:
                    node = Internal.getNodeFromName1(zp, Internal.__FlowSolutionNodes__)
                    Internal._rmNodesByName(node, var)

                ret.append(zp)
        except: raise

    return ret

def computeIndicatorField(octreeHexa, varName, nbTargetPts=-1, bodies=[],
                          refineFinestLevel=1, coarsenCoarsestLevel=1,
                          isAMR=False,valMin=0,valMax=1,isOnlySmallest=False):
    """Compute the indicator -1, 0 or 1 for each element of the HEXA octree
    (or quadtree) with respect to the indicatorValue variable varName located
    at element centers.
    The bodies set the indicator to 0 in the vicinity of bodies.
    nbTargetPts controls the number of points after adaptation.
    If refineFinestLevel=1, the finest levels are refined.
    If coarsenCoarsestLevel=1, the coarsest levels are coarsened wherever
    possible.
    Return the indicator field.
    Usage: computeIndicatorField(octreeHexa, indicVal, nbTargetPts, bodies,
    refineFinestLevel, coarsenCoarsestLevel)"""
    epsInf, epsSup = _computeIndicatorField(octreeHexa, varName, nbTargetPts, bodies,
                                            refineFinestLevel, coarsenCoarsestLevel,
                                            isAMR=isAMR,valMin=valMin,valMax=valMax,isOnlySmallest=isOnlySmallest)
    return octreeHexa, epsInf, epsSup

def _computeIndicatorField(octreeHexa, varName, nbTargetPts=-1, bodies=[],
                           refineFinestLevel=1, coarsenCoarsestLevel=1,
                           isAMR=False,valMin=0,valMax=1,isOnlySmallest=False):
    """Compute the indicator -1, 0 or 1 for each element of the HEXA octree
    (or quadtree) with respect to the indicatorValue variable varName located
    at element centers.
    The bodies set the indicator to 0 in the vicinity of bodies.
    nbTargetPts controls the number of points after adaptation.
    If refineFinestLevel=1, the finest levels are refined.
    If coarsenCoarsestLevel=1, the coarsest levels are coarsened wherever
    possible.
    Return the indicator field.
    Usage: computeIndicatorField(octreeHexa, indicVal, nbTargetPts, bodies,
    refineFinestLevel, coarsenCoarsestLevel)"""
    zvars = varName.split(':')
    hexa = C.getFields(Internal.__GridCoordinates__, octreeHexa)[0]
    bodiesA = []
    if bodies != []:
        bodiesA = C.getFields(Internal.__GridCoordinates__, bodies)
    fields = C.getField(varName, octreeHexa)[0]
    if zvars[0] != 'centers': fields = Converter.node2Center(fields)

    if isAMR:
        epsInf=valMin
        epsSup=valMax
        indicator = Post.computeIndicatorField_AMR(
            hexa, fields, nbTargetPts, bodiesA,
            refineFinestLevel, coarsenCoarsestLevel,valMin=valMin,valMax=valMax,isOnlySmallest=isOnlySmallest)
    else:
        indicator, epsInf, epsSup = Post.computeIndicatorField(
            hexa, fields, nbTargetPts, bodiesA,
            refineFinestLevel, coarsenCoarsestLevel)
    C.setFields([indicator], octreeHexa, 'centers')
    return epsInf, epsSup

def computeIndicatorValue(octreeHexa, t, varName):
    """Computes the indicator value on the octree mesh based on the maximum
    value of the absolute value of the field of name varName in t.
    The projected field is stored at elements of the octree mesh."""
    octreeHexa2 = Internal.copyRef(octreeHexa)
    _computeIndicatorValue(octreeHexa2, t, varName)
    return octreeHexa2

def _computeIndicatorValue(octreeHexa, t, varName):
    """Computes the indicator value on the octree mesh based on the maximum
    value of the absolute value of the field of name varName in t.
    The projected field is stored at elements of the octree mesh."""
    vars0 = varName.split(':')
    if vars0[0] == 'centers': loc = 'centers'
    else: loc = 'nodes'
    hexa = C.getFields(Internal.__GridCoordinates__, octreeHexa)[0]
    coords = C.getFields(Internal.__GridCoordinates__, t)
    field = C.getField(varName, t)
    if loc == 'centers': field = Converter.center2Node(field)
    indicValue = Post.computeIndicatorValue(hexa, coords, field)
    if loc == 'centers': indicValue = Converter.node2Center(indicValue)
    C.setFields([indicValue], octreeHexa, loc)
    return None

def sharpEdges(a, alphaRef=30.):
    """Detect sharp edges between adjacent cells of a surface. Angle out of
    [180-alpharef,180+alpharef] are considered as sharp.
    Usage: sharpEdges(a, alphaRef)"""
    arrays = C.getAllFields(a, 'nodes')
    res = Post.sharpEdges(arrays, alphaRef)
    zones = []
    for r in res: zones.append(C.convertArrays2ZoneNode('edge', [r]))
    return zones

def silhouette(a, vector):
    """Detect shape of an unstructured surface."""
    arrays = C.getAllFields(a, 'nodes')
    res = Post.silhouette(arrays, vector)
    zones = []
    for r in res: zones.append(C.convertArrays2ZoneNode('silhouette', [r]))
    return zones

#==============================================================================
# modify variable names
#==============================================================================
def renameVars(a, varsPrev, varsNew):
    """Rename variables names in varsPrev with names defined by varsNew."""
    t = Internal.copyRef(a)
    _renameVars(t, varsPrev, varsNew)
    return t

def _renameVars(a, varsPrev, varsNew):
    """Rename variables names in varsPrev with names defined by varsNew."""
    if len(varsPrev) != len(varsNew):
        raise ValueError("renameVars: lists of variables must be of same size.")

    fnodes = Internal.getNodesFromName3(a, Internal.__FlowSolutionNodes__)
    fcenters = Internal.getNodesFromName3(a, Internal.__FlowSolutionCenters__)
    for nov in range(len(varsPrev)):
        splP = varsPrev[nov].split(':')
        splN = varsNew[nov].split(':')
        if fcenters != []:
            if len(splP) != 1 and splP[0] == 'centers' and len(splN) != 1 and splN[0]=='centers':
                for fc in fcenters:
                    for j in fc[2]:
                        if j[0]==splP[1]: j[0] = splN[1]
        if fnodes != []:
            if len(splP) == 1 and len(splN) == 1:
                if fnodes != []:
                    for fn in fnodes:
                        for j in fn[2]:
                            if j[0]==splP[0]: j[0] = splN[0]
            elif len(splP) != 1 and splP[0] == 'nodes' and len(splN) != 1 and splN[0]=='nodes':
                if fnodes != []:
                    for fn in fnodes:
                        for j in fn[2]:
                            if j[0]==splP[1]: j[0] = splN[1]

# Statistique de cellN=1 par rapport au nombre total de pt des zones
def checkOccupancyCellN(lowerLimit, t):
    """Check cellN occupancy"""
    import Converter.Mpi as Cmpi
    total_cells        = 0

    list_zones         =[]
    list_occupancy     =[]
    list_zones_below   =[]

    t=Internal.rmGhostCells(t,t,2,adaptBCs=1) #want "real" cells

    for z in Internal.getZones(t):
        cells_zone = 0
        sol        = Internal.getNodeFromName(z, 'FlowSolution#Centers')
        celN       = Internal.getNodeFromName(sol, 'cellN')[1]

        total_cell_zone= C.getNCells(z)
        cells_zone_0   = numpy.count_nonzero(celN==0)
        cells_zone_1   = numpy.count_nonzero(celN==1)
        cells_zone_2   = numpy.count_nonzero(celN==2)

        cells_zone     = cells_zone_1

        total_cells  += cells_zone
        occupancy_rate=float(cells_zone)/float(total_cell_zone)*100

        list_zones.append(z[0])
        list_occupancy.append(occupancy_rate)
        if occupancy_rate < lowerLimit and cells_zone>0: list_zones_below.append(z[0])

    with open('CellN_occupancy_MPIrank_'+str(Cmpi.rank)+'.txt', 'w') as f:
        for i in range(len(list_zones)):
            string2write=list_zones[i]          + '   ' + \
                str(list_occupancy[i])+ ' \n'
            f.write(string2write)

        string2write='..........Below Threshold.......... \n'
        f.write(string2write)
        string2write='Total numberof zones ::'+str(len(list_zones_below))+ ' \n'
        f.write(string2write)
        string2write='List of zones='+str(list_zones_below) + ' \n'
        f.write(string2write)
        string2write='Min Occupancy='+str(numpy.min(list_occupancy)) + ' \n'
        f.write(string2write)
        string2write='Max Occupancy='+str(numpy.max(list_occupancy)) + '  \n'
        f.write(string2write)
        string2write='Mean Occupancy='+str(numpy.mean(list_occupancy))+'  \n'
        f.write(string2write)

    return list_zones,list_occupancy,list_zones_below

def printMinMaxAndErrors(t, printErrors=False, listVars=[]):
    """Check min/max of given fields"""
    if printErrors:
        errors = Internal.checkPyTree(t)
        print("Printing Errors")
        print(errors)

    if listVars:
        for z in Internal.getZones(t):
            print("______________________")
            #check max & min values of each zone
            for i in listVars: print('Zone=',z[0],'|',i,':', C.getMinValue(z, i), C.getMaxValue(z, i))
        print("______________________")
        #check max & min values
        for i in listVars: print('Tree|',i,':', C.getMinValue(t, i), C.getMaxValue(t, i))
    else:
        print("Warning: printMinMaxAndErrors: List of variables to print min & max values is EMPTY.")
    return None

# Ecrit la taille des zones de t
def printSizeZones(t, withGhost=False, isNCells=False, printFinal=False):
    """Print zone stats - Ncells per zone & full Tree"""
    import Converter.Mpi as Cmpi
    list_save_zones  =[]
    list_save_ncells =[]

    if not withGhost: Internal._rmGhostCells(t, t, 2, adaptBCs=0)

    total_cells = 0
    dict_cell_zone={}
    for z in Internal.getZones(t):
        if isNCells:var_loc = C.getNCells(z)
        else: var_loc = C.getNPts(z)
        total_cells += var_loc
        dict_cell_zone[z[0]] = var_loc
        list_save_zones.append(z[0])
        list_save_ncells.append(var_loc/1.e06)

    with open('Ncells_MPIrank_'+str(Cmpi.rank)+'.txt', 'w') as f:
        for i in range(len(list_save_zones)):
            string2write=list_save_zones[i]       + '   ' + \
                str(list_save_ncells[i])+ ' \n'
            f.write(string2write)

        if printFinal:
            string2write='Ncells Total='+str(total_cells/1.e06)
            f.write(string2write)
    return None

# probeLocations [CB:RR] - redondant peut-etre avec probe
def probeLocations(tprobe, tcase):
    import Connector.PyTree as X
    import Converter.Mpi as Cmpi
    ##Create list of Probes for tprobe
    interDict = X.getIntersectingDomains(tprobe,tcase)

    list_save_zones=[]
    list_save_i    =[]
    list_save_j    =[]
    list_save_k    =[]
    ##In intersection zone get i,j,k for smallest distance
    for p in Internal.getZones(tprobe):
        xnode=Internal.getNodeFromName(p,'CoordinateX')
        ynode=Internal.getNodeFromName(p,'CoordinateY')
        znode=Internal.getNodeFromName(p,'CoordinateZ')

        x_loc=Internal.getValue(xnode)
        y_loc=Internal.getValue(ynode)
        z_loc=Internal.getValue(znode)

        isave = 0
        jsave = 0
        ksave = 0
        z=interDict[p[0]]

        if z:
            z2=Internal.getNodeFromName(tcase, z[0])
            C._initVars(z2,'dist=sqrt(({CoordinateX}-%g)**2+({CoordinateY}-%g)**2+({CoordinateZ}-%g)**2)'%(x_loc,y_loc,z_loc))
            dist_array = Internal.getNodeByName(z2,'dist')[1]

            dim   = Internal.getZoneDim(z2)
            ni    = dim[1]
            nj    = dim[2]
            nk    = dim[3]
            pnt   = numpy.where(dist_array == numpy.amin(dist_array))

            list_save_zones.append(z[0])
            list_save_i.append(pnt[0][0])
            list_save_j.append(pnt[1][0])
            list_save_k.append(pnt[2][0])

    with open('probes_locations_ijk_MPIrank_'+str(Cmpi.rank)+'.txt', 'w') as f:
        for i in range(len(list_save_zones)):
            string2write=list_save_zones[i]  + '  '+ \
                str(list_save_i[i])+ '  '+ \
                str(list_save_j[i])+ '  '+ \
                str(list_save_k[i])+ '\n'
            f.write(string2write)

    return None


def cgns2tecplot(t,isRmGhost=False,isSet2Zero=False):
    listNameRemove=['CoordinateX', 'CoordinateY', 'CoordinateZ', 'TurbulentDistance', 'cellN']
    t = C.node2Center(t)
    if isRmGhost: t = Internal.rmGhostCells(t, t, 2, adaptBCs=0)
    if isSet2Zero:
        for z in Internal.getZones(t):
            sol        = Internal.getNodeFromName(z,'FlowSolution')
            varNames   = C.getVarNames(z)[0]
            for namem in listNameRemove:
                varNames.remove(namem)

            dist2walls = Internal.getNodeFromName(sol,'TurbulentDistance')[1]
            sh         = numpy.shape(dist2walls)
            for namem in varNames:
                if 'P1' in namem or 'M1' in namem:continue
                var       = Internal.getNodeFromName(sol,namem)[1]

                if len(sh) ==2:
                    for j in range(sh[1]):
                        for i in range(sh[0]):
                            if dist2walls[i,j]<0:
                                var[i,j]=0
                else:
                    for k in range(sh[2]):
                        for j in range(sh[1]):
                            for i in range(sh[0]):
                                if dist2walls[i,j]<0:
                                    var[i,j,k]=0
    return t
