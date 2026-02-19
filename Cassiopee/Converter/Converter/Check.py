# - check -
"""Check pyTree integrity."""""

from . import Internal
from . import PyTree as C
import numpy

# last is 78
CGNSTypes = {
    'CGNSTree_t':72,
    'CGNSLibraryVersion_t':1,
    'CGNSBase_t':0,
    'BaseIterativeData_t':20,
    'Zone_t':48,
    'Elements_t':52,
    'GridCoordinates_t':54,
    'FlowSolution_t':53,
    'ZoneGridConnectivity_t':64,
    'GridConnectivityProperty_t':65,
    'GridConnectivityType_t':66,
    'GridConnectivity1to1_t':67,
    'GridConnectivity_t':75,
    'Periodic_t':76,
    'OversetHoles_t':68,
    'ZoneIterativeData_t':69,
    'ZoneSubRegion_t':70,
    'ZoneType_t':71,

    'DataArray_t':3,
    'DataConversion_t':4,
    'DataClass_t':5,
    'Descriptor_t':6,
    'DimensionalExponents_t':7,
    'AdditionalExponents_t':8,
    'DimensionalUnits_t':9,
    'AdditionalUnits_t':10,
    'UserDefinedData_t':13,
    'DiscreteData_t':51,
    'Rind_t':50,
    'GridLocation_t':16,
    'Ordinal_t':17,
    'IndexArray_t':18,
    'IndexRange_t':19,
    '"int[IndexDimension]"':73,
    '"int"':74,

    'Axisymmetry_t':2,
    'AxisymmetryAxisVector_t':11,
    'AxisymmetryReferencePoint_t':12,

    'FamilyName_t':15,
    'AdditionalFamilyName_t':14,
    'Family_t':21,
    'FamilyBC_t':22,
    'FamilyBCDataSet_t':23,

    'GeometryReference_t':27,
    'GeometryEntity_t':28,
    'GeometryFile_t':29,
    'GeometryFormat_t':30,

    'ReferenceState_t':25,
    'RotatingCoordinates_t':31,
    'FlowEquationSet_t':32,
    'ChemicalKineticsModel_t':33,
    'EMConductivityModel_t':34,
    'EMElectricFieldModel_t':35,
    'GasModel_t':36,
    'GoverningEquations_t':37,
    'ThermalConductivityModel_t':38,
    'ThermalRelaxationModel_t':39,
    'TurbulenceClosure_t':40,
    'TurbulenceModel_t':41,
    'ViscosityModel_t':42,
    'ConvergenceHistory_t':43,
    'Gravity_t':44,
    'IntegralData_t':45,
    'ReferenceState_t':46,
    'SimulationType_t':47,
    'EquationDimension_t':78, # not SIDS (tolerate)

    'ArbitraryGridMotion_t':49,
    'RigidGridMotion_t':55,
    'TimeRigidMotion_t':77, # non standard

    'ZoneBC_t':56,
    'BC_t':57,
    'BCData_t':24,
    'BCDataSet_t':58,
    'BCProperty_t':59,
    'Area_t':60,
    'AreaType_t':61,
    'WallFunction_t':62,
    'WallFunctionType_t':63
}

#==============================================================================
# IN: t: pyTree to be checked
# IN: level: check level 0 (version node), 1 (node conformity),
# 2 (unique base name), 3 (unique zone name), 4 (unique BC name),
# 5 (BC ranges), 6 (BCMatch/NearMatch), 7 (FamilyZone et FamilyBCs),
# 8 (invalid CGNS Types), 9 (invalid connectivity),
# 10 (invalid field names), 11 NAN in fields, 12 name length
# if level=-n, perform check from 0 to n
# OUT: error = [noeud posant probleme, message]
#==============================================================================
def checkPyTree(t, level=-20):
    """Check different conformity in tree."""
    errors = []
    if level <= 0 or level == 0:
        # check version node
        errors += checkVersionNode(t)
    if level <= -1 or level == 1:
        # check nodes conformity
        errors += checkNodes(t)
    if level <= -2 or level == 2:
        # check unique base names
        errors += checkUniqueNames(t, 'CGNSBase_t')
    if level <= -3 or level == 3:
        # check unique zone names
        errors += checkUniqueNames(t, 'Zone_t')
    if level <= -4 or level == 4:
        # check unique BC names
        errors += checkUniqueNames(t, 'BC_t')
        # check unique BCMatch names
        errors += checkUniqueNames(t, 'GridConnectivity1to1_t')
        # check unique BCNearMatch/BCOverlap names
        errors += checkUniqueNames(t, 'GridConnectivity_t')
    if level <= -5 or level == 5:
        # check BC range
        errors += checkBCRanges(t, 'BC_t')
        # check BCMatch range
        errors += checkBCRanges(t, 'GridConnectivity1to1_t')
        # check BCMatch opposite range
        errors += checkDonorRanges(t, 'GridConnectivity1to1_t')
        # check BCNearMatch/Overlap range
        errors += checkBCRanges(t, 'GridConnectivity_t')
        # check BCNearMatch opposite range
        errors += checkDonorRanges(t, 'GridConnectivity_t')
    if level <= -6 or level == 6:
        # check BCMatch opposite range
        errors += checkOppositRanges(t, 'GridConnectivity1to1_t')
        # check BCNearMatch/Overlap opposite range
        errors += checkOppositRanges(t, 'GridConnectivity_t')
    if level <= -7 or level == 7:
        # check zone family
        errors += checkZoneFamily(t)
        # check BC family
        errors += checkBCFamily(t)
    if level <= -8 or level == 8:
        # check CGNS type for each node
        errors += checkCGNSType(t)
    if level <= -9 or level == 9:
        # check element nodes (connectivity)
        errors += checkElementNodes(t)
    if level <= -10 or level == 10:
        # check valid CGNS var name
        errors += checkCGNSVarNames(t)
        # check Coordinates in fields
        #errors += checkCoordinatesInFields(t)
        # check field dimension
        #errors += checkFieldConformity(t)
    if level <= -11 or level == 11:
        # check NAN or infinite in fields
        errors += checkNAN(t)
    if level <= -12 or level == 12:
        # check if some names are longer than 32 chars
        errors += checkNameLength(t)
    if level <= -13 or level == 13:
        errors += checkBaseZonesDim(t)
    # Ne retourne que le noeud et le message dans les erreurs
    retErrors = []
    for i in range(len(errors)//3): retErrors += [errors[3*i], errors[3*i+2]]
    return retErrors

#==============================================================================
# Correct pyTree
#==============================================================================
def correctPyTree(t, level=-20):
    """Correct non conformities in tree."""
    tp = Internal.copyRef(t)
    _correctPyTree(tp, level)
    return tp

def _correctPyTree(t, level=-20):
    """Correct non conformities in tree."""
    # Corrige le noeud version
    if level <= 0 or level == 0:
        _correctVersionNode(t)
    # Supprime les noeuds non conformes
    if level <= -1 or level == 1:
        _correctNodes(t)
    # Renomme les bases
    if level <= -2 or level == 2:
        _correctNames(t, 'CGNSBase_t')
    # Renomme les zones
    if level <= -3 or level == 3:
        _correctNames(t, 'Zone_t')
    # Renomme les BCs
    if level <= -4 or level == 4:
        _correctNames(t, 'BC_t')
        _correctNames(t, 'GridConnectivity1to1_t')
        _correctNames(t, 'GridConnectivity_t')
    # Supprime les BCs avec des ranges invalides
    if level <= -5 or level == 5:
        _correctBCRanges(t, 'BC_t')
        _correctBCRanges(t, 'GridConnectivity1to1_t')
        _correctDonorRanges(t, 'GridConnectivity1to1_t')
        _correctBCRanges(t, 'GridConnectivity_t')
        _correctDonorRanges(t, 'GridConnectivity_t')
    # Supprime les BCs avec des fenetres opposees invalides
    if level <= -6 or level == 6:
        _correctOppositRanges(t, 'GridConnectivity1to1_t')
        _correctOppositRanges(t, 'GridConnectivity_t')
    # Corrige les family oubliees
    if level <= -7 or level == 7:
        _correctZoneFamily(t)
        _correctBCFamily(t)
    # Corrige les noeud pas de bon type CGNS
    if level <= -8 or level == 8:
        _correctCGNSType(t)
    # Corrige les noeuds elements (connectivity)
    if level <= -9 or level == 9:
        _correctElementNodes(t)
    # Corrige les noms de variables non CGNS
    if level <= -10 or level == 10:
        _correctCGNSVarNames(t)
        #_correctCoordinatesInFields(t)
        #_correctFieldConformity(t)
    # Supprime les NAN dans les champs
    if level <= -11 or level == 11:
        _correctNAN(t)
    # Tronque les noms > 32 chars
    if level <= -12 or level == 12:
        _correctNameLength(t)
    # Corrige la dim de la base
    if level <= -13 or level == 13:
        _correctBaseZonesDim(t, splitBases=True)
    C.registerAllNames(t)

    return None

#==============================================================================
# Check version node
# Doit etre en premier dans l'arbre et non duplique.
#==============================================================================
def checkVersionNode(t):
    """Check version node."""
    errors = []
    if len(t) == 4 and t[3] == 'CGNSTree_t':
        version = Internal.getNodesFromType1(t, 'CGNSLibraryVersion_t')
        if version == []: errors += [None, t, 'Missing CGNS version node.']
        elif t[2][0] != version[0]: errors += [version[0], t, 'CGNS version node must be first in tree sons.']
        if len(version) > 1:
            for v in version[1:]: errors += [v, t, 'Only one CGNS version node is allowed.']
    return errors

#==============================================================================
# Correct version node, en cree un ou supprime ceux en trop
#==============================================================================
def _correctVersionNode(t):
    """Correct version node."""
    errors = checkVersionNode(t)
    le = len(errors)//3
    added = 0
    for e in range(le):
        node = errors[3*e]
        if node is None: t[2].insert(0, Internal.createCGNSVersionNode())
        else:
            c = 0
            for n in t[2]:
                if n[3] == 'CGNSLibraryVersion_t' and id(n) == id(node):
                    del t[2][c]; break
                c += 1
            if added == 0: t[2].insert(0, node); added = 1
    return None

#==============================================================================
# Check nodes renvoie la liste des noeuds non conformes
# Check nodes modifie aussi les noeuds ayant comme valeur des chaines
# de caracteres avec des blancs au debut ou a la fin
#==============================================================================
def checkNodes(node):
    """Check basic node conformity."""
    errors = []
    isStd = Internal.isStdNode(node)
    if isStd >= 0:
        for c in node[isStd:]: checkNode__(c, node, errors)
    else: checkNode__(node, node, errors)
    return errors

#==============================================================================
def checkNode__(node, parent, errors):
    sons = []
    # node doit etre une liste
    if isinstance(node, list):
        # node doit avoir pour longueur 4 [nom, valeur, [fils], type]
        if len(node) == 4:
            # si node[1] est une string -> strip
            if isinstance(node[1], str): node[1] = node[1].strip()
            if isinstance(node[1], numpy.ndarray) and (node[1].dtype.kind == 'S' or node[1].dtype.kind == 'a'):
                val = Internal.getValue(node)
                if isinstance(val, str):
                    val = val.strip()
                    Internal.setValue(node, val)

            # node[0] (nom) est une string ou None
            if not isinstance(node[0], str) and node[0] is None:
                errors += [node, parent, "Node[0] of node %s must be a string designing node name."%node[0]]

            # node[2] (liste des fils) doit etre une liste
            if not isinstance(node[2], list):
                errors += [node, parent, "Node[2] of node %s must be a list of sons."%node[0]]
            else: sons = node[2]

            # node[3] doit etre une string se terminant par _t ...
            if not isinstance(node[3], str):
                errors += [node, parent, "Node[3] of node %s must be a string designing the node type."%node[0]]
            #if node[3][-2:] != '_t' and node[3] != "int[IndexDimension]":
            #    errors += [node, parent, "Node[3] of node %s must be a string designing the node type."%node[0]]

        else: errors += [node, parent, "Node %s has a length != 4."%node[0]]
    else: errors += [node, parent, "Node is not a list."]

    for n in sons: checkNode__(n, node, errors)

#==============================================================================
# Delete les noeuds non conformes
#==============================================================================
def _correctNodes(t):
    """Delete non conform nodes."""
    errors = checkNodes(t)
    le = len(errors)//3
    for e in range(le):
        node = errors[3*e]
        parent = errors[3*e+1]
        c = Internal.getNodePosition(node, parent)
        del parent[2][c]
    return None

#==============================================================================
# Check name length
# Check if node[0] has less than 32 chars (legacy constraint for cgnsview)
#==============================================================================
def checkNameLength(node):
    """Check node[0] length. Must be < 32 chars."""
    errors = []
    isStd = Internal.isStdNode(node)
    if isStd >= 0:
        for c in node[isStd:]: checkNameLength__(c, node, errors)
    else: checkNameLength__(node, node, errors)
    return errors

#==============================================================================
def checkNameLength__(node, parent, errors):
    sons = []
    if len(node[0]) > 32:
        errors += [node, parent, "Node name %s has a length > 32."%node[0]]
    for n in sons: checkNameLength__(n, node, errors)

#==============================================================================
# Truncate > 32 node names
#==============================================================================
def _correctNameLength(t):
    """Truncate node names if necessary."""
    errors = checkNameLength(t)
    le = len(errors)//3
    for e in range(le):
        node = errors[3*e]
        parent = errors[3*e+1]
        name = node[0]
        if node[3] == 'Zone_t':
            newName = name[0:31]
            newName = C.getZoneName(newName)
            if len(newName) > 32: newName = newName[0:32]
            node[0] = newName
            Internal._renameNode(t, name, newName)
        elif node[3] == 'BC_t':
            newName = name[0:31]
            newName = C.getBCName(newName)
            if len(newName) > 32: newName = newName[0:32]
            node[0] = newName
            Internal._renameNode(t, name, newName)
        elif node[3] == 'CGNSBase_t':
            newName = name[0:31]
            newName = C.getBaseName(newName)
            if len(newName) > 32: newName = newName[0:32]
            node[0] = newName
            Internal._renameNode(t, name, newName)
        else:
            node[0] = name[0:32]

        print("INFO: replacing %s with %s."%(name, node[0]))
    return None

#==============================================================================
# Verifie que les noms des noeuds de type donne sont bien uniques
# Optimise pour Base, Zone, BC et GC
#==============================================================================
def checkUniqueNames(t, ntype):
    nameServer = {}
    errors = []

    if ntype == 'CGNSBase_t':
        nodes = Internal.getBases(t)
        for n in nodes:
            name = n[0]
            if name not in nameServer: nameServer[name] = 0
            else: errors += [n, t, "Base name %s is already used."%n[0]]

    elif ntype == 'Zone_t':
        bases = Internal.getBases(t)
        for b in bases:
            nodes = Internal.getNodesFromType1(b, 'Zone_t')
            for n in nodes:
                name = n[0]
                if name not in nameServer: nameServer[name] = 0
                else: errors += [n, b, "Zone name %s is already used."%n[0]]

    elif ntype == 'BC_t':
        zones = Internal.getZones(t)
        for z in zones:
            zbcs = Internal.getNodesFromType1(z, 'ZoneBC_t')
            for zbc in zbcs:
                nodes = Internal.getNodesFromType1(zbc, 'BC_t')
                for n in nodes:
                    name = n[0]
                    if name not in nameServer: nameServer[name] = 0
                    else: errors += [n, zbc, "BC name %s (zone %s) is already used."%(n[0],z[0])]

    elif ntype == 'GridConnectivity_1to1_t':
        zones = Internal.getZones(t)
        for z in zones:
            gcs = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
            for gc in gcs:
                nodes = Internal.getNodesFromType1(gc, 'GridConnectivity_1to1_t')
                for n in nodes:
                    name = n[0]
                    if name not in nameServer: nameServer[name] = 0
                    else: errors += [n, gc, "BCMatch name %s (zone %s) is already used."%(n[0],z[0])]

    elif ntype == 'GridConnectivity_t':
        zones = Internal.getZones(t)
        for z in zones:
            gcs = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
            for gc in gcs:
                nodes = Internal.getNodesFromType1(gc, 'GridConnectivity_t')
                for n in nodes:
                    name = n[0]
                    if name not in nameServer: nameServer[name] = 0
                    else: errors += [n, gc, "GridConnectivity name %s (zone %s) is already used."%(n[0],z[0])]
    else:
        nodes = Internal.getNodesFromType(t, ntype)
        for n in nodes:
            name = n[0]
            if name not in nameServer: nameServer[name] = 0
            else:
                (p,c) = Internal.getParentOfNode(t, n)
                errors += [n, p, "Node name %s is already used."%n[0]]
    return errors

#==============================================================================
# Change node name if already defined
#==============================================================================
def _correctNames(t, ntype):
    nameServer = {}
    if ntype == 'CGNSBase_t': nodes = Internal.getBases(t)
    elif ntype == 'Zone_t': nodes = Internal.getZones(t)
    else: nodes = Internal.getNodesFromType(t, ntype)
    zoneDonors = []
    for n in nodes:
        name = n[0]
        if name not in nameServer: nameServer[name] = 0
        else: # deja existant
            c = nameServer[name]; ret = 1
            while ret == 1:
                name2 = '%s.%d'%(name,c)
                if name2 not in nameServer: ret = 0
                else: ret = 1
                c += 1
            nameServer[name2] = 0
            nameServer[name] = c
            if n[3] == 'Zone_t': zoneDonors.append((n[0], name2))
            n[0] = name2

    # Modifie les zoneDonors
    _correctDonors(t, 'GridConnectivity1to1_t', zoneDonors)
    _correctDonors(t, 'GridConnectivity_t', zoneDonors)

    # Modifie les attachs
    return None

#==============================================================================
# Corrige les occurence d'un nom de zone dans les BCs
def _correctDonors(t, ntype, zoneDonors):
    zones = Internal.getZones(t)
    for z in zones:
        nodes = Internal.getNodesFromType2(z, ntype)
        for n in nodes:
            zdonorname = Internal.getValue(n)
            for zd in zoneDonors:
                if zd[0] == zdonorname: Internal._setValue(n, zd[1])
    return None

#==============================================================================
# Check BC ranges
# IN: t: arbre a verifier
# IN: ntype: type du noeud de BC a verifier (BC_t,...)
# Verifie que le range de la BC est contenu dans la grille
# Verifie que les faces de la BC sont contenues dans la grille
# Verifie que la fenetre n'est pas volumique
# Corrige la shape des ranges si celle-ci est en C
#==============================================================================
def checkBCRanges(t, ntype):
    errors = []
    if ntype == 'BC_t': ctype = 'ZoneBC_t'
    else: ctype = 'ZoneGridConnectivity_t'
    zones = Internal.getZones(t)
    for z in zones:
        dim = Internal.getZoneDim(z)
        if dim[0] != 'Structured': continue  # pas de check en non structure
        bcs = Internal.getNodesFromType1(z, ctype)
        for bc in bcs:
            nodes = Internal.getNodesFromType1(bc, ntype)
            for n in nodes:
                prange = Internal.getNodesFromName1(n, 'PointRange')
                for r in prange:
                    if r[1].shape == (2,3):
                        r[1] = numpy.reshape(r[1], (3,2), order='F')
                    if r[1].shape == (3,2):
                        win = Internal.range2Window(r[1])
                    else: win = None # pas de check en non structure
                    # Check structure uniquement pour l'instant
                    error = 0
                    if win is not None:
                        # Volumic window
                        if win[0] != win[1] and win[2] != win[3] and win[4] != win[5]: error = 1
                        # imin out of bounds
                        if win[0] < 0 or win[0] > dim[1]: error = 1
                        # imax out of bounds
                        if win[1] < 0 or win[1] > dim[1]: error = 1
                        # jmin out of bounds
                        if win[2] < 0 or win[2] > dim[2]: error = 1
                        # jmax out of bounds
                        if win[3] < 0 or win[3] > dim[2]: error = 1
                        # kmin out of bounds
                        if win[4] < 0 or win[4] > dim[3]: error = 1
                        # kmax out of bounds
                        if win[5] < 0 or win[5] > dim[3]: error = 1

                    if error == 1:
                        errors += [n, bc, "Range of BC %s is invalid for zone %s."%(n[0],z[0])]
    return errors

def checkBCFaces(t, ntype):
    errors = []
    if ntype == 'BC_t': ctype = 'ZoneBC_t'
    else: ctype = 'ZoneGridConnectivity_t'
    zones = Internal.getZones(t)
    for z in zones:
        bcs = Internal.getNodesFromType1(z, ctype)
        for bc in bcs:
            r = Internal.getElementRange(bc, type="NGON")
            if r is not None: nfaces = r[1]-r[0]+1
            else: nfaces = 0
            nodes = Internal.getNodesFromType2(z, ntype)
            for n in nodes:
                plist = Internal.getNodesFromName1(n, 'PointList')
                for r in plist:
                    faces = r[1]
                    faces1 = faces[faces > nfaces]
                    faces2 = faces[faces < 1]
                    if faces1.size > 0 or faces2.size > 0:
                        errors += [n, bc, "Faces of BC %s is invalid for zone %s."%(n[0],z[0])]
    return errors

#==============================================================================
# Check donor BC ranges
# IN: t: arbre a verifier
# IN: ntype: ntype du noeud de BC a verifier (BC_t,...)
# On verifie que le range du donneur est contenu dans la grille donneur
#==============================================================================
def checkDonorRanges(t, ntype):
    errors = []
    if ntype == 'BC_t': ctype = 'ZoneBC_t'
    else: ctype = 'ZoneGridConnectivity_t'
    zones = Internal.getZones(t)
    for z in zones:
        bcs = Internal.getNodesFromType1(z, ctype)
        for bc in bcs:
            nodes = Internal.getNodesFromType1(bc, ntype)
            for n in nodes:
                donorName = Internal.getValue(n)
                donors = Internal.getNodesFromName2(t, donorName)
                if donors != []:
                    if all([Internal.getType(d) == 'Zone_t' for d in donors]):
                        dim = Internal.getZoneDim(donors[0])
                        if dim[0] != 'Structured': continue  # pas de check en non structure
                        r = Internal.getElementRange(donors[0], type="NGON")
                        if r is not None: nfaces = r[1]-r[0]+1
                        else: nfaces = 0
                        prange = Internal.getNodesFromName1(n, 'PointRangeDonor')
                        for r in prange:
                            if r[1].shape == (2,3):
                                r[1] = numpy.reshape(r[1], (3,2), order='F')
                            if r[1].shape == (3,2):
                                win = Internal.range2Window(r[1])
                            else: win = [0,0,0,0,0,0] # pas de check en NS
                            error = 0
                            if win[0] != win[1] and win[2] != win[3] and win[4] != win[5]: error = 1
                            if win[0] < 0 or win[0] > dim[1]: error = 1
                            if win[1] < 0 or win[1] > dim[1]: error = 1
                            if win[2] < 0 or win[2] > dim[2]: error = 1
                            if win[3] < 0 or win[3] > dim[2]: error = 1
                            if win[4] < 0 or win[4] > dim[3]: error = 1
                            if win[5] < 0 or win[5] > dim[3]: error = 1

                            if error == 1:
                                errors += [n, bc, "Range of donor BC %s is invalid for zone %s."%(n[0],z[0])]

                    elif all([Internal.getType(d) == 'Family_t' for d in donors]):
                        if not all([Internal.getName(d) == donorName for d in donors]):
                            errors += [n, bc, "Type donor %s %s is invalid for zone %s."%(ntype,n[0],z[0])]
                    else:
                        errors += [n, bc, "Type donor %s %s is invalid for zone %s."%(ntype,n[0],z[0])]
    return errors

def checkDonorFaces(t, ntype):
    errors = []
    if ntype == 'BC_t': ctype = 'ZoneBC_t'
    else: ctype = 'ZoneGridConnectivity_t'
    zones = Internal.getZones(t)
    for z in zones:
        bcs = Internal.getNodesFromType1(z, ctype)
        for bc in bcs:
            r = Internal.getElementRange(bc, type="NGON")
            if r is not None: nfaces = r[1]-r[0]+1
            else: nfaces = 0
            nodes = Internal.getNodesFromType1(bc, ntype)
            for n in nodes:
                donorName = Internal.getValue(n)
                donors = Internal.getNodesFromName2(t, donorName)
                if donors != []:
                    plist = Internal.getNodesFromName1(n, 'PointListDonor')
                    for r in plist:
                        faces = r[1]
                        faces1 = faces[faces > nfaces]
                        faces2 = faces[faces < 1]
                        if faces1.size > 0 or faces2.size > 0:
                            errors += [n, bc, "Faces of donor BC %s is invalid for zone %s."%(n[0],z[0])]
    return errors

#==============================================================================
# Supprime les BCs de type donne avec des ranges invalides
#==============================================================================
def _correctBCRanges(t, ntype):
    errors = checkBCRanges(t, ntype)
    le = len(errors)//3
    for e in range(le):
        node = errors[3*e]
        parent = errors[3*e+1]
        c = Internal.getNodePosition(node, parent)
        del parent[2][c]
    return None

#==============================================================================
# Supprime les BCs avec des donor ranges invalides
#==============================================================================
def _correctDonorRanges(t, ntype):
    errors = checkDonorRanges(t, ntype)
    le = len(errors)//3
    for e in range(le):
        node = errors[3*e]
        parent = errors[3*e+1]
        c = Internal.getNodePosition(node, parent)
        del parent[2][c]
    return None

#==============================================================================
# modify PointRange min/max into max/min for 1to1 GC if Transform index
# is negative - to be compliant by the standard (connectMatch always do min/max)
#==============================================================================
def _reorderBCMatchPointRange(t):
    for z in Internal.getZones(t):
        for gc in Internal.getNodesFromType2(z, 'GridConnectivity1to1_t'):
            TR = Internal.getNodeFromName1(gc, 'Transform')
            TR = Internal.getValue(TR)
            trirac1 = TR[0]; trirac2 = TR[1]; trirac3 = TR[2]
            PRN = Internal.getNodeFromName1(gc, 'PointRangeDonor')
            win = Internal.range2Window(Internal.getValue(PRN))
            [imin, imax, jmin, jmax, kmin, kmax] = win
            if imin != imax and trirac1<0:
                win[0] = imax; win[1] = imin
            if jmin != jmax and trirac2<0:
                win[2] = jmax; win[3] = jmin
            if kmin != kmax and trirac3<0:
                win[4] = kmax; win[5] = kmin
            PR = Internal.window2Range(win)
            PRN[1] = PR
    return None

#==============================================================================
# Verifie les ranges des fenetres opposees
# Le donneur doit exister et les ranges etre coherents
# IN: ntype: GridConnectivity1to1_t ou GridConnectivity_t
#==============================================================================
def checkOppositRanges(t, ntype):
    errors = []
    delta = numpy.empty(3, Internal.E_NpyInt); deltaopp = numpy.empty(3, Internal.E_NpyInt)
    zones = Internal.getZones(t)
    for z in zones:
        zname = z[0]
        bcs = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
        for bc in bcs:
            nodes = Internal.getNodesFromType1(bc, ntype)
            dimZone = Internal.getZoneDim(z)[4]
            for n in nodes:
                prange = Internal.getNodesFromName1(n, 'PointRange')
                prangedonor = Internal.getNodesFromName2(n, 'PointRangeDonor')#NearMatch : necessaire d aller au niveau 2
                mtype = Internal.getNodeFromName1(n, 'GridConnectivityType')
                if mtype is not None: mtype = Internal.getValue(mtype)
                else: mtype = 'Match'
                zdonorname = Internal.getValue(n)
                zdonor = Internal.getNodesFromName2(t, zdonorname)
                if zdonor == []:
                    errors += [n, bc, "Donor zone %s of BC %s (zone %s) does not exist."%(zdonorname,n[0],z[0])]
                else:
                    if mtype != 'Overset' and prange != [] and prangedonor != []:
                        for z in zdonor:
                            if z[3] == 'Zone_t': zdonor = z; break
                        # Verifie que le donneur est dans la meme base (trop cher)
                        #(b1, c1) = Internal.getParentOfNode(t, z)
                        #(b2, c2) = Internal.getParentOfNode(t, zdonor)
                        #if (b1[0] != b2[0] and mtype == 'Match'):
                        #    errors += [n, "Donor zone %s of BC %s (zone %s) is not in the same base as zone %s."%(zdonorname,n[0],z[0],z[0])]
                        #else:
                        # Verifie que la paire (prange,prangedonor) existe bien ds la zone zdonor
                        nodesopp = Internal.getNodesFromType2(zdonor, ntype)
                        dim = Internal.getZoneDim(zdonor)
                        error = 1
                        for nopp in nodesopp:
                            if n is not nopp:
                                prangeopp = Internal.getNodesFromName1(nopp, 'PointRange')
                                prangedonoropp = Internal.getNodesFromName2(nopp, 'PointRangeDonor') # NearMatch: necessaire d'aller au niveau 2
                                mtypeopp = Internal.getNodeFromName1(nopp, 'GridConnectivityType')
                                if mtypeopp is not None: mtypeopp = Internal.getValue(mtypeopp)
                                else: mtypeopp = 'Match'
                                zoppdonorname = Internal.getValue(nopp)
                                if zoppdonorname == zname and mtype == mtypeopp:
                                    # current zone
                                    rangez = Internal.range2Window(prange[0][1])
                                    rangezd = Internal.range2Window(prangedonor[0][1])
                                    # donor zone
                                    rangezopp = Internal.range2Window(prangeopp[0][1])
                                    rangezoppd = Internal.range2Window(prangedonoropp[0][1])
                                    if rangez == rangezoppd and rangezd == rangezopp: error = 0
                        if error == 1:
                            errors += [n, bc, "Opposite window from zone %s of BC %s (zone %s) does not exist."%(zdonorname,n[0],zname)]
                        # Check des ranges
                        for ropp in prangedonor:
                            if ropp[1].shape == (2,3):
                                ropp[1] = numpy.reshape(ropp[1], (3,2), order='F')
                            if ropp[1].shape == (3,2):
                                winopp = Internal.range2Window(ropp[1])
                            else: winopp = [-1,0,0,0,0,0]
                            error = 0
                            if winopp[0] != winopp[1] and winopp[2] != winopp[3] and winopp[4] != winopp[5]: error = 1
                            if winopp[0] < 0 or winopp[0] > dim[1]: error = 1
                            if winopp[1] < 0 or winopp[1] > dim[1]: error = 1
                            if winopp[2] < 0 or winopp[2] > dim[2]: error = 1
                            if winopp[3] < 0 or winopp[3] > dim[2]: error = 1
                            if winopp[4] < 0 or winopp[4] > dim[3]: error = 1
                            if winopp[5] < 0 or winopp[5] > dim[3]: error = 1
                            if error == 1:
                                errors += [n, bc, "Donor range of BC %s is invalid for zone %s."%(n[0],z[0])]
                            else:
                                if ntype == 'GridConnectivity1to1_t': # BCMatch only
                                    # check consistency of current and donor windows in each direction, taking into account Transform
                                    transform = Internal.getNodesFromName1(n, 'Transform')
                                    for r in prange:
                                        if r[1].shape == (2,3):
                                            r[1] = numpy.reshape(r[1], (3,2), order='F')
                                        if r[1].shape == (3,2):
                                            win = Internal.range2Window(r[1])
                                        else: win = [0,0,0,0,0,0]
                                        if transform != []: # not mandatory, [+1,+2,+3] by default.
                                            transform = transform[0][1]
                                        else:
                                            transform = numpy.empty(3, Internal.E_NpyInt)
                                            transform[0] = 1; transform[1] = 2; transform[2] = 3
                                        delta[0] = abs(win[1] - win[0]) # delta i for win
                                        if dimZone > 1: delta[1] = abs(win[3] - win[2]) # delta j for win
                                        if dimZone == 3: delta[2] = abs(win[5] - win[4]) # delta k for win
                                        deltaopp[0] = abs(winopp[1] - winopp[0]) # delta i for winopp
                                        if dimZone > 1: deltaopp[1] = abs(winopp[3] - winopp[2]) # delta j for winopp
                                        if dimZone == 3: deltaopp[2] = abs(winopp[5] - winopp[4]) # delta k for winopp
                                        if dimZone == 3:
                                            if ((delta[0] != deltaopp[abs(transform[0])-1]) or
                                                (delta[1] != deltaopp[abs(transform[1])-1]) or
                                                    (delta[2] != deltaopp[abs(transform[2])-1])):
                                                errors += [n, bc, "window of BC %s for zone %s does not match with its opposite window."%(n[0],z[0])]
                                        elif dimZone == 2:
                                            if ((delta[0] != deltaopp[abs(transform[0])-1]) or
                                                    (delta[1] != deltaopp[abs(transform[1])-1])):
                                                errors += [n, bc, "window of BC %s for zone %s does not match with its opposite window."%(n[0],z[0])]

                                        elif dimZone == 1:
                                            if delta[0] != deltaopp[abs(transform[0])-1]:
                                                errors += [n, bc, "window of BC %s for zone %s does not match with its opposite window."%(n[0],z[0])]
    return errors

#==============================================================================
# Efface les GCs qui n'ont pas de donneur existant
#                qui ont des ranges non coherents
#                dont le donneur n'a pas le noeud reciproque
#==============================================================================
def _correctOppositRanges(t, ntype):
    errors = checkOppositRanges(t, ntype)
    le = len(errors)//3
    for e in range(le):
        node = errors[3*e]
        parent = errors[3*e+1]
        c = Internal.getNodePosition(node, parent)
        del parent[2][c]
    return None

#==============================================================================
# Si une zone est taggee avec une famille, la famille doit etre declaree
# dans sa base
#==============================================================================
def checkZoneFamily(t):
    errors = []
    bases = Internal.getBases(t)
    for b in bases:
        zones = Internal.getZones(b)
        for z in zones:
            fs = Internal.getNodesFromType1(z, 'FamilyName_t')
            fs += Internal.getNodesFromType1(z, 'AdditionalFamilyName_t')
            for f in fs: # zone taggee
                name = Internal.getValue(f)
                # Check for family
                ref = Internal.getNodeFromName1(b, name)
                if ref is None:
                    errors += [b, t, 'FamilyZone %s (referenced by zone %s) is not defined in base.'%(name,z[0])]
                elif ref[3] != 'Family_t':
                    errors += [b, t, 'FamilyZone %s (referenced by zone %s) is not defined in base.'%(name,z[0])]
    return errors

#==============================================================================
# Si une famille n'est pas declaree, on l'ajoute dans la base
#==============================================================================
def _correctZoneFamily(t):
    errors = checkZoneFamily(t)
    le = len(errors)//3
    for e in range(le):
        node = errors[3*e]
        name = errors[3*e+2]
        name = name.split(' '); name = name[1]
        C._addFamily2Base(node, name)
    return None

#==============================================================================
# Si une zone utilise une familyBC, la famille doit etre declaree
# dans sa base
#==============================================================================
def checkBCFamily(t):
    errors = []
    bases = Internal.getBases(t)
    for base in bases:
        zones = Internal.getZones(base)
        for z in zones:
            BCs = Internal.getNodesFromType2(z, 'BC_t')
            for b in BCs:
                f = Internal.getNodeFromType1(b, 'FamilyName_t')
                if f is not None: # zone avec BC family
                    name = Internal.getValue(f)
                    # Check for family
                    refs = Internal.getNodesFromName1(base, name)
                    for ref in refs:
                        if ref[3] != 'Family_t': refs.remove(ref)
                    if refs == []:
                        errors += [base, t, 'FamilyBC %s (referenced by zone %s) is not defined in base.'%(name,z[0])]
                    elif refs[0][3] != 'Family_t':
                        errors += [base, t, 'FamilyBC %s (referenced by zone %s) is not defined in base.'%(name,z[0])]
    return errors

#==============================================================================
# Si une famille BC n'est pas declaree, on l'ajoute dans la base
#==============================================================================
def _correctBCFamily(t):
    errors = checkBCFamily(t)
    le = len(errors)//3
    for e in range(le):
        node = errors[3*e]
        name = errors[3*e+2]
        name = name.split(' '); name = name[1]
        C._addFamily2Base(node, name, bndType='UserDefined')
    return None

#==============================================================================
# Verifie que dans une base, toutes les zones ont le meme cellDim que la base
#==============================================================================
def checkBaseZonesDim(t):
    errors = []
    bases = Internal.getBases(t)
    for b in bases:
        dimBase = b[1][0]
        zones = Internal.getNodesFromType1(b, 'Zone_t')
        for z in zones:
            dim = Internal.getZoneDim(z)
            if dim[4] != dimBase: errors += [b, b, "Zone %s cellDim (%d) is inconsistent with base %s cellDim (%d)."%(z[0],dim[4],b[0],dimBase)]
    return errors

#==============================================================================
# Update all bases with the max dim found in their zones
# if splitBases: enforce bases to have homogenous cellDims by splitting base
#==============================================================================
def _correctBaseZonesDim(t, splitBases=False):
    bases = Internal.getBases(t)
    for b in bases:
        zones = Internal.getZones(b)
        z1 = []; z2 = []; z3 = [] # zones de dim 1,2,3
        listOfAddBases = []
        listOfRmBases = []
        for z in zones:
            try:
                dim = Internal.getZoneDim(z)
                if dim[4] <= 1: z1.append(z)
                elif dim[4] == 2: z2.append(z)
                elif dim[4] == 3: z3.append(z)
            except: pass
        lz1 = len(z1); lz2 = len(z2); lz3 = len(z3)
        lzmax = max(lz1, lz2, lz3)
        # put max dim in base node
        if lz3 > 0: b[1][0] = 3
        elif lz2 > 0: b[1][0] = 2
        elif lz1 > 0: b[1][0] = 1
        else: pass

        if lz1+lz2+lz3 > lzmax and splitBases:
            listOfRmBases.append(b[0])
            if lz1 != 0:
                listOfAddBases.append([1, b[0]+'.1', z1])
            if lz2 != 0:
                listOfAddBases.append([2, b[0]+'.2', z2])
            if lz3 != 0:
                listOfAddBases.append([3, b[0]+'.3', z3])

    if splitBases:
        for baseName in listOfRmBases:
            Internal._rmNodeByPath(t, baseName)

        for b in listOfAddBases:
            cellDim, baseName, zones = b
            base = Internal.newCGNSBase(baseName, cellDim=cellDim, physDim=3, parent=t)
            base[2] += Internal.getZones(zones)

    return None

#===============================================================================
# check if the PointRanges for a zone z (ZoneBC ou ZoneGridConnectivity) are
# compatible with multigrid
#===============================================================================
def checkMGForBCRanges(z, ntype, multigrid, sizemin):
    puiss = 2**(multigrid)
    errors = []
    nodes = Internal.getNodesFromType2(z, ntype)
    for n in nodes:
        PRS = Internal.getNodesFromName1(n, 'PointRange')
        for PR in PRS:
            [imin,imax,jmin,jmax,kmin,kmax] = Internal.range2Window(PR[1])
            if imin != imax:
                di = (imax-imin)
                res = (di+1)/puiss
                if di%puiss != 0 : errors+=[n,"BC %s of zone %s is not multigrid of level %d in direction i."%(n[0],z[0],multigrid)]
                elif res < sizemin:
                    errors+=[n,"PointRange for BC %s of zone %s: not enough points on coarse grid in direction i."%(n[0],z[0])]
            if jmin != jmax:
                dj = (jmax-jmin)
                res = (dj+1)/puiss
                if dj%puiss != 0: errors+=[n,"BC %s of zone %s is not multigrid of level %d in direction j."%(n[0],z[0],multigrid)]
                elif res < sizemin:
                    errors+=[n,"PointRange for BC %s of zone %s: not enough points on coarse grid in direction j."%(n[0],z[0])]
            if kmin != kmax:
                dk = (kmax-kmin)
                res = (dk+1)/puiss
                if dk%puiss != 0 : errors+=[n,"BC %s of zone %s is not multigrid of level %d in direction k."%(n[0],z[0],multigrid)]
                elif res < sizemin:
                    errors+=[n,"PointRange for BC %s of zone %s: not enough points on coarse grid in direction k."%(n[0],z[0])]
    return errors

#===============================================================================
# check if the PointRangeDonor for a zone z (ZoneBC ou ZoneGridConnectivity)
# are compatible with multigrid
#===============================================================================
def checkMGForDonorBCRanges(z, ntype, multigrid, sizemin):
    puiss = 2**(multigrid)
    errors = []
    nodes = Internal.getNodesFromType2(z, ntype)
    for n in nodes:
        PRS = Internal.getNodesFromName1(n, 'PointRangeDonor')
        for PR in PRS:
            [imin,imax,jmin,jmax,kmin,kmax] = Internal.range2Window(PR[1])
            if imin != imax:
                di = (imax-imin)
                res = (di+1)/puiss
                if di%puiss != 0:
                    errors+=[n,"PointRangeDonor for GC %s of zone %s is not multigrid of level %d in direction i."%(n[0],z[0],multigrid)]
                elif res < sizemin:
                    errors+=[n,"PointRangeDonor for GC %s of zone %s: not enough points on coarse grid in direction i."%(n[0],z[0])]
            if jmin != jmax:
                dj = (jmax-jmin)
                res = (dj+1)/puiss
                if dj%puiss != 0:
                    errors+=[n,"PointRangeDonor for GC %s of zone %s is not multigrid of level %d in direction j."%(n[0],z[0],multigrid)]
                elif res < sizemin:
                    errors+=[n,"PointRangeDonor for GC %s of zone %s: not enough points on coarse grid in direction j."%(n[0],z[0])]
            if kmin != kmax:
                dk = (kmax-kmin)
                res = (dk+1)/puiss
                if dk%puiss != 0:
                    errors+=[n,"PointRangeDonor for GC %s of zone %s is not multigrid of level %d in direction k."%(n[0],z[0],multigrid)]
                elif res < sizemin: errors+=[n,"PointRangeDonor for GC %s of zone %s: not enough points on coarse grid in direction i."%(n[0],z[0])]
    return errors

#==============================================================================
# check if the tree is compatible with multigrid (zones, BC and connectivities)
# level is the MG level that must be ensured: N = 2^level+1
#==============================================================================
def checkMultigrid(t, level=1, nbMinCoarseB=5, nbMinCoarseW=3):
    """Check multigrid validity (zones, BC and connectivities)."""
    errors = []
    if level == 0: return errors
    puiss = 2**(level)

    for z in Internal.getZones(t):
        # check Dims
        dims = Internal.getZoneDim(z)
        if dims[0] == 'Structured':
            ni = dims[1]; nj = dims[2]; nk = dims[3]
            res = (ni-1)/puiss
            if (ni-1)%puiss != 0:
                errors+=[z,"Zone %s is not multigrid of level %d in direction i."%(z[0],level)]
            elif res < nbMinCoarseB:
                errors+=[z,"Zone %s: not enough points on coarse grid for level %d in direction i."%(z[0],level)]
            res = (nj-1)/puiss
            if (nj-1)%puiss != 0: errors+=[z,"Zone %s is not multigrid of level %d in direction j."%(z[0],level)]
            elif res < nbMinCoarseB:
                errors+=[z,"Zone %s: not enough points on coarse grid for level %d in direction j."%(z[0],level)]
            res = (nk-1)/puiss
            if (nk-1)%puiss != 0: errors+=[z,"Zone %s is not multigrid of level %d in direction k."%(z[0],level)]
            elif res < nbMinCoarseB:
                errors+=[z,"Zone %s: not enough points on coarse grid for level %d in direction k."%(z[0],level)]

            # check BC ranges (receptors)
            errors += checkMGForBCRanges(z,'BC_t',level,nbMinCoarseW)
            errors += checkMGForBCRanges(z,'GridConnectivity1to1_t',level,nbMinCoarseW)
            errors += checkMGForBCRanges(z,'GridConnectivity_t',level,nbMinCoarseW)
            # check BC ranges (donors)
            errors += checkMGForDonorBCRanges(z,'GridConnectivity1to1_t',level,nbMinCoarseW)
            errors += checkMGForDonorBCRanges(z,'GridConnectivity_t',level,nbMinCoarseW)
    return errors

#=============================================================================
# Check if the number of points of a zone does not exceed sizeMax
#=============================================================================
def checkSize(t, sizeMax=100000000):
    """Check if the number of points of zones dont exceed sizeMax."""
    errors = []
    for z in Internal.getZones(t):
        dims = Internal.getZoneDim(z)
        if dims[0] == 'Structured':
            npts = dims[1]*dims[2]*dims[3]
        else: npts = dims[1]
        if npts > sizeMax: errors += [z,"Zone %s exceeds the maximum number of points (Npts=%d)."%(z[0],npts)]
    return errors

#==============================================================================
# Verifie que le type des noeuds est dans CGNSTypes
#==============================================================================
def checkCGNSType(node):
    """Check if all node are of a valid CGNS type."""
    errors = []
    isStd = Internal.isStdNode(node)
    if isStd >= 0:
        for c in node[isStd:]: checkCGNSType__(c, node, errors)
    else: checkCGNSType__(node, node, errors)
    return errors

def checkCGNSType__(node, parent, errors):
    ntype = node[3]
    if ntype not in CGNSTypes:
        errors += [node, parent, 'Unknown CGNS type %s for node %s.\n'%(ntype, node[0])]
    sons = node[2]
    for n in sons: checkCGNSType__(n, node, errors)

#==============================================================================
# Delete les noeuds de type non valide
#==============================================================================
def _correctCGNSType(t):
    """Delete nodes of invalid CGNS types."""
    errors = checkCGNSType(t)
    le = len(errors)//3
    for e in range(le):
        node = errors[3*e]
        parent = errors[3*e+1]
        c = Internal.getNodePosition(node, parent)
        del parent[2][c]
    return None

#==============================================================================
# Check element nodes dans t
# Verifie:
# si une zone a NGON+PE et pas de NFace
# les vertex references existent pour element et ngon
# pas de faces negatives
#==============================================================================
def checkElementNodes(t):
    errors = []
    bases = Internal.getBases(t)
    for b in bases:
        zones = Internal.getZones(b)
        for z in zones:
            connects = Internal.getElementNodes(z)
            iBEs = []; iNGons = []; iNFaces = []; i = 0
            for c in connects:
                ctype = c[1][0]
                if ctype == 22: iNGons.append(i)
                elif ctype == 23: iNFaces.append(i)
                else: iBEs.append(i)
                i += 1

            if len(iNGons) > 1:
                print("Warning: CheckPyTree: multiple NGON connectivity for zone %s. Please merge NGON connectivities."%z[0])
            if len(iNFaces) > 1:
                print("Warning: CheckPyTree: multiple NFACE connectivity for zone %s. Please merge NFACE connectivities."%z[0])
            if len(iNGons) == 1 and len(iNFaces) == 0:
                errors += [z, b, 'NFace is missing for zone %s.'%z[0]]

            #if len(BE) > 1:
            #    errors += [z, b, 'Multiple BE connectivity for zone %s.'%z[0]]

            for iBE in iBEs: # check existing vertices
                c = Internal.getNodeFromName1(connects[iBE], 'ElementConnectivity')
                if c[1] is None:
                    print("Warning: CheckPyTree: ElementConnectivity is None (may not be loaded).")
                elif c[1].size > 0:
                    minv = numpy.min(c[1]); maxv = numpy.max(c[1])
                    npts = C.getNPts(z)
                    if minv < 1 or maxv > npts:
                        errors += [z, b, 'Element connectivity referenced unexisting vertices in zone %s.'%z[0]]

            for iNFace in iNFaces: # check positive faces
                c = Internal.getNodeFromName1(connects[iNFace], 'ElementConnectivity')
                if c[1] is None: print("Warning: CheckPyTree: ElementConnectivity is None (may not be loaded).")
                else:
                    minv = numpy.min(c[1])
                    if minv < 1: errors += [c, connects[iNFace], 'Negative NFace index']

            for iNGon in iNGons: # check existing vertices
                c = Internal.getNodeFromName1(connects[iNGon], 'ElementConnectivity')
                if c[1] is None:
                    print("Warning: CheckPyTree: ElementConnectivity is None (may not be loaded).")
                else:
                    minv = numpy.min(c[1]); maxv = numpy.max(c[1])
                    npts = C.getNPts(z)
                    if minv < 1 or maxv > npts:
                        errors += [z, b, 'NGON connectivity referenced unexisting vertices in zone %s.'%z[0]]
    return errors

#==============================================================================
# Fait un break connectivity pour les BE multiples
# Ajoute le noeud Face si on a un PE et pas de NFace
# Fait abs(index) pour les NFaces
#==============================================================================
def _correctElementNodes(t):
    _correctBCElementNodes(t)
    _cleanBEConnect(t)
    errors = checkElementNodes(t)
    le = len(errors)//3
    for e in range(le):
        zone = errors[3*e]
        parent = errors[3*e+1]
        msg = errors[3*e+2]
        if msg[0:8] == 'Negative':
            zone[1] = numpy.absolute(zone[1])
        if msg[0:8] == 'Multiple':
            zones = C.breakConnectivity(zone)
            c = Internal.getNodePosition(zone, parent)
            parent[2][c] = zones[0]; parent[2] += zones[1:]
        elif msg[0:5] == 'NFace':
            # Look for PE
            PE = Internal.getNodeFromName2(zone, 'ParentElements')
            if PE is not None:
                Internal._adaptPE2NFace(zone, remove=False)
    return None

#===============================================================================
# Corrige des boundary connectivity qui sont a zero (GE[1][1])
#===============================================================================
def _correctBCElementNodes(t):
    #_correctBCPL2ER(t)

    zones = Internal.getZones(t)
    for z in zones:
        GEl = Internal.getNodesFromType1(z, 'Elements_t')
        maxDim = 0 # dimension elements non BC
        minDim = 3
        for GE in GEl:
            itype = GE[1][0]
            ibc = GE[1][1]
            if ibc == 0:
                if itype >= 3 and itype == 4: maxDim = max(maxDim, 1); minDim = min(minDim, 1)
                if itype >= 5 and itype <= 9: maxDim = max(maxDim, 2); minDim = min(minDim, 2)
                if itype >= 10: maxDim = max(maxDim, 3); minDim = min(minDim, 3)
        if minDim != maxDim:
            for GE in GEl:
                itype = GE[1][0]
                ibc = GE[1][1]
                if ibc == 0:
                    if itype >= 3 and itype == 4 and maxDim > 1: GE[1][1] = 7
                    if itype >= 5 and itype <= 9 and maxDim > 2: GE[1][1] = 7
    return None

#==========================================================================================
#  In a Basic Element mesh, retrieve BCs with ElementRange refering a surface connectivity
#  1/ if BC is defined by a PointList= a set of vertices
#  2/ or if BC is defined by a PointList= a set of face centers
#==========================================================================================
def _correctBCPL2ER(t):
    import Post.PyTree as P
    for z in Internal.getZones(t):
        zdim = Internal.getZoneDim(z)
        if zdim[0] == 'Unstructured':
            ztype = zdim[3]
            if ztype != 'NGON':
                for zbc in Internal.getNodesFromType2(z, 'BC_t'):
                    bndName = Internal.getName(zbc)
                    #bndType = Internal.getValue(zbc)
                    gcl = Internal.getNodeFromType1(zbc, 'GridLocation_t')
                    PL = Internal.getNodeFromName1(zbc, 'PointList')
                    if PL is not None:
                        if gcl is None or Internal.getValue(gcl)=='Vertex':
                            C._initVars(z, 'tag', 0.)
                            for ind in Internal.getValue(PL)[0]:
                                C.setValue(z, 'tag', ind, 1.)
                            bc = P.exteriorFaces(z)
                            C._rmVars(z, ["tag"])
                            bc = P.selectCells(bc,"{tag}>0",strict=1)
                            Internal._rmNodesFromName(zbc, "PointList")
                            if bc != []:
                                bc[0] = C.getZoneName(bndName)
                                C._mergeConnectivity(z, bc, boundary=1)
                                bcn = Internal.getNodeFromName1(z, bc[0])
                                bcnr = Internal.getNodeFromName1(bcn, 'ElementRange')
                                ER = [bcnr[1][0], bcnr[1][1]]

                        elif Internal.getValue(gcl) == 'FaceCenter':
                            bc_gcname = 'Elements_%s'%(Internal.getName(zbc))
                            bc_gcnode = Internal.getNodeFromName(z, bc_gcname)
                            if bc_gcnode is not None:
                                ER = Internal.getNodeFromName(bc_gcnode, 'ElementRange')
                                if ER is not None: ER = Internal.getValue(ER)

                        if ER is not None:
                            r = numpy.empty((1,2), dtype=Internal.E_NpyInt, order='F')
                            r[0,0] = ER[0]
                            r[0,1] = ER[1]
                            zbc[2].append([Internal.__ELEMENTRANGE__, r, [], 'IndexRange_t'])
                            Internal._rmNodesFromName(zbc, 'GridLocation')
                            Internal._rmNodesFromName(zbc, 'PointList')
    return None

#===============================================================================
#  Dans un maillage NGON, enleve les connectivites BE
#===============================================================================
def _cleanBEConnect(t):
    for z in Internal.getZones(t):
        zdim = Internal.getZoneDim(z)
        if zdim[0] == 'Unstructured':
            if zdim[3] == 'NGON':
                for elt_t in Internal.getNodesFromType1(z, 'Elements_t'):
                    elt_no = Internal.getValue(elt_t)[0]
                    if elt_no != 22 and elt_no != 23:
                        Internal._rmNodesFromName(z, elt_t[0])
    return None

# Dans un ParentElement, shift les elements si pas deja le cas
# shift=1: add nface shift if possible
# shift=-1: sub nface shift if possible
def _shiftParentElements(t, shift=1):
    for z in Internal.getZones(t):
        ngon = Internal.getNGonNode(z)
        if ngon is not None:
            PE = Internal.getNodeFromName1(ngon, 'ParentElements')
            if PE is not None and PE[1] is not None:
                n = PE[1]
                if n is not None:
                    np = n[n > 0]
                    emin = numpy.min(np)
                    if emin == 1:
                        if shift == 1:
                            nf = n.size//2
                            np = n+nf
                            PE[1] = numpy.where(np==nf, 0, np)
                    else:
                        if shift == -1:
                            nf = n.size//2
                            np = n-nf
                            PE[1] = numpy.where(np==-nf, 0, np)

    return None

#===============================================================================
# Check non CGNS varnames in FlowSolution_t and BCDataSet
#===============================================================================
def checkCGNSVarNames(t):
    errors = []
    zones = Internal.getZones(t)
    for z in zones:
        # Change on container
        cont = Internal.getNodesFromType1(z, 'FlowSolution_t')
        for c in cont:
            datas = Internal.getNodesFromType1(c, 'DataArray_t')
            for d in datas:
                n = Internal.getCGNSName(d[0])
                if n != d[0]: errors += [d, c, '%s not a valid CGNS variable name for zone %s.'%(d[0],z[0])]
        # Change on BCDataSet (if any)
        bcs = Internal.getNodesFromType2(z, 'BC_t')
        for b in bcs:
            datas = Internal.getBCDataSet(z, b)
            for d in datas:
                n = Internal.getCGNSName(d[0])
                if n != d[0]: errors += [d, b, '%s not a valid CGNS variable name for BCDataSet in zone %s.'%(d[0],z[0])]
    return errors

def _correctCGNSVarNames(t):
    errors = checkCGNSVarNames(t)
    for e in range(len(errors)//3):
        node = errors[3*e]
        n = Internal.getCGNSName(node[0])
        node[0] = n
    return None

#===============================================================================
# Check if Coordinates are in fields, change name with + F
#===============================================================================
def checkCoordinatesInFields(t):
    errors = []
    zones = Internal.getZones(t)
    for z in zones:
        # Check GridCoordinates size
        npts = C.getNPts(z)
        cont = Internal.getNodesFromType1(z, 'GridCoordinates_t')
        for c in cont:
            datas = Internal.getNodesFromType1(c, 'DataArray_t')
            for d in datas:
                if d[1] is None: print("CheckPyTree: coordinates is None (may not be loaded).")
                elif d[1].shape != npts:
                    errors += [d, c, 'zone %s has coordinates of wrong size.'%(z[0])]

        # Check if coordinates are in a FlowSolution_t
        cont = Internal.getNodesFromType1(z, 'FlowSolution_t')
        for c in cont:
            datas = Internal.getNodesFromType1(c, 'DataArray_t')
            for d in datas:
                if d[0] == 'CoordinateX' or d[0] == 'CoordinateY' or d[0] == 'CoordinateZ':
                    errors += [d, c, 'zone %s contains Coordinates in Field.'%(z[0])]
        # Change on BCDataSet (if any)
        bcs = Internal.getNodesFromType2(z, 'BC_t')
        for b in bcs:
            datas = Internal.getBCDataSet(z, b)
            for d in datas:
                if d[0] == 'CoordinateX' or d[0] == 'CoordinateY' or d[0] == 'CoordinateZ':
                    errors += [d, b, 'zone %s contains Coordinates in BCDataSet.'%(z[0])]
    return errors

def _correctCoordinatesInFields(t):
    errors = checkCoordinatesInFields(t)
    for e in range(len(errors)//3):
        node = errors[3*e]
        n = node[0]+'F'
        node[0] = n
    return None

#===============================================================================
# Check field conformity and location
#===============================================================================
def checkFieldConformity(t):
    errors = []
    zones = Internal.getZones(t)
    for z in zones:
        npts = C.getNPts(z)
        ncells = C.getNCells(z)
        # Change on container
        cont = Internal.getNodesFromType1(z, 'FlowSolution_t')
        for c in cont:
            loc = Internal.getNodeFromType1(c, 'GridLocation_t')
            if loc is None: loc = npts
            elif Internal.getValue(loc) == 'Centers': loc = ncells
            elif Internal.getValue(loc) == 'Vertex': loc = npts
            datas = Internal.getNodesFromType1(c, 'DataArray_t')
            for d in datas:
                if d[1].size != loc:
                    errors += [d, c, '%s has not the right size for field in zone %s.'%(d[0],z[0])]
    return errors

def _correctFieldConformity(t):
    errors = checkFieldConformity(t)
    for e in range(len(errors)//3):
        node = errors[3*e]
        parent = errors[3*e+1]
        c = Internal.getNodePosition(node, parent)
        del parent[2][c]
    return None

#===============================================================================
# Check NAN in fields
#===============================================================================
def checkNAN(t):
    errors = []
    zones = Internal.getZones(t)
    for z in zones:
        vars = C.getVarNames(z)[0]
        for v in vars:
            isFinite = C.isFinite(z, v)
            if not isFinite:
                d = Internal.getZoneName2(z, v)
                errors += [d, v, 'Field %s of zone %s has NAN.'%(v,z[0])]
    return errors

def _correctNAN(t):
    errors = checkNAN(t)
    for e in range(len(errors)//3):
        node = errors[3*e]
        field = errors[3*e+1]
        C._setNANValuesAt(node, field, 0.)
