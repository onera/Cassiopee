"""Conversion module for Cassiopee package.
"""
#from numpy import *
__version__ = '4.1'
__author__ = "Stephanie Peron, Christophe Benoit, Gaelle Jeanfaivre, Pascal Raud, Benoit Rodriguez, Simon Verley, Bruno Maugars, Thomas Renaud"
#
# Python Interface for conversion between array / file / CGNS
#
import numpy
import os.path
try: from . import converter
except: import converter
import KCore

# INT size for numpys connectivities
from KCore.Dist import EDOUBLEINT
if EDOUBLEINT: E_NpyInt = numpy.int64
else: E_NpyInt = numpy.int32

__all__ = [
    'array', 'getApi', 'addVars', '_addVars', 'addVars2',
    'center2ExtCenter', 'center2Node', 'conformizeNGon',
    'convertArray2Hexa', 'convertArray2NGon', 'convertArray2Node',
    'convertArray2Tetra', 'convertBAR2Struct', 'convertTri2Quad',
    'convertStrand2Penta', 'convertPenta2Strand',
    'convertArrays2File', 'convertFile2Arrays', 'copy',
    'createGlobalHook', 'createHook',
    'createGlobalIndex', '_createGlobalIndex',
    'recoverGlobalIndex', '_recoverGlobalIndex',
    'createSockets', 'diffArrays', 'diffArraysGeom', 'isFinite',
    'setNANValuesAt', 'extCenter2Node', 'extractVars', 'getIndexField',
    'freeHook', 'getArgMax', 'getArgMin', 'getMaxValue', 'getMeanRangeValue',
    'getMeanValue', 'getMinValue', 'getNCells', 'getNPts', 'getValue',
    'getVarNames', 'identifyElements', 'identifyFaces', 'identifyNodes',
    'identifySolutions', 'initVars', '_initVars', 'isNamePresent', 'listen',
    'magnitude', 'nearestElements', 'nearestFaces', 'nearestNodes',
    'node2Center', 'node2ExtCenter', 'normL0', 'normL2',
    'normalize', '_normalize', 'randomizeVar', 'rmVars',
    'send', 'setPartialFields', 'setValue', 'addGhostCellsNGon',
    'checkFileType', 'convertHO2LO', 'convertLO2HO', 'convertExt2Format__',
    'mergeConnectivity','adaptSurfaceNGon',
    '_signNGonFaces', '_unsignNGonFaces', 'makeParentElements'
]

# -- Create an array --
# Les champs sont mis a zero, sauf si pour les champs cellN et cellNF
def array(vars, n1, n2, sub, api=1):
    """Create a structured or unstructured array.
    Usage: array(vars, ni, nj, nk)
    or array(vars, npoints, nelts, eltType)"""
    if isinstance(sub, str):
        return arrayNS(vars, n1, n2, sub, api)
    else:
        return arrayS(vars, n1, n2, sub, api)

# -- Create a structured array --
def arrayS(vars, ni, nj, nk, api=1):
    """Create a structured array.
    Usage: array(vars, ni, nj, nk)"""
    ni = int(ni); nj = int(nj); nk = int(nk)
    vars = vars.replace(' ', '')
    if vars.startswith(',') or vars.endswith(','):
        print("Warning: array: your var string is suspicious.")
    vl = vars.split(','); v = len(vl)
    if api == 1:
        a = numpy.zeros((v, ni*nj*nk), numpy.float64)
        for i in range(v):
            if vl[i] == 'cellN' or vl[i] == 'cellNF': a[i,:] = 1.
    else:
        a = []
        for i in range(v):
            if vl[i] == 'cellN' or vl[i] == 'cellNF':
                a.append(numpy.ones((ni*nj*nk), numpy.float64))
            else:
                a.append(numpy.zeros((ni*nj*nk), numpy.float64))
    return [vars, a, ni, nj, nk]

# -- Create an unstructured array --
def arrayNS(vars, npoints, nelts, eltType, api=1):
    """Create a unstructured array.
    Usage: array(vars, npoints, nelts, eltType)"""
    npoints = int(npoints); nelts = int(nelts)
    vars = vars.replace(' ', '')
    if vars.startswith(',') or vars.endswith(','):
        print("Warning: array: your var string is suspicious.")
    vl = vars.split(','); v = len(vl)
    eltTypeN = eltType.replace('*', '')
    if eltTypeN == 'NODE': nt = 1
    elif eltTypeN == 'BAR': nt = 2
    elif eltTypeN == 'TRI': nt = 3
    elif eltTypeN == 'QUAD': nt = 4
    elif eltTypeN == 'TETRA': nt = 4
    elif eltTypeN == 'PYRA': nt = 5
    elif eltTypeN == 'PENTA': nt = 6
    elif eltTypeN == 'HEXA': nt = 8
    else:
        raise ValueError("arrayNS: this function doesnt work for %s."%eltType)

    if eltType[-1] == '*':
        if api == 1:
            a = numpy.zeros((v, nelts), numpy.float64)
            for i in range(v):
                if vl[i] == 'cellN' or vl[i] == 'cellNF': a[i,:] = 1.
        else:
            a = []
            for i in range(v):
                if vl[i] == 'cellN' or vl[i] == 'cellNF':
                    a.append(numpy.ones((nelts), numpy.float64))
                else:
                    a.append(numpy.zeros((nelts), numpy.float64))
    else:
        if api == 1:
            a = numpy.zeros((v, npoints), numpy.float64)
            for i in range(v):
                if vl[i] == 'cellN' or vl[i] == 'cellNF': a[i,:] = 1.
        else:
            a = []
            for i in range(v):
                if vl[i] == 'cellN' or vl[i] == 'cellNF':
                    a.append(numpy.ones((npoints), numpy.float64))
                else:
                    a.append(numpy.zeros((npoints), numpy.float64))

    if api == 1: c = numpy.ones((nt, nelts), E_NpyInt)
    else: c = [numpy.ones((nelts, nt), E_NpyInt)]

    return [vars, a, c, eltType]

# -- get array api (1,2,3) --
# IN: array
# OUT: -1 (not an array), 1,2,3 following api
# OUT: 1 (array1)
# OUT: 2 (array2) or array3 (structured, mono element because identical)
# OUT: 3 (array3) (ME, NGONv4)
def getApi(a):
    if len(a) == 5: # structure
        if not isinstance(a[0], str): return -1
        if isinstance(a[1], numpy.ndarray): return 1
        if isinstance(a[1], list): return 2 # compat. 3
        return -1
    elif len(a) == 4: # non structure
        if not isinstance(a[0], str): return -1
        if not isinstance(a[3], str): return -1
        eltType = a[3]
        if isinstance(a[1], numpy.ndarray) and isinstance(a[2], numpy.ndarray): return 1
        if isinstance(a[1], list):
            if not isinstance(a[2], list): return -1
            if eltType == "NGON" or eltType == "NGON*":
                if len(a[2]) == 2: return 2
                if len(a[2]) == 4 or len(a[2]) == 5:
                    if a[2][3][-1] == a[2][1].size: return 3
                    else: return 2
                else: return -1
            else: # ELTS
                if len(a[2]) == 1: return 2 # compat. 3
                else: return 3
            if isinstance(a[2], list): return 3
            return -1
        return -1
    else:
        return -1 # not a valid array

# -- Get the number of points of an array --
def getNPts(a):
    """Return the total number of points."""
    if isinstance(a[0], list):
        b = 0
        for i in a: b += getNPts__(i)
        return b
    else: return getNPts__(a)

def getNPts__(a):
    npts = 0
    if len(a) == 5: # structure
        npts = a[2]*a[3]*a[4]
    elif len(a) == 4: # non structure
        if a[3][-1] != '*':
            if isinstance(a[1], list): npts = a[1][0].size
            else: npts = a[1].shape[1]
        else: npts = 0 # don't know
    return npts

# -- Get the number of cells of an array --
def getNCells(a):
    """Return the total number of cells."""
    if isinstance(a[0], list):
        b = 0
        for i in a: b += getNCells__(i)
        return b
    else: return getNCells__(a)

def getNCells__(a):
    ncells = 0
    if len(a) == 5: # structure
        nic = a[2]-1; njc = a[3]-1; nkc = a[4]-1
        nic = max(nic, 1); njc = max(njc, 1); nkc = max(nkc, 1)
        ncells = nic*njc*nkc
    elif len(a) == 4: # non structure
        if isinstance(a[1], list):
            if a[3] == 'NGON' or a[3] == 'NGON*':
                ncells = a[2][3].size
            else: ncells = a[2][0].shape[0]
        else:
            if a[3] == 'NGON' or a[3] == 'NGON*':
                offset = a[2][0,1]
                ncells = a[2][0,offset+2]
            else: ncells = a[2].shape[1]
    return ncells

# -- Copy an array --
def copy(array):
    """Copy  an  array.
    Usage: a2 = copy(a1)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(converter.copy(i))
        return b
    else:
        return converter.copy(array)

def setPartialFields(a, ad, listIndices):
    """Set values defined in a 1D array for points in listIndices
    Usage: b=setPartialFields(a,ad,listIndices)"""
    return converter.setPartialFields(a,ad,listIndices)

# -- rmVars --
def rmVars(a, var):
    """Remove variables."""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(rmVars__(i, var))
        return b
    else:
        return rmVars__(a, var)

def rmVars__(a, var):
    varstring = a[0]
    vs = varstring.split(',')
    eVars = []
    if isinstance(var, list):
        for v in vs:
            present = False
            for w in var:
                if v == w: present = True; break
            if not present: eVars.append(v)
    else:
        for v in vs:
            if v != var: eVars.append(v)
    return extractVars(a, eVars)

# -- Extract a list of variables from an array --
# only for array1
def extractVars(array, vars):
    """Extract variables from a.
    Usage: a = extractVars(array, ['x','ro'])"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            varsi = []
            if isinstance(vars, list):
                for j in vars:
                    p = KCore.isNamePresent(i, j)
                    varsi.append(p+1)
            else:
                p = KCore.isNamePresent(i, vars)
                varsi.append(p+1)
            b.append(converter.extractVars(i, varsi))
        return b
    else:
        varsi = []
        for i in vars:
            p = KCore.isNamePresent(array, i)
            varsi.append(p+1)
        return converter.extractVars(array, varsi)

# -- Add variables --
def addVars(array, add=None):
    """Add variables to an array.
    Usage: a = addVars(array, 'varString')
    or: a = addVars([a1, a2, a3])"""
    if add is not None:
        if isinstance(array[0], list):
            b = []
            for i in array: b.append(addVar(i, add))
            return b
        else:
            return addVar(array, add)
    else:
        if isinstance(array[0][0], list):
            l = 0
            for i in array:
                if l == 0: l = len(i)
                elif len(i) != l:
                    raise ValueError("addVars: can not add variables.")
            c = 0; b = []
            for i in array[0]:
                l = []
                for i in array: l.append(i[c])
                c += 1
                b.append(addVars1(l))
            return b
        else:
            return addVars1(array)

# -- Add variables defined by a string --
def addVar(array, vars):
    """Add variables to an array.
    Usage: a = addVar(array, vars)"""
    if isinstance(vars, list):
        a = array
        for i in vars: a = converter.addVar(a, i) # super cher
        return a
    else:
        return converter.addVar(array, vars)

# -- Add variables defined by arrays --
def addVars1(arrays):
    """Put all variables in an array. arrays must be a list of arrays
    defined on the same grid.
    Usage: a = addVars(arrays)"""
    return converter.addVars(arrays)

# -- Nouvelle version (array1/array2) --
def addVars2(a, varNames=None):
    if varNames is None: # concat une liste
        b = list(a)
        b[0] = copy(b[0])
    else: b = copy(a)
    _addVars(b, varNames)
    if varNames is None: return b[0]
    else: return b

def _addVars(a, varNames=None):
    if varNames is None and isinstance(a[0][0], list):
        c = 0
        for i in a[0]:
            _addVars2__([i, a[1][c]]); c += 1
    elif varNames is None: _addVars2__(a)
    else:
        if isinstance(a[0], list):
            for i in a: _addVars__(i, varNames)
        else: _addVars__(a, varNames)
    return None

# add a list of varnames
def _addVars__(a, varNames):
    if isinstance(varNames, str): varNames = [varNames]
    if isinstance(a[1], list): # array2
        for v in varNames:
            if KCore.isNamePresent(a, v) == -1:
                a[0] += ','+v
                s = a[1][0].shape
                a[1].append(numpy.zeros(s, dtype=numpy.float64, order='F'))
    else: # array1
        lg = 0
        for v in varNames:
            if KCore.isNamePresent(a, v) == -1:
                a[0] += ','+v
                lg += 1
        s = a[1].shape[1]
        nfld = a[1].shape[0]
        n = numpy.zeros((nfld+lg, s), dtype=numpy.float64)
        for i in range(nfld): n[i,:] = a[1][i,:]
        a[1] = n
    return None

# a est une liste d'arrays qu'il faut concatener dans le premier
# modifie le premier array in place
def _addVars2__(a):
    if isinstance(a[0][1], list): # array2/3
        a0 = a[0]
        rest = a[1:]
        for i in rest: a0[0] += ','+i[0]
        for i in rest: a0[1] += i[1]
    else: # array1
        a0 = a[0]
        rest = a[1:]
        for i in rest: a0[0] += ','+i[0]
        nfld = a0[1].shape[0]
        for i in rest: nfld += i[1].shape[0]
        n = numpy.empty((nfld, a0[1].shape[1]), dtype=numpy.float64)
        nfld = a0[1].shape[0]
        for j in range(nfld):
            n[j,:] = a0[1][j,:]
        for i in rest:
            nfld2 = i[1].shape[0]
            for j in range(nfld2):
                n[nfld+j,:] = i[1][j,:]
            nfld += nfld2
        a0[1] = n
    return None

# -- randomize a variable --
def randomizeVar(array, var, deltaMin, deltaMax):
    """Randomize a field defined by var within a range [a-deltaMin, a+deltaMax].
    Usage: a = randomizeVar(array, var, deltaMin, deltaMax)"""
    if isinstance(array[0], list):
        out = []
        for i in array: out.append(converter.randomizeVar(i, var, deltaMin, deltaMax))
        return out
    else:
        return converter.randomizeVar(array, var, deltaMin, deltaMax)

# -- Init variables --
def initVars(a, var, v1=None, v2=None, mode=0, isVectorized=False):
    """Initialize a variable by a value, a function or a formula."""
    b = copy(a)
    _initVars(b, var, v1, v2, mode, isVectorized)
    return b

def _initVars(a, var, v1=None, v2=None, mode=0, isVectorized=False):
    if isinstance(a[0], list):
        for i in a: _initVars__(i, var, v1, v2, mode, isVectorized)
    else: _initVars__(a, var, v1, v2, mode, isVectorized)
    return None

def _initVars__(a, var, v1, v2, mode=0, isVectorized=False):
    if v1 is None:
        if mode == 0: _initVarByEq__(a, var) # numpy eval
        else: _initVarByEq2__(a, var) # expression eval
    elif callable(v1):
        _initVarByFunction__(a, var, v1, v2, isVectorized)
    else:
        _initVarByConst__(a, var, v1)
    return None

def _initVarByConst__(a, var, val):
    # Init one or several vars by a constant value.
    if not isinstance(var, list): var = [var]
    posvars = []
    for v in var:
        posvar = KCore.isNamePresent(a, v)
        if posvar == -1:
            _addVars(a, var)
            posvar = KCore.isNamePresent(a, v)
        posvars.append(posvar)

    if not isinstance(a[1], list): # array1
        for posvar in posvars: a[1][posvar,:] = val
    else:
        for posvar in posvars: a[1][posvar][:] = val
    return None

def _initVarByFunction__(a, var, F, fargs=[], isVectorized=False):
    # Init one or several vars with a function F.
    isResMultivars = lambda res: isinstance(res, tuple)
    isResSinglevar = lambda res: not isResMultivars(res)

    errorMsgVar = "The number of arguments returned by the function ({}) is "\
                  "different from the number of variables to initialise ({})."
    errorMsgName = "_initVarByFunction__: can't find variables {} in array."

    def isResValid(funcRes, numVarsToInit, errorMsgVar):
        singlevar = (numVarsToInit == 1)
        if singlevar:
            if isResMultivars(funcRes):
                raise IndexError(errorMsgVar.format(len(funcRes), 1))
        elif isResSinglevar(funcRes):
            raise IndexError(errorMsgVar.format(1, numVarsToInit))

    if not isinstance(var, list): var = [var]

    # Find position of vars
    posvars = []
    for v in var:
        posvar = KCore.isNamePresent(a, v)
        if posvar == -1:
            _addVars(a, var)
            posvar = KCore.isNamePresent(a, v)
        posvars.append(posvar)
    ninitvars = len(posvars)
    singlevar = (ninitvars == 1)

    # Find positions of the function's arguments
    posargs = numpy.array([KCore.isNamePresent(a, i) for i in fargs], dtype=int)
    posvarsNotFound = (posargs == -1).nonzero()[0]
    if posvarsNotFound.size:
        raise NameError(errorMsgName.format(', '.join(fargs[i] for i in posvarsNotFound)))

    # Evaluate variable(s)
    # Manage the following cases:
    #   - is the node a[1] a numpy.ndarray (array1) or a list (array3)
    #   - is the function F vectorized
    #   - does the function F have arguments
    #   - is there one or several lhs variables to initialise
    if isVectorized:
        # F is vectorized: manipulate arrays
        if isinstance(a[1], list): # array3
            if len(fargs):
                res = F(*[a[1][posarg][:] for posarg in posargs])
                if singlevar:
                    # Case A1
                    if isResMultivars(res):
                        raise IndexError(errorMsgVar.format(len(res), 1))
                    a[1][posvars[0]][:] = res
                else:
                    # Case A2
                    if isResSinglevar(res):
                        raise IndexError(errorMsgVar.format(1, ninitvars))
                    for i, posvar in enumerate(posvars):
                        a[1][posvar][:] = res[i]
            else:
                res = F()
                if singlevar:
                    # Case B1
                    if isResMultivars(res):
                        raise IndexError(errorMsgVar.format(len(res), 1))
                    a[1][posvars[0]][:] = res
                else:
                    # Case B2
                    if isResSinglevar(res):
                        raise IndexError(errorMsgVar.format(1, ninitvars))
                    for i, posvar in enumerate(posvars):
                        a[1][posvar][:] = res[i]
        else: # array1
            if len(fargs):
                # Case C
                res = F(*[a[1][posarg,:] for posarg in posargs])
                isResValid(res, ninitvars, errorMsgVar)
                a[1][posvars,:] = res
            else:
                res = F()
                isResValid(res, ninitvars, errorMsgVar)
                if singlevar:
                    # Case D1
                    a[1][posvars[0],:] = res
                else:
                    # Case D2
                    for i, posvar in enumerate(posvars):
                        a[1][posvar,:] = res[i]
    else:
        # F is not vectorized: loop over all elements of the list/numpy.ndarray
        if isinstance(a[1], list):
            if singlevar:
                varFlatten = a[1][posvars[0]].ravel(order='K')
                if len(fargs):
                    # Case E1a
                    nelems = a[1][posvars[0]].size
                    fargsFlatten = numpy.array([a[1][posarg].ravel(order='K') for posarg in posargs])
                    res = F(*fargsFlatten[:,0])
                    isResValid(res, ninitvars, errorMsgVar)
                    varFlatten[0] = res
                    for j in range(1, nelems):
                        res = F(*fargsFlatten[:,j])
                        varFlatten[j] = res
                else:
                    # Case E1b
                    res = F()
                    isResValid(res, ninitvars, errorMsgVar)
                    varFlatten[:] = res

            else:
                if len(fargs):
                    # Case E2a
                    nelems = a[1][posvars[0]].size
                    varsFlatten = [a[1][posvar].ravel(order='K') for posvar in posvars]
                    fargsFlatten = numpy.array([a[1][posarg].ravel(order='K') for posarg in posargs])
                    res = F(*fargsFlatten[:,0])
                    isResValid(res, ninitvars, errorMsgVar)
                    for j in range(nelems):
                        res = F(*fargsFlatten[:,j])
                        for i, posvar in enumerate(posvars):
                            varsFlatten[i][j] = res[i]
                else:
                    # Case E2b
                    res = F()
                    isResValid(res, ninitvars, errorMsgVar)
                    for i, posvar in enumerate(posvars):
                        varFlatten = a[1][posvar].ravel(order='K')
                        varFlatten[:] = res[i]
        else:
            if len(fargs):
                nelems = a[1].shape[1]
                res = F(*[a[1][posarg,0] for posarg in posargs])
                isResValid(res, ninitvars, errorMsgVar)
                if singlevar:
                    # Case F1
                    posvar = posvars[0]
                    a[1][posvar,0] = res
                    for j in range(1, nelems):
                        a[1][posvar,j] = F(*[a[1][posarg,j] for posarg in posargs])
                else:
                    # Case F2
                    for j in range(nelems):
                        res = F(*[a[1][posarg,j] for posarg in posargs])
                        for i, posvar in enumerate(posvars):
                            a[1][posvar,j] = res[i]
            else:
                res = F()
                isResValid(res, ninitvars, errorMsgVar)
                if singlevar:
                    # Case G1
                    posvar = posvars[0]
                    a[1][posvar,:] = res
                else:
                    # Case G2
                    for i, posvar in enumerate(posvars):
                        a[1][posvar,:] = res[i]
    return None

# Initialisation par une formule par numpy
def _initVarByEq__(a, eq):
    import re
    varstring = a[0]
    vars = varstring.split(',')

    replacements = {
        'minimum(': 'numpy.minimum(',
        'maximum(': 'numpy.maximum(',
        'cos(': 'numpy.cos(',
        'cosh(': 'numpy.cosh(',
        'sin(': 'numpy.sin(',
        'sinh(': 'numpy.sinh(',
        'sqrt(': 'numpy.sqrt(',
        'log(': 'numpy.log(',
        'tan(': 'numpy.tan(',
        'atan(': 'numpy.arctan(',
        'arctan(': 'numpy.arctan(',
        'arctan2(': 'numpy.arctan2(',
        'exp(': 'numpy.exp(',
        'degrees(': 'numpy.degrees(',
        'logical_and(': 'numpy.logical_and(',
        'pi': 'numpy.pi'
    }

    eq = eq.replace('centers:', '')
    eq = eq.replace('nodes:', '')
    # Only substitute dict key by its value if the key is a prefix and if it is
    # not preceded by a full stop
    for instr in sorted(replacements, key=len, reverse=True):
        pattern = re.compile(rf'(?<!\.)\b{re.escape(instr)}')
        eq = pattern.sub(replacements[instr], eq)

    # Split suivant ; si plusieurs formules sont definies
    eq = eq.split(';')

    for eq0 in eq:
        # Extrait la variable a initialiser de eq
        s = eq0.split('=', 1)
        #if len(s) != 2:
        #    print('Error: initVars: equation is incorrect.'); return None

        var = s[0]; var = var.replace('{', ''); var = var.replace('}', '')
        var = var.lstrip(); var = var.rstrip()
        varp = KCore.isNamePresent(a, var)
        if varp == -1:
            _addVars(a, var); varp = KCore.isNamePresent(a, var)

        # Initialisation de la variable
        if not isinstance(a[1], list): # array1
            loc = s[1]; c = 0
            for v in vars:
                loc = loc.replace('{%s}'%v, 'ap[%d,:]'%c); c += 1
            ap = a[1]
            ap[varp,:] = eval(loc)
        else: # array2/3
            loc = s[1]; ap = a[1]
            ap1 = [ ap[c].ravel(order='K') for c in range(len(ap)) ]
            c = 0
            for v in vars:
                loc = loc.replace('{%s}'%v, 'ap1[%d][:]'%c); c += 1
            ap1[varp][:] = eval(loc)
    return None

# Initialisation par une formule avec expression.ast
def _initVarByEq2__(a, eq):
    from . import expression as expr
    # Extrait les variables de a
    varstring = a[0]
    vars = varstring.split(',')

    eq = eq.replace('centers:', '')
    eq = eq.replace('nodes:', '')

    # Split suivant ; si plusieurs formules sont definies
    eq = eq.split(';')

    for eq0 in eq:
        ast_eq = expr.ast(eq0)
        # Extrait la variable a initialiser de eq
        s = eq0.split('=', 1)

        var = s[0]; var = var.replace('{', ''); var = var.replace('}', '')
        var = var.lstrip(); var = var.rstrip()
        varp = KCore.isNamePresent(a, var)
        if varp == -1:
            _addVars(a, var); varp = KCore.isNamePresent(a, var)

        ast_eq.run(a)
    return None

# Get index field
def getIndexField__(a):
    npts = getNPts(a)
    n = numpy.arange(0, npts, dtype=numpy.float64)
    # structure/non structure/array1/array2
    if len(a) == 5: # structure
        if isinstance(a[1], list): return ['index', [n], a[2], a[3], a[4]] # array2
        else: return ['index', n.reshape(1,npts), a[2], a[3], a[4]] # array1
    else:
        if isinstance(a[1], list): return ['index', [n], a[2], a[3]] # array2
        else: return ['index', n.reshape(1,npts), a[2], a[3]] # array1

def getIndexField(array):
    """Return the index field in an array.
    Usage: getIndexField(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(getIndexField__(i))
        return b
    else:
        return getIndexField__(array)

# Convert a file extension into a file format
def convertExt2Format__(fileName):
    """Convert a file extension into a file format."""
    ext2Format = {
        '.plt': 'bin_tp',
        '.dat': 'fmt_tp',
        '.tp': 'fmt_tp',
        '.v3d': 'bin_v3d',
        '.fv3d': 'fmt_v3d',
        '.vtk': 'bin_vtk',
        '.mesh': 'fmt_mesh',
        '.msh': 'fmt_gmsh',
        '.stl': 'fmt_stl',
        '.fstl': 'fmt_stl',
        '.bstl': 'bin_stl',
        '.selig': 'fmt_selig',
        '.gltf': 'bin_gltf',
        '.glb': 'bin_gltf',
        '.fig': 'fmt_xfig',
        '.svg': 'fmt_svg',
        '.pov': 'fmt_pov',
        '.cgns': 'bin_cgns',
        '.adf': 'bin_adf',
        '.hdf': 'bin_hdf',
        '.mod': 'bin_hdf',
        '.grid': 'bin_tau',
        '.h5': 'bin_fsdm',
        '.pickle': 'bin_pickle',
        '.df3': 'bin_df3',
        '.3ds': 'bin_3ds',
        '.ply': 'bin_ply',
        '.obj': 'fmt_obj',
        '.gts': 'fmt_gts',
        '.png': 'bin_png',
        '.jpg': 'bin_jpg',
        '.jpeg': 'bin_jpg',
        '.d': 'fmt_cedre',
        '.su2': 'fmt_su2',
        '.foam': 'fmt_foam',
        '.gbin': 'bin_plot3d',
        '.gfmt': 'fmt_plot3d',
        '.arc': 'bin_arc',
        '.iges': 'fmt_iges',
        '.igs': 'fmt_iges',
        '.stp': 'fmt_step',
        '.step': 'fmt_step',
        '.ref': 'bin_pickle',
        '.ref1': 'bin_pickle',
        '.ref2': 'bin_pickle'
    }
    extension = os.path.splitext(fileName)[-1]
    fmt = ext2Format.get(extension.lower(), 'unknown')
    return fmt

def convertFile2Arrays(fileName, format=None, nptsCurve=20, nptsLine=2,
                       density=-1., zoneNames=None, BCFaces=None, BCFields=None,
                       hmax=0.0, hausd=1., grow=0.0, mergeTol=-1, occAlgo=0,
                       centerArrays=None, api=1):
    """Read file and return arrays containing file data.
    Usage: a = convertFile2Arrays(fileName, options)"""
    try: import locale; locale.setlocale(locale.LC_NUMERIC, 'C') # force .
    except: pass
    if format is None: format = convertExt2Format__(fileName); autoTry = True
    else: autoTry = False
    exists = os.path.exists(fileName)
    if not exists: raise IOError("convertFile2Arrays: file %s not found."%fileName)

    if format == 'bin_pickle':
        import pickle
        print('Reading \''+fileName+'\'...', end="")
        try:
            file = open(fileName, 'rb')
            oldData = False
            if oldData: a = pickle.load(file, encoding='latin1')
            else: a = pickle.load(file)
            file.close()
        except:
            raise TypeError("convertFile2Arrays: file %s can not be read."%fileName)
        else:
            # convert i8/i4 types if necessary
            if isinstance(a[0], list):
                for i in a:
                    if len(i) == 4:
                        if isinstance(i[2], list): # array2/3
                            for c, k in enumerate(i[2]):
                                if k.dtype != E_NpyInt: i[2][c] = k.astype(E_NpyInt, order='K')
                        else: # array1
                            if i[2].dtype != E_NpyInt: i[2] = i[2].astype(E_NpyInt, order='K')
            else:
                if len(a) == 4:
                    if isinstance(a[2], list): # array2/3
                        for c, k in enumerate(a[2]):
                            if k.dtype != E_NpyInt: a[2][c] = k.astype(E_NpyInt, order='K')
                    else: # array1
                        if a[2].dtype != E_NpyInt: a[2] = a[2].astype(E_NpyInt, order='K')

            print('done.')
            return a
    elif format == 'fmt_iges' or format == 'fmt_step':
        try: import OCC
        except: raise ImportError("convertFile2Arrays: CAD readers requires OCC module.")
        a = OCC.convertCAD2Arrays(fileName, format=format, h=hmax,
                                  chordal_err=hausd, growth_ratio=grow,
                                  merge_tol=mergeTol, algo=occAlgo)
        if zoneNames is not None:
            for c in range(len(a)): zoneNames.append('zone%d'%c)
        return a
    elif format == 'fmt_free':
        print('Reading '+fileName+' (fmt_free)...', end="")
        try:
            file = open(fileName, 'r')
            f = file.read()
            file.close()
            f = f.split('\n')
            np = len(f) # nbre d'enregistrements
            if f[np-1] == '' or f[np-1] == ' ': np = np-1
            nv = 1 # nbre de variables
            if np > 2: line = f[2]
            elif np > 1: line = f[1]
            elif np > 0: line = f[0]
            else: line = ''
            line = line.split(' ')
            nv = len(line)
            # varstring
            varString = ''
            for i in range(nv-1): varString += 'v%d,'%(i+1)
            if nv > 0: varString += 'v%d'%nv
            a = array(varString, np, 1, 1)
            pt = a[1]
            for i in range(np):
                line = f[i]; line = line.split(' ')
                for n in range(len(line)): pt[n,i] = float(line[n])
            print('done.')
            zoneNames.append('zone0')
            return [a]
        except:
            raise TypeError("convertFile2Arrays: file %s can not be read."%fileName)
    else:
        try:
            return converter.convertFile2Arrays(fileName, format, nptsCurve, nptsLine, density, zoneNames, BCFaces, BCFields, centerArrays, api)
        except:
            if not autoTry: raise
            else: pass

            format = checkFileType(fileName)
            try:
                return converter.convertFile2Arrays(fileName, format, nptsCurve, nptsLine, density, zoneNames, BCFaces, centerArrays, api)
            except:
                FORMATS = [
                    'bin_ply', 'fmt_tp', 'fmt_v3d',
                    'bin_tp', 'bin_v3d', 'bin_vtk', 'fmt_mesh',
                    'fmt_gmsh', 'bin_gmsh', 'fmt_stl',
                    'bin_stl', 'bin_gltf',
                    'fmt_xfig', 'fmt_svg', 'bin_3ds',
                    'fmt_obj', 'fmt_gts' , 'fmt_pov', 'bin_arc'
                ]
                for fmt in FORMATS:
                    try:
                        a = converter.convertFile2Arrays(fileName, fmt, nptsCurve, nptsLine, density, zoneNames, BCFaces, BCFields, centerArrays, api)
                        return a
                    except:
                        return converter.convertFile2Arrays(fileName, format, nptsCurve, nptsLine, density, zoneNames, BCFaces, BCFields, centerArrays, api)

def convertArrays2File(arrays, fileName, format=None, isize=8, rsize=8,
                       endian='big', colormap=0, dataFormat='%.9e ',
                       zoneNames=[], BCFaces=[]):
    """Write arrays to output file.
    Usage: convertArrays2File(arrays, fileName, format, options)"""
    try: import locale; locale.setlocale(locale.LC_NUMERIC, 'C') # force .
    except: pass
    if arrays != [] and not isinstance(arrays[0], list): arrays = [arrays]
    znames = []
    if len(zoneNames) == 0:
        for iz in range(len(arrays)): znames.append('Zone%d'%(iz+1))
    else:
        znames = zoneNames
        if len(zoneNames) != len(arrays):
            raise ValueError("convertArrays2File: zoneNames list, if not empty, must have the same length of arrays list.")
    if format is None:
        format = convertExt2Format__(fileName)
    if format == 'bin_pickle':
        import pickle
        file = open(fileName, 'wb')
        print('Writing \''+fileName+'\'...', end="")
        pickle.dump(arrays, file, protocol=pickle.HIGHEST_PROTOCOL); file.close()
        print('done.')
    elif format == 'fmt_free':
        file = open(fileName, 'w')
        print('Writing %s (fmt_free)...'%fileName),
        for a in arrays:
            out = ''
            pt = a[1]
            nv = pt.shape[0]
            s = pt.shape[1]
            for i in range(s):
                for n in range(nv-1):
                    out += dataFormat%(pt[n,i])
                    out += ' '
                out += dataFormat%(pt[nv-1,i])
                out += '\n'
            file.write(out)
        file.close()
        print('done.')
    else:
        converter.convertArrays2File(arrays, fileName, format, isize, rsize,
                                     endian, colormap, dataFormat,
                                     znames, BCFaces)

def diffArrays(arrays1, arrays2, arrays3=[], atol=1.e11, rtol=0.):
    """Diff arrays defining solutions. Return the delta field."""
    if arrays3 != []:
        return converter.diffArrays(arrays1, arrays2, arrays3)
    else:
        return converter.diffArrays(arrays1, arrays2, atol, rtol)

def diffArraysGeom(array1, array2, atol=1.e-10, rtol=0.):
    """Diff arrays defining solutions, geometrically. Return the delta field."""
    if isinstance(array1[0], list):
        ret = []
        for c, a in enumerate(array1):
            res = diffArraysGeom__(a, array2[c], atol, rtol)
            if res is None: return None
            ret.append(res[0])
        return ret
    else:
        return diffArraysGeom__(array1, array2, atol, rtol)

def diffArraysGeom__(array1, array2, atol=1.e-10, rtol=0.):
    hook = createHook(array1, 'nodes')
    ids = identifyNodes(hook, array2, tol=atol)
    if numpy.any(ids == -1):
        print("ids", ids)
        freeHook(hook)
        return None # one array is different on coordinates
    ids[:] -= 1
    pt = array1[1]
    # renumerote a in place
    if isinstance(pt, list): # array2/3
        pt2 = []
        for c, p in enumerate(pt):
            p2 = numpy.copy(p, order='K')
            pr2 = p2.ravel(order='K')
            pr = pt[c].ravel(order='K')
            pr2[:] = pr[ids[:]]
            pt2.append(p2)
    else: # array1
        pt2 = numpy.copy(pt, order='K')
        pt2[:,:] = pt[:,ids[:]]

    array1[1] = pt2
    ret2 = diffArrays([array1], [array2], atol, rtol)
    freeHook(hook)
    return ret2

def isFinite__(a, var=None):
    nfld = a[1].shape[0]
    vars = getVarNames(a)
    ret = True
    for c, v in enumerate(vars):
        if var is None or v == var:
            ptr = a[1][c]
            ptr = ptr.ravel(order="K")
            #b = numpy.isfinite(ptr)
            #res = numpy.all(b)
            res = converter.isFinite(ptr)
            if res > 0:
                ret = False
                print('Warning: NAN or INF value in field (%s)'%v)
    return ret

def isFinite(array, var=None):
    """Return true if all fields have no NAN or INF values."""
    if isinstance(array[0], list):
        ret = True
        for a in array:
            ret1 = isFinite__(a, var)
            if not ret1: ret = False
        return ret
    else: return isFinite__(array, var)

def _setNANValuesAt__(a, var=None, value=0.):
    nfld = a[1].shape[0]
    vars = getVarNames(a)
    for c, v in enumerate(vars):
        if var is None or v == var:
            ptr = a[1][c]
            ptr = ptr.ravel(order="K")
            converter.setNANValuesAt(ptr, value)
    return None

def _setNANValuesAt(array, var=None, value=0.):
    """Set NAN values at value."""
    if isinstance(array[0], list):
        for a in array: _setNANValuesAt__(a, var)
    else: return _setNANValuesAt__(array, var)

def setNANValuesAt(array, var=None, value=0.):
    """Set NAN values at value."""
    b = copy(array)
    return _setNANValuesAt(b, var, value)

def getValue(array, ind):
    """Return the values of an array for a point of index ind or (i,j,k)...
    Usage: getValue(array, ind)"""
    if isinstance(array[0], list):
        raise TypeError("getValue: only for one array.")

    if isinstance(ind, tuple):
        if len(array) != 5: # structure
            raise TypeError("getValue: (i,j,k) indexing is only for structured array.")
        ni = array[2]; nj = array[3]
        if len(ind) == 3: index = (ind[0]-1)+(ind[1]-1)*ni+(ind[2]-1)*ni*nj
        elif len(ind) == 2: index = (ind[0]-1)+(ind[1]-1)*ni
        else: raise ValueError("getValue: too much values in index tuple.")
    else: index = ind

    n = array[1]
    if isinstance(n, list): # array3
        nfld = len(n)
        v = []
        for nf in range(nfld): v.append(n[nf].ravel('k')[index])
        return v
    else: # array1
        nfld = n.shape[0]
        v = []
        for nf in range(nfld): v.append(n[nf, index])
        return v

def setValue(array, ind, values):
    """Set the values in an array for a point of index ind or (i,j,k)...
    Usage: setValue(array, ind, values)"""
    if isinstance(array[0], list): raise TypeError("setValue: only for one array.")

    if isinstance(ind, tuple):
        if len(array) != 5: # structure
            raise TypeError("setValue: (i,j,k) indexing is only for structured array.")
        ni = array[2]; nj = array[3]
        if len(ind) == 3: index = (ind[0]-1)+(ind[1]-1)*ni+(ind[2]-1)*ni*nj
        elif len(ind) == 2: index = (ind[0]-1)+(ind[1]-1)*ni
        else: raise ValueError("setValue: too much values in index tuple.")
    else: index = ind
    ar = array[1]
    nf2 = len(values)
    if isinstance(ar, list): # array3
        nf = len(ar)
        if nf2 != nf: raise ValueError("setValue: values is badly dimensioned.")
        for v in range(nf): ar[v].ravel('k')[index] = values[v]
    else:
        nf = ar.shape[0]
        if nf2 != nf: raise ValueError("setValue: values is badly dimensioned.")
        ar[:, index] = values[:]
    return None

def getArgMin(array, varName):
    """Get array value where the variable defined by varName is minimum.
    Usage: getArgMin(array, varName)"""
    if isinstance(array[0], list):
        pos1 = KCore.isNamePresent(array[0], varName)
        b = converter.getArgMin(array[0], varName)
        for i in array:
            pos2 = KCore.isNamePresent(i, varName)
            a = converter.getArgMin(i, varName)
            if a[pos2] < b[pos1]:
                b = a
                pos1 = pos2
        return b
    else:
        return converter.getArgMin(array, varName)

def getMinValue(array, varName):
    """Get the minimum value of variable defined by varName in array.
    Usage: getMinValue(array, varName)"""
    if not isinstance(varName, list): varNames = [varName]
    else: varNames = varName
    out = []
    for v in varNames:
        if isinstance(array[0], list): pos = KCore.isNamePresent(array[0], v)
        else: pos = KCore.isNamePresent(array, v)
        out.append(getArgMin(array, varName)[pos])
    if len(out) == 1: return out[0]
    else: return out

def getArgMax(array, varName):
    """Get array value where the variable defined by varName is maximum.
    Usage: getArgMax(array, varName)"""
    if isinstance(array[0], list):
        pos1 = KCore.isNamePresent(array[0], varName)
        b = converter.getArgMax(array[0], varName)
        for i in array:
            pos2 = KCore.isNamePresent(i, varName)
            a = converter.getArgMax(i, varName)
            if a[pos2] > b[pos1]:
                b = a
                pos1 = pos2
        return b
    else:
        return converter.getArgMax(array, varName)

def getMaxValue(array, varName):
    """Get the maximum value of variable defined by varName in array.
    Usage: getMaxValue(array, varName)"""
    if not isinstance(varName, list): varNames = [varName]
    else: varNames = varName
    out = []
    for v in varNames:
        if isinstance(array[0], list): pos = KCore.isNamePresent(array[0], v)
        else: pos = KCore.isNamePresent(array, v)
        out.append(getArgMax(array, varName)[pos])
    if len(out) == 1: return out[0]
    else: return out

def getMeanValue(array, varName):
    """Get the mean value of the variable defined by varName in an array.
    Usage: getMeanValue(array, varName)"""
    if isinstance(array[0], list):
        pos1 = KCore.isNamePresent(array[0], varName)
        b = 0.
        lb = 0.
        for i in array:
            pos2 = KCore.isNamePresent(i, varName)
            a = converter.getMeanValue(i, varName)
            if isinstance(i[1], list): # array2/3
                la = i[1][0].size * 1.
            else: la = i[1].shape[1] * 1.
            b = lb*b + la*a
            lb = lb + la
            b = b / lb
        return b
    else:
        return converter.getMeanValue(array, varName)

def getMeanRangeValue(array, varName, rmin, rmax):
    """Get the mean value of the variable defined by varName for a sorted range in an array.
    Usage: getMeanRangeValue(array, varName, rmin, rmax)"""
    if isinstance(array[0], list):
        pos1 = KCore.isNamePresent(array[0], varName)
        b = 0.
        lb = 0.
        for i in array:
            pos2 = KCore.isNamePresent(i, varName)
            a = converter.getMeanRangeValue(i, varName, rmin, rmax)
            if isinstance(i[1], list): # array2/3
                la = i[1][0].size*1.*(rmax-rmin)
            else: la = i[1].shape[1]*1.*(rmax-rmin)
            b = lb*b + la*a
            lb = lb + la
            b = b / lb
        return b
    else:
        return converter.getMeanRangeValue(array, varName, rmin, rmax)

def normL0(array, varName):
    """Get the L0 norm of the field defined by varName in the array.
    If celln exists in the array, the norm for blanked points is not computed.
    Usage: normL0(array, varName)"""
    if isinstance(array[0], list):
        norm = 0
        for i in array: norm = max(converter.normL0(i, varName), norm)
        return norm
    else: return converter.normL0(array, varName)

def normL2(array, varName):
    """Get the L2 norm of the field defined by varName in the array.
    If celln exists in the array, the norm for blanked points is not computed.
    Usage: normL2(array, varName)"""
    if isinstance(array[0], list):
        norm = 0.; nsize = 0
        for i in array:
            if isinstance(i[1], list): # array2
                nsize0 = i[1][0].size
            else: nsize0 = i[1].shape[1]
            nsize += nsize0
            norm0 = converter.normL2(i, varName)
            norm += norm0*norm0*nsize0
        if nsize > 0: norm = numpy.sqrt(norm/nsize)
        return norm
    else:
        return converter.normL2(array, varName)

def normalize(a, vars):
    """Get the normalisation of the fields defined by vars in the array.
    Usage: normalize(a, vars)"""
    b = copy(a)
    _normalize(b, vars)
    return b

def _normalize(a, vars):
    if isinstance(a[0], list):
        for i in a: converter.normalize(i, vars)
    else: converter.normalize(a, vars)
    return None

def magnitude(array, vars):
    """Get the magnitude of the fields defined by vars in the array.
    Usage: magnitude(array, vars)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(converter.magnitude(i, vars))
        return b
    else:
        return converter.magnitude(array, vars)

def convertHexa2Struct(array):
    """Convert an hexa array in a struct array.
    Usage: convertHexa2Struct(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(converter.convertHexa2Struct(i))
        return b
    else:
        return converter.convertHexa2Struct(array)

# Convert a BAR array without ramifications, closed into an i-array
def convertBAR2Struct(array):
    """Convert a BAR array without ramifications, closed into an i-array.
    Usage: convertBAR2Struct(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            if len(i) == 4: b.append(converter.convertBAR2Struct(i))
            else: b.append(i)
        return b
    else:
        if len(array) == 4: return converter.convertBAR2Struct(array)
        else: return array

def convertTri2Quad(array, alpha=30.):
    """Convert a TRI array to a QUAD array.
    Usage: convertTri2Quad(array, alpha)"""
    if isinstance(array[0], list):
        b = []; c = []
        for i in array:
            res = converter.convertTri2Quad(i, alpha)
            b += res[0]; c += res[1]
        return b, c
    else:
        return converter.convertTri2Quad(array, alpha)

def conformizeNGon(array, tol=1.e-6):
    """Conformize topologically a NGON array.
    Usage: conformizeNGon(array, tol)"""
    if len(array) == 4 and array[3] == 'NGON':
        return converter.conformizeNGon(array, tol)
    else: return array

def adaptSurfaceNGon(array):
    """Convert a surface NGon from one type (A: NGON=bars, NFACE=polygon)
    to another (B: NGON=polygon, NFACE=NULL).
    """
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(converter.adaptSurfaceNGon(i))
        return b
    else:
        return converter.adaptSurfaceNGon(array)

# -- interne --
def convertArray2Tetra1__(array, arrayC=[], split='simple'):
    try: sub = array[3]
    except: raise TypeError("convertArray2Tetra: arg must be an array.")

    if isinstance(sub, str): t = sub
    else: t = 'STRUCT'
    if split == 'simple': # no points added
        if t == 'STRUCT':
            return converter.convertStruct2Tetra(array)
        elif t == 'HEXA' or t == 'QUAD':
            return converter.convertHexa2Tetra(array)
        elif t == 'PENTA':
            return converter.convertPenta2Tetra(array)
        elif t == 'PYRA':
            return converter.convertPyra2Tetra(array)
        elif t == 'TETRA' or t == 'TRI' or t == 'NODE' or t == 'BAR':
            return array
        elif t == 'NGON':
            try: import Transform as T
            except: raise ImportError("convertArray2Tetra: requires Transform for NGONs.")
            tmp = T.breakElements(array)
            brd = []
            for i in tmp:
                if i[3] != 'NGON': brd.append(convertArray2Tetra1__(i))
                else:
                    brd.append(convertArray2Tetra1__(i, split="withBarycenters"))
            brd = T.join(brd)
            return brd
        elif t == 'HEXA*' or t == 'QUAD*':
            tmp = center2Node(array)
            tmp = converter.convertHexa2Tetra(tmp)
            return node2Center(tmp)
        elif t == 'PENTA*':
            tmp = center2Node(array)
            tmp = converter.convertPenta2Tetra(tmp)
            return node2Center(tmp)
        elif t == 'PYRA*':
            tmp = center2Node(array)
            tmp = converter.convertPyra2Tetra(tmp)
            return node2Center(tmp)
        elif t == 'TETRA*' or t == 'TRI*' or t == 'NODE*' or t == 'BAR*':
            return array
        else:
            raise TypeError("convertArray2Tetra: type of element is not taken into account: %s."%t)

    elif split == 'withBarycenters': # new points added at centers of elements and faces
        if t == 'STRUCT':
            if arrayC == []:
                return converter.convertStruct2TetraBary(array)
            else:
                return converter.convertStruct2TetraBaryBoth(array, arrayC)

        elif t == 'HEXA' or t == 'QUAD' or t == 'PENTA' or t == 'PYRA' or t == 'BAR':
            if arrayC == []:
                return converter.convertArray2TetraBary(array)
            else:
                return converter.convertArray2TetraBaryBoth(array, arrayC)

        elif t == 'TETRA' or t == 'TRI' or t == 'NODE':
            if arrayC == []: return array
            else: return [array,arrayC]

        elif t == 'NGON':
            if arrayC == []: return converter.convertNGon2TetraBary(array)
            else:
                return converter.convertNGon2TetraBaryBoth(array, arrayC)
        elif t == 'HEXA*' or t == 'PENTA*' or t == 'PYRA*' or t == 'BAR*':
            tmp = center2Node(array)
            tmp =  converter.convertArray2TetraBary(tmp)
            return node2Center(tmp)
        elif t == 'TETRA*' or t == 'TRI*' or t == 'NODE*':
            return array
        else:
            raise TypeError("convertArray2Tetra: with split='withBarycenters', type of element is not taken into account: %s."%t)
    else:
        raise TypeError("convertArray2Tetra: type of split unknown: %s."%split)

# -- Convert array(s) to tetra
def convertArray2Tetra(array, split='simple'):
    """Convert a array in an unstructured tetra array.
    Unstructured array is triangular in 2D and tetrahedral in 3D.
    Usage: convertArray2Tetra(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(convertArray2Tetra1__(i, split=split))
        return b
    else: return convertArray2Tetra1__(array, split=split)

# -- interne --
def convertArray2Hexa1__(array):
    try: sub = array[3]
    except: raise TypeError("convertArray2Hexa: arg must be an array.")
    if isinstance(sub, str): t = sub
    else: t = 'STRUCT'
    if t == 'STRUCT': return converter.convertStruct2Hexa(array)
    elif t == 'HEXA' or t == 'QUAD' or t == 'BAR' or t == 'NODE':
        return array
    elif t == 'HEXA*' or t == 'QUAD*' or t == 'BAR*' or t == 'NODE*':
        return array
    elif t == 'TRI' or t == 'TETRA' or t == 'PENTA':
        return converter.convertUnstruct2Hexa(array)
    elif t == 'NGON':
        try: import Transform as T
        except: raise ImportError("convertArray2Hexa: requires Transform for NGONs.")
        tmp = T.breakElements(array)
        brd = []
        for i in tmp:
            if i[3] != 'NGON': brd.append(convertArray2Hexa1__(i))
        brd = T.join(brd)
        return brd
    elif t == 'TRI*' or t == 'TETRA*' or t == 'PENTA*':
        return converter.convertUnstruct2Hexa(array)
    else: raise TypeError("convertArray2Hexa: type of element is not taken into account: %s."%t)

# -- convert arrays(s) to hexa
def convertArray2Hexa(array):
    """Convert a array in an unstructured hexa array.
    Unstructured array can be quad in 2D and hexa in 3D.
    Usage: convertArray2Hexa(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(convertArray2Hexa1__(i))
        return b
    else: return convertArray2Hexa1__(array)

def convertArray2NGon__(array, api=1):
    try: sub = array[3]
    except: raise TypeError("convertArray2NGon: arg must be an array.")
    if isinstance(sub, str): t = sub
    else: t = 'STRUCT'
    if t == 'STRUCT': return converter.convertStruct2NGon(array, api)
    elif t == 'NGON': return array
    else: return converter.convertUnstruct2NGon(array, api)

def convertArray2NGon(array, api=1):
    """Convert a array in a NGON array.
    Usage: convertArray2NGon(array, api)"""
    if isinstance(array[0], list):
        b = []
        for i in array: b.append(convertArray2NGon__(i, api))
        return b
    else: return convertArray2NGon__(array, api)

def convertPenta2Strand(array):
    """Convert a PENTA array to a STRAND array."""
    if isinstance(array[0], list):
        b = []
        for i in array: b.append(converter.convertPenta2Strand(i))
        return b
    else: return converter.convertPenta2Strand(array)

def convertStrand2Penta(array):
    """Convert a STRAND array to a PENTA array."""
    if isinstance(array[0], list):
        b = []
        for i in array: b.append(converter.convertStrand2Penta(i))
        return b
    else: return converter.convertStrand2Penta(array)

def node2Center(array, accurate=0):
    """Convert array defined on nodes to array defined on centers.
    Usage: node2Center(array, accurate=0)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(converter.node2Center(i, accurate))
        return b
    else:
        return converter.node2Center(array, accurate)

def center2Node(array, cellNType=0, BCFields=None):
    """Convert array defined on centers to array defined on nodes.
    Usage: center2Node(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(converter.center2Node(i, cellNType, BCFields))
        return b
    else:
        b = converter.center2Node(array, cellNType, BCFields)
        return b

def node2ExtCenter(array):
    """Convert array defined on nodes to an array defined on extended centers.
    Usage: node2ExtCenter(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(converter.node2ExtCenter(i))
        return b
    else:
        return converter.node2ExtCenter(array)

def extCenter2Node(array):
    """Convert array defined on extended centers to an array defined on nodes.
    Usage: extCenter2Node(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(converter.extCenter2Node(i))
        return b
    else:
        return converter.extCenter2Node(array)

def center2ExtCenter(array):
    """Convert array defined for centers to an array defined for extended centers.
    Usage: center2ExtCenter(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(converter.center2ExtCenter(i))
        return b
    else:
        return converter.center2ExtCenter(array)

# -- convert arrays(s) to node array(s).
def convertArray2Node(array):
    """Convert an array in an unstructured node array.
    Usage: convertArray2Node(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(converter.convertArray2Node(i))
        return b
    else:
        return converter.convertArray2Node(array)

#==============================================================================
def addGhostCellsNGon(arrayN, arrayC=[],depth=2):
    """Add depth layers of ghost cells in an NGON array.
    Usage: addGhostCellsNGon(arrayN,arrayC,depth)"""
    if isinstance(arrayN[0], list):
        b = []
        nzones = len(arrayN); nzonesC = len(arrayC)
        if nzones == nzonesC:
            for noz in range(nzones):
                if arrayN[noz] != [] and arrayC[noz] == []:
                    res = converter.addGhostCellsNGonNodes(arrayN[noz],depth)
                    b.append(res)
                elif arrayN[noz] == [] and arrayC[noz] != []:
                    res = converter.addGhostCellsNGonCenters(arrayC[noz],depth)
                    b.append(res)
                else:
                    res = converter.addGhostCellsNGonBoth(arrayN[noz], arrayC[noz],depth)
                    b.append(res)
        else:
            for noz in range(nzones):
                b.append(converter.addGhostCellsNGon(arrayN[noz],depth))
        return b
    else:
        if arrayN != [] and arrayC==[]: return converter.addGhostCellsNGonNodes(arrayN,depth)
        elif arrayN == [] and arrayC != []: return converter.addGhostCellsNGonCenters(arrayC,depth)
        else: return converter.addGhostCellsNGonBoth(arrayN, arrayC, depth)

def rmGhostCellsNGon(arrayN, arrayC=[], depth=2):
    """Delete depth layers of cells at exterior borders of a NGON mesh.
    Usage: rmGhostCellsNGon(arrayN, arrayC, depth)"""
    if isinstance(arrayN[0], list):
        b = []
        nzones = len(arrayN); nzonesC = len(arrayC)
        if nzones == nzonesC:
            for noz in range(nzones):
                if arrayN[noz] != [] and arrayC[noz] == []:
                    res = converter.rmGhostCellsNGonNodes(arrayN[noz],depth)
                    b.append(res)
                elif arrayN[noz] == [] and arrayC[noz] != []:
                    res = converter.rmGhostCellsNGonCenters(arrayC[noz],depth)
                    b.append(res)
                else:
                    res = converter.rmGhostCellsNGonBoth(arrayN[noz], arrayC[noz],depth)
                    b.append(res)
        else:
            for noz in range(nzones):
                b.append(converter.rmGhostCellsNGon(arrayN[noz],depth))
        return b
    else:
        if arrayN != [] and arrayC==[]: return converter.rmGhostCellsNGonNodes(arrayN,depth)
        elif arrayN == [] and arrayC != []: return converter.rmGhostCellsNGonCenters(arrayC,depth)
        else: return converter.rmGhostCellsNGonBoth(arrayN, arrayC, depth)

#==============================================================================
# Retourne une liste des noms de variables presents dans l'arrays
#==============================================================================
def getVarNames(a):
    """Get variable names.
    Usage: getVarNames(a)"""
    if isinstance(a[0], list):
        allvars = []
        for i in a:
            v = i[0].split(",")
            allvars.append(v)
    else: allvars = a[0].split(",")
    return allvars

#==============================================================================
# Fonctions de preconditionement (hook)
# IN: function: le nom de la fonction qui va utiliser le hook
#==============================================================================
def createHook(a, function='None'):
    """Create a hook for a given function.
    Usage: hook = createHook(a, function)"""
    if function == 'None': return None
    elif function == 'faceCenters': # 0
        # Retourne un KDT pour les centres des faces
        if isinstance(a[0], list):
            b = []
            for i in a:
                b.append(converter.registerFaces(i))
            return b
        else:
            return converter.registerFaces(a)
    elif function == 'nodes': # 2
        # Retourne un KDT pour les noeuds
        if isinstance(a[0], list):
            b = []
            for i in a:
                b.append(converter.registerNodes(i))
            return b
        else:
            return converter.registerNodes(a)
    elif function == 'elementCenters': # 3
        # Retourne un KDT pour les centres des elements
        if isinstance(a[0], list):
            b = []
            for i in a:
                b.append(converter.registerElements(i))
                #b.append(converter.registerElements(convertArray2NGon(i)))
            return b
        else:
            return converter.registerElements(a)
            #return converter.registerElements(convertArray2NGon(a))
    elif function == 'extractMesh': # 1
        # Retourne un ADT pour les elements
        try: import Post as P
        except: inl = a
        else: inl, modified = P.growOfEps__(a, 1.e-6, nlayers=2, planarity=False)
        return converter.registerCells(inl, None, None, 0, 0.)

    elif function == 'adt': # 1 ADT pour les interpolations
        return converter.registerCells(a, None, None, 0, 0.)

    else: raise ValueError("function %s is invalid."%function)

def createHookAdtCyl(a, center=(0,0,0), axis=(0,0,1), depth=0, thetaShift=0.):
    """Create a hook for cylindrical adt."""
    return converter.registerCells(a, center, axis, depth, thetaShift)

#===============================================================================
# Fonctions de preconditionement (hook)
# IN: function: le nom de la fonction qui va utiliser le hook
# IN: a: liste des zones utilisees pour faire le preconditionnement global
# IN: indir: 0: pas d'indirection sur les zones
#            1: sauvegarde de indirZones:
#               indirection sur le no de la zone pour chq pt du hook
# OUT: hook(,indirZones)
#===============================================================================
def createGlobalHook(a, function='None', indir=0):
    """Create a hook for a set of zones and for a given function.
    Usage: hook = createGlobalHook(a, function)"""
    if function == 'None': return None
    elif function == 'faceCenters': # 0
        # Retourne un KDT pour les centres des faces
        if not isinstance(a[0],list): return converter.registerAllFaces([a], indir)
        else: return converter.registerAllFaces(a, indir)
    elif function == 'nodes': # 2
        # Retourne un KDT pour les noeuds
        if not isinstance(a[0],list): return converter.registerAllNodes([a], indir)
        else: return converter.registerAllNodes(a, indir)
    elif function == 'elementCenters': # 3
        # Retourne un KDT pour les centres des elements
        if not isinstance(a[0],list): return converter.registerAllElements(convertArray2NGon([a]), indir)
        else: return converter.registerAllElements(convertArray2NGon(a), indir)
    elif function == 'extractMesh': # 1
        raise ValueError('function=extractMesh not implemented for global hook.')
    else: raise ValueError("function is invalid.")

#==============================================================================
def freeHook(hook):
    """Free hook.
    Usage: freeHook(hook)"""
    if isinstance(hook, list):
        for i in hook:
            if i is not None: converter.freeHook(i)
    else:
        if hook is not None: converter.freeHook(hook)

#==============================================================================
# Fonctions d'identification geometrique
#==============================================================================
def identifyNodes(hook, a, tol=1.e-11):
    """Identify nodes of a in KDT. return identified indices.
    Usage: identifyNodes(hook, a)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(converter.identifyNodes(hook, i, tol))
        return b
    else:
        return converter.identifyNodes(hook, a, tol)

def identifyFaces(hook, a, tol=1.e-11):
    """Identify face centers of a in KDT. return identified indices.
    Usage: identifyFaces(hook, a)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(converter.identifyFaces(hook, convertArray2NGon(i), tol))
        return b
    else:
        return converter.identifyFaces(hook, convertArray2NGon(a), tol)

def identifyElements(hook, a, tol=1.e-11):
    """Identify element centers of a in KDT. return identified indices.
    Usage: identifyElement(hook, a)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            # if len(i) == 5: i = convertArray2Hexa(i)
            b.append(converter.identifyElements(hook, i, tol))
        return b
    else:
        # if len(a) == 5: a = convertArray2Hexa(a)
        return converter.identifyElements(hook, a, tol)

#=============================================================================
def identifySolutions(coordsRcv, solDnr, hookDnr, vars=[], tol=1.e6):
    """Identify points in a hook to mesh points and set the solution if donor
    and receptor points are distant from tol.
    Usage: identifySolutions(coordsRcv, solDnr, hookDnr, vars, tol)"""
    if vars != []: solDnr = extractVars(solDnr, vars)
    res = None
    if isinstance(coordsRcv[0], list): # receptor is a list of zones
        if isinstance(solDnr[0], list):
            res = converter.identifySolutions(hookDnr, solDnr, coordsRcv, tol)
        elif not isinstance(solDnr[0], list): # une seule zone
            res = converter.identifySolutions(hookDnr, [solDnr], coordsRcv, tol)

    else: # receptor is a single zone
        if isinstance(solDnr[0], list):
            res = converter.identifySolutions(hookDnr, solDnr, [coordsRcv], tol)[0]
        elif not isinstance(solDnr[0], list): # une seule zone
            res = converter.identifySolutions(hookDnr, [solDnr], [coordsRcv], tol)[0]
    return res

#==============================================================================
# Fonctions d'identification du noeud/element/face le plus proche
#==============================================================================
def nearestNodes(hook, a):
    """Find in KDT nearest points to nodes of a. return identified indices.
    Usage: nearestNodes(hook, a)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(converter.nearestNodes(hook, i))
        return b
    else:
        return converter.nearestNodes(hook, a)

def nearestFaces(hook, a):
    """Find in KDT nearest points to face centers of a. return identified indices.
    Usage: nearestFaces(hook, a)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(converter.nearestFaces(hook, convertArray2NGon(i)))
        return b
    else:
        return converter.nearestFaces(hook, convertArray2NGon(a))

def nearestElements(hook, a):
    """Find in KDT nearest points to element centers of a. return identified indices.
    Usage: nearestElement(hook, a)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            if len(i) == 5: i = convertArray2Hexa(i)
            b.append(converter.nearestElements(hook, i))
        return b
    else:
        if len(a)==5: a = convertArray2Hexa(a)
        return converter.nearestElements(hook, a)

# Create global index
def createGlobalIndex(a, start=0):
    """Create the global index field."""
    b = copy(a)
    _createGlobalIndex(b, start)
    return b

def _createGlobalIndex(a, start=0):
    """Create the global index field."""
    _initVars(a, 'globalIndex', 0)
    if isinstance(a[0], list):
        for i in a: converter.createGlobalIndex(i, start)
        return a
    else:
        converter.createGlobalIndex(a, start)
        return a

def recoverGlobalIndex(a, b):
    """Recover fields of b in a following the global index field."""
    c = copy(a)
    _recoverGlobalIndex(c, b)
    return c

def _recoverGlobalIndex(a, b):
    """Recover fields of b in a following the global index field."""
    if isinstance(b[0], list):
        for i in b:
            variables = getVarNames(i); _addVars(a, variables)
    else:
        variables = getVarNames(b); _addVars(a, variables)
    if isinstance(b[0], list):
        if isinstance(a[0], list):
            for bi in b:
                for ai in a:
                    converter.recoverGlobalIndex(bi, ai)
        else:
            for bi in b:
                converter.recoverGlobalIndex(bi, a)
    else:
        if isinstance(a[0], list):
            for ai in a:
                converter.recoverGlobalIndex(b, ai)
        else:
            converter.recoverGlobalIndex(b, a)
    return None

# mergeConnectivity: merge deux connectivites en elements basiques
# en un seul maillage multiconnectivite (array3 uniquement)
# concatenation simple
def mergeConnectivity(a1, a2):
    if a1[0] != a2[0]: raise ValueError('mergeConnectivity: only for same fields.')
    if len(a1) != 4 or len(a2) != 4: raise ValueError('mergeConnectivity: only for unstructured arrays.')
    if a1[3] == 'NGON' or a2[3] == 'NGON': raise ValueError('mergeConnectivity: only for element arrays.')

    # fields
    f1 = a1[1]; f2 = a2[1]
    fo = []
    for c, f in enumerate(f1): fo.append(numpy.concatenate((f, f2[c])))
    npts = len(f1[0])

    # connectivity
    c1 = a1[2]; c2 = a2[2]
    for cp in c2: cp[:] = cp[:]+npts
    co = c1+c2
    return [a1[0], fo, co, a1[3]+','+a2[3]]


# Retourne -1: la variable n'est presente dans aucun array
# Retourne 0: la variable est presente dans au moins un array
# Retourne 1: la variable est presente dans tous les arrays
def isNamePresent(a, varname):
    """Test if varName is present in a."""
    if isinstance(a[0], list):
        one = 0
        for i in a:
            p = KCore.isNamePresent(i, varname)
            if p != -1: one += 1
        if one == len(a): return 1
        elif one == 0: return -1
        else: return 0
    else: # un seul array
        p = KCore.isNamePresent(a, varname)
        if p == -1: return -1
        else: return 1

def _signNGonFaces(a):
    """Return a consistently oriented pyTree with signed faces.
    Usage: _signNGonFaces(a)"""
    if isinstance(a[0], list):
        for i in a: converter.signNGonFaces(i)
    else: converter.signNGonFaces(a)
    return None

def _unsignNGonFaces(a):
    """Unsign NFACE connectivity"""
    if isinstance(a[0], list):
        isSigned = 1
        for i in a:
            if isSigned == 1: isSigned = converter.unsignNGonFaces(i)
            else: break # stop if unsigned
    else:
        isSigned = converter.unsignNGonFaces(a)
    return isSigned

def makeParentElements(a):
    PEs = []
    if isinstance(a[0], list):
        for i in a:
            PEs.append(converter.makeParentElements(i))
    else:
        PEs.append(converter.makeParentElements(a))
    return PEs


# convert to low order mesh
def convertHO2LO(a, mode=0):
    """Convert a HO mesh to a low order mesh.
    Usage: convertHO2LO(a, mode)"""
    if isinstance(a[0], list):
        out = []
        for i in a:
            out.append(converter.convertHO2LO(i, mode))
        return out
    else:
        b = converter.convertHO2LO(a, mode)
        return b

# convert to high order mesh
def convertLO2HO(a, mode=0, order=2):
    """Convert a LO mesh to a high order mesh.
    Usage: convertLO2HO(a, mode, order)"""
    if isinstance(a[0], list):
        out = []
        for i in a:
            out.append(converter.convertLO2HO(i, mode, order))
        return out
    else:
        b = converter.convertLO2HO(a, mode, order)
        return b

#==============================================================================
# Client/Server - send
#==============================================================================
def send(data, host='localhost', rank=0, port=15555):
    """Send data to socket."""
    import socket; import Compressor
    port = port+rank; sizeBuf = 1024

    #print('connecting to port', port)
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect((host, port))
    except:
        print('Send: can not connect to %s [%d]. Nothing sent.'%(host,port))
        return

    data = Compressor.pack(data, method=0) # serialize
    size = len(data)

    # Blocks
    header = (size,sizeBuf)
    #print('sending', header)
    header = Compressor.pack(header, method=0)
    header = header.ljust(255)

    s.send(header)
    nbytes = 0
    while nbytes < size:
        if nbytes+sizeBuf > size:
            nbytes += s.send(data[nbytes:nbytes+sizeBuf])
        else:
            nbytes += s.send(data[nbytes:])
        #print('send', nbytes, size)
    s.close()

#==============================================================================
# Client/server - createSockets
#==============================================================================
def createSockets(nprocs=1, port=15555):
    """Create sockets for communication."""
    import socket
    sockets = []
    for i in range(nprocs):
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.bind(('', port+i))
        s.listen(5)
        sockets.append(s)
    return sockets

#==============================================================================
# Client/server - listen
#==============================================================================
def listen(s):
    """Listen for sends."""
    import socket
    import Compressor
    while True:
        #s.listen(5)
        s.setblocking(1)
        client, address = s.accept()
        data = None
        nb = client.recv(255)
        if nb != b"":
            nb = nb.rstrip()
            (size,sizeBuf) = Compressor.unpack(nb, method=0)
            data = b''
            nbytes = 0
            while nbytes < size:
                if nbytes+sizeBuf < size:
                    received = client.recv(sizeBuf)
                else:
                    received = client.recv(size-nbytes)
                data += received; nbytes += len(received)
            data = Compressor.unpack(data, method=0)
            client.close()
            return data

#=============================================================================
# checkfile type
#=============================================================================
def checkFileType(fileName):
    """Find file type."""
    import os
    try: file = os.open(fileName, os.O_RDONLY)
    except: raise IOError("checkFileType: file %s not found."%fileName)
    #header = file.read(512)  # lecture des 512 premiers octets
    header = os.read(file, 512)
    os.close(file)

    if header[1:4] == b'HDF': return 'bin_hdf'
    if header[4:7] == b'ADF': return 'bin_adf'
    if header[0:5] == b'#!TDV': return 'bin_tp'
    if header[0:5] == b'TITLE' or header[0:5] == b'title' or header[0:9] == b"VARIABLES" or header[0:9] == b"variables" or header[0:8] == b"FILETYPE" or header[0:8] == b"filetype":
        return 'fmt_tp'
    if header.find(b"MeshVersionUnformatted") != -1: return 'bin_mesh'
    if header.find(b"MeshVersionFormatted") != -1: return 'fmt_mesh'
    if header.find(b"NDIME=") != -1: return 'fmt_su2'
    if header.find(b"DONNEES GENERALES") != -1 or header.find(b"--------") != -1: return 'fmt_cedre'
    if header[0:8] == b'CEDRE_IO': return 'bin_cedre'
    if header.find(b"solid") == 0: return 'fmt_stl'
    if header.find(b"v ") != -1: return 'fmt_obj'
    if header.find(b"$MeshFormat") != -1:
        EndMesh = header.find(b"$EndMeshFormat")
        if EndMesh == 20: return 'fmt_gmsh'
        elif EndMesh == 25: return 'bin_gmsh'
    if header.find(b"ply") != -1: return 'bin_ply'
    if header.find(b"#FIG") == 0: return 'fmt_xfig'
    if header.find(b"<svg") != -1: return 'fmt_svg'
    if header.find(b"PNG") == 1: return 'bin_png'
    import binascii as b
    beader = b.b2a_hex(header)
    eol = b"0a"
    if (beader[0:8] == b"04000000" or beader[0:8] == b"00000004") and (header[16:18] == b'va' or header[16:17] == b'x' or header[16:17] == b'y' or header[16:17] == b'z' or header[16:18] == b'VA' or header[16:17] == b'X' or header[16:17] == b'Y' or header[16:17] == b'Z'):
        return 'bin_v3d'
    if (beader[0:8] == b"08000000" or beader[0:8] == b"00000008") and (header[20:22] == b'va' or header[20:21] == b'x' or header[20:21] == b'y' or header[20:21] == b'z' or header[20:22] == b'VA' or header[20:21] == b'X' or header[20:21] == b'Y' or header[20:21] == b'Z'):
        return 'bin_v3d'
    if (beader[10:12] == eol and (beader[50:52] == b"78" or beader[50:52] == b"79" or beader[50:52] == b"80" or  beader[12:14] == b"78" or   beader[12:14] == b"79" or beader[12:14] == b"80"  or beader[50:52] == b"58" or beader[50:52] == b"59" or beader[50:52] == b"60" or  beader[12:14] == b"58" or beader[12:14] == b"59" or beader[12:14] == b"60" or beader[46:52] == b"766172" or beader[46:52] == b"564152" or beader[12:18] == b"766172" or beader[12:18] == b"564152")):
        return 'fmt_v3d'
    if beader.find(b"4d4d") == 0: return 'bin_3ds'
    if beader[0:4] == b"d8ff": return 'bin_jpg'
    dt = numpy.dtype('<i4')
    ieader = numpy.fromfile(fileName, dtype=dt, count=128, sep="")
    if ieader[0] == 4:
        ninjnk = ieader[1] * 3 * 4   # 3 pour ni,nj,nk, 4 pour 4 octets cas 3D multibloc
        ninj   = ieader[1] * 2 * 4 # 3 pour ni,nj, 4 pour 4 octets cas 2D multibloc
        if ieader[3] == ninjnk or ieader[3] == ninj:
            return 'bin_plot3d'
    if ieader[0] == 12 or ieader[0] == 8:  # cas 2D ou 3d monobloc
        return 'bin_plot3d'

    fileSize = os.path.getsize(fileName)
    try: ntri = header[80:82]; ntri = int(ntri)
    except: ntri = 0
    sizet = ntri*50+84  #format bin_stl 80 octets d'entete/nombre de triangles/50 octets par triangles
    if fileSize == sizet: return 'bin_stl'

    eolx1 = beader.find(eol, 0)
    eolx2 = beader.find(eol, eolx1 + 1)
    eol1 = (eolx1 +1)//2
    if eol1 != 0:
        file = open(fileName, 'r')
        i = 0
        ligne0=[]; ligne1=[]
        ninjnk_size = 0; npi = 0
        for line in file:
            if i == 0:
                ligne0 = line.split()
                npi = int(ligne0[0])
            else:
                newline = line.split()
                sfloat = newline[0].count(".")
                if sfloat == 0: ligne1.extend(newline)
                else: break
            i += 1
        file.close()
        ninjnk_size += len(ligne1)
        if ninjnk_size == 2 * npi or ninjnk_size == 3 * npi: return 'fmt_plot3d'
    return 'unknown'
