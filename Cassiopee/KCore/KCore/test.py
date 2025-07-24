"""Testing module for Cassiopee validation."""
# Systeme de validation des modules Cassiopee
# arrays: testA, [outA], stdTestA, [outTestA]
# pyTrees: testT, [outT], stdTestT, [outTestT]
# objects: testO
from __future__ import print_function
import numpy, sys, os

# global tolerance on float fields
TOLERANCE = 1.e-11

# whether to diffArrays geometrically or topologically. default is topologically.
GEOMETRIC_DIFF = False

# Data directory to store references
from KCore.Dist import getDataFolderName
DATA = getDataFolderName()

#=============================================================================
# Retourne la variable VALIDLOCAL si elle existe dans l'environnement
# Cette variable est utilisee dans les cas tests pour ecrire
# le resultat du test localement
#=============================================================================
def getLocal():
    """Return the local directory name to write test data."""
    a = os.getenv('VALIDLOCAL')
    if a is None: a = '.'
    return a

#=============================================================================
# Verifie que arrays est egal au arrays de reference stocke dans un fichier
# number est le no du test dans le script
#=============================================================================
def testA(arrays, number=1):
    """Test arrays."""
    import Converter as C
    if not isinstance(arrays[0], list): arrays = [arrays]

    # Check Data directory
    a = os.access(DATA, os.F_OK)
    if not a:
        print("{} directory doesn't exist. Created.".format(DATA))
        os.mkdir(DATA)

    # Construit le nom du fichier de reference
    fileName = sys.argv[0]
    baseName = os.path.basename(fileName)
    dirName = os.path.dirname(fileName)
    fileName = os.path.splitext(baseName)[0]

    if dirName == '': reference = '%s/%s.ref%d'%(DATA, fileName, number)
    else: reference = '%s/%s/%s.ref%d'%(dirName, DATA, fileName, number)
    a = os.access(reference, os.R_OK)
    if not a:
        print("Warning: reference file %s has been created."%reference)
        C.convertArrays2File(arrays, reference, 'bin_pickle')
        return True
    else:
        old = C.convertFile2Arrays(reference, 'bin_pickle')
        if GEOMETRIC_DIFF:
            # geometrical check
            if all(coord in oldArr[0].split(',') for oldArr in old for coord in 'xyz'):
                ret = C.diffArrayGeom(arrays, old, tol=TOLERANCE)
                if ret is None:
                    print('DIFF: geometrical diff, 1-to-1 match in '
                          'identifyNodes failed, topological diff performed '
                          'instead.')
                    ret = C.diffArrays(arrays, old)
            else:
                print("Warning: missing coordinates for geometrical diff., "
                      "topological diff performed instead.")
                ret = C.diffArrays(arrays, old)
        else:
            # topological check
            ret = C.diffArrays(arrays, old)

        isSuccessful = True
        varNames = list(dict.fromkeys(','.join(i[0] for i in ret).split(',')))
        nvarNames = len(varNames)
        l0 = [0. for _ in range(nvarNames)]
        l2 = [0. for _ in range(nvarNames)]
        for i in ret:
            for v in i[0].split(','):
                vidx = varNames.index(v)
                l0[vidx] = max(l0[vidx], C.normL0(i, v))
                l2[vidx] = max(l2[vidx], C.normL2(i, v))

        for vidx, v in enumerate(varNames):
            if l0[vidx] > TOLERANCE:
                print('DIFF: Variable=%s, L0=%.12f, L2=%.12f'%(v, l0[vidx], l2[vidx]))
                isSuccessful = False

        return isSuccessful

# idem testA avec ecriture fichier
def outA(arrays, number=1):
    """Test and write arrays."""
    import Converter as C
    C.convertArrays2File(arrays, 'out%d.plt'%number)
    testA(arrays)

#=============================================================================
# Verifie que le pyTree est egal au pyTree de reference stocke dans un fichier
# number est le no du test dans le script
#=============================================================================
def testT(t, number=1):
    """Test pyTrees."""
    import Converter.PyTree as C
    import Converter.Internal as Internal

    # Transforme t en pyTree, pour pouvoir relire la reference
    t, ntype = Internal.node2PyTree(t)

    # Verifie la compatibilite avec la CGNS lib
    #checkCGNSlib(t, number)

    # Check OWNDATA / copy
    C._ownNumpyArrays(t)

    # Check Data directory
    a = os.access(DATA, os.F_OK)
    if not a:
        print("{} directory doesn't exist. Created.".format(DATA))
        os.mkdir(DATA)

    # Construit le nom du fichier de reference
    fileName = sys.argv[0]
    baseName = os.path.basename(fileName)
    dirName = os.path.dirname(fileName)
    fileName = os.path.splitext(baseName)[0]

    if dirName == '': reference = '%s/%s.ref%d'%(DATA, fileName, number)
    else: reference = '%s/%s/%s.ref%d'%(dirName, DATA, fileName, number)
    a = os.access(reference, os.R_OK)
    if not a:
        print("Warning: reference file %s has been created."%reference)
        C.convertPyTree2File(t, reference, 'bin_pickle')
        return True
    else:
        old = C.convertFile2PyTree(reference, 'bin_pickle')
        checkTree(t, old)
        if GEOMETRIC_DIFF:
            # geometrical check
            nXYZ = [Internal.getNodesFromName(old, "Coordinate" + ax) for ax in 'XYZ']
            nXYZ = [arr[0] if len(arr) else [] for arr in nXYZ]
            if all(len(arr) for arr in nXYZ) and all(arr[1] is not None for arr in nXYZ):
                ret = C.diffArrayGeom(t, old, tol=TOLERANCE)
                if ret is None:
                    print('DIFF: geometrical diff, 1-to-1 match in '
                          'identifyNodes failed, topological diff performed '
                          'instead.')
                    ret = C.diffArrays(t, old)
            else:
                print("Warning: missing coordinates for geometrical diff., "
                      "topological diff performed instead.")
                ret = C.diffArrays(t, old)
        else:
            # topological check
            ret = C.diffArrays(t, old)
        C._fillMissingVariables(ret)
        mvars = C.getVarNames(ret)
        if len(mvars) > 0: mvars = mvars[0]
        else: mvars = []

        isSuccessful = True
        for v in mvars:
            l0 = C.normL0(ret, v)
            l2 = C.normL2(ret, v)
            if l0 > TOLERANCE:
                print('DIFF: Variable=%s, L0=%.12f, L2=%.12f'%(v,l0,l2))
                isSuccessful = False
        return isSuccessful

def outT(t, number=1):
    """Test and write pyTrees."""
    import Converter.PyTree as C
    C.convertPyTree2File(t, 'out%d.cgns'%number)
    testT(t, number)

#=============================================================================
# Verifie que le fichier est identique au fichier de reference
# Diff byte to byte
#=============================================================================
def testF(infile, number=1):
    # Chek infile
    a = os.access(infile, os.F_OK)
    if not a:
        print("DIFF: file "+infile+' doesnt exist.')

    # Check Data directory
    a = os.access(DATA, os.F_OK)
    if not a:
        print("{} directory doesn't exist. Created.".format(DATA))
        os.mkdir(DATA)
    fileName = sys.argv[0]
    baseName = os.path.basename(fileName)
    dirName = os.path.dirname(fileName)
    fileName = os.path.splitext(baseName)[0]

    if dirName == '': reference = '%s/%s.ref%d'%(DATA, fileName, number)
    else: reference = '%s/%s/%s.ref%d'%(dirName, DATA, fileName, number)
    a = os.access(reference, os.R_OK)
    if not a:
        print("Can not open file %s for reading."%reference)
        print("Reference file %s has been created."%reference)
        os.system("cp "+infile+" "+reference)
        return True
    else:
        print("Diffing with '"+reference+"'... done.")
        import filecmp
        ret = filecmp.cmp(reference, infile, shallow=False)
        #ret = os.system("diff "+reference+" "+infile)
        if not ret:
            print("DIFF: with file "+reference+'.')
            return False
        else: return True

def checkObject_(objet, refObjet, reference):
    # tests sur les types
    if type(refObjet) != type(objet):
        if isinstance(refObjet, (int, float)): pass
        elif (isinstance(refObjet, (numpy.int32, numpy.int64, numpy.intc)) and
              isinstance(objet, (numpy.int32, numpy.int64, numpy.intc))): pass
        elif (isinstance(refObjet, (numpy.float32, numpy.float64)) and
              isinstance(objet, (numpy.float32, numpy.float64))): pass
        else:
            print("DIFF: object type differs from "+reference+'.')
            return False
    # autres tests
    if isinstance(refObjet, bool):
        if refObjet != objet:
            print("DIFF: object value differs from %s (diff=%g,%g)."%(reference, objet, refObjet))
            return False
    elif isinstance(refObjet, (int, numpy.int32, numpy.int64, numpy.intc)):
        diff = abs(refObjet-objet)
        if diff > 0:
            print("DIFF: object value differs from %s (diff=%g)."%(reference, diff))
            return False
    elif isinstance(refObjet, (float, numpy.float32, numpy.float64)):
        diff = abs(refObjet-objet)
        if diff > TOLERANCE:
            print("DIFF: object value differs from %s (diff=%g)."%(reference, diff))
            return False
    elif isinstance(refObjet, dict):
        for k in refObjet.keys():
            v1 = refObjet[k]
            v2 = objet.get(k, None)
            ret = checkObject_(v1, v2, reference)
            if not ret: return False
    elif isinstance(refObjet, numpy.ndarray): # array
        if refObjet.shape != objet.shape:
            print("DIFF: object shape differs from "+reference+'.')
            return False
        if refObjet.dtype == 'S1' or objet.dtype == 'S1':
            if refObjet.dtype != objet.dtype:
                print("DIFF: object type differs from "+reference+'.')
            if not numpy.all(refObjet == objet):
                print("DIFF: object string differs from "+reference+'.')
                return False
        else:
            diff = numpy.abs(refObjet-objet)
            diff = (diff < TOLERANCE)
            if not diff.all():
                print("DIFF: object value differs from "+reference+'.')
                return False
    elif isinstance(refObjet, list): # liste
        for i, ai in enumerate(refObjet):
            if not checkObject_(objet[i], ai, reference): return False
    elif refObjet != objet: # autre objet
        print("DIFF: object differs from "+reference+'.')
        return False
    return True

#=============================================================================
# Verifie que l'objet python est identique a celui stocke dans le
# fichier de reference
#=============================================================================
def testO(objet, number=1):
    """Test python object."""
    if isinstance(objet, dict):
        # perform some sort on dict to be predictible
        from collections import OrderedDict
        objet = OrderedDict(sorted(objet.items(), key=lambda t: t[0]))

    # Check Data directory
    a = os.access(DATA, os.F_OK)
    if not a:
        print("{} directory doesn't exist. Created.".format(DATA))
        os.mkdir(DATA)
    fileName = sys.argv[0]
    baseName = os.path.basename(fileName)
    dirName = os.path.dirname(fileName)
    fileName = os.path.splitext(baseName)[0]

    if dirName == '': reference = '%s/%s.ref%d'%(DATA, fileName, number)
    else: reference = '%s/%s/%s.ref%d'%(dirName, DATA, fileName, number)
    a = os.access(reference, os.R_OK)

    # OWNDATA check / copy
    if isinstance(objet, numpy.ndarray) and not objet.flags['OWNDATA']:
        objet = numpy.copy(objet)

    if not a:
        print("Can not open file "+reference+" for reading.")
        print("Reference file "+reference+" has been created.")
        import pickle as pickle
        file = open(reference, 'wb')
        try: pickle.dump(objet, file, protocol=pickle.HIGHEST_PROTOCOL)
        except: pickle.dump('Undumpable object', file, protocol=pickle.HIGHEST_PROTOCOL)
        file.close()
        return True
    else:
        try: import pickle
        except ImportError: import cPickle as pickle
        file = open(reference, 'rb')
        oldData = False
        if oldData: a = pickle.load(file, encoding='latin1')
        else: a = pickle.load(file)
        file.close()
        print("Reading '"+reference+"'... done.")
        if isinstance(a, str) and a == 'Undumpable object': return True
        return checkObject_(objet, a, reference)

#=============================================================================
# Verifie que les arbres t1 et t2 sont identiques
# t1: courant; t2: reference
# Attention: seule la structure et les champs entiers sont testes
# Pour tester les champs reels, utiliser diffArrays
# Retourne 1 si identiques
# Retourne 0 si differents
#=============================================================================
def checkTree(t1, t2):
    """Check that pyTree t1 and t2 are identical."""
    dict1 = {}
    buildDict__('.', dict1, t1)
    dict1.pop('./CGNSTree/CGNSLibraryVersion', None) # avoid comparison when version change
    dict2 = {}
    buildDict__('.', dict2, t2)
    dict2.pop('./CGNSTree/CGNSLibraryVersion', None)

    for k in dict2.keys():
        node2 = dict2[k]
        # cherche le noeud equivalent dans t1
        if k not in dict1:
            print('DIFF: node %s existe dans reference mais pas dans courant.'%k)
        else:
            node1 = dict1[k]
            checkTree__(node1, node2)

def buildDict__(curr, mdict, node):
    d = curr+'/'+node[0]
    mdict[d] = node
    for i in node[2]: buildDict__(d, mdict, i)

def checkTree__(node1, node2):
    if node1[0] != node2[0]: # nom du noeud
        print('DIFF: nom des noeuds differents:')
        print('DIFF: reference: %s.'%node2[0])
        print('DIFF: courant: %s.'%node1[0])
        return 0
    if node1[3] != node2[3]: # type du noeud
        print('DIFF: type de noeud differents pour le noeud: %s.'%node1[0])
        print('DIFF: reference: %s.'%node2[3])
        print('DIFF: courant: %s.'%node1[3])
        return 0
    if GEOMETRIC_DIFF and node1[0] in ['NGonElements', 'NFaceElements']:
        return 1
    if len(node1[2]) != len(node2[2]):
        childNames1 = [n[0] for n in node1[2]]
        childNames2 = [n[0] for n in node2[2]]
        childNamesSet1 = set(childNames1)
        childNamesSet2 = set(childNames2)
        diffSet12 = childNamesSet1 - childNamesSet2
        diffSet21 = childNamesSet2 - childNamesSet1
        print('DIFF: longueur des fils differente pour le noeud: %s.'%node1[0])
        if len(diffSet12) > 0:
            print('  - Noms des noeuds de courant qui ne sont pas dans '\
                  'ref:\n{}.'.format(', '.join(f'{i}' for i in diffSet12)))
        if len(diffSet21) > 0:
            print('  - Noms des noeuds de ref qui ne sont pas dans '\
                  'courant:\n{}.'.format(', '.join(f'{i}' for i in diffSet21)))
        if len(diffSet12) == 0 and len(diffSet21) == 0:
            from collections import Counter
            if len(childNamesSet1) != len(childNames1):
                doublons = set(i for i, count in Counter(childNames1).items() if count > 1)
                print('  - Noeuds doublons detectes dans courant: {}.'.format(
                      ', '.join(f'{i}' for i in doublons)))
            if len(childNamesSet2) != len(childNames2):
                doublons = set(i for i, count in Counter(childNames2).items() if count > 1)
                print('  - Noeuds doublons detectes dans ref: {}.'.format(
                      ', '.join(f'{i}' for i in doublons)))
        return 0
    if GEOMETRIC_DIFF and node1[0] in ['ElementRange', 'ElementConnectivity']:
        return 1
    val1 = node1[1]; val2 = node2[1]
    if isinstance(val1, str):
        if not isinstance(val2, str):
            print('DIFF: types de valeurs differents pour le noeud: %s.'%node1[0])
            print('DIFF: reference: {}'.format(type(val2)))
            print('DIFF: courant: str')
            return 0
        if val1 != val2:
            print('DIFF: valeurs differentes pour le noeud: %s.'%node1[0])
            print('DIFF: reference:'+val2)
            print('DIFF: courant:'+val1)
            return 0
    elif isinstance(val1, float):
        if not isinstance(val2, float):
            print('DIFF: types de valeurs differents pour le noeud: %s.'%node1[0])
            print('DIFF: reference: {}'.format(type(val2)))
            print('DIFF: courant: float')
            return 0
        if val1 != val2:
            print('DIFF: valeurs differentes pour le noeud: %s.'%node1[0])
            print('DIFF: reference: %f'%val2)
            print('DIFF: courant: %f'%val1)
            return 0
    elif isinstance(val1, int):
        if not isinstance(val2, int):
            print('DIFF: types de valeurs differents pour le noeud: %s.'%node1[0])
            print('DIFF: reference: {}'.format(type(val2)))
            print('DIFF: courant: int')
            return 0
        if val1 != val2:
            print('DIFF: valeurs differentes au noeud:'%node1[0])
            print('DIFF: reference: %d'%val2)
            print('DIFF: courant: %d'%val1)
            return 0
    elif isinstance(val1, numpy.ndarray):
        if not isinstance(val2, numpy.ndarray):
            print('DIFF: types numpy de valeurs differents pour le noeud: %s.'%node1[0])
            print('DIFF: reference:', val2.dtype)
            print('DIFF: courant:', val1.dtype)
            print('DIFF: reference:', val2)
            print('DIFF: courant:', val1)
            return 0
        if val1.dtype == numpy.int32 or val1.dtype == numpy.int64:
            if val2.dtype != numpy.int32 and val2.dtype != numpy.int64:
                print('DIFF: types numpy de valeurs differents pour le noeud: %s.'%node1[0])
                print('DIFF: reference:', val2)
                print('DIFF: courant:', val1)
                return 0
            if val1.shape != val2.shape:
                print('DIFF: shape differentes pour le noeud: %s.'%node1[0])
                print('DIFF: reference:', val2.shape)
                print('DIFF: courant:', val1.shape)
                print('DIFF: reference:', val2)
                print('DIFF: courant:', val1)
                return 0
            if ((val1 == val2).all()) == False:
                print('DIFF: valeurs differentes pour le noeud: %s.'%node1[0])
                print('DIFF: reference:', val2)
                print('DIFF: courant:', val1)
                return 0
        if val1.dtype == numpy.float64:
            if val2.dtype != numpy.float64:
                print('DIFF: types numpy de valeurs differents pour le noeud: %s.'%node1[0])
                print('DIFF: reference:', val2.dtype)
                print('DIFF: courant:', val1.dtype)
                print('DIFF: reference:', val2)
                print('DIFF: courant:', val1)
                return 0
            if val1.shape != val2.shape:
                print('DIFF: shape differentes pour le noeud: %s.'%node1[0])
                print('DIFF: reference:', val2.shape)
                print('DIFF: courant:', val1.shape)
                print('DIFF: reference:', val2)
                print('DIFF: courant:', val1)
                return 0
        if val1.dtype == numpy.float32:
            if val2.dtype != numpy.float32:
                print('DIFF: types numpy de valeurs differents pour le noeud: %s.'%node1[0])
                print('DIFF: reference:', val2.dtype)
                print('DIFF: courant:', val1.dtype)
                print('DIFF: reference:', val2)
                print('DIFF: courant:', val1)
                return 0
            if val1.shape != val2.shape:
                print('DIFF: shape differentes pour le noeud: %s.'%node1[0])
                print('DIFF: reference:', val2.shape)
                print('DIFF: courant:', val1.shape)
                print('DIFF: reference:', val2)
                print('DIFF: courant:', val1)
                return 0

        #     if (numpy.abs(val1 -val2)<1.e-6).all() == False:
        #         print('DIFF: valeurs differentes pour le noeud: %s.'%node1[0])
        #         delta = numpy.max(numpy.abs(val1 -val2))
        #         print('DIFF: ', delta)
        #         return 0
    return 1

#==============================================================================
# tests standards
#==============================================================================
# Batterie de tests standards pour les arrays
#==============================================================================
def stdTestA(F, *keywords):
    ntype = 'std'
    if len(sys.argv) >= 2: ntype = sys.argv[1]
    if ntype == 'std': stdTest1__(0, 0, 0, F, *keywords)
    elif ntype == 'out': stdTest1__(1, 0, 0, F, *keywords)
    elif ntype == 'mem': stdTest1__(0, 10, 0, F, *keywords)
    elif ntype == 'heavy': stdTest1__(0, 0, 1, F, *keywords)

# Std tests avec ecriture de fichiers
def outTestA(F, *keywords):
    stdTest1__(1, 0, 0, F, *keywords)
# Std tests avec tests memoire
def memTestA(loop, F, *keywords):
    stdTest1__(0, loop, 0, F, *keywords)
# Std tests avec heavy memory load
def heavyTestA(F, *keywords):
    stdTest1__(0, 0, 1, F, *keywords)

#==============================================================================
# retourne 0: a est un array
# retourne 1: a est un arrays
# retourne 2: sinon
#==============================================================================
def checkType__(a):
    if isinstance(a, list):
        l = len(a)
        if l == 0: return 2
        if isinstance(a[0], str) and (l == 4 or l == 5):
            return 0
        else:
            b = a[0]
            if not isinstance(b, list): return 2
            if len(b) == 0: return 2
            if isinstance(b[0], str) and (len(b) == 4 or len(b) == 5):
                return 1
            else: return 2
    else: return 2

#==============================================================================
# IN: output=1: ecrit les fichiers resultats
# IN: memory>0: test les fuites memoire en bouclant sur les fonctions
# IN: heavy=1: tests lourds (gros maillages, plein de variables)
#==============================================================================
def stdTest1__(output, memory, heavy, F, *keywords):
    import Converter as C
    import Generator as G

    testName = sys.argv[0]
    coverage = 0

    # 1- Structure 1D
    try:
        a = G.cart((0,0,0), (1,1,1), (10,1,1))
        C._initVars(a, '{F}={x}+{y}+{z}')
        b = F(a, *keywords)
        res = checkType__(b)
        if output == 1 and res == 0: C.convertArrays2File([b], 'out1.plt')
        elif output == 1 and res == 1: C.convertArrays2File(b, 'out1.plt')
        elif output == 1 and res == 2: print('STRUCT1D', b)
        if res == 0: testA([b], 1)
        elif res == 1: testA(b, 1)
        elif res == 2: testO(b, 1)
        for i in range(memory): b = F(a, *keywords)
        if heavy == 1:
            a = G.cart((0,0,0), (1,1,1), (10000,1,1) )
            for i in range(100): C._addVars(a, 'F'+str(i))
            b = F(a, *keywords)
        coverage += 1
    except TypeError:
        if output == 1: print('STRUCT1D: uncovered.')
    except: print('%s: Structure 1D: fails.'%testName); raise

    # 2- Structure 2D
    try:
        a = G.cart((0,0,0), (1,1,1), (10,10,1))
        C._initVars(a, '{F}={x}+{y}+{z}')
        b = F(a, *keywords)
        res = checkType__(b)
        if output == 1 and res == 0: C.convertArrays2File([b], 'out2.plt')
        elif output == 1 and res == 1: C.convertArrays2File(b, 'out2.plt')
        elif output == 1 and res == 2: print('STRUCT2D', b)
        if res == 0: testA([b], 2)
        elif res == 1: testA(b, 2)
        elif res == 2: testO(b, 2)
        for i in range(memory): b = F(a, *keywords)
        if heavy == 1:
            a = G.cart((0,0,0), (1,1,1), (1000,1000,1) )
            for i in range(100): C._addVars(a, 'F'+str(i))
            b = F(a, *keywords)
        coverage += 1
    except TypeError:
        if output == 1: print('STRUCT2D: uncovered.')
    except: print('%s: Structure 2D: fails.'%testName); raise

    # 3- Structure 3D
    try:
        a = G.cart((0,0,0), (1,1,1), (10,10,10) )
        C._initVars(a, '{F}={x}+{y}+{z}')
        b = F(a, *keywords)
        res = checkType__(b)
        if output == 1 and res == 0: C.convertArrays2File([b], 'out3.plt')
        elif output == 1 and res == 1: C.convertArrays2File(b, 'out3.plt')
        elif output == 1 and res == 2: print('STRUCT3D', b)
        if res == 0: testA([b], 3)
        elif res == 1: testA(b, 3)
        elif res == 2: testO(b, 3)
        for i in range(memory): b = F(a, *keywords)
        if heavy == 1:
            a = G.cart((0,0,0), (1,1,1), (1000,1000,100) )
            for i in range(100): C._addVars(a, 'F'+str(i))
            b = F(a, *keywords)
        coverage += 1
    except TypeError:
        if output == 1: print('STRUCT3D: uncovered.')
    except: print('%s: Structure 3D: fails.'%testName); raise

    # 4- BAR
    try:
        a = G.cartTetra((0,0,0), (1,1,1), (10,1,1))
        C._initVars(a, '{F}={x}+{y}+{z}')
        b = F(a, *keywords)
        res = checkType__(b)
        if output == 1 and res == 0: C.convertArrays2File([b], 'out4.plt')
        elif output == 1 and res == 1: C.convertArrays2File(b, 'out4.plt')
        elif output == 1 and res == 2: print('BAR', b)
        if res == 0: testA([b], 4)
        elif res == 1: testA(b, 4)
        elif res == 2: testO(b, 4)
        for i in range(memory): b = F(a, *keywords)
        if heavy == 1:
            a = G.cartTetra((0,0,0), (1,1,1), (10000,1,1) )
            for i in range(100): C._addVars(a, 'F'+str(i))
            b = F(a, *keywords)
        coverage += 1
    except TypeError:
        if output == 1: print('BAR: uncovered.')
    except: print('%s: BAR: fails.'%testName); raise

    # 5- TRI
    try:
        a = G.cartTetra((0,0,0), (1,1,1), (10,10,1))
        C._initVars(a, '{F}={x}+{y}+{z}')
        b = F(a, *keywords)
        res = checkType__(b)
        if output == 1 and res == 0: C.convertArrays2File([b], 'out5.plt')
        elif output == 1 and res == 1: C.convertArrays2File(b, 'out5.plt')
        elif output == 1 and res == 2: print('TRI', b)
        if res == 0: testA([b], 5)
        elif res == 1: testA(b, 5)
        elif res == 2: testO(b, 5)
        for i in range(memory): b = F(a, *keywords)
        if heavy == 1:
            a = G.cartTetra((0,0,0), (1,1,1), (1000,1000,1))
            for i in range(100): C._addVars(a, 'F'+str(i))
            b = F(a, *keywords)
        coverage += 1
    except TypeError:
        if output == 1: print('TRI: uncovered.')
    except: print('%s: TRI: fails.'%testName); raise

    # 6- QUAD
    try:
        a = G.cartHexa((0,0,0), (1,1,1), (10,10,1))
        C._initVars(a, '{F}={x}+{y}+{z}')
        b = F(a, *keywords)
        res = checkType__(b)
        if output == 1 and res == 0: C.convertArrays2File([b], 'out6.plt')
        elif output == 1 and res == 1: C.convertArrays2File(b, 'out6.plt')
        elif output == 1 and res == 2: print('QUAD', b)
        if res == 0: testA([b], 6)
        elif res == 1: testA(b, 6)
        elif res == 2: testO(b, 6)
        for i in range(memory): b = F(a, *keywords)
        if heavy == 1:
            a = G.cartHexa((0,0,0), (1,1,1), (1000,1000,100))
            for i in range(100): C._addVars(a, 'F'+str(i))
            b = F(a, *keywords)
        coverage += 1
    except TypeError:
        if output == 1: print('QUAD: uncovered.')
    except: print('%s: QUAD: fails.'%testName); raise

    # 7- TETRA
    try:
        a = G.cartTetra((0,0,0), (1,1,1), (10,10,10))
        C._initVars(a, '{F}={x}+{y}+{z}')
        b = F(a, *keywords)
        res = checkType__(b)
        if output == 1 and res == 0: C.convertArrays2File([b], 'out7.plt')
        elif output == 1 and res == 1: C.convertArrays2File(b, 'out7.plt')
        elif output == 1 and res == 2: print('TETRA', b)
        if res == 0: testA([b], 7)
        elif res == 1: testA(b, 7)
        elif res == 2: testO(b, 7)
        for i in range(memory): b = F(a, *keywords)
        if heavy == 1:
            a = G.cartTetra((0,0,0), (1,1,1), (1000,1000,100) )
            for i in range(100): C._addVars(a, 'F'+str(i))
            b = F(a, *keywords)
        coverage += 1
    except TypeError:
        if output == 1: print('TETRA: uncovered.')
    except: print('%s: TETRA: fails.'%testName); raise

    # 8- HEXA
    try:
        a = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
        C._initVars(a, '{F}={x}+{y}+{z}')
        b = F(a, *keywords)
        res = checkType__(b)
        if output == 1 and res == 0: C.convertArrays2File([b], 'out8.plt')
        elif output == 1 and res == 1: C.convertArrays2File(b, 'out8.plt')
        elif output == 1 and res == 2: print('HEXA', b)
        if res == 0: testA([b], 8)
        elif res == 1: testA(b, 8)
        elif res == 2: testO(b, 8)
        for i in range(memory): b = F(a, *keywords)
        if heavy == 1:
            a = G.cartHexa((0,0,0), (1,1,1), (1000,1000,100) )
            for i in range(100): C._addVars(a, 'F'+str(i))
            b = F(a, *keywords)
        coverage += 1
    except TypeError:
        if output == 1: print('HEXA: uncovered.')
    except: print('%s: HEXA: fails.'%testName); raise

    # 9- PENTA
    try:
        a = G.cartPenta((0,0,0), (1,1,1), (10,10,10))
        C._initVars(a, '{F}={x}+{y}+{z}')
        b = F(a, *keywords)
        res = checkType__(b)
        if output == 1 and res == 0: C.convertArrays2File([b], 'out9.plt')
        elif output == 1 and res == 1: C.convertArrays2File(b, 'out9.plt')
        elif output == 1 and res == 2: print('PENTA', b)
        if res == 0: testA([b], 9)
        elif res == 1: testA(b, 9)
        elif res == 2: testO(b, 9)
        for i in range(memory): b = F(a, *keywords)
        if heavy == 1:
            a = G.cartPenta((0,0,0), (1,1,1), (1000,1000,100) )
            for i in range(100): C._addVars(a, 'F'+str(i))
            b = F(a, *keywords)
        coverage += 1
    except TypeError:
        if output == 1: print('PENTA: uncovered.')
    except: print('%s: PENTA: fails.'%testName); raise

    # 10- PYRA
    try:
        a = G.cartPyra( (0,0,0), (1,1,1), (10,10,10) )
        C._initVars(a, '{F}={x}+{y}+{z}')
        b = F(a, *keywords)
        res = checkType__(b)
        if output == 1 and res == 0: C.convertArrays2File([b], 'out10.plt')
        elif output == 1 and res == 1: C.convertArrays2File(b, 'out10.plt')
        elif output == 1 and res == 2: print('PYRA', b)
        if res == 0: testA([b], 10)
        elif res == 1: testA(b, 10)
        elif res == 2: testO(b, 10)
        for i in range(memory): b = F(a, *keywords)
        if heavy == 1:
            a = G.cartPyra( (0,0,0), (1,1,1), (1000,1000,100) )
            for i in range(100): C._addVars(a, 'F'+str(i))
            b = F(a, *keywords)
        coverage += 1
    except TypeError:
        if output == 1: print('PYRA: uncovered.')
    except: print('%s: PYRA: fails.'%testName); raise

    # 11- NGON 1D
    try:
        a = G.cartNGon((0,0,0), (1,1,1), (10,1,1))
        C._initVars(a, '{F}={x}+{y}+{z}')
        b = F(a, *keywords)
        res = checkType__(b)
        if output == 1 and res == 0: C.convertArrays2File([b], 'out11.tp')
        elif output == 1 and res == 1: C.convertArrays2File(b, 'out11.tp')
        elif output == 1 and res == 2: print('NGON1D', b)
        if res == 0: testA([b], 11)
        elif res == 1: testA(b, 11)
        elif res == 2: testO(b, 11)
        for i in range(memory): b = F(a, *keywords)
        if heavy == 1:
            a = G.cartNGon((0,0,0), (1,1,1), (10000,1,1) )
            for i in range(100): C._addVars(a, 'F'+str(i))
            b = F(a, *keywords)
        coverage += 1
    except TypeError:
        if output == 1: print('NGON1D: uncovered.')
    except: print('%s: NGON 1D: fails.'%testName); raise

    # 12- NGON 2D
    try:
        a = G.cartNGon((0,0,0), (1,1,1), (10,10,1))
        C._initVars(a, '{F}={x}+{y}+{z}')
        b = F(a, *keywords)
        res = checkType__(b)
        if output == 1 and res == 0: C.convertArrays2File([b], 'out12.tp')
        elif output == 1 and res == 1: C.convertArrays2File(b, 'out12.tp')
        elif output == 1 and res == 2: print('NGON2D', b)
        if res == 0: testA([b], 12)
        elif res == 1: testA(b, 12)
        elif res == 2: testO(b, 12)
        for i in range(memory): b = F(a, *keywords)
        if heavy == 1:
            a = G.cartNGon((0,0,0), (1,1,1), (1000,1000,1) )
            for i in range(100): C._addVars(a, 'F'+str(i))
            b = F(a, *keywords)
        coverage += 1
    except TypeError:
        if output == 1: print('NGON2D: uncovered.')
    except: print('%s: NGON 2D: fails.'%testName); raise

    # 13- NGON 3D
    try:
        a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
        C._initVars(a, '{F}={x}+{y}+{z}')
        b = F(a, *keywords)
        res = checkType__(b)
        if output == 1 and res == 0: C.convertArrays2File([b], 'out13.plt')
        elif output == 1 and res == 1: C.convertArrays2File(b, 'out13.plt')
        elif output == 1 and res == 2: print('NGON3D', b)
        if res == 0: testA([b], 13)
        elif res == 1: testA(b, 13)
        elif res == 2: testO(b, 13)
        for i in range(memory): b = F(a, *keywords)
        if heavy == 1:
            a = G.cartNGon((0,0,0), (1,1,1), (1000,1000,100) )
            for i in range(100): C._addVars(a, 'F'+str(i))
            b = F(a, *keywords)
        coverage += 1
    except TypeError:
        if output == 1: print('NGON3D: uncovered.')
    except: print('%s: NGON 3D: fails.'%testName); raise

    # 14- liste d'arrays
    try:
        a = G.cart((0,0,0), (1,1,1), (10,10,10))
        C._initVars(a, '{F}={x}+{y}+{z}')
        b = G.cartTetra((0,0,0), (1,1,1), (10,10,10))
        C._initVars(b, '{F}={x}+{y}+{z}')
        A = F([a,b], *keywords)
        res = checkType__(A)
        if output == 1 and res == 0: C.convertArrays2File([A], 'out14.plt')
        elif output == 1 and res == 1: C.convertArrays2File(A, 'out14.plt')
        elif output == 1 and res == 2: print('List of arrays', b)
        if res == 0: testA([A], 14)
        elif res == 1: testA(A, 14)
        elif res == 2: testO(A, 14)
        coverage += 1
    except TypeError:
        if output == 1: print('Array list: uncovered.')
    except: print('%s: Array list: fails.'%testName); raise

    # Write coverage
    writeCoverage(coverage/14.*100)

#==============================================================================
def stdTestT(F, *keywords):
    type = 'std'
    if len(sys.argv) >= 2: type = sys.argv[1]
    if type == 'std': stdTestT__(0, F, *keywords)
    elif type == 'out': stdTestT__(1, F, *keywords)

# standard tests avec ecriture
def outTestT(F, *keywords):
    stdTestT__(1, F, *keywords)

#==============================================================================
def stdTestT__(output, F, *keywords):
    import Converter.PyTree as C
    import Generator.PyTree as G
    import Converter.Internal as Internal

    testName = sys.argv[0]
    coverage = 0

    # 1- Une zone
    try:
        a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
        C._initVars(a, '{F}={CoordinateX}+{CoordinateY}+{CoordinateZ}')
        C._initVars(a, '{centers:G}={centers:CoordinateX}+{centers:CoordinateY}+{centers:CoordinateZ}')
        C._addBC2Zone(a, 'wall', 'BCWall', 'imin')
        C._addBC2Zone(a, 'overlap', 'BCOverlap', 'jmin')
        b = F(a, *keywords)
        res = Internal.isStdNode(b)
        if output == 1 and res != -2:
            t = C.newPyTree(['Base']); t[2][1][2] += [b]
            C.convertPyTree2File(t, 'out1.cgns')
        if res != -2: testT(b, 1)
        else: testO(b, 1)
        coverage += 1
    except TypeError: pass #
    except: print('%s: One zone: fails.'%testName); raise

    # 2- Une liste de zones
    try:
        a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
        b = G.cartTetra( (0,0,0), (1,1,1), (10,10,10) )
        B = F([a,b], *keywords)
        res = Internal.isStdNode(B)
        if output == 1 and res != -2:
            t = C.newPyTree(['Base']); t[2][1][2] += B
            C.convertPyTree2File(t, 'out2.cgns')
        if res != -2: testT(B, 2)
        else: testO(B, 2)
        coverage += 1
    except TypeError: pass
    except: print('%s: Zone list: fails.'%testName); raise

    # 3- Un arbre
    try:
        a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
        t = C.newPyTree(['Base'])
        C._addState(t[2][1], 'EquationDimension', 3)
        t[2][1][2] += [a]
        t = F(t, *keywords)
        res = Internal.isStdNode(t)
        if output == 1 and res != -2: C.convertPyTree2File(t, 'out3.cgns')
        if res != -2: testT(t, 3)
        else: testO(t, 3)
        coverage += 1
    except TypeError: pass
    except: print('%s: Full tree: fails.'%testName); raise

    # 4- Une base
    try:
        a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
        t = C.newPyTree(['Base'])
        C._addState(t[2][1], 'EquationDimension', 3)
        t[2][1][2] += [a]
        b = t[2][1]
        b = F(b, *keywords)
        t[2][1] = b
        res = Internal.isStdNode(t)
        if output == 1 and res != -2: C.convertPyTree2File(t, 'out4.cgns')
        if res != -2: testT(b, 4)
        else: testO(b, 4)
        coverage += 1
    except TypeError: pass #
    except: print('%s: base: fails.'%testName); raise

##     # 5- Une liste de bases
##     try:
##         a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
##         b = G.cart( (12,0,0), (1,1,1), (10,10,10) )
##         t = C.newPyTree(['Base', 'Base2'])
##         C._addState(t[2][1], 'EquationDimension', 3)
##         t[2][1][2] += [a]; t[2][2][2] += [b]
##         bases = [ t[2][1], t[2][2] ]
##         #bases = F(bases, *keywords)
##         t = C.newPyTree(); t[2] += bases
##         if (output == 1): C.convertPyTree2File(t, 'out5.cgns')
##         testT(bases, 5)
##         coverage += 1
##     except TypeError: pass
##     except: print('%s: list of bases: fails.'%testName); raise

    # Write coverage
    writeCoverage(coverage/4.*100)

#==============================================================================
# Ecrit la couverture du test a l'ecran
#==============================================================================
def writeCoverage(coverage):
    testName = sys.argv[0]
    print('%s: coverage=%d'%(testName, coverage)+'%')

#==============================================================================
# Verifie la validite du fichier avec la cgnslib
# Retourne 0 (FAIL), 1 (SUCCESS)
#==============================================================================
def checkCGNSlib(t, number=1):
    import Converter.PyTree as C
    C.convertPyTree2File(t, '.test%d.cgns'%number)
    import subprocess
    cmd = "cgnscheck .test%d.cgns"%number
    try:
        s = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        # Look for errors in output
        i1 = s.find("ERROR:")
        if i1 != -1: print("FAILED: CGNSlib check."); return 0
        return 1
    except: print("FAILED: CGNSlib check."); return 0

#==============================================================================
# Ecrit la memoire prise par le process
# IN: msg: message a ecrire en meme temps que la memoire
# Return the memory in kB
#==============================================================================
def printMem(msg, waitTime=0.1):
    """Write process memory to stdout."""
    import time
    pid = os.getpid()
    time.sleep(waitTime)
    try: f = open("/proc/{}/smaps".format(pid))
    except:
        #f = open("/proc/{}/status".format(pid))
        return 0.
    s = f.readlines()
    f.close()

    tot = 0.
    found = False
    for ts in s:
        if found:
            tot += int(ts[5:-3])
            found = False
        if ts.find("heap") >= 0:
            found = True
    if tot > 1.e6:
        print('{:<40} : {} GB '.format(msg,tot/1.e6))
    elif tot > 1000.:
        print('{:<40} : {} MB '.format(msg,tot/1000.))
    else: print('{:<40} : {} kB '.format(msg,tot))
    sys.stdout.flush()
    return tot
