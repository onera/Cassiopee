import re as RE
import os
import shutil
exp = RE.compile("#include")

MODS = ['TKBool', 'TKBRep', 'TKernel',
        'TKG2d', 'TKG3d', 'TKGeomAlgo',
        'GeomBase', 'TKIGES', 'TKMath', 'TKPrim',
        'TKShHealing', 'TKTopAlgo', 'TKXSBase']

DIRS = ['.', 'occ_inc']
for i in MODS: DIRS += ['occ_src/'+i]

# Le dictionnaire qui garde trace des fichiers deja visites
# file:0 (a ouvrir), file:1 (deja visite)
dic = {}

# cherche les includes: ajoute le fichier cxx + le fichier include
def checkDepend(file, dic):
    #print('Checking %s'%file)
    dic[file] = 1
    f = None
    for i in DIRS:
        try:
            f = open(i+'/'+file, "r")
            print('Found %s.'%(i+'/'+file))
            break
        except: pass
    if f is None:
        print('File %s not found.'%file); return

    lines = f.readlines()
    for i in lines:
        r = exp.search(i)
        if r is not None:
            i = i.replace("#include", "")
            i = i.replace("\"","")
            i = i.replace(" ","")
            i = i.replace("\n","")
            i = i.replace("\r","")
            i = i.replace("<","")
            i = i.replace(">","")
            if i not in dic: dic[i] = 0 # add include file

            c = i.replace("hxx", "cxx")
            if c not in dic: dic[c] = 0 # add cxx file
            d = i.replace("hxx", "dxx")
            if c not in dic: dic[c] = 0 # add cxx file
    f.close()

#==============================================================================
# Reconstruit un repertoire identique avec uniquement les bons fichiers
#==============================================================================
def rebuild(dic):
    keys = dic.keys()
    names = []
    for k in keys:
        # find orig
        name = None
        for i in DIRS:
            if os.access(i+'/'+k, os.R_OK):
                #print('Found %s.'%(i+'/'+k))
                name = i+'/'+k
                names.append(name)
                break
    for n in names:
        target = n.replace('occ_src', 'OCC_SRC')
        target = target.replace('occ_inc', 'OCC_INC')
        base = n.split('/')[0]
        if (base != '.'):
            shutil.copyfile(n, target);
            print('copy %s -> %s'%(n, target))

#==============================================================================
# Ecrit les sources (.cxx) par modules
#==============================================================================
def writeSources(dic):
    keys = dic.keys()
    print('== Sources (%d) =='%len(keys))
    sources = {}
    for i in MODS: sources[i] = []

    names = []
    for k in keys:
        # find orig
        name = None
        for i in DIRS:
            if os.access(i+'/'+k, os.R_OK):
                #print('Found %s.'%(i+'/'+k))
                name = i+'/'+k
                names.append(name)
                break

    for n in names:
        spl = n.split('/')
        last = spl[-1]
        last = last.split('.')
        if len(last) == 2 and last[1] == 'cxx':
            for i in MODS:
                if spl[1] == i: sources[i] += [n]

    for i in MODS:
        print('%s_srcs = ['%i)
        for k in sources[i]: print('%s,'%k.replace('occ_src', 'OCC_SRC'))
        print(' ]')
#==============================================================================

checkDepend('import_OCC_CAD_wrapper.cpp', dic)
checkDepend('CADviaOCC.cpp', dic)
checkDepend('OCCSurface.cpp', dic)

end = False
while end == False:
    keys = dic.keys() ; l = len(keys) ; update = 0
    for k in keys:
        if (dic[k] == 0): checkDepend(k, dic) ; update = 1
    if update == 0: end = True

rebuild(dic)
writeSources(dic)
