# ice mesher 2D
import matplotlib
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
matplotlib.rcParams['font.size'] = 10

import os
os.environ['MPLCONFIGDIR'] = '/tmp'

try:
    import CPlot.Decorator as _Dec

    def _safe_createSubPlot(*args, **kwargs):
        import matplotlib.pyplot as plt
        figsize = kwargs.get('figsize', (10, 6))
        fig, ax = plt.subplots(figsize=figsize)
        return fig, ax
    _Dec.createSubPlot = _safe_createSubPlot

    def _safe_savefig(path, *args, **kwargs):
        import matplotlib.pyplot as plt
        try:
            plt.savefig(path, dpi=150, bbox_inches='tight')
        except Exception as e:
            print(f"  [savefig] ⚠️ {e}")
    _Dec.savefig = _safe_savefig

    def _safe_closeAll(*args, **kwargs):
        import matplotlib.pyplot as plt
        plt.close('all')
    _Dec.closeAll = _safe_closeAll

    print("[PATCH] CPlot.Decorator patché avec succès")
except ImportError:
    print("[PATCH] CPlot.Decorator non disponible — ignoré")

import Converter.PyTree as C
import Transform.PyTree as T
import Generator.PyTree as G
import Post.PyTree as P
import Geom.PyTree as D
import Converter.Internal as Internal
import Dist2Walls.PyTree as D2W

import matplotlib.pyplot as plt
import numpy as np

factorThreshold = dict(
    minF = 0.1,
    maxF = 6.0,
    minSweeps = 0,
    )

DEFAULT_MESH_QUALITY = dict(
    maxAcrit = 10,     # no bad-angle cells
    minVBL = 1e-10,
    minVT3= 1e-8,
    maxVcrit = 0,     # no bad-volume cells
    minAngleBL = 0.0,  # minimum angle in BL quad mesh [°]
    minAngleT3 = 5.0,   # minimum angle in T3 mesh [°]
    maxRBL = 25.0,
    maxRT3 = 25.0
    )

def defaultGeoQuality () :
    return dict(
        maxInward  = 0.005,   # 0.5 % corde
        maxOutward = 0.005,
        stdDev     = 0.005,
        meanDev    = 0.005
    )

nptsName = 2500
SWEEPSname = 50

def exportName (name, nptsNameGiven):
    global SWEEPSname, nptsName
    nptsName = nptsNameGiven
    SWEEPSname = name

#===============================================================================================
# Updates the mesh threshold for BL and T3 vmin according to hf and spacing between points
#================================================================================================
def updateMeshThreshold (BAR1, hf, dz):
    global DEFAULT_MESH_QUALITY

    secu = 0.01
    x, y, _ = getCoords(C.convertBAR2Struct(BAR1))
    
    # Calcul des espacements entre chaque point consécutif
    dl = np.sqrt(np.diff(x)**2 + np.diff(y)**2)   
    mindl = np.min(dl)
    # print (f'mindl = {mindl}')
    
    vminBL = hf * mindl * secu
    vminT3 = 0.433 * mindl**2 * dz * secu

    oldBL = DEFAULT_MESH_QUALITY['minVBL'] 
    oldT3 = DEFAULT_MESH_QUALITY['minVT3']

    if vminBL > 0:
        vminBL = 10 ** np.floor(np.log10(vminBL))
        DEFAULT_MESH_QUALITY['minVBL'] = vminBL

    if vminT3 > 0:
        vminT3 = 10 ** np.floor(np.log10(vminT3))   
        DEFAULT_MESH_QUALITY['minVT3'] = vminT3

    print (f"new vmin= {DEFAULT_MESH_QUALITY['minVBL']} (vs old {oldBL}) et vminT3 = {DEFAULT_MESH_QUALITY['minVT3']} (vs old {oldT3})")

#===============================================================================================
# checks the profile order (1D STRUCT)
# returns 1 if normal outside and 0 otherwise
#================================================================================================
def checkOrder(profile):
    distrib = D.line((0,0,0), (0.01,0,0), N=2)
    lay = G.addNormalLayers(profile, distrib, niter=0, check=False)
    lay1 = T.subzone(lay, (1,-1,1), (-1,-1,1))
    l0 = D.getLength(profile)
    l1 = D.getLength(lay1)
    if l1 > l0: return 1
    else: return 0

#===============================================================================================
# generaters a TXT file with the BC markers to copy paste into solver (fsui)
#===============================================================================================
def genBCMarkers (m2, profileName):
        symBCs = [bc for bc in Internal.getNodesFromType(m2, 'BC_t') 
            if Internal.getValue(bc) == 'BCSymmetryPlane']
        nSym = len(symBCs)

        with open(f'Solveur/Résolution/{profileName}/0_{nptsName}/info/bc_markers.txt', 'w') as f:
            f.write(f"# Généré automatiquement - {nSym} sym détectées\n")
            f.write("[Marker]\n")
            f.write("    aircraft   = 1\n")
            f.write("    farfield   = 2\n")
            for i in range(nSym):
                f.write(f"    sym{i:<8}= {i+3}\n")
            f.write("[Families]\n")
            f.write("    Aircraft   = aircraft\n")
            f.write("    Farfield   = farfield\n")
            for i in range(nSym):
                f.write(f"    Symmetry{i:<3}= sym{i}\n")
            f.write("[BoundaryConditions]\n")
            f.write("    [[Aircraft]]\n        FlowSolverBCType = BCWallViscousIsothermal\n        [[[input data specification]]]\n            type = scalar\n            temperature = 280 [K]\n")
            f.write("    [[Farfield]]\n        FlowSolverBCType = BCFarfield\n")
            for i in range(nSym):
                f.write(f"    [[Symmetry{i}]]\n        FlowSolverBCType = BCSymmetryPlane\n")

        print(f"BC block generated for {nSym} symmetries in bc_markers.txt")

#===============================================================================================
# IN: BAR1: profile with ice, correctly meshed, in XY plane, simple loop, no branching
# IN: ht: total height of boundary layer mesh
# IN: hf: height of first wall cell
# OUT: meshBL
#===============================================================================================
def mesherBL(smoothed1D, ht, hf, extruder):

    BAR1 = C.convertArray2Tetra(smoothed1D)
    profile = C.initVars(BAR1, 'CoordinateZ', 0.)
    profile = C.convertBAR2Struct(profile)
    if checkOrder(profile) == 0:
        T._reorder(profile, (-1, 2, 3))

    nLayersRatio = 1.15
    nLayers = int(np.log(ht / hf) / np.log(nLayersRatio)) + 1
    nLayers = max(nLayers, 10)
    nLayers = min(nLayers, 100)

    if extruder == 0:
        niter = 100
        distrib = D.line((0, 0, 0), (ht, 0, 0), N=nLayers)
        goOn = True
        while goOn:
            bl = G.addNormalLayers(profile, distrib, niter=niter, check=True)
            nk = Internal.getZoneDim(bl)[2]
            if nk < nLayers:
                print('addNL was stopped %d out of %d'%(nk, nLayers))
                niter += 200
                if niter > 1000:
                    goOn = False
            else:
                goOn = False
        T._reorder(bl, (-1, 2, 3))
        print('addNL has reached %d out of %d'%(nk, nLayers))
    else:
        zoneDim = Internal.getZoneDim(profile)
        distrib = G.cart((0, 0, 0), (1, ht / (nLayers - 1), 1),
                         (zoneDim[1], nLayers, 1))
        profile2 = T.reorder(profile, (-1, 2, 3))
        bl = G.hyper2D(profile2, distrib, "O")

    C.convertPyTree2File(bl, 'monitoring/bl.cgns')

    d = G.cartr2((0, 0, 0), (hf, 1, 1), (1.1, 1.1, 1.1), (ht, 1, 1))
    d = T.subzone(d, (1, 1, 1), (-1, 1, 1))
    C._initVars(d, '{CoordinateX}={CoordinateX}/%g' % ht)
    bl = G.map(bl, d, dir=2)
    C.convertPyTree2File(bl, 'monitoring/bl2.cgns')
    bl = C.convertArray2Hexa(bl)
    bl = G.close(bl)

    return BAR1, bl, nk, nLayers


#===============================================================================================
# IN: BAR1: profile with ice, correctly meshed, in XY plane, simple loop, no branching
# IN: BAR3: outer domain boundary, correctly meshed, in XY plane, simple curve with no branching
# IN: bl: boundary layer mesh generated in mesherBL
# IN: dz: fictionnal thickness for 3D simulation
# OUT: mesh: ready for coda
#===============================================================================================
def mesherT3(BAR1, BAR3, bl, dz):

    BAR3 = C.convertBAR2Struct(BAR3)
    BAR3 = C.initVars(BAR3, 'CoordinateZ', 0.)
    if checkOrder(BAR3) == 0:
        T._reorder(BAR3, (-1, 2, 3))
    BAR3 = C.convertArray2Tetra(BAR3)
    C.convertPyTree2File(BAR3, 'monitoring/exterior.cgns')

    #=========================
    # exterior of bl mesh - OK
    #=========================
    ext = P.exteriorFaces(bl)
    ext = T.splitConnexity(ext)
    if len(ext) < 2:
        raise RuntimeError("BL exterior has < 2 components — geometry too rough")
    # Select outermost component
    ext = ext[1] # always true
    T._reorder(ext, (-1,))

    meshValid = checkBL(ext)

    borders = T.join(ext, BAR3)
    m2 = G.T3mesher2D(borders, triangulateOnly=0, grading=1.05, metricInterpType=0)
    T._reorder(m2, (1,))
    C.convertPyTree2File(m2, 'monitoring/T3mesh.cgns')

    # Merge and BCs (simplified — call full mesher() for production)
    m = [bl, m2]
    T._addkplane(m)
    T._contract(m, (0, 0, 0), (1, 0, 0), (0, 1, 0), dz)
    merged = C.mergeConnectivity(m[0], m[1])

    #==========
    # BCs - OK
    #==========
    e0 = P.exteriorFaces(m)
    e0 = T.breakElements(e0)

    # recover BAR1
    e1 = T.addkplane(BAR1)
    T._contract(e1, (0,0,0), (1,0,0), (0,1,0), dz)
    e1[0] = 'wall'
    C._addBC2Zone(merged, 'wall', 'BCWall', subzone=e1)

    # recover BAR3
    e3 = T.addkplane(BAR3)
    T._contract(e3, (0,0,0), (1,0,0), (0,1,0), dz)
    e3[0] = 'far'
    C._addBC2Zone(merged, 'far', 'BCFarfield', subzone=e3)

    # recover symmetry planes
    es1 = P.selectCells(e0, '{CoordinateZ}<%g'%(0.5*dz), strict=True)
    es2 = P.selectCells(e0, '{CoordinateZ}>%g'%(dz-0.5*dz), strict=True)

    c = 0
    for e in es1+es2:
        if C.getNCells(e) > 0:
            e[0] = 'sym%d'%c
            C._addBC2Zone(merged, 'sym%d'%c, 'BCSymmetryPlane', subzone=e); c += 1

    # Reorientation axes (Y↔Z swap)
    for z in Internal.getZones(merged):
        yp = Internal.getNodeFromName(z, 'CoordinateY')
        zp = Internal.getNodeFromName(z, 'CoordinateZ')
        temp = zp[1]
        zp[1] = yp[1]
        yp[1] = -temp

    return bl, m2, merged

#===============================================================================================
# Normalizes the point density of a profile (coarsens or refines)
# Protects the trailing edge (TE) from excessive coarsening to preserve geometric accuracy.
# Calculates an adaptive splitting sensitivity angle based on the final mesh spacing.

# IN: profile: 1D contour profile to be normalized
# IN: chord length (used for density calculation)
# IN: targetDensity: number of points per unit chord to check if coarsens or refines
# IN: factor: used in T.oneovern if the profile is too dense
# IN: te_protect_frac: fraction of the chord at the TE protected from coarsening (default: 0.02)
# OUT: profile: the normalized (coarsened or refined) profile
# OUT: splitSensib: adaptive splitting angle threshold (in degrees, bounded between 5 and 20)
#===============================================================================================
def normaliseDensity(profile, chord, targetDensity, factor=None, te_protect_frac=0.02):

    
    N = getNpts(profile)
    density = N / chord
    ratio = density / targetDensity

    if factor == None:
        factor = ratio 
    else :
        factor = factor + 1

    print(f"  [normaliseDensity] N={N}  chord={chord:.4f}"
          f"  density={density:.1f} pts/chord  target={targetDensity:.0f}")

    type = 'none'
    aBefore = checkArea(profile)

    if ratio > 1.2:
        # Identifier les indices proches du TE (x > chord*(1 - te_protect_frac))
        p = C.convertBAR2Struct(profile)
        x, y, _ = getCoords(p)
        te_x_threshold = chord * (1.0 - te_protect_frac)
        
        te_indices = np.where(x >= te_x_threshold)[0]
        print(f"    → {len(te_indices)} points protégés au TE (x >= {te_x_threshold:.4f})")

        # Coarsening global via oneovern
        # print(f"    → Coarsening {N}/{factor} via T.oneovern")
        # p_coarse = T.oneovern(p, (factor, 1, 1)) premiere stratégie

        # Coarsening global via refine
        targetN = int(targetDensity * chord)    
        print(f"    → Coarsening {N}→{targetN} via refine")
        p_coarse = D.refine(p, N=targetN)
        
        # Réinjecter les points TE supprimés si nécessaire
        xc, yc, _ = getCoords(p_coarse)
        te_x_coarse = np.where(xc >= te_x_threshold)[0]
        
        if len(te_x_coarse) < 2:
            # Le TE a été trop coarsen — on reconstruit avec les points originaux
            print(f"    → TE sous-résolu après coarsening, réinjection des points TE")
            
            # Points non-TE du profil coarsen
            non_te_mask = xc < te_x_threshold
            x_new = np.concatenate([xc[non_te_mask], x[te_indices]])
            y_new = np.concatenate([yc[non_te_mask], y[te_indices]])
            
            # Reconstruction de la zone
            npts_new = len(x_new)
            rebuilt = D.line(
                (float(x_new[0]),  float(y_new[0]),  0.),
                (float(x_new[-1]), float(y_new[-1]), 0.),
                N=npts_new
            )
            Internal.getNodeFromName(rebuilt, 'CoordinateX')[1].flat[:] = x_new
            Internal.getNodeFromName(rebuilt, 'CoordinateY')[1].flat[:] = y_new
            Internal.getNodeFromName(rebuilt, 'CoordinateZ')[1].flat[:] = np.zeros(npts_new)
            profile = C.convertArray2Tetra(rebuilt)
        else:
            profile = C.convertArray2Tetra(p_coarse)

        type = 'coarser'
       
    elif ratio < 0.8:
        nTarget = max(10, int(round(targetDensity * chord)))
        print(f"    → Under-dense (×{ratio:.1f}): refine to {nTarget} pts")
        profile = D.refine(profile, N=nTarget)
        type = 'refinement'
    else:
        print(f"    → Density OK : no change")
    
    aAfter = checkArea(profile)
    if aBefore > 1e-12:
        areaError = (aAfter - aBefore) / aBefore * 100
        print(f"\n[Conservation] Area error due to {type}: {areaError:.8f} % ")
        
    return profile

# Suppose qu'on rafine tjrs. 
def getDensity(profile, chord, targetDensity, factor, te_protect_frac=0.02):
    
    factor = factor + 1
    N = getNpts(profile)
    density = N / chord

    print(f"  [getDensity] N={N}  chord={chord:.4f}"
          f"  density={density:.1f} pts/chord ")

    aBefore = checkArea(profile)

    # Identifier les indices proches du TE (x > chord*(1 - te_protect_frac))
    p = C.convertBAR2Struct(profile)
    x, y, _ = getCoords(p)
    te_x_threshold = chord * (1.0 - te_protect_frac)
    
    te_indices = np.where(x >= te_x_threshold)[0]
    print(f"    → {len(te_indices)} points protégés au TE (x >= {te_x_threshold:.4f})")

    # Coarsening global via oneovern
    print(f"    → Coarsening ×{factor} via T.oneovern")
    p_coarse = T.oneovern(p, (factor, 1, 1))
    
    # Réinjecter les points TE supprimés si nécessaire
    xc, yc, _ = getCoords(p_coarse)
    te_x_coarse = np.where(xc >= te_x_threshold)[0]
    
    if len(te_x_coarse) < 2:
        # Le TE a été trop coarsen — on reconstruit avec les points originaux
        print(f"    → TE sous-résolu après coarsening, réinjection des points TE")
        
        # Points non-TE du profil coarsen
        non_te_mask = xc < te_x_threshold
        x_new = np.concatenate([xc[non_te_mask], x[te_indices]])
        y_new = np.concatenate([yc[non_te_mask], y[te_indices]])
        
        # Reconstruction de la zone
        npts_new = len(x_new)
        rebuilt = D.line(
            (float(x_new[0]),  float(y_new[0]),  0.),
            (float(x_new[-1]), float(y_new[-1]), 0.),
            N=npts_new
        )
        Internal.getNodeFromName(rebuilt, 'CoordinateX')[1].flat[:] = x_new
        Internal.getNodeFromName(rebuilt, 'CoordinateY')[1].flat[:] = y_new
        Internal.getNodeFromName(rebuilt, 'CoordinateZ')[1].flat[:] = np.zeros(npts_new)
        profile = C.convertArray2Tetra(rebuilt)
    else:
        profile = C.convertArray2Tetra(p_coarse)
       
    N_final  = getNpts(profile)
    ds_mean  = chord / N_final
    # Loi empirique : entre 5° (très fin) et 20° (très grossier)
    splitSensib = float(np.clip(10.0 * np.sqrt(ds_mean / chord) * np.sqrt(targetDensity), 5.0, 20.0))

    print(f"    → N_final={N_final}  ds_mean={ds_mean:.4e}"
          f"  splitSensib adaptatif={splitSensib:.2f}°")
    
    aAfter = checkArea(profile)
    if aBefore > 1e-12:
        areaError = (aAfter - aBefore) / aBefore * 100
        print(f"\n[Conservation] Area error due to coarser: {areaError:.8f} % ")

    return profile, splitSensib

#==========================================================
# Merge segments shorter than factor * hf using G.close
#==========================================================
def repairShortSegments(curve, hf, factor=5.0):

    try:
        p = C.convertBAR2Struct(curve)
        x, y, z = getCoords(p)
        x, y, z = x.copy(), y.copy(), z.copy()

        ds = np.sqrt(np.diff(x)**2 + np.diff(y)**2)
        tol_phys = factor * hf
        tol_stat = 1e-2 * float(np.mean(ds))
        tol = max(tol_phys, tol_stat)
        print(f"tol = {tol}")

        n_before = len(x)
        keep = np.ones(len(x), dtype=bool)

        i = 0
        while i < len(x) - 1:
            if ds[i] < tol:
                # Supprimer le point i+1 (point d'arrivée du segment court)
                # sauf si c'est le premier ou dernier point (endpoints fixes)
                if i + 1 < len(x) - 1:
                    keep[i + 1] = False
                    # Recalculer ds pour le segment suivant
                    ds_new = np.sqrt((x[i] - x[i+2])**2 + (y[i] - y[i+2])**2)
                    # Mise à jour locale pour éviter de sauter un autre court segment
                    if i + 2 < len(ds) + 1:
                        ds = np.concatenate([ds[:i],[ds_new],ds[i+2:]])
                        x = x[keep[:len(x)]]  # pas encore — on filtre à la fin
                        y = y[keep[:len(y)]]
                        z = z[keep[:len(z)]]
                        keep = np.ones(len(x), dtype=bool)
                        ds = np.sqrt(np.diff(x)**2 + np.diff(y)**2)
                        # Repart du même i pour vérifier le nouveau segment
                        continue
            i += 1

        n_after = len(x)
        removed = n_before - n_after
        print(f"  [repairShortSegments] {removed} nœud(s) supprimé(s) "
              f"(tol={tol:.2e} m) : {n_before} → {n_after} pts")

        if removed == 0:
            return curve

        # Reconstruction de la zone structurée
        import Converter.Internal as Internal
        import Geom.PyTree as D
        new_curve = D.line((float(x[0]), float(y[0]), 0.),
                           (float(x[-1]), float(y[-1]), 0.), N=len(x))
        Internal.getNodeFromName(new_curve, 'CoordinateX')[1].flat[:] = x
        Internal.getNodeFromName(new_curve, 'CoordinateY')[1].flat[:] = y
        Internal.getNodeFromName(new_curve, 'CoordinateZ')[1].flat[:] = z
        return new_curve

    except Exception as e:
        print(f"  [repairShortSegments] failed ({e}) — curve unchanged")
        return curve
    
#===========================================================================
# Classifies a sub-curve based on its junction angles and internal rugosity
#===========================================================================
def getZoneType(rug, angleLeft, signLeft, angleRight, signRight, roughThreshold, protectAngle):
    
    if angleLeft <= angleRight:
        sharpAngle, sharpSign = angleLeft, signLeft
    else:
        sharpAngle, sharpSign = angleRight, signRight

    isSharp = sharpAngle < protectAngle
    
    if isSharp and sharpSign == -1: return 'p'
    if isSharp and sharpSign == +1: return 'v'
    if rug > roughThreshold: return 'r'
    
    return 'n'

def scoreZone(rug, alphaInternal, angleLeft, signLeft, angleRight, signRight,
              length, chord, roughThreshold, protectAngle):
    """
    Retourne un dict de scores [0-1] pour chaque type, et le type dominant.
    """
    # --- Score PIC (convexe) ---
    # Jonction la plus critique
    if angleLeft <= angleRight:
        sharpAngle, sharpSign = angleLeft, signLeft
    else:
        sharpAngle, sharpSign = angleRight, signRight

    score_p, score_v = 0.0, 0.0

    if sharpSign == -1:  # convexe
        # Score entre 0 et 1 selon la sévérité de l'angle
        score_p = float(np.clip((protectAngle - sharpAngle) / protectAngle, 0, 1))
    elif sharpSign == +1:  # concave
        score_v = float(np.clip((protectAngle - sharpAngle) / protectAngle, 0, 1))

    # --- courbure interne (alphaInternal = tableau des angles intérieurs) ---
    # Une zone pic a ses angles INTERNES aussi tournés vers l'extérieur
    if len(alphaInternal) > 2:
        meanAlpha = float(np.mean(alphaInternal[1:-1]))
        # 180° = plat, <180° = courbé
        curvatureScore = float(np.clip((180.0 - meanAlpha) / 60.0, 0, 1))
        # Renforce le score pic ou valley selon le signe moyen
        signs = np.sign(180.0 - alphaInternal[1:-1])
        dominantSign = float(np.mean(signs))
        if dominantSign < 0:   # majorité convexe
            score_p = max(score_p, 0.4 * curvatureScore)
        elif dominantSign > 0: # majorité concave
            score_v = max(score_v, 0.4 * curvatureScore)

    # --- Score RUGUEUX ---
    score_r = float(np.clip(rug / (roughThreshold * 3.0), 0, 1))
    consistency = curvatureConsistency(alphaInternal)
    # Zone chaotique → renforce rugosité
    score_r = max(score_r, (1.0 - consistency) * 0.8)
    # Zone uniforme → renforce pic ou creux (× 1.2, borné à 1)
    if consistency > 0.8:
        score_p = min(score_p * 1.2, 1.0)
        score_v = min(score_v * 1.2, 1.0)

    # --- Score NEUTRE (résiduel) ---
    score_n = 1.0 - max(score_p, score_v, score_r)
    score_n = max(score_n, 0.0)

    scores = {'p': score_p, 'v': score_v, 'r': score_r, 'n': score_n}

    # Type dominant (avec seuil minimal pour éviter le bruit)
    if score_p > 0 and score_p >= score_v:
        zType = 'p'
    elif score_v > 0 and score_v > score_p:
        zType = 'v'
    elif score_r > 0.3:   # seuil minimal pour éviter le bruit sur 'r'
        zType = 'r'
    else:
        zType = 'n'

    return zType, scores

def splitting (profile, splitSensib, roughThreshold, chord, maxSegFrac=0.05):
    """
    Passe 1 : coarse split pour trouver les grandes transitions (pics, TE).
    Passe 2 : fine split uniquement sur les zones rugueuses.
    """
    # Passe 1 : angle large → grandes zones homogènes
    subCoarse = T.splitCurvatureAngle(profile, splitSensib * 2.0)

    # Passe 2 : sur les zones rugueuses seulement, re-split fin
    subFinal = []
    for sub in subCoarse:
        alphas = getSharpestAngles(sub)
        rug = rugosity(alphas)
        if rug > roughThreshold:
            # Re-split local plus fin
            subLocal = T.splitCurvatureAngle(sub, splitSensib * 0.5)
            subFinal.extend(subLocal)
        else:
            subFinal.append(sub)

    # Passe 3 : split forcé sur les segments trop longs
    if chord is not None and maxSegFrac > 0:
        maxLen = maxSegFrac * chord
        subFinal2 = []
        for sub in subFinal:
            L = D.getLength(sub)
            if L > maxLen and getNpts(sub) >= 6:
                # Nombre de morceaux nécessaires
                nParts = int(np.ceil(L / maxLen))
                nPts = getNpts(sub)
                # Découpe régulière en nParts morceaux par sous-zonage
                pts_per_part = max(3, nPts // nParts)
                x, y, z = getCoords(sub)
                start = 0
                parts = []
                for k in range(nParts):
                    end = start + pts_per_part if k < nParts - 1 else nPts
                    end = min(end, nPts)
                    if end - start < 2:
                        break
                    part = D.line(
                        (float(x[start]), float(y[start]), 0.),
                        (float(x[end-1]), float(y[end-1]), 0.),
                        N=end - start
                    )
                    Internal.getNodeFromName(part, 'CoordinateX')[1].flat[:] = x[start:end]
                    Internal.getNodeFromName(part, 'CoordinateY')[1].flat[:] = y[start:end]
                    Internal.getNodeFromName(part, 'CoordinateZ')[1].flat[:] = z[start:end]
                    parts.append(part)
                    start = end - 1  # overlap d'un point pour la continuité
                if parts:
                    print(f"  [splitting] Segment long ({L:.4f} > {maxLen:.4f}) → {len(parts)} morceaux")
                    subFinal2.extend(parts)
                else:
                    subFinal2.append(sub)
            else:
                subFinal2.append(sub)
        subFinal = subFinal2

    print(f"  [splitting] {len(subCoarse)} (coarse) → {len(subFinal)} zones (fine + longueur)")
    return subFinal

def curvatureConsistency(alphas):
    """
    Retourne un score [0-1] : 1 = courbure uniforme (même signe),
    0 = courbure chaotique (alternance de signes = rugosité géométrique).
    """
    if len(alphas) < 4:
        return 1.0
    deviations = 180.0 - alphas[1:-1]
    signs = np.sign(deviations)
    # Nombre de changements de signe
    nFlips = int(np.sum(np.abs(np.diff(signs)) > 0))
    maxFlips = len(signs) - 1
    if maxFlips == 0:
        return 1.0
    return float(1.0 - nFlips / maxFlips)

def computeSplitSensib(profile, targetDensity, chord, percentile=15.0, minVal=3.0, maxVal=20.0):
    """
    Calcule splitSensib à partir de la distribution réelle des angles du profil.
    On cible le percentile bas des déviations angulaires — ainsi les zones
    à courbure modérée mais non nulle seront quand même coupées.
    """
    alphas = getSharpestAngles(profile)
    deviations = 180.0 - alphas  # 0 = plat, grand = angle vif

    # On ignore les noeuds parfaitement plats (déviation < 0.01°)
    significant = deviations[deviations > 0.01]

    if len(significant) < 5:
        # Profil trop lisse ou trop peu de points : fallback empirique
        ds_mean = chord / getNpts(profile)
        return float(np.clip(10.0 * np.sqrt(ds_mean / chord) * np.sqrt(targetDensity), minVal, maxVal))

    # splitSensib = percentile bas des déviations significatives
    # → on coupe dès qu'on dépasse ce seuil, ce qui capture les zones à courbure modérée
    raw = float(np.percentile(significant, percentile))
    result = float(np.clip(raw, minVal, maxVal))

    print(f"  [computeSplitSensib] déviation médiane={np.median(significant):.2f}°"
          f"  p{percentile:.0f}={raw:.2f}°  → splitSensib={result:.2f}°")
    return result
#===========================================================================
# Detects topology
#===========================================================================
def detectTopology (BAR1, profileBody, distThreshold, splitSensib,
                   roughThreshold, protectAngle, chord=1.0, borrowFrac=0.0):
    
    _SIGN = {+1: "PIC▲", -1: "creux▼", 0: "flat"}
    profile = C.initVars(BAR1, 'CoordinateZ', 0.)
    profile = C.convertBAR2Struct(profile)
    if isClockwise(profile):
        T._reorder(profile, (-1, 2, 3))
    
    # subCurves = T.splitCurvatureAngle(profile, splitSensib)
    subCurves = splitting (profile, splitSensib, roughThreshold, chord=chord, maxSegFrac=0.05)

    # ====================================================================
    # PURGE DES MICRO-COURBES (< 3 points)
    # ====================================================================
    idx = 0
    while idx < len(subCurves):
        if getNpts(subCurves[idx]) < 3 and len(subCurves) > 1:
            if idx < len(subCurves) - 1:
                subCurves[idx+1] = T.join([subCurves[idx], subCurves[idx+1]])
                subCurves.pop(idx)
            else:
                subCurves[idx-1] = T.join([subCurves[idx-1], subCurves[idx]])
                subCurves.pop(idx)
        else:
            idx += 1

    n = len(subCurves)
    print(f"\n[remesh]  {n} sub-curve(s)")

    # =====================================
    # Definition of distThreshold and chord
    # =====================================
    if profileBody is not None:
        if isClockwise(profileBody):
            T._reorder(profileBody, (-1, 2, 3))
        chord = getChord(profileBody) 
        thr = distThreshold if distThreshold is not None else 0.001 * chord
        print(f"  Mode DISTANCE  corde={chord:.4f}  seuil={thr:.4f}")
        subCurves = D2W.distance2Walls(subCurves, [profileBody], type='ortho', loc='nodes', signed=0)
    else:
        xN, _, _ = getCoords(profile)
        chord = getChord(profile) 

    # ================================
    # Junction calculation
    # ================================
    isOpen = closeGeometryRelative(profile, toClose=False)
    junctionData = [(180.0, 0)] * (n + 1)
    for i in range(1, n):
        junctionData[i] = junctionAngleSigned(subCurves[i-1], subCurves[i])
    if n > 1 and not isOpen:
        jd = junctionAngleSigned(subCurves[-1], subCurves[0])
        junctionData[0] = junctionData[n] = jd

    # ================================
    # Pre-compute protection
    # ================================
    isProtList = []
    pinfoList = []
    alphaCache = [getSharpestAngles(sub) for sub in subCurves]
    rugCache = [rugosity(alpha) for alpha in alphaCache]

    for i, sub in enumerate(subCurves):
        p, inf = isProtected(sub, rugCache[i], alphaCache[i], profileBody=profileBody,
                             distThreshold=distThreshold,
                             roughThreshold=roughThreshold,
                             protectAngle=protectAngle,
                             maxCornerLength=0.0)

        isProtList.append(p)
        pinfoList.append(inf)

    # ================================
    # Cleaning of the protective islands to create 2 zones (protected and unprotected)
    # ================================
    nSub = len(isProtList)

    # If there is a mix of protected and unprotected segments
    if nSub > 0 and (False in isProtList) and (True in isProtList):
        # Start scanning from a free segment to avoid splitting a continuous block
        startIdx = isProtList.index(False)
        # Contiguous blocks of protected sub-curves (True)
        blocks = []
        currentBlock = []
        for offset in range(nSub):
            idx = (startIdx + offset) % nSub
            if isProtList[idx]:
                currentBlock.append(idx)
            else:
                if len(currentBlock) > 0:
                    blocks.append(currentBlock)
                    currentBlock = []
        if len(currentBlock) > 0:
            blocks.append(currentBlock)

        # If more than one protected block exists, we clean
        if len(blocks) > 1:
            
            # trailing edge is protected in priority and protects adjacents, fall back to the heaviest block.
            te_idx = max( range(nSub), key=lambda i: np.mean((getCoords(subCurves[i])[0])))
            print(f"  [Topology] 📍 Trailing edge auto-detected at sub-curve [{te_idx+1:02d}]")

            bestBlockIdx = next(
                (i for i, block in enumerate(blocks) if te_idx in block),
                None
            )

            if bestBlockIdx is None:
                # Fallback: keep the heaviest block (legacy behaviour)
                score = [sum(getNpts(subCurves[i]) for i in block) for block in blocks]
                bestBlockIdx = int(np.argmax(score))
                print("  [Topology] ⚠️  Trailing edge sub-curve not in any protected block — "
                    "falling back to largest block.")
            # ── END AUTO-DETECT ────────────────────────────────────────────────────

            # Unprotect all other blocks (the "false pearls" / ice chunks)
            for i, block in enumerate(blocks):
                if i != bestBlockIdx:
                    for idx in block:
                        isProtList[idx] = False
                        print(f"  [Topology] 🧹 Sub-curve [{idx+1:02d}] unprotected "
                            f"(merged into ice)")
                        
    # ====================================================================
    # PRE-GROUPING: Merge adjacent zones of the SAME TYPE
    # ====================================================================
    # Identification
    zoneTypes = []
    scoresMap = {}
    for i in range(n):
        if not isProtList[i]:
            aL, sL = junctionData[i]
            aR, sR = junctionData[(i + 1) % n]
            zType, scores = scoreZone(
                rugCache[i], alphaCache[i],
                aL, sL, aR, sR,
                D.getLength(subCurves[i]), chord,
                roughThreshold, protectAngle
            )
            zoneTypes.append(zType)
            scoresMap[i] = scores
            print (f"[{i}] {zType} : p={scores.get('p',0):.2f} v={scores.get('v',0):.2f} r={scores.get('r',0):.2f}")
        else:
            zoneTypes.append('PROT')
            scoresMap[i] = {}

    #Regrouping
    startIdx = 0
    if 'PROT' in zoneTypes:
        startIdx = zoneTypes.index('PROT')

    blocks = []
    currentBlock = [startIdx]

    for offset in range(1, nSub):
        idx = (startIdx + offset) % nSub
        
        # If same type, add to current block
        if zoneTypes[idx] == zoneTypes[currentBlock[-1]]:
            currentBlock.append(idx)
        else:
            blocks.append(currentBlock)
            currentBlock = [idx]
            
    blocks.append(currentBlock)

    mergedCurves = []
    mergedProt = []
    mergedTypes = []
    mergedScores = {}
    for idx_block, block in enumerate(blocks):
        if len(block) == 1:
            mergedCurves.append(subCurves[block[0]])
        else:
            curvesToJoin = [subCurves[i] for i in block]
            mergedCurves.append(T.join(curvesToJoin))
            print(f"  [Grouping] 🔗 Merged {len(block)} adjacent '{zoneTypes[block[0]]}' zones into one macro-block.")
            
        mergedProt.append(isProtList[block[0]])
        mergedTypes.append(zoneTypes[block[0]])

        mergedScores[idx_block] = scoresMap[block[0]]
    
    # Replace old lists with the new merged ones
    subCurves = mergedCurves
    isProtList = mergedProt
    macroTypes = mergedTypes
    scoresMap = mergedScores
    n = len(subCurves)

    for i in range(n):
        print (f"[{i}] {macroTypes[i]} : p={scoresMap[i].get('p',0):.2f} v={scoresMap[i].get('v',0):.2f} r={scoresMap[i].get('r',0):.2f}")
        

    alphaCache = [getSharpestAngles(sub) for sub in subCurves]
    rugCache   = [rugosity(alpha) for alpha in alphaCache]

    # Recompute junction data for the new macro-curves
    junctionData = [(180.0, 0)] * (n + 1)
    for i in range(1, n):
        junctionData[i] = junctionAngleSigned(subCurves[i-1], subCurves[i])
    if n > 1 and not isOpen:
        jd = junctionAngleSigned(subCurves[-1], subCurves[0])
        junctionData[0] = junctionData[n] = jd

    prevSubCurves = list(subCurves)
    for i, sub in enumerate(prevSubCurves):

        N = getNpts(sub)
        L = D.getLength(sub)
        aL, sL = junctionData[i]
        aR, sR = junctionData[i + 1]

        protected = isProtList[i]
        pinfo = pinfoList[i]
        xA, yA, zA = getCoords(sub)

        #  Point borrowing at protected/free boundaries 
        if not protected:
            targetDist = L * borrowFrac
            if i > 0 and isProtList[i - 1]:
                xP, yP, zP = getCoords(prevSubCurves[i - 1])
                k = min(getKFromDist(xP, yP, targetDist, fromEnd=True), len(xP) - 3)
                if k > 0:
                    xA = np.concatenate([xP[-k-1:-1], xA])
                    yA = np.concatenate([yP[-k-1:-1], yA])
                    zA = np.concatenate([zP[-k-1:-1], zA])
            if i < n - 1 and isProtList[i + 1]:
                xN2, yN2, zN2 = getCoords(prevSubCurves[i + 1])
                k = min(getKFromDist(xN2, yN2, targetDist, fromEnd=False), len(xN2) - 3)
                if k > 0:
                    xA = np.concatenate([xA, xN2[1:k+1]])
                    yA = np.concatenate([yA, yN2[1:k+1]])
                    zA = np.concatenate([zA, zN2[1:k+1]])
        else:
            if i > 0 and not isProtList[i - 1]:
                prevL = D.getLength(prevSubCurves[i - 1])
                k = min(getKFromDist(xA, yA, prevL * borrowFrac, fromEnd=False), len(xA) - 3)
                if k > 0:
                    xA, yA, zA = xA[k:], yA[k:], zA[k:]
            if i < n - 1 and not isProtList[i + 1]:
                nextL = D.getLength(prevSubCurves[i + 1])
                k = min(getKFromDist(xA, yA, nextL * borrowFrac, fromEnd=True), len(xA) - 3)
                if k > 0:
                    xA, yA, zA = xA[:-k], yA[:-k], zA[:-k]

        if len(xA) < 2:
            xA, yA, zA = getCoords(sub)   # revert

        seg = D.line((float(xA[0]), float(yA[0]), 0.),
                     (float(xA[-1]), float(yA[-1]), 0.), N=len(xA))
        Internal.getNodeFromName(seg, 'CoordinateX')[1].flat[:] = xA
        Internal.getNodeFromName(seg, 'CoordinateY')[1].flat[:] = yA
        Internal.getNodeFromName(seg, 'CoordinateZ')[1].flat[:] = zA

        subCurves[i] = seg
      
    # Si besoin de print pour voir ce qu'il se passe Rafine protege ? 
    # for i, sub in enumerate(subCurves):

    #     npts = getNpts(sub)
    #     xA2, yA2, zA2 = getCoords(sub)
    #     pStart = (xA2[0], yA2[0], zA2[0])
    #     pEnd = (xA2[-1], yA2[-1], zA2[-1])

    #     protected = isProtList[i]
    #     if protected:
    #         L = D.getLength (sub) 
    #         n = max(3, int(L * 200))
    #         if not (L< 1e-10 or npts < 3) : 
    #             sub = refineProfile (sub, n, i)
    #         subCurves[i] = forceEndpoints(sub, pStart, pEnd)

    # n = len(subCurves)
    # scoresMap = {}
    # for i in range(n):  # n = len(subCurves) après merge
    #     if not isProtList[i]:
    #         aL, sL = junctionData[i]
    #         aR, sR = junctionData[(i + 1) % n]
    #         zType, scores = scoreZone(
    #             rugCache[i], alphaCache[i],
    #             aL, sL, aR, sR,
    #             D.getLength(subCurves[i]), chord,
    #             roughThreshold, protectAngle
    #         )
    #         print (f"[{i}] {zType} : p={scores.get('p',0):.2f} v={scores.get('v',0):.2f} r={scores.get('r',0):.2f}")
    #         scoresMap[i] = scores
    #         macroTypes[i] = zType
    #     else:
    #         scoresMap[i] = {}

        
    return subCurves, isProtList, macroTypes, junctionData, profile, scoresMap
    
def refineProfile (curve, npts = 300, i = 0):

    C.convertPyTree2File(curve, f'monitoring/bef{i}.plt')
    # x, y, _ = getCoords(curve)
    # h1 = np.sqrt((x[1]-x[0])**2 + (y[1]-y[0])**2)
    # h2 = np.sqrt((x[-1]-x[-2])**2 + (y[-1]-y[-2])**2)   
    # d = D.distrib2(curve, h1, h2)

    h = D.getLength(curve)/(npts-1)
    d = D.distrib1(curve, h)

    C.convertPyTree2File(d, f'monitoring/line{i}.plt')
    curve = G.map(curve, d)
    C.convertPyTree2File(curve, f'monitoring/{i}.plt')
   
    return curve

#============================================================================
# IN: BAR1: 1D struct or BAR (input geometry)
# IN: splitSensib: [°] cuts if angle < (180 - splitSensib) 
# IN: peakFactor: refine ration for rough zones
# IN: valleyFactor: coarse ratio for smoothed or valley zones
# IN: roughThreshold: rough zone to smooth
# IN: protectAngle: [°] corner below = protected singularity
# OUT: processedFree : list of free 1D-STRUCT zones (remeshed, not smoothed)
# OUT: processedProtected : list of protected 1D-STRUCT zones (untouched)
# OUT: shapeLabels: list of str, same order as processed
# OUT: orderMap : list of ('free'|'protected', index_in_its_list)
# OUT: to reconstruct the original order for T.join
# OUT: profile : the structured profile used for splitCurvatureAngle
# OUT: chord : chord length
#============================================================================
def mapper(subCurves, isProtList, macroTypes,
                     peakFactor=1.5, valleyFactor=0.5,
                     roughFactor=1.4):
    
    processedFree, processedProtected, shapeLabels = [], [], []
    orderMap = []
    for i, sub in enumerate(subCurves):

        xA2, yA2, zA2 = getCoords(sub)
        pStart = (xA2[0], yA2[0], zA2[0])
        pEnd = (xA2[-1], yA2[-1], zA2[-1])

        npts = getNpts(sub)

        protected = isProtList[i]
        if not protected:
            currentType = macroTypes[i] 
            sub, shape = adaptiveRemesh(sub, zType=currentType,
                             peakFactor=peakFactor, valleyFactor=valleyFactor,
                             roughFactor=roughFactor)
                        
            sub = forceEndpoints(sub, pStart, pEnd)
            processedFree.append(Internal.copyTree(sub))
            orderMap.append(('free', len(processedFree) - 1))
            shapeLabels.append(shape)
        else:
            # L = D.getLength (sub) 
            # n = max(3, int(L * 200))
            # if not (L< 1e-10 or npts < 2) : 
            #     sub = refineProfile (sub, n, i)
            # sub = forceEndpoints(sub, pStart, pEnd)
            # sub = T.oneovern(sub, (2, 1, 1))
            # sub = forceEndpoints(sub, pStart, pEnd)

            processedProtected.append(Internal.copyTree(sub))
            orderMap.append(('protected', len(processedProtected) - 1))
            shapeLabels.append('protected')
        
        print(f"  [{i+1:02d}] {'PROT' if protected else 'FREE'} ({npts} → {getNpts(sub)} pts)")
        
    return processedFree, processedProtected, shapeLabels, orderMap

#===========================================================================
# Reconstruct the full ordered profile from the two separate lists.
# Uses orderMap to interleave free and protected zones in correct order.
# Repairs endpoints between consecutive zones.
#===========================================================================
def joinInOrder(processedFree, processedProtected, orderMap, mergeTol=1e-2):

    ordered = []
    for kind, idx in orderMap:
        if kind == 'free':
            ordered.append(Internal.copyTree(processedFree[idx]))
        else:
            ordered.append(Internal.copyTree(processedProtected[idx]))

    # =========================================================================
    # Repair shared endpoints (exact averaging)
    # =========================================================================
    n = len(ordered)
    for i in range(n):
        next_i = (i + 1) % n
        kind_i    = orderMap[i][0]
        kind_next = orderMap[next_i][0]

        xA, yA, zA = getCoords(ordered[i])
        xB, yB, zB = getCoords(ordered[next_i])

        gap = np.sqrt((float(xA[-1]) - float(xB[0]))**2 +
                    (float(yA[-1]) - float(yB[0]))**2)

        if gap > 1e-12:
            print ("il y a une rupture")
            # Le point de référence = le côté protégé s'il existe
            if kind_i == 'protected':
                # La zone i est protégée → on force la zone suivante à rejoindre i
                rx, ry, rz = float(xA[-1]), float(yA[-1]), float(zA[-1])
                Internal.getNodeFromName(ordered[next_i], 'CoordinateX')[1].flat[0] = rx
                Internal.getNodeFromName(ordered[next_i], 'CoordinateY')[1].flat[0] = ry
                Internal.getNodeFromName(ordered[next_i], 'CoordinateZ')[1].flat[0] = rz

            elif kind_next == 'protected':
                # La zone suivante est protégée → on force la zone i à rejoindre next
                rx, ry, rz = float(xB[0]), float(yB[0]), float(zB[0])
                Internal.getNodeFromName(ordered[i], 'CoordinateX')[1].flat[-1] = rx
                Internal.getNodeFromName(ordered[i], 'CoordinateY')[1].flat[-1] = ry
                Internal.getNodeFromName(ordered[i], 'CoordinateZ')[1].flat[-1] = rz

            else:
                # Deux zones libres → moyenne acceptable
                mx = 0.5 * (float(xA[-1]) + float(xB[0]))
                my = 0.5 * (float(yA[-1]) + float(yB[0]))
                mz = 0.5 * (float(zA[-1]) + float(zB[0]))
                Internal.getNodeFromName(ordered[i],      'CoordinateX')[1].flat[-1] = mx
                Internal.getNodeFromName(ordered[i],      'CoordinateY')[1].flat[-1] = my
                Internal.getNodeFromName(ordered[i],      'CoordinateZ')[1].flat[-1] = mz
                Internal.getNodeFromName(ordered[next_i], 'CoordinateX')[1].flat[0]  = mx
                Internal.getNodeFromName(ordered[next_i], 'CoordinateY')[1].flat[0]  = my
                Internal.getNodeFromName(ordered[next_i], 'CoordinateZ')[1].flat[0]  = mz

        xA2, yA2, _ = getCoords(ordered[i])
        xB2, yB2, _ = getCoords(ordered[next_i])
        print(f"  [joinInOrder] jonction {i}→{next_i}: "
            f"{orderMap[i][0]} | {orderMap[next_i][0]} | gap_after={gap:.2e}")
        
    # =========================================================================
    # Join (endpoints are aligned)
    # =========================================================================
    # for k, zone in enumerate(ordered):
    #     x, y, _ = getCoords(zone)
    #     print(f"  [joinInOrder DEBUG] zone {k} ({orderMap[k][0]}): "
    #         f"start=({x[0]:.6f},{y[0]:.6f}) end=({x[-1]:.6f},{y[-1]:.6f}) N={len(x)}")
    joined = T.join(ordered)

    C.convertPyTree2File(ordered, f'monitoring/joined.cgns')

    # Remove micro-segments 
    if mergeTol > 0:
        try:
            print(f"  [joinInOrder] Cleaning needed for micro-segment(s)")
            try:
                joined_struct = C.convertBAR2Struct(joined)
            except Exception:
                joined_struct = joined
            x, y, _ = getCoords(joined_struct)
            ds = np.sqrt(np.diff(x)**2 + np.diff(y)**2)
            # Check for very short segments (< 0.1% of mean spacing)
            meanDs = float(np.mean(ds)) if len(ds) > 0 else 1.

            # G.close merges nodes within tol of each other
            tol = 1e-2 * meanDs
            joined_clean = G.close(joined, tol)     # joined_clean = G.close(joined)
            nBefore = getNpts(joined)
            nAfter = getNpts(joined_clean)
            if nBefore != nAfter:
                print(f"  [joinInOrder] cleaned {nBefore - nAfter} micro-segment(s)"
                      f" (tol={tol:.2e})")
            print(f"  [joinInOrder] cleaned {nBefore - nAfter} micro-segment(s)"
                      f" (tol={tol:.2e})")
            
            return joined_clean
        except Exception as e:
            print(f"  [joinInOrder] micro-segment cleanup failed ({e}) — using raw join")

    return joined

#===============================================================================================
# Main iterative loop for mesh smoothing and generation. 
# 1. Normalizes the point density of the profile and repairs short segments.
# 2. Analyzes the topology to separate free zones (to be smoothed) from protected zones.
# 3. Iteratively applies smoothing (T.smooth) to free zones and remeshes the boundary layer (BL) and far-field (T3).
#   → Adapts smoothing sweeps and geometry mapping parameters (peak/valley factors) if mesh quality checks fail.
# Stops when geometry deviation and mesh quality (volume, orthogonality) criteria are met or max iterations reached.

# IN: BAR1: Inner profile boundary
# IN: BAR3: Outer domain boundary
# IN: profileBody: Baseline geometry used as reference for distance checks and projections
# IN: ht, hf, dz
# IN: distThreshold: Optional distance threshold for topology detection
# IN: splitSensib: Sensibility angle for topological splitting (degrees, default 15.0)
# IN: protectAngle: Angle threshold above which nodes are strictly protected (default 160.0)
# IN: borrowFrac: Fraction of protected points to borrow for transitions (default 0.0)
# IN: targetDensity: Desired profile point density per chord unit (default 325.0)

# IN: peakFactor/valleyFactor/roughFactor: Control parameters for the geometry mapper
# IN: roughThreshold: Angle threshold to identify rough geometric zones (default 4.0)

# IN: sweepsInit/sweepsMax/addSweeps: Initial, maximum, and step values for the number of smoothing sweeps
# IN: twoWays, step: T.smooth parameter and step size parameter
# IN: maxIter: Maximum number of mesh generation attempts (default 50)
# IN: meshQuality: Dictionary of thresholds for mesh quality evaluation (default None -> uses global default)
# IN: extruder: boundary layer generation algorithm (0 or 1)

# IN: profileName: Name of the current profile (used for exporting files and folders)
# IN: limZoom: List containing [xmin, xmax, ymin, ymax] used to crop visual exports

# OUT: bestMesh: The final valid 2D combined mesh (BL + T3), or None if it failed entirely
# OUT: bestSmooth: The final 1D smoothed profile contour
# OUT: bestReport: Dictionary containing the mesh quality metrics from the final iteration
#===============================================================================================
def smoothLoop(BAR1, BAR3, profileBody, ht, hf, dz = 1,
               # remesh params (used once)
               distThreshold=None, splitSensib=15.0,
               protectAngle=160.0, borrowFrac=0.0,
               targetDensity=325.,
               peakFactor=1.5, valleyFactor=0.5,
               roughFactor=1.4, roughThreshold=4.0,
               # smoothing iteration
               sweepsInit=5, sweepsMax=1000, addSweeps=50.0,
               twoWays=1, step=1,
               maxIter=50,
               # mesh quality thresholds
               meshQuality=None,
               extruder=0, profileName='horn', limZoom=[0,0,0,0]):
    
    nptsName = getNpts(BAR1)
    exportName (SWEEPSname, nptsName)
    currPF = peakFactor
    currVF = valleyFactor
    currRF = roughFactor
    currSweeps = sweepsInit

    oldPF = peakFactor
    oldVF = valleyFactor
    oldRF = roughFactor
    oldSweeps = sweepsInit
    BLpb = False

    history = {
        'iter': [], 'pf': [], 'vf': [], 'rf': [], 'sweeps': [], 'density': [], 
        'rmaxBL': [], 'vminBL': [], 'rmaxT3': [], 'vminT3': [],
        'nk': [], 'nLayers': [], 'npts': []
    }

    bestMesh = None
    bestSmooth = None
    bestReport = {}
    
    os.makedirs(f"photos/{profileName}", exist_ok=True)
    os.makedirs("monitoring", exist_ok=True)
    exportPathInfo = f"Solveur/Résolution/{profileName}/0_{nptsName}/info"
    exportPathPhoto = f"Solveur/Résolution/{profileName}/0_{nptsName}/photos"

    if meshQuality is None:
        meshQuality = DEFAULT_MESH_QUALITY

    chord = getChord(BAR1)

    if closeGeometryRelative(BAR1, toClose=False):
        BAR1 = closeGeometryRelative(BAR1, toClose=True)

    aBefore = checkArea(BAR1)
    currfactor = 0
    BAR1norm = normaliseDensity(BAR1, chord, targetDensity, currfactor)
    if splitSensib is not None:
        currSplitSensib = splitSensib
    else:
        currSplitSensib = computeSplitSensib(BAR1norm, targetDensity, chord)

    print (f"Chosen split : {currSplitSensib}")
    BAR1norm = repairShortSegments(BAR1norm, hf, factor=5.0)
    currN = getNpts(BAR1norm)
    oldN = currN
    C.convertPyTree2File([BAR1, BAR1norm], f'Solveur/Résolution/{profileName}/0_{nptsName}/profile.plt')   

    subCurves, isProtList, macroTypes, junctionData, profile, scoresMap = \
        detectTopology (BAR1norm, profileBody, distThreshold, currSplitSensib,
                   roughThreshold, protectAngle, chord=chord, borrowFrac=borrowFrac)

    #  STEP 3: smoothing iteration 
    for it in range(1, maxIter + 1):

        print(f"\n{'━'*60}")
        print(f"  SMOOTH ITER {it}/{maxIter}  sweeps={currSweeps}")
        print(f"{'━'*60}")

        pf_str = f"pf={currPF:.3f} (old {oldPF:.3f})" if currPF != oldPF else f"pf={currPF:.3f}"
        vf_str = f"vf={currVF:.3f} (old {oldVF:.3f})" if currVF != oldVF else f"vf={currVF:.3f}"
        rf_str = f"rf={currRF:.3f} (old {oldRF:.3f})" if currRF != oldRF else f"rf={currRF:.3f}"
        npts_str = f"npts={currN} (old {oldN})" if currN != oldN else f"npts={currN}"
        sweeps_str = f"sweeps={currSweeps} (old {oldSweeps})" if currSweeps != oldSweeps else f"sweeps={currSweeps}"

        print(f"\n[smoothLoop Iter] State : {pf_str}, {vf_str}, {rf_str}, {npts_str}, {sweeps_str}")

        finalAreaError = 0.0
        profile = BAR1norm

        # if BLpb == True :
        #     BAR1norm, currSplitSensib = getDensity (BAR1, chord, targetDensity, currfactor)
        #     BAR1norm = repairShortSegments(BAR1norm, hf, factor=5.0)

        #     C.convertPyTree2File([BAR1, BAR1norm], f'Solveur/Résolution/{profileName}/{SWEEPSname}_Sweeps/profile.plt')

        #     subCurves, isProtList, macroTypes, junctionData, profile = detectTopology (BAR1norm, profileBody, distThreshold, currSplitSensib,
        #             roughThreshold, protectAngle, borrowFrac=borrowFrac)

        procFree, procProt, shapeLabels, currentOrderMap = mapper(
                                                            subCurves, isProtList, macroTypes, 
                                                            peakFactor=currPF, valleyFactor=currVF, roughFactor=currRF,
                                                            )      

        processed_ordered = []
        for kind, idx in currentOrderMap:
            if kind == 'free':
                processed_ordered.append(procFree[idx])
            else:
                processed_ordered.append(procProt[idx])

        subcurves_with_shapes = [
            (curve, shape, scoresMap.get(map_pos, {}))
            for map_pos, ((curve, shape), (kind, orig_idx)) in enumerate(zip(
                zip(processed_ordered, shapeLabels), currentOrderMap
            ))
        ]
        visualizeRemeshing(subcurves_with_shapes, figName="remeshing", exportPath=exportPathPhoto, limZoom=limZoom)

        currentFree = procFree
        currentProt = procProt
        joined_before = joinInOrder(currentFree, currentProt, currentOrderMap)
        C.convertPyTree2File(joined_before, f'Solveur/Résolution/{profileName}/0_{nptsName}/join.plt')

        dNpts = getNpts(joined_before) - oldN
        print (f'We end up with npts = {oldN} + {dNpts} = {getNpts(joined_before)}')

        aAfter = checkArea(joined_before)
        if aBefore > 1e-12:
            areaError = (aAfter - aBefore) / aBefore * 100
            print(f"\n[Conservation] Area error due to mapper in smoothLoop: {areaError:.8f} % ")

        # Smooth each free contiguous block together
        # (join adjacent free zones, smooth, split back — keeps endpoints movable)
        currentFree, currentProt, currentOrderMap = smoothFreeBlocks(
            currentFree, currentProt, currentOrderMap,
            currSweeps, twoWays, step
        )
        
        # Reconstruct full profile in order
        smoothed1D = joinInOrder(currentFree, currentProt, currentOrderMap)
        
        smoothed1D = C.convertBAR2Struct(smoothed1D)

        aAfter = checkArea(smoothed1D)
        if aBefore > 1e-12:
            areaError = (aAfter - aBefore) / aBefore * 100
            print(f"\n[Conservation] Area error due to script: {areaError:.8f} % ")
            finalAreaError = areaError

        print ("[GeometryCheck] at end of remeshing and joining")

        checkGeometryIntegrity(smoothed1D, label=f'iter {it}', exportPath=exportPathInfo)
        oldSweeps = currSweeps
        oldPF = currPF
        oldVF = currVF 
        oldRF = currRF
        oldN = currN
        currN = getNpts(smoothed1D)

        history['iter'].append(it)
        history['pf'].append(currPF)
        history['vf'].append(currVF)
        history['rf'].append(currRF)
        history['sweeps'].append(currSweeps)

        history['rmaxT3'].append(float('nan'))
        history['vminT3'].append(float('nan'))

        history['rmaxBL'].append(float('nan'))
        history['vminBL'].append(float('nan'))
        history['nk'].append(float('nan'))
        history['nLayers'].append(float('nan'))
        history['npts'].append(currN)


        updateMeshThreshold (smoothed1D, hf, dz)
        passedDev, reportDev = checkGeoQuality(profile, smoothed1D)
        checkMeshQuality (curve=smoothed1D)
        if passedDev:
                print(f"  ✅ Geometry deviation GOOD at iter {it} (sweeps={currSweeps})")
        else :
            print(f"  ⛔ Geometry deviation BAD at iter {it} (sweeps={currSweeps})")
            currSweeps, currPF, currVF, currRF = getGeometryDeviation (reportDev, currSweeps, currPF, currVF, currRF)
            continue 

        # Reaffine
        # print(f"  GEOMETRY CHECK WITH DEMESHING")
        # L = D.getLength(smoothed1D)
        # currN = getNpts(smoothed1D)
        # # newN = max(10, int(currN * 0.95))
        # newN = currN
        # # cleanD = D.line((0., 0., 0.), (1, 0., 0.), N=newN)
        # h1 = getStepSize(BAR1norm, end="first")
        # h2 = getStepSize(BAR1norm, end="last")

        # cleanD = D.distrib2(smoothed1D, h1, h2, normalized=True, algo=1)
        # smoothed1D = G.map(smoothed1D, cleanD)

        # aAfter = checkArea(smoothed1D)
        # checkGeometryIntegrity(smoothed1D, label=f'iter {it}')


    # return bestMesh, bestSmooth, bestReport
        print(f"  Building mesh...")
        
        # BL Mesher (crash critique → on skip l'itération)
        try:
            BAR1BL, bl, nk, nLayers = mesherBL(Internal.copyTree(smoothed1D), ht, hf, extruder)
            history['nk'][-1] = nk
            history['nLayers'][-1] = nLayers
        except Exception as e:
            print(f"  ⚠️  BL Mesh build failed: {e}")
            currSweeps = int(currSweeps + addSweeps)
            continue

        passedBL, reportBL = checkMeshQuality(bl=bl, thresholds=meshQuality)

        if nk < nLayers:
            passedBL = False

        history['rmaxBL'][-1] = reportBL.get('rmaxBL', float('nan'))
        history['vminBL'][-1] = reportBL.get('vminBL', float('nan'))

        try:
            plotHistoryBL(history, exportPath=exportPathInfo)
        except Exception as e:
            print(f"  ⚠️  Export history failed (non-blocking): {e}")

        if not passedBL:
            print(f"  ⛔ BL Mesh quality BAD at iter {it} (sweeps={currSweeps})")

            currSweeps = currSweeps + 10
            currVF = max(factorThreshold['minF'], currVF * 0.9)
            
            # BLpb = True
            # currfactor = currfactor + 1
            continue
        
        BLpb = False    
        print(f"  ✅ BL Mesh quality GOOD at iter {it} (sweeps={currSweeps})")

        # T3 Mesher (crash → skip iteration)
        try:
            bl, m2, fullMesh = mesherT3(BAR1BL, BAR3, bl, dz)
            genBCMarkers (fullMesh, profileName=profileName)
        except Exception as e:
            print(f"  ⚠️  T3 Mesh build failed: {e}")
            currSweeps, currPF, currVF, currRF = getMeshQuality({}, currSweeps, currPF, currVF, currRF)
            continue

        passedT3, reportT3 = checkMeshQuality(m2=m2, thresholds=meshQuality)

        history['rmaxT3'][-1] = reportT3.get('rmaxT3', float('nan'))
        history['vminT3'][-1] = reportT3.get('vminT3', float('nan'))

        # Export historique 
        try:
            plotHistory(history, exportPath=exportPathInfo)
            exportHistory(history, bestReport={**reportBL, **reportT3},
                        reportDev=reportDev, finalAreaError=finalAreaError,
                        exportPath=exportPathInfo)
        except Exception as e:
            print(f"  ⚠️  Export history failed (non-blocking): {e}")

        if not passedT3:
            print(f"  ⛔ T3 Mesh quality BAD at iter {it} (sweeps={currSweeps})")
            currSweeps, currPF, currVF, currRF = getMeshQuality(reportT3, currSweeps, currPF, currVF, currRF)
            continue

        print(f"  ✅ T3 Mesh quality GOOD at iter {it} (sweeps={currSweeps})")

        bestMesh = fullMesh
        bestSmooth = smoothed1D
        bestReport = {**reportBL, **reportT3} 

        if passedT3:
            print(f"  ✅ Mesh quality OK at iter {it} (sweeps={currSweeps})")
            break

        if it == maxIter:
            print(f"  ⚠️  Max iters reached — returning best mesh so far")
            break

        currSweeps = min(int(currSweeps + addSweeps), sweepsMax)

    if bestSmooth is None:
        print("  ⚠️  No mesh was built — returning remeshed (unsmoothed) profile")
        bestSmooth = joinInOrder(currentFree, currentProt, currentOrderMap)
        bestMesh = None
        bestReport = {'passed': False}
    
    return bestMesh, bestSmooth, bestReport

#===========================================================================
# Join contiguous free sub-curves into blocks, smooth each block together
#===========================================================================
def smoothFreeBlocks(procFree, processedProtected, orderMap, sweeps, twoWays, step):
    
    n = len(orderMap)
    
    newFreeList = []      # Nouvelles zones libres (possiblement fusionnées)
    newProtList = list(processedProtected)  # Zones protégées inchangées
    newOrderMap = []      # Nouvelle carte d'ordre
    
    i = 0
    while i < n:
        kind, idx = orderMap[i]
        
        # =====================================================================
        # CAS 1 : Zone protégée → copier tel quel
        # =====================================================================
        if kind == 'protected':
            newOrderMap.append(('protected', idx))
            i += 1
            continue
        
        # =====================================================================
        # CAS 2 : Zone(s) libre(s) contiguë(s) → fusionner + lisser
        # =====================================================================
        block_idxs = []   # indices dans orderMap
        free_idxs = []    # indices dans procFree
        j = i
        
        # Trouver toutes les zones libres consécutives
        while j < n and orderMap[j][0] == 'free':
            block_idxs.append(j)
            free_idxs.append(orderMap[j][1])
            j += 1
        
        block_zones = [Internal.copyTree(procFree[fi]) for fi in free_idxs]
        
        # Si une seule zone → lisser directement
        if len(block_zones) == 1:
            smoothed_block = T.consSmooth(block_zones[0], sweeps, twoWays, step)
            newFreeList.append(smoothed_block)
            newOrderMap.append(('free', len(newFreeList) - 1))
        
        # Si plusieurs zones → fusionner, lisser, GARDER FUSIONNÉ
        else:
            # Réparer les endpoints avant fusion
            for k in range(len(block_zones) - 1):
                xA, yA, zA = getCoords(block_zones[k])
                xB, yB, zB = getCoords(block_zones[k + 1])
                
                mx = 0.5 * (float(xA[-1]) + float(xB[0]))
                my = 0.5 * (float(yA[-1]) + float(yB[0]))
                mz = 0.5 * (float(zA[-1]) + float(zB[0]))
                
                Internal.getNodeFromName(block_zones[k],     'CoordinateX')[1].flat[-1] = mx
                Internal.getNodeFromName(block_zones[k],     'CoordinateY')[1].flat[-1] = my
                Internal.getNodeFromName(block_zones[k],     'CoordinateZ')[1].flat[-1] = mz
                
                Internal.getNodeFromName(block_zones[k + 1], 'CoordinateX')[1].flat[0]  = mx
                Internal.getNodeFromName(block_zones[k + 1], 'CoordinateY')[1].flat[0]  = my
                Internal.getNodeFromName(block_zones[k + 1], 'CoordinateZ')[1].flat[0]  = mz
            
            # Fusionner
            joined = T.join(block_zones)

            # Trouver les zones PROT adjacentes dans newOrderMap
            prev_map_idx = block_idxs[0] - 1
            next_map_idx = block_idxs[-1] + 1  # = j dans ta boucle

            # Forcer endpoint début → dernière coord de la PROT précédente
            if prev_map_idx >= 0 and orderMap[prev_map_idx][0] == 'protected':
                xP, yP, zP = getCoords(processedProtected[orderMap[prev_map_idx][1]])
                Internal.getNodeFromName(joined, 'CoordinateX')[1].flat[0]  = float(xP[-1])
                Internal.getNodeFromName(joined, 'CoordinateY')[1].flat[0]  = float(yP[-1])
                Internal.getNodeFromName(joined, 'CoordinateZ')[1].flat[0]  = float(zP[-1])

            # Forcer endpoint fin → première coord de la PROT suivante  
            if next_map_idx < n and orderMap[next_map_idx][0] == 'protected':
                xN2, yN2, zN2 = getCoords(processedProtected[orderMap[next_map_idx][1]])
                Internal.getNodeFromName(joined, 'CoordinateX')[1].flat[-1] = float(xN2[0])
                Internal.getNodeFromName(joined, 'CoordinateY')[1].flat[-1] = float(yN2[0])
                Internal.getNodeFromName(joined, 'CoordinateZ')[1].flat[-1] = float(zN2[0])

            smoothed_block = T.consSmooth(joined, sweeps, twoWays, step)
            
            # ✅ GARDER LE BLOC FUSIONNÉ (ne pas split !)
            newFreeList.append(smoothed_block)
            newOrderMap.append(('free', len(newFreeList) - 1))
        
        i = j  # Passer au bloc suivant
    
    return newFreeList, newProtList, newOrderMap

#===========================================================================
# Return (x, y, z) numpy arrays from a PyTree zone.
#===========================================================================
def getCoords(curve):
    x = Internal.getNodeFromName(curve, 'CoordinateX')[1].ravel()
    y = Internal.getNodeFromName(curve, 'CoordinateY')[1].ravel()
    z = Internal.getNodeFromName(curve, 'CoordinateZ')[1].ravel()
    return x, y, z

#===========================================================================
# Number of nodes in a structured 1-D zone
#===========================================================================
def getNpts(curve):
    return Internal.getZoneDim(curve)[1]

#===========================================================================
# Sharpest angle at each node [°] 
# 180° = flat, small = sharp corner.
#===========================================================================
def getSharpestAngles(curve):
    a = D.getSharpestAngle(curve)
    return Internal.getNodeFromName(a, 'alpha')[1].ravel()

#===========================================================================
# True if closed curve is clockwise : vectorized Shoelace formula.
#===========================================================================
def isClockwise(curve):
    x, y, z = getCoords(curve)
    area = np.sum((np.roll(x, -1) - x) * (np.roll(y, -1) + y))
    return area > 0

#===========================================================================
# Standard deviation of (180° - α) at ​​interior nodes.
# 0 = perfectly smooth, large = irregular.
#===========================================================================
def rugosity(alphas):
    # alphas = getSharpestAngles(curve)
    if len(alphas) > 2:
        alphas = alphas[1:-1]
    return float (np.std(180.0 - alphas))

#===========================================================================
# Finds how many points 'k' correspond to a physical 'target Dist'
# If fromEnd=True, counts from the end of the curve to the beginning.
#===========================================================================
def getKFromDist(x, y, targetDist, fromEnd=True):
    distAcc = 0.0
    k = 0
    nPts = len(x)
    if targetDist <= 1e-12: 
        return 0
    
    if fromEnd:
        for j in range(nPts - 1, 0, -1):
            distAcc += np.sqrt((x[j]-x[j-1])**2 + (y[j]-y[j-1])**2)
            k += 1
            if distAcc >= targetDist: break
    else:
        for j in range(0, nPts - 1):
            distAcc += np.sqrt((x[j+1]-x[j])**2 + (y[j+1]-y[j])**2)
            k += 1
            if distAcc >= targetDist: break
            
    # On sécurise : au moins 2 points, max la moitié de la courbe
    return max(2, min(k, nPts // 2))
#===========================================================================
# Calculate the physical distance (the step size 'h') of the first or last segment of a curve.
#===========================================================================
def getStepSize(curve, end="first"):
    x, y, z = getCoords(curve)
    if len(x) < 2: return 1e-5
    
    if end == "first":
        return float(np.sqrt((x[1]-x[0])**2 + (y[1]-y[0])**2 + (z[1]-z[0])**2))
    else:
        return float(np.sqrt((x[-1]-x[-2])**2 + (y[-1]-y[-2])**2 + (z[-1]-z[-2])**2))

#===========================================================================
# Verify the 1D profile is:
# - closed (first == last point within tolerance)
# - has no duplicate consecutive points
# - forms a single connected loop (no dangling branches)
# Prints a report. Returns True if OK, False if issues found.
#===========================================================================
def checkGeometryIntegrity(profile, label='profile', exportPath='photos'):
    import os
    import matplotlib.pyplot as plt
    
    ok = True
 
    # Convert to consistent format
    p = C.initVars(profile, 'CoordinateZ', 0.)
    try:
        p = C.convertBAR2Struct(p)
    except Exception:
        pass
 
    x, y, _ = getCoords(p)
    N = len(x)
 
    # Check closure
    dist2 = (float(x[0]) - float(x[-1]))**2 + (float(y[0]) - float(y[-1]))**2
    ref2 = (float(x[0]) - float(x[1]))**2  + (float(y[0]) - float(y[1]))**2 if N >= 2 else 1.
    isClosed = dist2 < 1e-12 * max(ref2, 1e-30)
 
    # Check for duplicate consecutive points (zero-length segments)
    ds = np.sqrt(np.diff(x)**2 + np.diff(y)**2)
    nDuplicates = int(np.sum(ds < 1e-14))
 
    # Check for very short segments (< 0.1% of mean spacing)
    meanDs = float(np.mean(ds)) if len(ds) > 0 else 1.
    
    # ON IDENTIFIE LES MAUVAIS SEGMENTS
    dup_mask = ds < 1e-14
    nDuplicates = int(np.sum(dup_mask))
    
    bad_mask = ds < 1e-2 * meanDs

    print ("tol = ", 1e-2 * meanDs)
    nVeryShort = int(np.sum(bad_mask))
 
    print(f"\n[checkGeometryIntegrity] {label}")
    print(f"  N = {N}  |  closed = {isClosed}  |  duplicates = {nDuplicates}"
          f"  |  very-short segments = {nVeryShort}")
 
    if not isClosed:
        gap = float(np.sqrt(dist2))
        print(f"  ⚠️  NOT CLOSED — gap = {gap:.2e}  (should be < {1e-6*float(np.sqrt(ref2)):.2e})")
        ok = False
        
    if nDuplicates > 0:
        print(f"  ⚠️  {nDuplicates} duplicate consecutive point(s) — may cause mesher issues")
        ok = False
        
    if nVeryShort > 0:
        print(f"  ⚠️  {nVeryShort} very-short segment(s) — may cause T3mesher2D failures")
        ok = False
        
    # ====================================================================
    # GÉNÉRATION DE L'IMAGE DEBUG (ZOOM SUR LES DÉFAUTS)
    # ====================================================================
    if nDuplicates > 0 or nVeryShort > 0:
        try:
            dup_indices = np.where(dup_mask)[0]
            bad_indices = np.where(bad_mask)[0]
            
            os.makedirs(exportPath, exist_ok=True)
            fig, ax = plt.subplots(figsize=(8, 6))
            
            # Tracer tout le profil en gris fin
            ax.plot(x, y, color='gray', linestyle='-', linewidth=0.5, label='Profil complet', alpha=0.7)
            
            # Mettre en évidence les points en DOUBLE (en gros et en ORANGE)
            for i in dup_indices:
                ax.plot(x[i:i+2], y[i:i+2], color='orange', linestyle='-', marker='X', markersize=8, linewidth=2)
                ax.text(x[i], y[i], f" Dup {i}", color='orange', fontsize=9, ha='right', va='top')

            # Mettre en évidence les MICRO-SEGMENTS (en gros et en ROUGE)
            for i in bad_indices:
                ax.plot(x[i:i+2], y[i:i+2], color='red', linestyle='-', marker='o', markersize=6, linewidth=2)
                ax.text(x[i], y[i], f" Pt {i}", color='red', fontsize=9, ha='left', va='bottom')
            
            # Légende dynamique
            if nDuplicates > 0:
                ax.plot([], [], color='orange', marker='X', markersize=8, linewidth=2, label='Doublons (dist ≈ 0)')
            if nVeryShort > 0:
                ax.plot([], [], color='red', marker='o', markersize=6, linewidth=2, label='Micro-segments')
            
            ax.set_aspect('equal')
            ax.legend(loc='best')
            ax.set_title(f"DEBUG: {nDuplicates} Doublon(s) | {nVeryShort} Micro-segment(s) - {label}")
            
            # On zoom autour du TOUT PREMIER défaut détecté pour le voir clairement
            all_bad = np.concatenate((dup_indices, bad_indices))
            if len(all_bad) > 0:
                first_bad = np.min(all_bad)
                marge = 0.05 * (np.max(x) - np.min(x)) # Marge de 5% de la corde
                ax.set_xlim(x[first_bad] - marge, x[first_bad] + marge)
                ax.set_ylim(y[first_bad] - marge, y[first_bad] + marge)

            safe_label = str(label).replace(" ", "_")
            outPath = os.path.join(exportPath, f"DEBUG_geometrie_{safe_label}.png")
            #plt.savefig(outPath, dpi=200, bbox_inches='tight')
            plt.close(fig)
            print(f"  [Visualisation] 📸 Photo de debug générée : {outPath} (Regarde le zoom !)")
        except Exception as e:
            print(f"  [Visualisation] Échec de la génération de l'image de debug : {e}")

    if ok:
        print(f"  ✓ Geometry OK")
 
    return ok

#===========================================================================
# Signed interior angle [°] at the junction between two sub-curves.
# sign +1 = convex peak (outward curve)
# sign -1 = concave trough (inward curve)
# OUT: magnitude, sign
#===========================================================================
def junctionAngleSigned(subA, subB):
    xA, yA, zA = getCoords(subA)
    xB, yB, zB = getCoords(subB)
    if len(xA) < 2 or len(xB) < 2:
        return 180.0, 0

    vA = np.array([xA[-1] - xA[-2], yA[-1] - yA[-2]], dtype=float)
    vB = np.array([xB[1]  - xB[0],  yB[1]  - yB[0]],  dtype=float)

    nA, nB = np.linalg.norm(vA), np.linalg.norm(vB)
    if nA < 1e-14 or nB < 1e-14:
        return 180.0, 0

    cosT = np.clip(np.dot(vA, vB) / (nA * nB), -1.0, 1.0)
    magnitude = 180.0 - float(np.degrees(np.arccos(cosT)))
    cross = vA[0]*vB[1] - vA[1]*vB[0]

    if abs(cross) < 1e-14:
        return magnitude, 0
    sign = -1 if cross > 0 else +1
    return float(magnitude), int(sign)

#===========================================================================
# Calculates the orthogonal (unsigned) distance of each node of the curve to the reference profile `profileBody`.
# OUT: numpy array 1-D of nodal distances
#===========================================================================
def getWallDistances(curve, profileBody):

    distNode = Internal.getNodeFromName(curve, 'TurbulentDistance')
    if distNode is not None:
        return distNode[1].ravel()

    # Fallback to compute if missing
    result = D2W.distance2Walls(curve, [profileBody],
                                type='ortho',
                                loc='nodes',
                                signed=0)
    
    dist = Internal.getNodeFromName(result, 'TurbulentDistance')[1].ravel()
    return dist

#===========================================================================
# Calculates the chord from a profile (or profileBody)
#===========================================================================
def getChord(profile):

    xN, _, _ = getCoords(profile)
    chord = np.max(xN) #float(np.max(xN) - np.min(xN))

    # print(f"\n[Loop]  chord = {chord:.4f}")

    return chord

#===========================================================================
# Adjusts the number of points on a subcurve based on its pre-calculated type
# Convex peak ('p'): ×peakFactor
# Concave trough ('v'): ×valleyFactor
# Rough ('r'): ×roughFactor
#===========================================================================
def adaptiveRemesh(curve, zType, peakFactor=1.5, valleyFactor=0.5, roughFactor=1.4):
        
    info = {}
    L = D.getLength(curve)
    N = getNpts(curve)
    grading = 1.0

    if zType == 'p':
        factor = peakFactor
        label = f"MACRO-PIC convexe → ×{factor:.3f}"
        info['shape'] = 'p'
    elif zType == 'v':
        factor = valleyFactor
        label = f"MACRO-CREUX concave → ×{factor:.3f}"
        info['shape'] = 'v'
    elif zType == 'r':
        factor = roughFactor
        label = f"MACRO-RUGOSITÉ → ×{factor:.3f}"
        info['shape'] = 'r'
    else:
        # factor = 1.0
        # label = f"MACRO-NEUTRE → inchangé"
        # info['shape'] = 'n'
        factor = peakFactor
        label = f"MACRO-NEUTRE → MACRO-PIC convexe → ×{factor:.3f}"
        info['shape'] = 'p'

    nNew = max(4, int(N * factor))

    h1, h2 = getH(L, nNew, grading)   
    print(f"        Remesh: {N:4d} → {nNew:4d} pts  |  {label}")

    # Remaillage uniforme sur tout le macro-bloc
    # ref = D.line((0., 0., 0.), (1., 0., 0.), N=nNew)
    d = D.distrib2(curve, h1, h2, normalized=True, algo=1)
    
    if nNew != N:
        return G.map(curve, d), info["shape"]  
    else :
        return curve, info["shape"]
#===========================================================================
# GARDÉE pour si besoin d'une autre stratégie de remaillage 1D (plus complex)
#===========================================================================
def adaptiveRemesh2(curve, angleLeft=180.0, signLeft=0,
                          angleRight=180.0, signRight=0,
                          peakFactor=1.5, valleyFactor=0.5, roughFactor = 1.4,
                          roughThreshold=5.0, protectAngle=120.0):
        
    info = {}
    L = D.getLength(curve)
    N = getNpts(curve)
    rug = rugosity(curve)
    nNew = N

    # Jonction la plus critique + côté où elle se trouve
    if angleLeft <= angleRight:
        sharpAngle, sharpSign = angleLeft,  signLeft
        sharpIsLeft = True
    else:
        sharpAngle, sharpSign = angleRight, signRight
        sharpIsLeft = False

    isSharp = sharpAngle < protectAngle
    isPeak = isSharp and sharpSign == -1
    isValley = isSharp and sharpSign == +1

    if isPeak:
        factor = peakFactor
        label = f"PIC convexe (α={sharpAngle:.1f}°) → ×{factor:.1f}"
        grading = 1.25
        info['shape'] = 'p'

    elif isValley:
        factor = valleyFactor
        label = f"CREUX concave (α={sharpAngle:.1f}°) → ×{factor:.1f}"
        grading = 1.0
        info['shape'] = 'v'

    elif rug > roughThreshold:
        factor = roughFactor
        label = f"RUGOSITÉ (σ={rug:.1f}°) → ×{factor:.1f}"
        grading = 1.0
        sharpIsLeft = False  # uniforme, pas d'orientation à inverser
        info['shape'] = 'r'

    else:
        factor = 1.0
        label = f"neutre (σ={rug:.1f}°, α={sharpAngle:.1f}°)"
        grading = 1.0
        sharpIsLeft = False
        info['shape'] = 'n'

    nNew = max(3, int(N * factor))
    h1, h2 = getH(L, nNew, grading)   
    
    if not sharpIsLeft:
        h1, h2 = h2, h1  

    print(f"      Remesh: {N:4d} → {nNew:4d} pts  |  {label}")
    d = D.distrib2(curve, h1, h2, normalized=True, algo=1)
    if nNew != N:
        return G.map(curve, d), info["shape"]  
    else :
        return curve, info["shape"]

#===========================================================================
# Calculates the h1, h2 using the sum of geometric series with a grading r
#===========================================================================
def getH (L, N, r=1.1):
    n = N - 1
    if n <= 0: return L, L
    
    # If the ratio is 1.0 (no grading, uniform mesh)
    if abs(r - 1.0) < 1e-6:
        h = L / n
        return h, h
        
    # Calcul exact de la première maille (h1)
    h1 = L * (1.0 - r) / (1.0 - r**n)
    
    # Calcul exact de la dernière maille (h2)
    h2 = h1 * (r**(n - 1))
    
    return h1, h2

#===========================================================================
# Protected sub-curve (leading edge, tip) if:
# 1 - minimum angle < protectAngle
# 2 - (optional) curvilinear length < maxCornerLength (0 = disabled) (to develop to protected profile profil)
#===========================================================================
def isProtected(curve, rug, alphaCache, profileBody=None, distThreshold=None, roughThreshold = 4.0, protectAngle=60.0, maxCornerLength=0.0):
    
    info = {}

    info['rugosity'] = False
    
    # Si la rugosité est forte, on veut protéger la zone du lissage
    isRugous = (rug > roughThreshold)

    if profileBody is not None:

        dist = getWallDistances(curve, profileBody)
        distMax = float(np.max(dist))
        distMean = float(np.mean(dist))
        info['mode'] = 'distance'
        info['distMax'] = distMax
        info['distMean'] = distMean

        # Estimation automatique du seuil si non fourni
        if distThreshold is None:
            # Corde ≈ étendue en X du profile (normalisé → corde ≈ 1)
            chord = getChord(profileBody)# float(np.max(xN) - np.min(xN))
            distThreshold = 0.001 * chord           # 0.1 % de la corde
            info['distThreshold_auto'] = distThreshold

        protected = bool(distMax <= float(distThreshold))
        reason = f"distMax ({distMax:.4e}) <= seuil ({distThreshold:.4e})" if protected else f"distMax ({distMax:.4e}) > seuil ({distThreshold:.4e})"

        if isRugous:
            protected = False
            info['rugosity'] = True
            reason = f"Rugosité forte ({rug:.2f} > {roughThreshold})"

        xN, _, _ = getCoords(curve)
        if np.mean(xN) >= chord * 0.2: 
            protected = True
            info['forced'] = True
            info['rugosity'] = False
            reason = f"Forcé (x_mean {np.mean(xN):.4f} >= 20% corde {chord*0.2:.4f})"

        length = D.getLength(curve)
        if length > 0.2 * chord:
            protected = True
            reason = f"Forcé (longueur {length:.4f} > 20% corde {chord*0.2:.4f})"

        alphas = alphaCache
        if len(alphas) > 2:
            alphas = alphas[1:-1]
        info['minAngle'] = float(np.min(alphas))

        print(f"  [isProtected - Distance] -> {protected} | Raison: {reason}")
        return protected, info
    
    alphas = alphaCache
    if len(alphas) > 2:
        alphas = alphas[1:-1]
    minAngle = float(np.min(alphas))
    hasSharp = minAngle < protectAngle
    if maxCornerLength > 0:
        protected = hasSharp and D.getLength(curve) < maxCornerLength
    else:
        protected = hasSharp

    info['mode'] = 'angle'
    info['minAngle'] = minAngle

    if isRugous:
        protected = False
        info['rugosity'] = True

    return protected, info

#==================================================================
# Force first and last nodes to exact junction coordinates
#==================================================================
def forceEndpoints(curve, pStart, pEnd):
    for coord, idxStart, idxEnd in [
        ('CoordinateX', pStart[0], pEnd[0]),
        ('CoordinateY', pStart[1], pEnd[1]),
        ('CoordinateZ', pStart[2], pEnd[2]),
    ]:
        node = Internal.getNodeFromName(curve, coord)
        if node is not None:
            node[1].flat[0] = idxStart
            node[1].flat[-1] = idxEnd
    return curve

#===========================================================================
#                     VALIDATION / INDICATOR RESULTS
#===========================================================================
    
#==============================
# check quality of output mesh
#==============================
def checkQuality(m):
    # check positive volumes
    G._getVolumeMap(m)
    volmin = C.getMinValue(m, 'centers:vol')
    volmax = C.getMaxValue(m, 'centers:vol')
    volmean = C.getMeanValue(m, 'centers:vol')
    if volmin > 0: print("INFO: all positive volumes.")
    else: print("Warning: negative volume detected in final mesh.")

    return 0

#===============================================
# check deviation of profile from input profile
#===============================================
def checkDeviation(rough, smoothed, refLen, method = "ortho"):

    rIce = Internal.copyTree(rough) #rough Ice
    sIce = Internal.copyTree(smoothed) #smooth ice

    tol = 1e-6

    rIce = C.initVars(rIce, 'CoordinateZ', 0.)
    sIce = C.initVars(sIce, 'CoordinateZ', 0.)

    try:
        rIce = C.convertBAR2Struct(rIce)
        sIce = C.convertBAR2Struct(sIce)
    except Exception as e:
        print(f"[checkDeviation] convertBAR2Struct skipped: {e}")

    try:
        rIce = T.join(rIce)
        sIce = T.join(sIce)
    except Exception as e:
        raise RuntimeError(f"[checkDeviation] T.join failed — cannot compute distances: {e}")

    rIceCellN = C.initVars(rIce, 'cellN', 1.)
    distances = D2W.distance2Walls(sIce, [rIceCellN], signed=1, type=method, loc='nodes')
    values = Internal.getNodeFromName(distances, 'TurbulentDistance')[1].ravel()

    maxInward = float(np.min(values)) / refLen
    maxOutward = float(np.max(values)) / refLen
    stdDev = float(np.std(values)) / refLen
    meanDev = float(np.mean(values)) / refLen
    
    aSmooth = checkArea(sIce)
    aRough = checkArea(rIce)
    if aRough > 1e-12:
        areaError = (aSmooth - aRough) / aRough * 100
        print(f"\n[Conservation] Area error : {areaError:.8f} % ")

    sumInward = float(np.sum(np.abs(values[values < -tol])))
    sumOutward = float(np.sum(values[values > tol]))
    sumStable = float(np.sum(np.abs(values[np.abs(values) <= tol])))
    total = sumInward + sumOutward + sumStable

    return dict(maxInward=maxInward, maxOutward=maxOutward,
                stdDev=stdDev, meanDev=meanDev)

#======================================================================
# Check if a PyTree curve is open or closed using a relative tolerance.
# If open, joins the end to the start.
#======================================================================
def closeGeometryRelative(curve, toClose=False):

    dim = Internal.getZoneDim(curve)
    if dim[0] == 'Unstructured':
        # Conversion coûteuse seulement si c'est non-structuré
        curve2Check = C.convertBAR2Struct(curve)
    else:
        curve2Check = curve

    cx, cy, cz = getCoords(curve2Check)
    n1, n1Ngb, n2 = 0, 1, -1

    refL2 = (cx[n1] - cx[n1Ngb])**2 + (cy[n1] - cy[n1Ngb])**2 + (cz[n1] - cz[n1Ngb])**2
    dist2 = (cx[n1] - cx[n2])**2    + (cy[n1] - cy[n2])**2    + (cz[n1] - cz[n2])**2
    
    isOpen = (dist2 >= 1e-12 * refL2)
    if not toClose:
        return isOpen
    
    if not isOpen:
        print(f"Closed geometry (Actual edge points: {n1} and {n2}).")
    else:
        print(f"Open Geometry (Actual edge points: {n1} and {n2}). Creation of the connection...")
        start = (cx[n1], cy[n1], cz[n1])
        end = (cx[n2], cy[n2], cz[n2])
        segment = D.line(end, start, N=5)
        curve2Check = T.join(curve2Check, segment)
        print("-> Connection successful")
    
    return curve2Check

#==========================================
# check area of profile from input profile 
#==========================================
def checkArea(curve):

    profile = Internal.copyTree(curve)

    try:
        profile = C.initVars(profile, 'CoordinateZ', 0.)

        # Ensure BAR (unstructured) for T3mesher2D
        try:
            profile = C.convertBAR2Struct(profile)
        except Exception:
            pass
        profile = C.convertArray2Tetra(profile)


        # Close if open
        if closeGeometryRelative(profile, toClose=False):
            profile = closeGeometryRelative(profile, toClose=True)  

        profile = C.convertArray2Tetra(profile)

        meshed = G.T3mesher2D(profile, triangulateOnly=1, metricInterpType=0)

        # PyTree integ with node field (avoids node2Center warning)
        C._initVars(meshed, 'field', 1.)
        area = float(P.integ(meshed, var='field')[-1])

    except Exception as e:
        print(f"  [checkArea] failed: {e}")
        area = 0.0

    # print(f"  Area = {area:.8f}")
    return area

# =============================================================================
# Check quality of BL (quad) and T3 meshes using G.checkMesh.
#   Returns (passed, report_dict).
#   passed = True if all thresholds are met.
# =============================================================================
def checkMeshQuality(bl=None, m2=None, curve=None, thresholds=None):
    if thresholds is None:
        thresholds = DEFAULT_MESH_QUALITY

    report = {
        'passed': True,
        'failures': {} 
    }
    passed = True

    if curve is not None:

        # Extraction Boundary Layer (C)
        infoC = G.checkMesh(curve)
        
        vcritC = int(infoC.get('vcrit', 0))
        vminC, vmaxC, vmeanC = float(infoC.get('vmin', 0.)), float(infoC.get('vmax', 0.)), float(infoC.get('vmean', 0.))
        
        acritC = int(infoC.get('acrit', 0))
        aminC, amaxC, ameanC = float(infoC.get('amin', 0.)), float(infoC.get('amax', 0.)), float(infoC.get('amean', 0.))
        
        rcritC = int(infoC.get('rcrit', 0))
        rminC, rmaxC, rmeanC = float(infoC.get('rmin', 0.)), float(infoC.get('rmax', 0.)), float(infoC.get('rmean', 0.))

        # # Évaluations
        # okAcritC = acritC <= thresholds['maxAcrit']
        # okVcritC = vcritC <= thresholds['maxVcrit']
        # okVminC = vminC >= thresholds['minVC']
        # okAminC = aminC  >= thresholds['minAngleC']
        # okRmaxC = rmaxC  <= thresholds['maxRC']

        # # Enregistre
        # if not okAcritC: report['failures']['acrit'] = acritC
        # if not okVcritC: report['failures']['vcrit'] = vcritC
        # if not okVminC:  report['failures']['vmin'] = vminC
        # if not okAminC:  report['failures']['amin'] = aminC
        # if not okRmaxC:  report['failures']['rmax'] = rmaxC

        # report.update({
        #     'vcritC': vcritC, 'vminC': vminC, 'vmaxC': vmaxC, 'vmeanC': vmeanC,
        #     'acritC': acritC, 'aminC': aminC, 'amaxC': amaxC, 'ameanC': ameanC,
        #     'rcritC': rcritC, 'rminC': rminC, 'rmaxC': rmaxC, 'rmeanC': rmeanC
        # })

        # passedC = okAcritC and okVcritC and okRmaxC and okVminC and okAminC
        # passed = passed and passedC

        print(f"\n  C Mesh Quality Check:")
        print(f"      Volume : vcrit = {vcritC:<8}  | vmin = {vminC:.2e}  | vmax = {vmaxC:.2e} | vmean = {vmeanC:.2e}")
        print(f"      Angle : acrit = {acritC:<8}  | amin = {aminC:>5.1f}° | amax = {amaxC:>5.1f}° | amean = {ameanC:>5.1f}°")
        print(f"      Regularity : rcrit = {rcritC:<8} (info)       | rmin = {rminC:.2e} | rmax = {rmaxC:.2e} | rmean = {rmeanC:.2e}")



    if bl is not None:

        # Extraction Boundary Layer (BL)
        infoBL = G.checkMesh(bl)
        
        vcritBL = int(infoBL.get('vcrit', 0))
        vminBL, vmaxBL, vmeanBL = float(infoBL.get('vmin', 0.)), float(infoBL.get('vmax', 0.)), float(infoBL.get('vmean', 0.))
        
        acritBL = int(infoBL.get('acrit', 0))
        aminBL, amaxBL, ameanBL = float(infoBL.get('amin', 0.)), float(infoBL.get('amax', 0.)), float(infoBL.get('amean', 0.))
        
        rcritBL = int(infoBL.get('rcrit', 0))
        rminBL, rmaxBL, rmeanBL = float(infoBL.get('rmin', 0.)), float(infoBL.get('rmax', 0.)), float(infoBL.get('rmean', 0.))

        # Évaluations
        okAcritBL = acritBL <= thresholds['maxAcrit']
        okVcritBL = vcritBL <= thresholds['maxVcrit']
        okVminBL = vminBL >= thresholds['minVBL']
        okAminBL = aminBL  >= thresholds['minAngleBL']
        okRmaxBL = rmaxBL  <= thresholds['maxRBL']

        # Enregistre
        if not okAcritBL: report['failures']['acrit'] = acritBL
        if not okVcritBL: report['failures']['vcrit'] = vcritBL
        if not okVminBL:  report['failures']['vmin'] = vminBL
        if not okAminBL:  report['failures']['amin'] = aminBL
        if not okRmaxBL:  report['failures']['rmax'] = rmaxBL

        report.update({
            'vcritBL': vcritBL, 'vminBL': vminBL, 'vmaxBL': vmaxBL, 'vmeanBL': vmeanBL,
            'acritBL': acritBL, 'aminBL': aminBL, 'amaxBL': amaxBL, 'ameanBL': ameanBL,
            'rcritBL': rcritBL, 'rminBL': rminBL, 'rmaxBL': rmaxBL, 'rmeanBL': rmeanBL
        })

        # passedBL = okVcritBL and okRmaxBL and okVminBL and okAminBL and okAcritBL
        passedBL = okVminBL and okRmaxBL
        passed = passed and passedBL

        print(f"\n  BL Mesh Quality Check:")
        print(f"      Volume : vcrit = {vcritBL:<8} {'✓' if okVcritBL else '✗'} (tol <= {thresholds['maxVcrit']}) | vmin = {vminBL:.2e} {'✓' if okVminBL else '✗'} (tol >= {thresholds['minVBL']}) | vmax = {vmaxBL:.2e} | vmean = {vmeanBL:.2e}")
        print(f"      Angle : acrit = {acritBL:<8} {'✓' if okAcritBL else '✗'} (tol <= {thresholds['maxAcrit']:<2}) | amin = {aminBL:>5.1f}° | amax = {amaxBL:>5.1f}° | amean = {ameanBL:>5.1f}°")
        print(f"      Regularity : rcrit = {rcritBL:<8} (info)       | rmin = {rminBL:.2e} | rmax = {rmaxBL:.2e} {'✓' if okRmaxBL else '✗'} (tol <= {thresholds['maxRBL']})| rmean = {rmeanBL:.2e}")

        regMap = G.getRegularityMap(bl)
        a = C.getMaxValue(regMap, 'centers:regularity')
        print (f"Comparaison avec régularitymap : rmax = {rmaxBL} vs rmaxReg = {a}")

    if m2 is not None:

        # Extraction Unstructured (T3)
        infoT3 = G.checkMesh(m2)

        vcritT3 = int(infoT3.get('vcrit', 0))
        vminT3, vmaxT3, vmeanT3 = float(infoT3.get('vmin', 0.)), float(infoT3.get('vmax', 0.)), float(infoT3.get('vmean', 0.))
        
        acritT3 = int(infoT3.get('acrit', 0))
        aminT3, amaxT3, ameanT3 = float(infoT3.get('amin', 0.)), float(infoT3.get('amax', 0.)), float(infoT3.get('amean', 0.))
        
        rcritT3 = int(infoT3.get('rcrit', 0))
        rminT3, rmaxT3, rmeanT3 = float(infoT3.get('rmin', 0.)), float(infoT3.get('rmax', 0.)), float(infoT3.get('rmean', 0.))

        okAcritT3 = acritT3 <= thresholds['maxAcrit']
        okVcritT3 = vcritT3 <= thresholds['maxVcrit']
        okVminT3 = vminT3 >= thresholds['minVT3']
        okAminT3 = aminT3  >= thresholds['minAngleT3']
        okRmaxT3 = rmaxT3  <= thresholds['maxRT3']

        if not okAcritT3: report['failures']['acrit'] = acritT3
        if not okVcritT3: report['failures']['vcrit'] = vcritT3
        if not okVminT3:  report['failures']['vmin'] = vminT3
        if not okAminT3:  report['failures']['amin'] = aminT3
        if not okRmaxT3:  report['failures']['rmax'] = rmaxT3

        report.update({
            'vcritT3': vcritT3, 'vminT3': vminT3, 'vmaxT3': vmaxT3, 'vmeanT3': vmeanT3,
            'acritT3': acritT3, 'aminT3': aminT3, 'amaxT3': amaxT3, 'ameanT3': ameanT3,
            'rcritT3': rcritT3, 'rminT3': rminT3, 'rmaxT3': rmaxT3, 'rmeanT3': rmeanT3
        })

        regMap = G.getRegularityMap(m2)
        a = C.getMaxValue(regMap, 'centers:regularity')
        print (f"Comparaison avec régularitymap : rmax = {rmaxT3} vs rmaxReg = {a}")

        print ('On trouve ', report.get('rmaxT3'))
        # passedT3 = okVcritT3 and okRmaxT3 and okVminT3 
        passedT3 = okRmaxT3 and okVminT3 
        passed = passed and passedT3

        print(f"\n  Mesh Quality Check:")
        print(f"    --- Unstructured (T3) ---")
        print(f"      Volume : vcrit = {vcritT3:<8} {'✓' if okVcritT3 else '✗'} (tol <= {thresholds['maxVcrit']}) | vmin = {vminT3:.2e} {'✓' if okVminT3 else '✗'} (tol >= {thresholds['minVT3']}) | vmax = {vmaxT3:.2e} | vmean = {vmeanT3:.2e}")
        print(f"      Angle : acrit = {acritT3:<8} {'✓' if okAcritT3 else '✗'} (tol <= {thresholds['maxAcrit']:<2}) | amin = {aminT3:>5.1f}° | amax = {amaxT3:>5.1f}° | amean = {ameanT3:>5.1f}°")
        print(f"      Regularity : rcrit = {rcritT3:<8} (info)       | rmin = {rminT3:.2e} | rmax = {rmaxT3:.2e} {'✓' if okRmaxT3 else '✗'} (tol <= {thresholds['maxRT3']}) | rmean = {rmeanT3:.2e}")

    report['passed'] = passed
    return passed, report
# =============================================================================
# Check quality of BL (quad) and T3 meshes using G.checkMesh.
#   Returns (passed, report_dict).
#   passed = True if all thresholds are met.
# =============================================================================
def checkGeoQuality(profile, smoothed1D, devThresholds=None):

    if devThresholds is None:
        devThresholds = defaultGeoQuality()

    report = {
        'passed': True,
        'failures': {} 
    }
    
    chord = getChord(smoothed1D)
    refLen = chord
    infoDev = checkDeviation(profile,smoothed1D, refLen)

    maxInward, maxOutward, stdDev, meanDev = float(infoDev.get('maxInward', 0.)), float(infoDev.get('maxOutward', 0.)), float(infoDev.get('stdDev', 0.)), float(infoDev.get('meanDev', 0.))
    
    # Évaluations
    okInward = abs(maxInward) <= devThresholds['maxInward']
    okOutward = maxOutward <= devThresholds['maxOutward']
    okStd = stdDev <= devThresholds['stdDev']
    okMean = abs(meanDev) <= devThresholds['meanDev']

    # Enregistre
    if not okInward: report['failures']['maxInward'] = maxInward
    if not okOutward: report['failures']['maxOutward'] = maxOutward
    if not okStd: report['failures']['stdDev'] = stdDev
    if not okMean: report['failures']['meanDev'] = meanDev

    report.update({'maxInward': maxInward, 'maxOutward': maxOutward, 'stdDev': stdDev, 'meanDev': meanDev})

    passedDev = okInward and okOutward and okStd and okMean
    
    print(f"\n  Geometry Deviation Check (Smoothing vs Rough Ice):")
    print(f"      Inward : max = {abs(maxInward):>+7.4f} {'✓' if okInward else '✗'} (tol <= {devThresholds['maxInward']:>+7.4f})")
    print(f"      Outward : max = {maxOutward:>+7.4f} {'✓' if okOutward else '✗'} (tol <= {devThresholds['maxOutward']:>+7.4f})")
    print(f"      Dispersion : std = {stdDev:>7.4f} {'✓' if okStd else '✗'} (tol <= {devThresholds['stdDev']:>7.4f})")
    print(f"      Mean deviation : mean = {meanDev:>7.4f} {'✓' if okMean else '✗'} (tol <= {devThresholds['stdDev']:>7.4f})")

    report['passed'] = passedDev
    return passedDev, report

#===========================================================================
# Visualize subcurves colored by remeshing decision (before smoothing)
#   'p' (convex peak)    → green
#   'v' (concave trough) → blue
#   'r' (rough)          → orange
#   'n' / protected      → red
#===========================================================================
def visualizeRemeshing(subcurves_with_shapes,  figName="remeshing", exportPath="photos", limZoom=[0,1,0,1]):
    import os
    import numpy as np
    os.makedirs(exportPath, exist_ok=True)

    COLOR_MAP = {
        'p':         ('green',  'Convexe (pic)'),
        'v':         ('blue',   'Concave (creux)'),
        'r':         ('orange', 'Rugueux'),
        'n':         ('gray',   'Neutre'),
        'protected': ('red',    'Protégé / Non modifié'),
    }

    try:
        import CPlot.Decorator as Decorator
        _use_decorator = True
    except ImportError:
        import matplotlib.pyplot as plt
        _use_decorator = False

    # Extraction des coordonnées et Bounding Box globale 
    allX, allY = [], []
    for curve, _, _ in subcurves_with_shapes:
        x, y, _ = getCoords(curve)
        allX.extend(x)
        allY.extend(y)
        
    allX = np.array(allX)
    allY = np.array(allY)

    # LE FILTRE ULTIME ANTI-NAN
    # On crée un masque qui ne garde que les points qui SONT des nombres
    valid_mask = ~(np.isnan(allX) | np.isnan(allY))
    
    # Si tout est corrompu (cas extrême), on annule le dessin proprement
    if not np.any(valid_mask):
        print("  [visualizeRemeshing] ⚠️ Image ignorée : données corrompues (NaN).")
        return

    # On applique le masque pour nettoyer les tableaux
    allX_clean = allX[valid_mask]
    allY_clean = allY[valid_mask]

    # On calcule les limites sur les données saines
    xmin, xmax = np.min(allX_clean), np.max(allX_clean)
    ymin, ymax = np.min(allY_clean), np.max(allY_clean)

    # Calcul des limites (Global et Zoom) 

    # Vue Globale
    mX = (xmax - xmin) * 0.05
    mY = (ymax - ymin) * 0.10
    xlimGlobal = [xmin - mX, xmax + mX]
    ylimGlobal = [ymin - mY, ymax + mY]

    # Vue Zoom (Bord d'attaque)
    # On trouve automatiquement le point avec le X minimum
    indiceBordAttaque = allX.argmin()
    x_BA = allX[indiceBordAttaque]
    y_BA = allY[indiceBordAttaque]

    xlimZoom = [limZoom[0], limZoom[1]]
    ylimZoom = [limZoom[2], limZoom[3]]

    def draw_and_save(xlim, ylim, suffix):
        dx = xlim[1] - xlim[0]
        dy = ylim[1] - ylim[0]

        if dx < 1e-6: dx = 1e-6
        if dy < 1e-6: dy = 1e-6

        if dy < 0.05 * dx:
            mid_y = (ylim[0] + ylim[1]) / 2.0
            dy = 0.05 * dx
            ylim = [mid_y - dy/2.0, mid_y + dy/2.0]

        if dx < 0.05 * dy:
            mid_x = (xlim[0] + xlim[1]) / 2.0
            dx = 0.05 * dy
            xlim = [mid_x - dx/2.0, mid_x + dx/2.0]

        ratio = dx / dy
        hauteur_calculee = max(3.0, min(10.0 / ratio, 15.0))
        taille_figure = (10, hauteur_calculee)

        if _use_decorator:
            fig, ax = Decorator.createSubPlot(box=False, figsize=taille_figure)
        else:
            fig, ax = plt.subplots(figsize=taille_figure)

        legend_seen = set()

        # ↓ la boucle est ICI dans draw_and_save, pas à l'extérieur
        for curve, shape, scores in subcurves_with_shapes:
            x, y, _ = getCoords(curve)   # ← x et y extraits ici, pour ce segment
            color, label = COLOR_MAP.get(shape, ('gray', f'inconnu ({shape})'))

            full_label = label

            lbl = full_label if full_label not in legend_seen else '_nolegend_'
            legend_seen.add(full_label)

            ax.plot(x, y, color=color, linestyle='-', marker='o',
                    markersize=2, linewidth=1.5, label=lbl, zorder=2)

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_aspect('equal')

        titre = "Remeshing decision" + (" (Zoom)" if suffix == "zoom" else " (Global)")
        ax.set_title(titre, fontsize=9)
        ax.legend(loc='best', fontsize=7, framealpha=0.9)

        outPath = os.path.join(exportPath, f"{figName}_remeshing_{suffix}.png")
        try:
            if _use_decorator:
                Decorator.savefig(outPath, pad=0.0, transparent=True)
                Decorator.closeAll()
            else:
                plt.savefig(outPath, dpi=150, bbox_inches='tight', transparent=True)
                plt.close(fig)
            print(f"  [visualizeRemeshing] → {outPath}")
        except Exception as e:
            print(f"  [visualizeRemeshing] ⚠️ Impossible de sauvegarder '{suffix}': {e}")
            if _use_decorator:
                Decorator.closeAll()
            else:
                plt.close('all')

    #  Génération des deux figures 
    draw_and_save(xlimGlobal, ylimGlobal, "global")
    draw_and_save(xlimZoom, ylimZoom, "zoom")

def checkBL(ext_outer):

    print("   -> Test du contour extérieur de la couche limite...")
    try:
        test_contour = G.T3mesher2D(ext_outer, triangulateOnly=1)
        print("   ✅ Contour valide, aucune queue d'aronde !")
        return True 
    except Exception as e:
        print("   ❌ Échec : La bordure extérieure se croise (Queue d'aronde détectée).")
        return False # Le test a échoué
    

def getMeshQuality (report, sweeps, pf, vf, rf):

    failures = report.get('failures', {})
    
    if 'vmin' in failures:
        vf = max(factorThreshold['minF'], vf * 0.9)   # ← moins de points dans les creux
        rf = max(factorThreshold['minF'], rf * 0.9)
        # pf = pf * (1 - 0.1) 
        # if sweeps < 1000 :
        #     sweeps = int(sweeps * (1.10))

    if 'rmax' in failures:
        vf = min(factorThreshold['maxF'], vf * 1.05)   # ← plus de points dans les creux
        # rf = rf * (1 + 0.05)
        pf = min(factorThreshold['maxF'], pf * 1.05)    # ← plus de points sur les pics
    
    if 'acrit' in failures:
        pf = min(factorThreshold['maxF'], pf * 1.05)
        vf = min(factorThreshold['maxF'], vf * 1.05)

    return sweeps, pf, vf, rf

def getGeometryDeviation (report, sweeps, pf, vf, rf):

    failures = report.get('failures', {})

    if 'maxInward' in failures:
        vf = min(factorThreshold['maxF'], vf * 1.05)   # ← moins de points dans les creux
        rf = min(factorThreshold['maxF'], rf * 1.05)
        pf = min(factorThreshold['maxF'], pf * 1.05)
        if sweeps < 1000 :
            sweeps = int(sweeps * (0.90))
    
    if 'stdDev' in failures:
        # Dispersion hétérogène → raffiner la résolution du profil
        pf = min(factorThreshold['maxF'], pf * 1.05)

    return sweeps, pf, vf, rf

def plotHistoryBL(history, exportPath="photos", figName="optimizationHistoryBL.png"):
    import os
    import matplotlib.pyplot as plt
    import numpy as np

    os.makedirs(exportPath, exist_ok=True)
    
    iters = np.array(history['iter'])
    
    # Création d'une figure avec 2 sous-graphes empilés qui partagent l'axe X
    fig, axs = plt.subplots(2, 1, figsize=(10, 12), sharex=True)
    fig.suptitle("Évolution des Paramètres et Critères au fil des Itérations", fontsize=16)

    # ==========================================
    # ZONE 1 : Facteurs (gauche) + Sweeps & Points (droite)
    # ==========================================
    line_pf = axs[0].plot(iters, history['pf'], label='Peak Factor (pf)', marker='o', color='green', linewidth=2)
    line_vf = axs[0].plot(iters, history['vf'], label='Valley Factor (vf)', marker='s', color='blue', linewidth=2)
    line_rf = axs[0].plot(iters, history['rf'], label='Rough Factor (rf)', marker='^', color='orange', linewidth=2)

    ax0_twin1 = axs[0].twinx()
    ax0_twin2 = axs[0].twinx()
    ax0_twin2.spines['right'].set_position(('outward', 60))
    line_sw = ax0_twin1.plot(iters, history['sweeps'], label='Number of Sweeps', marker='v', color='purple', linewidth=2, linestyle='--')
    line_pts = ax0_twin2.plot(iters, history['npts'], label='Nb Points Profil', marker='d', color='brown', linewidth=2, linestyle='--')
    
    axs[0].set_ylabel('Valeur des Facteurs')
    ax0_twin1.set_ylabel('Sweeps', color='purple')
    ax0_twin2.set_ylabel('Nb Points', color='brown')
    ax0_twin1.tick_params(axis='y', colors='purple')
    ax0_twin2.tick_params(axis='y', colors='brown')
    
    axs[0].grid(True, linestyle=':', alpha=0.7)
    
    # Regrouper les légendes du graphe du haut
    lines_top = line_pf + line_vf + line_rf + line_sw + line_pts
    labels_top = [l.get_label() for l in lines_top]
    axs[0].legend(lines_top, labels_top, loc='upper left')

    # ==========================================
    # ZONE 2 : Couches BL (gauche) et Vmin (droite)
    # ==========================================
    ax3 = axs[1]
    ax3_twin = ax3.twinx()

    # On filtre les itérations où le maillage a échoué (NaN)
    valid_idx = ~np.isnan(history['nk'])
    valid_iters = iters[valid_idx]
    valid_bl = np.array(history['nk'])[valid_idx]
    valid_vmin = np.array(history['vminBL'])[valid_idx]

    if len(valid_iters) > 0:
        # Axe gauche : Couches BL créées
        line1 = ax3.plot(valid_iters, valid_bl, label='Couches BL extrudées', marker='x', color='red', linewidth=2)
        ax3.set_ylabel('Nombre de Couches BL', color='red')
        ax3.tick_params(axis='y', labelcolor='red')
        ax3.axhline(y=history['nLayers'][0], color='red', linestyle='--', alpha=0.5, label='nLayers target')
        

        # Axe droite : Vmin (Échelle Logarithmique)
        line2 = ax3_twin.plot(valid_iters, valid_vmin, label='Vmin BL', marker='.', color='black', linewidth=2)
        ax3_twin.set_yscale('log')
        ax3_twin.set_ylabel('Vmin (Échelle Log)', color='black')
        ax3_twin.tick_params(axis='y', labelcolor='black')
        ax3_twin.axhline(y=DEFAULT_MESH_QUALITY['minVBL'], color='black', linestyle='--', alpha=0.5, label='Tolérance vminBL')
        print (DEFAULT_MESH_QUALITY['minVBL'])
        
        # Regrouper les légendes du graphe du bas
        lines_bot = line1 + line2
        labels_bot = [l.get_label() for l in lines_bot]
        ax3.legend(lines_bot, labels_bot, loc='upper left')
    
    ax3.set_xlabel('Itérations')
    ax3.grid(True, linestyle=':', alpha=0.7)

    # Ajustement et sauvegarde
    plt.tight_layout()
    plt.subplots_adjust(top=0.95) 
    
    outPath = os.path.join(exportPath, figName)
    plt.savefig(outPath, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"  [Graphe] 📈 Historique BL d'optimisation généré : {outPath}")

def plotHistory(history, exportPath="photos", figName="optimizationHistory.png"):

    os.makedirs(exportPath, exist_ok=True)
    
    iters = np.array(history['iter'])
    
    # Création d'une figure avec 3 sous-graphes empilés qui partagent l'axe X (Itérations)
    fig, axs = plt.subplots(2, 1, figsize=(10, 12), sharex=True)
    fig.suptitle("Évolution des Paramètres et Critères au fil des Itérations", fontsize=16)

    # ==========================================
    # ZONE 1 : Les Facteurs (PF, VF, RF)
    # ==========================================
    axs[0].plot(iters, history['pf'], label='Peak Factor (pf)', marker='o', color='green', linewidth=2)
    axs[0].plot(iters, history['vf'], label='Valley Factor (vf)', marker='s', color='blue', linewidth=2)
    axs[0].plot(iters, history['rf'], label='Rough Factor (rf)', marker='^', color='orange', linewidth=2)

    ax0_twin = axs[0].twinx()
    ax0_twin.plot(iters, history['sweeps'], label='Number of Sweeps', marker='v', color='purple', linewidth=2)
    
    axs[0].set_ylabel('Valeur des Facteurs')
    ax0_twin.set_ylabel('Sweeps', color='purple')
    ax0_twin.tick_params(axis='y', labelcolor='purple')
    axs[0].grid(True, linestyle=':', alpha=0.7)
    axs[0].legend(loc='upper left')

    # ==========================================
    # ZONE 2 : Les Critères de Maillage (Rmax et Vmin)
    # ==========================================
    ax3 = axs[1]
    ax3_twin = ax3.twinx() # Création d'un axe secondaire à droite pour Vmin (très petit)

    valid_idx = ~np.isnan(history['rmaxT3'])
    valid_iters = iters[valid_idx]
    valid_rmax = np.array(history['rmaxT3'])[valid_idx]
    valid_vmin = np.array(history['vminT3'])[valid_idx]

    if len(valid_iters) > 0:
    # Rmax (Axe de gauche, linéaire)
        line1 = ax3.plot(valid_iters, valid_rmax, label='Rmax T3', marker='x', color='red', linewidth=2)
        line2 = ax3_twin.plot(valid_iters, valid_vmin, label='Vmin T3', marker='.', color='black', linewidth=2)
        
        ax3.set_ylabel('Rmax', color='red')
        ax3.tick_params(axis='y', labelcolor='red')
        ax3.axhline(y=DEFAULT_MESH_QUALITY['maxRT3'], color='red', linestyle='--', alpha=0.5, label='Tolérance RmaxT3') # Ligne de tolérance

        # Vmin (Axe de droite, Logarithmique car très petit)
        ax3_twin.set_yscale('log')
        ax3_twin.set_ylabel('Vmin (Échelle Log)', color='black')
        ax3_twin.tick_params(axis='y', labelcolor='black')
        ax3_twin.axhline(y=DEFAULT_MESH_QUALITY['minVT3'], color='black', linestyle='--', alpha=0.5, label='Tolérance vminT3') # Ligne de tolérance
        
        # Rassembler les légendes des deux axes
        lines = line1 + line2
        labels = [l.get_label() for l in lines]
        ax3.legend(lines, labels, loc='upper left')
    
    ax3.set_xlabel('Itérations')
    ax3.grid(True, linestyle=':', alpha=0.7)

    # Ajustement et sauvegarde
    plt.tight_layout()
    plt.subplots_adjust(top=0.95) # Laisse de la place pour le titre principal
    
    outPath = os.path.join(exportPath, figName)
    plt.savefig(outPath, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"  [Graphe] 📈 Historique d'optimisation généré : {outPath}")


def exportHistory(history, bestReport=None, reportDev=None, finalAreaError=None, devThresholds = None, tol = None, exportPath="photos", fileName="optimizationHistory.txt"):
    """Exporte l'historique d'optimisation dans un fichier texte bien formaté,
       incluant les derniers rapports de qualité à la fin."""
    import os
    import numpy as np
    os.makedirs(exportPath, exist_ok=True)
    outPath = os.path.join(exportPath, fileName)

    with open(outPath, 'w', encoding='utf-8') as f:
        # En-tête du tableau
        f.write(f"{'Iter':>5} | {'PF':>6} | {'VF':>6} | {'RF':>6} | {'Sweeps':>6} | {'RmaxT3':>9} | {'VminT3':>10}\n")
        f.write("-" * 75 + "\n")

        # Remplissage ligne par ligne
        for i in range(len(history['iter'])):
            it = history['iter'][i]
            pf = history['pf'][i]
            vf = history['vf'][i]
            rf = history['rf'][i]
            sw = history['sweeps'][i]
            
            rmaxT3 = f"{history['rmaxT3'][i]:9.2f}" if not np.isnan(history['rmaxT3'][i]) else "      NaN"
            vminT3 = f"{history['vminT3'][i]:10.2e}" if not np.isnan(history['vminT3'][i]) else "       NaN"

            f.write(f"{it:5d} | {pf:6.3f} | {vf:6.3f} | {rf:6.3f} | {sw:6d} | {rmaxT3} | {vminT3}\n")

        # =========================================================================
        # Ajout des rapports finaux (Format Console Exact)
        # =========================================================================
        f.write("\n\n" + "=" * 75 + "\n")
        f.write("BILAN FINAL DU DERNIER MAILLAGE (BEST MESH)\n")
        f.write("=" * 75 + "\n")

        if devThresholds is None:
            devThresholds = defaultGeoQuality()

        if finalAreaError is not None:
            f.write(f"\n[Conservation] Area error : {finalAreaError:.8f} %\n")

        if reportDev is not None:
            i_max = reportDev.get('maxInward', 0)
            o_max = reportDev.get('maxOutward', 0)
            std = reportDev.get('stdDev', 0)
            mean = reportDev.get('meanDev', 0)

            f.write(f"\n  Geometry Deviation Check (Smoothing vs Rough Ice):\n")
            f.write(f"      Inward : max = {abs(i_max):>+7.4f} {'✓' if abs(i_max)<=devThresholds['maxInward'] else '✗'} (tol <= {devThresholds['maxInward']:>+7.4f})\n")
            f.write(f"      Outward : max = {o_max:>+7.4f} {'✓' if o_max<=devThresholds['maxOutward'] else '✗'} (tol <= {devThresholds['maxOutward']:>+7.4f})\n")
            f.write(f"      Dispersion : std = {std:>7.4f} {'✓' if std<=devThresholds['stdDev'] else '✗'} (tol <= {devThresholds['stdDev']:>7.4f})\n")
            f.write(f"      Mean deviation : mean = {mean:>7.4f} {'✓' if abs(mean)<=devThresholds['meanDev'] else '✗'} (tol <= {devThresholds['meanDev']:>7.4f})\n")

        if bestReport is not None:
            if tol is None:
                tol = DEFAULT_MESH_QUALITY
            
            # Couche Limite (BL)
            if 'vminBL' in bestReport:
                f.write(f"\n  BL Mesh Quality Check:\n")
                f.write(f"      Volume : vcrit = {bestReport['vcritBL']:<8} {'✓' if bestReport['vcritBL']<=tol['maxVcrit'] else '✗'} (tol <= {tol['maxVcrit']}) | vmin = {bestReport['vminBL']:.2e} {'✓' if bestReport['vminBL']>=tol['minVBL'] else '✗'} (tol >= {tol['minVBL']}) | vmax = {bestReport['vmaxBL']:.2e} | vmean = {bestReport['vmeanBL']:.2e}\n")
                f.write(f"      Angle : acrit = {bestReport['acritBL']:<8} {'✓' if bestReport['acritBL']<=tol['maxAcrit'] else '✗'} (tol <= {tol['maxAcrit']}) | amin = {bestReport['aminBL']:>5.1f}° | amax = {bestReport['amaxBL']:>5.1f}° | amean = {bestReport['ameanBL']:>5.1f}°\n")
                f.write(f"      Regularity : rcrit = {bestReport['rcritBL']:<8} (info)       | rmin = {bestReport['rminBL']:.2e} | rmax = {bestReport['rmaxBL']:.2e} {'✓' if bestReport['rmaxBL']<=tol['maxRBL'] else '✗'} (tol <= {tol['maxRBL']})| rmean = {bestReport['rmeanBL']:.2e}\n")

            # Non Structuré (T3)
            if 'vminT3' in bestReport:
                f.write(f"\n  Mesh Quality Check:\n")
                f.write(f"    --- Unstructured (T3) ---\n")
                f.write(f"      Volume : vcrit = {bestReport['vcritT3']:<8} {'✓' if bestReport['vcritT3']<=tol['maxVcrit'] else '✗'} (tol <= {tol['maxVcrit']}) | vmin = {bestReport['vminT3']:.2e} {'✓' if bestReport['vminT3']>=tol['minVT3'] else '✗'} (tol >= {tol['minVT3']}) | vmax = {bestReport['vmaxT3']:.2e} | vmean = {bestReport['vmeanT3']:.2e}\n")
                f.write(f"      Angle : acrit = {bestReport['acritT3']:<8} {'✓' if bestReport['acritT3']<=tol['maxAcrit'] else '✗'} (tol <= {tol['maxAcrit']}) | amin = {bestReport['aminT3']:>5.1f}° | amax = {bestReport['amaxT3']:>5.1f}° | amean = {bestReport['ameanT3']:>5.1f}°\n")
                f.write(f"      Regularity : rcrit = {bestReport['rcritT3']:<8} (info)       | rmin = {bestReport['rminT3']:.2e} | rmax = {bestReport['rmaxT3']:.2e} {'✓' if bestReport['rmaxT3']<=tol['maxRT3'] else '✗'} (tol <= {tol['maxRT3']}) | rmean = {bestReport['rmeanT3']:.2e}\n")

    print(f"  [Export] 💾 Historique et diagnostics enregistrés : {outPath}")