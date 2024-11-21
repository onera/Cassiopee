"""Fluid state variables for different reference dimensions."""

import math

#==============================================================================
# [1] Retourne un etat de reference adimensionne correspondant a un
# adimensionnement par la densite, la vitesse du son et la temperature
# Les grandeurs retournees sont des grandeurs adimensionnees
# Il faut utiliser un maillage adimensionne
# IN: MInf: Mach infini de l'ecoulement incident
# IN: alphaZ: angle de l'ecoulement incident par rapport a z
# IN: alphaY: angle de l'ecoulement incident par rapport a y
# IN: ReInf: Reynolds infini
# IN: MutSMuInf: ratio mut/mu a l'infini (0.2 pour les ecoulements externes)
# jusqu'a 100 pour les ecoulements internes
# IN: TurbLevelInf: niveau de turbulence a l'infini
# IN: si Mtip est different de None, c'est Mtip qui sert a definir la turbulence
# (pour k-w uniquement)
#==============================================================================
def adim1(MInf=0.5, alphaZ=0., alphaY=0., ReInf=1.e8, MutSMuInf=0.2,
          TurbLevelInf=1.e-4, Mtip=None):
    """Return a state corresponding to adimensioning by density, sound speed and temperature."""
    Gamma    = 1.4   # constante des gaz parfait (sans dimension)
    DRgp     = 287.053 # R gaz parfait (dimensionne)
    DCs      = 110.4 # Cs air (dimensionne)
    DTs      = 288.15 # Temperature de reference dans Sutherland (dimensionne)
    DTInf    = DTs # On choisit que l'etat de reference de sutherland
                   # soit la temperature a l'infini
    DMus     = 1.78938e-5 # Mu dimensionne a la temperature de reference
    DMuInf   = DMus  # toujours identification des grandeurs de reference et infini

    Cs = DCs / DTInf # Cs adimensionne
    Ts = DTs / DTInf # Ts adimensionne, correspond alors a mus=1.78938e-5

    RoInf   = 1. # adimensionne
    TInf    = 1. # adimensionne
    AInf    = 1. # adimensionne

    UInf = MInf

    alz = alphaZ * math.pi / 180.
    aly = alphaY * math.pi / 180.

    PInf    = 1. / (Gamma)
    RouInf1 = RoInf * UInf * math.cos(alz)
    RovInf1 = RoInf * UInf * math.sin(alz)
    RowInf1 = 0.
    RouInf  = RouInf1 * math.cos(aly)
    RovInf  = RovInf1
    RowInf  = RouInf1 * math.sin(aly)
    RoEInf  = PInf/(Gamma-1.) + 0.5*RoInf*UInf*UInf
    cvInf = (RoEInf/RoInf - 0.5*UInf*UInf)/TInf
    cpInf = Gamma*cvInf

    # si Mtip existe, les grandeurs visqueuses et turbulentes
    # sont prises par rapport a Mtip
    if Mtip is not None: UInf = Mtip*AInf

    MuInf = RoInf*UInf / max(ReInf, 1.e-10)
    Mus = MuInf # identification

    Pr = 0.70951 # Prandtl (sans dimension)

    # Pour k-omega
    RokInf = 1.5 * RoInf * (TurbLevelInf * UInf)**2
    MutInf = MutSMuInf * MuInf
    RoomegaInf = RoInf * RokInf / max(MutInf, 1.e-10)

    # Pour spalart
    RonutildeInf = MutInf

    return [RoInf, RouInf, RovInf, RowInf, RoEInf, PInf, TInf, cvInf, MInf,
            ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
            Mus, Cs, Ts, Pr]

#==============================================================================
# [2] Retourne un etat de reference adimensionne correspondant a un
# adimensionnement par la densite, la vitesse et la temperature
# Les grandeurs retournees sont des grandeurs adimensionnees
# IN: MInf: Mach infini de l'ecoulement incident
# IN: alphaZ: angle de l'ecoulement incident par rapport a z
# IN: alphaY: angle de l'ecoulement incident par rapport a y
# IN: ReInf: Reynolds infini
# IN: MutSMuInf: ratio mut/mu a l'infini (0.2 pour les ecoulements externes)
# jusqu'a 100 pour les ecoulements internes)
# IN: TurbLevelInf: taux de turbulence a l'infini en pourcentage de UInf
# (pour k-w uniquement)
# Cet adimensionnement ne doit pas etre utilise pour Minf=0.
#==============================================================================
def adim2(MInf=0.5, alphaZ=0., alphaY=0., ReInf=1.e8, MutSMuInf=0.2,
          TurbLevelInf=1.e-4):
    """Return a state corresponding to adimensioning by density, velocity and temperature."""
    Gamma    = 1.4  # constante des gaz parfait (sans dimension)
    DRgp     = 287.053 # R gaz parfait (dimensionne)
    DCs      = 110.4 # Cs air (dimensionne)
    DTs      = 288.15 # Temperature de reference dans Sutherland (dimensionne)
    DTInf    = DTs # On choisit que l'etat de reference de sutherland
                   # soit la temperature a l'infini
    DMus     = 1.78938e-5 # Mu dimensionne a la temperature de reference
    DMuInf   = DMus  # toujours identification des grandeurs de reference et infini

    Cs = DCs / DTInf # Cs adimensionne
    Ts = DTs / DTInf # Ts adimensionne, correspond alors a mus=1.78938e-5

    RoInf   = 1. # adimensionne
    TInf    = 1. # adimensionne
    UInf    = 1. # adimensionne

    AInf = UInf / MInf

    alz = alphaZ * math.pi / 180.
    aly = alphaY * math.pi / 180.

    PInf    = (AInf*AInf) / (Gamma)
    RouInf1 = RoInf * UInf * math.cos(alz)
    RovInf1 = RoInf * UInf * math.sin(alz)
    RowInf1 = 0.
    RouInf  = RouInf1 * math.cos(aly)
    RovInf  = RovInf1
    RowInf  = RouInf1 * math.sin(aly)
    RoEInf  = PInf/(Gamma-1.) + 0.5*RoInf*UInf*UInf
    cvInf = (RoEInf / RoInf - 0.5*UInf*UInf)/ TInf
    cpInf = Gamma*cvInf

    MuInf = RoInf*UInf / max(ReInf, 1.e-10)
    Mus = MuInf # identification

    Pr = 0.70951 # Prandtl (sans dimension)

    # Pour k-omega
    RokInf = 1.5 * RoInf * (TurbLevelInf * UInf)**2
    MuInf = RoInf*UInf / max(ReInf, 1.e-10)
    MutInf = MutSMuInf * MuInf
    RoomegaInf = RoInf * RokInf / max(MutInf, 1.e-10)

    # Pour spalart
    RonutildeInf = MutInf

    return [RoInf, RouInf, RovInf, RowInf, RoEInf, PInf, TInf, cvInf, MInf,
            ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
            Mus, Cs, Ts, Pr]

#==============================================================================
# [3] Retourne un etat de reference adimensionne correspondant a un
# adimensionnement par la densite, la vitesse du son et la temperature
# Les grandeurs retournees sont des grandeurs adimensionnees
# Il faut utiliser un maillage adimensionne
# IN: MInf: Mach infini de l'ecoulement incident
# IN: alphaZ: angle de l'ecoulement incident par rapport a z
# IN: alphaY: angle de l'ecoulement incident par rapport a y
# IN: ReInf: Reynolds infini
# IN: LInf: longeur de reference (en m), sert seulement a calculer la viscosite
# IN: MutSMuInf: ratio mut/mu a l'infini (0.2 pour les ecoulements externes)
# jusqu'a 100 pour les ecoulements internes
# IN: TurbLevelInf: taux de turbulence a l'infini en pourcentage de Uinf
# (pour k-w uniquement)
# IN: Mtip: si Mtip different de None, il sert pour la definition de la turbulence
#==============================================================================
def adim3(MInf=0.5, alphaZ=0., alphaY=0., ReInf=1.e8, LInf=1., MutSMuInf=0.2,
          TurbLevelInf=1.e-4, Mtip=None):
    """Return a state corresponding to adimensioning by density, sound speed and temperature."""
    Gamma    = 1.4   # constante des gaz parfait (sans dimension)
    DRgp     = 287.053 # R gaz parfait (dimensionne)
    DCs      = 110.4 # Cs air (dimensionne)
    DTs      = 288.15 # Temperature de reference dans Sutherland (dimensionne)
    DTInf    = DTs # On choisit que l'etat de reference de sutherland
                   # soit la temperature a l'infini
    DMus     = 1.78938e-5 # Mu dimensionne a la temperature de reference
    DMuInf   = DMus  # toujours identification des grandeurs de reference et infini

    Cs = DCs / DTInf # Cs adimensionne
    Ts = DTs / DTInf # Ts adimensionne, correspond alors a mus=1.78938e-5

    RoInf   = 1. # adimensionne
    TInf    = 1. # adimensionne
    AInf    = 1. # adimensionne

    UInf = MInf

    alz = alphaZ * math.pi / 180.
    aly = alphaY * math.pi / 180.

    PInf    = 1. / (Gamma)
    RouInf1 = RoInf * UInf * math.cos(alz)
    RovInf1 = RoInf * UInf * math.sin(alz)
    #RowInf1 = 0.
    RouInf  = RouInf1 * math.cos(aly)
    RovInf  = RovInf1
    RowInf  = RouInf1 * math.sin(aly)
    RoEInf  = PInf/(Gamma-1.) + 0.5*RoInf*UInf*UInf
    cvInf = (RoEInf/RoInf - 0.5*UInf*UInf)/TInf
    cpInf = Gamma*cvInf

    # si Mtip existe, les grandeurs visqueuses et turbulentes
    # sont prises par rapport a Mtip
    if Mtip is not None: UInf = Mtip*AInf

    MuInf = RoInf*UInf*LInf / max(ReInf, 1.e-10)
    Mus = MuInf # identification

    Pr = 0.70951 # Prandtl (sans dimension)

    # Pour k-omega
    RokInf = 1.5 * RoInf * (TurbLevelInf * UInf)**2
    MutInf = MutSMuInf * MuInf
    RoomegaInf = RoInf * RokInf / max(MutInf, 1.e-10)

    # Pour spalart
    RonutildeInf = MutInf

    return [RoInf, RouInf, RovInf, RowInf, RoEInf, PInf, TInf, cvInf, MInf,
            ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
            Mus, Cs, Ts, Pr]

#==============================================================================
# [4] Retourne un etat de reference dimensionne correspondant a
# l'air sec au niveau du sol, considere comme un gaz parfait
# Retourne des valeurs physiques dimensionnees (USI)
# IN: UInf: vitesse en m/s (default 10km/h)
# IN: TInf: temperature en K (0 degree=273.15K) (default 25 degres)
# IN: PInf: pression en Pa (defaut 1 atm)
# IN: LInf: longeur de reference (en m), sert seulement a calculer le Re
# IN: Mtip: si Mtip different de None, il sert pour la definition de la turbulence
#==============================================================================
def dim1(UInf=2.7777, TInf=298.15, PInf=101325., LInf=1.,
         alphaZ=0., alphaY=0., MutSMuInf=0.2,
         TurbLevelInf=1.e-4, Mtip=None, MuMultiplier=1.):
    """Return a dimensioned state specifying velocity, temperature and pressure."""
    alz = alphaZ * math.pi / 180.
    aly = alphaY * math.pi / 180.

    Gamma = 1.4 # constante des gaz parfait (sans dimension)
    Rgp = 287.053 # dimensionne

    # lois approx air
    RoInf = 1.292 * 273.15 / TInf # in kg/m3, loi air sec
    # ou loi des gaz parfaits
    RoInf = PInf/(Rgp*TInf)

    # mu en kg/m/s par loi directe
    MuInf = 8.8848e-15*TInf**3 -3.2398e-11*TInf**2 + 6.2657e-8*TInf+2.3543e-6
    # ou mu par Sutherland
    Cs      = 110.4      # dimensionne
    Ts      = 288.15     # temperature correspondante dimensionnee
    Mus     = 1.78938e-5 # mu correspondant dimensionne
    MuInf   = Mus*math.sqrt(TInf/Ts)*(1.+Cs/Ts)/(1.+Cs/TInf)
    MuInf = MuInf * MuMultiplier

    # nu en m2/s
    NuInf = MuInf / RoInf

    # diffusivite thermique en W/m/K
    lambdaInf = 1.5207e-11*TInf**3 -4.857e-8*TInf**2+1.0184e-4*TInf-3.9333e-4

    #print('RoInf=', RoInf)
    #print('TInf=', TInf)
    #print('MuInf=', MuInf)
    #print('NuInf=', NuInf)
    #print('Pinf=', PInf)
    #print('lambdaInf=', lambdaInf)

    # deduit
    RouInf1 = RoInf * UInf * math.cos(alz)
    RovInf1 = RoInf * UInf * math.sin(alz)
    #RowInf1 = 0.
    RouInf  = RouInf1 * math.cos(aly)
    RovInf  = RovInf1
    RowInf  = RouInf1 * math.sin(aly)

    RoeInf = PInf/(Gamma-1.)
    RoEInf = RoeInf+0.5*RoInf*UInf*UInf
    aInf = math.sqrt(Gamma*PInf/RoInf)
    MInf = UInf/aInf
    #print('aInf=', aInf)
    #print('MInf=', MInf)

    eInf = RoeInf/RoInf
    cvInf = eInf/TInf
    cpInf = Gamma*cvInf
    #print('cpInf=', cpInf)

    # si Mtip existe, les grandeurs visqueuses et turbulentes
    # sont prises par rapport a Mtip
    if Mtip is not None: UInf = Mtip*aInf

    ReInf = UInf*LInf/NuInf
    #print('ReInf=', ReInf)

    Pr = MuInf*cpInf/lambdaInf # Prandtl (sans dimension)
    #print('Pr=', Pr)

    # Cs
    #print('Cs=', Cs)

    # Pour k-omega
    RokInf = 1.5 * RoInf * (TurbLevelInf * UInf)**2
    MutInf = MutSMuInf * MuInf
    RoomegaInf = RoInf * RokInf / max(MutInf, 1.e-10)

    # Pour spalart
    RonutildeInf = MutInf

    return [RoInf, RouInf, RovInf, RowInf, RoEInf, PInf, TInf, cvInf, MInf,
            ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
            Mus, Cs, Ts, Pr]

#==============================================================================
# [5] Retourne un etat de reference dimmensionne correspondant a
# l'air sec au niveau du sol, considere comme un gaz parfait
# Different input than dim1
#==============================================================================
def dim2(UInf=2.7777, TInf=298.15, RoInf=1.225, LInf=1.,
         alphaZ=0., alphaY=0., MutSMuInf=0.2,
         TurbLevelInf=1.e-4, Mtip=None, MuMultiplier=1):
    """Return a dimensioned state specifying velocity, temperature and density."""
    Rgp = 287.053
    PInf = RoInf*Rgp*TInf
    return dim1(UInf=UInf, TInf=TInf, PInf=PInf, LInf=LInf,
                alphaZ=alphaZ, alphaY=alphaY, MutSMuInf=MutSMuInf,
                TurbLevelInf=TurbLevelInf, Mtip=Mtip, MuMultiplier=MuMultiplier)

#==============================================================================
# [6] Retourne un etat de reference dimmensionne correspondant a
# l'air sec au niveau du sol, considere comme un gaz parfait
# Different input
#==============================================================================
def dim3(UInf=2.7777, PInf=101325., RoInf=1.2, LInf=1.,
         alphaZ=0., alphaY=0., MutSMuInf=0.2,
         TurbLevelInf=1.e-4, Mtip=None, MuMultiplier=1):
    """Return a dimensioned state specifying velocity, pressure and density."""
    Rgp = 287.053
    TInf = PInf/(RoInf*Rgp)
    return dim1(UInf=UInf, TInf=TInf, PInf=PInf, LInf=LInf,
                alphaZ=alphaZ, alphaY=alphaY, MutSMuInf=MutSMuInf,
                TurbLevelInf=TurbLevelInf, Mtip=Mtip, MuMultiplier=MuMultiplier)

#==============================================================================
# [7] Retourne un etat de reference dimensionne correspondant a
# l'air sec au niveau du sol, considere comme un gaz parfait
# Retourne des valeurs physiques dimensionnees (USI)
# IN: UInf: vitesse en m/s (default 10km/h)
# IN: TInf: temperature en K (0 degree=273.15K) (default 25 degres)
# IN: PInf: pression en Pa (defaut 1 atm)
# IN: LInf: longeur de reference (en m), sert seulement a calculer le Re
# IN: Mus: viscosite de Sutherland, sert seulement a modifier le Re
#==============================================================================
def dim4(UInf=2.7777, TInf=298.15, PInf=101325., LInf=1.,
         alphaZ=0., alphaY=0., Mus=1.78938e-5, MutSMuInf=0.2,
         TurbLevelInf=1.e-4, Mtip=None):
    """Dimensional state 4."""
    alz = alphaZ * math.pi / 180.
    aly = alphaY * math.pi / 180.

    Gamma = 1.4 # constante des gaz parfait (sans dimension)
    Rgp = 287.053 # dimensionne

    # lois approx air
    RoInf = 1.292 * 273.15 / TInf # in kg/m3, loi air sec
    # ou loi des gaz parfaits
    RoInf = PInf/(Rgp*TInf)

    # mu en kg/m/s par loi directe
    MuInf = 8.8848e-15*TInf**3 -3.2398e-11*TInf**2 + 6.2657e-8*TInf+2.3543e-6
    # ou mu par Sutherland
    Cs      = 110.4 # dimensionne
    Ts      = 288.15 # temperature correspondante dimensionnee
    #Mus     = 1.78938e-5 # mu correspondant dimensionne
    MuInf = Mus*math.sqrt(TInf/Ts)*(1.+Cs/Ts)/(1.+Cs/TInf)

    # nu en m2/s
    NuInf = MuInf / RoInf

    # diffusivite thermique en W/m/K
    lambdaInf = 1.5207e-11*TInf**3 -4.857e-8*TInf**2+1.0184e-4*TInf-3.9333e-4

    # deduit
    RouInf1 = RoInf * UInf * math.cos(alz)
    RovInf1 = RoInf * UInf * math.sin(alz)
    #RowInf1 = 0.
    RouInf  = RouInf1 * math.cos(aly)
    RovInf  = RovInf1
    RowInf  = RouInf1 * math.sin(aly)

    RoeInf = PInf/(Gamma-1.)
    RoEInf = RoeInf+0.5*RoInf*UInf*UInf
    aInf = math.sqrt(Gamma*PInf/RoInf)
    MInf = UInf/aInf

    eInf = RoeInf/RoInf
    cvInf = eInf/TInf
    cpInf = Gamma*cvInf

    # si Mtip existe, les grandeurs visqueuses et turbulentes
    # sont prises par rapport a Mtip
    if Mtip is not None: UInf = Mtip*aInf 

    ReInf = UInf*LInf/NuInf

    Pr = MuInf*cpInf/lambdaInf # Prandtl (sans dimension)

    # Cs

    # Pour k-omega
    RokInf = 1.5 * RoInf * (TurbLevelInf * UInf)**2
    MutInf = MutSMuInf * MuInf
    RoomegaInf = RoInf * RokInf / max(MutInf, 1.e-10)

    # Pour spalart
    RonutildeInf = MutInf

    return [RoInf, RouInf, RovInf, RowInf, RoEInf, PInf, TInf, cvInf, MInf,
            ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
            Mus, Cs, Ts, Pr]
