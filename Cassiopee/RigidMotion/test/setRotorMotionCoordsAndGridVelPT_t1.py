# - setRotorMotionAndGridVel (pyTree) -
# Rotor motion
import RigidMotion.PyTree as R
import Generator.PyTree as G
import KCore.test as test

time0 = 0.01
a = G.cart((0.2,-0.075,0), (0.01,0.01,0.1), (131,11,1))
# Mettre tous les parametres
RotorMotion={'Motion_Blade1':{'initial_angles' : [0.,0],#PSI0,PSI0_b
                              'alp0': -12.013,'alp_pnt' : [0.,0.,0.], 'alp_vct':[0.,1.,0.],
                              'rot_pnt' : [0.,0.,0.],'rot_vct':[0.,0.,1.],'rot_omg':104.71,
                              'span_vct' : [1.,0.,0.],
                              'pre_lag_pnt' : [0.075,0.,0.],'pre_lag_vct' : [0.,0.,1.],'pre_lag_ang' : -4.,
                              'pre_con_pnt' : [0.,0.,0.],'pre_con_vct' : [0.,1.,0.],'pre_con_ang' : 0.,
                              'del_pnt' : [0.075,0.,0.],'del_vct' : [0.,0.,1.],'del0' : -0.34190,
                              'del1c' : 0.48992E-01 , 'del1s': -0.95018E-01,
                              'bet_pnt' : [0.076,0.,0.],'bet_vct' : [0.,1.,0.],'bet0' : -2.0890,
                              'bet1c' : 3.4534, 'bet1s' : 0.0,
                              'tet_pnt' : [0.156,0.,0.],'tet_vct' : [1.,0.,0.],'tet0' : 12.807,
                              'tet1c' : 1.5450, 'tet1s' : -3.4534}}

dictBlade = RotorMotion["Motion_Blade1"]
init_angles = dictBlade["initial_angles"]
psi0 = init_angles[0]; psi0_b = init_angles[1]
transl_speed = (-87.9592,0.,0.)
alp_pnt = dictBlade["alp_pnt"]
alp_vct = dictBlade["alp_vct"]
alp0 = dictBlade["alp0"]
rot_pnt = dictBlade["rot_pnt"]
rot_vct = dictBlade["rot_vct"]
rot_omg = dictBlade["rot_omg"]
del_pnt = dictBlade["del_pnt"]
del_vct = dictBlade["del_vct"]
del0 = dictBlade["del0"]
delc = (dictBlade["del1c"],)
dels = (dictBlade["del1s"],)
bet_pnt = dictBlade["bet_pnt"]
bet_vct = dictBlade["bet_vct"]
bet0 = dictBlade["bet0"]
betc = (dictBlade["bet1c"],)
bets = (dictBlade["bet1s"],)
tet_pnt = dictBlade["tet_pnt"]
tet_vct = dictBlade["tet_vct"]
tet0 = dictBlade["tet0"]
tetc = (dictBlade["tet1c"],)
tets = (dictBlade["tet1s"],)
span_vct = dictBlade['span_vct']
pre_lag_pnt = dictBlade["pre_lag_pnt"]
pre_lag_vct = dictBlade["pre_lag_vct"]
pre_lag_ang = dictBlade["pre_lag_ang"]
pre_con_pnt = dictBlade["pre_con_pnt"]
pre_con_vct = dictBlade["pre_con_vct"]
pre_con_ang = dictBlade["pre_con_ang"]
R._setPrescribedMotion2(a, 'Motion_Blade1', transl_speed=transl_speed,
                        psi0=psi0, psi0_b=psi0_b,
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

R._setRotorMotionCoordsAndGridVel(a, time=time0)
test.testT(a)
