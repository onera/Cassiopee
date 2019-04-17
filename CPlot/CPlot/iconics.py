# iconics Photos
try: import Tkinter as TK
except: import tkinter as TK
PHOTO = []

saveImg = TK.PhotoImage(data="""
R0lGODlhFAAUAOfPAAEBAQcJCAoFCQsLCw4PEw8SCw4SExQMEhYUDhkSDhMTExQWGRYYFBUZGhkV
EhoVHBoYFhsbGxccIxcfKh4eIBcjIRYjKR4gIhggKxogMCYVDSAXEiEcGiEgHCUlJSYmKCcoKigo
JioqKiUuNS0uMC8wMjMzMzY2ODQ9NTA5Pjs5NDk5OQQzXQcsYQUuawMvcAU1YgAzbQY4Ywc8bgk1
Ygo4YwA1cgA3egA5dAA9egk7eD4+QABBfgxEfxxMeilEYUFAPEJCQkdIS0pKS0xNUU9QU1FSVFRV
WFdYWlpaXF5eYF5gZGJiZGZmaGtra3Nzc3N0eXt7ewE+gABEgwBGiABIhQBKiwpHhgBOkABQjwBS
lABVmQBYlwBZnRtChQBcogBeqABgpQBjqgBorgBmsQBqsgFtuCVcjSZekSlakitcmDFVhS1jjyVh
lCZlmSdpnS5hkitnmixqnDBllzBlmSdtoSttoSdyqSxxpSt1qyh5ril4sElpmE92n1x2kX5/gWl/
pHeChoSEhIyNjY6SlZOTlJKWmZeYnJ2enpqco56uvZ+xu6KenaOjo6OlrKWor62traapsKqts62w
tq6xuLGxsbO1urW4vrq6vZOswra5wLu9wr7Axb7ByL7JwcHCxcTFycXIxsbIy8jFw8jGy8nIxsrK
zcbJ0M3N0c7Q09DMzdDO09DQztLS1dXV2dbY29jW1djX29zY19ra3Mfa4N3d4N7g39/g4+Pd1ODd
3ODe4+Dg3uLi5Obm6OXo5efo6ejm6ejq5+vr7O7u8O7w7u7x8PHu7fLr8PLx7fPz8/b2+Pb49/f6
+f379f7+/v//////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////yH5BAEKAP8ALAAAAAAUABQA
AAj+AP8JHEiwoMETJkSIMDGr4SxeEGvZUrHChCCBJjyE+PBhVayPHjiGqGUioQmMFDp48JDK1ccI
FFbGUqhQoAgPFzxQMNWqZUwPHVrRFCHQQ4SjFDadctUK5lFXHD2AEPghJIiNIT+AWMl15YeiHKpy
DXGhKlasHgSGAMC2rdu3AIj+C+GMWTNld5UpS8Y3mTBhyZyl/SfCmTJme/kK+yUM2K9fvXo5kysC
sV7FwoIBi1wLIjPKd+8OQyZsWGSIqJNRtuwX2bBgqGfVapjs5D8TfpMN+8WrF69Zr5iyYpUKmG0T
xZKR9s2r1qtWrVaZMqWKl+0TuYXFhtWKFapQpkoqWReI/e8vXbNndZ8OKlSpWbZVxMCBg4eWLV3C
lNkvposWLVTYZtCABgUEADs=
""")

undoImg = TK.PhotoImage(data="""
R0lGODlhFAAUAMZrACZj0zNjuyln1Spo1ixp1TBp1Cxr1y5r0zdo0TZp0TZp0jZq0jZq0zxpyzZr
0zZr1DtqzTdr1DZs1DVt0jds1Ddt1C9w2jdt1Txuwzht1Dhuzzht1Tdu1Tdu1jhu1Dhu1TZv1kVs
vzhv1Tlv1Tlv1jhw1jRy2Dlw1jlw10Juyj9vzUBvzD1wzzpx0kBvzUFvyzpx1Dlx10FvzDpx1kJv
yz1w1Dpx1zxyzjpy1jpy1zty1z1y0z5y0ztz1ztz2Dxz2EFy0EBz0zx02UBz1Tx12UB01T122T52
1z122kB30EN10z931EZ1zj542kV20kd2zT561D552z963EJ71T993EF93j1/4UF+3jyB4j2B4kKB
4ESB4ECH5U6HykKI5kuK1UaN6EmO6E2S5UyS6U2U61Wa4lGa7Vqh7Vej816l71up9f//////////
/////////////////////////////////////////////////////////////////////////yH+
FUNyZWF0ZWQgd2l0aCBUaGUgR0lNUAAh+QQBCgB/ACwAAAAAFAAUAAAHsIB/goOEhYaHTEBBRTws
hiGFT0MeIyczOT0/O4QpDYNOGR8kMzo+QkRGSEeDNAgQf0oVojY+PRotTVFSUkuCLwkKNRQbJDY4
AYQ3VVVXGH8yCgwRxDEwh0lbW1N/KwsOFCIoOIeCVFZa3AwSHCU5E+R/UFhZfy4OFx0gBAMH8F5c
f1Q8uFAAgAADFkyQCwMGHrwxZBySMyNG4iE0XywaUqOx0JkyHQmlCUmoC8mTDgMBADs=
""")

deleteImg = TK.PhotoImage(data="""
R0lGODlhFAAUAOeEALUFBbcGBrUHB7cHB7cICLkICLgJCboJCboKCrsKCrwKCrwLC7sMDL0MDL4M
DL8NDcAODrcSEsAPD8EPD8AQEMEQEMIQEMAREbwTE8IREbQXF8MREbMYGL8UFLYYGLMaGrQaGsUT
E8YTE7gaGrocHMgWFskWFsAaGrodHbseHrweHswZGb0gIMEiIs8cHMYhIcIjI8AlJdMfH8wmJtYi
IskoKNcjI8QtLcgtLbg0NMYuLrszM8cuLr4yMrg1NcYvL9smJr0zM8UwMMYwMNIrK88tLcYxMcQy
Mrk3N8UyMsUzM8Q0NN4pKcU0NMQ1Nd4qKsE3N8Q2NsI3N8A4ON8rK8Q3N8Y4OMI6OsA7O8Q6Osc5
OdUzM9kyMuMuLsU7O8A+PsU8PMM9PeQvL8Q9PcI+PsU+PuYxMcY/P8s9PecyMsJCQto4OMNCQsZB
QekzM9s6OuA4OMdDQ+03N+I8PMRKSu44OPA5OcdOTvM8PMpOTvVAQMxTU/ZBQc1UVMtVVc1VVcxW
VsxXV89kZNBmZv//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////yH+FUNyZWF0ZWQgd2l0aCBU
aGUgR0lNUAAh+QQBCgD/ACwAAAAAFAAUAAAI+QD/CRxIsKBBgYLoHBzoB0nBQXniqFkYKEoUHwMh
xmlT5otBQE2OHFGSQ+CdjWXAZMFC8E8SITCF7BjIJmWWKk6mCOxj5McPHTqCFCRz08mSI1D2DNHB
gweOHgevGD1yo8WJDhcoVECx8J+UGBgOLHAAocKGFF0FRhDrQEKGECzSCvSSoO2GECVgyD3DoAGE
DSJMrHBRo+sZBAseWMA7WAaNIgfLGDjgdwOLFy5k2ADCZEvBMQQK1J2gQuAMGkCodDGzZmAYAQEK
KHhAgiCRJ2LcyKnzRuAHALEXjDDIJY0cPHqsDOQwoICHhXDs8EFTEISGtHO0yN0eEAA7
""")

fitImg = TK.PhotoImage(data="""
R0lGODlhFAAVAOeSADo6Ojs7Ozw8PD09PT4+PkBAQEFBQUREREZGRklJSUpKSktLS0xMTP8AAE1N
TU9PT1BQUFFRUVJSUlRUVFVVVVZWVldXV1hYWFlZWVpaWltbW15eXl9fX2BgYGFhYWJiYmNjY2Rk
ZGhjY2VlZWpkZGZmZmdnZ21mZmlpaW9nZ2pqamxsbG1tbW5ubm9vb3BwcHFxcXJycnR0dHV1dXZ2
dnt7e3x8fH19fX5+foCAgIGBgYKCgoODg4SEhIWFhYaGhouLi+Bqat5ra9xsbOBsbN9tbedqapCQ
kJOTk5SUlJaWlpeXl5mZmZqampubm5+fn6GhoaKioqSkpKWlpaampqenp6ioqKmpqaqqqqurq62t
ra6urq+vr7GxsbKysrOzs7a2trm5ubq6uru7u7y8vL29vb6+vsDAwMPDw8TExMbGxsjIyMrKysvL
y83Nzc7Ozs/Pz9DQ0NHR0dLR0dPT09TU1NXV1dbW1tfX19jY2Nra2tvb29zc3N3d3d7d3d7e3uDg
4OHh4eLi4uPj4+Tk5OXl5ebm5urq6uvr6+zs7O3t7e/v7/Dw8PLy8vT09PX19fv7+/7+/v//////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////yH5BAEKAP8ALAAAAAAUABUA
AAj+AP8JHEjISo4RI25IATSw4cBBNgQAACBAIoAAMvo4FFjmAIAVYhCJ8ZKIzAsABr44DDMgARpJ
d5iUALGkjqQ1DAJoGXinAINAimZMHAoABqJCDwbEETgigBtFFABMoNIkyRULACQcoiNgw782AG5I
igGghqMuM7c80vFRkg8AaYAAeFMHQIZHLYiieNQBgJs7AHR4MABJCQAvWACIGCLkBIAoYgAciYSg
goEFYUIAgNK3SIMGRgBooAKgQxgIAAgoAKP5CQcARD6HziClNJgHAC4QbgLgyxUAJIQESbF5DAAk
lCOwhaNH9KMVRE088gCgDh4ANs4AwCHJRthHXEZ0fMjyiAcAGpJ+ABjzD0MAOI02AKhwxYmSLBhE
M7IzwIJANgI4IIgjPVg0UQA5NGIIBAKsMVBtCqghCSFXvMBCFYJIwoYDmzlkBQEAuECGIiMtYgZZ
A0yx0T92qDAUASBOZMJSKw6ExAoWPJDBDmnU6KOPAQEAOw==
""")

mainImg = TK.PhotoImage(data="""
R0lGODlhFAAUAMZ8ABA5ihE9jwc/pwpEqhdEkh1HfSJGiR9Pkh9OnB5QmCFSiyNTgx9YfyNUmBlY
oxNasxNbsiNYmh1XsCBZoytXpxlfsBxfrxNjvS9boihhlh9iqypgmBxlqRtktRRmvixfryJrphpt
ujJlqDBnpyFxpTFntilssh10txp2xyRzvEZtgCR4rjl0mRp8xB17xCp7uC18tTt6nzR/rSyDtzGB
yzeLsTuKuieQ0TOQtjiPs0KLsjWRuSiU1ESMt0iIxzOVwUSQsz+UxFeSozGc1kmYw0yWxz2d0UCc
1k+ctUOey0mdwk6gtV2Ztlmcuz6n2ziq3Tqq2k+jzkin3Vejyl6juD+u3lSpx0au1z+w4E6r2FKs
0FerymKl00mx4Gmlyk+y0W6oumKt0laz2mysynerv3CvzVq64mq12W+003eyzGC831fB4nK6ymS/
5lvL6nvB4HnC53PG54rD1o3F343G33/P64nQ5IbR7o7S7ZTP7ojV74nW7P///////////////yH+
FUNyZWF0ZWQgd2l0aCBUaGUgR0lNUAAh+QQBCgB/ACwAAAAAFAAUAAAHtYB/goOEhX9phoRshXN5
ZYmFcmCCb3BhKpCDVHV4ZHRnTUxeXEWZQmptaGNTUkc9LIZ2e3p3cWJdWVFEOjGJbmtfSFZaSk5G
QTYvND6FS4NbWDlJQzMuKCkZHRUWhmZVOEA8KycgHB4XCoIHhFdPO383JIIhfwsQD4lQP38yLYMm
gwokMFRDEAwGgzQ4yJRpQwSGkEZMgJhIRAOKhkogwFhIAgGOgz5QMADyD4YBAkoKCgAAUiAAOw==
""")

eyeImg = TK.PhotoImage(data="""
R0lGODlhFAAUAOetAAAAAAEBAQICAgMDAwUFBgYGBgcHBwkJCQwMDAwMDRAQEBMTExUVFRUWGBYW
FhgYGBoaGhsbGxwcHBwcHR0dHR4eHx8fHwkmOiEhISIiIiMjIxApOSQkJCUlJScnJykpKR8zQS8v
LxM6VjMzMxE9Wx48UTY2Njk5OSRDVz09PUNDQ0VFRUZGRkdHR09PTxxeilFRURJlnldXVxlvqFxc
XBVxtl1dXRdythJ2uhB2wBJ1wGFhYht2sg95xidzpQh9zA98yBV7xGVlZB55vhJ+yCd7uhKCzmpq
aSZ/vCp+uSWAvyiAvSqAuCqAuSWBwmpucCqBuheG021ubSiCwDF/uSGFyCuCvhqH0m9vbzKBux2I
0W5xczaBuSGJ0CKJznJycjiCvSiLzj6EuimMziqMzXZ2dnh4eHl5eXt7e3x8fEiOxEmQxoKDg4OD
goODg4SEhFmWxF6Vv4qKin2PnGOYwWqbwZSUlJiYmJycnHypzoeuzqioqKqqqqurq6ysrK2tra+v
r5u2zpm61LS0tLi4uLy8vKnE26zG3MLCwsbGxsjIyL/O3LzQ4rzR48LU49LS0snW4sbX5tTU1NXV
1dfX19jY2NnZ2dzc3N3d3d/f3+Dg4OHh4eLi4uPj4+3t7fLy8vPz8/P2+vX4+/n5+Pr6+vv7+/z8
/Pz9/f39/f39/v3+/v7+/v7//v//////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////yH5BAEKAP8ALAAAAAAUABQA
AAj+AP8JHEiwoMGDCBMKtPTnzh0/lBT+62TGAoCLGCV8yXQwzYEAGFwc2fJEyAoGAQyUIbhphAk7
nVYdEnOjBhhDqCS98RAC079LNhCtatVKjxIyUJiMCZKHqClALSwlGkrU0RIkGwYUANGkByOirUjx
AUuUTpgLCjSc8DTHyxpWZMmqosIjwQcKJyzt8KHjVFywopLEiAABgAAWDF7kCPWXqCouOCYsQGAB
T4UZQ1KRBUW2zhURGR48cFBCCxy4RDVNegQWkpUfJBoQQEHESSSiqwZZ+ldIzqdWqAJNiWKkShcg
goZSKrNbYCUYaSqVWhQnSxE1jUYRomHmYJ8UHGQYnGnDRoqKDlh8JuS0x00ZNHcUSZxPX2JAADs=
""")

selectallImg = TK.PhotoImage(data="""
R0lGODlhFAAUAMZ8ADVHYU1dc2t4i215jHWBk4OOnpCbqZOdq5ynt56puKOrt6+6xbC6xrG7xrG8
x7O9yLS+ybW/yrrAybLE0LjCzbnDzrrDzrbH07jJ1LrL1r3N18PL1L7N2L7O2MHQ2sfO18LQ2sTS
28TS3MXS3MXT3MrR2szR18bT3cbU3c3T2c7T2s3U3snW38/U287V3srX38vX4MzY4NDX4dLX3c3Z
4c7Z4dHY4dPY3tPY39TY38/a4s/b4tTZ4NDb49Hc49bb4tTe5dbe59Xf5tje5Nbf5tre49bg5tbg
59fg59fh59rg5tng6djh6N3g5drh6drh6tni6N3h5tvi6trj6dvj6d/i5tvk6tzk6t3k6t3l6t7l
69/l7N7m69/m7ODm7eHm7eDn7OHn7eHo7eLp7ePp7uPq7uXp7+Tq7uXr7+br8Obs8Ofs8Oft8eju
8enu8uru8urv8uzw8+3x9O7y9O7y9fD09vH09vH09/T19/P2+PX3+ff5+v///////////////yH+
FUNyZWF0ZWQgd2l0aCBUaGUgR0lNUAAh+QQBCgB/ACwAAAAAFAAUAAAH/oABBYOEhYaEgnyKi4yN
i4N8ElWNJiaNVRJ8kAAKjQICjQoAmgV8eE0bHyUuNisyTlJPS0FmeKSKOVBcZGxxdnp7e3p1cWGK
kHw4TFliaXBzd3l5dnNwX8elfDdHVF1la3BydHRycGte2IozQFBYYGdqbW9ubWpnW+kDBjtCUFlg
Y86gQXNmDBgEA24ROPDCBxEoVrR0AQOmi5YrCQjc4tNiBIweQJJAmUKFyhQoSJSk49hBBAsaPoAI
MWJECBAfQ1aqwKABxAkWMWro0FEjBosfK1NMuJCBg4cQJFCcIBHCA490eIowcPAAQoQFDSpUsFCB
QhRbmzox+hRqFCRJEJQsMcK00ZHdR4IO6S0UIBAAOw==
""")

copyImg = TK.PhotoImage(data="""
R0lGODlhFAAUAOfwAC9HbDBIbTJKbjNLcTRNcDRNczROdjVPdDZPdjVQdzZQdjdQdjdRdzhRdzhS
eDlSeTpTejlUfDpUejpUezpVfTtWfzxXgD1Xfj1Zgj5Zgz9agT9agz9ahD9bhD9bhUBbhUBchkJc
hEFdh0JeiERehkNfiUNfikRhi0Vhi0ZhiUVijEVijUZijUZjjUdjjkdkj0lkjEhkj0hlkEhlkUhm
kUlmkUtmj0pnk0tok0xpj0tplE1pkkxplU1qlk5rl09rlU5smE9tmVBtmlBumlFumFFvm1JwnFNw
mlJwnVRxm1NxnlRynVRyn1VznlVzoFZ0n1Z1old1oFd1o1h2oVd2pFh3pFh3pVl3pFl4oll5o1p5
p1p6pFp6pVt6p1t6qFt7plx7qVx8ql18ql59rF5+rGCArmGBr2GBsGKCsWODsmOEs2SEs2SFtGWF
tWWGtWWGtmaHtmaHt2eIt2eIuGeJuGeJuWiKuWiKumiLummLumqMvGqNvHqMn2uOvmyPv2yQv22R
wW6Twn+SpY+hspKktJSmtZantpeot5iot5ipuJmpuJuruaGvvKGzwaixvKayv6ezwKe4xay3w6+5
xa28yLK8x7S+yLHAy7fAyrbCzbbDzrrCzLnFz7zEzbfI1LzH0b3I0r/J08DJ07vL1r7N2MTN1cHQ
2sLQ2sjP18PR28TS28TS3MXT3MbT3cvS2cnV3snW38/U28rX38vX38zY4M3Y4c3Z4dDb49Xa4NHc
49Lc5NPd5Njc4dTe5dXe5tnd4trd4tbf5tfg59jh59ji6Nni6Nni6dri6Nvj6dvk6tzk6t3k6t3l
6t7l697m69/m69/m7N/n7ODn7OHn7OHo7eLp7ePp7uPq7uTq7uTr7+Xr7+bs8Ofs8Ojt8eju8enu
8eru8urv8uvv8+zw8+zw9Ozx8+zx9O3x9O/y9e/z9fDz9vH09vP1+PT2+PX4+fb4+v//////////
/////////////////////////////////////////////////////yH+FUNyZWF0ZWQgd2l0aCBU
aGUgR0lNUAAh+QQBCgD/ACwAAAAAFAAUAAAI/gBzcMky5cmSIz9spAhxwYGAfxAjfgn0Z48dOGrI
XEHCo8UHBhEjbgHkR0+dN2nGVDGiY4WHBSEhYumzB8+cNmjESCmC4wQHBTH/Rclzh84bNWW6OAlS
w0SGA0GbyInjZs2ZMFSQ9HghwkKBoEnYrEljRowVJkFuoOhAYUBQImbKkAFTxYmQGCUyRDAAIOgO
L1qsQFEyJMYjUaVQuYrFC9gvXrgcwXDCREkRHzhKhIKGzdu5de3euWOXrhMJID564JjBYgOoZtW2
hTunbt06dOQ2aaAhw4UKEyAofEoW7Ro3cOPMlRP3DdOEESNEeMBQIQGnYsykWdPGrVu3bdksexGQ
AOFBAwQFAmQKdkzZM2rVrFmrNq1SUIiacgUjhmxZs2fONKPMJPf9c8kst/QizDDFGFNMMcNIUiAl
q8Biyy279OKLL73oAkmBkZBySiuvzFKLLbbQIgsjBTbiySikmJLKKqywosopixQoyCCEFGIIIoko
okgih/AREAA7
""")

getImg = TK.PhotoImage(data="""
R0lGODlhDwARAIQXAMfHx7a2turq6uPj4/j4+Lu7u9zc3JCQkLS0tImJiYGBgb29vaOjo6qqqvHx
8XFxcaurq7CwsJiYmNTU1KGhocHBwcTExAAAAMzMzP///wAAAAAAAAAAAAAAAAAAAAAAACH+EUNy
ZWF0ZWQgd2l0aCBHSU1QACH5BAEAAB8ALAAAAAAPABEAAAVC4Cd+FzCeKFmm6ZVdLOq+Zkxm+BrP
ed3iQBiLF/SNiMEfcCk8Ll+XaPP2ykmGr8dserwkFttCKyFCbG2fxhl9iI5CADs=
""")

pinImg = TK.PhotoImage(data="""
R0lGODlhGgAaAOevAFBQUFRWVzR+nTyAnzOHqm96fjmHrDWJrlyCjDqKqTOLtTKNtTmMsjSOtjiN
tjeOsjaQt0OOsDuQvIGChTyStkOQtj+RuTuSu0aRujyUuz2Uvj+UvECVvjyYvT2Yv0OYxUaZw0Ob
vUqZv0mZw0CdwIyMjEibyU+bxk2cxkWgwkiiw1KfyEakxk6kwEulyVikzFWl0Eqpzk6px5iZmpqa
m0+v0VOz01S00li11F+31Vi511y922q92VzB4GLE4GTE3V/F42vD22TI5G7F53HG4mXK5m3O5mrQ
6nXN6GvS623X83LY8HnX7HDa9HDa9Xrb83Pd93be93je9svKy4Hf93zh+X3i+n7i+Yng83rk+3/j
+oDj+oPj+YHk+onj+IXk+oLl+4Pl+4Dm/ITl+4bl+obl+4jl+4nl+obm+4fm/Ifm/Ynm+4nm/Ivm
/Izm+47m+ovn+4fo/ZTl94zn+5Hm+ZXm+Y7o/I/o/JDo/JHo/JTo/JXo+5bo/Jfo+pPp/JPp/ZXp
/Jbp+4/r/pfp/Jjp/JXq/Zbq/KDo+Jfr/Znr/Zfs/5/r/KDr+5zs/Zzs/p3s/Z7s/qDs/aPs/aTs
/KDt/qbs/KXt/aft/aru/aru/q3u/avv/q/u/azv/q7v/rLv/bDw/bHw/rfx/rjx/rnx/bnx/rrx
/bvx/rzx/rvy/r/y/sf0/8j0/tD2/9j3////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////yH+EUNyZWF0ZWQgd2l0aCBH
SU1QACH5BAEKAP8ALAAAAAAaABoAAAj+AP8JHEiwoMGDCBMqPNio0J01YBYmlEQp0Z85aMicqSKx
IKNNkBDlgVPGTR8vHQcuUrXJkSE8bMaYeUMlpcBTpjw9KmQnTRctVZzY/BRJj5w6igTFEZMFis1/
gLgs8ZGDBxEkQ55OwvQliQ4WEBSMePqPk6QwQlJk0GAzkyhQfEj5kXLjgQSbhFK1WjVJk5ojKhig
SPnlkitXqA4FurKDgoUVKfeMYlUKEx0mRmQcqGAgZSVVnSwNKiPlR4sEAgY4SKkqVKE2W6IUsRGC
QAQRH1JeIpTGShMgNUg02PDCxFMsT5T0iNFhAQcYGJ4WQBAEhwsPF0CcIFsCANmDMwIOfC84ZcL4
gjTOq19/MCAAOw==
""")

pin2Img = TK.PhotoImage(data="""
R0lGODlhGgAaAOfDAGRlZ2VlZ1iBjXl7fXp7fTSLtHt8fjWTwESauT6byDucyzycyTycyjydxjyd
yz6dyz2fyj+fyT+gyECgyECgykifvUGhykOjykSk1ESnz0iozUipzEOr1negq0mq10es3Uat2E6r
2JucnUev25ycnVis1E6u4Vus10ey3U6431O54mC56lq93FW+41m94lm942G83Vu+3VnA4l3A3FvB
4VvB42m961vC41zC42u97FvD41rE517G5mvC7G3C7XbA63HB8HLB73TD713L7mXK62fK7WHM7WnO
6mXQ72bS8XDR73HR7m7S8mnU8WrU8mzY9GzZ9MXFxm7b9njZ82vd/HHc93jb9XXe9nXf+Hbf+X7d
9nTg+HXg+Xbg+czLzXjg+Xfh+Xfh+njh+oHf93vh+Xvi+nnj+n3j+n7j+n7j+4Hj+YHk+nbn/4Tk
+oTl+9LS04Lm/Ibl+4vk+Yfl+4Tm+4jl+ojm+43m+4rn/Ivn+4zn/I7n+47n/Izo/I3o/I3o/Y/o
/JHo/JLo/JLo/Y/p/JTo/JHp/JHp/ZPp/ZTp/Zbp/ZXq/Jbq/Zrp/Jnq/Jjr/Znr/Zrr/Zvr/Zrs
/Z7r/J/r/Jzs/Z3s/Z7s/Z7s/qHs/Z/t/aHt/qPt/KXt/abt/aTu/qjt/aXu/qfu/qnu/q7t+6vu
/a/t/Kjv/6zu/qvv/rHu/a7v/a/v/rXu+7Dw/rPw/rbw/rPx/rbx/7fx/rvx/bry/rzy/r3y/7/z
/sDz/8Hz/sHz/8Lz/sT0/8b1/tD3/9j4/uP5////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////yH+EUNyZWF0ZWQgd2l0aCBH
SU1QACH5BAEKAP8ALAAAAAAaABoAAAj+AP8JHEiwoMGDCBMqNBjLVKM9bcpsebKwoKxcnyAFwuPm
TJgqTYxUFGgLFaZHiPzMWVOGCxQkI/+pQsXJEqNDedykMYOFScxUvmyNsqTI0B87Vor4iMkKGK1K
he6MWfJixIcVMf+dMkWIjhgnO1B4MJH1H5hSruRc4ZFhwY+yAl/xGtRFRgQFQuD+S7VLUJYaExgE
gUuJVN8vNyg4AKIhKyxNnnoJEoPDwoMcCC5wqDhLVCRNtybrkJCgRIUNICq2ApUpkiQ+aI6wiDGD
RgsVFVfh2jRpESA7ash0kZKESI+KjoIJq6UpUZ84a9BMgRHCRsVOv3SFoiSojhYlLiAvHMBwIuau
S3rgsKEyJEWDAnoFdhAQv6CIAPUJRhmQf+AbAv0J5IUBAf5DAgAVBQQAOw==
""")

refreshIconImg = TK.PhotoImage(data="""
R0lGODlhFAAUAPdoAAAAAAEBAQICAgMDAwQEBAcHBwkJCQwMDA0NDQ4ODg8PDxAQEBQUFBsbGx0d
HR4eHh8fHyAgICIiIiMjIyQkJCgoKCkpKSoqKjIyMjY2Njg4ODk5OTo6Ojs7Ozw8PD09PT4+PkFB
QUJCQklJSUtLS0xMTFBQUFJSUldXV1hYWFlZWWBgYGJiYmVlZWxsbHBwcHV1dXd3d3l5eXp6ent7
e39/f4KCgoiIiJGRkZWVlZaWlpubm6SkpKampqenp6ioqKqqqqysrK2tra6urq+vr7CwsLS0tL6+
vsLCwsTExMbGxs3Nzc7Ozs/Pz9DQ0NLS0tXV1dfX19nZ2dvb2+Tk5OXl5ejo6Orq6uzs7O3t7e/v
7/Dw8PHx8fLy8vT09PX19ff39/j4+Pn5+fr6+vv7+/z8/P39/f7+/v///wAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACH5BAHoA2gALAAAAAAUABQA
hwAAAAEBAQICAgMDAwQEBAcHBwkJCQwMDA0NDQ4ODg8PDxAQEBQUFBsbGx0dHR4eHh8fHyAgICIi
IiMjIyQkJCgoKCkpKSoqKjIyMjY2Njg4ODk5OTo6Ojs7Ozw8PD09PT4+PkFBQUJCQklJSUtLS0xM
TFBQUFJSUldXV1hYWFlZWWBgYGJiYmVlZWxsbHBwcHV1dXd3d3l5eXp6ent7e39/f4KCgoiIiJGR
kZWVlZaWlpubm6SkpKampqenp6ioqKqqqqysrK2tra6urq+vr7CwsLS0tL6+vsLCwsTExMbGxs3N
zc7Ozs/Pz9DQ0NLS0tXV1dfX19nZ2dvb2+Tk5OXl5ejo6Orq6uzs7O3t7e/v7/Dw8PHx8fLy8vT0
9PX19ff39/j4+Pn5+fr6+vv7+/z8/P39/f7+/v///wAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjmANGg2XICgAAACA0CgNBEoEOHNA70GBJEiJAg
RXgwcADlocMUITw6vABAghORKkCIFBgBYYQnHlOuRLMBAAEADxo6lLlSCowXNRCseMhTYI4VLpQg
aSHwQwqiKtEwmQCgwAADBSgIzKACqsANAXag0YGwgkANT3dGPZJEYIkGDUZ8TSuw6Ew0HOiisTuz
g14VHwSCQKCAhEAjQAT6JRoYjQWDONDcABDVw9/GFBIsUIhhStO/UVksiWIjxo+HIEDfdQii68eQ
q9GIQPFwxgEfRCpa3B2EyI8DMh5qMaEwofGDJrQIDAgAOw==
""")

invertSelectionImg = TK.PhotoImage(data="""
R0lGODlhGAAYAPIHAAEBAQICAgMDAwQEBAUFBQYGBgcHBwAAACH5BAEAAAcAIf8LSW1hZ2VNYWdp
Y2sNZ2FtbWE9MC40NTQ1NQAsAAAAABgAGAAAA194utz+MDpB6ag4t8w7ZV44HAIYcuORKt0krOXC
QStpVnSs1Fj+yjcd44LbBR3Eok2mZM1Yx0WyZ2R+dh3YbTCaUqHWivcKDKemNbDzJxma0m0huV1m
Xoh4y2tO7/shCQA7
""")

PHOTO += [saveImg, undoImg, deleteImg, copyImg, fitImg, selectallImg, eyeImg,
          mainImg, getImg, pinImg, pin2Img, refreshIconImg, invertSelectionImg]
