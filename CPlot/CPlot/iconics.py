# iconics Photos
try: import tkinter as TK
except: import Tkinter as TK

# Class to store and provide load on demand of PhotoImages
class Photo:
    def __init__(self):
        self._plist = {}
    def __getitem__(self, key):
        from .Ttk import ICONMODE

        if key in self._plist: return self._plist[key]

        if ICONMODE == 0: #lightMode
            if key == 0: img = createSaveImgLightMode()
            elif key == 1: img = createUndoImgLightMode()
            elif key == 2: img = createDeleteImgLightMode()
            elif key == 3: img = createExportImgIconImgLightMode()
            elif key == 4: img = createFitImgLightMode()
            elif key == 5: img = createSelectallImgLightMode()
            elif key == 6: img = createEyeImgLightMode()
            elif key == 7: img = createMainImgLightMode()
            elif key == 8: img = createGetImgLightMode()
            elif key == 9: img = createPinImgLightMode()
            elif key == 10: img = createPin2ImgLightMode()
            elif key == 11: img = createRefreshIconImgLightMode()
            elif key == 12: img = createInvertSelectionImgLightMode()
            elif key == 13: img = createLogoCassiopeeImg()
        else: #darkMode
            if key == 0: img = createSaveImgDarkMode()
            elif key == 1: img = createUndoImgDarkMode()
            elif key == 2: img = createDeleteImgDarkMode()
            elif key == 3: img = createExportImgIconImgDarkMode()
            elif key == 4: img = createFitImgDarkMode()
            elif key == 5: img = createSelectallImgDarkMode()
            elif key == 6: img = createEyeImgDarkMode()
            elif key == 7: img = createMainImgDarkMode()
            elif key == 8: img = createGetImgDarkMode()
            elif key == 9: img = createPinImgDarkMode()
            elif key == 10: img = createPin2ImgDarkMode()
            elif key == 11: img = createRefreshIconImgDarkMode()
            elif key == 12: img = createInvertSelectionImgDarkMode()
            elif key == 13: img = createLogoCassiopeeImg()
        self._plist[key] = img
        return img

PHOTO = Photo()

def createSaveImgDarkMode():
    saveImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7MK6iAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAi0
lEQVR4nO3SQQ6CQAxA0d5DvAnsIVzKJQcyHgi8xyesHInJlGkJKv3rdl7SjEh0uoAeeJKUmX+bSdZu
W+EpRQ3wNpwPGeCl4ShYh7MPjDus2S+GgVqJNt6wOQl4nf2of3DqO3DJLr72K+DhAVdaNHnjaoalsI
AlTv1Nn2vCv1EDd874CLSll4t+txluQ70WDrFYpQAAAABJRU5ErkJggg==
""")
    return saveImg

def createSaveImgLightMode():
    saveImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7MK6iAAAA
CXBIWXMAAAsTAAALEwEAmpwYAAAAh0lEQVR4nO3UwQmAMAyF
4beHdRO9Ky7l0YHEgdQ9Ij0I0oNNQopI80MOhZYPcijg1dgE
4ARAj3krvXOfZyl8JKgWFuNkCMdZvoLZOBWAqQTMeU/ahx0T
7a1hMphsDlN1q14BNOAXAGwWcIC81gLW5nA2X7U28eaOAp/H
zoFHYzyig3pv3m+7AKg9J4Kk0mDnAAAAAElFTkSuQmCC
""")
    return saveImg

####

def createUndoImgDarkMode():
    undoImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7MK6iAAAACXBIWXMAAAsTAAALEwEAmpwYAAAA0El
EQVR4nO3UIQ7CMBSA4XoIdyBwCVDbHNyHYLkFcAoIdnAMmAMxHCDI8D+pmHlpQgtdi9iv23x5TfOUav
ungEUUFCAKqouC6qKgFr2AAlgDWSjU1AEYxIB1d2Dsis/x0wMY/jT5h7MdIAE2BnzvBEvc4c7MgKdf4
453tgJeOcM17ng+FXChQgR0BVyFgnsCfoaC5VOfQsE7AS9DoKblkzT5mTLDpLrcFzIFLth1A/q+4KsD
OvKC6oDSAs29TVoHTICzgCrgqFdjYx+prU1Z9gbU8vhYAZrzyAAAAABJRU5ErkJggg==
""")
    return undoImg

def createUndoImgLightMode():
    undoImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7MK
6iAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAyUlEQVR4
nO3UMQ4BQRSA4b8n7iBcgmp3O+4jWrfAKYh2OQbbUd
ChEPoRiW1eJszY2RnF/Mmrv7zJ5EHsz5qEQlUoVIVC
VShUfZkHUABzIPOFKs1sgE4IWAEXoG+Ljx3hV6Bbdf
NPNYAEWGjwtS0scdNGGjytgtu0FPDsF7jEbUoFXOCp
poDvvuCWgG++YPnUO1/wSsBTH6ju+CR1fqZMs+lrcl
fIEDgYnswz0HYFnyzQHg47GqC5y03LBsBecyC279NY
20eKxTDpCdOoqIxxt+JEAAAAAElFTkSuQmCC
""")
    return undoImg

####

def createDeleteImgDarkMode():
    deleteImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7MK6iAAAACXB
IWXMAAAsTAAALEwEAmpwYAAAAtUlEQVR4nO2WSwqEMBBEvYSi9z
9OXI0Lj/NkUJDJwv6kBmHGgqy0eN1tUqbrHv2lgAGYjzU0+IvbD
/TAi1MrMAWgOT97pbUWT+VHp+93a5Us2IRfQN3gvhoV1tgcntEE
R+EyaAQuhzq/22I8Cx/DDxldaTttgOugAbgeGgC7E+4bo9bAyW2
uNjh3HCcc4ZCJ10tFEkkGJxGDEjg3/hbnKNQBL9mrizuRmvzslZ
fQZU3of/Qb2gDIApP3EFwGCQAAAABJRU5ErkJggg==
""")
    return deleteImg

def createDeleteImgLightMode():
    deleteImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7MK6iAAAACXBIWXMAAAsT
AAALEwEAmpwYAAAAtElEQVR4nO2WwQ6DMAxD/ROJtv//HHqCwz4HLp00EAMn
Naq0YSknar02UBPg1r/KAZRa3uAfIn4DMAKYa70APAPQtL98mN41kTv3unbr
H7JgBv4NSoNt06qZaNuZ58GAo3AZNAKXQ5n3Np08y1zDlY5OJT9pFi6FsvBL
oCw4knDyVsvglvi4muHe4zoZEQ6ZeD1UJJFkcEvEoAReev0WS69BwHZGl0gi
Nfm97jI0rAn9t35AC0Zz72ONUPbuAAAAAElFTkSuQmCC
""")
    return deleteImg

####

def createFitImgDarkMode():
    fitImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7MK6iAAAACXBIWXMAAAsTAAALEwEAmpwYAAACUk
lEQVR4nNWWPU9UQRSGtxH/gWghFGIjoUewwE6iITYifwFM+PglgJBAQlZbZclqh4qVJIaa2GmFH60f
jYRl85hzeUdPxrv3zi4E9E1uMnPOnHnnfMyZW6n8bwB2gLdnQZzhNIguAnPAc+ADf2DjZ8As0H2ShJ
eAx0CDctiaqh3yuKR3gR/a9CfwBBgHrjiyPsnWgX3JvgNjnZI+AJra6ClwuSzHQA9Qk8pspzrxtAkc
AtNOPgAsA++cxzZeMp1bNyPbZrLnHOU0hDcjBc6LMEQgD6Z7CHTJxorN8C2p6FRIWXgd6WtXPCvACP
AC2ARuAquu+LYc+YZk1ZQr01CRZDmVp4ZPwPUC2yHgs9YuStarvQ6AC0XEc5G3AwrhQRGpsx/WwS2/
/ZKFYpspMrTmYLgXebtSRur2WIq8vq95vcgodKQ+zUP1jkTrbiv0H4HRSDcom13Nr2r+Po9wm2JsRu
uNMGAv0lnRFWH7rIjfHCfUoyLfA251HOqc4hqPCmW1kghXkAuaT2heWFyh26xH18muyFAC6Q13na5F
TWS6yLDbNZAeyawNouYwXEL6RWvn22ogBuCRjGuadwGvJGso/INqo1saL7uW+RI4J9u6ZGuVMqht2n
v6u9uIfFEhbAXTzTvS0AW/lnrryMfcszjr5P06wK4jtPFCyKkjPdR3J4nUGU+5Z3Aj5NzpM0SyXhde
I51sizTy3N5TVCQ1XY9wN9F4Qofbd+Ftz9MWlV5N/Nmz6l1LzmkKbDP90tRzfm/r0p0cYSvk5fhUwN
HD8nfD/9fxC4rui3E7qA1hAAAAAElFTkSuQmCC
""")
    return fitImg

def createFitImgLightMode():
    fitImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7MK6iAAAACXBIWXMAAAsTAAALEwE
AmpwYAAACSklEQVR4nNWWPU9UQRSGnwb8B7AUSoEVSI8fBXQQzcZG+A1igvBLQJdEEr
LQwkoWOhCooDC2WkolYqvSuFnvZs1J3klObu7eOywQ5E0muXPunHnnfM7ALcQn4ONNE
Lc1rh0lYAHYAU4csX1vA/NA/1USDgDrwF9H1mnYmqoOeSk8B8616R9gA5gGhhzZfclq
QEOy30C5W9JXQEsbbQJ3I2J8D3gvuenOdmNpC0iAOScfBSrAF0ds32+BB27da+m2LmL
5gHNvIL0DLDsPZI2WDtArnXnJf8Um3bpzbyA9dMnzDhgH9oBdYAJYccl34Mi3JKsWkZ
a0QcPFdFnK34GHObqPgDOtfeNibns1gb484oWUtaNyYbOANOCxDm7xHZGspj0t7h2xo
0UvNK9obu6NRSVl9bTm9TylE1ebuOy1mHo8letPganUvzGX7WivNvA1i/C4oCNZEnmc
un/fUv/2CvY6vi7i3YK9ji7j6imRG+lkB1d/jnF1OrksIVBDaKtOYxHKb0nzmZjkCt3
GSgC1wUQlYnVahCeunIYlC717Lk+x3zUQK35v9ZnqNI/0h9YuSjYY20AMa1K2k6L2t+
9aZkVxPFR7HJN7Q8v8APRIty7ZKhEo6T713aZXDSHJydZElvakuuDPGGsDyu5atLgHj
OgA6WtxycU0kCYaz2JJA2bdNbjlYp73EBh07jXSl3SJsu7TtpKkpvIItRlqfkY50XDu
vbClWZlejXzsNZVI0TGNQZ+SrZ7xvK3r35US3viDnhTsYsls+P81/gFbbA+hcXORfQA
AAABJRU5ErkJggg==
""")
    return fitImg

####

def createMainImgDarkMode():
    mainImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7MK6iAAAACXB
IWXMAAAsTAAALEwEAmpwYAAAAb0lEQVR4nO3VQQrAIAxEUY/Xcz
f0IvYgUywUmoWKmAml5q2jHwUxpfBFAA7UCTPcFGEz6GCGpdHda
eE3+inXCAPYAJyYl8teI+GywEoeCd88ZpQIu111Wv2jF8aMUnuI
jBklwg+Pq5beB281E/7vAhcGfwmdhfhLAAAAAElFTkSuQmCC"
""")
    return mainImg

def createMainImgLightMode():
    mainImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7MK6iAAAACXBIWXMAAAsTAAALEwE
AmpwYAAAAbUlEQVR4nO2VQQ6AMAgE93m+W+JH6kPw0oMctByAWLuTcNt00g0JAPkoBw
B9GMkU62Aonr9qeZHuKEKzf7mGeANwOpZpNK2/5aYFSO/y0BqjMgaKy6pe+9BLUsbgu
TpRGQPFWlW1OA58VIb8nAuWjtuCVaO6OAAAAABJRU5ErkJggg==
""")
    return mainImg

####

def createEyeImgDarkMode():
    eyeImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7MK6iAAAACXBIWXMAAAsTAAALEwEAmpwYAAABl0l
EQVR4nO2Wv0oDQRCHY6FBNFErfQCxNDFYi71iJxa24kskWpigEUSUPISSYCfYWIldIsSHiPljKRKrfL
IwgbjZ2btogk0+2GbvN7+5nZ3du0hkzJj/BlgHssAjUAfaMuoydwKkhplwB3ghPGVg+y8JV4AHfs89s
Dxo0j3g02HWBNJAEpiRsQZk5JnNB7AbJuEEcAp0HCZFIOaJjQElR5zxyhlvX9KCUraiGtjv4UpuuHZ6
ABe4afpW6vCJAy3F69wWH6KTtrRTxgB4A2pA3sxZmiOP30FXlAC+PMKEZWoS2eQtjWk4jTawakQV/Pw
os6zSpuZoNB8VI6qOILHZZx9VI9pUjk+XZIhSnw1Q6g6w0RVeeYQZR3PlZeVacx17/C5ts+chHac54F
3xegIm7YAF4FUJKA1wgdwpHlVgXgtckq+LljwesFItaRlYDHrraeBWMWjJ5ZACZmWkZE+18t4Yz6Bq9
ZZsXzk6YamJR+AW9SEXQc6zGq0q2UEaUgWIAlvy9TLd3+j59WnIXEE0Ud1pzJjI6PkGDLCOslXUiekA
AAAASUVORK5CYII="
""")
    return eyeImg

def createEyeImgLightMode():
    eyeImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7MK6iAAAACXBIWXMAAAsTAAALEwEAmpwYAA
ABkklEQVR4nO2WvUoDQRSFPwsjoomk8wHEUpNgL9aKnVjYii/hT2FEIwgS8hBKgp1gY6ddYpGH
WPNjKRKrRAZuYBj3zu5iRIs9cJs7557D3Duzs5AixT/AGlAGHoEOMJDoSO4UKE3ScBt4AUYxow
ls/cRwGXhIYOjGPbCU1HQX+AgR6wGHQAGYkygCR7Lm8t+BnTiGU8A5MAwRqQNZT61Za4TUGa0z
0VZNa0rb6r5CRyPMfARUNY0rpaAXsVMXOaCvaF265APPITEztZERgVcgACqSs3Hs0dsfk1aBTw
/RrNuohHBMzkbRozcAVgypFXEt3DYHIRyTs5GN0GwZUvsXjHMRmm1D2lCuzzgKMVp9kaDVQ2B9
TKx6iObjYCMj5oHncJ149K5dsecJXacF4E3RegKm3YK8Z96NBB+QO89c81rhorwumnkuYqeaaV
O0vZgFbhWBvnwczLs7L1GSmWrtvRHNWDAt21OuTtwIRCPOiL4hKy+LthutK+WEB1LFDLApr5c5
/V3r16cruZpwDDdFCv4MXw/YJO5+W1zLAAAAAElFTkSuQmCC
""")
    return eyeImg

####

def createSelectallImgDarkMode():
    selectallImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7MK6iAAAACXBIWXMAAAsT
AAALEwEAmpwYAAAA2ElEQVR4nO3WLw+CQBiA8StWDUZNJjcxaTQYZSbn5ufx
c+o0gJsDq/1xuDcgwjiO14DwNBjHb/dnDGO6WhewBW64FwK+CxxSv8AFVunv
4QhYAKssbL3nuKGejF3mwHZ7TrViYC7jpsC9AEYTjopQDXgP7IBnDjqT5z25
RhOeyL11Ci9FNeATMJL7yYm9pA7S1/Jqwlm8VzZTTTjpDIxtZqoNJ12BI/DA
Ik24Uo2CQwX3/YkEBlVgvyYeABt518Ea1goYpidgfh3Ql5l+rJpp0q+Qr7Xn
Xe3qBSbm7FwdlDawAAAAAElFTkSuQmCC"
""")
    return selectallImg

def createSelectallImgLightMode():
    selectallImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7MK6iAAAACXBIWXMAAAsTAAALEwEAmpw
YAAAA4ElEQVR4nO3WsQrCMBCA4X9x1cHVQZwE66Sjg6PFSQSfx/cUFawg1dU9UrhCKVGT9B
xqe3BD06YflztKoY0mxhq4AiYwEyAOgZMKaJ6XENgo5X/DKTADFhbYuecmAI1k79wCO/fce
OQdmMq+MXB7AxtNOP2AVoa3wAZ4WtCJPB/JtdGER7K2LOAuaGX4DAxlPZvYY2GQbMerBhvg
AAzkXsehUtXhOhcq/1ap+lSfgD3wcJx6Ndh4Zn3gRAHNP5E9HziuiGfoSt6184G1ol8q4Of
RlUrLp1afX6FYsedtNChe6YxhwyFtFlAAAAAASUVORK5CYII=
""")
    return selectallImg

####

def createRefreshIconImgDarkMode():
    refreshIconImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7MK6iAAAACXBIWXMAAAsTA
AALEwEAmpwYAAABLElEQVR4nO3WP0oDQRTH8ek1N9BSLFJoLEwMOUDES6SwSA
QPEzR1+hC8TS6QP6VgxGIt9CsPX2AJM0tm562Q4A8Glpn3+LCzs8s69x/nHHA
HLNk9GTBI3jxgQVy+dTykwrG5B74Uf/wz2P32pOOUgE1wSsLJOAmw9vdzB25Q
GQx0CvCsSjhqR/YKzoAh0ASOdMj1E/BZFbwELgp6LkOf2RQ426DAKTAF3nW8A
Oe61vDdeQo8zKGvnnWZO9GaZ0v4WuflTkOZaE3LEj7WednaUNZaU7OAe561Xg
B+s4Q/gHpuvq5zvky15sYClsxy76xcEzhcZ9o7soIlYx3bWeuB26BX1q+TL+2
tHkFXvkJruCOnXZ/pqMpPZunsBRz7e1uUeQx8a4TPge7O8EHnBwTF9IXD9sZB
AAAAAElFTkSuQmCC"
""")
    return refreshIconImg

def createRefreshIconImgLightMode():
    refreshIconImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7MK6iAAAACXBIWXMAAAsTAAALEwEAmpwYAAA
BK0lEQVR4nO3VvUoDQRDA8X9vfAMtJUUKPwo1wQdQfIkUFkbwYURTp5cjb5MXSGIZiMHiUsSVwQ
kcYW+5zc2KBgcGjtldfuzO3R78x3fcAlPAVcwc6GEQkwhU8lPzoS7sIvMOWCn++JMwVrjbAjbB3
ZZwbdzVgCXuCy9cLyV8FcDzlHDsifwdOAeegAtgT1Oen4FlKngKHAfWnASu2crh2+kaPQQyYKE5
BJo6dlqy88qxuVCOd43OPONSO9A5L5bwudazQB9fdc6lJdzQ+iIAv+ucfQu46xnrlsBzS/gDaBX
qLa35YGmDRNvqqEeFb3ZUgsrLdaRr+5Y9Hmj6+poV0DPrz8l5srOxRtC3FDeX8/yRGtrTfsor09
XI3w9PDNFxDHxjhAt6HQPvbnwBMTWbmOKXHaoAAAAASUVORK5CYII=
""")
    return refreshIconImg

####

def createInvertSelectionImgDarkMode():
    invertSelectionImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7MK6iAAAACXBIWXMAAAsT
AAALEwEAmpwYAAAAiElEQVR4nO3VTQ5AMBBA4TkecSiWztm6x7OxEKGdqfET+
tbT+aIRRGq/C+iACX0RaD3giL3gARdVYXMUZpc2VdihADSSy0Paw+UqOLdP3g
xHz6u0wK0RT748ajiXddF7YZQLvw0D41MwW/wKOHDcoJyzfzKBJrF0UM6V/yT
WLQd7uTtWT1qTk82PRqkd8x0TsgAAAABJRU5ErkJggg==
""")
    return invertSelectionImg

def createInvertSelectionImgLightMode():
    invertSelectionImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7MK6iAAAACXBIWXMAAAsTA
AALEwEAmpwYAAAAiElEQVR4nO2UUQrAIAxDczxlh5qfO6fuHhn72N9wUQoVZ6
B/tQ/aRGDpj9oAnAAoVgEQLcClAfpUtgCzsxa4WVyrxuDmYiVqwQMs55zG5pJ
PQS9wMV6lDI6N8C/zmLte1bhgigPnBh9eYL7AzcG5Epkk9nV9maEyNIl9PTl/
1f1wh4OSBxRT6gJMAiJ68qE2fgAAAABJRU5ErkJggg==
""")
    return invertSelectionImg

####

def createExportImgIconImgDarkMode():
    exportImgIconImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7MK6iAAAACXBIWXMAAAsTAAALEwEAmpwYAAAA1
UlEQVR4nO3WQQqCQBTG8ekUEdlNWrQsWnSJNgbtW3qEILpJ2GVaWVcI3NQ/pBfMIsWZNwmaHwzI+N
785KGgMX2MMcASuFE/ORCrhwdccctT1kYLu2YNPATfNgabd48exwMOguMJq3EUsPTH1gsX/wwGphV
47gufgDEQAanHQ+ELR9b+5C/gVMZcoOcmYXV6uDWjPgIjWcX1Jxf5zofAITScfKlJrPv7kn0VnFTU
7aTmDgzKcC+4Ru0KmGnO0DV1Bg6RVsCuv7dVyVzgRSA8A+a14U7nBbDgh0RUNOYPAAAAAElFTkSuQ
mCC"
""")
    return exportImgIconImg

def createExportImgIconImgLightMode():
    exportImgIconImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7MK6iAAAACXBIWXMAAAsTAAALEwEAmp
wYAAAA0UlEQVR4nO3UTQqCUBSG4bdVRFQ7adCwCFpFE4PmDV1CEO4kbDONsi0ETurEpSM4
KPH+ZKB+8E2u5/qIXIU+76yAGyA1mwMRAZJZoKZP7dYXFstugIfiuyZhQuHiAAfBxRH2xs
UDRk94ceCiX8KzCjx3hU/AGJgAqcNDiStswCLTTsCp4gY9NwlLgNZOD8s/X3UCjLRJaf2i
3/kQOIaG4w8zcen64cu6FxxXzO115g4MKvDasdm0Buae9/Db1BpYOvHnygKiVxt4GQg36M
IGbm9eXP2K+Xnrja4AAAAASUVORK5CYII=
""")
    return exportImgIconImg

####

def createPinImgDarkMode():
    getImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAACNiR0
NAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAyElEQVR4nN3
RIQ7CQBCF4YIqAi5A74AhAYEn4QQcgeBag0GC4Ay9BwG
uUlUJCcWi+MmGbVKa2bJsKwhPzma+zMx63t8FmAMpMGg
C6wNXXsmAYR2sBex5T+aMAssS5o4CfmFVKWlVc09PEwN
bYKzrM+AuYDf1UVXgqdTwANYG9AgEn9YLgESYZKffp/p
uK6BtezMTmk/qW0ECehbWH3kuASINXEropg4WCevH32J
hjhVqRXThgoXCm0IPQNcWm5gwpwAdhTaC/XSeQduqIrT
hI8EAAAAASUVORK5CYII=
""")
    return getImg

def createPinImgLightMode():
    getImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAACNiR0N
AAAACXBIWXMAAAsTAAALEwEAmpwYAAAAwklEQVR4nN3R
PwrCMBTH8a9OOugF7B1cBB3cBU/gEcStXVw66uAZvIeo
V+nUUUFdnawEnhBCUtukU3+QJX8+vPcCbcwKyIFxE9gI
eAAF8AQmIVgHOAlWNIFuDCwI7Wmt2lZe9ngo1RyBPTCT
/SXwtmAv+ShnrsaDD5A60AsQ/WtPXcgslRzkfCFz2wLd
qjNzoak209qJgJul/SmeSQS4G+guBEss7avfr5VYw37R
0bUPFlvOFHoGBlWxeQnmlb6gLc8XZJVPb3UgrVUAAAAA
SUVORK5CYII=
""")
    return getImg

####

def createPin2ImgDarkMode():
    pin2Img = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAACNiR0NA
AAACXBIWXMAAAsTAAALEwEAmpwYAAAAzklEQVR4nK3UPW
oCQRgG4G1CQvAMySkCUXOWFNYW0UrvIHgMLRUhkMJG8ED
B3if4ExiXdXZ2mbfb2Z2Hb2e+mQJdDItcwRBHjHKioxQU
M3zjOQuKlUu2qejXdcK04t3AbXboNKl0HIy94lAC21eKd
QXWutKl+vzEoCe8B2hdVpjHwF642xH0F5+1v1paw0np+T
8bvCRh9/oyqHRRtAke7qA3LZWU68R9FlTFEQzG+o1QkfO
Mj1NLpXx7Dt6aXmNRFI8nNBVLQtsmQLPf/N0/T4QjDE1c
ciIAAAAASUVORK5CYII=
""")
    return pin2Img

def createPin2ImgLightMode():
    pin2Img = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAACNiR0N
AAAACXBIWXMAAAsTAAALEwEAmpwYAAAAyElEQVR4nK2U
MQoCMRBFXyOKeAY9heCqZ7Gw3sKx0jsIHkPLFUGw2Ebw
QGLvipLAEpKYxHyYIjOZx5DMDEABlGRUCbwAyQmVQOgO
uAD9XNAKaIA6FLpSCVtLbKli2m7AIKbSdcs3Ap4G8K9K
TxZYcqVHD6xRdvWBesCkBf0Fq4C9Dzg1ftsFfQALAqXf
cGOctZ2BIZESR6UHEtVxQM2WCtIn8Z4LKpYR1L5ZLFQ8
8zxXLRVy96txwhoTX05XQWMVuvqSoNk3f/EGNJZisCcQ
dnEAAAAASUVORK5CYII=
""")
    return pin2Img

####

def createGetImgDarkMode():
    getImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAACNiR0NAAAA
CXBIWXMAAAsTAAALEwEAmpwYAAAAnklEQVR4nNXUOw5BQRSA
4Sk0olRbgV1YgYZGLENiCSxBYwFKiQUotHagsQuNT25McT2u
x9xJxJ+cZE7z5VQTwl+GOdZo5gJXrm3RygnmQd2C9VGPYD3U
czAdVQ2moV6D36Peg0WzuuASg9J0UsATJvG9/xioAM8Yx30X
0V5ICX2MSvswgpsk8D40cIxXd0OOMI1XLnKBbRyKry0L+LMu
qm5wD0wN0zsAAAAASUVORK5CYII=
""")
    return getImg

def createGetImgLightMode():
    getImg = TK.PhotoImage(data="""
iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAACNiR0NAAAA
CXBIWXMAAAsTAAALEwEAmpwYAAAAnUlEQVR4nNXUsQmDUBDG
8X+RJlimdoJs4QRptJGMIWSEZASbDHCl4AAWabOBTbaw0SBc
8Qgviu9d4wcH9zU/rjrYax5AAxytQAEmoAMSS3CyQsUBTVD5
AaNR8YBRqPwBg1FZAINQWQHnuceCTyB3Jg0BB6DS/b0F8IEj
cNX+UjQjMBegdHqhYItRDsBHrz5boTe9srYCT0Cvr23H+QKl
lEgpqPve+gAAAABJRU5ErkJggg==
""")
    return getImg

## OLD ##

def createLogoCassiopeeImg():
    logoCassiopeeImg = TK.PhotoImage(data="""
R0lGODlhEAAQAMZUAAAAAAABAAADABUPFR8SKCQUMCkTPioaOjcbUUsTf0MaWEkbZkocZXUD/lMb
eUAkWnMF91IcgHML5YMG/4AH/34J/28S2n0O+mYdsnYT6YkM/0s4AoQR/1ssdYcU/34Z6YsU/30Z
/2I0RWooqUxBC1g8NnIqoXYl15gh/5ci/YYuyIwp6n02upY122hZEpFAoZ00/KE5+INiCIdoB4Nq
AIVtFJNtCYt1CZRzGZF4Bp9ompR6DKZ8Bb93q7aTD6+aG8OdN8agLcmgM8ujONOyH+CxLem0AO23
BNa2YuK6G+65C9S4bt+2k+PBZ//GBNjGkeHHff/KCP/TCvngaP//////////////////////////
////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////yH+
EUNyZWF0ZWQgd2l0aCBHSU1QACH5BAEKAH8ALAAAAAAQABAAAAeJgH+Cg4SFhoYhJ4RNhyCCCoRQ
G4YjAH8UJoNPMoMVEX8pKn8sloNDJRoLf5YeAoIShEs8NS49OoUJpYJFR1FSTjQNDoIohS84Nz5K
Nh0BghkFhC2WEw9CDB+HhBgATCIHhxwWBH8rgjOEP4QDgheuQIJTRIQwgwiCSEZJOYQxEAaFguwg
oa3goEAAOw==
""")
    return logoCassiopeeImg