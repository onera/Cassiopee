# iconics Photos
import tkinter as TK

# ICONFORMAT = 'PNG' or 'GIF'
ICONFORMAT = 'PNG'
try: createSaveImgLightMode()
except: ICONFORMAT = 'GIF'

# Class to store and provide load on demand of PhotoImages
class Photo:
    def __init__(self):
        self._plist = {}
    def __del__(self):
        for p in self._plist: self._plist[p] = None
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
            elif key == 14: img = createCopyImgLightMode()
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
            elif key == 14: img = createCopyImgDarkMode()
        self._plist[key] = img
        return img

PHOTO = Photo()

if ICONFORMAT == 'GIF':

    def createSaveImgDarkMode():
        saveImg = TK.PhotoImage(data="""
        R0lGODlhHgAeAPAAAP///wAAACH5BAEAAAEALAAAAAAeAB4AAAJYjI+py+0Po5wJWLCsue9q9XmOh4El2aDmpprcGpCl+6VtJY43Ivf+jPsJgbyhkGE8ZpI+ZE73okmeQaKtVp1GsxFqEQu1ssBXcWwbdjJl5bWZAo/L53RFAQA7
        """)
        return saveImg

    def createSaveImgLightMode():
        saveImg = TK.PhotoImage(data="""
        R0lGODlhHgAeAPAAAAAAAAAAACH5BAEAAAEALAAAAAAeAB4AAAJYjI+py+0Po5wJWLCsue9q9XmOh4El2aDmpprcGpCl+6VtJY43Ivf+jPsJgbyhkGE8ZpI+ZE73okmeQaKtVp1GsxFqEQu1ssBXcWwbdjJl5bWZAo/L53RFAQA7
        """)
        return saveImg

    ####

    def createUndoImgDarkMode():
        undoImg = TK.PhotoImage(data="""
        R0lGODlhHgAeAPAAAP///wAAACH5BAEAAAEALAAAAAAeAB4AAAJDjI+py+0PIwNSAlrfxXnu7m3iyGnkKULoekUsqn6h7NCzWy6jG7egbvsdUsIEsWjYIZPBHwn0grKcq+Jpic1qt9xFAQA7
        """)
        return undoImg

    def createUndoImgLightMode():
        undoImg = TK.PhotoImage(data="""
        R0lGODlhHgAeAPAAAAAAAAAAACH5BAEAAAEALAAAAAAeAB4AAAJDjI+py+0PIwNSAlrfxXnu7m3iyGnkKULoekUsqn6h7NCzWy6jG7egbvsdUsIEsWjYIZPBHwn0grKcq+Jpic1qt9xFAQA7
        """)
        return undoImg

    ####

    def createDeleteImgDarkMode():
        deleteImg = TK.PhotoImage(data="""
        R0lGODlhHgAeAPAAAP///wAAACH5BAEAAAEALAAAAAAeAB4AAAJQjI+py+0Po5wKAGmfzXq37X3gtYxcZYojujomabyRHNDzm1I4SB073KvZgsIc8Xf8AW/KpSs3VPF8Rka0OC1VEddYy/plOcUdzJiITqvXkgIAOw==
        """)
        return deleteImg

    def createDeleteImgLightMode():
        deleteImg = TK.PhotoImage(data="""
        R0lGODlhHgAeAPAAAAAAAAAAACH5BAEAAAEALAAAAAAeAB4AAAJQjI+py+0Po5wKAGmfzXq37X3gtYxcZYojujomabyRHNDzm1I4SB073KvZgsIc8Xf8AW/KpSs3VPF8Rka0OC1VEddYy/plOcUdzJiITqvXkgIAOw==
        """)
        return deleteImg

    ####

    def createCopyImgDarkMode():
        deleteImg = TK.PhotoImage(data="""
        R0lGODlhHgAeAPAAAP///wAAACH5BAEAAAEALAAAAAAeAB4AAAJYjI+py+0Po5wL2IuzZrp7W33iRTlYaYLGKB5kwI6rGrcwXXszsHe9f+MFgbLgL3NEGocbZvPV5DxxpwaRev3sklBbV+dUqsJfcCJnRqCBinUVBY/L53RDAQA7
        """)
        return deleteImg

    def createCopyImgLightMode():
        deleteImg = TK.PhotoImage(data="""
        R0lGODlhHgAeAPAAAAAAAAAAACH5BAEAAAEALAAAAAAeAB4AAAJYjI+py+0Po5wL2IuzZrp7W33iRTlYaYLGKB5kwI6rGrcwXXszsHe9f+MFgbLgL3NEGocbZvPV5DxxpwaRev3sklBbV+dUqsJfcCJnRqCBinUVBY/L53RDAQA7
        """)
        return deleteImg

    ####

    def createFitImgDarkMode():
        fitImg = TK.PhotoImage(data="""
        R0lGODlhHgAeAPAAAP///wAAACH5BAEAAAEALAAAAAAeAB4AAAJwjI+pywzQYnuynoehXbnjfX0BNYqWWSLopBme2jqtl4WsVlN1GcMuqRABE6RhENf7/ZIw5FHpiw5pwB3vde2EqlRpz1hMhplbKNHJMp+XMvQRqYPsVtGsVs2JYfHtOg80B0d2QzVYsQcik3houLhYAAA7
        """)
        return fitImg

    def createFitImgLightMode():
        fitImg = TK.PhotoImage(data="""
        R0lGODlhHgAeAPAAAAAAAAAAACH5BAEAAAEALAAAAAAeAB4AAAJwjI+pywzQYnuynoehXbnjfX0BNYqWWSLopBme2jqtl4WsVlN1GcMuqRABE6RhENf7/ZIw5FHpiw5pwB3vde2EqlRpz1hMhplbKNHJMp+XMvQRqYPsVtGsVs2JYfHtOg80B0d2QzVYsQcik3houLhYAAA7
        """)
        return fitImg

    ####

    def createMainImgDarkMode():
        mainImg = TK.PhotoImage(data="""
        R0lGODlhHgAeAPAAAP///wAAACH5BAEAAAEALAAAAAAeAB4AAAJVjI+py+0PowO0UsmsxUtXrngXiIgAWYpoeK5d6yYjadZwZtfTfczp/go2fD2egSgzBpBLZezJ8RQ3Rw2L2sRKk1prFTv9fMVZcni2LSPTbC/0DY+jCgA7
        """)
        return mainImg

    def createMainImgLightMode():
        mainImg = TK.PhotoImage(data="""
        R0lGODlhHgAeAPAAAAAAAAAAACH5BAEAAAEALAAAAAAeAB4AAAJVjI+py+0PowO0UsmsxUtXrngXiIgAWYpoeK5d6yYjadZwZtfTfczp/go2fD2egSgzBpBLZezJ8RQ3Rw2L2sRKk1prFTv9fMVZcni2LSPTbC/0DY+jCgA7
        """)
        return mainImg

    ####

    def createEyeImgDarkMode():
        eyeImg = TK.PhotoImage(data="""
        R0lGODlhHgAeAPAAAP///wAAACH5BAEAAAEALAAAAAAeAB4AAAJbjI+py+0P4wGUwmoVzm4D/UkJN4ViCZ4IabDttrixV6kvian5WEs7KANyhj8G7GY6Gn/M4qwXIEI701W1oSxdscWsjxbkPbeu8KNs0qVfNjbPOZZ54+26/Q4pAAA7
        """)
        return eyeImg

    def createEyeImgLightMode():
        eyeImg = TK.PhotoImage(data="""
        R0lGODlhHgAeAPAAAAAAAAAAACH5BAEAAAEALAAAAAAeAB4AAAJcjI+py+0P4wGUwmoVzm4D/UkJN4ViCZ4IabDttrixV6kvian5WEs7KANyhj8G7GY6Gn/M4qwXIEI701W1oSxdsUWnjxbkPbeu8KNs0qVfNjbPK5Rll9u2/Y5PFAAAOw==
        """)
        return eyeImg

    ####

    def createSelectallImgDarkMode():
        selectallImg = TK.PhotoImage(data="""
        R0lGODlhHgAeAPAAAP///wAAACH5BAEAAAEALAAAAAAeAB4AAAJdjI+py+0Po5wJ2IsziLpn6IXbI3pgaRlpMC5ou64KelxqW4m1zeII/7kBfTuYTHNjYGrFZU/pbDpluWgo6XpNiULtiPrzxrg9MTar3VFI4HW4PROTqy+3/Y7P6w0FADs=
        """)
        return selectallImg

    def createSelectallImgLightMode():
        selectallImg = TK.PhotoImage(data="""
        R0lGODlhHgAeAPAAAAAAAAAAACH5BAEAAAEALAAAAAAeAB4AAAJdjI+py+0Po5wJ2IsziLpn6IXbI3pgaRlpMC5ou64KelxqW4m1zeII/7kBfTuYTHNjYGrFZU/pbDpluWgo6XpNiULtiPrzxrg9MTar3VFI4HW4PROTqy+3/Y7P6w0FADs=
        """)
        return selectallImg

    ####

    def createRefreshIconImgDarkMode():
        refreshIconImg = TK.PhotoImage(data="""
        R0lGODlhHgAeAPAAAP///wAAACH5BAEAAAEALAAAAAAeAB4AAAJgjI+py+0PopzxKYqlRZmGun2dJoYjWVonaLDNCiTuAl/xC2dOruMY12P8aB4fa3IoCpEippII6rRQ0Ju0eVvanNQtsMhcRofP7SgZrpaHYva1zZvJ4mkzHULPhvb8fqAAADs=
        """)
        return refreshIconImg

    def createRefreshIconImgLightMode():
        refreshIconImg = TK.PhotoImage(data="""
        R0lGODlhHgAeAPAAAAAAAAAAACH5BAEAAAEALAAAAAAeAB4AAAJgjI+py+0PopzxKYqlRZmGun2dJoYjWVonaLDNCiTuAl/xC2dOruMY12P8aB4fa3IoCpEippII6rRQ0Ju0eVvanNQtsMhcRofP7SgZrpaHYva1zZvJ4mkzHULPhvb8fqAAADs=
        """)
        return refreshIconImg

    ####

    def createInvertSelectionImgDarkMode():
        invertSelectionImg = TK.PhotoImage(data="""
        R0lGODlhHgAeAPAAAP///wAAACH5BAEAAAEALAAAAAAeAB4AAAJQjI+py+0Po5wJ2IstDblf6nlgmI3kN52YqbKnq5ZM/C40oGg2neN73GsAgynfzBjRHYnMXROhdFaMyOiTQ806rNiO4fZbfVsP5GZ2TqvXiAIAOw==
        """)
        return invertSelectionImg

    def createInvertSelectionImgLightMode():
        invertSelectionImg = TK.PhotoImage(data="""
        R0lGODlhHgAeAPAAAAAAAAAAACH5BAEAAAEALAAAAAAeAB4AAAJQjI+py+0Po5wJ2IstDblf6nlgmI3kN52YqbKnq5ZM/C40oGg2neN73GsAgynfzBjRHYnMXROhdFaMyOiTQ806rNiO4fZbfVsP5GZ2TqvXiAIAOw==
        """)
        return invertSelectionImg

    ####

    def createExportImgIconImgDarkMode():
        exportImgIconImg = TK.PhotoImage(data="""
        R0lGODlhHgAeAPAAAP///wAAACH5BAEAAAEALAAAAAAeAB4AAAJZjI+py+0PopzxKYqlRZmGun2dJoYjWVonaLDNCiTuAl/xS8Kz7OoOKlr9cj6czmM8Am3KJa+54yAnLSTDGrxhmTfegWo0Rb9j5rALKWebZ2ib/S6G5vS6oQAAOw==
        """)
        return exportImgIconImg

    def createExportImgIconImgLightMode():
        exportImgIconImg = TK.PhotoImage(data="""
        R0lGODlhHgAeAPAAAAAAAAAAACH5BAEAAAEALAAAAAAeAB4AAAJZjI+py+0PopzxKYqlRZmGun2dJoYjWVonaLDNCiTuAl/xS8Kz7OoOKlr9cj6czmM8Am3KJa+54yAnLSTDGrxhmTfegWo0Rb9j5rALKWebZ2ib/S6G5vS6oQAAOw==
        """)
        return exportImgIconImg

    ####

    def createPinImgDarkMode():
        getImg = TK.PhotoImage(data="""
        R0lGODlhFAAUAPAAAP///wAAACH5BAEAAAEALAAAAAAUABQAAAIojI+py83gHoARzVvPxXVvqXmLaHySSTlkFqyd21CmaqWxYo/szvdsAQA7
        """)
        return getImg

    def createPinImgLightMode():
        getImg = TK.PhotoImage(data="""
        R0lGODlhFAAUAPAAAAAAAAAAACH5BAEAAAEALAAAAAAUABQAAAIojI+py83gHoARzVvPxXVvqXmLaHySSTlkFqyd21CmaqWxYo/szvdsAQA7
        """)
        return getImg

    ####

    def createPin2ImgDarkMode():
        pin2Img = TK.PhotoImage(data="""
        R0lGODlhFAAUAPAAAP///wAAACH5BAEAAAEALAAAAAAUABQAAAIyjI+poN38oJMIWGdotjpRzmEBBIJiiZ5ouBzp6bLiGM+Va9P6Pt0Z1vF5bMFckZhzFAAAOw==
        """)
        return pin2Img

    def createPin2ImgLightMode():
        pin2Img = TK.PhotoImage(data="""
        R0lGODlhFAAUAPAAAAAAAAAAACH5BAEAAAEALAAAAAAUABQAAAIyjI+poN38oJMIWGdotjpRzmEBBIJiiZ5ouBzp6bLiGM+Va9P6Pt0Z1vF5bMFckZhzFAAAOw==
        """)
        return pin2Img

    ####

    def createGetImgDarkMode():
        getImg = TK.PhotoImage(data="""
        R0lGODlhFAAUAPAAAP///wAAACH5BAEAAAEALAAAAAAUABQAAAIojI+pyw3Q3IvrTYrsxUFv6n1SiIUWaIoJyZ0qc70lxGV0fciRjvdLAQA7
        """)
        return getImg

    def createGetImgLightMode():
        getImg = TK.PhotoImage(data="""
        R0lGODlhFAAUAPAAAAAAAAAAACH5BAEAAAEALAAAAAAUABQAAAIojI+pyw3Q3IvrTYrsxUFv6n1SiIUWaIoJyZ0qc70lxGV0fciRjvdLAQA7
        """)
        return getImg

else:
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

    def createCopyImgDarkMode():
        deleteImg = TK.PhotoImage(data="""
        iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7MK6iAAAACX
        BIWXMAAAsTAAALEwEAmpwYAAAAhElEQVR4nO3WTQqAIBBAYQ9S
        FxHbJ923A1n3eBG4CRRNGvvzrYUPxkFUqvXbgAlYubYFsCl4Py
        SRS8FiqRoBHTBXhz3eB+Eai0QEFl8kIrBYr4c1YKrDqmBCDc4d
        tc54HIzEqE0GPEjAn1oudwfsgPE07BdpKIGDIdhjvzlWAD/eZ+
        t3bTevaj5hdmBpAAAAAElFTkSuQmCC
        """)
        return deleteImg

    def createCopyImgLightMode():
        deleteImg = TK.PhotoImage(data="""
        iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7MK6iAAAACXB
        IWXMAAAsTAAALEwEAmpwYAAAAhklEQVR4nO3WTQqAIBBA4XeQuo
        jYPum+Hci6x0RQqzQl/Cnywaw/0BGE1p+bgBWQhLMAJgQviVE5x
        oZgyThF6oC5BrzX++ASiyQuuMQiiQsusUjySVgBugZ81mBJfdSK
        cDoHrCPgoS0X90dtazwnC4xPYB1xpz7YmWScd35zTAb8cp+tf7U
        BsEzPhL9359QAAAAASUVORK5CYII=
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
        iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7
        MK6iAAAACXBIWXMAAAsTAAALEwEAmpwYAAABt0lE
        QVR4nO2WTyuEURTG31Espmym+AIsFJZoqLHxBViI
        FJmF+QQWdobEmrD0MWRDiWTMwkb+RGEoNkpYiZ9u
        neH0dueO+5pQPHUWc85znme6973n3iD4x08AqAHm
        gGs+sBNBZ1f1G61Zo+1qMIQwtiMY71h0ZlwNV0Lq
        9DVzaHaJZsFFehFSrILGMdF8+ZXGW649BdqBJWAf
        uJfYl1xbmT3fjPKvG4BVymMdaPI2sAHoA56U+CEw
        DDRKjABHqv4I9AZfAZBR+26wBsQtvLjUijA9Y1FN
        B0Kmd0BC1atNqN8J4Wjzfl/T1tDyGsyr+pDUzbIO
        qvxCqMdwWnyMc5YPZ1TVb1X+RuVHLX17PsZ5i0Ba
        1U9V/kTl0181TgLPIYHF0BjMSyRV3pxnDaPR8Wlj
        EZkIiZhhUefg1wtHY9zLtAhgOSS0AdRaeLVSs66Q
        N2TOZoFXJXgm57tZIiO5Igx3siLzHkgBB5SH4aSC
        iC+HXIlaFdADrADHwIPEseRMrapEb975kvmT93FB
        SP77VFqzWzQvXaRpyweTi2Bmm3xZV4O5caaAC9Xg
        /XKQl0wR53K83m+yfwTfiTeflAtCzSughgAAAABJ
        RU5ErkJggg==
        """)
        return eyeImg

    def createEyeImgLightMode():
        eyeImg = TK.PhotoImage(data="""
        iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAYAAAA7M
        K6iAAAACXBIWXMAAAsTAAALEwEAmpwYAAABr0lEQV
        R4nO2WTytEYRTGf6NYqNkovgALhSVC2PgCLESKzMJ
        8AgsrhsTa36WPIRtKo8mYhY38icJQbJSw0qC3zq3T
        2zt33Ds3lHnq1O2c5zzP7b7vPe8LFfwSaoBl4B74l
        MiE0DlU/UZrSbSLYkk1eHEQwjjj0Fn0a7gTUjfRoU
        c0836kgpBiERrHRLPwJ43TJda0A9gEjoFnCfO8AbS
        XWPP9MG/dCGw7Nowdu0AzEWEIeFPip8A40CQxAZyp+
        iswWK5pUq27iR2g1sGrlZrHMz1TYU1HLNMnoE7Vqy
        U81AlHmw8HNW2zPq+JFVUfk7r5rKMqv2r1GE5rEOO
        sY+NMqvqjyj+o/KSj7yiIcc4hkFD1S5W/UPlEucZd
        wLslsG6NwZyE4XrYsHqMRmcQY4MZS8QMi3offoNwd
        M80IbFpCe0BcQcvLjXNXaPMOZsCPpTglfzfLRJJyX
        l1w52Nat73ASffGJmG00vIm0O2SK0KGAC2ZCe/SJx
        LbkA4LuRK3WT+33mcF5JZ06jQL5q3fqQFx4YptuZB
        J1/Kr8GcNvPAjWoIc3NIq/5rYM46ySrgx/AF1c/MR
        0JNiqgAAAAASUVORK5CYII=
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