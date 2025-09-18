"""Panels for tkCassiopee."""

import tkinter as TK
try: import Tk as CTK
except ImportError: from . import Tk as CTK
import Converter
import Converter.Internal as Internal
import Converter.PyTree as C
from . import PyTree as CPlot
from . import Ttk as TTK
from . import iconics
import time

#==============================================================================
# LoadPanel
LOADPANEL = None
# RenderPanel
RENDERPANEL = None
# ErrorWindow
ERRORWINDOW = None
# Materials (render)
MATERIALS = ['Solid', 'Flat', 'Glass', 'Chrome',
             'Metal', 'Wood', 'Marble', 'Granite', 'Brick', 'XRay',
             'Cloud', 'Gooch', 'Sphere', 'Texmat']
# global VARS for RENDER panel
VARS = []
# global VARS for ACTIVATION panel
AVARS = []
# global VARS for OPEN panel
OVARS = []
# data for  panel
mailData = {}
# data for Document panel
docData = {}
# Widgets dict
WIDGETS = {}

#==============================================================================
# Affiche le about et les informations d'install
#==============================================================================
def about():
    """About panel."""
    import KCore.buildInfo
    logoImg = TK.PhotoImage(data="""
R0lGODlhLAFWAOf8AAoGDAkLBw8MEBIQFBUQDhESDxITGhUTFhcVGBUWHRYXFRoWHhgXGhoYGxwa
IhocGR0bHiEfIh8hHiUjJiknKi4sLzAxLzMxNDY1Mzc1OTo4PDo8OT48QABLpABMnkFDQURCRklH
SgBWmgBVpkhKRwBXogBXsABYqgBZpBRTpk1LThRUoBVVnABbnwBbqwlbmgZbpwBdpwddoghcsABf
qgBfsVJRVABipxFeqgBjrhJgpVdUWABmqwBnphVioRRjnFVXVFpYXBxkqgRsqw1qr11bXhZpsyZq
q19hXhhus2JgZCRspRRxtwd1uhxxsGZkaGRmYzBwpC5yrGhpZyR2tmtsaih6szB3vG9scDZ4sy16
u29wbTJ9sC6Aw0F+p3N1cjWCvHh1eXZ4dUCCtkSCvkuBtziHu0SFunp8eX99gFaEtUiKuH+BfkyJ
xUuLwVKKuoOFgleMw06QvoiGilqQv4iJhlKWxYyOi4+NkWGYyY+RjpKQlGiYxGGbxWmcwZOVkpaU
mG+bzpqYnJialm+ix3Ghzp6coHmjyZyem22n0YKl0KCin3eq0KShpoCq0X+rzKSmo6ulpJGp1Kmn
q6mrqISy0oew3pSw1I+z1a2vrLCuso+23rGzsJK71aC31rS2s5660ZG93py82Lq4vLi6t6HC35rE
4by+u7+9wa7B2qnC4qHG26nF3MHDwMTCxqvK6LHN5bzL38jKx7jM5czJzbjO3rbT7L7S68/RztLQ
1L3V58fT7sXV6MnV49PV0sTY5NfZ1sPd8NvZ3cnd983d8czf5dHd7Mvf69bd49vd2tnd7djh6dHi
9djg9dzg8N/h3uLg5Nfj8dHl8uDl6Nfo7+Tm4+Tl7+fl6d/n8Ojq59zu++Tt9e3r7+zs9+ju8Ont
/e3v6/Dy7/Px9ev0/eX2/en39/P18vT0//H2+fj1+u/59Pz29fb49f/3/f/58fb7/vn7+P36//r8
+fT++f/87f/7+fL///v9+vj/9Pb/+//9+//+9f/9//n///z/+/7//CwAAAAALAFWAAAI/gD/CRxI
sKDBgwgTKlzIsKHDhxAjSpxIsaLFixgzatzIsaPHjyBDihxJsqTJkyhTqlzJsqXLlzBjypxJs6bN
mzhz6tzJs6fPn0CDCh1KtKjPewiR/pOX1CE8h0oNRh3ILx5Bpv/uIY36biBWhk8PKp36cKpVskbT
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
ClbQghfEYAY1uEEOdtCDHwRhCEUIloAAADs=
""")
    winl = TTK.Toplevel(border=0)
    winl.title('About Cassiopee')
    winl.columnconfigure(0, weight=1)
    winl.rowconfigure(0, weight=1)
    winl.minsize(320, 300)
    # position de la fenetre parent
    xpos = winl.master.winfo_rootx()+45
    ypos = winl.master.winfo_rooty()+45
    winl.geometry("%+d%+d" % (xpos, ypos))
    scrollbar = TTK.Scrollbar(winl, orient=TK.VERTICAL, width=10)
    scrollbar.grid(sticky=TK.NSEW, row=0, column=1)

    textWidget = TK.Text(winl, yscrollcommand=scrollbar.set,
                         width=43, height=20, background='White')
    textWidget.tag_config('title', justify=TK.CENTER, foreground='blue')
    textWidget.grid(sticky=TK.NSEW, row=0, column=0)
    scrollbar.config(command=textWidget.yview)
    myText = "- Cassiopee %s - \n"%Converter.__version__
    textWidget.insert(TK.END, myText, 'title')
    myText = "A CFD pre- and post-processing tool"
    textWidget.insert(TK.END, myText, 'title')
    textWidget.insert(TK.END, '\n\n')

    textWidget.image_create(TK.INSERT, image=logoImg)
    textWidget.image = logoImg

    myText = "\n\n Licensed under GPL3.\n\n"
    textWidget.insert(TK.END, myText)
    authors = KCore.__allAuthors__
    authors = authors.split(',')
    myText = " -- Authors: --\n "
    col = 0; naut = len(authors)
    for i in authors:
        col += 1
        if col > 0: myText += '\n'; col = 0
        if naut > 1: myText += i+','
        else: myText += i+'.\n\n'
        naut -= 1

    textWidget.insert(TK.END, myText)
    textWidget.insert(TK.END, "\n")
    buildInfo = KCore.buildInfo.buildDict
    myText = " Built on %s.\n"%buildInfo['date']
    textWidget.insert(TK.END, myText)
    myText = " Python: %s.\n"%buildInfo['python']
    textWidget.insert(TK.END, myText)
    myText = " Numpy: %s.\n"%buildInfo['numpy']
    textWidget.insert(TK.END, myText)

    myText = " Converter: libhdf5: "
    if 'hdf' in buildInfo and buildInfo['hdf'] != "None":
        myText += 'present.\n'
    else: myText += 'not present (HDF format unavailable).\n'
    textWidget.insert(TK.END, myText)

    #myText = " CPlot: png: "
    #if 'png' in buildInfo and buildInfo['png'] != "None":
    #    myText += 'present.\n'
    #else: myText += 'not present (png export unavailable).\n'
    #textWidget.insert(TK.END, myText)
    myText = " CPlot: mpeg: "
    if 'mpeg' in buildInfo and buildInfo['mpeg'] != "None":
        myText += 'present.\n'
    else: myText += 'not present (mpeg export unavailable).\n'
    textWidget.insert(TK.END, myText)

    myText = " Connector: libmpi: "
    if 'mpi' in buildInfo and buildInfo['mpi'] != "None":
        myText += 'present.\n'
    else: myText += 'not present (direct MPI communications unavailable).\n'
    textWidget.insert(TK.END, myText)

    myText = "\n\n"
    myText += " Cassiopee uses third party software:\n"
    myText += " - freeglut (see CPlot/GLUT/LICENSE).\n"
    myText += " - glew (see CPlot/GLEW/LICENSE).\n"
    myText += " - netgen (see Generator/Netgen/LICENSE).\n"
    myText += " - tetgen (see Generator/Tetgen/LICENSE).\n"
    myText += " - metis (see KCore/Metis/LICENSE).\n"
    myText += " - scotch (see XCore/scotch/LICENSE).\n"
    myText += " - MMGs (see Generator/MMGS/LICENSE).\n"
    myText += " - opencascade (see opencascade/LICENSE_LGPL_21.txt).\n"
    myText += " - xatlas (see Geom/xatlas/LICENCE).\n"
    myText += " - zstd (see Compressor/zstd/LICENCE).\n"
    myText += " - fpc (see Compressor/fpc.cpp).\n"
    myText += " - zfp (see Compressor/zfp/LICENCE).\n"
    myText += " - png (see KCore/Images/png/LICENCE).\n"
    myText += " - libjpeg (see KCore/Images/libjpeg/README).\n"
    myText += " - Icons by Icons8.\n"

    textWidget.insert(TK.END, myText)
    return

#==============================================================================
# Affiche la fenetre pour la cle d'activation
# Ecrit simplement un fichier avec la cle
#==============================================================================
def activation():
    global AVARS
    AVARS = []
    winl = TTK.Toplevel(border=0)
    winl.title('Activation key')
    winl.columnconfigure(0, weight=1)
    winl.columnconfigure(1, weight=1)
    winl.rowconfigure(0, weight=1)
    AVARS.append(winl)
    # position de la fenetre parent
    xpos = winl.master.winfo_rootx()+45
    ypos = winl.master.winfo_rooty()+45
    winl.geometry("%+d%+d" % (xpos, ypos))
    textWidget = TK.Text(winl, width=50, height=6, background='White')
    d = checkKey()
    if d == 0: version = 'open source'
    else:
        year = d/12 ; month = d - 12*year
        version = 'full (until %2d/%4d)'%(month,year)
    myText = 'You are using the %s version of Cassiopee.\n\n'%version
    myText += 'To switch to full version, you must request an activation key.\n'
    myText += 'Please mail to cassiopee@onera.fr or enter\n'
    myText += 'an activation key below:'
    textWidget.insert(TK.END, myText)
    textWidget.grid(sticky=TK.NSEW, row=0, column=0, columnspan=2)
    V = TK.StringVar(winl); V.set('')
    AVARS.append(V)
    entry = TK.Entry(winl, textvariable=V, background='White')
    entry.grid(row=1, column=0, columnspan=2, sticky=TK.EW)
    B = TTK.Button(winl, text="Submit", command=submitKey)
    BB = CTK.infoBulle(parent=B, text='Submit key.')
    B.grid(row=2, column=0, columnspan=1, sticky=TK.EW)
    B = TTK.Button(winl, text="Cancel", command=cancelKey)
    B.grid(row=2, column=1, columnspan=1, sticky=TK.EW)

def checkKey():
    import KCore
    return KCore.kcore.activation('0')

def readKeyFile(file):
    d = {}
    try:
        f = open(file, 'r')
        ret = f.readlines()
        f.close()
        for i in ret:
            i = i.replace(' ','')
            i = i.replace('\n', '')
            s = i.split(':')
            if len(s) == 2: d[s[0]] = s[1]
            if len(s) == 1: d['0'] = s[0]
    except: pass
    return d

def submitKey(event=None):
    key = AVARS[1].get()
    key = key.split(':')
    if len(key) == 2: name = key[0]; key = key[1]
    else: name = '0'; key = key[0]

    import KCore.installPath
    path = KCore.installPath.libPath
    # Essai dans installPath/.CassiopeKey
    file = path+'/.CassiopeeKey'
    d = readKeyFile(file)
    d[name] = key
    fail = False
    try:
        f = open(path+'/.CassiopeeKey', 'w')
        for k in d:
            f.write(k+':'+d[k]+'\n')
        f.close()
        CTK.TXT.insert('START', 'Key submitted.\n')
    except: fail = True

    if not fail:
        AVARS[0].destroy()
        return

    # Essai dans home/.CassiopeeKey
    import os.path
    path = os.path.expanduser('~')
    file = path+'/.CassiopeeKey'
    d = readKeyFile(file)
    d[name] = key
    try:
        f = open(path+'/.CassiopeeKey', 'w')
        for k in d:
            f.write(k+':'+d[k]+'\n')
        f.close()
        CTK.TXT.insert('START', 'Key submitted.\n')
    except:
        CTK.TXT.insert('START', 'Can not write key (check permission).\n')
        CTK.TXT.insert('START', 'Error: ', 'Error') ; return
    AVARS[0].destroy()

def cancelKey(event=None):
    AVARS[0].destroy()

#==============================================================================
# Error panel
#==============================================================================
def _deleteErrorWindow():
    try: ERRORWINDOW.destroy()
    except: pass

def _destroyErrorWindow(event):
    global ERRORWINDOW
    ERRORWINDOW = None

#==============================================================================
# IN: errors: [A,B,C,D,...]: display B, C,
# A is intended to be error number or garbage
#==============================================================================
def displayErrors(errors, header=''):
    if len(errors) == 0: return
    global ERRORWINDOW
    if ERRORWINDOW is not None:
        try: ERRORWINDOW.withdraw()
        except: ERRORWINDOW = None

    if ERRORWINDOW is None:
        ERRORWINDOW = TTK.Toplevel()
        ERRORWINDOW.columnconfigure(0, weight=1)
        ERRORWINDOW.rowconfigure(0, weight=1)
        ERRORWINDOW.title("Errors...")
        # position de la fenetre parent
        xpos = ERRORWINDOW.master.winfo_rootx()+45
        ypos = ERRORWINDOW.master.winfo_rooty()+45
        ERRORWINDOW.geometry("%+d%+d" % (xpos, ypos))
        #ERRORWINDOW.protocol("WM_DELETE_WINDOW", _deleteErrorWindow)
        #ERRORWINDOW.bind("<Destroy>", _destroyErrorWindow)
        scrollbar = TTK.Scrollbar(ERRORWINDOW, orient=TK.VERTICAL, width=10)
        scrollbar.grid(sticky=TK.NSEW, row=0, column=1)
        myText = TK.Text(ERRORWINDOW, yscrollcommand=scrollbar.set,
                         width=40, height=20, background='white')
        myText.tag_config("Title", foreground="blue")
        myText.mark_set('START', TK.INSERT)
        myText.mark_gravity('START', TK.LEFT)
        myText.grid(sticky=TK.NSEW, row=0, column=0)
        scrollbar.config(command=myText.yview)
    else:
        # trick pour avoir la fenetre d'erreur au premier plan
        ERRORWINDOW.withdraw(); ERRORWINDOW.deiconify(); ERRORWINDOW.focus_set()
        myText = ERRORWINDOW.winfo_children()[1] # text
        # myText.delete(1.0, TK.END)
    # Errors
    nerr = len(errors)//2; allText = ''
    for l in range(nerr): allText += ' - '+errors[2*l+1]+'\n'
    allText += '\n'
    myText.delete('1.0', TK.END) # clear previous errors
    myText.insert('START', allText)
    # Header
    ti = time.localtime()
    sti = '['+str(ti.tm_hour)+':'+str(ti.tm_min)+':'+str(ti.tm_sec)+']'
    myText.insert('START', '* %s %s\n'%(header,sti), 'Title')
    # Keep it to a reasonable size
    myText.delete(2000.0, TK.END)

#=============================================================================
# Mail panel
#=============================================================================
def destroyMailWindow(event=None):
    mailData['mailWindow'].destroy()

def cancelMail(event=None):
    global mailData
    mailData['mailFriends'] = None
    mailData['mailWindow'].destroy()
    mailData.pop('mailWindow')

def createMail():
    try:
        import smtplib
        from email.mime.image import MIMEImage
        from email.mime.text import MIMEText
        from email.mime.multipart import MIMEMultipart
    except:
        CTK.TXT.insert('START', 'Mail module unavailable.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    import os, os.path

    # Dump
    try: exportResolution = CTK.PREFS['exportResolution']
    except: exportResolution = '1280x720'
    CPlot.setState(exportResolution=exportResolution)
    CPlot.setState(export='.tmp001203.png')

    friendText = mailData['mailFriends']
    if friendText is None:
        os.remove('.tmp001203.png')
        return # cancelMail
    me = mailData['mailMe']
    bugReport = mailData['bugReport']
    if bugReport: friends = ['cassiopee@onera.fr']
    else: friends = friendText.split(';')
    titleText = mailData['titleText']
    if bugReport: titleText = '[BUG]'+titleText
    messageText = mailData['messageText']
    # ajoute la location si possible
    if CTK.FILE != '':
        p = os.path.abspath(CTK.FILE)
        messageText += '\nLink: %s\n'%p

    # Create the container (outer) email message.
    msg = MIMEMultipart()
    msg['Subject'] = '[Cassiopee] '+titleText
    msg['From'] = me
    COMMASPACE = ', '
    msg['To'] = COMMASPACE.join(friends)
    msg.preamble = 'Send by Cassiopee.'
    if messageText != '':
        msg.attach(MIMEText(messageText, 'plain'))

    # Wait for dump
    p = False
    while not p:
        p = os.path.exists('.tmp001203.png')
        if p: time.sleep(0.5); p = os.path.exists('.tmp001203.png')

    # Assume we know that the image files are all in PNG format
    fp = open('.tmp001203.png', 'rb')
    img = MIMEImage(fp.read())
    fp.close()
    msg.attach(img)

    # Send the email via our own SMTP server.
    try:
        s = smtplib.SMTP('localhost')
    except:
        os.remove('.tmp001203.png')
        CTK.TXT.insert('START', 'No valid SMTP server on your machine.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
        # essai de passer par gmail
        #s = smtplib.SMTP('smtp.gmail.com', 587)
        #s = smtplib.SMTP('mailhost.onera')

    s.sendmail(me, friends, msg.as_string())
    s.quit()
    os.remove('.tmp001203.png')
    CTK.TXT.insert('START', 'Mail sends.\n')

def sendMail(event=None):
    global mailData
    if CTK.t == []: return
    mailData['mailMe'] = mailData['meWidget'].get("1.0", TK.END)
    mailData['mailFriends'] = mailData['friendWidget'].get("1.0", TK.END)
    txt = mailData['titleWidget'].get("1.0", TK.END)
    #mailData['titleText'] = txt.encode('utf-8')
    mailData['titleText'] = txt
    txt = mailData['messageWidget'].get("1.0", TK.END)
    #mailData['messageText'] = txt.encode('utf-8')
    mailData['messageText'] = txt
    mailData['bugReport'] = False
    CTK.PREFS['mailMe'] = mailData['mailMe']
    if 'mailFriends' in CTK.PREFS:
        CTK.PREFS['mailFriends1'] = CTK.PREFS['mailFriends']
    CTK.PREFS['mailFriends'] = mailData['mailFriends']
    CTK.PREFS['mailTitle'] = mailData['titleText']
    CTK.savePrefFile()
    createMail()
    #mailData['mailWindow'].destroy()

def sendBug(event=None):
    global mailData
    if CTK.t == []: return
    mailData['mailMe'] = mailData['meWidget'].get("1.0", TK.END)
    mailData['mailFriends'] = mailData['friendWidget'].get("1.0", TK.END)
    txt = mailData['titleWidget'].get("1.0", TK.END)
    #mailData['titleText'] = txt.encode('utf-8')
    mailData['titleText'] = txt
    txt = mailData['messageWidget'].get("1.0", TK.END)
    #mailData['messageText'] = txt.encode('utf-8')
    mailData['messageText'] = txt
    mailData['bugReport'] = True
    CTK.PREFS['mailMe'] = mailData['mailMe']
    CTK.PREFS['mailFriends'] = mailData['mailFriends']
    CTK.PREFS['mailTitle'] = mailData['titleText']
    CTK.savePrefFile()
    createMail()
    #mailData['mailWindow'].destroy()

def openMailWindow():
    global mailData
    if 'mailWindow' in mailData:
        try: mailData['mailWindow'].withdraw()
        except: mailData.pop('mailWindow')

    if 'mailWindow' not in mailData:
        MAILWINDOW = TTK.Toplevel()
        mailData['mailWindow'] = MAILWINDOW
        MAILWINDOW.columnconfigure(0, weight=0)
        MAILWINDOW.columnconfigure(1, weight=1)
        MAILWINDOW.columnconfigure(2, weight=1)
        MAILWINDOW.rowconfigure(0, weight=0)
        MAILWINDOW.rowconfigure(1, weight=0)
        MAILWINDOW.rowconfigure(2, weight=0)
        MAILWINDOW.rowconfigure(3, weight=1)
        MAILWINDOW.rowconfigure(4, weight=0)
        MAILWINDOW.title("Mail image...")
        # position de la fenetre parent
        xpos = MAILWINDOW.master.winfo_rootx()+45
        ypos = MAILWINDOW.master.winfo_rooty()+45
        MAILWINDOW.geometry("%+d%+d" % (xpos, ypos))
        #MAILWINDOW.protocol("WM_DELETE_WINDOW", cancelMail)
        #MAILWINDOW.bind("<Destroy>", destroyMailWindow)
        B = TK.Label(MAILWINDOW, text="From:")
        B.grid(row=0, column=0, sticky=TK.EW)
        B = TK.Text(MAILWINDOW, width=40, height=1, background='White')
        BB = CTK.infoBulle(parent=B, text='Your email-adress: toto@tata.com')
        B.grid(row=0, column=1, columnspan=2, sticky=TK.EW)
        if 'mailMe' in CTK.PREFS: B.insert(TK.END, CTK.PREFS['mailMe'])
        mailData['meWidget'] = B
        B = TK.Label(MAILWINDOW, text="To:")
        B.grid(row=1, column=0, sticky=TK.EW)
        B = TK.Text(MAILWINDOW, width=40, height=2, background='White')
        BB = CTK.infoBulle(parent=B, text='Destinaries email-adress: toto@tata.com; titi@tata.com\nCassiopee support: cassiopee@onera.fr')
        B.grid(row=1, column=1, columnspan=2, sticky=TK.EW)
        if 'mailFriends' in CTK.PREFS: B.insert(TK.END, CTK.PREFS['mailFriends'])
        mailData['friendWidget'] = B
        B = TK.Label(MAILWINDOW, text="Title:")
        B.grid(row=2, column=0, sticky=TK.EW)
        B = TK.Text(MAILWINDOW, width=40, height=1, background='White')
        B.grid(row=2, column=1, columnspan=2, sticky=TK.EW)
        if 'mailTitle' in CTK.PREFS: B.insert(TK.END, CTK.PREFS['mailTitle'])
        mailData['titleWidget'] = B
        B = TK.Label(MAILWINDOW, text="Message:")
        B.grid(row=3, column=0, sticky=TK.EW)
        B = TK.Text(MAILWINDOW, width=40, height=5, background='White')
        B.grid(row=3, column=1, columnspan=2, sticky=TK.EW)
        mailData['messageWidget'] = B
        B = TTK.Button(MAILWINDOW, text="Cancel", command=cancelMail)
        B.grid(row=4, column=0, sticky=TK.EW)
        B = TTK.Button(MAILWINDOW, text="E-mail image", command=sendMail)
        BB = CTK.infoBulle(parent=B, text='Send an e-mail with your current screenshot.')
        B.grid(row=4, column=1, sticky=TK.EW)
        B = TTK.Button(MAILWINDOW, text="Report bug", command=sendBug)
        BB = CTK.infoBulle(parent=B, text='Send an e-mail to cassiopee team with \ninformations on the software state.')
        B.grid(row=4, column=2, sticky=TK.EW)
    else:
        MAILWINDOW = mailData['mailWindow']
        MAILWINDOW.withdraw(); MAILWINDOW.deiconify(); MAILWINDOW.focus_set()

#==============================================================================
# Document panel
#==============================================================================
def destroyDocumentWindow(event=None):
    docData['docWindow'].destroy()

def cancelDocument(event=None):
    global docData
    docData['docText'] = None
    docData['docWindow'].destroy()
    docData.pop('docWindow')

def createDoc():
    try:
        from odf.opendocument import OpenDocumentText, load
        from odf import text
        from odf.text import P
        from odf.draw import Frame, Image
        from odf.style import Style, MasterPage, PageLayout, PageLayoutProperties, TextProperties, GraphicProperties, ParagraphProperties, DrawingPageProperties
    except:
        CTK.TXT.insert('START', 'odf module unavailable.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    import os, os.path

    # Dump
    try: exportResolution = CTK.PREFS['exportResolution']
    except: exportResolution = '1920x1080'
    CPlot.setState(exportResolution=exportResolution)
    CPlot.setState(export='.tmp001204.png')

    text = docData['docText']
    if text is None:
        os.remove('.tmp001204.png')
        return # cancel

    blog = docData['blog'] # blog or not blog

    # Wait for dump (main image)
    p = False
    while not p:
        p = os.path.exists('.tmp001204.png')
        if p: time.sleep(0.5); p = os.path.exists('.tmp001204.png')

    # Dump 1D plot
    wrotePlot = False; sizePlot = None
    if CTK.TKPLOTXY is not None:
        wrotePlot = CTK.TKPLOTXY.DESKTOP.export('.tmp001205.png')
        if wrotePlot:
            sizePlot = CTK.TKPLOTXY.DESKTOP.getActiveGraphFigSize()

    # Write document
    docName = docData['docName']
    docName = docName.strip()
    if docName[-4:] == '.odt': docName = docName[0:-4]
    docNameOdt = docName+'.odt'
    if os.path.exists(docNameOdt): doc = load(docNameOdt)
    else: doc = OpenDocumentText()

    btext = '===============================================================\n'
    if blog:
        import getpass
        userName = getpass.getuser()
        date = time.strftime("%H:%M:%S (%d/%m/%Y)")
        btext += 'Blogged by %s at %s.\n'%(userName, date)
        if CTK.FILE != '':
            p = os.path.abspath(CTK.FILE)
            btext += '\nLink: %s\n'%p

    text = btext+text
    text = text.split('\n')
    for tx in text:
        p = P(text=tx); doc.text.addElement(p)

    # Main image
    (w, h) = CPlot.getState('win')
    coeff = h*1./w
    imgFrame = Frame(width="%fcm"%15., height="%fcm"%(15.*coeff), x="%fcm"%0.1, y="56pt")
    doc.text.addElement(imgFrame)
    href = doc.addPicture('.tmp001204.png')
    imgFrame.addElement(Image(href=href))

    # Plot image (if any)
    if wrotePlot:
        coeff = sizePlot[1]*1./sizePlot[0]
        imgFrame = Frame(width="%fcm"%15., height="%fcm"%(15.*coeff), x="%fcm"%0.1, y="56pt")
        doc.text.addElement(imgFrame)
        href = doc.addPicture('.tmp001205.png')
        imgFrame.addElement(Image(href=href))

    p = P()
    p.addElement(imgFrame)
    doc.text.addElement(p)

    doc.save(docName, True)
    os.remove('.tmp001204.png')
    CTK.TXT.insert('START', 'Image added to %s.\n'%(docNameOdt))

def writeDocument(event=None):
    global docData
    docData['docName'] = docData['docWidget'].get("1.0", TK.END)
    docData['docText'] = docData['textWidget'].get("1.0", TK.END)
    docData['blog'] = False
    CTK.PREFS['docName'] = docData['docName']
    #s = docData['docText'].encode('base64', 'strict')
    s = docData['docText']
    s = s.split('\n'); s = "".join(s)
    CTK.PREFS['docText'] = s
    CTK.savePrefFile()
    createDoc()
    #docData['docWindow'].destroy()

def writeBlog(event=None):
    global docData
    docData['docName'] = docData['docWidget'].get("1.0", TK.END)
    docData['docText'] = docData['textWidget'].get("1.0", TK.END)
    docData['blog'] = True
    CTK.PREFS['docName'] = docData['docName']
    #s = docData['docText'].encode('base64', 'strict')
    s = docData['docText']
    s = s.split('\n'); s = "".join(s)
    CTK.PREFS['docText'] = s
    CTK.savePrefFile()
    createDoc()
    #docData['docWindow'].destroy()

def openDocFile(event=None):
    import tkinter.filedialog as tkFileDialog
    initFile = docData['docWidget'].get("1.0", TK.END)
    initFile = ''
    file = tkFileDialog.asksaveasfilename(
        filetypes=[('Open document', '*.odt'), ('All files', '*.*')], initialfile=initFile)
    if (file == '' or file is None or file == ()): # user cancel
        return
    file = CTK.fixFileString2__(file)
    docData['docWidget'].delete("1.0", TK.END)
    docData['docWidget'].insert(TK.END, file)

def openDocWindow():
    global docData
    if 'docWindow' in docData:
        try: docData['docWindow'].withdraw()
        except: docData.pop('docWindow')
    if 'docWindow' not in docData:
        DOCWINDOW = TTK.Toplevel()
        docData['docWindow'] = DOCWINDOW
        DOCWINDOW.columnconfigure(0, weight=0)
        DOCWINDOW.columnconfigure(1, weight=1)
        DOCWINDOW.columnconfigure(2, weight=1)
        DOCWINDOW.rowconfigure(0, weight=0)
        DOCWINDOW.rowconfigure(1, weight=1)
        DOCWINDOW.rowconfigure(2, weight=0)
        DOCWINDOW.title("Save to document...")
        # position de la fenetre parent
        xpos = DOCWINDOW.master.winfo_rootx()+45
        ypos = DOCWINDOW.master.winfo_rooty()+45
        DOCWINDOW.geometry("%+d%+d" % (xpos, ypos))
        #DOCWINDOW.protocol("WM_DELETE_WINDOW", cancelDocument)
        #DOCWINDOW.bind("<Destroy>", destroyDocumentWindow)
        B = TTK.Button(DOCWINDOW, text="File:", command=openDocFile)
        BB = CTK.infoBulle(parent=B, text='Open an existing document file.')
        B.grid(row=0, column=0, sticky=TK.EW)
        B = TK.Text(DOCWINDOW, width=40, height=1, background='White')
        B.grid(row=0, column=1, columnspan=2, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Put here your document name (mydoc.odt).')
        if 'docName' in CTK.PREFS: B.insert(TK.END, CTK.PREFS['docName'])
        docData['docWidget'] = B
        B = TK.Label(DOCWINDOW, text="Text:")
        B.grid(row=1, column=0, sticky=TK.EW)
        B = TK.Text(DOCWINDOW, width=40, height=5, background='White')
        B.grid(row=1, column=1, columnspan=2, sticky=TK.EW)
        docData['textWidget'] = B
        if 'docText' in CTK.PREFS:
            #s = CTK.PREFS['docText'].decode('base64', 'strict')
            s = CTK.PREFS['docText']
            B.insert(TK.END, s)
        B = TTK.Button(DOCWINDOW, text="Cancel", command=cancelDocument)
        B.grid(row=2, column=0, sticky=TK.EW)
        B = TTK.Button(DOCWINDOW, text="Add image", command=writeDocument)
        BB = CTK.infoBulle(parent=B, text='Add current screenshot and text to \nyour document.')
        B.grid(row=2, column=1, sticky=TK.EW)
        B = TTK.Button(DOCWINDOW, text="Blog it", command=writeBlog)
        BB = CTK.infoBulle(parent=B, text='Add current screenshot and text to \nyour document with date and time.')
        B.grid(row=2, column=2, sticky=TK.EW)
    else:
        DOCWINDOW = docData['docWindow']
        DOCWINDOW.withdraw(); DOCWINDOW.deiconify(); DOCWINDOW.focus_set()

#==============================================================================
# Render panel
#==============================================================================
def _deleteRenderWindow():
    try: RENDERPANEL.destroy()
    except: pass

def _destroyRenderWindow(event):
    global RENDERPANEL
    RENDERPANEL = None

# update listboxes du render panel
def updateRenderPanel():
    if RENDERPANEL is None: return
    filter = VARS[9].get()
    selzones = set()
    listz = WIDGETS['myLists'][0]
    for l in WIDGETS['myLists']:
        cur = l.curselection()
        for s in cur:
            zname = listz.get(s)
            selzones.add(zname)
    for l in WIDGETS['myLists']: l.delete(0, TK.END)

    bases = Internal.getBases(CTK.t)
    for b in bases:
        zones = Internal.getNodesFromType1(b, 'Zone_t')
        for z in zones:
            name = b[0]+'/'+z[0]
            zoneName = name
            ri = Internal.getNodeFromName1(z, '.RenderInfo')
            if ri is not None:
                rt = Internal.getNodeFromName1(ri, 'Material')
                if rt is not None: material = Internal.getValue(rt)
                else: material = 'None'
                rt = Internal.getNodeFromName1(ri, 'Color')
                if rt is not None: color = Internal.getValue(rt)
                else: color = 'None'
                rt = Internal.getNodeFromName1(ri, 'Blending')
                if rt is not None: blending = "%5.2f"%(Internal.getValue(rt))
                else: blending = "%5.2f"%(1.)
                rt = Internal.getNodeFromName1(ri, 'MeshOverlay')
                if rt is not None: meshOverlay = str(Internal.getValue(rt))
                else: meshOverlay = '0'
                rt = Internal.getNodeFromName1(ri, 'ShaderParameters')
                if rt is not None:
                    v = rt[1]
                    param1 = "%5.2f"%(v[0])
                    param2 = "%5.2f"%(v[1])
                else:
                    param1 = "%5.2f"%(1.); param2 = "%5.2f"%(1.)
            else:
                material = 'None'; color = 'None'
                blending = "%5.2f"%(1.); meshOverlay = '0'
                param1 = "%5.2f"%(1.); param2 = "%5.2f"%(1.)

            data = [zoneName, material, color, blending, meshOverlay, param1, param2]
            ret = False
            if CTK.matchString(filter, zoneName):
                for i, d in enumerate(data):
                    WIDGETS['myLists'][i].insert(TK.END, d)

    for sel in selzones:
        for i in range(listz.size()):
            if sel == listz.get(i):
                listz.see(i)
                listz.selection_set(i)

# appele quand quelque chose est selectionne dans n'importe quelle listbox
def renderSelect(event=None):
    myset = set() # selected lines
    for l in WIDGETS['myLists']:
        sel = l.curselection()
        for i in sel: myset.add(i)

    # setter values
    material = None; color = None; blending = None;
    meshOverlay = None; shader1 = None; shader2 = None

    # select les zones dans CPlot
    CPlot.unselectAllZones()
    selected = []
    myList = WIDGETS['myLists'][0]
    for i in myset:
        name = myList.get(i)
        name = name.strip(); name = name.split('/')
        baseName = name[0]; zoneName = name[1]
        noz = CPlot.getCPlotNumber(CTK.t, baseName, zoneName)
        selected.append( (noz, 1) )
        if len(myset) == 1:
            z = Internal.getNodeFromPath(CTK.t, baseName+'/'+zoneName)
            ri = Internal.getNodeFromName1(z, '.RenderInfo')
            if ri is not None:
                #rt = Internal.getNodeFromName1(ri, 'Material')
                #if rt is not None: material = Internal.getValue(rt)
                #rt = Internal.getNodeFromName1(ri, 'Color')
                #if rt is not None: color = Internal.getValue(rt)
                #rt = Internal.getNodeFromName1(ri, 'Blending')
                #if rt is not None: blending = Internal.getValue(rt)
                rt = Internal.getNodeFromName1(ri, 'MeshOverlay')
                if rt is not None: meshOverlay = Internal.getValue(rt)
                #rt = Internal.getNodeFromName1(ri, 'ShaderParameters')
                #if rt is not None: shader1 = rt[1][0]; shader2 = rt[1][1]

    CPlot.setSelectedZones(selected)

    # set les datas dans les setters
    if len(myset) == 1: # uniquement si une seule zone selectionnee et pour le mesh overlay
        if meshOverlay == 0 or meshOverlay is None: VARS[3].set(0)
        elif meshOverlay == 1: VARS[3].set(1)
        #if blending is not None: WIDGETS['blending'].set(blending*100)
        #else: WIDGETS['blending'].set(100)
        #if shader1 is not None: WIDGETS['shader1'].set(shader1*50)
        #else: WIDGETS['shader1'].set(50)
        #if shader2 is not None: WIDGETS['shader2'].set(shader2*50)
        #else: WIDGETS['shader2'].set(50)

# called when double click on zone name (listbox 0)
def fitZone(i, event=None):
    myList = WIDGETS['myLists'][i]
    list0 = WIDGETS['myLists'][0]
    sel = myList.curselection()
    if len(sel) == 0: return
    sel = sel[0]
    name = list0.get(sel)
    name = name.strip(); name = name.split('/')
    baseName = name[0]; zoneName = name[1]
    noz = CPlot.getCPlotNumber(CTK.t, baseName, zoneName)
    selected = [(noz, 1)]
    CPlot.setSelectedZones(selected)
    CPlot.lookFor()

# called when right click on zone name (listbox 0)
def deactivateZone(event=None):
    myList = WIDGETS['myLists'][0]
    sel = myList.curselection()
    activated = []
    for s in sel:
        name = myList.get(s)
        name = name.strip(); name = name.split('/')
        baseName = name[0]; zoneName = name[1]
        noz = CPlot.getCPlotNumber(CTK.t, baseName, zoneName)
        sp = CPlot.getActiveStatus(noz)
        if sp: activated.append((noz, 0))
        else: activated.append((noz, 1))
    if activated:
        CPlot.setActiveZones(activated)

# called when double left click in color listbox
def openScalar(event=None):
    myList = WIDGETS['myLists'][2]
    cur = myList.curselection()
    varName = ''
    for c in cur:
        varName = myList.get(c); break
    if varName[0:4] == 'Iso:':
        # open tkView applet + focus + select render mode
        # open scalar in render widget
        # set the right scalar variable
        tkView = CTK.getModule('tkView')
        tkView.showApp()
        tkView.VARS[6].set('Render')
        tkView.setMode()
        tkView.WIDGETS['scalar'].grid(in_=tkView.WIDGETS['render'], row=1, column=0, columnspan=2, sticky=TK.EW)
        tkView.VARS[18].set(varName[4:])
        tkView.displayField()
        CPlot.setState(mode='render')

def selectAll(event=None):
    myList = RENDERPANEL.winfo_children()[1]
    myList.selection_set(0, TK.END)

def getSelection(event=None):
    updateRenderPanel() # pour forcer l'update
    myList = WIDGETS['myLists'][0]
    for l in WIDGETS['myLists']: l.selection_clear(0, TK.END)
    nzs = CPlot.getSelectedZones()
    if nzs == []: return
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        baseName = CTK.t[2][nob][0]
        name = baseName+'/'+CTK.t[2][nob][2][noz][0]
        for c in range(myList.size()):
            name2 = myList.get(c)
            if name2.strip() == name:
                myList.see(c)
                myList.selection_set(c)
                break

# reselect la list box i items from CPlot selected zones
def reselect(i):
    zList = WIDGETS['myLists'][0]
    myList = WIDGETS['myLists'][i]
    items = zList.get(0, TK.END)
    items = [it.strip() for it in items]
    nzs = CPlot.getSelectedZones()
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        base = CTK.t[2][nob]
        zone = CTK.t[2][nob][2][noz]
        zoneName = base[0]+'/'+zone[0]
        for c, item in enumerate(items):
            if item == zoneName:
                zList.see(c)
                zList.selection_set(c); break

def setColorVar(l):
    if l == 'Custom>':
        import tkinter.colorchooser as tkColorChooser
        ret = tkColorChooser.askcolor()
        l = ret[1]
    VARS[1].set(l)

def updateVarNameList(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ == 0 or nzs == []:
        vars = C.getVarNames(CTK.t, excludeXYZ=True)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        vars = C.getVarNames(CTK.t[2][nob][2][noz], excludeXYZ=True)
    m = WIDGETS['colors'].children['menu']
    m.delete(0, TK.END)
    allvars = ['White', 'Black', 'Grey', 'Blue', 'Red', 'Green', 'Yellow',
               'Orange', 'Brown', 'Magenta', 'Custom>']
    if len(vars) > 0:
        for v in vars[0]: allvars.append('Iso:'+v)
    for i in allvars:
        m.add_command(label=i, command=lambda v=VARS[1],l=i:setColorVar(l))

def updateVarNameList2(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        zvars = C.getVarNames(CTK.dt, excludeXYZ=True, mode=1)
    else:
        zvars = C.getVarNames(CTK.t, excludeXYZ=True, mode=1)
    allvars = ['White', 'Black', 'Grey', 'Blue', 'Red', 'Green', 'Yellow',
               'Orange', 'Brown', 'Magenta', 'Custom>']
    if len(zvars) > 0:
        for v in zvars[0]: allvars.append('Iso:'+v)

    if 'colors' in WIDGETS:
        WIDGETS['colors']['values'] = allvars

def setMaterial(event=None):
    nzs = CPlot.getSelectedZones()
    material = VARS[0].get()
    if nzs == []: return
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        a = CPlot.addRender2Zone(CTK.t[2][nob][2][noz], material=material)
        CTK.replace(CTK.t, nob, noz, a)
    CTK.TKTREE.updateApp()
    updateRenderPanel()
    CPlot.render()

def setColor(event=None):
    nzs = CPlot.getSelectedZones()
    color = VARS[1].get()
    if color == 'Custom>':
        import tkinter.colorchooser as tkColorChooser
        ret = tkColorChooser.askcolor()
        color = ret[1]
    VARS[1].set(color)

    if nzs == []: return
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        a = CPlot.addRender2Zone(CTK.t[2][nob][2][noz], color=color)
        CTK.replace(CTK.t, nob, noz, a)
    CTK.TKTREE.updateApp()
    updateRenderPanel()
    reselect(0)
    CPlot.render()

def setBlend(event=None):
    VARS[6].set('Blending [%.2f]'%(WIDGETS['blending'].get() / 100.))
    nzs = CPlot.getSelectedZones()
    blend = WIDGETS['blending'].get() / 100.
    blend = float(blend)
    if nzs == []: return
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        a = CPlot.addRender2Zone(CTK.t[2][nob][2][noz], blending=blend)
        CTK.replace(CTK.t, nob, noz, a)
    CTK.TKTREE.updateApp()
    updateRenderPanel()
    CPlot.render()

def setMesh(event=None):
    nzs = CPlot.getSelectedZones()
    if nzs == []: return
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        zone = CTK.t[2][nob][2][noz]
        ri = Internal.getNodeFromName1(zone, '.RenderInfo')
        meshOverlay = 0
        if ri is not None:
            rt = Internal.getNodeFromName1(ri, 'MeshOverlay')
            if rt is not None: meshOverlay = Internal.getValue(rt)
        if meshOverlay == 0: meshOverlay = 1
        else: meshOverlay = 0
        a = CPlot.addRender2Zone(CTK.t[2][nob][2][noz], meshOverlay=meshOverlay)
        CTK.replace(CTK.t, nob, noz, a)
    CTK.TKTREE.updateApp()
    updateRenderPanel()
    CPlot.render()

def setShaderParameter1(event=None):
    VARS[7].set('Shader1 [%.2f]'%(WIDGETS['shader1'].get() / 50.))
    nzs = CPlot.getSelectedZones()
    p1 = WIDGETS['shader1'].get() / 50.
    p1 = float(p1)
    if nzs == []: return
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        r = Internal.getNodeFromName(CTK.t[2][nob][2][noz], 'ShaderParameters')
        if r is not None: p2 = r[1][1]
        else: p2 = 1.
        a = CPlot.addRender2Zone(CTK.t[2][nob][2][noz], shaderParameters=[p1,p2])
        CTK.replace(CTK.t, nob, noz, a)
    CTK.TKTREE.updateApp()
    updateRenderPanel()
    CPlot.render()

def setShaderParameter2(event=None):
    VARS[8].set('Shader2 [%.2f]'%(WIDGETS['shader2'].get() / 50.))
    nzs = CPlot.getSelectedZones()
    p2 = WIDGETS['shader2'].get() / 50.
    p2 = float(p2)
    if nzs == []: return
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        r = Internal.getNodeFromName(CTK.t[2][nob][2][noz], 'ShaderParameters')
        if r is not None: p1 = r[1][0]
        else: p1 = 1.
        a = CPlot.addRender2Zone(CTK.t[2][nob][2][noz], shaderParameters=[p1,p2])
        CTK.replace(CTK.t, nob, noz, a)
    CTK.TKTREE.updateApp()
    updateRenderPanel()
    CPlot.render()

# fonction de scroll globale
def updateScrollLists(*args):
    for l in WIDGETS['myLists']:
        l.yview(*args)

# fonctions de scroll pour chaque listbox
# chacun se base sur i
def yscrolli(i, *args):
    l0 = WIDGETS['myLists'][i]
    for c, l in enumerate(WIDGETS['myLists']):
        if c != i:
            if l.yview() != l0.yview():
                l.yview_moveto(args[0])
    WIDGETS['scrollbar'].set(*args)

# called when filtering
def filterZones(event=None):
    updateRenderPanel()

def openRenderPanel():
    global RENDERPANEL
    if RENDERPANEL is not None:
        try: RENDERPANEL.withdraw()
        except: RENDERPANEL = None
    if RENDERPANEL is None:
        ttk = CTK.importTtk()
        RENDERPANEL = TTK.Toplevel(CTK.WIDGETS['masterWin'])
        RENDERPANEL.columnconfigure(0, weight=1) # zoneName
        RENDERPANEL.columnconfigure(1, weight=1) # material
        RENDERPANEL.columnconfigure(2, weight=1) # Color
        RENDERPANEL.columnconfigure(3, weight=0) # blend
        RENDERPANEL.columnconfigure(4, weight=0) # Mesh
        RENDERPANEL.columnconfigure(5, weight=0) # shader param1
        RENDERPANEL.columnconfigure(6, weight=0) # shader param2
        RENDERPANEL.columnconfigure(7, weight=0) # scrollbar
        RENDERPANEL.rowconfigure(0, weight=0)
        RENDERPANEL.rowconfigure(1, weight=0)
        RENDERPANEL.rowconfigure(2, weight=1)
        RENDERPANEL.rowconfigure(3, weight=0)

        RENDERPANEL.title("Render panel")
        # position de la fenetre parent
        xpos = RENDERPANEL.master.winfo_rootx()+45
        ypos = RENDERPANEL.master.winfo_rooty()+45
        RENDERPANEL.geometry("%+d%+d" % (xpos, ypos))
        #RENDERPANEL.protocol("WM_DELETE_WINDOW", _deleteRenderWindow)
        #RENDERPANEL.bind("<Destroy>", _destroyRenderWindow)
        scrollbar = TTK.Scrollbar(RENDERPANEL, orient=TK.VERTICAL, width=10)
        scrollbar.grid(sticky=TK.NSEW, row=2, column=7)
        WIDGETS['scrollbar'] = scrollbar
        myLists = []
        for i in range(7):
            if i == 0: # I cant find another way of setting constant in lambda functions
                # ZoneName listbox
                myList = TTK.Listbox(RENDERPANEL, selectmode=TK.EXTENDED,
                                     yscrollcommand=lambda *args: yscrolli(0,*args),
                                     #width=30, height=20,
                                     background='white', exportselection=1)
                myList.bind('<Double-Button>', lambda event: fitZone(0, event))
                myList.bind('<Button-3>', deactivateZone)

            elif i == 1:
                # Material listbox
                myList = TTK.Listbox(RENDERPANEL, selectmode=TK.EXTENDED,
                                     yscrollcommand=lambda *args: yscrolli(1,*args),
                                     #width=30, height=20,
                                     background='white', exportselection=1)
                myList.bind('<Double-Button>', lambda event: fitZone(1, event))
            elif i == 2:
                # Color listbox
                myList = TTK.Listbox(RENDERPANEL, selectmode=TK.EXTENDED,
                                     yscrollcommand=lambda *args: yscrolli(2,*args),
                                     #width=30, height=20,
                                     background='white', exportselection=1)
                myList.bind('<Double-Button>', openScalar)
            elif i == 3:
                # Blend listbox
                myList = TTK.Listbox(RENDERPANEL, selectmode=TK.EXTENDED,
                                     yscrollcommand=lambda *args: yscrolli(3,*args),
                                     width=4,
                                     background='white', exportselection=1)
                myList.bind('<Double-Button>', lambda event: fitZone(3, event))
            elif i == 4:
                # Mesh listbox
                myList = TTK.Listbox(RENDERPANEL, selectmode=TK.EXTENDED,
                                     yscrollcommand=lambda *args: yscrolli(4,*args),
                                     width=4,
                                     background='white', exportselection=1)
                myList.bind('<Double-Button>', lambda event: fitZone(4, event))
            elif i == 5:
                # Shader param1 listbox
                myList = TTK.Listbox(RENDERPANEL, selectmode=TK.EXTENDED,
                                     yscrollcommand=lambda *args: yscrolli(5,*args),
                                     width=4,
                                     background='white', exportselection=1)
                myList.bind('<Double-Button>', lambda event: fitZone(5, event))
            elif i == 6:
                # shader param2 listbox
                myList = TTK.Listbox(RENDERPANEL, selectmode=TK.EXTENDED,
                                     yscrollcommand=lambda *args: yscrolli(6,*args),
                                     width=4,
                                     background='white', exportselection=1)
                myList.bind('<Double-Button>', lambda event: fitZone(6, event))

            myList.grid(sticky=TK.NSEW, row=2, column=i)
            myList.bind('<<ListboxSelect>>', renderSelect)

            myLists.append(myList)
        WIDGETS['myLists'] = myLists
        scrollbar.config(command=updateScrollLists)

        # -0- Material
        V = TK.StringVar(RENDERPANEL); V.set('Solid'); VARS.append(V)
        # -1- Color
        V = TK.StringVar(RENDERPANEL); V.set('White'); VARS.append(V)
        # -2- Blend
        V = TK.StringVar(RENDERPANEL); V.set('1.'); VARS.append(V)
        # -3- Mesh
        V = TK.IntVar(RENDERPANEL); V.set('0'); VARS.append(V)
        # -4- Shader Param1
        V = TK.StringVar(RENDERPANEL); V.set('1.'); VARS.append(V)
        # -5- Shader Param2
        V = TK.StringVar(RENDERPANEL); V.set('1.'); VARS.append(V)
        # -6- Blending info bulle
        V = TK.StringVar(RENDERPANEL); V.set('Blending.'); VARS.append(V)
        # -7- Shader parameter1 info bulle
        V = TK.StringVar(RENDERPANEL); V.set('Shader1.'); VARS.append(V)
        # -8- Shader parameter2 info bulle
        V = TK.StringVar(RENDERPANEL); V.set('Shader2.'); VARS.append(V)
        # -9- Filter zone name string
        V = TK.StringVar(RENDERPANEL); V.set(''); VARS.append(V)

        # -- labels --
        label1 = TTK.Label(RENDERPANEL, text='Base/Zone')
        label1.grid(row=0, column=0, sticky=TK.EW)
        label2 = TTK.Label(RENDERPANEL, text='Material')
        label2.grid(row=0, column=1, sticky=TK.EW)
        label3 = TTK.Label(RENDERPANEL, text='Color')
        label3.grid(row=0, column=2, sticky=TK.EW)
        label4 = TTK.Label(RENDERPANEL, text='Blend')
        label4.grid(row=0, column=3, sticky=TK.EW)
        label5 = TTK.Label(RENDERPANEL, text='Mesh')
        label5.grid(row=0, column=4, sticky=TK.EW)
        label6 = TTK.Label(RENDERPANEL, text='Shader 1')
        label6.grid(row=0, column=5, sticky=TK.EW)
        label7 = TTK.Label(RENDERPANEL, text='Shader 2')
        label7.grid(row=0, column=6, sticky=TK.EW)

        # -- Filters --
        B = TK.Entry(RENDERPANEL, textvariable=VARS[9], background='White', width=40)
        B.bind('<KeyRelease>', filterZones)
        B.grid(row=3, column=0, columnspan=7, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Filter zones by this regexp.')

        # -- Setters --

        # zone get selection from CPlot
        B = TTK.Button(RENDERPANEL, command=getSelection, text='Get',
                       image=iconics.PHOTO[8], padx=0, pady=0)
        BB = CTK.infoBulle(parent=B, text='Get selection from plotter.')
        B.grid(row=1, column=0, sticky=TK.EW)

        # material setter
        B = TTK.OptionMenu(RENDERPANEL, VARS[0], *MATERIALS, command=setMaterial)
        B.grid(row=1, column=1, sticky=TK.EW)

        # color setter
        F = TTK.Frame(RENDERPANEL, borderwidth=0)
        F.columnconfigure(0, weight=1)
        if ttk is None:
            B = TK.OptionMenu(F, VARS[1], '', command=setColor)
            B.grid(sticky=TK.EW)
            F.bind('<Enter>', updateVarNameList)
            F.grid(row=1, column=2, sticky=TK.EW)
            WIDGETS['colors'] = B
        else:
            B = ttk.Combobox(F, textvariable=VARS[1],
                             values=[], state='readonly', height=11)
            B.bind('<<ComboboxSelected>>', setColor)
            B.grid(sticky=TK.EW)
            F.bind('<Enter>', updateVarNameList2)
            F.grid(row=1, column=2, sticky=TK.EW)
            WIDGETS['colors'] = B

        # blending setter
        B = TTK.Scale(RENDERPANEL, from_=0, to=100, orient=TK.HORIZONTAL,
                      command=setBlend, showvalue=0, borderwidth=1, value=100)
        WIDGETS['blending'] = B
        B.grid(row=1, column=3, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, textVariable=VARS[6])

        # Mesh setter (toggle)
        #B = TTK.Button(RENDERPANEL, text="Toggle", command=setMesh)
        B = TTK.Checkbutton(RENDERPANEL, variable=VARS[3], command=setMesh)
        B.grid(row=1, column=4, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text="Toggle mesh overlay")

        # shader param1 setter
        B = TTK.Scale(RENDERPANEL, from_=0, to=100, orient=TK.HORIZONTAL,
                      command=setShaderParameter1, showvalue=0, borderwidth=1, value=100)
        WIDGETS['shader1'] = B
        WIDGETS['shader1'].set(50)
        B.grid(row=1, column=5, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, textVariable=VARS[7])

        # shaderparam2 setter
        B = TTK.Scale(RENDERPANEL, from_=0, to=100, orient=TK.HORIZONTAL,
                      command=setShaderParameter2, showvalue=0, borderwidth=1, value=100)
        WIDGETS['shader2'] = B
        WIDGETS['shader2'].set(50)
        B.grid(row=1, column=6, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, textVariable=VARS[8])

    else:
        # trick pour avoir la fenetre au premier plan
        RENDERPANEL.withdraw(); RENDERPANEL.deiconify(); RENDERPANEL.focus_set()
    updateRenderPanel()

#====================================================================================
# Load panel: panel for partial load of files
#====================================================================================
def openLoadFileDialog(event=None):
    import tkinter.filedialog as tkFileDialog
    files = tkFileDialog.askopenfilenames(
        filetypes=[('CGNS files', '*.cgns'), ('CGNS files', '*.adf'), ('CGNS files', '*.hdf'), ('CGNS/ADF files', '*.adf'), ('CGNS/HDF files', '*.hdf'), ('All files', '*.*')], initialfile=CTK.FILE, multiple=0)
    if files == '' or files is None or files == (): # user cancel
        return
    files = CTK.fixFileString__(files, CTK.FILE)
    CTK.tkLoadFile(files, mode='partial')
    updateLoadPanel()

# IN: wname: string dans le widget
# OUT: retourne: le nom (/Base/Zone) et le nom tagge (/Base/Zone [X])
def ripTag(wname):
    l = len(wname)
    if wname[l-4:l] == ' [X]':
        pname = wname[0:l-4]
        tname = wname
    else:
        pname = wname
        tname = wname+ ' [X]'
    return pname, tname

# Get the number of the element in list equal to e or et
def getNumber(l, e, et):
    for i, li in enumerate(l):
        if li == e: return i
        if li == et: return i
    return -1

# Met a jour les listes en fonction du handle et de CTK.t (si deja loade)
def updateLoadPanel():
    if CTK.HANDLE is None: return
    OVARS[0].set(CTK.HANDLE.fileName)
    vars = CTK.HANDLE.fileVars
    znp = CTK.HANDLE.znp

    # VARS
    lb = WIDGETS['LBVARS']
    lb.delete(0, TK.END)
    for i, value in enumerate(vars):
        lb.insert(i, value)
        OVARS[3].append(value)

    # ZONES
    lb = WIDGETS['LBZONES']
    lb.delete(0, TK.END)
    for i, value in enumerate(znp):
        lb.insert(i, value)
        OVARS[4].append(value)

def loadVars(event=None):
    if CTK.HANDLE is None: return
    CTK.setCursor(2, WIDGETS['loadVars'])
    # Recupere les variables selectionnees
    selection = WIDGETS['LBVARS'].curselection()
    varList = []
    for s in selection:
        v = WIDGETS['LBVARS'].get(s)
        v, tname = ripTag(v)
        WIDGETS['LBVARS'].delete(s)
        WIDGETS['LBVARS'].insert(s, tname)
        i = getNumber(OVARS[3], v, tname)
        OVARS[3][i] = tname
        varList.append(v)
    # Load les variables donnees pour toutes zones existants dans t
    CTK.HANDLE._loadVariables(CTK.t, varList)

    CTK.t = CTK.upgradeTree(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    if 'tkContainers' in CTK.TKMODULES: CTK.TKMODULES['tkContainers'].updateApp()
    if CTK.TKPLOTXY is not None: CTK.TKPLOTXY.updateApp()
    CTK.display(CTK.t)
    CTK.setCursor(0, WIDGETS['loadVars'])

def unloadVars(event=None):
    if CTK.HANDLE is None: return
    # Recupere les variables selectionnees
    selection = WIDGETS['LBVARS'].curselection()
    varList = []
    for s in selection:
        v = WIDGETS['LBVARS'].get(s)
        v, tname = ripTag(v)
        WIDGETS['LBVARS'].delete(s)
        WIDGETS['LBVARS'].insert(s, v)
        i = getNumber(OVARS[3], v, tname)
        OVARS[3][i] = v
        varList.append(v)
    # Enleve les variables selectionnees pour toutes les zones existants dans t
    C._rmVars(CTK.t, varList)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    if CTK.TKPLOTXY is not None: CTK.TKPLOTXY.updateApp()
    CTK.display(CTK.t)

def filterVarList(event=None):
    filter = OVARS[1].get()
    listbox = WIDGETS['LBVARS']
    listbox.delete(0, TK.END)
    for s in OVARS[3]:
        if CTK.matchString(filter, s):
            listbox.insert(TK.END, s)
    listbox.config(yscrollcommand=WIDGETS['SCROLLVARS'].set)
    WIDGETS['SCROLLVARS'].config(command=listbox.yview)
    return True

def loadZones(event=None):
    if CTK.HANDLE is None: return
    CTK.setCursor(2, WIDGETS['loadZones'])
    import Converter.Filter as Filter
    # First load
    if len(Internal.getZones(CTK.t)) == 0: firstLoad = True
    else: firstLoad = False
    # Recupere les zones selectionnees
    selection = WIDGETS['LBZONES'].curselection()
    zList = []
    for s in selection:
        v = WIDGETS['LBZONES'].get(s)
        v, tname = ripTag(v)
        WIDGETS['LBZONES'].delete(s)
        WIDGETS['LBZONES'].insert(s, tname)
        i = getNumber(OVARS[4], v, tname)
        #v = v.encode('utf-8') # Cedre!!
        OVARS[4][i] = tname
        zList.append(v)
    # Charge les GC+GC+BC pour les zones selectionnees + variables deja dans t
    CTK.HANDLE._loadZonesWoVars(CTK.t, zList)
    vars = C.getVarNames(CTK.t, excludeXYZ=True)
    if len(vars)>1:
        vars = vars[0]
        if len(vars)>0:
            Filter._loadVariables(CTK.t, CTK.HANDLE.fileName, zList, vars, CTK.HANDLE.format)

    Filter._convert2PartialTree(CTK.t)
    CTK.t = CTK.upgradeTree(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    if CTK.TKPLOTXY is not None: CTK.TKPLOTXY.updateApp()
    CTK.display(CTK.t)
    if firstLoad: CPlot.fitView(); module = CTK.getModule('tkContainers'); module.updateApp()
    CTK.setCursor(0, WIDGETS['loadZones'])

def unloadZones(event=None):
    if CTK.HANDLE is None: return
    # Decharge les zones selectionnees
    selection = WIDGETS['LBZONES'].curselection()
    zList = []
    for s in selection:
        v = WIDGETS['LBZONES'].get(s)
        v, tname = ripTag(v)
        WIDGETS['LBZONES'].delete(s)
        WIDGETS['LBZONES'].insert(s, v)
        i = getNumber(OVARS[4], v, tname)
        #v = v.encode('utf-8') # Cedre!!
        OVARS[4][i] = v
        zList.append(v)
    for p in zList:
        Internal._rmNodeByPath(CTK.t, p)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    if CTK.TKPLOTXY is not None: CTK.TKPLOTXY.updateApp()
    CTK.display(CTK.t)

def filterZoneList(event=None):
    filter = OVARS[2].get()
    listbox = WIDGETS['LBZONES']
    listbox.delete(0, TK.END)
    for s in OVARS[4]:
        if CTK.matchString(filter, s):
            listbox.insert(TK.END, s)
    listbox.config(yscrollcommand=WIDGETS['SCROLLZONES'].set)
    WIDGETS['SCROLLZONES'].config(command=listbox.yview)
    return True

def openLoadPanel(event=None):
    global LOADPANEL
    if LOADPANEL is not None:
        try: LOADPANEL.withdraw()
        except: LOADPANEL = None
    if LOADPANEL is None:
        LOADPANEL = TTK.Toplevel()
        LOADPANEL.columnconfigure(0, weight=1)
        #LOADPANEL.columnconfigure(1, weight=0)
        #LOADPANEL.columnconfigure(2, weight=1)
        #LOADPANEL.columnconfigure(3, weight=0)
        LOADPANEL.rowconfigure(0, weight=1)
        LOADPANEL.title("File load panel")
        LOADPANEL.bind_all("<Control-o>", openLoadFileDialog)

        # position de la fenetre parent
        xpos = LOADPANEL.master.winfo_rootx()+45
        ypos = LOADPANEL.master.winfo_rooty()+45
        LOADPANEL.geometry("%+d%+d" % (xpos, ypos))
        F = TK.Frame(LOADPANEL, border=0)
        F.columnconfigure(0, weight=1)
        F.columnconfigure(1, weight=1)
        F.columnconfigure(2, weight=0)
        F.columnconfigure(3, weight=1)
        F.columnconfigure(4, weight=1)
        F.columnconfigure(5, weight=0)
        F.rowconfigure(0, weight=0)
        F.rowconfigure(1, weight=1)
        F.rowconfigure(2, weight=0)
        F.rowconfigure(3, weight=0)
        F.grid(row=0, column=0, columnspan=1, sticky=TK.NSEW)

        # -0- filename
        V = TK.StringVar(F); V.set(CTK.FILE); OVARS.append(V)
        # -1- Vars filter
        V = TK.StringVar(F); V.set(''); OVARS.append(V)
        # -2- Zone filter
        V = TK.StringVar(F); V.set(''); OVARS.append(V)
        # -3- Store Vars list
        OVARS.append([])
        # -4- Store Zone list
        OVARS.append([])

        # Entete
        B = TTK.Entry(F, textvariable=OVARS[0], background='White')
        BB = CTK.infoBulle(parent=B, text='File for reading.')
        B.grid(row=0, column=0, columnspan=5, sticky=TK.EW)
        B = TTK.Button(F, text='...', command=openLoadFileDialog)
        BB = CTK.infoBulle(parent=B, text='Select file.')
        B.grid(row=0, column=4, columnspan=1, sticky=TK.EW)

        # Zones
        scrollbar = TTK.Scrollbar(F, orient=TK.VERTICAL, width=10)
        scrollbar.grid(sticky=TK.NSEW, row=1, column=2)
        myList = TK.Listbox(F, selectmode=TK.EXTENDED,
                            yscrollcommand=scrollbar.set,
                            width=40, height=20, background='white')
        myList.grid(sticky=TK.NSEW, row=1, column=0, columnspan=2)
        scrollbar.config(command=myList.yview)
        WIDGETS['LBZONES'] = myList
        WIDGETS['SCROLLZONES'] = scrollbar

        B = TTK.Button(F, text="Load", command=loadZones)
        B.grid(row=2, column=0, sticky=TK.EW)
        WIDGETS['loadZones'] = B
        BB = CTK.infoBulle(parent=B, text='Load zones.')

        B = TTK.Button(F, text="Unload", command=unloadZones)
        B.grid(row=2, column=1, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Unload zones.')

        B = TK.Entry(F, textvariable=OVARS[2], background='White', width=40)
        B.bind('<KeyRelease>', filterZoneList)
        B.grid(row=3, column=0, columnspan=2, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Filter zones by this regexp.')

        # Vars
        scrollbar = TTK.Scrollbar(F, orient=TK.VERTICAL, width=10)
        scrollbar.grid(sticky=TK.NSEW, row=1, column=5)
        myList = TK.Listbox(F, selectmode=TK.EXTENDED,
                            yscrollcommand=scrollbar.set,
                            width=40, height=20, background='white')
        myList.grid(sticky=TK.NSEW, row=1, column=3, columnspan=2)
        scrollbar.config(command=myList.yview)
        WIDGETS['LBVARS'] = myList
        WIDGETS['SCROLLVARS'] = scrollbar

        B = TTK.Button(F, text="Load", command=loadVars)
        B.grid(row=2, column=3, sticky=TK.EW)
        WIDGETS['loadVars'] = B
        BB = CTK.infoBulle(parent=B, text='Load vars.')

        B = TTK.Button(F, text="Unload", command=unloadVars)
        B.grid(row=2, column=4, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Unload vars.')

        B = TK.Entry(F, textvariable=OVARS[1], background='White', width=40)
        B.bind('<KeyRelease>', filterVarList)
        B.grid(row=3, column=3, columnspan=2, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Filter vars by this regexp.')

    else:
        # trick pour avoir la fenetre au premier plan
        LOADPANEL.withdraw(); LOADPANEL.deiconify(); LOADPANEL.focus_set()
    updateLoadPanel()
