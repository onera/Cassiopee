# matplotlib decorators for CPlot
# strongly inspired by Luis Bernardos ideas

from . import PyTree as CPlot
from . import ColorMaps
import matplotlib.pyplot as plt
import matplotlib
import numpy

#import matplotlib.cbook as cbook
#import matplotlib.patches as patches
#import matplotlib.colors as mplcolors
#from matplotlib import cm
#from matplotlib.colors import ListedColormap, LinearSegmentedColormap

dpi = 100.

#==========================================================
def setBatch(batch=True):
    """Set batch mode for matplotlib."""
    if batch: matplotlib.use('Agg') # sans serveur X
    else: matplotlib.use('TkAgg') # avec serveur X

#==========================================================
def createColormap(type='Blue2Red'):
    """Create colormap."""
    if type == 'Blue2Red':
        colors = [(0.00,[0,0,1]),
                  (0.25,[0,1,1]),
                  (0.50,[0,1,0]),
                  (0.75,[1,1,0]),
                  (1.00,[1,0,0])]
    elif type == 'BiColorRGB':
        c1 = CPlot.getState('colormap1')
        c2 = CPlot.getState('colormap2')
        colors = [(0.00,[c1[0],c1[1],c1[2]]),
                  (1.00,[c2[0],c2[1],c2[2]])]
    elif type == 'BiColorHSV':
        c1 = CPlot.getState('colormap1')
        c2 = CPlot.getState('colormap2')
        colors = [(0.00,[c1[0],c1[1],c1[2]]),
                  (1.00,[c2[0],c2[1],c2[2]])]
    elif type == 'TriColorRGB':
        c1 = CPlot.getState('colormap1')
        c2 = CPlot.getState('colormap2')
        c3 = CPlot.getState('colormap3')
        colors = [(0.00,[c1[0],c1[1],c1[2]]),
                  (0.50,[c3[0],c3[1],c3[2]]),
                  (1.00,[c2[0],c2[1],c2[2]])]
    elif type == 'TriColorHSV':
        c1 = CPlot.getState('colormap1')
        c2 = CPlot.getState('colormap2')
        c3 = CPlot.getState('colormap3')
        colors = [(0.00,[c1[0],c1[1],c1[2]]),
                  (0.50,[c3[0],c3[1],c3[2]]),
                  (1.00,[c2[0],c2[1],c2[2]])]
    elif type == 'MultiColorRGB':
        colors = ColorMaps.export2MatplotLib2(CPlot.getState('colormapC'))
    elif type == 'MultiColorHSV':
        colors = ColorMaps.export2MatplotLib2(CPlot.getState('colormapC'))
    elif type == 'Diverging':
        colors = ColorMaps.export2MatplotLib2(ColorMaps.Diverging)
    else: # Blue2Red par defaut
        colors = [(0.00,[0,0,1]),
                  (0.25,[0,1,1]),
                  (0.50,[0,1,0]),
                  (0.75,[1,1,0]),
                  (1.00,[1,0,0])]

    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Blue2Red', colors)
    return cmap

#==========================================================
# IN: fig: the figure
# IN: xmin,xmax, ymin, ymax
# IN: levels: list or numpy of levels
# IN: cmap: color map
# IN: levels: level of colormaps
# IN: title: color bar title
# IN: discrete: True if discrete color map
#==========================================================
def createColorBar(fig, ax, levels=None, title=None, cmap=None, valueFormat='%0.3f', discrete=True,
                   fontSize=20, color="black", location="right", pad=0.):
    """Create a color bar."""
    if cmap is None: 
        cmap = CPlot.getState('colormap')
        if cmap == 0 or cmap == 1: # primaire
            cmap = createColormap('Blue2Red')
        elif cmap == 2 or cmap == 3: # primaire
            cmap = createColormap('BiColorRGB')
        elif cmap == 4 or cmap == 5: # primaire
            cmap = createColormap('BiColorHSV')
        elif cmap == 6 or cmap == 7: # primaire
            cmap = createColormap('TriColorRGB')
        elif cmap == 8 or cmap == 9: # primaire
            cmap = createColormap('TriColorHSV')
        elif cmap == 10 or cmap == 11: # primaire
            cmap = createColormap('MultiColorRGB')
        elif cmap == 12 or cmap == 13: # primaire
            cmap = createColormap('MultiColorHSV')
        elif cmap == 14 or cmap == 15: # primaire
            cmap = createColormap('Diverging')
        else: # default
            cmap = createColormap('Blue2Red')
    elif isinstance(cmap, str): cmap = createColormap(cmap)

    if levels is None:
        levels = CPlot.getState('isoScale')
        if levels[0] <= 2: levels[0] = 3

    # Nombre de separation sur la barre
    Ns = min(levels[0], 10)
    cbar_ticks = numpy.linspace(levels[1], levels[2], levels[0])
    nlevels = numpy.linspace(levels[1], levels[2], levels[0])

    if discrete:
        norm = matplotlib.colors.BoundaryNorm(nlevels, cmap.N)
        cset = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    else:
        norm = matplotlib.colors.Normalize(levels[1], levels[2], clip=True)
        cset = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

    fontdict = {'fontsize':fontSize}
    from mpl_toolkits.axes_grid1 import make_axes_locatable    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(location, size="2%", pad=pad)
    if location == "top" or location == "bottom": orientation = "horizontal"
    else: orientation = "vertical"
    cbar = fig.colorbar(cset, format=valueFormat, ticks=cbar_ticks, cax=cax, orientation=orientation)

    # Regle les marqueurs de la barre
    cbar.ax.tick_params(which='major', length=4., width=0.5, color=color, labelsize=fontSize)
    cbar.ax.tick_params(labelcolor=color)

    # Titre de la barre
    if title is not None:
        cbar.ax.set_title(title, fontdict=fontdict, color=color, pad=0.15*dpi)

    fig.subplots_adjust(wspace=0.)    
    return cbar

#==========================================================
# Return a numpy image from a file
def getImage(fileName):
    """Read image from file."""
    img = plt.imread(fileName)
    return img

#==========================================================
# Create subplot with image in it
# IN: img: image or image file name
# OUT: fig, ax: figure and axes 
#==========================================================
def createSubPlot(img='.decorator.png', title=None, box=False):
    """Create a sub plot figure."""
    if isinstance(img, str): img = getImage(img)
    sh = img.shape
    win = (sh[1], sh[0])
    fig, ax = plt.subplots(figsize=(win[0]/float(dpi), win[1]/float(dpi)), dpi=dpi)
    ax.imshow(img, animated=True, aspect=1)
    # Axes vides
    #ax.plot([], [], 'o', ms=8., mfc='None')
    if box: ax.set(xticks=[], yticks=[])
    else: ax.set_axis_off()
    if title is not None: ax.set_title(title)
    return fig, ax

#==========================================================
# Create a text 
#==========================================================
def createText(ax, posx=0, posy=0, text='', size=20, color="black",
               box=False, boxColor="black", boxBackColor="white", **kwargs):
    """Create text."""
    if not box:
        ax.text(posx, posy, text, size=size, ha='left', va='bottom', color=color, transform=ax.transAxes, **kwargs)
    else:
        ax.text(posx, posy, text, size=size, ha='left', va='bottom', color=color, transform=ax.transAxes,
                bbox=dict(boxstyle="round, pad=0.2, rounding_size=0.02", ec=boxColor, fc=boxBackColor), **kwargs)
    return ax

#==========================================================
# Save current figure to fileName
#==========================================================
def savefig(fileName, pad=0.):
    """Save current figure in a file."""
    print("Write %s"%fileName)
    plt.savefig(fileName, dpi=dpi, bbox_inches='tight', pad_inches=pad)

#==========================================================
# Show current fig
#==========================================================
def show():
    """Show current figure."""
    plt.show()

#==============================================================
# Change xyz (3D) to image position (written by Luis Bernardos)
#==============================================================
def xyz2Pixel__(points, win, posCam, posEye, dirCam, viewAngle):
    """Return the two-component image-pixel positions of a set of points located in the 3D world of CPlot."""

    # ------------------------------- #
    # BUILD FRENET UNIT FRAME (b,n,c) #
    # ------------------------------- #
    # <c> is the Camera axes unit vector
    c = numpy.array(posCam) - numpy.array(posEye)
    R = numpy.sqrt(c.dot(c)) # R is distance between posCam and posEye
    c /= R

    # <b> is binormal
    b = numpy.cross(numpy.array(dirCam),c)
    b /= numpy.sqrt(b.dot(b))

    # <n> is normal
    n = numpy.cross(c,b)
    n /= numpy.sqrt(b.dot(b))

    # <h> is the constant total height of the curvilinear window
    va = numpy.deg2rad(viewAngle)
    h = R * va
    h = 2 * R * numpy.tan(va/2.)

    # used to transform curvilinear unit to pixel
    crv2Pixel = float(win[1]) / h

    pixels = []

    # The window plane is defined as a set of three points (p0, p1, p2)
    p0 = numpy.array(posEye)
    p1 = p0+b
    p2 = p0+n
    p01 = p1 - p0 # segment
    p02 = p2 - p0 # segment

    for point in points:
        # ----------------------------------- #
        # COMPUTE pixel-position of point <p> #
        # ----------------------------------- #
        p = numpy.array(point)

        # Shall compute the intersection of the view of point <p> with the window plane

        # Such line is defined through two points (la, lb) as
        la, lb = numpy.array(posCam), p
        lab = lb - la # segment

        # Intersection is given by equation x = la + lab*t
        den = -lab.dot(numpy.cross(p01,p02))

        # Only for information (computation not required):
        # t = numpy.cross(p01,p02).dot(la-p0) / den
        # x = la + lab*t

        # parametric components (u, v) are actually what we look for
        u = numpy.cross(p02,-lab).dot(la-p0) / den
        v = numpy.cross(-lab,p01).dot(la-p0) / den

        # width and height in pixels are expressed in terms of (u, v)
        # Pixel components relative to Figure origin (upper left)
        pxP_w =  u*crv2Pixel + 0.5*float(win[0])
        pxP_h = -v*crv2Pixel + 0.5*float(win[1])

        pixels += [[pxP_w, pxP_h]]
    return pixels

def xyz2Pixel(Xs):
    """Transform 3D coordinates in pixel image coordinates for a set of points."""
    posCam = CPlot.getState("posCam")
    posEye = CPlot.getState("posEye")
    dirCam = CPlot.getState("dirCam")
    viewAngle = CPlot.getState("viewAngle")
    win = CPlot.getState("win")
    return xyz2Pixel__(Xs, win, posCam, posEye, dirCam, viewAngle)

#==========================================================
# Draw an arrow. X2 is the arrow head.
#==========================================================
def createArrow(ax, X1, X2, width=0.001, text=None, textSize=10, shiftText=(0,0), **kwargs):
    """Draw an arrow on figure."""
    poss = xyz2Pixel([X1,X2])
    pos1x, pos1y = poss[0]; pos2x, pos2y = poss[1]
    ax.arrow(pos1x, pos1y, pos2x-pos1x, pos2y-pos1y, width=width, length_includes_head=True, **kwargs)
    if text is not None:
        sx = shiftText[0]
        sy = shiftText[1]
        if 'color' in kwargs:
            color = kwargs['color']
        else: color = 'black'
        ax.text(pos1x+sx, pos1y+sy, s=text, size=textSize, color=color)
    return ax