# matplotlib decorators for CPlot
# strongly inspired by Luis Bernardos ideas

from . import PyTree as CPlot
from . import ColorMaps
import matplotlib.pyplot as plt
import matplotlib
import numpy

#==========================================================
def setBatch(batch=True):
    """Set batch mode for matplotlib."""
    if batch: matplotlib.use('Agg') # sans serveur X
    else: matplotlib.use('TkAgg') # avec serveur X

#==========================================================
def getInfo2DMode(xlim, ylim, zplane, ppw):
    """Toggle 2D mode for CPlot display."""
    xmin, xmax = xlim
    ymin, ymax = ylim
    cx = (xmax+xmin)/2.
    cy = (ymax+ymin)/2.
    width = xmax-xmin
    height = ymax-ymin
    pph = int((ppw*height)/width) # pixels per height

    D = max(width, height) # focal length = distance between posCam and posEye
    viewAngle = 2*numpy.arctan(height/(2*D)) # camera's vertical angle of view based on focal length
    viewAngle *= 180./numpy.pi

    posCam = (cx, cy, zplane+D)
    posEye = (cx, cy, zplane)
    dirCam = (0., 1., 0.)
    exportResolution="%dx%d"%(ppw, pph)

    return posCam, posEye, dirCam, viewAngle, exportResolution

#==========================================================
def createColormap(ctype='Blue2Red', colormapC1=None, colormapC2=None, colormapC3=None, colormapC=None, colors=None):
    """Create colormap."""

    if ctype in ['Jet', 'Magma', 'Inferno', 'Viridis', 'Plasma', 'Diverging', 'Greys', 'Greens']:
        colors = ColorMaps.export2MatplotLib2(ColorMaps.cmapDict[ctype])
    elif ctype == 'NiceBlue':
        c1 = matplotlib.colors.to_rgba('#000000')
        c2 = matplotlib.colors.to_rgba('#FFFFFF')
        c3 = matplotlib.colors.to_rgba('#0061A5')
        colors = [(0.00,[c1[0],c1[1],c1[2]]),
                  (0.50,[c3[0],c3[1],c3[2]]),
                  (1.00,[c2[0],c2[1],c2[2]])]
    elif ctype == 'BiColor':
        c1 = colormapC1
        c2 = colormapC2
        colors = [(0.00,[c1[0],c1[1],c1[2]]),
                  (1.00,[c2[0],c2[1],c2[2]])]
    elif ctype == 'TriColor':
        c1 = colormapC1
        c2 = colormapC2
        c3 = colormapC3
        colors = [(0.00,[c1[0],c1[1],c1[2]]),
                  (0.50,[c3[0],c3[1],c3[2]]),
                  (1.00,[c2[0],c2[1],c2[2]])]
    elif ctype == 'MultiColor':
        colors = colormapC
    else: # Blue2Red par defaut
        ctype = 'Blue2Red'
        colors = [(0.00,[0,0,1]),
                  (0.25,[0,1,1]),
                  (0.50,[0,1,0]),
                  (0.75,[1,1,0]),
                  (1.00,[1,0,0])]

    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(ctype, colors)
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
def createColorBar(fig, ax, levels=None, title=None, cmap=None, discrete=False, location='right', size='3%', pad=0.5, **kwargs):
    """Create a color bar."""

    ctype = 'Blue2Red'
    colormapC = None
    colormapC1 = kwargs.get('colormapC1', '#000000')
    if isinstance(colormapC1, str):
        colormapC1 = matplotlib.colors.to_rgba(colormapC1)
    colormapC2 = kwargs.get('colormapC2', '#FFFFFF')
    if isinstance(colormapC2, str):
        colormapC2 = matplotlib.colors.to_rgba(colormapC2)
    colormapC3 = kwargs.get('colormapC3', '#777777')
    if isinstance(colormapC3, str):
        colormapC3 = matplotlib.colors.to_rgba(colormapC3)

    if cmap is None:
        cmap = CPlot.getState('colormap')
        if cmap in [2,3,4,5]:
            ctype = 'BiColor' # BiColorRGB & BiColorHSV
            colormapC1 = CPlot.getState('colormapC1')
            colormapC2 = CPlot.getState('colormapC2')
        elif cmap in [6,7,8,9]:
            ctype = 'TriColor' # TriColorRGB & TriColorHSV
            colormapC1 = CPlot.getState('colormapC1')
            colormapC2 = CPlot.getState('colormapC2')
            colormapC3 = CPlot.getState('colormapC3')
        elif cmap in [10,11,12,13]:
            ctype = 'MultiColor' # MultiColorRGB & MultiColorHSV
            colormapC = CPlot.getState('colormapC')
            colormapC = ColorMaps.export2MatplotLib2(colormapC)
        elif cmap in [14,15]: ctype = 'Diverging'
        else: ctype = 'Blue2Red'
    elif isinstance(cmap, str):
        ctype = cmap
    elif isinstance(cmap, list):
        ctype = 'MultiColor'
        colormapC = cmap

    cmap = createColormap(ctype, colormapC1=colormapC1, colormapC2=colormapC2, colormapC3=colormapC3, colormapC=colormapC)

    if levels is None:
        levels = CPlot.getState('isoScale')
        if levels[0] <= 2: levels[0] = 3

    # kwargs arguments
    extend = kwargs.get('extend', 'neither') # 'neither', 'both', 'min', 'max'
    extendrect = kwargs.get('extendrect', True)
    fontsize = kwargs.get('fontsize', 12)
    labelsize = kwargs.get('labelsize', 10)
    titlepad = kwargs.get('titlepad', 12)
    labelFormat = kwargs.get('labelFormat', '%.2f')
    nticks = kwargs.get('nticks', 5)

    if discrete: # adapt levels[0] if necessary
        while (levels[0]-1)%(nticks-1) != 0: levels[0] += 1

    # Nombre de separation sur la barre
    nlevels = numpy.linspace(levels[1], levels[2], levels[0])
    cticks = numpy.linspace(levels[1], levels[2], nticks)

    if discrete:
        norm = matplotlib.colors.BoundaryNorm(nlevels, cmap.N)
        cset = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    else:
        norm = matplotlib.colors.Normalize(levels[1], levels[2], clip=True)
        cset = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(location, size=size, pad=pad)
    if location == "top" or location == "bottom": orientation = "horizontal"
    else: orientation = "vertical"
    cbar = fig.colorbar(cset, format=labelFormat, ticks=cticks, cax=cax, orientation=orientation, extendrect=extendrect, extend=extend)

    # Regle les marqueurs de la barre
    cbar.ax.tick_params(which='major', labelsize=labelsize)

    # Titre de la barre
    if title is not None: cbar.ax.set_title(title, fontsize=fontsize, pad=titlepad)

    fig.subplots_adjust(wspace=0.)
    return cbar

#==========================================================
# Return a numpy image from a file
#==========================================================
def getImage(fileName):
    """Read image from file."""
    img = plt.imread(fileName)
    return img

#==========================================================
# Create subplot with image in it
# IN: img: image or image file name
# OUT: fig, ax: figure and axes
#==========================================================
def createSubPlot(img='.decorator.png', figsize=None, dpi=None, title=None, box=False, xlim=[], ylim=[], **kwargs):
    """Create a subplot figure."""
    if isinstance(img, str): img = getImage(img)
    sh = img.shape
    winx, winy = sh[1], sh[0]

    # kwargs arguments
    fontsize = kwargs.get('fontsize', 12.)
    labelsize = kwargs.get('labelsize', 10.)
    interpolation = kwargs.get('interpolation', 'none')

    extent = None
    if xlim != [] and ylim != []: extent = xlim+ylim

    if figsize is None:
        figsize = (8, 8*winy/winx) # automatic figsize setting based on image ratio

    if dpi is None:
        dpi = winx/figsize[0] # automatic dpi setting based on the image actual size

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.imshow(img, animated=True, aspect=1., extent=extent, interpolation=interpolation)
    if extent is not None: ax.set(xlim=xlim, ylim=ylim)

    if box:
        if extent is not None: # xlim and ylim are defined
            xlabel = kwargs.get('xlabel', 'x')
            xlabelpad = kwargs.get('xlabelpad', 11.)
            ax.set_xlabel(xlabel, fontsize=fontsize, labelpad=xlabelpad)

            ylabel = kwargs.get('ylabel', 'y')
            ylabelpad = kwargs.get('ylabelpad', 11.)
            ax.set_ylabel(ylabel, fontsize=fontsize, labelpad=ylabelpad)

            ax.tick_params(which='major', labelsize=labelsize)
        else:
            boxColor = kwargs.get('boxColor', 'black')
            boxLineWidth = kwargs.get('boxLineWidth', 1.)
            ax.set(xticks=[], yticks=[])
            for spine in ax.spines.values():
                spine.set_linewidth(boxLineWidth)
                spine.set_color(boxColor)
    else:
        ax.set_axis_off()

    if title is not None:
        titlepad = kwargs.get('titlepad', 12.)
        ax.set_title(title, fontsize=fontsize, pad=titlepad)

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
def savefig(fileName, pad=0., dpi='figure', **kwargs):
    """Save current figure in a file."""
    print("Write %s"%fileName)
    plt.savefig(fileName, dpi=dpi, bbox_inches='tight', pad_inches=pad, **kwargs)

#==========================================================
# Show current fig
#==========================================================
def show():
    """Show current figure."""
    plt.show()

#==========================================================
# Close all figures
#==========================================================
def closeAll():
    """Close all figures."""
    plt.close('all')

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