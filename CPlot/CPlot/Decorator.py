# matplotlib decorators for CPlot
import CPlot.PyTree as CPlot
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
    elif type == 'Black2White':
        colors = [(0.00,[0.,0.,0.]),
                  (1.00,[1.,1.,1.])]
    elif type == 'Viridis':
        colors = [(0.00,[253./255.,231./255.,37./255.]),
                  (0.50,[33./255.,145./255.,140./255.]),
                  (1.00,[68./255.,1./255.,84./255.])]
    elif type == 'Inferno':
        colors = [(0.00,[252./255.,255./255.,164./255.]),
                  (0.50,[188./255.,55./255.,84./255.]),
                  (1.00,[0./255.,0./255.,4./255.])]
    elif type == 'Magma':
        colors = [(0.00,[252./255.,253./255.,191./255.]),
                  (0.50,[183./255.,55./255.,121./255.]),
                  (1.00,[0./255.,0./255.,4./255.])]
    elif type == 'Plasma':
        colors = [(0.00,[240./255.,249./255.,33./255.]),
                  (0.50,[204./255.,71./255.,120./255.]),
                  (1.00,[13./255.,8./255.,135./255.])]
    elif type == 'NiceBlue':
        colors = [(0.00,[0./255.,0./255.,0./255.]),
                  (0.50,[255./255.,255./255.,255./255.]),
                  (1.00,[0./255.,97./255.,165./255.])]
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
def createColorBar(fig, ax, levels=None, title=None, cmap=None, valueFormat='%0.3f', discrete=True):
    """Create a color bar."""
    if cmap is None: 
        cmap = CPlot.getState('colormap')
        if cmap == 0 or cmap == 1:
            cmap = createColormap('Blue2Red')
        elif cmap == 2 or cmap == 3:
            cmap = createColormap('Green2Red')
        elif cmap == 4 or cmap == 5:
            cmap = createColormap('BiColorRGB')
        elif cmap == 6 or cmap == 7:
            cmap = createColormap('BiColorHSV')
        elif cmap == 8 or cmap == 9:
            cmap = createColormap('Diverging')
        elif cmap == 10 or cmap == 11:
            cmap = createColormap('TriColorRGB')
        elif cmap == 12 or cmap == 13:
            cmap = createColormap('TriColorHSV')
        elif cmap == 14 or cmap == 15: 
            cmap = createColormap('Black2White')
        elif cmap == 16 or cmap == 17: 
            cmap = createColormap('Veridis')
        elif cmap == 18 or cmap == 19: 
            cmap = createColormap('Inferno')
        elif cmap == 20 or cmap == 21: 
            cmap = createColormap('Magma')
        elif cmap == 22 or cmap == 23: 
            cmap = createColormap('Plasma')
        elif cmap == 24 or cmap == 25: 
            cmap = createColormap('NiceBlue')
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

    fontdict = {'fontsize':25}
    from mpl_toolkits.axes_grid1 import make_axes_locatable    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.5)
    cbar = fig.colorbar(cset, format=valueFormat, ticks=cbar_ticks, cax=cax)

    # Regle les marqueurs de la barre
    cbar.ax.tick_params(which='major', length=4., width=0.5, color='black', labelsize=20)

    # Titre de la barre
    if title is not None:
        cbar.ax.set_title(title, fontdict=fontdict, color='black', pad=0.15*dpi)
    
    fig.subplots_adjust(right=1.)
    return cbar

#==========================================================
# Return a numpy image size
def getImage(fileName):
    """Read image from file."""
    img = plt.imread(fileName)
    return img

#==========================================================
# Create subplot with image in it
# IN: img: image or image file name
# OUT: fig, ax: figure and axes 
#==========================================================
def createSubPlot(img='.decorator.png', title=None):
    """Create a sub plot figure."""
    if isinstance(img, str): img = getImage(img)
    sh = img.shape
    win = (sh[1], sh[0])
    fig, ax = plt.subplots(figsize=(win[0]/float(dpi), win[1]/float(dpi)), dpi=dpi)
    ax.imshow(img, animated=True, aspect=1)
    # Axes vides
    #ax.plot([], [], 'o', ms=8., mfc='None')
    # ax.set_axis_off()
    ax.set(xticks=[], yticks=[])
    if title is not None: ax.set_title(title)
    return fig, ax

#==========================================================
# Create a text 
#==========================================================
def createText(ax, posx=0, posy=0, text='', size=20, box=False):
    """Create text."""
    if not box:
        ax.text(posx, posy, text, size=size, ha='left', va='bottom', transform=ax.transAxes)
    else:
        ax.text(posx, posy, text, size=size, ha='left', va='bottom', transform=ax.transAxes,
                bbox=dict(boxstyle="round, pad=0.2, rounding_size=0.02", ec="black", fc="white"))
    return ax

#==========================================================
# Save current figure to fileName
#==========================================================
def savefig(fileName, pad=0.2):
    """Save current figure."""
    print("Write %s"%fileName)
    plt.savefig(fileName, dpi=dpi, bbox_inches='tight', pad_inches=pad)

#==========================================================
# Show current fig
#==========================================================
def show():
    """Show current figure."""
    plt.show()
