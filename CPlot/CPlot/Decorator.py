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
    """Set batch mode."""
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
def createColorBar(fig, ax, levels, title=None, cmap=None, discrete=True):
    if cmap is None: 
        cmap = CPlot.getState('colormap')
        if cmap == 0 or cmap == 1:
            cmap = createColormap('Blue2Red')
        else: # default
            cmap = createColormap('Blue2Red')
    elif isinstance(cmap, str): cmap = createColormap(cmap)

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

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.5)
    cbar = fig.colorbar(cset, format='%0.2f', ticks=cbar_ticks, cax=cax)
    # Regle les marqueurs de la barre
    cbar.ax.tick_params(which='major', length=4., width=0.5, color='black')
    # Titre de la barre
    if title is not None:
        fontdict = {'fontsize':20}
        cbar.ax.set_title(title, fontdict=fontdict, color='black', pad=0.15*dpi)
    
    fig.subplots_adjust(right=1.)
    # plt.tight_layout()

#==========================================================
# Return a numpy image size
def getImage(fileName):
    img = plt.imread(fileName)
    return img

#==========================================================
# Create subplot with image in it
# IN: img: image or image file name
# OUT: fig, ax: figure and axes 
#==========================================================
def createSubPlot(img='.decorator.png', title=None):
    if isinstance(img, str): img = getImage(img)
    sh = img.shape
    win = (sh[1], sh[0])
    fig, ax = plt.subplots(figsize=(win[0]/float(dpi), win[1]/float(dpi)), dpi=dpi)
    ax.imshow(img, animated=True, aspect=1)
    # Axes vides
    ax.plot([], [], 'o', ms=8., mfc='None')
    # ax.set_axis_off()
    ax.set(xticks=[], yticks=[])
    if title is not None: ax.set_title(title)
    return fig, ax

#==========================================================
# Save current figure to fileName
#==========================================================
def saveFig(fileName, pad=0.2):
    print("Write %s."%fileName)
    plt.savefig(fileName, dpi=dpi, bbox_inches='tight', pad_inches=pad)
