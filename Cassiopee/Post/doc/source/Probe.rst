.. Probe documentation master file

:tocdepth: 2


Post.Probe: Probe extraction
===========================================================

Specific probe extraction module.

This page describes how the object-oriented Probe module operates. The various methods presented below should be used after instantiating a Probe object of the same-named class. 

All of these methods and functions can be executed in both sequential and parallel contexts.

.. py:module:: Post.Probe

List of functions
#################


**-- Probe creation**

.. autosummary::
    :nosignatures:

    Post.Probe.Probe

**-- Probe methods**

.. autosummary::
    :nosignatures:

    Probe.printInfo
    Probe.prepare
    Probe.extract
    Probe.flush
    Probe.read

Contents
########

Probe creation
---------------------------------------------------------------------------

.. py:function:: Post.Probe.Probe(fileName, t=None, X=None, ind=None, blockName=None, tPermeable=None, fields=None, append=False, bufferSize=100)

    Create a Python object from the Probe class. Five modes are possible:
    
    + **mode 0**: if t and X are provided, the probe will extract the fields from t at the absolute position X=(x,y,z).
    + **mode 1**: if t and (ind, blockName) are provided, the probe will extract the fields from t at the given index of the specified block.
    + **mode 2**: if t is provided without a specific location (index or position), the probe will store all zones of t with both coordinates and fields. t must be 1D or 2D.
    + **mode 3**: if t and tPermeable are provided, the probe will interpolate the fields on tPermeable using t as the donor tree. tPermeable must be 1D or 2D.
    + **mode 4**: if t is not provided, the probe will simply store specified values over time.

    Data are periodically flushed to probe file when the buffer size exceeds the user input bufferSize.

    For Probe objects in modes 0, 1, or 4, all data are stored in 1D numpy arrays of size bufferSize.

    For Probe objects in modes 2 or 3, the data are stored in 2D/3D numpy arrays of size (bufferSize,ni,nj) where (ni,nj) are the original dimensions of the receiver tree zones (t for mode 2 or tPermeable for mode 3).

    .. note:: Modes 2 and 3 only operates with structured meshes.

    :param fileName: name of the probe file
    :type fileName: string
    :param t: pyTree containing the flow solution
    :type t: pyTree
    :param X: absolute position of a single probe (mode 0 only) 
    :type  X: tuple of 3 floats (x,y,z)
    :param ind: index of a single probe located in blockName (mode 1 only)
    :type ind: integer or tuple of 3 integers (i,j,k)
    :param blockName: name of the block containing the probe index (mode 1 only)
    :type blockName: string
    :param tPermeable: interpolated tree (mode 3 only)
    :type tPermeable: pyTree, zone or list of zones
    :param fields: list of fields to extract (located at the mesh nodes OR the mesh centers)
    :type fields: list of strings 
    :param append: if True, append result to existing file
    :type append: Boolean
    :param bufferSize: size of internal buffer
    :type bufferSize: int
    
    :rtype: probe instance

    *Example of use:*

    * `Probe - mode 0 (pyTree) <Examples/Post/probeMode0PT.py>`_:

    .. literalinclude:: ../build/Examples/Post/probeMode0PT.py

    * `Probe - mode 1 (pyTree) <Examples/Post/probeMode1PT.py>`_:

    .. literalinclude:: ../build/Examples/Post/probeMode1PT.py

    * `Probe - mode 2 (pyTree) <Examples/Post/probeMode2PT.py>`_:

    .. literalinclude:: ../build/Examples/Post/probeMode2PT.py

    * `Probe - mode 3 (pyTree) <Examples/Post/probeMode3PT.py>`_:

    .. literalinclude:: ../build/Examples/Post/probeMode3PT.py

    * `Probe - mode 4 (pyTree) <Examples/Post/probeMode4PT.py>`_:

    .. literalinclude:: ../build/Examples/Post/probeMode4PT.py

---------------------------------------

Probe methods
---------------------------------------------------------------------------

.. py:function:: Probe.printInfo()

    **printInfo** is a method of the Probe class.

    Print all the information about a probe object.

    *Example of use:*

    * `Print information about a given probe (pyTree) <Examples/Post/probePrintInfoPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/probePrintInfoPT.py

---------------------------------------

.. py:function:: Probe.prepare(t, loc='nodes', extrap=1, nature=1, penalty=1, verbose=2)

    **prepare** is a method of the Probe class.

    Only for mode 3. Prepare the interpolation data for probe extraction from a donor computational tree to a 1D or 2D receiver mesh.

    The arguments are mostly the same as those for Connector.setInterpData2(). Donor tree flow fields must be located at the nodes. Receiver tree flow fields can be located at the nodes (loc='nodes') or at the centers (loc='centers').

    If the cellN field (cell nature field) is not found in the donor tree, the function assumes that all cells are computed cells (cellN = 1). Similarly, if the cellN field is not found in the receiver tree, the function assumes that all cells are interpolated cells (cellN = 2).

    extrap options:

    + **0**: extrapolation is disabled. Extrapolated cells are flagged as orphans.
    + **1**: extrapolation is enabled.

    nature options:

    + **0**: candidate donors can be everything but blanked cells (cellN = 0).
    + **1**: candidate donors can only be computed cells (cellN = 1).

    penalty options: 

    + **0**: all candidates have the same weight.
    + **1**: donor cell candidates located at a zone border are penalized against interior cells.

    verbose options:

    + **0**: no information is printed.
    + **1**: print only the summary of the interpolation data search.
    + **2**: print the summary and indices of the orphan points (if any).

    :param t: pyTree containing the flow solution
    :type t: pyTree
    :param loc: extraction location
    :type loc: string ('centers' or 'nodes')
    :param extrap: extrapolation parameter
    :type extrap: integer (0 or 1)
    :param nature: donor cell nature parameter 
    :type nature: integer (0 or 1)
    :param penalty: penalization parameter
    :type penalty: integer (0 or 1)
    :param verbose: print parameter
    :type verbose: integer (0, 1, or 2)

    *Example of use:*

    * `Probe (mode 3) prepare operation for receiver flow fields located at centers (pyTree) <Examples/Post/probePrepareCenterPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/probePrepareCenterPT.py

    * `Probe (mode 3) prepare operation for receiver flow fields located at nodes (pyTree) <Examples/Post/probePrepareNodePT.py>`_:

    .. literalinclude:: ../build/Examples/Post/probePrepareNodePT.py

---------------------------------------

.. py:function:: Probe.extract(t=None, time=-1., value=[], onlyTransfer=False)

    **extract** is a method of the Probe class.

    Extract the probe information (coordinates and field) at a given time from t.

    Flow solution and grid coordinate extractions from the donor tree t are based on the information located in the Internal.__FlowSolutionNodes__, Internal.__FlowSolutionCenters, and Internal.__FlowSolutionCenters__ containers. Users can modify these container names by editing the constant names in the Internal module. New containers are created each time the local buffer size exceeds the bufferSize limits specified by the user.

    Probe information data is then stored in the following fixed container names: GridCoordinates#i, FlowSolution#i, and FlowSolution#Centers#i, where i is an iterator that is automatically set and handled by the Probe class.

    To extract the receiver tree interpolated data at a given time without stacking the information (coordinates and field) in the probe, set onlyTransfer to True for mode 3. This option is useful when users only want to monitor integral quantities instead of monitoring all solution fields over time.

    :param t: pyTree containing the flow solution
    :type t: pyTree
    :param time: extraction time
    :type time: float
    :param value: list of values to be stored in probe (only for mode 4)
    :type value: list of floats
    :param onlyTransfer: deactivate information stacking in the probe object (only for mode 3)
    :type onlyTransfer: boolean

    *Example of use:*

    .. note:: for simpler examples of each mode, see the Post.Probe.Probe function definition.

    * `Probe extraction using both mode 3 and mode 4 to monitor integral quantities (pyTree) <Examples/Post/probeExtractPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/probeExtractPT.py
    
---------------------------------------

.. py:function:: Probe.flush()

    **flush** is a method of the Probe class.

    Force buffered data to be written in the associated probe file before reaching the user input bufferSize limit.

    *Example of use:*

    * `Flush unsaved probe data to disk. (pyTree) <Examples/Post/probeFlushPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/probeFlushPT.py

---------------------------------------

.. py:function:: Probe.read(cont=None, ind=None, probeName=None)

    Read the data stored in the probe file and return either a zone or a list of zones.

    The Probe object must already be instantiated with the desired probe file name. Otherwise, use Post.Probe.Probe(fileName) to create a new probe object from an existing probe file for reading purposes only.

    Can be used in two ways:

    + **extract all points**: if cont is provided, this function extracts all data from the given time container. 
    + **extract all times**: if ind and probeName are provided, this function now extracts the given index value(s) of the given zone at all times.

    If ind is provided but probeName is not, the function will automatically selects the first zone in the probe.

    When extracting all points, the output is a list of zones. Each zone corresponds to a probe zone at a single time included in the loaded container.

    When extracting all times, the output is a single zone. This zone contains numpy arrays of size (ntime,npoints) depending on the number of point indices provided.

    .. note:: For modes 0, 1, and 4, probes consist of only one zone named 'probe'. For modes 2 and 3, probes contain the same number of zones as the original tree (t for mode 2 or tPermeable for mode 3).

    :param cont: time container number
    :type cont: integer
    :param ind: index (or indices) of probe point(s) located in probeName
    :type ind: integer, tuple of 3 integers (i,j,k), list of integers, or list of tuples of 3 integers [(i1,j1,k1), (i2,j2,k2), ...]
    :param probeName: only used when ind is provided. Name or number of the probe zone to extract
    :type probeName: string or int

    :rtype: zone or list of zones

    *Example of use:*

    * `Read data from probe file. (pyTree) <Examples/Post/probeReadPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/probeReadPT.py

---------------------------------------
