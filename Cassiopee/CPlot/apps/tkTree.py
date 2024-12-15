# - tkTree -
"""View a pyTree in a tree widget."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import tkNodeEdit # for node inspector
import Converter.Internal as Internal
try: import tkinter.dnd as Tkdnd # drag and drop
except: import Tkdnd
import numpy

try: range = xrange
except: pass

# local widgets list
STATUS = 0
WIDGETS = {}; VARS = []

# buffer for copy/paste
BUFFER = None

def strFormat(v):
    if isinstance(v, int): return "%ld"%v
    elif isinstance(v, float): return "%g"%v
    else: return "%s"%v

#==============================================================================
def report_callback_exception():
    """report exception on sys.stderr."""
    import traceback
    import sys
    sys.stderr.write("Exception in Tree control callback.\n")
    traceback.print_exc()

#------------------------------------------------------------------------------
class Struct:
    """Helper object for add_node() method"""
    def __init__(self):
        pass

#------------------------------------------------------------------------------
class Node:
    """Tree helper class that's instantiated for each element in the tree.  It
    has several useful attributes:
    parent_node     - immediate parent node
    id              - id assigned at creation
    expanded_icon   - image displayed when folder is expanded to display
                      children
    collapsed_icon  - image displayed when node is not a folder or folder is
                      collapsed.
    parent_widget   - reference to tree widget that contains node.
    expandable_flag - is true when node is a folder that may be expanded or
                      collapsed.
    expanded_flag   - true to indicate node is currently expanded.
    h_line          - canvas line to left of node image.
    v_line          - canvas line below node image that connects children.
    indic           - expand/collapse canvas image.
    label           - canvas text label
    symbol          - current canvas image

    Please note that methods prefixed PVT_* are not meant to be used by
    client programs."""

    def __init__(self, parent_node, id, collapsed_icon, x, y,
                 parent_widget=None, expanded_icon=None, label=None,
                 expandable_flag=0):
        """Create node and initialize it.  This also displays the node at the
        given position on the canvas, and binds mouseclicks."""
        # immediate parent node
        self.parent_node = parent_node
        # internal name used to manipulate things
        self.id = id
        # bitmaps to be displayed
        self.expanded_icon = expanded_icon
        self.collapsed_icon = collapsed_icon
        # tree widget we belong to
        if parent_widget: self.widget = parent_widget
        else: self.widget = parent_node.widget
        # for speed
        sw = self.widget
        # our list of child nodes
        self.child_nodes = []
        # flag that node can be expanded
        self.expandable_flag = expandable_flag
        self.expanded_flag = 0
        # add line
        if parent_node and sw.line_flag:
            self.h_line = sw.create_line(x, y, x-sw.dist_x, y)
        else:
            self.h_line = None
        self.v_line = None
        # draw appropriate image
        self.symbol = sw.create_image(x, y, image=self.collapsed_icon)
        # add expand/collapse indicator
        self.indic = None
        if expandable_flag and sw.line_flag and sw.plus_icon and sw.minus_icon:
            self.indic=sw.create_image(x-sw.dist_x, y, image=sw.plus_icon)
        # add label
        self.label=sw.create_text(x+sw.text_offset, y, text=label,
                                  anchor='w', font=CTK.TEXTFONT)
        # single-click to expand/collapse
        if self.indic:
            sw.tag_bind(self.indic, '<Button-1>', self.PVT_click)
            sw.tag_bind(self.symbol, '<Button-1>', self.PVT_click)
            sw.tag_bind(self.label, '<Button-1>', self.PVT_clickSelect)
            sw.tag_bind(self.label, '<Control-Shift-Button-1>',
                        self.PVT_clickMultipleSelect)
            sw.tag_bind(self.label, '<Double-Button-1>',
                        self.PVT_clickEdit)
            sw.tag_bind(self.label, '<Button-3>', self.PVT_clickRight)
        else:
            sw.tag_bind(self.symbol, '<Button-1>', self.PVT_clickSelect)
            sw.tag_bind(self.label, '<Button-1>', self.PVT_clickSelect)
            sw.tag_bind(self.label, '<Double-Button-1>', self.PVT_clickEdit)
            sw.tag_bind(self.label, '<Button-3>', self.PVT_clickRight)
        # For editing
        sw.tag_bind(self.label, '<Key>', self.PVT_key)
        # for drag'n'drop target detection
        sw.tag_bind(self.symbol, '<Any-Enter>', self.PVT_enter)
        sw.tag_bind(self.label, '<Any-Enter>', self.PVT_enter)

    # ----- PUBLIC METHODS -----
    def set_collapsed_icon(self, icon):
        """Set node's collapsed image"""
        self.collapsed_icon = icon
        if not self.expanded_flag:
            self.widget.itemconfig(self.symbol, image=icon)

    def set_expanded_icon(self, icon):
        """Set node's expanded image"""
        self.expanded_icon = icon
        if self.expanded_flag:
            self.widget.itemconfig(self.symbol, image=icon)

    def parent(self):
        """Return node's parent node"""
        return self.parent_node

    def prev_sib(self):
        """Return node's previous sibling (the child immediately above it)"""
        i = self.parent_node.child_nodes.index(self)-1
        if i >= 0:
            return self.parent_node.child_nodes[i]
        else: return None

    def next_sib(self):
        """Return node's next sibling (the child immediately below it)"""
        i = self.parent_node.child_nodes.index(self)+1
        if i < len(self.parent_node.child_nodes):
            return self.parent_node.child_nodes[i]
        else:
            return None

    def next_visible(self):
        """Return next lower visible node"""
        n = self
        if n.child_nodes:
            # if you can go right, do so
            return n.child_nodes[0]
        while n.parent_node:
            # move to next sibling
            i = n.parent_node.child_nodes.index(n)+1
            if i < len(n.parent_node.child_nodes):
                return n.parent_node.child_nodes[i]
            # if no siblings, move to parent's sibling
            n = n.parent_node
        # we're at bottom
        return self

    def prev_visible(self):
        """Return next higher visible node"""
        n = self
        if n.parent_node:
            i = n.parent_node.child_nodes.index(n)-1
            if i < 0:
                return n.parent_node
            else:
                j = n.parent_node.child_nodes[i]
                return j.PVT_last()
        else:
            return n

    def children(self):
        """Return list of node's children"""
        return self.child_nodes[:]

    def get_label(self):
        """Return string containing text of current label"""
        return self.widget.itemcget(self.label, 'text')

    def set_label(self, label):
        """Set current text label"""
        self.widget.itemconfig(self.label, text=label)

    def expanded(self):
        """Returns true if node is currently expanded, false otherwise"""
        return self.expanded_flag

    def expandable(self):
        """Returns true if node can be expanded (i.e. if it's a folder)"""
        return self.expandable_flag

    def full_id(self):
        """Return list of IDs of all parents and node ID"""
        if self.parent_node:
            return self.parent_node.full_id()+(self.id,)
        else: return (self.id,)

    def expand(self):
        """Expand node if possible"""
        if not self.expanded_flag: self.PVT_set_state(1)

    def collapse(self):
        """Collapse node if possible"""
        if self.expanded_flag:
            self.PVT_set_state(0)

    def delete(self, me_too=1):
        """Delete node from tree. ("me_too" is a hack not to be used by
        external code, please!)"""
        sw = self.widget
        if not self.parent_node and me_too:
            # can't delete the root node
            raise ValueError("can't delete root node.")
        self.PVT_delete_subtree()
        # move everything up so that distance to next subnode is correct
        n = self.next_visible()
        x1, y1 = sw.coords(self.symbol)
        x2, y2 = sw.coords(n.symbol)
        if me_too: dist = y2-y1
        else: dist = y2-y1-sw.dist_y
        self.PVT_tag_move(-dist)
        n = self
        if me_too:
            if sw.pos == self:
                # move cursor if it points to current node
                sw.move_cursor(self.parent_node)
            self.PVT_unbind_all()
            sw.delete(self.symbol)
            sw.delete(self.label)
            sw.delete(self.h_line)
            sw.delete(self.v_line)
            sw.delete(self.indic)
            self.parent_node.child_nodes.remove(self)
            # break circular ref now, so parent may be GC'ed later
            n = self.parent_node
            self.parent_node = None
        n.PVT_cleanup_lines()
        n.PVT_update_scrollregion()

    def insert_before(self, nodes):
        """Insert list of nodes as siblings before this node.  Call parent
        node's add_node() function to generate the list of nodes."""
        i = self.parent_node.child_nodes.index(self)
        self.parent_node.PVT_insert(nodes, i, self.prev_visible())

    def insert_after(self, nodes):
        """Insert list of nodes as siblings after this node.  Call parent
        node's add_node() function to generate the list of nodes."""
        i=self.parent_node.child_nodes.index(self)+1
        self.parent_node.PVT_insert(nodes, i, self.PVT_last())

    def insert_children(self, nodes):
        """Insert list of nodes as children of this node.  Call node's
        add_node() function to generate the list of nodes."""
        self.PVT_insert(nodes, 0, self)

    def toggle_state(self):
        """Toggle node's state between expanded and collapsed, if possible"""
        if self.expandable_flag:
            if self.expanded_flag: self.PVT_set_state(0)
            else: self.PVT_set_state(1)

    # ----- functions for drag'n'drop support -----
    def PVT_enter(self, event):
        """detect mouse hover for drag'n'drop"""
        self.widget.target = self

    def dnd_end(self, target, event):
        """Notification that dnd processing has been ended. It DOES NOT imply
        that we've been dropped somewhere useful, we could have just been
        dropped into deep space and nothing happened to any data structures,
        or it could have been just a plain mouse-click w/o any dragging."""
        if not self.widget.drag:
            # if there's been no dragging, it was just a mouse click
            self.widget.move_cursor(self)
            self.toggle_state()
        self.widget.drag = 0

    # ----- PRIVATE METHODS (prefixed with "PVT_") -----
    # these methods are subject to change, so please try not to use them
    def PVT_last(self):
        """Return bottom-most node in subtree"""
        n = self
        while n.child_nodes:
            n = n.child_nodes[-1]
        return n

    def PVT_find(self, search):
        """Used by searching functions"""
        if self.id != search[0]:
            # this actually only goes tilt if root doesn't match
            return None
        if len(search) == 1:
            return self
        # get list of children IDs
        i = map(lambda x: x.id, self.child_nodes)
        # if there is a child that matches, search it
        try:
            return self.child_nodes[i.index(search[1])].PVT_find(search[1:])
        except:
            return None

    def PVT_insert(self, nodes, pos, below):
        """Create and insert new children. "nodes" is list previously created
        via calls to add_list(). "pos" is index in the list of children where
        the new nodes are inserted. "below" is node which new children should
        appear immediately below."""
        if not self.expandable_flag: return
        #raise TypeError("not an expandable node")
        # for speed
        sw = self.widget
        # expand and insert children
        children = []
        self.expanded_flag = 1
        sw.itemconfig(self.symbol, image=self.expanded_icon)
        if sw.minus_icon and sw.line_flag:
            sw.itemconfig(self.indic, image=sw.minus_icon)
        if len(nodes):
            # move stuff to make room
            below.PVT_tag_move(sw.dist_y*len(nodes))
            # get position of first new child
            xp, dummy = sw.coords(self.symbol)
            dummy, yp = sw.coords(below.symbol)
            xp = xp+sw.dist_x
            yp = yp+sw.dist_y
            # create vertical line
            if sw.line_flag and not self.v_line:
                self.v_line=sw.create_line(
                    xp, yp,
                    xp, yp+sw.dist_y*len(nodes))
                sw.tag_lower(self.v_line, self.symbol)
            n = sw.node_class
            for i in nodes:
                # add new subnodes, they'll draw themselves
                # this is a very expensive call
                children.append(
                    n(parent_node=self, expandable_flag=i.flag, label=i.name,
                      id=i.id, collapsed_icon=i.collapsed_icon,
                      expanded_icon=i.expanded_icon, x=xp, y=yp))
                yp = yp+sw.dist_y
            self.child_nodes[pos:pos] = children
            self.PVT_cleanup_lines()
            self.PVT_update_scrollregion()
            sw.move_cursor(sw.pos)

    def PVT_set_state(self, state):
        """Common code forexpanding/collapsing folders. It's not re-entrant,
        and there are certain cases in which we can be called again before
        we're done, so we use a mutex."""
        while self.widget.spinlock:
            pass
        self.widget.spinlock=1
        # expand & draw our subtrees
        if state:
            self.child_nodes = []
            self.widget.new_nodes = []
            if self.widget.get_contents_callback:
                # this callback needs to make multiple calls to add_node()
                try:
                    self.widget.get_contents_callback(self)
                except:
                    report_callback_exception()
            self.PVT_insert(self.widget.new_nodes, 0, self)
        # collapse and delete subtrees
        else:
            self.expanded_flag=0
            self.widget.itemconfig(self.symbol, image=self.collapsed_icon)
            if self.indic:
                self.widget.itemconfig(self.indic, image=self.widget.plus_icon)
            self.delete(0)
        # release mutex
        self.widget.spinlock=0

    def PVT_cleanup_lines(self):
        """Resize connecting lines"""
        if self.widget.line_flag:
            n = self
            while n:
                if n.child_nodes:
                    x1, y1 = self.widget.coords(n.symbol)
                    x2, y2 = self.widget.coords(n.child_nodes[-1].symbol)
                    self.widget.coords(n.v_line, x1, y1, x1, y2)
                n = n.parent_node

    def PVT_update_scrollregion(self):
        """Update scroll region for new size"""
        x1, y1, x2, y2 = self.widget.bbox('all')
        self.widget.configure(scrollregion=(x1, y1, x2+5, y2+5))

    def PVT_delete_subtree(self):
        """Recursively delete subtree & clean up cyclic references to make
        garbage collection happy"""
        sw = self.widget
        sw.delete(self.v_line)
        self.v_line = None
        for i in self.child_nodes:
            # delete node's subtree, if any
            i.PVT_delete_subtree()
            i.PVT_unbind_all()
            # delete widgets from canvas
            sw.delete(i.symbol)
            sw.delete(i.label)
            sw.delete(i.h_line)
            sw.delete(i.v_line)
            sw.delete(i.indic)
            # break circular reference
            i.parent_node = None
        # move cursor if it's in deleted subtree
        if sw.pos in self.child_nodes:
            sw.move_cursor(self)
        # now subnodes will be properly garbage collected
        self.child_nodes = []

    def PVT_unbind_all(self):
        """Unbind callbacks so node gets garbage-collected. This wasn't easy
        to figure out the proper way to do this.  See also tag_bind() for the
        Tree widget itself."""
        for j in (self.symbol, self.label, self.indic, self.h_line,
                  self.v_line):
            for k in self.widget.bindings.get(j, ()):
                self.widget.tag_unbind(j, k[0], k[1])

    def PVT_tag_move(self, dist):
        """Move everything below current icon, to make room for subtree using
        the Disney magic of item tags.  This is the secret of making
        everything as fast as it is."""
        # mark everything below current node as movable
        bbox1=self.widget.bbox(self.widget.root.symbol, self.label)
        bbox2=self.widget.bbox('all')
        self.widget.dtag('move')
        self.widget.addtag('move', 'overlapping',
                           bbox2[0], bbox1[3], bbox2[2], bbox2[3])
        # untag cursor & node so they don't get moved too
        self.widget.dtag(self.widget.cursor_box, 'move')
        self.widget.dtag(self.symbol, 'move')
        self.widget.dtag(self.label, 'move')
        # now do the move of all the tagged objects
        self.widget.move('move', 0, dist)

    def PVT_click(self, event):
        """Handle mouse clicks by kicking off possible drag'n'drop
        processing"""
        if self.widget.drop_callback:
            if Tkdnd.dnd_start(self, event):
                x1, y1, x2, y2=self.widget.bbox(self.symbol)
                self.x_off=(x1-x2)/2
                self.y_off=(y1-y2)/2
        else:
            # no callback, don't bother with drag'n'drop
            self.widget.drag = 0
            self.dnd_end(None, None)

    def PVT_highlight(self, item):
        # mark focused item.  note that this code recreates the
        # rectangle for each update, but that's fast enough for
        # this case.
        sw = self.widget
        bbox = sw.bbox(item)
        sw.delete("highlight")
        if bbox:
            i = sw.create_rectangle(bbox, fill="white", tag="highlight")
            sw.lower(i, item)

    def PVT_key(self, event):
        sw = self.widget
        item = sw.focus()
        if not item: return
        insert = sw.index(item, TK.INSERT)

        if event.char >= " ":
            # printable character
            if sw.tk.call(sw._w, 'select', 'item'):
                sw.dchars(item, TK.SEL_FIRST, TK.SEL_LAST)
                sw.select_clear()
            sw.insert(item, "insert", event.char)
            self.PVT_highlight(item)
        elif event.keysym == "BackSpace":
            if sw.tk.call(sw._w, 'select', 'item'):
                sw.dchars(item, TK.SEL_FIRST, TK.SEL_LAST)
                sw.select_clear()
            else:
                if insert > 0: sw.dchars(item, insert-1, insert)
            self.PVT_highlight(item)
        elif event.keysym == "Return":
            pid = self.id
            # Change the node name
            pid[0] = sw.itemcget(item, "text")
            try: tkNodeEdit.updateNode(pid)
            except: pass
            sw.focus('')
            sw.select_clear()
            sw.delete("highlight")
            #x1, y1, x2, y2 = sw.bbox(self.symbol, self.label)
            #sw.coords(sw.cursor_box, x1-1, y1-1, x2+1, y2+1)
            x1, y1, x2, y2 = sw.bbox(self.label)
            sw.coords(sw.cursor_box, x1, y1, x2, y2)

            if pid[3] == 'Zone_t':
                if CTK.__MAINTREE__ <= 0: CTK.display(CTK.t)
                bases = CTK.t[2][1:]
                active = []
                for b in bases:
                    baseName = b[0]
                    for z in b[2]:
                        if z[3] == 'Zone_t':
                            if id(z) == id(pid):
                                zoneName = baseName+Internal.SEP1+pid[0]
                                i = CPlot.getCPlotNumber(CTK.t, baseName, pid[0])
                                active.append((i, zoneName))
                CPlot.setZoneNames(active)
                (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)

            elif pid[3] == 'CGNSBase_t':
                if CTK.__MAINTREE__ <= 0: CTK.display(CTK.t)
                zones = Internal.getNodesFromType1(pid, 'Zone_t')
                bases = CTK.t[2][1:]
                active = []
                dnz = CPlot.updateCPlotGlobalNumbering(CTK.t)
                for z in zones:
                    zoneName = pid[0]+Internal.SEP1+z[0]
                    #i = CPlot.getCPlotNumber(CTK.t, pid[0], z[0])
                    i = dnz[pid[0]][z[0]]
                    active.append((i, zoneName))
                CPlot.setZoneNames(active)
                (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        elif event.keysym == "Right":
            sw.icursor(item, insert+1)
            sw.select_clear()
        elif event.keysym == "Left":
            sw.icursor(item, insert-1)
            sw.select_clear()

    def PVT_clickEdit(self, event):
        self.widget.move_cursor(self)
        if self.widget.type(TK.CURRENT) != "text": return
        self.PVT_highlight(TK.CURRENT)
        # move focus to item
        self.widget.focus_set() # move focus to canvas
        self.widget.focus(TK.CURRENT) # set focus to text item
        self.widget.select_from(TK.CURRENT, 0)
        self.widget.select_to(TK.CURRENT, TK.END)

    def PVT_clickRight(self, event):
        if CTK.t == []: return
        self.widget.focus('')
        self.widget.select_clear()
        self.widget.delete("highlight")

        self.widget.move_cursor(self)
        # Toggle activate status of a zone
        pid = self.id
        if pid[3] == 'CGNSTree_t':
            if CTK.__MAINTREE__ <= 0: CTK.display(CTK.t)
            bases = Internal.getBases(pid)
            activated = []
            s = -1
            for b in bases:
                baseName = b[0]
                nodes = Internal.getNodesFromType1(b, 'Zone_t')
                if nodes != []:
                    noz = CPlot.getCPlotNumber(CTK.t, baseName, nodes[0][0])
                    if s == -1:
                        sp = CPlot.getActiveStatus(noz)
                        if sp == 0: s = 1
                        else: s = 0
                        break
            nodes = Internal.getZones(CTK.t)
            for no in range(len(nodes)): activated.append((no, s))
            if s == 0: CTK.TXT.insert('START', 'Tree deactivated.\n')
            elif s == 1: CTK.TXT.insert('START', 'Tree activated.\n')
            CPlot.setActiveZones(activated)
        elif pid[3] == 'Zone_t':
            if CTK.__MAINTREE__ <= 0: CTK.display(CTK.t)
            bases = Internal.getBases(CTK.t)
            for b in bases:
                zones = Internal.getNodesFromType1(b, 'Zone_t')
                for z in zones:
                    if id(z) == id(pid): ret = b; break
            noz = CPlot.getCPlotNumber(CTK.t, ret[0], pid[0])
            active = CPlot.getActiveStatus(noz)
            if active == 1:
                CPlot.setActiveZones([(noz,0)])
                CTK.TXT.insert('START', 'Zone '+pid[0]+' deactivated.\n')
            else:
                CPlot.setActiveZones([(noz,1)])
                CTK.TXT.insert('START', 'Zone '+pid[0]+' activated.\n')
        elif pid[3] == 'CGNSBase_t':
            if CTK.__MAINTREE__ <= 0: CTK.display(CTK.t)
            baseName = pid[0]
            nodes = Internal.getNodesFromType1(pid, 'Zone_t')
            activated = []
            s = -1
            dnz = CPlot.updateCPlotGlobalNumbering(CTK.t)
            if nodes != []:
                #noz = CPlot.getCPlotNumber(CTK.t, baseName, nodes[0][0])
                noz = dnz[baseName][nodes[0][0]]
                s = CPlot.getActiveStatus(noz)
                if s == 0: s = 1
                else: s = 0
            for z in nodes:
                #noz = CPlot.getCPlotNumber(CTK.t, baseName, z[0])
                noz = dnz[baseName][z[0]]
                activated.append( (noz, s) )
            if s == 0:
                CTK.TXT.insert('START', 'Base '+baseName+' deactivated.\n')
            elif s == 1:
                CTK.TXT.insert('START', 'Base '+baseName+' activated.\n')
            CPlot.setActiveZones(activated)
        elif pid[3] == 'Family_t':
            if CTK.__MAINTREE__ <= 0: CTK.display(CTK.t)
            isFamilyBC = Internal.getNodeFromType1(pid, 'FamilyBC_t')
            if isFamilyBC is not None: return
            nodes = C.getFamilyZones(CTK.t, pid[0])
            zones = []
            for n in nodes: # ne garde que les zones
                if n[3] == 'Zone_t': zones.append(n)

            activated = []
            active = -2
            dnz = CPlot.updateCPlotGlobalNumbering(CTK.t)
            for z in zones:
                base, c = Internal.getParentOfNode(CTK.t, z)
                #noz = CPlot.getCPlotNumber(CTK.t, base[0], z[0])
                noz = dnz[base[0]][z[0]]
                if active == -2:
                    active = CPlot.getActiveStatus(noz)
                if active == 1: activated.append( (noz, 0) )
                else: activated.append( (noz, 1) )
            if active == 1:
                CTK.TXT.insert('START', 'Family '+pid[0]+' deactivated.\n')
            elif active == 0:
                CTK.TXT.insert('START', 'Family '+pid[0]+' activated.\n')
            CPlot.setActiveZones(activated)

    def PVT_clickMultipleSelect(self, event):
        self.widget.move_cursor(self)
        self.PVT_displayNode(False)

    def PVT_clickSelect(self, event):
        if CTK.t == []: return
        self.widget.focus('')
        self.widget.select_clear()
        self.widget.delete("highlight")
        self.widget.move_cursor(self)
        self.PVT_displayNode(True)

    def PVT_displayNode(self, clear=False):
        pid = self.id
        try: tkNodeEdit.updateNode(pid)
        except: pass

        if pid[3] == 'CGNSLibraryVersion_t':
            v = pid[1]
            if isinstance(v, float): v = strFormat(v)
            if isinstance(v, numpy.ndarray): v = str(v[0])
            CTK.TXT.insert('START', v+'\n')

        elif pid[3] == 'ZoneType_t':
            v = Internal.getValue(pid)
            CTK.TXT.insert('START', v+'\n')

        elif pid[3] == 'CGNSTree_t':
            if CTK.__MAINTREE__ <= 0: CTK.display(CTK.t)
            bases = Internal.getBases(pid)
            s = -1

            for b in bases:
                baseName = b[0]
                nodes = Internal.getNodesFromType1(b, 'Zone_t')
                if nodes != []:
                    noz = CPlot.getCPlotNumber(CTK.t, baseName, nodes[0][0])
                    if s == -1:
                        sp = CPlot.getSelectedStatus(noz)
                        if sp == 0: s = 1
                        else: s = 0
                        break

            nodes = Internal.getZones(CTK.t)
            if clear: CPlot.unselectAllZones(); s = 1 # force select
            selected = []
            for no in range(len(nodes)): selected.append((no, s))

            CPlot.setSelectedZones(selected)
            if s == 1: CTK.TXT.insert('START', 'Tree selected.\n')
            else: CTK.TXT.insert('START', 'Tree unselected.\n')

        elif pid[3] == 'CGNSBase_t':
            if CTK.__MAINTREE__ <= 0: CTK.display(CTK.t)
            baseName = pid[0]
            nodes = Internal.getNodesFromType1(pid, 'Zone_t')
            s = 1
            dnz = CPlot.updateCPlotGlobalNumbering(CTK.t)
            if nodes != []:
                #noz = CPlot.getCPlotNumber(CTK.t, baseName, nodes[0][0])
                noz = dnz[baseName][nodes[0][0]]
                s = CPlot.getSelectedStatus(noz)
                if s == 0: s = 1
                else: s = 0
            if clear: CPlot.unselectAllZones(); s = 1 # force select
            selected = []
            for z in nodes:
                #noz = CPlot.getCPlotNumber(CTK.t, baseName, z[0])
                noz = dnz[baseName][z[0]]
                selected.append((noz, s))
            CPlot.setSelectedZones(selected)
            if s == 1:
                CTK.TXT.insert('START', 'Base '+baseName+' selected.\n')
            else:
                CTK.TXT.insert('START', 'Base '+baseName+' unselected.\n')

        elif pid[3] == 'IndexRange_t':
            if pid[1].shape == (3,2):
                win = Internal.range2Window(pid[1])
                string = ''
                for i in win: string += '%d '%i
            else: string = str(pid[1])
            string += '\n'
            CTK.TXT.insert('START', string)

        elif pid[3] == '"int"':
            string = str(Internal.getValue(pid))+'\n'
            CTK.TXT.insert('START', string)

        elif pid[3] == 'GridLocation_t':
            v = Internal.getValue(pid)
            CTK.TXT.insert('START', v+'\n')

        elif pid[3] == 'GridConnectivityType_t':
            v = Internal.getValue(pid)
            CTK.TXT.insert('START', v+'\n')

        elif pid[3] == 'GridConnectivity_t':
            v = Internal.getValue(pid)
            CTK.TXT.insert('START', v+'\n')

        elif pid[3] == 'Zone_t':
            if CTK.__MAINTREE__ <= 0: CTK.display(CTK.t)
            ret = Internal.getParentOfNode(CTK.t, pid)
            noz = CPlot.getCPlotNumber(CTK.t, ret[0][0], pid[0])
            s = CPlot.getSelectedStatus(noz)
            if clear: CPlot.unselectAllZones(); s = 0 # force select
            if s == 0:
                CPlot.setSelectedZones([(noz,1)])
                txt = 'Zone '+pid[0]+' selected ('
            else:
                CPlot.setSelectedZones([(noz,0)])
                txt = 'Zone '+pid[0]+' unselected ('
            dim = pid[1]
            if dim.shape[0] == 3:
                txt += str(dim[0,0])+';'+str(dim[1,0])+';'+str(dim[2,0])+')\n'
            elif dim.shape[0] == 2:
                txt += str(dim[0,0])+';'+str(dim[1,0])+')\n'
            else:
                txt += str(dim[0,0])+')\n'
            CTK.TXT.insert('START', txt)

        elif pid[3] == 'DataArray_t':
            v = pid[1]
            txt = ''
            if isinstance(v, numpy.ndarray):
                if v.dtype == 'c': txt = Internal.getValue(pid)
                else:
                    pt = v.ravel('k'); size = pt.size
                    txt += str(v.shape)+': '
                    if size > 0: txt += strFormat(pt[0])
                    if size > 1: txt += ' ' + strFormat(pt[1])
                    if size > 2: txt += ' ' + strFormat(pt[2])
                    if size > 3: txt += ' ' + strFormat(pt[3])
                    if size > 4: txt += ' ' + strFormat(pt[4])
                    if size > 5: txt += ' ' + strFormat(pt[5])
                    if size > 6: txt += '...'
            else: txt += str(v)
            CTK.TXT.insert('START', txt+'\n')
            ret = Internal.getParentOfNode(CTK.t, pid)
            cont = ret[0]
            if cont[3] == 'FlowSolution_t':
                gp = Internal.getNodeFromType1(cont, 'GridLocation_t')
                field = None
                if cont[0] == Internal.__FlowSolutionNodes__:
                    field = pid[0]
                elif cont[0] == Internal.__FlowSolutionCenters__:
                    field = 'centers:'+pid[0]
                elif gp is not None and Internal.getValue(gp) == 'CellCenter':
                    Internal.__FlowSolutionCenters__ = cont[0]
                    field = 'centers:'+pid[0]
                    if 'tkContainers' in CTK.TKMODULES: CTK.TKMODULES['tkContainers'].updateApp()
                else:
                    Internal.__FlowSolutionNode__ = cont[0]
                    field = pid[0]
                    if 'tkContainers' in CTK.TKMODULES: CTK.TKMODULES['tkContainers'].updateApp()

                zvars = C.getVarNames(CTK.t)[0]
                ifield = 0; lenvars = 0
                for i in zvars:
                    if i != 'CoordinateX' and i != 'CoordinateY' and i != 'CoordinateZ':
                        lenvars += 1
                for i in zvars:
                    if i == field: break
                    if i != 'CoordinateX' and i != 'CoordinateY' and i != 'CoordinateZ':
                        ifield += 1
                CPlot.setState(mode=3, scalarField=ifield)

        elif pid[3] == 'BC_t':
            v = Internal.getValue(pid)
            CTK.TXT.insert('START', 'BC type: %s\n'%v)

        elif pid[3] == 'GoverningEquations_t':
            v = Internal.getValue(pid)
            CTK.TXT.insert('START', v+'\n')

        elif pid[3] == 'TurbulenceModel_t':
            v = Internal.getValue(pid)
            CTK.TXT.insert('START', v+'\n')

        elif pid[3] == 'GridConnectivity1to1_t':
            v = Internal.getValue(pid)
            CTK.TXT.insert('START', '1to1 with %s\n'%v)

        elif pid[3] == 'GridConnectivity_t':
            v = Internal.getValue(pid)
            CTK.TXT.insert('START', 'Cnx with %s\n'%v)

        elif pid[3] == 'IndexArray_t':
            v = pid[1]
            txt = ''
            if isinstance(v, numpy.ndarray):
                txt += str(v.shape)+': '
                p = v.ravel('k'); ln = p.size
                if ln > 0: txt += strFormat(p[0])
                if ln > 1: txt += ' ' + strFormat(p[1])
                if ln > 2: txt += ' ' + strFormat(p[2])
                if ln > 3: txt += '...'
            else: txt += str(v)
            CTK.TXT.insert('START', txt+'\n')

        elif pid[3] == '"int[IndexDimension]"':
            CTK.TXT.insert('START', str(pid[1])+'\n')

        elif pid[3] == 'GridCoordinates_t':
            CTK.TXT.insert('START', 'Displaying mesh\n')
            CPlot.setMode(0)

        elif pid[3] == 'FlowSolution_t':
            CTK.TXT.insert('START', 'Displaying '+pid[0]+'\n')
            contName = pid[0]
            CPlot.setMode(3)
            gp = Internal.getNodeFromType1(pid, 'GridLocation_t')
            if gp is not None:
                if Internal.getValue(gp) == 'CellCenter' and Internal.__FlowSolutionCenters__ != contName:
                    Internal.__FlowSolutionCenters__ = contName
                    if 'tkContainers' in CTK.TKMODULES: CTK.TKMODULES['tkContainers'].updateApp()
                    CTK.display(CTK.t)
                elif Internal.__FlowSolutionNodes__ != contName:
                    Internal.__FlowSolutionNodes__ = contName
                    if 'tkContainers' in CTK.TKMODULES: CTK.TKMODULES['tkContainers'].updateApp()
                    CTK.display(CTK.t)

        elif pid[3] == 'FamilyName_t' or pid[3] == 'AdditionalFamilyName_t':
            v = Internal.getValue(pid)
            CTK.TXT.insert('START', v+'\n')

        elif pid[3] == 'FamilyBC_t':
            v = Internal.getValue(pid)
            CTK.TXT.insert('START', v+'\n')

        elif pid[3] == 'Family_t':
            isFamilyBC = Internal.getNodeFromType1(pid, 'FamilyBC_t')
            if isFamilyBC is not None:
                CTK.TXT.insert('START', 'Family '+pid[0]+' is a familyBC.\n')
            else: # Zone family
                if CTK.__MAINTREE__ <= 0: CTK.display(CTK.t)
                nodes = C.getFamilyZones(CTK.t, pid[0])
                zones = []
                for n in nodes: # ne garde que les zones
                    if n[3] == 'Zone_t': zones.append(n)

                selected = []
                s = 1
                dnz = CPlot.updateCPlotGlobalNumbering(CTK.t)
                if zones != []:
                    z = zones[0]
                    base, c = Internal.getParentOfNode(CTK.t, z)
                    #noz = CPlot.getCPlotNumber(CTK.t, base[0], z[0])
                    noz = dnz[base[0]][z[0]]
                    s = CPlot.getSelectedStatus(noz)
                    if s == 0: s = 1
                    else: s = 0
                if clear: s = 1 # force select
                for z in zones:
                    base, c = Internal.getParentOfNode(CTK.t, z)
                    #noz = CPlot.getCPlotNumber(CTK.t, base[0], z[0])
                    noz = dnz[base[0]][z[0]]
                    selected.append((noz, s))
                if clear: CPlot.unselectAllZones()
                CPlot.setSelectedZones(selected)
                if s == 1:
                    CTK.TXT.insert('START', 'Family of zones '+pid[0]+' selected.\n')
                else:
                    CTK.TXT.insert('START', 'Family of zones '+pid[0]+' unselected.\n')

        elif pid[3] == 'Descriptor_t':
            v = Internal.getValue(pid)
            CTK.TXT.insert('START', v+'\n')

        elif pid[3] == 'Elements_t':
            val = pid[1]
            connectType = val[0]; boundary = val[1]
            v = 'Connectivity of type '
            name, nnodes = Internal.eltNo2EltName(connectType)
            v += name

            if boundary > 0: v += ' (boundary)'
            CTK.TXT.insert('START', v+'.\n')

        else:
            v = 'Node type: %s'%pid[3]
            CTK.TXT.insert('START', v+'\n')

#------------------------------------------------------------------------------
class Tree(TK.Canvas):
    def __init__(self, master, root_id, root_label='',
                 get_contents_callback=None, dist_x=8+CTK.FONTSIZE,
                 dist_y=8+CTK.FONTSIZE, text_offset=12,
                 line_flag=1, expanded_icon=None,
                 collapsed_icon=None, regular_icon=None, plus_icon=None,
                 minus_icon=None, node_class=Node, drop_callback=None,
                 *args, **kw_args):
        # pass args to superclass
        TK.Canvas.__init__(self, master, *args, **kw_args)
        # this allows to subclass Node and pass our class in
        self.node_class = node_class
        # keep track of node bindings
        self.bindings = {}
        # cheap mutex spinlock
        self.spinlock = 0
        # flag to see if there's been any d&d dragging
        self.drag = 0
        # default images (BASE64-encoded GIF files)
        if expanded_icon is None:
            self.expanded_icon=TK.PhotoImage(
                data='R0lGODlhEAANAKIAAAAAAMDAwICAgP//////ADAwMAAAAAAA' \
                'ACH5BAEAAAEALAAAAAAQAA0AAAM6GCrM+jCIQamIbw6ybXNSx3GVB' \
                'YRiygnA534Eq5UlO8jUqLYsquuy0+SXap1CxBHr+HoBjoGndDpNAAA7')
        else:
            self.expanded_icon = expanded_icon
        if collapsed_icon is None:
            self.collapsed_icon = TK.PhotoImage(
                data='R0lGODlhDwANAKIAAAAAAMDAwICAgP//////ADAwMAAAAAAA' \
                'ACH5BAEAAAEALAAAAAAPAA0AAAMyGCHM+lAMMoeAT9Jtm5NDKI4Wo' \
                'FXcJphhipanq7Kvu8b1dLc5tcuom2foAQQAyKRSmQAAOw==')
        else:
            self.collapsed_icon = collapsed_icon
        if regular_icon is None:
            self.regular_icon = TK.PhotoImage(
                data='R0lGODlhCwAOAJEAAAAAAICAgP///8DAwCH5BAEAAAMALAAA' \
                'AAALAA4AAAIphA+jA+JuVgtUtMQePJlWCgSN9oSTV5lkKQpo2q5W+' \
                'wbzuJrIHgw1WgAAOw==')
        else:
            self.regular_icon = regular_icon
        if plus_icon is None:
            self.plus_icon = TK.PhotoImage(
                data='R0lGODdhCQAJAPEAAAAAAH9/f////wAAACwAAAAACQAJAAAC' \
                'FIyPoiu2sJyCyoF7W3hxz850CFIA\nADs=')
        else:
            self.plus_icon = plus_icon
        if minus_icon is None:
            self.minus_icon = TK.PhotoImage(
                data='R0lGODdhCQAJAPEAAAAAAH9/f////wAAACwAAAAACQAJAAAC' \
                'EYyPoivG614LAlg7ZZbxoR8UADs=')
        else:
            self.minus_icon = minus_icon

        self.expanded_zone = TK.PhotoImage(data="""
R0lGODlhEAANAMIEAAAAAICAgJKC1e/s+f///////////////yH+EUNyZWF0ZWQgd2l0aCBHSU1Q
ACH5BAEAAAQALAAAAAAQAA0AAAM6SBrM+pCEQWmIbw6xbXNSx3GVBYRiygXA534Cq5UlO8jUqLYs
quuy0+SXap1CxBHr+HoBjoSndDpNAAA7
""")

        self.collapsed_zone = TK.PhotoImage(data="""
R0lGODlhDwANAMIEAAAAAICAgJKC1ebi9f///////////////yH+EUNyZWF0ZWQgd2l0aCBHSU1Q
ACH5BAEAAAQALAAAAAAPAA0AAAMySBTM+jAMMUeAT9Jtm5NDKI4WoFXcFphhipanq7Kvu8b1dLc5
tcuom2foARAAyKRSmQAAOw==
""")

        self.expanded_base = TK.PhotoImage(data="""
R0lGODlhEAANAMIEAAAAAICAgKq4R8HqwP///////////////yH+EUNyZWF0ZWQgd2l0aCBHSU1Q
ACH5BAEAAAQALAAAAAAQAA0AAAM6SBrM+pCEQWmIbw6xbXNSx3GVBYRiygXA534Cq5UlO8jUqLYs
quuy0+SXap1CxBHr+HoBjoSndDpNAAA7
""")

        self.collapsed_base = TK.PhotoImage(data="""
R0lGODlhDwANAMIEAAAAAICAgECiNpHZjv///////////////yH+EUNyZWF0ZWQgd2l0aCBHSU1Q
ACH5BAEAAAQALAAAAAAPAA0AAAMySBTM+jAMMUeAT9Jtm5NDKI4WoFXcFphhipanq7Kvu8b1dLc5
tcuom2foARAAyKRSmQAAOw==
""")

        self.expanded_family = TK.PhotoImage(data="""
R0lGODlhEAANAMIEAAAAAMRPT4CAgOWurv///////////////yH+EUNyZWF0ZWQgd2l0aCBHSU1Q
ACH5BAEAAAQALAAAAAAQAA0AAAM6SCrM+pCIQamIb46wbXNSx3GVBYRiygnA534Bq5UlO8jUqLYs
quuy0+SXap1CxBHr+HoBjoSndDpNAAA7
""")

        self.collapsed_family = TK.PhotoImage(data="""
R0lGODlhDwANAMIEAAAAAMRPT4CAgOa0tP///////////////yH+EUNyZWF0ZWQgd2l0aCBHSU1Q
ACH5BAEAAAQALAAAAAAPAA0AAAMySCTM+lCMMIeAT9Jtm5NDKI4WoFXcJphhipanq7Kvu8b1dLc5
tcuom2foARAAyKRSmQAAOw==
""")

        self.expanded_user = TK.PhotoImage(data="""
R0lGODlhEAANAMIEAAAAAICAgP/Gwf/l4////////////////yH+EUNyZWF0ZWQgd2l0aCBHSU1Q
ACH5BAEAAAQALAAAAAAQAA0AAAM6SBrM+pCEQWmIbw6xbXNSx3GVBYRiygXA534Cq5UlO8jUqLYs
quuy0+SXap1CxBHr+HoBjoSndDpNAAA7
""")

        self.collapsed_user = TK.PhotoImage(data="""
R0lGODlhDwANAMIEAAAAAICAgP+pof7Rzf///////////////yH+EUNyZWF0ZWQgd2l0aCBHSU1Q
ACH5BAEAAAQALAAAAAAPAA0AAAMySBTM+jAMMUeAT9Jtm5NDKI4WoFXcFphhipanq7Kvu8b1dLc5
tcuom2foARAAyKRSmQAAOw==
""")

        # horizontal distance that subtrees are indented
        self.dist_x = dist_x
        # vertical distance between rows
        self.dist_y = dist_y
        # how far to offset text label
        self.text_offset = text_offset
        # flag controlling connecting line display
        self.line_flag = line_flag
        # called just before subtree expand/collapse
        self.get_contents_callback = get_contents_callback
        # called after drag'n'drop
        self.drop_callback = drop_callback
        # create root node to get the ball rolling
        self.root = node_class(parent_node=None, label=root_label,
                               id=root_id, expandable_flag=1,
                               collapsed_icon=self.collapsed_icon,
                               expanded_icon=self.expanded_icon,
                               x=dist_x, y=dist_y, parent_widget=self)
        # configure for scrollbar(s)
        x1, y1, x2, y2 = self.bbox('all')
        self.configure(scrollregion=(x1, y1, x2+5, y2+5))
        # add a cursor
        self.cursor_box = self.create_rectangle(0, 0, 0, 0,
                                                outline='sky blue')
        #fill='sky blue')
        self.move_cursor(self.root)
        # make it easy to point to control
        self.bind('<Enter>', self.PVT_mousefocus)
        # totally arbitrary yet hopefully intuitive default keybindings
        # stole 'em from ones used by microsoft tree control
        # page-up/page-down
        self.bind('<Next>', self.pagedown)
        self.bind('<Prior>', self.pageup)
        # arrow-up/arrow-down
        self.bind('<Down>', self.next)
        self.bind('<Up>', self.prev)
        # arrow-left/arrow-right
        #self.bind('<Left>', self.ascend)
        # (hold this down and you expand the entire tree)
        #self.bind('<Right>', self.descend)
        # home/end
        self.bind('<Home>', self.first)
        self.bind('<End>', self.last)
        # space bar
        self.bind('<Key-space>', self.toggle)
        # erase
        self.bind('<Delete>', self.PVT_erase)
        # Resize canvas
        self.bind('<Control-e>', self.expandCanvas)
        self.bind('<Control-r>', self.shrinkCanvas)

    # ----- PRIVATE METHODS (prefixed with "PVT_") -----
    # these methods are subject to change, so please try not to use them
    def PVT_mousefocus(self, event):
        """Soak up event argument when moused-over"""
        self.focus_set()

    def PVT_erase(self, event):
        self.move_cursor(self.pos)
        # Erase node
        CTK.saveTree()
        id = self.pos.id

        if id[3] == 'CGNSTree_t': # efface tout l'arbre
            del id[2]
            CTK.t = [id[0], None, [], 'CGNSTree_t']
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.display(CTK.t); return
        elif id[3] == 'CGNSBase_t': # efface une base
            ret = Internal.getParentOfNode(CTK.t, id)
            deletedZoneNames = []
            for z in id[2]:
                if z[3] == 'Zone_t': deletedZoneNames.append(id[0]+Internal.SEP1+z[0])
            if Internal.isStdNode(ret[0]) >= 0: del ret[0][ret[1]]
            else: del ret[0][2][ret[1]]
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            if deletedZoneNames != []:
                CPlot.delete(deletedZoneNames)
                CPlot.render()
        elif id[3] == 'Zone_t': # efface une zone
            ret = Internal.getParentOfNode(CTK.t, id)
            deletedZoneNames = [ret[0][0]+Internal.SEP1+id[0]]
            del ret[0][2][ret[1]]
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CPlot.delete(deletedZoneNames)
            CPlot.render()
        else: # efface un autre type de noeuds
            ret = Internal.getParentOfNode(CTK.t, id)
            if Internal.isStdNode(ret[0]) >= 0: del ret[0][ret[1]]
            else: del ret[0][2][ret[1]]
        sav = self.pos.prev_visible()
        self.pos.delete()
        self.move_cursor(sav)

    # ----- PUBLIC METHODS -----
    def tag_bind(self, tag, seq, *args, **kw_args):
        """Keep track of callback bindings so we can delete them later. I
        shouldn't have to do this!!!!"""
        # pass args to superclass
        aargs = (self, tag, seq)+args
        func_id = TK.Canvas.tag_bind(*aargs, **kw_args)
        # save references
        self.bindings[tag] = self.bindings.get(tag, [])+[(seq, func_id)]

    def add_list(self, list=None, name=None, id=None, flag=0,
                 expanded_icon=None, collapsed_icon=None):
        """Add node construction info to list"""
        n = Struct()
        n.name = name
        n.id = id
        n.flag = flag
        if collapsed_icon: n.collapsed_icon = collapsed_icon
        else:
            if flag:
                # it's expandable, use closed folder icon
                n.collapsed_icon = self.collapsed_icon
            else:
                # it's not expandable, use regular file icon
                n.collapsed_icon = self.regular_icon
        if flag:
            if expanded_icon:
                n.expanded_icon = expanded_icon
            else:
                n.expanded_icon = self.expanded_icon
        else:
            # not expandable, don't need an icon
            n.expanded_icon = None
        if list is None: list = []
        list.append(n)
        return list

    def add_node(self, name=None, id=None, flag=0, expanded_icon=None,
                 collapsed_icon=None):
        """Add a node during get_contents_callback()"""
        if id is not None and len(id) >= 4:
            if id[3] == 'Zone_t':
                collapsed_icon = self.collapsed_zone
                expanded_icon = self.expanded_zone
            elif id[3] == 'CGNSBase_t':
                collapsed_icon = self.collapsed_base
                expanded_icon = self.expanded_base
            elif id[3] == 'Family_t' or id[3] == 'FamilyName_t' or id[3] == 'AdditionalFamilyName_t':
                collapsed_icon = self.collapsed_family
                expanded_icon = self.expanded_family
            elif id[3] == 'UserDefinedData_t' and len(id[2])>0:
                collapsed_icon = self.collapsed_user
                expanded_icon = self.expanded_user
        self.add_list(self.new_nodes, name, id, flag, expanded_icon,
                      collapsed_icon)

    def find_full_id(self, search):
        """Search for a node"""
        return self.root.PVT_find(search)

    def cursor_node(self):
        """Return node under cursor"""
        return self.pos

    def see(self, *items):
        """Scroll (in a series of nudges) so items are visible"""
        x1, y1, x2, y2 = self.bbox(*items)

        while x2 > self.canvasx(0)+self.winfo_width():
            old = self.canvasx(0)
            self.xview('scroll', 1, 'units')
            # avoid endless loop if we can't scroll
            if old == self.canvasx(0): break

        while y2 > self.canvasy(0)+self.winfo_height():
            old = self.canvasy(0)
            self.yview('scroll', 1, 'units')
            if old == self.canvasy(0): break

        # done in this order to ensure upper-left of object is visible
        while x1 < self.canvasx(0):
            old = self.canvasx(0)
            self.xview('scroll', -1, 'units')
            if old == self.canvasx(0): break
        while y1 < self.canvasy(0):
            old = self.canvasy(0)
            self.yview('scroll', -1, 'units')
            if old == self.canvasy(0): break

    def move_cursor(self, node):
        """Move cursor to node"""
        self.pos = node
        #x1, y1, x2, y2 = self.bbox(node.symbol, node.label)
        #self.coords(self.cursor_box, x1-1, y1-1, x2+1, y2+1)
        x1, y1, x2, y2 = self.bbox(node.label)
        self.coords(self.cursor_box, x1, y1, x2, y2)
        self.see(node.symbol, node.label)

    def toggle(self, event=None):
        """Expand/collapse subtree"""
        self.pos.toggle_state()

    def next(self, event=None):
        """Move to next lower visible node"""
        self.move_cursor(self.pos.next_visible())

    def prev(self, event=None):
        """Move to next higher visible node"""
        self.move_cursor(self.pos.prev_visible())

    def ascend(self, event=None):
        """Move to immediate parent"""
        if self.pos.parent_node:
            # move to parent
            self.move_cursor(self.pos.parent_node)

    def descend(self, event=None):
        """Move right, expanding as we go"""
        if self.pos.expandable_flag:
            self.pos.expand()
            if self.pos.child_nodes:
                # move to first subnode
                self.move_cursor(self.pos.child_nodes[0])
                return
        # if no subnodes, move to next sibling
        self.next()

    def first(self, event=None):
        """Go to root node"""
        # move to root node
        self.move_cursor(self.root)

    def last(self, event=None):
        """Go to last visible node"""
        # move to bottom-most node
        self.move_cursor(self.root.PVT_last())

    def pageup(self, event=None):
        """Previous page"""
        n = self.pos
        j = self.winfo_height()//self.dist_y
        for i in range(j-3):
            n = n.prev_visible()
        self.yview('scroll', -1, 'pages')
        self.move_cursor(n)

    def pagedown(self, event=None):
        """Next page"""
        n = self.pos
        j = self.winfo_height()//self.dist_y
        for i in range(j-3):
            n = n.next_visible()
        self.yview('scroll', 1, 'pages')
        self.move_cursor(n)

    def clear(self):
        for cur_node in self.root.children():
            cur_node.delete()
        self.root.collapse()

    def expandCanvas(self, event=None):
        aw = int(self.cget("width"))+10 # 230
        ah = int(self.cget("height"))+10 # 210
        self.config(width=aw, height=ah)

    def shrinkCanvas(self, event=None):
        aw = max(int(self.cget("width"))-10,230) # 230
        ah = max(int(self.cget("height"))-10,210) # 210
        self.config(width=aw, height=ah)

    # ----- functions for drag'n'drop support -----
    def where(self, event):
        """Determine drag location in canvas coordinates. event.x & event.y
        don't seem to be what we want."""
        # where the corner of the canvas is relative to the screen:
        x_org = self.winfo_rootx()
        y_org = self.winfo_rooty()
        # where the pointer is relative to the canvas widget,
        # including scrolling
        x = self.canvasx(event.x_root-x_org)
        y = self.canvasy(event.y_root-y_org)
        return x, y

    def dnd_accept(self, source, event):
        """Accept dnd messages, i.e. we're a legit drop target, and we do
        implement d&d functions."""
        self.target = None
        return self

    def dnd_enter(self, source, event):
        """Get ready to drag or drag has entered widget (create drag
        object)"""
        # this flag lets us know there's been drag motion
        self.drag = 1
        x, y = self.where(event)
        x1, y1, x2, y2=source.widget.bbox(source.symbol, source.label)
        dx, dy = x2-x1, y2-y1
        # create dragging icon
        if source.expanded_flag:
            self.dnd_symbol = self.create_image(x, y,
                                                image=source.expanded_icon)
        else:
            self.dnd_symbol = self.create_image(x, y,
                                                image=source.collapsed_icon)
        self.dnd_label = self.create_text(x+self.text_offset, y,
                                          text=source.get_label(),
                                          justify='left',
                                          anchor='w',
                                          font=CTK.TEXTFONT)

    def dnd_motion(self, source, event):
        """Move drag icon"""
        self.drag = 1
        x, y = self.where(event)
        x1, y1, x2, y2 = self.bbox(self.dnd_symbol, self.dnd_label)
        self.move(self.dnd_symbol, x-x1+source.x_off, y-y1+source.y_off)
        self.move(self.dnd_label, x-x1+source.x_off, y-y1+source.y_off)

    def dnd_leave(self, source, event):
        """Finish dragging or drag has left widget (destroy drag object)"""
        self.delete(self.dnd_symbol)
        self.delete(self.dnd_label)

    def dnd_commit(self, source, event):
        """Object has been dropped here"""
        # call our own dnd_leave() to clean up
        self.dnd_leave(source, event)
        # process pending events to detect target node
        # update_idletasks() doesn't do the trick if source & target are
        # on  different widgets
        self.update()
        if not self.target:
            # no target node
            return
        # we must update data structures based on the drop
        if self.drop_callback:
            try:
                # called with dragged node and target node
                # this is where a file manager would move the actual file
                # it must also move the nodes around as it wishes
                self.drop_callback(source, self.target)
            except:
                report_callback_exception()

    def focusOnGivenZone(self, baseName, zoneName=None, bcName=None):
        """Focus on given zone."""
        node = self.root
        if not node.expanded(): # root is not expanded
            node.expand()
        children = node.children()
        # Look for base node
        baseNode = None
        for c in children: # bases
            if c.id[0] == baseName: baseNode = c; break
        if baseNode is None: return
        if zoneName is None:
            baseNode.widget.move_cursor(baseNode)
            baseNode.PVT_highlight(TK.CURRENT)
            return
        if not baseNode.expanded(): # baseNode is not expanded
            baseNode.expand()
        children = baseNode.children()
        # Look for zone node
        zoneNode = None
        for c in children: # zones
            if c.id[0] == zoneName: zoneNode = c; break
        if zoneNode is None:
            baseNode.widget.move_cursor(baseNode)
            baseNode.PVT_highlight(TK.CURRENT)
            return
        if bcName is None:
            zoneNode.widget.move_cursor(zoneNode)
            zoneNode.PVT_highlight(TK.CURRENT)
            return
        if not zoneNode.expanded(): # zoneNode is not expanded
            zoneNode.expand()
        children = zoneNode.children()
        # Look for ZoneBC and ZoneGC
        ZBCNode = None; ZGCNode = None
        for c in children:
            if c.id[0] == 'ZoneBC': ZBCNode = c
            if c.id[0] == 'ZoneGridConnectivity': ZGCNode = c
        # Look for bcNode in ZoneBC
        bcNode = None
        if ZBCNode is not None:
            if not ZBCNode.expanded(): ZBCNode.expand()
            children = ZBCNode.children()
            for c in children:
                if c.id[0] == bcName: bcNode = c; break
            if bcNode is not None:
                if not ZBCNode.expanded(): # zoneBC is not expanded
                    ZBCNode.expand()
                bcNode.widget.move_cursor(bcNode)
                bcNode.PVT_highlight(TK.CURRENT)
                return
        # Look for bcNode in ZoneGridConnectivity
        if ZGCNode is not None:
            if not ZGCNode.expanded(): ZGCNode.expand()
            children = ZGCNode.children()
            for c in children:
                if c.id[0] == bcName: bcNode = c; break
            if bcNode is not None:
                if not ZGCNode.expanded(): # zoneGC is not expanded
                    ZGCNode.expand()
                bcNode.widget.move_cursor(bcNode)
                bcNode.PVT_highlight(TK.CURRENT)
                return
        zoneNode.widget.move_cursor(zoneNode)
        zoneNode.PVT_highlight(TK.CURRENT)

#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, border=2, relief=CTK.FRAMESTYLE,
                           text='tkTree  [ + ]  ', font=CTK.FRAMEFONT)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.columnconfigure(0, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Expand', accelerator='Ctrl+e',
                          command=expandCanvas)
    FrameMenu.add_command(label='Shrink', accelerator='Ctrl+r',
                          command=shrinkCanvas)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    WIDGETS['frameMenu'] = FrameMenu

    # - Frame sunken -
    Frame2 = TK.Frame(Frame, border=1, relief=TK.SUNKEN)
    Frame2.columnconfigure(0, weight=1)
    Frame2.grid(sticky=TK.NSEW)

    aw = 230; ah = 210
    if 'tkTreeWidth' in CTK.PREFS: aw = int(CTK.PREFS['tkTreeWidth'])
    if 'tkTreeHeight' in CTK.PREFS: ah = int(CTK.PREFS['tkTreeHeight'])

    # - Tree -
    B = Tree(master=Frame2,
             root_id=[],
             root_label='tree',
             get_contents_callback=get_contents,
             width=aw, height=ah,
             background='White')
    B.bind('<MouseWheel>', onMouseWheel)
    B.bind('<Button-4>', onMouseWheel)
    B.bind('<Button-5>', onMouseWheel)
    B.bind('<Control-c>', onCopy)
    B.bind('<Control-x>', onCut)
    B.bind('<Control-v>', onPaste)
    B.grid(row=0, column=0, sticky=TK.NSEW)
    WIDGETS['tree'] = B
    sb = TTK.Scrollbar(Frame, width=10)
    sb.grid(row=0, column=2, sticky=TK.NS)
    B.configure(yscrollcommand=sb.set)
    sb.configure(command=B.yview)
    sb = TTK.Scrollbar(Frame, orient=TK.HORIZONTAL, width=10)
    sb.grid(row=1, column=0, sticky=TK.EW)
    B.configure(xscrollcommand=sb.set)
    sb.configure(command=B.xview)

    # - update button -
    #B = TK.Button(Frame, text="Update tree", command=updateApp)
    #B.grid(row=2, column=0, sticky=TK.EW)

#==============================================================================
def onMouseWheel(event):
    tree = WIDGETS['tree']
    if event.num == 5 or event.delta == -120:
        tree.yview('scroll', +1, 'units')
    elif event.num == 4 or event.delta == 120:
        tree.yview('scroll', -1, 'units')

#==============================================================================
def onCopy(event):
    global BUFFER
    node = getCurrentSelectedNode()
    BUFFER = node

#==============================================================================
def onCut(event):
    CTK.saveTree()
    global BUFFER
    CTK.saveTree()
    node = getCurrentSelectedNode()
    sw = WIDGETS['tree'].cursor_node()
    sw.widget.move_cursor(sw.parent_node)
    BUFFER = node
    (p, c) = Internal.getParentOfNode(CTK.t, node)
    del p[2][c]
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.display(CTK.t)
    updateApp()

#==============================================================================
def onPaste(event):
    CTK.saveTree()
    nodep = Internal.copyTree(BUFFER)
    node = getCurrentSelectedNode()
    node[2].append(nodep)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.display(CTK.t)
    updateApp()
    sw = WIDGETS['tree'].cursor_node()
    for n in sw.children():
        if n.id[0] == nodep[0]: n.widget.move_cursor(n)

#==============================================================================
def showApp():
    WIDGETS['frame'].grid(sticky=TK.NSEW, column=1); updateApp()

#==============================================================================
def hideApp(event=None):
    WIDGETS['frame'].grid_forget()

#==============================================================================
def updateApp():
    if CTK.t != []:
        W = WIDGETS['tree']
        oldTree = W
        oldRoot = W.root
        oldChildren = oldRoot.children()
        OPENTREE = False
        OPENBASES = {} # Dict des noeuds ouverts (niv 0)
        OPENZONES = [] # Liste de dict des noeuds ouverts (niv 1)
        OPENNODES = [] # Liste de liste de dict des noeuds ouverts (niv 2)
        SELECTED = None # nom du noeud selectionne
        LEVSEL = -1 # niveau du noeud selectionne
        if len(oldChildren) > 0: OPENTREE = True
        for c in oldChildren:
            children = c.children()
            if len(children) > 0: OPENBASES[c.get_label()] = 1
            if c.widget.pos == c: SELECTED = c.get_label(); LEVSEL = 0
            ZONES = {}; NODES = []
            for d in children:
                children2 = d.children()
                if len(children2) > 0: ZONES[d.get_label()] = 1
                if d.widget.pos == d: SELECTED = d.get_label(); LEVSEL = 1
                NODESL = {}
                for e in children2:
                    children3 = e.children()
                    if len(children3) > 0: NODESL[e.get_label()] = 1
                    if e.widget.pos == e:
                        SELECTED = e.get_label(); LEVSEL = 2
                    for f in children3:
                        if f.widget.pos == f:
                            SELECTED = f.get_label(); LEVSEL = 3
                if len(children2) > 0: NODES.append(NODESL)
            if len(children) > 0: OPENNODES.append(NODES)
            if len(children) > 0: OPENZONES.append(ZONES)

        oldTree.clear()
        oldTree.root.widget.delete(oldTree.root.symbol)
        oldTree.root.id = None

        W.root = W.node_class(
            parent_node=None, label='tree',
            id=CTK.t, expandable_flag=1,
            collapsed_icon=W.collapsed_icon,
            expanded_icon=W.expanded_icon,
            x=W.dist_x, y=W.dist_y,
            parent_widget=W)

        newRoot = W.root
        if OPENTREE: newRoot.expand()
        children = newRoot.children()
        bcount = 0
        for c in children: # Bases
            if c.get_label() in OPENBASES: c.expand()
            if LEVSEL == 0 and c.get_label() == SELECTED:
                c.widget.move_cursor(c)

            children2 = c.children()
            zcount = 0
            for d in children2: # Zones
                dlabel = d.get_label()
                if dlabel in OPENZONES[bcount]: d.expand()
                if LEVSEL == 1 and dlabel == SELECTED: d.widget.move_cursor(d)
                children3 = d.children()
                for e in children3:
                    elabel = e.get_label()
                    if elabel in OPENNODES[bcount][zcount]:
                        e.expand()
                    if LEVSEL == 2 and elabel == SELECTED:
                        e.widget.move_cursor(e)
                    children4 = e.children()
                    for f in children4:
                        if LEVSEL == 3 and f.get_label() == SELECTED:
                            f.widget.move_cursor(f)
                if len(children3) > 0: zcount += 1
            if len(children2) > 0: bcount += 1

#==============================================================================
def saveApp():
    CTK.PREFS['tkTreeWidth'] = WIDGETS['tree'].cget("width")
    CTK.PREFS['tkTreeHeight'] = WIDGETS['tree'].cget("height")
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    w = 230; h = 210
    CTK.PREFS['tkTreeWidth'] = str(w)
    CTK.PREFS['tkTreeHeight'] = str(h)
    WIDGETS['tree'].config(width=w, height=h)
    CTK.savePrefFile()

#==============================================================================
def displayFrameMenu(event=None):
    WIDGETS['frameMenu'].tk_popup(event.x_root+50, event.y_root, 0)

#==============================================================================
# Renvoit le noeud selectionne dans le tktree
# Retourne le noeud sous forme de l'arbre CGNS
#==============================================================================
def getCurrentSelectedNode():
    tree = WIDGETS['tree']
    node = tree.cursor_node()
    return node.id

#==============================================================================
# Appele pour savoir ce qu'il y a dans un noeud de l'arbre widget
#==============================================================================
def get_contents(node):
    # Add displayed nodes
    id = node.id
    if len(id) != 4: return
    for i in id[2]:
        if i[2] == []: folder = 0
        else: folder = 1
        node.widget.add_node(name=i[0], id=i, flag=folder)

#==============================================================================
def expandCanvas():
    WIDGETS['tree'].expandCanvas()

def shrinkCanvas():
    WIDGETS['tree'].shrinkCanvas()

#==============================================================================
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 2:
        CTK.FILE = sys.argv[1]
        try:
            CTK.t = C.convertFile2PyTree(CTK.FILE)
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.display(CTK.t)
        except: pass

    # Main window
    (win, menu, file, tools) = CTK.minimal('tkTree '+C.__version__)

    createApp(win); updateApp(); showApp()

    # - Main loop -
    win.mainloop()
