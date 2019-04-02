# Created by Steven Miller
# Columbia University
#
# Released under the MIT License
#
# uses Python 3.6.5
# uses gdspy version 1.3.1

import gdspy
import numpy as np
import copy
from random import SystemRandom as sysrand

_twopi = 2*np.pi

# A few functions to be used
def getID():
    return sysrand().randint(1,1e10)

# Useful function to split list (in half by default)
def split_list(alist, wanted_parts=2):
	length = len(alist)
	if type(alist) is not list:
		return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts].tolist() 
				 for i in range(wanted_parts) ]
	else:
		return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] 
				 for i in range(wanted_parts) ]

def direction_str_to_num(direction):
    if direction == '+x':
        return 0
    elif direction == '-x':
        return -np.pi
    elif direction == '+y':
        return np.pi/2
    elif direction == '-y':
        return -np.pi/2
    elif direction == None:
        return None
    elif isinstance(direction, float) or isinstance(direction, int):
        return direction
    else:
        return direction

def convert_Cell(cell, port_list=None):
    '''
    Patch for converting gdspy Cell to Vgdspy Cell (with ports and obj_id)
    '''
    gdspy.current_library.cell_dict.pop(cell.name)
    c = Cell(cell.name)
    c.elements = cell.elements 
    c.labels = cell.labels
    if port_list:
        c.port_list.extend(port_list)
    return c

def port_check(cell, depth=None, terminals=False, verbose=False):
    """
    This function goes through the port_list inside a cell and places a circle 
    around each port. Identifies connected ports and end ports. 
    
    TODO: To be implemented: 
        - check of each port connection for port attributes using port.is_connected()
        - integration with the traverse_portlist() function which does a 
          depth-first-search (DFS) of the portlist to verify all ports are connected
    
    Parameters
    ----------
    cell : number
    
    depth : integer or ``None``
        If not ``None``, defines from how many reference levels to check ports.
        If None, all reference levels are included.
    terminals : boolean
        if True, only prints terminals. If False, all are printed. 


    """
    obj_ID_list = []
    flag = None
    layers = {'io_o': 101, 
             'io_e': 100, 
             'connected_o': 98,
             'connected_e': 97,
             'other': 99,
             'error_layer': 110,
             'error_width': 111,
             'error_direction': 112,
             'error_origin': 113}
    if cell.port_list:
        for p in cell.port_list:
            if p.port_type == 'o':
                if terminals and not p.terminal:
                    continue
#                if depth is not None and p.depth > depth + 1:
#                    continue
                # Track object IDs in a list
                if p.obj_id not in obj_ID_list:
                    obj_ID_list.append(p.obj_id)
                # Check connection
                if isinstance(p.prev_port, Port):
                    flag = p.is_connected(p.prev_port, verbose=False)
                elif isinstance(p.next_port, Port):
                    flag = p.is_connected(p.next_port, verbose=False)
                # Set Layer 
                if p.next_port is None or p.prev_port is None: layer=layers['io_' + p.port_type]
                elif isinstance(p.prev_port, Port) or isinstance(p.next_port, Port): layer=layers[flag]
                else: layer=layers['other']
                # Place port marker in cell
                port_marker = Polygon([(p.origin[0], p.origin[1] + p.width/2), 
                                       (p.origin[0], p.origin[1] - p.width/2), 
                                       (p.origin[0] + p.width/2, p.origin[1])], layer=layer)
                if p.next_port == -1:
                    port_marker.rotate(p.direction, center=p.origin)
                elif p.prev_port == -1:
                    port_marker.rotate(p.direction-np.pi, center=p.origin)
                cell.add(port_marker)
                # Print
                if verbose: print(p.depth)
                if verbose: print((p.obj_id))
            elif p.port_type == 'e':
                # Track object IDs in a list
                if p.obj_id not in obj_ID_list:
                    obj_ID_list.append(p.obj_id)
                # Set Circle Width
                W = 2*p.width if 2*p.width < 2.5 else 2.5
                # Set Layer 
                if p.next_port is None or p.prev_port is None: layer=layers['io_' + p.port_type]
                elif isinstance(p.prev_port, Port) or isinstance(p.next_port, Port): layer=layers['connected_' + p.port_type]
                else: layer=layers['other']
                # Place port marker in cell
                if p.terminal or p.next_port is None or p.prev_port is None:
                    cell.add(Round(p.origin, W, inner_radius=W-0.1, layer=layer))
                else:
                    port_marker = Polygon([(p.origin[0], p.origin[1] + p.width/2), 
                                           (p.origin[0], p.origin[1] - p.width/2), 
                                           (p.origin[0] + p.width/2, p.origin[1])], layer=layer)
                    if p.next_port == -1:
                        port_marker.rotate(p.direction, center=p.origin)
                    elif p.prev_port == -1:
                        port_marker.rotate(p.direction-np.pi, center=p.origin)
                    cell.add(port_marker)
                # Print
                if verbose: print(p.depth)
                if verbose: print((p.obj_id))
    if verbose: print(obj_ID_list)

def traverse_portlist(port_list, start_port=None, port_traversal_list=[], depth=None):
    """
    - This should be recursive to allow for branching (multiple input/output ports)
    - Also can be recurvsive for number of levels down to go check on ports. 
    - Make this sensitive to hierarchy depth so we skip ones below a certain depth
    """
    if not start_port:
        a = [p for p in port_list if p.prev_port is None]
        start_port = a[0]
    p = start_port
    while True:
        if p not in port_traversal_list:
            port_traversal_list.append(p)
            # can insert here code from port_check above placing circles, etc...
            print(port_traversal_list[-1])
        if isinstance(p.next_port, Port):
            p = p.next_port
        else:
            # Check for other ports of same object
            other_ports = [pt for pt in port_list if pt.obj_id == p.obj_id and pt not in port_traversal_list]
            if other_ports:
                for op in other_ports:
                    traverse_portlist(port_list, start_port=op, port_traversal_list=port_traversal_list)
            else:
#                TODO: somehow this recursion is going extra times. I artificially fixed it by only 
#                appending the list if the port isn't in it yet, but it still loops thru, so needs to be 
#                debugged at some point.
                remaining_ports = [p for p in port_list if p not in port_traversal_list]
                return port_traversal_list, remaining_ports
                break


class Port:
    """
    Port object contains origin, direction, width, name, and input/output flag 
    for the inputs and outputs of a device. is_connected() method can be used
    to determine if two ports are successfully connected.

    Parameters
    ----------
    origin : number
        Origin of the Port.
    direction : {'+x', '-x', '+y', '-y'} or number
        Direction or angle (in *radians*) the path points to. None if port 
        does not require a direction.
    width : number
        width of path at port connection.
    layer : integer
        The GDSII layer number for the input path.
    name : string
        Name of the Port.
    port_type : {'o', 'opt', 'optical', 'e', 'elec', 'electrical', None}
        String indicating type of path.
    terminal : bool
        Flags port as output terminal for outermost cell hierarchy.
    prev_port : Port, or -1, or None
        Previous connected port. If None, port has no previous connections
        (either an input port or a floating port). If -1, port has a previous 
        port connection internal to the object.
    next_port : Port, or -1, or None
        Next connected port. If None, port has no next connections (either
        an output port or a floating port). If -1, port has a next port 
        connection internal to the object.
        
    Attributes
    ----------
    obj_id : list of integers
        List of unique object idendifiers attached to the port, in hierarchical
        order from outer to inner. Each object ID in the list is a random 10-digit 
        integer from function getID(). Object IDs are automatically added to this 
        obj_id list in the Cell.add() function.
    depth : positive integer
        Hierarchy depth of the port. Default as 1, deeper in hierarchy incrementally
        adds when port is part of an object/cell within a larger object/cell.
    terminal : bool
        If True, port is a terminal for an object.
    prev_port : Port, or -1, or None
        Previous connected port. If None, port has no previous connections
        (either an input port or a floating port). If -1, port has a previous 
        port connection internal to the object.
    next_port : Port, or -1, or None
        Next connected port. If None, port has no next connections (either
        an output port or a floating port). If -1, port has a next port 
        connection internal to the object.
    """
    def __init__(self,
                 origin,
                 direction,
                 width,
                 layer,
                 name,
                 port_type='o',
                 terminal=False,
                 prev_port=None,
                 next_port=None,
                 obj_id=[]):
        self.origin = origin
        self.direction = direction
        self.width = width
        self.layer = layer
        self.name = name
        if port_type in ['o', 'opt', 'optical']:
            self.port_type = 'o'
        elif port_type in ['e', 'elec', 'electrical']:
            self.port_type = 'e'
        elif port_type is None:
            self.port_type == port_type
        else:
            self.port_type = None
            raise Exception('port_type unsupported!')
        
        self._epsilon = .001
        
        self.direction = direction_str_to_num(direction)
        self.direction = self.direction % _twopi    # Wrap direction by 2pi
        
        self.obj_id = obj_id
        if isinstance(self.obj_id, int):
            self.obj_id = [self.obj_id]
        self.depth = len(self.obj_id)
        self.terminal = terminal
        self.prev_port = prev_port
        self.next_port = next_port
        
    def __str__(self):
        message = 'Port ' + str(self.name) + ', obj ID: ' + str(self.obj_id) + ', origin: ' + \
                str((round(self.origin[0], 3), round(self.origin[1], 3))) + \
                ', dir: ' + str(round(self.direction/np.pi)) + '*pi, width: ' + str(self.width) + \
                ', layer: ' + str(self.layer)
        if isinstance(self.prev_port, Port):
            message += ', prev connected to ' + str(self.prev_port.name)
            message += ', next connected to ' + str(self.next_port)
        elif isinstance(self.next_port, Port):
            message += ', next connected to ' + str(self.next_port.name)
            message += ', prev connected to ' + str(self.prev_port)
        else:
            message += ', prev connected to ' + str(self.prev_port)
            message += ', next connected to ' + str(self.next_port)
        return message
    
    def connect_to(self, port):
        if self.port_type == port.port_type:
            self.prev_port = port
            port.next_port = self
#        TODO: How can a port have more than 1 "next_port" ??
        # TODO: Check if all port attributes are connected
    
    def add_objectID(self, obj):
        if obj.obj_id not in self.obj_id:
            self.obj_id.insert(0, obj.obj_id)
        self.depth = len(self.obj_id)
    
    def reverse_direction(self):
        self.next_port, self.prev_port = self.prev_port, self.next_port
        self.direction -= np.pi
        
#    TODO: lists of ports local to cells? port connections between cells/paths is different than inside a path
#        (they're not necessarily geometrically linked but they are logically). Ports inside a cell could be 
#        absracted away to a logical connection between input and output such that inside the cell all the 
#        internals are checked, but once they're checked then the ports are treated logically as just inputs 
#        and outputs. How do I do that? Could flag the ones that have next or prev ports as None, those would
#        be the inputs and outputs, so those get carried with the cell, or all get carried with the cell but
#        those get put in a separate category. How do we construct a port checker that traverses all the connections?
#        Have a global list (or dict or set) of ports with each port having a next and prev saved in there.
#        Then have a routine that puts them all in a list and checks for orphaned ones that don't show up
#        in the list? But branches can't have just a 1D list, so it'd have to be something slightly different. 
#        How can a port have more than 1 "next_port"?
#        
#        There could be a hierarchical thing where inside a cell you have a checker that has access to top level
#        ports of the cells it contains, then inside each of those cells they have their own internal checker, 
#        but each doesn't have access to deeper levels.

    def port_check(self, port_list):
        pass
#        TODO: create this method

    def is_connected(self, port, verbose=True):
        """
        Check if this port is connected to another.
        """
        flag = 'connected_' + port.port_type
        if not isinstance(port,Port):
            raise TypeError('port argument not a Port!')
            
#        # Port Position Check
        if (abs(self.origin[0] - port.origin[0]) > self._epsilon) or \
           (abs(self.origin[1] - port.origin[1]) > self._epsilon):
            if verbose:
                print('Warning: Position Mismatch! port ' + str(self.name) + ' at ' + str(self.origin) + 
                      ' with port ' + str(port.name) + ' at ' + str(port.origin))
            flag = 'error_origin'
        
        # Port Direction Check
        if abs(self.direction%_twopi - port.direction%_twopi) > self._epsilon or \
           abs(self.direction%_twopi - port.direction%_twopi) - np.pi > self._epsilon:
            if verbose:
                print('Warning: Direction Mismatch! port ' + str(self.name) + ' with angle ' + str(round(self.direction, 3)) + 
                      ' with port ' + str(port.name) + ' with angle ' + str(round(port.direction, 3)))
            flag = 'error_direction'

        # Port Width Check
        if not self.width == port.width:
            if verbose:
                print('Warning: Width Mismatch! port ' + str(self.name) + ' of width ' + str(self.width) + 
                      ' with port ' + str(port.name) + ' of width ' + str(port.width))
            flag = 'error_width'

        # Port Layer Check
        if not self.layer == port.layer:
            if verbose:
                print('Warning: Layer Mismatch! port ' + str(self.name) + ' in layer ' + str(self.layer) + 
                      ' with port ' + str(port.name) + ' in layer ' + str(port.layer))
            flag = 'error_layer'
        
        return flag

    def translate(self, dx, dy):
        """
        Move port by dx, dy.
        """
        vec = np.array((dx, dy))
        po = np.array(self.origin)
        self.origin = po + vec
        self.origin = (self.origin[0], self.origin[1])

    def rotate(self, angle, center=(0,0)):
        """
        Rotate port.
        """
        ca = np.cos(angle)
        sa = np.sin(angle)
        sa = np.array((-sa, sa))
        c0 = np.array(center)
        po = np.array(self.origin)
        self.origin = (po - c0) * ca + (po - c0)[::-1] * sa + c0
        self.origin = (self.origin[0], self.origin[1])
        self.direction += angle


#class Component(gdspy.PolygonSet):
#    """
#    This component class is an extension of PolygonSet with added Port handling 
#    capabilities. It also includes get_cell() which returns a singleton cell to
#    be used multiple times (checks if cell already exists and doesn't make duplicates).
#    
#    Parameters
#    ----------
#
#        
#    Attributes
#    ----------
#
#    """
#    def __init__(self):
#        pass
#    
#    def get_cell(self, name):
#        """
#        Returns a cell containing element with origin at (0,0). Checks if cell name
#        already exists and if so returns that cell, otherwise creates new cell
#        by translating element back to origin (0,0) and placing it in new cell.
#
#        Parameters
#        ----------
#        element: Component
#            object to move, containing port_list
#        name : string
#            The name of the cell.
#            
#        Returns
#        -------
#        out : ``Cell``
#            Cell containing element with origin at (0,0).
#        """
#        if name in gdspy.current_library.cell_dict:
#            self.cell = gdspy.current_library.cell_dict.get(name)
#        else:
#            self.cell = Cell(name)
#            self.translate(-self.origin[0], -self.origin[1])
#            self.cell.add(self, port_list=self.port_list)
#        return self.cell
#    
##    The problem with the cell is a "singleton" cell needs to have origin at (0,0)
##    and the CellReference will place it in position. But the PolygonSet must 
##    have origin at wherever it will be going, so first making the PolygonSet
##    and then putting it into a cell doesn't work unless I translate it back to (0,0)
##    which might be cumbersome. Also where do the ports get redefined?? For Cell
##    you would need CellReference to take care of the port translation. I could 
##    make a Component_CellReference class that inherits everything from CellReference
##    but adds in port definition
##
#
#    def put(self, cell, rotation=0, magnification=None, x_reflection=False):
#        """
#        Rotation in *radians*
#        """
#        ref = CellReference(self.cell, port_list=self.port_list, origin=self.origin, 
#                               rotation=rotation*180/np.pi, magnification=magnification, x_reflection=x_reflection)
#        cell.add(ref, port_list=ref.port_list)
#        
#        for p in self.port_list:
#            p.origin = (p.origin[0] + self.origin[0], p.origin[1] + self.origin[1])
#            po = np.array(p.origin)
#
#            ca = np.cos(rotation)
#            sa = np.sin(rotation)
#            sa = np.array((-sa, sa))
#            c0 = np.array(self.origin)
#            p.origin = (po - c0) * ca + (po - c0)[::-1] * sa + c0
#            p.origin = (p.origin[0], p.origin[1])
#            p.direction += rotation
#
#        return self
#
#    def rotate(self, angle, center=(0,0)):
#        """
#        Rotate the component. Update component port list.
#        """
#        super().rotate(angle, center)
#        ca = np.cos(angle)
#        sa = np.sin(angle)
#        sa = np.array((-sa, sa))
#        c0 = np.array(center)
#        for p in self.port_list:
#            po = np.array(p.origin)
#            p.origin = (po - c0) * ca + (po - c0)[::-1] * sa + c0
#            p.origin = (p.origin[0], p.origin[1])
#            p.direction += angle
##            Make sure port direction rotation isn't broken for all components. This might only work if 
##            origin is actually attached to object
#            
#    def scale(self, scalex, scaley=None, center=(0, 0)):
#        """
#        Scale the component. Update component port list.
#        """
#        super().scale(scalex, scaley, center)
#        c0 = np.array(center)
#        s = scalex if scaley is None else np.array((scalex, scaley))
#        for p in self.port_list:
#            po = np.array(p.origin)
#            p.origin = (po - c0) * s + c0
#            p.origin = (p.origin[0], p.origin[1])
#            if scaley:
#                if s[0] < 0:
#                    p.direction = np.pi - p.direction
#                if s[1] < 0:
#                    p.direction = - p.direction
#            else:
#                if s < 0:
#                    p.direction = p.direction - np.pi
#
#    def translate(self, dx, dy):
#        """
#        Move the polygons from one place to another. Update component port list.
#        """
#        super().translate(dx, dy)
#        vec = np.array((dx, dy))
#        for p in self.port_list:
#            po = np.array(p.origin)
#            p.origin = po + vec
#            p.origin = (p.origin[0], p.origin[1])
        

class Cell(gdspy.Cell):
    """
    Collection of elements, both geometric objects and references to
    other cells. Every cell has a port_list with all ports of all hierarchies. 
    Each port carries anobj_id which is a list of objects it belongs to, and 
    length of that list isdepth of hierarchy, carried in port.depth. The cell 
    has terminalswhich areflagged in the port.terminal. Ports carry a port_type 
    which is optical or electrical.

    Parameters
    ----------
    name : string
        The name of the cell.
    exclude_from_current : bool
        If ``True``, the cell will not be automatically included in the
        current library.

    Attributes
    ----------
    name : string
        The name of this cell.
    elements : list
        List of cell elements (``PolygonSet``, ``CellReference``,
        ``CellArray``).
    labels : list
        List of ``Label``.
    """

    def __init__(self, name, exclude_from_current=False):
        super().__init__(name, exclude_from_current)
        self.port_list = []
        self.object_list = []
        self.obj_id = getID()
    
    def add(self, element):
        """
        Add a new element or list of elements to this cell. Add port_list of the
        element to the Cell's port_list. Add Cell object ID to the each port
        in the element's port_list.

        Parameters
        ----------
        element : object, list
            The element or list of elements to be inserted in this cell.

        Returns
        -------
        out : ``Cell``
            This cell.

        """
        super().add(element)
        if isinstance(element, list):
            for e in element:
                if e.port_list:
                    for p in e.port_list:
                        p.add_objectID(self)
                    self.port_list.extend(e.port_list)
        else:
            if element.port_list:
                for p in element.port_list:
                    p.add_objectID(self)
                self.port_list.extend(element.port_list)

    def put(self, cell, input_port=None, origin=(0,0), rotation=0, magnification=None, x_reflection=False):
        """
        Description
        
        Parameters
        ----------
        cell : 
            .
        input_port : 
            
        origin : 
            
        rotation : 
            Rotation in *radians*
        magnification : 
            
        x_reflection : 
            

        Returns
        -------
        out : ``CellReference``
            Cell reference.

        """
        if input_port:
            origin = input_port.origin
            rotation = input_port.direction
        ref = CellReference(self, port_list=self.port_list, origin=origin, 
                               rotation=rotation, magnification=magnification, x_reflection=x_reflection)
        cell.add(ref)
        if input_port:
            ref.connect_to(input_port)
        return ref
    
    def get_terminal(self, num=0, terminal_type='o'):
        terminal_list = [p for p in self.port_list if p.terminal is True and p.port_type == terminal_type]
        if not terminal_list:
            raise Exception('No terminals available for object ID ' + self.obj_id + '. Check if set_terminals() was called')
        if num < 0:
            return terminal_list
        else:
            return terminal_list[num]
    
    def set_terminals(self):
        """
        Set all ending ports of a Cell as cell terminals by making port.terminal = True. 
        Does not include starting ports. Call this function at the end of creating a component Cell
        """
        for p in self.port_list:
            if p.port_type == 'o':
                p.terminal = False
                if p.next_port is None:
                    p.terminal = True
            elif p.port_type == 'e':
                p.terminal = False
                if p.next_port is None:
                    p.terminal = True
                if p.prev_port is None:
                    p.reverse_direction()
                    p.terminal = True

class CellReference(gdspy.CellReference):
    """
    Simple reference to an existing cell. Includes updating of port list.

    Parameters
    ----------
    ref_cell : ``Cell`` or string
        The referenced cell or its name.
    port_list : list of Ports
        List of Ports
    origin : array-like[2]
        Position where the reference is inserted.
    rotation : number
        Angle of rotation of the reference (in *radians*).
    magnification : number
        Magnification factor for the reference.
    x_reflection : bool
        If ``True`` the reference is reflected parallel to the x
        direction before being rotated.
    ignore_missing : bool
        If ``False`` a warning is issued when the referenced cell is not
        found.
    """
    
    def __init__(self,
                 ref_cell,
                 port_list=None,
                 origin=(0, 0),
                 rotation=None,
                 magnification=None,
                 x_reflection=False,
                 ignore_missing=False):
        
        # Rotation in degrees for gdsopy.CellReference()
        super().__init__(ref_cell, origin, rotation*180/np.pi, magnification, x_reflection, ignore_missing)
        if port_list:
            self.port_list = copy.deepcopy(port_list)
        else:
            self.port_list = copy.deepcopy(ref_cell.port_list)
        self.obj_id = getID()
        if not rotation:
            rotation = 0
        # Update Port
        for p in self.port_list:
            p.obj_id[0] = self.obj_id
            p.origin = (p.origin[0] + origin[0], p.origin[1] + origin[1])
            if rotation is not 0:
                ca = np.cos(rotation)
                sa = np.sin(rotation)
                sa = np.array((-sa, sa))
                c0 = np.array(origin)
                po = np.array(p.origin)
                p.origin = (po - c0) * ca + (po - c0)[::-1] * sa + c0
                p.origin= tuple(p.origin)
                p.direction += rotation

    def translate(self, dx, dy):
        super().translate(dx, dy)
        # Update Port
        if self.port_list:
            for p in self.port_list:
                p.origin = (p.origin[0] + dx, p.origin[1] + dy)

    def connect_to(self, port):
        p0 = [p for p in self.port_list if p.prev_port is None or p.next_port is None]
        if not p0:
            print('Unable to Connect!')
        else:
            p0 = p0[0]
            p0.connect_to(port)

    def get_terminal(self, num=0, terminal_type='o'):
        terminal_list = [p for p in self.port_list if p.terminal is True and p.port_type == terminal_type]
        if not terminal_list:
            raise Exception('No terminals available for object ID ' + str(self.obj_id) + '. Check if set_terminals() was called')
        if num < 0:
            return terminal_list
        else:
            return terminal_list[num]
    
#    I don't think set_terminals() function is needed here


class Path(gdspy.Path):
    """
    Series of geometric objects that form a path or a collection of
    parallel paths.

    Parameters
    ----------
    width : number
        The width of each path.
    initial_point : array-like[2]
        Starting position of the path.
    number_of_paths : positive integer
        Number of parallel paths to create simultaneously.
    distance : number
        Distance between the centers of adjacent paths.
    connection_port : Port
        Port this path should be connected to
    path_type : {'o', 'opt', 'optical', 'e', 'elec', 'electrical', None}
        String indicating type of path. This sets port_type automatically.
        
    Attributes
    ----------
    x : number
        Current position of the path in the x direction.
    y : number
        Current position of the path in the y direction.
    w : number
        *Half*-width of each path.
    n : integer
        Number of parallel paths.
    direction : {'+x', '-x', '+y', '-y'} or number
        Direction or angle (in *radians*) the path points to.
    distance : number
        Distance between the centers of adjacent paths.
    length : number
        Length of the central path axis.  If only one path is created,
        this is the real length of the path.
    port_list : list of Ports
        List of Ports
    start_port : Port
        Port at the beginning of the Path
    end_port : Port
        Port at the end of the Path
    """
    def __init__(self,
                 width,
                 initial_point=(0, 0),
                 number_of_paths=1,
                 distance=0,
                 connection_port=None,
                 path_type='optical'):
        if connection_port:
            initial_point = connection_port.origin
        super().__init__(width, initial_point, number_of_paths, distance)
        self.port_list = []
#        self.object_list = None    # not needed
        self.obj_id = getID()
        self.start_port = None
        self.end_port = None
        self.connection_port = connection_port
        self.path_type = path_type
        self.exclude_ports = False   # Turn to False for adiabatic bends to exclude merged segment ports
#        if type(self.connection_port) == Port:
#            self.port_list.append(self.connection_port)
        
###################### Maybe add this method:
#    def add_ports(self, ports_to_add, port_list, portset):
#        if not portset:
#            port_list[0] = port0
#        port_list[1] = port1
#        portset.update([port0, port1])
#        


    def segment(self,
                length,
                direction=None,
                final_width=None,
                final_distance=None,
                axis_offset=0,
                layer=0,
                datatype=0,
                verbose=True,
                exclude_ports=False):
        self.exclude_ports = exclude_ports
#        if length == 0:
#            self.exclude_ports = True
        direction = direction_str_to_num(direction)
        if direction is None:
            direction = self.direction
        port0 = Port((self.x, self.y), direction, self.w*2, layer, 'a_seg', obj_id=[self.obj_id], port_type=self.path_type)
        super().segment(length, direction, final_width, final_distance, axis_offset, layer, datatype)
        port1 = Port((self.x, self.y), direction, self.w*2, layer, 'b_seg', obj_id=[self.obj_id], port_type=self.path_type)

        if not self.port_list:
            if self.connection_port:
                port0.connect_to(self.connection_port)
        else:
            port0.connect_to(self.port_list[-1])
        port0.next_port = -1
        port1.prev_port = -1
        if not self.exclude_ports:
            self.port_list.extend([port0, port1])
            self.start_port = self.port_list[0]
            self.end_port = self.port_list[-1]
        
        self.exclude_ports = False

    def arc(self,
            radius,
            initial_angle,
            final_angle,
            number_of_points=0.1,
            max_points=199,
            final_width=None,
            final_distance=None,
            layer=0,
            datatype=0,
            verbose=True,
            exclude_ports=False):
        if exclude_ports:
            self.exclude_ports = exclude_ports
        
        if final_angle > initial_angle:
            _arc_direction = initial_angle + np.pi/2
        else:
            _arc_direction = initial_angle - np.pi/2
        port0 = Port((self.x, self.y), _arc_direction, self.w*2, layer, 'a_arc', 
                     obj_id=[self.obj_id], port_type=self.path_type)
        super().arc(radius, initial_angle, final_angle, number_of_points, max_points, final_width, final_distance, layer, datatype)
        port1 = Port((self.x, self.y), _arc_direction + (final_angle - initial_angle), 
                     self.w*2, layer, 'b_arc', obj_id=[self.obj_id], port_type=self.path_type)
        
        if not self.port_list:
            if self.connection_port:
                port0.connect_to(self.connection_port)
        else:
            port0.connect_to(self.port_list[-1])
        port0.next_port = -1
        port1.prev_port = -1
        if not self.exclude_ports:
            self.port_list.extend([port0, port1])
            self.start_port = self.port_list[0]
            self.end_port = self.port_list[-1]
            
        self.exclude_ports = False
        
    def arc_adiabatic(self,
                      R0,
                      Rf,
                      initial_angle,
                      final_angle,
                      n_sections = 100,
                      transition='linear',
                      number_of_points=200,
                      max_points=199,
                      final_width=None, 
                      final_distance=None, 
                      layer=0,
                      datatype=0,
                      verbose=True,
                      exclude_ports=False):
        """
        Add an adiabatic curved section to the path.

        Parameters
        ----------
        R0 : number
            Initial radius of curvature for the section.
        Rf : number
            Final radius of curvature for the section.
        initial_angle : number
            Initial angle of the curve (in *radians*).
        final_angle : number
            Final angle of the curve (in *radians*).
        n_sections : number
            Number of discrete sections of constant radius (merged later)
        transition : {'linear', 'sqrt', 'sq', 'tanh_in', 'tanh_out', 'erf_in', 'erf_out'}
            Parmetric function for adiabatic transition. 
        number_of_points : integer or float
            If integer: number of vertices that form the object
            (polygonal approximation).  If float: approximate curvature
            resolution.  The actual number of points is automatically
            calculated.
        max_points : integer
            if ``number_of_points > max_points``, the element will be
            fractured in smaller polygons with at most ``max_points``
            each.
        final_width : number
            If set, the paths of this segment will have their widths
            linearly changed from their current value to this one.
        final_distance : number
            If set, the distance between paths is linearly change from
            its current value to this one along this segment.
        layer : integer, list
            The GDSII layer numbers for the elements of each path.  If
            the number of layers in the list is less than the number of
            paths, the list is repeated.
        datatype : integer, list
            The GDSII datatype for the elements of each path (between 0
            and 255).  If the number of datatypes in the list is less
            than the number of paths, the list is repeated.
        verbose : boolean
            If True, port warnings are printed.

        Returns
        -------
        out : ``Path``
            This object.

        Notes
        -----
        The GDSII specification supports only a maximum of 199 vertices
        per polygon.
        
        n_sections=1 is NOT supported
        """
        delta_theta = (final_angle - initial_angle)/n_sections
        npoints = number_of_points // n_sections
        mpoints = max_points // n_sections
        if npoints < 4: npoints = 4
        if mpoints < 4: mpoints = 4

        if final_angle > initial_angle:
            _arc_direction = initial_angle + np.pi/2
        else:
            _arc_direction = initial_angle - np.pi/2
        port0 = Port((self.x, self.y), _arc_direction, self.w*2, layer, 'a_arc', 
                     obj_id=[self.obj_id], port_type=self.path_type)
        
        for i in range(n_sections):
    
            arg = 1.0 * i / (1.0*(n_sections-1))
            if transition is 'linear':	# Radius changes linearly with angle
                transition_fn = 1.0 * arg
            elif transition is 'sqrt':
                transition_fn = 1.0 * (arg)**(.25)
            elif transition is 'sq':
                transition_fn = 1.0 * (arg)**4.0
            elif transition is 'tanh_in':
                transition_fn = 1.0 * np.tanh(np.pi * arg)
            elif transition is 'tanh_out':
                transition_fn = 1.0 * np.tanh(np.pi * arg - np.pi)+ 1.0
            elif transition is 'erf_in':
                transition_fn = 1.0 * np.erf(2.0 * arg)
            elif transition is 'erf_out':
                transition_fn = 1.0 * np.erf(2.0 * arg - 2.0)+ 1.0
            else:
                transition_fn = 0
                print('CAREFUL: "transition" is UNDEFINED')
    
            R = R0 + transition_fn * (Rf - R0)
            # Don't include intermediate ports for tiny segments in port_list
            self.exclude_ports = True
            self.turn(R, delta_theta, number_of_points=npoints, max_points=mpoints, final_width=final_width, 
                      final_distance=final_distance, layer=layer, datatype=datatype, exclude_ports=self.exclude_ports)
    
        port1 = Port((self.x, self.y), _arc_direction + (final_angle - initial_angle), 
                     self.w*2, layer, 'b_arc', obj_id=[self.obj_id], port_type=self.path_type)
        self.exclude_ports = False
        
        if not self.port_list:
            if self.connection_port:
                port0.connect_to(self.connection_port)
        else:
            port0.connect_to(self.port_list[-1])
        port0.next_port = -1
        port1.prev_port = -1
        if not self.exclude_ports:
            self.port_list.extend([port0, port1])
            self.start_port = self.port_list[0]
            self.end_port = self.port_list[-1]
        
        self.merge(layer, n_sections=n_sections)
        self.exclude_ports = False
        
    def turn(self,
             radius,
             angle,
             number_of_points=0.1,
             max_points=199,
             final_width=None,
             final_distance=None,
             layer=0,
             datatype=0,
             exclude_ports=False):
        self.exclude_ports = exclude_ports
        super().turn(radius, angle, number_of_points, max_points, final_width, final_distance, layer, datatype)

    def turn_adiabatic(self,
                      R0,
                      Rf,
                      angle,
                      n_sections = 100,
                      transition='linear',
                      number_of_points=200,
                      max_points=199,
                      final_width=None, 
                      final_distance=None, 
                      layer=0,
                      datatype=0,
                      verbose=True,
                      exclude_ports=False):
        """
        Add an adiabatic curved section to the path.

        Parameters
        ----------
        R0 : number
            Initial radius of curvature for the section.
        Rf : number
            Final radius of curvature for the section.
        angle : {'r', 'l', 'rr', 'll'} or number
            Angle (in *radians*) of rotation of the path.  The values
            'r' and 'l' represent 90-degree turns cw and ccw,
            respectively; the values 'rr' and 'll' represent analogous
            180-degree turns.
        n_sections : number
            Number of discrete sections of constant radius (merged later)
        transition : {'linear', 'sqrt', 'sq', 'tanh_in', 'tanh_out', 'erf_in', 'erf_out'}
            Parmetric function for adiabatic transition. 
        number_of_points : integer or float
            If integer: number of vertices that form the object
            (polygonal approximation).  If float: approximate curvature
            resolution.  The actual number of points is automatically
            calculated.
        max_points : integer
            if ``number_of_points > max_points``, the element will be
            fractured in smaller polygons with at most ``max_points``
            each.
        final_width : number
            If set, the paths of this segment will have their widths
            linearly changed from their current value to this one.
        final_distance : number
            If set, the distance between paths is linearly change from
            its current value to this one along this segment.
        layer : integer, list
            The GDSII layer numbers for the elements of each path.  If
            the number of layers in the list is less than the number of
            paths, the list is repeated.
        datatype : integer, list
            The GDSII datatype for the elements of each path (between 0
            and 255).  If the number of datatypes in the list is less
            than the number of paths, the list is repeated.
        verbose : boolean
            If True, port warnings are printed.
        Returns
        -------
        out : ``Path``
            This object.

        Notes
        -----
        The GDSII specification supports only a maximum of 199 vertices
        per polygon.
        
        n_sections=1 is NOT supported
        """
        exact = True
        if angle == 'r':
            delta_i = np.pi/2
            delta_f = 0
        elif angle == 'rr':
            delta_i = np.pi/2
            delta_f = -delta_i
        elif angle == 'l':
            delta_i = -np.pi/2
            delta_f = 0
        elif angle == 'll':
            delta_i = -np.pi/2
            delta_f = -delta_i
        elif angle < 0:
            exact = False
            delta_i = np.pi/2
            delta_f = delta_i + angle
        else:
            exact = False
            delta_i = -np.pi/2
            delta_f = delta_i + angle
        self.direction = direction_str_to_num(self.direction)
        if exact:
            exact = False
        self.arc_adiabatic(R0, Rf, self.direction + delta_i, self.direction + delta_f, 
                           n_sections, transition, number_of_points, max_points, 
                           final_width, final_distance, layer, datatype, verbose)
    #    if exact:
    #        self.direction = ['+x', '+y', '-x', '-y'][int(round(path.direction / (0.5 * np.pi))) % 4]

        self.exclude_ports = False
        
    def parametric(self,
                   curve_function,
                   curve_derivative=None,
                   number_of_evaluations=99,
                   max_points=199,
                   final_width=None,
                   final_distance=None,
                   layer=0,
                   datatype=0,
                   verbose=True,
                   exclude_ports=False):
        self.exclude_ports = exclude_ports
        port0 = Port((self.x, self.y), self.direction, self.w*2, layer, 'a', obj_id=[self.obj_id], port_type=self.path_type)
#        TODO: Not sure how to set the input direction of parametric - currently this won't catch an error properly
        super().parametric(curve_function, curve_derivative, number_of_evaluations, max_points, final_width, final_distance, layer, datatype)
        port1 = Port((self.x, self.y), self.direction, self.w*2, layer, 'b', obj_id=[self.obj_id], port_type=self.path_type)

        if not self.port_list:
            if self.connection_port:
                port0.connect_to(self.connection_port)
        else:
            port0.connect_to(self.port_list[-1])
        port0.next_port = -1
        port1.prev_port = -1
        if not self.exclude_ports:
            self.port_list.extend([port0, port1])
            self.start_port = self.port_list[0]
            self.end_port = self.port_list[-1]

        self.exclude_ports = False
        
    def rotate(self, angle, center=(0, 0)):
        super().rotate(angle, center)
        ca = np.cos(angle)
        sa = np.sin(angle)
        sa = np.array((-sa, sa))
        c0 = np.array(center)
        for p in self.port_list:
            po = np.array(p.origin)
            p.origin = (po - c0) * ca + (po - c0)[::-1] * sa + c0
            p.origin = tuple(p.origin)
            p.direction += angle

    def merge(self, layer, n_sections=None, datatype=0):
        # Merge polygons in the path
        if not self.polygons:
            raise Exception('No paths present to merge!')
        if not n_sections:
            n_sections = len(self.polygons)
        # Initialize values
        new_right, new_left = [], []
        for p in self.polygons[-n_sections:]:
            [right, left] = split_list(p)
            last_right = right[-1]
            new_right = new_right + right[:-1]	# Leave out last one to avoid redundancy
            first_left = left[0]
            new_left = left[1:] + new_left	# Leave out first one to avoid redundancy
    
        new_poly = new_right+[last_right] + [first_left]+new_left
        new_poly = np.array(new_poly)	# Need to cast it to ndarray
        del self.polygons[-n_sections:]
        del self.layers[-n_sections:]
        del self.datatypes[-n_sections:]

        self.polygons.append(new_poly)
        self.layers.append(layer)
        self.datatypes.append(datatype)

class PolyPath(gdspy.PolyPath):
    """
    Series of geometric objects that form a polygonal path or a
    collection of parallel polygonal paths.

    Parameters
    ----------
    points : array-like[N][2]
        Endpoints of each path segment.
    width : number or array-like[N]
        Width of the path.  If an array is given, width at each
        endpoint.
    number_of_paths : positive integer
        Number of parallel paths to create simultaneously.
    distance : number or array-like[N]
        Distance between the centers of adjacent paths.  If an array is
        given, distance at each endpoint.
    corners : positive integer
        Type of joins: 0 - miter join, 1 - bevel join
    ends : positive integer
        Type of path ends: 0 - no extension, 1 - rounded, 2 - extended
        by half width
    max_points : integer
        The paths will be fractured in polygons with at most
        ``max_points``.
    layer : integer, list
        The GDSII layer numbers for the elements of each path.  If the
        number of layers in the list is less than the number of paths,
        the list is repeated.
    datatype : integer, list
        The GDSII datatype for the elements of each path (between 0 and
        255).  If the number of datatypes in the list is less than the
        number of paths, the list is repeated.
    connection_port : Port
        Port this path should be connected to

    Returns
    -------
    out : ``PolyPath``
        This object.

    Notes
    -----
    The bevel join will give  strange results if the number of paths is
    greater than 1.    
    """
    def __init__(self,
                 points,
                 width,
                 number_of_paths=1,
                 distance=0,
                 corners=0,
                 ends=0,
                 max_points=199,
                 layer=0,
                 datatype=0,
                 connection_port=None):
        super().__init__(points, width, number_of_paths, distance, corners, ends, max_points, layer, datatype)
        self.object_list = None
#        TODO: Port handling for PolyPath
#        self.start_port = None
#        self.end_port = None
#        self.port_list = [self.start_port, self.end_port]
#        if type(connection_port) == Port:
#            self.port_list.append(connection_port)


class Round(gdspy.Round):
    """
    """

    def __init__(self,
                 center,
                 radius,
                 inner_radius=0,
                 initial_angle=0,
                 final_angle=0,
                 number_of_points=0.01,
                 max_points=199,
                 layer=0,
                 datatype=0,
                 **kwargs):
        super().__init__(center, radius, inner_radius, initial_angle, final_angle, number_of_points, max_points, layer, datatype)
        self.port_list = []
        self.object_list = None

class Rectangle(gdspy.Rectangle):
    """
    """

    def __init__(self,
                 point1, 
                 point2, 
                 layer=0, 
                 datatype=0,
                 **kwargs):
        super().__init__(point1, point2, layer, datatype)
        self.port_list = []
        self.object_list = None

class Text(gdspy.Text):
    """
    """

    def __init__(self,
                 text,
                 size,
                 position=(0, 0),
                 horizontal=True,
                 angle=0,
                 layer=0,
                 datatype=0,
                 **kwargs):
        super().__init__(text, size, position, horizontal, angle, layer, datatype)
        self.port_list = []
        self.object_list = None

class Polygon(gdspy.Polygon):
    """
    """

    def __init__(self,
                 points, 
                 layer=0, 
                 datatype=0, 
                 verbose=True,
                 **kwargs):
        super().__init__(points, layer, datatype, verbose)
        self.port_list = []
        self.object_list = []
        self.obj_id = [getID()]

    def add_port_list(self, port_list, **kwargs):
        self.port_list.extend(port_list)
    

    
###############################################################################
if __name__ == '__main__':
    try:

        Port_test = True
        Path_test = False
        

        # Port class tests
        if Port_test:
            Port_test_cell = Cell('PORT_TEST_CELL')
            path1 = Path(0.45)
            path1.segment(10, '+y', layer=1)
            path1.segment(10, '+y', layer=1)
            path1.w = 0.5
            path1.segment(10, '+y', layer=1)
            path1.segment(10, '+x', layer=1)
            path1.segment(10, '+x', layer=2)
            path1.y += 1
            path1.segment(10, '+x', layer=2)
            
            Port_test_cell.object_list.append(path1)

            
            for obj in Port_test_cell.object_list:
                Port_test_cell.add(obj)
            
            port_check(Port_test_cell)

        # Path class tests
        if Path_test:
            Path_test_cell = Cell('PATH_TEST_CELL')
            path1 = Path(0.45)
            path1.segment(10, '+y', layer=1)
            path1.turn(10, np.pi/2, layer=1)
            path1.arc(10, np.pi, np.pi/2, layer=1)
            Path_test_cell.object_list.append(path1)
            print(path1.obj_id)
            path2 = Path(0.45, initial_point=path1.end_port.origin, connection_port=path1.end_port)
            path2.segment(10, '+x', layer=2)
            Path_test_cell.object_list.append(path2)
            print(path2.obj_id)
            path4 = Path(0.45, initial_point=path1.start_port.origin, connection_port=path2.end_port)
            path4.arc_adiabatic(50, 5, -np.pi/2, -np.pi/4, transition='tanh_in', layer=4)
            path4.arc_adiabatic(5, 50, -np.pi/4, 0, transition='tanh_out', layer=4)
            path4.turn_adiabatic(50, 5, -np.pi/2, transition='tanh_in', layer=4)
            path4.segment(10, layer=4)
            Path_test_cell.object_list.append(path4)
            print(path4.obj_id)

            for obj in Path_test_cell.object_list:
                Path_test_cell.add(obj)
            
            port_check(Path_test_cell)
#            traversed_ports, remaining_obj = traverse_objectlist(Path_test_cell.object_list)
#            print(traversed_ports)
#            print(remaining_obj)
            
    

        # Export
        print('Exporting...')
        filename = __file__.split('.')[0] + '.gds'
        gdspy.write_gds(filename, unit=1.0e-6, precision=1.0e-9)
        print("Output gds to ", filename)
    
        # Clear cells in library
        gdspy.current_library.cell_dict.clear()

    finally:
        # Clear cells in library
        gdspy.current_library.cell_dict.clear()
