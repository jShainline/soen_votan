from __future__ import division, print_function, absolute_import
import numpy as np
import phidl, phidl.geometry as pg, phidl.routing as pr
from phidl import Device, Layer, LayerSet, Port, quickplot2 as qp
from phidl import make_device
import pprint
import gdspy
import itertools
from itertools import groupby
import warnings
from collections import OrderedDict
from functools import lru_cache
import copy
#==============================================================================
# Importing things from other files
#==============================================================================

from nc_constants import lys
from nc_constants import *
#from nc_constant_devices import *

# Import the property files defined in the PDK
#try:
#    import klayout_technology as tech
#except ImportError:
#    print('Warning: klayout_technology did not load.')
#else:
#    # Then add some extra waveguides to the tech definition
#    # These waveguides need to be defined by program, not an XML, because they are complicated
#    from calc.doped_waveguides import doped_wg_dict, add_contacts
#    tech.PROPERTIES.WAVEGUIDES.update(doped_wg_dict())
#
#
#pp = pprint.PrettyPrinter(indent=0, width = 1)
#np.set_printoptions(precision=3)
#NEXT_LID=1
#PHYSLABEL_UID = 0
#import pdb
#==============================================================================
#
# Utility functions
#
#==============================================================================

def make_die(D, physlabel = 'my favorite die name 12/31/2050', maintain_hierarchy=False):
    ''' Takes care of floorplan, labeling, and float rounding error healing.
        Maintaining the hierarchy does not flatten or fix nanometer gaps.
        It reduces build time, write time, and file size. It is NOT recommended for final designs.
    '''
    D.center = (0,0)
    layer_list = [lys['wg_deep'],lys['m5_wiring'], lys['m2_nw']]
    die = D << pg.basic_die(
              size = DIE_SIZE,
              street_width = CLEAVE_STREET_WIDTH,
              street_length = CLEAVE_STREET_LENGTH,
              die_name = physlabel,
              text_size = DIE_LABEL_SIZE,
              text_location = 'S',
              layer = layer_list,
              draw_bbox = True,
              bbox_layer = lys['FLOORPLAN'],
              )
    die.center = (0,0)

    if not maintain_hierarchy:
        # flatten
        D.flatten()
        # eliminate nm gap errors
        combine_floating_points(D)
    return D

def add_label(D, label):
    if type(label) is dict:
       for p in D.get_ports(depth = None):
            p.info.update(label)
    else:
        for p in D.get_ports(depth = None):
            p.info.update({'lid':label})

def get_WB_list(D):
    wbs = {}
    for r in D.references:
        wbs = {**wbs, **get_WB_list(r.parent)}
        if 'found_wirebond' in wbs:
            return {D.ports[1]:[D]}
    for device in wbs:
        wbs[device] += [D]
    return wbs

def set_WB_info(D,depth = None):
    device_list = []
    new_depth = depth
    if depth is not None:
        new_depth -= 1
    if depth is None or depth > 0:
        for r in D.references:
            if 'is_wirebond' in r.parent.info and r.parent.info['is_wirebond'] is True:
                   r.parent.info.update(pull_info(D))
                   device_list +=[D]
            device_list+=set_WB_info(r.parent,new_depth)
    return device_list

def get_devices(D,depth = None):
    device_list = []
    new_depth = depth
    if depth is not None:
        new_depth-=1
    if depth is None or depth > 0:
        for r in D.references:
            device_list += get_devices(r.parent,new_depth)
    return (device_list + [D])

def pull_info(D):
    for r in D.references:
        D.info.update(pull_info(r.parent))
    return D.info


def add_auto_physlabel(D, offset = (0,0), position = 'N',
                       spacer = 5,layer = None):
    global NEXT_LID
    sp = spacer

    if layer is None:
        ValueError('add_auto_physlabel() layer must be defined')

    physlabel=''
    physlabel += '%0.2i (%s)' % (NEXT_LID, D.info['identifier'])
    pl = nc_text(text = physlabel, size = PHYSLABEL_SIZE, layer = layer)

    if position is 'N':
        pl.ymin, pl.x =    D.ymax + sp, D.x
    elif position is 'NE':
        pl.ymin, pl.xmax = D.ymax + sp, D.xmax
    elif position is 'EN':
        pl.ymax, pl.xmin = D.ymax     , D.xmax + sp

    elif position is 'E':
        pl.y, pl.xmin =    D.y        , D.xmax + sp
    elif position is 'ES':
        pl.ymin, pl.xmin = D.ymin     , D.xmax + sp
    elif position is 'SE':
        pl.ymax, pl.xmax = D.ymin - sp, D.xmax

    elif position is 'S':
         pl.ymax, pl.x =   D.ymin - sp, D.x
    elif position is 'SW':
        pl.ymax, pl.xmin = D.ymin - sp, D.xmin
    elif position is 'WS':
        pl.ymin, pl.xmax = D.ymin     , D.xmin - sp

    elif position is 'W':
        pl.y, pl.xmax =    D.y, D.xmin - sp
    elif position is 'WN':
        pl.ymax, pl.xmax = D.ymax     , D.xmin - sp
    elif position is 'NW':
        pl.ymin, pl.xmin = D.ymax + sp, D.xmin


    pl.move(list(offset))
    pl = D << pl

    add_label(D,NEXT_LID)
    NEXT_LID += 1
    return pl


def promote_device_to_cascade_ports(D, name_map):
    '''
        Device port names can be anything.
        To string them together in a circuit,
        one approach is to use a naming convention.

        Cascade is very simple circuit where every stage has
        one input, one output, and some number of external/biasing contacts.
        Let's also say the input is on the left for now.

        name_map is keyed by device port names and valued by special circuit port names
            * "in" input to the stage
            * "out" output to the stage
            * "gnd*" will be connected to a common ground
            * "ext_*" will be double-promoted, so it will be visible outside of the cascade circuit
            * anything else will not be promoted at all
    '''
    # Rough checks for adherence
    def quickTally(prefix):
        filt = filter(lambda ck: ck.startswith(prefix), name_map.values())
        return len(list(filt))
    assert quickTally('in') == 1
    assert quickTally('out') == 1
    assert quickTally('gnd') <= 1
    promoted_ports = dict()
    for device_pname, port in D.ports.items():
        try:
            circuit_pname = name_map[device_pname]
            assert any(circuit_pname.startswith(prefix) for prefix in ['in', 'out', 'gnd', 'ext_'])
        except (KeyError, AssertionError):
            # do not add a promoted one
            pass
        else:
            port.name = circuit_pname
            promoted_ports[circuit_pname] = port
    D.ports = promoted_ports


def parameter_combinations(parameters_dict):
    """ Creates parameter combinations from a dict filled with list values, e.g.
    parameters_dict = {
        'width' : [0.1,0.2],
        'length' : [4,5,6],
        'height' : [22]
        }
    Will produce a list of dictionaries, each of which can be used as kwargs input:
        [{'height': 22, 'length': 4, 'width': 0.1},
         {'height': 22, 'length': 5, 'width': 0.1},
         {'height': 22, 'length': 6, 'width': 0.1},
         {'height': 22, 'length': 4, 'width': 0.2},
         {'height': 22, 'length': 5, 'width': 0.2},
         {'height': 22, 'length': 6, 'width': 0.2}]
"""
    value_combinations = list(itertools.product(*parameters_dict.values()))
    keys = list(parameters_dict.keys())
    return [{keys[n]:values[n] for n in range(len(keys))} for values in value_combinations]


def gen_param_variations(
    function,
    param_defaults = None,
    param_override = None,
    param_variations = {'channel_width' : [1,5,6,7]},
    permute_parameters = False,
    device_rotation = 0,
    label_layer = 255,
    multicore = False, verbose = True,
    ):
    """ Takes e.g.

    param_variations = {
            'channel_width' : [1,2,3]
            'gate_width' : [0.1,0.2,0.4],
            }
       or equivalently
    param_variations = dict(
            channel_width = [1,2,3],
            gate_width = [0.1,0.2,0.4],
            )
    """
    if param_defaults is None:
        param_defaults = dict()
    if param_override is None:
        param_override = dict()
    if permute_parameters is True:
        parameter_list = parameter_combinations(param_variations)
    else:
        parameter_list = [dict(zip(param_variations,t)) for t in zip(*param_variations.values())]
    N_devices = len(parameter_list)

    if not multicore:
        D_list = []
        for n, params in enumerate(parameter_list):
            D_new = make_device(function, config = param_defaults, **params, **param_override)
            label_text = ''
            for name, value in params.items():
                if type(value) is float:
                    label_text += ('%s=%0.2f' % (name,value)) + '\n'
                else:
                    label_text += ('%s=%s' % (name,value)) + '\n'
            D_new.label(text = label_text, position = D_new.center, layer = label_layer)
            D_new.rotate(device_rotation)

            D_list.append(D_new)
    else:
        from nc_parallel import parallel_call
        common_kwargs = param_defaults.copy()
        common_kwargs.update(param_override)
        D_list = parallel_call(function, common_kwargs=common_kwargs, kwarg_list=parameter_list, verbose=verbose)
    return D_list


def arrange_array(
    D_list,
    x_spacing = 50,
    y_spacing = 50,
    align = 'ymax xmax',
    max_row_width = None,
    max_row_devices = None,
    ):
    ''' Update: spacings can be list/arrays to indicate non-rectilinear movements.
        This is backwards compatible with scalar arguments.

        Todo: rename x_spacing to i_vector or something like that because it is no longer purely in x
    '''
    D = Device('array')

    if np.isscalar(x_spacing):
        x_spacing = [x_spacing, 0]
    if np.isscalar(y_spacing):
        y_spacing = [0, y_spacing]

    # Add all devices
    ref_list = []
    for D_new in D_list:
        ref = D << D_new
        ref_list.append(ref)

    ref_sizes = np.array([ref.size for ref in ref_list])
    x_sizes = ref_sizes[:,0]
    y_sizes = ref_sizes[:,1]
    x_cumulative = np.cumsum(x_sizes + x_spacing[0])
    y_cumulative = np.cumsum(y_sizes + y_spacing[1])

    x = 0
    y = 0
    col = 0
    row = 0
    row_max_yheight = 0
    for n, ref in enumerate(ref_list):
        col += 1
        if 'xmax' in align: ref.xmax = x + y_spacing[0] * row
        elif 'xmin' in align: ref.xmin = x + y_spacing[0] * row
        else: ref.x = x + y_spacing[0] * row
        if 'ymax' in align: ref.ymax = y + x_spacing[1] * col
        elif 'ymin' in align: ref.ymin = y + x_spacing[1] * col
        else: ref.y = y + x_spacing[1] * col
        if y_sizes[n] > row_max_yheight: row_max_yheight = y_sizes[n]
        x = x_cumulative[n]
        # Check to see if we should move to the next row
        if max_row_width is not None:
            start_next_row = (x >= max_row_width)
        elif max_row_devices is not None:
            start_next_row = (col >= max_row_devices)
        else:
            start_next_row = False
        # If so, reset coordinates
        if start_next_row:
            x = 0
            col = 0
            row += 1
            y -= (row_max_yheight + y_spacing[1])
            row_max_yheight = 0
            x_cumulative -= x_cumulative[n]
    return D


def autoarray(
    function,
    param_defaults = {},
    param_override = {},
    param_variations = {'channel_width' : [1,5,6,7]},
    permute_parameters = False,
    device_rotation = 0,
    x_spacing = 50,
    y_spacing = 50,
    align = 'ymax',
    max_row_width = 1000,
    label_layer = 255
    ):

    D_list = gen_param_variations(
        function = function,
        param_defaults = param_defaults,
        param_override = param_override,
        param_variations = param_variations,
        permute_parameters = permute_parameters,
        device_rotation = device_rotation,
        label_layer = label_layer)

    D = arrange_array(
        D_list = D_list,
        x_spacing = x_spacing,
        y_spacing = y_spacing,
        align = align,
        max_row_width = max_row_width,
        )
    label_text = {}
    label_text.update(param_defaults)
    label_text.update(param_override)
    D.label(text = pp.pformat(label_text), position = (D.xmin, D.ymin), layer = 255)
    return D



def autoarray_xy(
    function,
    param_defaults = {},
    param_override = {},
    param_variations_x = {'width' : [1,5,6,7]},
    param_variations_y = {'length' : [1.1,2,70]},
    # permute_parameters_x = False,
    # permute_parameters_y = False,
    device_rotation = 0,
    x_spacing = 50,
    y_spacing = 50,
    align = 'ymax',
    multicore = False,
    label_layer = 255,
    ):

    param_variations = OrderedDict()
    param_variations.update(param_variations_y)
    param_variations.update(param_variations_x)

    D_list = gen_param_variations(
        function = function,
        param_defaults = param_defaults,
        param_override = param_override,
        param_variations = param_variations,
        permute_parameters = True,
        device_rotation = device_rotation,
        multicore = multicore,
        label_layer = label_layer)

    num_xvals = 1
    for val in param_variations_x.values():
        num_xvals *= len(val)
    D = arrange_array(
        D_list = D_list,
        x_spacing = x_spacing,
        y_spacing = y_spacing,
        align = align,
        max_row_devices = num_xvals,
        )

    label_text = {}
    label_text.update(param_defaults)
    label_text.update(param_override)
    D.label(text = pp.pformat(label_text), position = (D.xmin, D.ymin), layer = 255)

    return D



def packer(
        D_list,
        box_size = (9500,9500),
        spacing = 100,
        ):
    """ Takes in a list of Devices (D_list), packs them as tightly as possible
    with a minimum of spacing between them, and returns a list of Die devices
    """
    import rectpack

    bins = [box_size]*10
    spacing = spacing

    packer = rectpack.newPacker(rotation=False,
                                pack_algo = rectpack.MaxRectsBlsf,
                                sort_algo = rectpack.SORT_AREA,
                                bin_algo = rectpack.PackingBin.Global)
    # Add the rectangles to packing queue
    for n, D in enumerate(D_list):
        w,h = D.size + spacing
        packer.add_rect(w,h, rid = D)
    # Add the bins where the rectangles will be placed
    for b in bins:
        packer.add_bin(*b)
    # Start packing
    packer.pack()

    nbins = len(packer)
    Die_list = []
    for n in range(nbins):
        Die_list.append(Device())

    all_rects = packer.rect_list()
    for rect in all_rects:
        b, x, y, w, h, D = rect
        xcenter = x + w/2 + spacing/2
        ycenter = y + h/2 + spacing/2
        d = Die_list[b].add_ref(D)
        d.center = (xcenter,ycenter)

    return Die_list


def crammer(
        D_list,
        spacing = 10,
        box_size = (1,1)
        ):

    # Convert Devices to rectangles
    rect_list = []
    total_area = 0
    for n, D in enumerate(D_list):
        w,h = D.size + spacing
        total_area += w*h
        rid = D
        rect_list.append({'width':w,'height':h,'rid':rid})

    box_area = box_size[0]*box_size[1]
    box_size_start = np.array(box_size)/box_area*total_area
    for f in np.linspace(100, 1, 100):
        fit_all_bins, all_rects = _rect_pack(rect_list, box_size = f*box_size_start)
        if fit_all_bins:
            best_all_rects = all_rects
#            print(f*box_size_start)

    Crammed_D = Device()
    for rect in best_all_rects:
        b, x, y, w, h, D = rect
        xcenter = x + w/2 + spacing/2
        ycenter = y + h/2 + spacing/2
        d = Crammed_D.add_ref(D)
        d.center = (xcenter,ycenter)

    return Crammed_D


def _rect_pack(rect_list, box_size = (100,100)):
    import rectpack
    bins = [box_size]

    packer = rectpack.newPacker(rotation=False,
                                pack_algo = rectpack.MaxRectsBlsf,
                                sort_algo = rectpack.SORT_AREA,
                                bin_algo = rectpack.PackingBin.Global)

    for rect in rect_list:
        packer.add_rect(**rect)
    # Add the bins where the rectangles will be placed
    for b in bins:
        packer.add_bin(*b)
    # Start packing
    packer.pack()

    try:
        packer[0]
        fit_all_bins = (len(packer[0]) == len(rect_list))
    except:
        fit_all_bins = False
    all_rects = packer.rect_list()

    return fit_all_bins, all_rects





def _create_floating_point_merge_map(data, tol = 1e-6):
    """ Creates a dictionary which maps nearby floating points to a single
    point.  So if I give it
    x_data = [-2,-1,0,1.00001,1.0002,1.0003,4,5, 5.003, 6,7,8]
    create_floating_point_merge_map(data = x_data, tol = 1e-3)
    it returns:
    {1.00001: 1.0002,
     1.0002:  1.0002,
     1.0003:  1.0002} """
    data = np.unique(data)
    data = np.sort(data)
    indices = np.diff(data) < tol

    data_map_list = []
    n = 0
    groups = [list(j) for i, j in itertools.groupby(indices)]
    for g in groups:
        if g[0] == True:
            data_map_list += [data[n:n+len(g)+1]]
        n += len(g)

    data_map_dict = {}
    for dm in data_map_list:
        target = np.median(dm)
        for d in dm:
            data_map_dict[d] = target
    return data_map_dict

def combine_floating_points(D, tol = 1e-6):
    """ Takes a flattened Device, and merges all nearby floating point values
    together to eliminate nanometer-size errors"""
    polygons = D.get_polygons(by_spec = False)
    all_points = np.vstack(polygons)

    x_correction = _create_floating_point_merge_map(all_points[:,0], tol = tol)
    y_correction = _create_floating_point_merge_map(all_points[:,1], tol = tol)

    for poly in D.polygons:
        if hasattr(poly, 'points'):
            polygons = [poly.points]
        else:
            polygons = poly.polygons
        for polygon in polygons:
            for point in polygon:
                if point[0] in x_correction:
                    point[0] = x_correction[point[0]]
                if point[1] in y_correction:
                    point[1] = y_correction[point[1]]

def drc_exclude(D):
        D << pg.bbox(D.bbox, layer = lys['DRC_exclude'])
        return D


# visual debugging
try:
    from lyipc.client import kqp
except:
    pass

#==============================================================================
#
# Helper geometry functions
#    - These are functions which may be reused in more than one test device
#    - For instance, a function which creates a resistance-measurement geometry
#      may be used on multiple layers
#    - Layers should be called out explicitly -- pass individual layers with
#      meaningful names as arguments (so they can be reused), e.g.
#        def waveguide(width = 5, deep_layer_etch = 2, shallow_layer_etch = 4)
#    - Set default layer #s to the 240-250 range so they're easy to spot
#
#==============================================================================

def nc_text(text = 'abcd', size = 10, position=(0, 0), justify = 'left', layer = 255):
    # generates text with the DRC exclude layer around it
    D = Device()

    t = D << pg.text(text = text, size = size, justify = justify, layer = layer)
    t.center = position

    drc_exclude = D << pg.compass(size = t.size, layer = lys['DRC_exclude'])
    drc_exclude.center = t.center

    return D

def wire_basic(p1, p2, directions = 'xxyx', width = 5, layer = 240):
    D = Device(name = 'wire_basic')

    if isinstance(p1, phidl.device_layout.Port): p1 = p1.midpoint
    if isinstance(p2, phidl.device_layout.Port): p2 = p2.midpoint

    directions = directions.lower()
    num_x = sum([xy == 'x' for xy in directions])
    num_y = sum([xy == 'y' for xy in directions])
    distance = np.asarray(p2)-p1

    points = [p1]
    for xy in directions:
        if xy == 'x':
            travel = np.array([distance[0]/num_x, 0])
            new_point = points[-1] + travel
        elif xy == 'y':
            travel = np.array([0, distance[1]/num_y])
            new_point = points[-1] + travel
        else:
            raise ValueError('[PHIDL] wire_basic() directions argument must be string with only "x" or "y" characters')
        if not np.all(np.abs(travel) < width/1e3): # Only add point if traveling some significant distance
            points.append(new_point)

    route_path = gdspy.PolyPath(points, width = width)
    D.add_polygon(route_path.polygons, layer = layer)
    return D


def resistor(
    Pad = None,
    res_num_squares = 2,
    res_wire_width = RES_WIDTH,
    ):

    D = Device()

    pad_size = Pad.size
    res_layer = lys['m3_res']
    pad_layer = lys['m5_wiring']

    res_length = res_num_squares*res_wire_width
    Rp = pg.rectangle(size = pad_size, layer = [Layer(res_layer), Layer(pad_layer)])
    Rr = pg.rectangle(size = [res_length, res_wire_width], layer = res_layer)
    rp1 = D << Rp
    rp2 = D << Rp
    rr = D << Rr
    rp2.xmin = rp1.xmax + res_length
    rr.center = rp1.center
    rr.xmin = rp1.xmax

    D.add_port(name = 1, midpoint = rp1.center, width = pad_size[1], orientation = -180)
    D.add_port(name = 2, midpoint = rp2.center, width = pad_size[1], orientation = 0)
    D.flatten()
    return D


def padify(port, pad_size = WIREBOND_SIZE, pad_offset = [150,25],
           layer = 0, inset_layer_distance = 5, inset_layer = None,
           inset_include_taper = False,
           ):
    """ Creates a contact pad with an optional inset layer, offset
    a certain distance from a port.  Useful for creating pad structures.
    """

    D = Device('padify')
    midp = port.midpoint
    Pad = pg.rectangle(size = pad_size, layer = layer)
    Pad.add_port(name = 1, midpoint = [Pad.xmin,Pad.y],
                 width = pad_size[1], orientation = 180)
    pad = D.add_ref(Pad)
    pad.center = midp
    pad.move(pad_offset)
    pad.rotate(port.orientation, center = midp)
#    if override_port_width is not None:
#        new_port = port._copy()
#        new_port.width = override_port_width
    eps1 = port.endpoints
    eps2 = pad.ports[1].endpoints
    pts = [eps1[0], eps1[1], eps2[0], eps2[1]]
    taper = D.add_polygon(pts, layer = layer)

    if inset_layer is not None:
        if inset_include_taper: inset_items = [pad, taper]
        else:                   inset_items = [pad]
        D.add_ref(pg.offset(inset_items, distance = -inset_layer_distance, layer = inset_layer))
        D.add_ref(pg.copy_layer(pad, layer = layer, new_layer = inset_layer))

    D.add_port(name = 1, midpoint = pad.center, width = pad_size[1], orientation = port.orientation)
    D.flatten()
    return D

def port2pad(port, Pad = None, pad_offset = [150,25],
    position_override = None, wiring_layer = 0, pad_info = {'test':True}):
    D = Device('port2pad')
    MPad = pg.extract(Pad, layers = metal_layers)
    pad = D.add_ref(Pad)
    wiring_pad = D << pg.compass(MPad.size, layer = wiring_layer)
    wiring_pad.center = pad.center
    midp = port.midpoint
    center = pad.center
    orientation = port.orientation
    pad.move(center, midp + pad_offset)
    pad.rotate(orientation, center = midp)
    wiring_pad.move(center, midp + pad_offset)
    wiring_pad.rotate(orientation, center = midp)
    port_pts = port.endpoints
    pad_pts = wiring_pad.ports['W'].endpoints

    if position_override is not None:
        pad.center = position_override
        wiring_pad.center = position_override
        pad_pts = [wiring_pad.ports['E'].endpoints[1], wiring_pad.ports['E'].endpoints[0]]

        #pad_pts_shift = np.array([np.sin(port.orientation)*

    pts = [pad_pts[0], pad_pts[1], port_pts[0], port_pts[1]]
    D.add_polygon(pts, layer = wiring_layer)
    new_pad_info= {**pad_info,**{'has_info':True}}
    D.add_port(pad.ports[1])
    D.ports[1].info.update(new_pad_info)

    return D

def pad2pad(port, Pad = None, directions = 'xxyx', width = 5,
            pad_offset = [150,25], wiring_layer = 0):
    D = Device('pad2pad')
    pad = D.add_ref(Pad)
    midp = port.midpoint
    center = pad.center
    orientation = port.orientation
    pad.move(center, midp + pad_offset)
    pad.rotate(orientation, center = midp)
    D << wire_basic(midp, pad.center, directions = directions,
        width = width, layer = wiring_layer)
    D.add_port(pad.ports[1])
    return D


# def inductorify(port = None, pad_size = VIA_SIZE, wire_width = 0.1, num_squares = 1000,
#                 inductor_layer = 241, pad_layer = 242):
#     D = Device('inductor')

#     wire_pitch = wire_width/MEANDER_FILL_FACTOR
#     area_per_sq = (wire_pitch)*wire_width
#     side_length = np.sqrt(area_per_sq*num_squares)

#     L = snspd_with_pads(
#         wire_width = wire_width,
#         size = [side_length,side_length],
#         num_squares = None,
#         pad_size = pad_size,
#         snspd_layer = inductor_layer,
#         pad_layer = pad_layer,
#         gnd_layer = pad_layer,
#         )

#     l = D.add_ref(L)
#     if port is not None:
#         l.connect(1, port)
#     D.add_port(l.ports[1])
#     D.add_port(l.ports[2])

#     D.flatten()

#     return D


def inductor_with_pads(
    Pad = None,
    wire_width = 0.1,
    num_squares = 1000,
    ):
    return snspd_with_pads(Pad = Pad,
            num_squares = num_squares,
            wire_width = wire_width,
            )



def ytron_with_pads(
        Pad = None,
        rho = 0.05,
        arm_widths = [1,1],
        ytron_layer = 240
        ):

    D = Device('ytron_with_pads')

    Ytron = pg.ytron_round(
            rho = rho,
            arm_widths = arm_widths,
            arm_lengths = 2*np.array(arm_widths),
            source_length = sum(arm_widths),
            theta = 5,
            theta_resolution = 10,
            layer = ytron_layer)

    y = D << Ytron

    # Build optimal steps
    step_ratio = 2
    step_l = D <<  pg.optimal_step(start_width = arm_widths[0],
                                 end_width = arm_widths[0]*step_ratio,
                                 num_pts = 25, width_tol = 1e-3,
                                 anticrowding_factor = 1.2,
                                 layer = ytron_layer)
    step_r = D <<  pg.optimal_step(start_width = arm_widths[1],
                                 end_width = arm_widths[1]*step_ratio,
                                 num_pts = 25, width_tol = 1e-3,
                                 anticrowding_factor = 1.2,
                                 layer = ytron_layer).reflect()
    step_l.connect(1, y.ports['left'])
    step_r.connect(1, y.ports['right'])

    pad_width = max(Pad.size)
    pad_l = D << port2pad(step_l.ports[2], Pad, pad_offset = [pad_width*1.5,pad_width*0.75], wiring_layer = ytron_layer)
    pad_r = D << port2pad(step_r.ports[2], Pad, pad_offset = [pad_width*1.5,-pad_width*0.75], wiring_layer = ytron_layer)
    pad_s = D << port2pad(y.ports['source'], Pad, pad_offset = [pad_width*1.5,0], wiring_layer = ytron_layer)

    D.add_port(name = 'wiring1', port = pad_l.ports[1])
    D.add_port(name = 'wiring2', port = pad_r.ports[1])
    D.add_port(name = 'wiring3', port = pad_s.ports[1])

    D.info = y.info

    return D


def via_test_defaults(Pad = None, pad_layer = 240,
    via_width = 100,
    wiring1_layer = 241, wiring2_layer = 242, via_layer = 243):
    pad_size = Pad.size
    D = pg.test_via(num_vias = 40, wire_width = WIRING_WIDTH, via_width = via_width,
        via_spacing = max(WIRING_WIDTH, min(VIA_SIZE))*4, pad_size = pad_size, min_pad_spacing = 0,
        pad_layer = pad_layer, wiring1_layer = wiring1_layer,
        wiring2_layer =wiring2_layer, via_layer = via_layer)
    [[xmin,ymin], [xmax, ymax]] = D.bbox
    pad1 = D.add_ref(Pad)
    pad2 = D.add_ref(Pad)
    pad1.xmin, pad1.ymin = xmin, ymin
    pad2.xmax, pad2.ymax = xmax, ymax

    return D


def comb_isolation_defaults(Pad = None, comb_layer = 240, overlap_zigzag_layer = 241,
                            comb_pad_layer = 242, overlap_pad_layer = 243):
    pad_size = Pad.size
    D = pg.test_comb(pad_size = pad_size, wire_width = 4, wire_gap = MIN_WIRING_GAP,
              comb_layer = comb_layer, overlap_zigzag_layer = overlap_zigzag_layer,
              comb_pad_layer = comb_pad_layer, comb_gnd_layer = comb_pad_layer,
              overlap_pad_layer = overlap_pad_layer)
    [[xmin,ymin], [xmax, ymax]] = D.bbox
    x = (xmax+xmin)/2
    y =(ymax+ymin)/2
    pad1 = D.add_ref(Pad)
    pad2 = D.add_ref(Pad)
    pad3 = D.add_ref(Pad)
    pad4 = D.add_ref(Pad)
    pad1.xmin, pad1.y = xmin, y
    pad2.xmax, pad2.y = xmax, y
    pad3.x, pad3.ymax = x, ymax
    pad4.x, pad4.ymin = x, ymin
    return D


def sheet_res_defaults(Pad = None, res_layer = 240, pad_layer = 241,
                            gnd_layer = 242):
    pad_size = Pad.size
    D = pg.test_res(pad_size = pad_size,
                     num_squares = 1000,
                     width = 5,
                     res_layer = res_layer,
                     pad_layer = pad_layer,
                     gnd_layer = gnd_layer)
    [[xmin,ymin], [xmax, ymax]] = D.bbox
    pad1 = D.add_ref(Pad)
    pad2 = D.add_ref(Pad)
    pad1.xmin, pad1.ymin = xmin, ymin
    pad2.xmax, pad2.ymax = xmax, ymax
    return D


def ntron(channel_width = 1, gate_width = 0.1, gate_length = 0.05,
          width_terminals = 2, layer = 240):
    D = Device('ntron')
    S = pg.optimal_step(start_width = channel_width, end_width = width_terminals, num_pts = 50, width_tol = 1e-3,
                     anticrowding_factor = 1.2, layer = layer)
    s_bot = D.add_ref(S).rotate(-90)
    s_top = D.add_ref(S).reflect()
    s_top.connect(1, s_bot.ports[1])

    g = D.add_polygon(gdspy.PolyPath(points = [[0,0],
                                               [width_terminals,0],
                                               [width_terminals*2,0],
                                               [width_terminals*2+gate_length,0],
                                               ],
                                     width = [width_terminals,width_terminals,gate_width,gate_width]), layer = layer)
    g = g[0]
    g.y = s_bot.ymax
    g.xmax = s_bot.xmin

    D.add_port(name = 'drain', port = s_top.ports[2])
    D.add_port(name = 'source', port = s_bot.ports[2])
    D.add_port(name = 'gate', midpoint = (g.xmin, g.y), orientation = -180, width = width_terminals)

    D.info['gate_width'] = gate_width
    D.info['channel_width'] = channel_width
    return D


def ntron_with_pads(
    Pad = None,
    channel_width = 1,
    gate_width = 0.1,
    gate_length = 0.05,
    width_terminals = 2,
    ):

    ntron_layer = lys['m2_nw']
#    parameters_dict = locals() # save the parameters input into the function
    Ntron = ntron(channel_width = channel_width, gate_width = gate_width, gate_length = gate_length,
                      width_terminals = width_terminals, layer = ntron_layer)
#    Pad = pg.compass(pad_size, layer = ntron_layer)

    D = Device('ntron_with_pads')
    nt = D.add_ref(Ntron)

#    Pads = Device('ntron_pads')
    pad_width = max(Pad.size)
    pad_d = D << port2pad(nt.ports['drain'], Pad, pad_offset = [pad_width*1.5,0], wiring_layer = ntron_layer)
    pad_g = D << port2pad(nt.ports['gate'], Pad, pad_offset = [pad_width*1.5,0], wiring_layer = ntron_layer)
    pad_s = D << port2pad(nt.ports['source'], Pad, pad_offset = [pad_width*1.5,0], wiring_layer = ntron_layer)

#    pad_g = D.add_ref(port2pad(nt.ports['gate'], pad_size = pad_size, pad_offset = [max(pad_size), 0],  layer = ntron_layer,
#                     inset_layer = pad_layer, inset_layer_distance = PAD_OVERHANG))
#    pad_s = D.add_ref(port2pad(nt.ports['source'], pad_size = pad_size, pad_offset = [max(pad_size), 0],  layer = ntron_layer,
#                     inset_layer = gnd_layer, inset_layer_distance = PAD_OVERHANG))

    D.add_port(name = 'drain', port = pad_d.ports[1])
    D.add_port(name = 'gate', port = pad_g.ports[1])
    D.add_port(name = 'source', port = pad_s.ports[1])

    D.info = nt.info
    # D.info.update(parameters_dict)

    return D


def ntron_with_inductor(
    Pad = None,
    channel_width = 1,
    gate_width = 0.1,
    gate_length = 0.05,
    width_terminals = 2,
    L_ntron_num_squares=1000,
    ):
    D = Device()

    Ntron = ntron_with_pads(
            Pad = Pad,
            channel_width = channel_width,
            gate_width = gate_width,
            gate_length = gate_length,
            width_terminals = width_terminals,
            )

    L = inductor_with_pads(Pad = Pad,
                    wire_width = channel_width*INDUCTOR_WIDTH_RATIO,
                    num_squares = L_ntron_num_squares,
                    )

    nt = D.add_ref(Ntron)
    l = D.add_ref(L)
    l.move(l.ports[2], nt.ports['drain'])
#
    D.add_port(port = l.ports[1], name = 'wiring1')
    D.add_port(port = nt.ports['gate'], name = 'wiring2')
    D.add_port(port = nt.ports['drain'], name = 'wiring3')
    D.add_port(port = nt.ports['source'], name = 'wiring4')
    return D


def ktron(
    nanowire_width =  3,
    nanowire_height =  None,
    heater_width = 2,
    heater_height = 2,
    heater_num_squares = None,
    heater_offset = 0, # Offset of heater from left edge of nanowire
    nanowire_layer = 240,
    resistor_layer = 241,
    resistor_pad_layer = 242,
    ):

    # Create blank device
    D = Device(name = 'ktron')


    if nanowire_height is None:
        nanowire_height = nanowire_width

    if heater_num_squares is not None:
        heater_height = heater_width*heater_num_squares
    else:
        heater_num_squares = heater_height/heater_width

    # Create components
    res_wiring_layer = [Layer(resistor_layer), Layer(resistor_pad_layer)]
    heater_size = [heater_width, heater_height]
    nanowire_size = [nanowire_width, nanowire_height]
    Nanowire = pg.compass(size = nanowire_size, layer = nanowire_layer)
    Heater = pg.compass(size = heater_size, layer = resistor_layer)
    Wire = pg.compass(size = [heater_width, heater_width], layer = res_wiring_layer)

    # Add references to components
    nw = D.add_ref(Nanowire)
    h = D.add_ref(Heater)
    w1 = D.add_ref(Wire)
    w2 = D.add_ref(Wire)
    h.xmin = nw.xmin + heater_offset
    h.y = nw.y
    w1.x = w2.x = h.x
    w1.ymin = h.ymax
    w2.ymax = h.ymin

    ## Record meta-information
    D.info['nanowire_width'] = nanowire_width
    D.info['heater_size'] = np.round(heater_size,2).tolist()
    D.info['heater_area'] = np.round(heater_size[0]*heater_size[1],2)
    D.info['heater_num_squares'] = np.round(heater_num_squares,2)

    D.add_port(name = 1, port = w1.ports['W'])
    D.add_port(name = 2, port = w2.ports['W'])
    D.add_port(name = 3, port = nw.ports['N'])
    D.add_port(name = 4, port = nw.ports['S'])

    D.flatten()

    return D

def ktron_with_pads(
    Pad = None,
    nanowire_width =  3,
    nanowire_height =  None,
    heater_width = 1,
    heater_height = 2,
    heater_num_squares = None,
    heater_offset = 0,
    pad_offset = 50,
    ):

    # Create blank device
    D = Device('ktron_with_pads')

    nanowire_layer = lys['m2_nw']
    resistor_layer = lys['m3_res']
    resistor_pad_layer = lys['m5_wiring']

    # Create components
    H = ktron(
    nanowire_width =  nanowire_width,
    nanowire_height =  nanowire_height,
    heater_width = heater_width,
    heater_height = heater_height,
    heater_num_squares = heater_num_squares,
    heater_offset = heater_offset,
    nanowire_layer = nanowire_layer,
    resistor_layer = resistor_layer,
    resistor_pad_layer = resistor_pad_layer,
    )


    # Add references to components meander_pad_layer
    h = D.add_ref(H)

    pad_width = max(Pad.size)
    res_wiring_layer = [Layer(resistor_layer), Layer(resistor_pad_layer)]
    pad_n = D << port2pad(h.ports[1], Pad, pad_offset = [pad_width*2,-pad_width*0.75], wiring_layer = res_wiring_layer)
    pad_s = D << port2pad(h.ports[2], Pad, pad_offset = [pad_width*2,pad_width*0.75], wiring_layer = res_wiring_layer)
    pad_w = D << port2pad(h.ports[3], Pad, pad_offset = [pad_offset,0], wiring_layer = nanowire_layer)
    pad_e = D << port2pad(h.ports[4], Pad, pad_offset = [pad_offset,0], wiring_layer = nanowire_layer)


    D.add_port(name = 1, port = pad_n.ports[1])
    D.add_port(name = 2, port = pad_s.ports[1])
    D.add_port(name = 3, port = pad_w.ports[1])
    D.add_port(name = 4, port = pad_e.ports[1])

    # Record information
    D.info = H.info
    # D.info.update(parameters_dict)


    return D



def ktron_with_inductor(
    Pad = None,
    nanowire_width =  3,
    nanowire_height =  None,
    heater_width = 1,
    heater_height = 2,
    heater_num_squares = None,
    heater_offset = 0,
    pad_offset = 50,
    L_num_squares=1000,
    L_width_ratio = INDUCTOR_WIDTH_RATIO,
    ):
    D = Device('ktron_with_L')

    Ktron = ktron_with_pads(
        Pad = Pad,
        nanowire_width =  nanowire_width,
        nanowire_height =  nanowire_height,
        heater_width = heater_width,
        heater_num_squares = heater_num_squares,
        heater_offset = heater_offset,
        pad_offset = pad_offset,
            )

    L = inductor_with_pads(Pad = Pad,
                    wire_width = nanowire_width*L_width_ratio,
                    num_squares = L_num_squares,
                    )

    kt = D.add_ref(Ktron)
    l = D.add_ref(L)
    l.move(l.ports[2], kt.ports[3])
#
    D.add_port(port = l.ports[1], name = 'wiring3')
    D.add_port(port = kt.ports[2], name = 'wiring2')
    D.add_port(port = kt.ports[1], name = 'wiring1')
    D.add_port(port = kt.ports[4], name = 'wiring4')
    return D


def htron(
    nanowire_width =  0.15,
    nanowire_spacing = 0.1,
    meander_num_squares = 1000,
    heater_num_squares = 1,
    meander_layer = 240,
    resistor_layer = 241):
    # Create blank device
    D = Device(name = 'htron')

    # Basic calculations
    extra_meander_width = 2
    extra_meander_height = 1
    area_per_meander_sq = (nanowire_width+nanowire_spacing)*nanowire_width
    meander_area = area_per_meander_sq*meander_num_squares
    meander_total_width = np.sqrt(meander_area/heater_num_squares)
    meander_total_height = heater_num_squares*meander_total_width
    meander_size = np.array([meander_total_width, meander_total_height])
    heater_size = meander_size

    meander_size = meander_size + [extra_meander_width,extra_meander_height]

    # meander_size = heater_size + np.array([meander_extra_width,0])
    meander_pitch = nanowire_width + nanowire_spacing
    # heater_standoff_y = 1

    # Create components
    Meander = pg.snspd_expanded(wire_width = nanowire_width, wire_pitch = meander_pitch, size = meander_size,
               terminals_same_side = False, connector_width = nanowire_width*4, layer = meander_layer)
    # heater_size_actual = heater_size + np.array([0, heater_standoff_y])
    Heater = pg.compass(size = heater_size, layer = resistor_layer)

    # Add references to components
    m = D.add_ref(Meander)
    h = D.add_ref(Heater)
    h.center = m.center


    # Record meta-information
    heater_area = heater_size[0]*heater_size[1]
    D.info['nanowire_width'] = nanowire_width
    D.info['nanowire_pitch'] = nanowire_width + nanowire_spacing
    D.info['meander_num_squares'] = np.round(m.info['num_squares'],2)
    D.info['meander_size'] = np.round((m.xsize, m.ysize),2).tolist()
    D.info['heater_size'] = np.round(heater_size,2).tolist()
    D.info['heater_area'] = np.round(heater_size[0]*heater_size[1],2)
    D.info['heater_num_squares'] = np.round(heater_num_squares,2)
    D.info['overlap_area'] = np.round(m.ysize*heater_size[0],1)
    D.info['overlap_num_squares'] = np.round(heater_area/area_per_meander_sq,1)

    D.add_port(name = 1, port = h.ports['N'])
    D.add_port(name = 2, port = h.ports['S'])
    D.add_port(name = 3, port = m.ports[1])
    D.add_port(name = 4, port = m.ports[2])

    # D.label(pp.pformat(D.info), D.center) # "Pretty-print" the info dict
    D.flatten()
    return D

#    qp(D)
def htron_with_pads(
    Pad = None,
    nanowire_width =  0.15,
    nanowire_spacing = 0.1,
    meander_num_squares = 1000,
    heater_num_squares = 1,
    meander_layer = 240,
    resistor_layer = 242,
    resistor_pad_layer = 243,
    ):

    # parameters_dict = locals() # save the parameters input into the function

    # Create blank device
    D = Device('htron_with_pads')

    # Create components
    H = htron(
    nanowire_width =  nanowire_width,
    nanowire_spacing = nanowire_spacing,
    meander_num_squares = meander_num_squares,
    heater_num_squares = heater_num_squares,
    meander_layer = meander_layer,
    resistor_layer = resistor_layer,
    )


    # Add references to components meander_pad_layer
    h = D.add_ref(H)

    pad_width = max(Pad.size)
    res_wiring_layer = [Layer(resistor_layer), Layer(resistor_pad_layer)]
    pad_n = D << port2pad(h.ports[1], Pad, pad_offset = [pad_width*1.5,0], wiring_layer = res_wiring_layer)
    pad_s = D << port2pad(h.ports[2], Pad, pad_offset = [pad_width*1.5,0], wiring_layer = res_wiring_layer)
    pad_w = D << port2pad(h.ports[3], Pad, pad_offset = [pad_width*1.5,0], wiring_layer = meander_layer)
    pad_e = D << port2pad(h.ports[4], Pad, pad_offset = [pad_width*1.5,0], wiring_layer = meander_layer)


    D.add_port(name = 1, port = pad_n.ports[1])
    D.add_port(name = 2, port = pad_s.ports[1])
    D.add_port(name = 3, port = pad_w.ports[1])
    D.add_port(name = 4, port = pad_e.ports[1])

    # Record information
    D.info = H.info
    # D.info.update(parameters_dict)


    return D

def snspd_with_pads(
        Pad = None,
        wire_width = 0.15,
        num_squares = 6000,
        ):
#    parameters_dict = locals() # save the parameters input into the function
    D = Device('snspd_with_pads')

    snspd_layer = lys['m2_nw']

    wire_pitch = wire_width/MEANDER_FILL_FACTOR
    area_per_sq = (wire_pitch)*wire_width
    side_length = np.sqrt(area_per_sq*num_squares)

    s = D.add_ref(pg.snspd_expanded(wire_width = wire_width, wire_pitch = wire_pitch,
                              size = [side_length,side_length], num_squares = None, terminals_same_side = False,
                              connector_width = wire_width*3, layer = snspd_layer)).rotate(-90)
    pad_width = max(Pad.size)
    pad1 = D << port2pad(s.ports[1], Pad, pad_offset = [pad_width*1.5,0], wiring_layer = snspd_layer)
    pad2 = D << port2pad(s.ports[2], Pad, pad_offset = [pad_width*1.5,0], wiring_layer = snspd_layer)

    D.add_port(name = 1, port = pad1.ports[1])
    D.add_port(name = 2, port = pad2.ports[1])

    D.info = s.info

    D.info['snspd_layer'] = snspd_layer
    D.info['wire_width'] = wire_width
    D.info['num_squares'] = num_squares

    D.info['expected_resistance'] = num_squares*EXPECTED_RSQ_WSI
    D.info['length'] = num_squares/wire_width

    return D


def snspd_with_inductor(
        Pad = None,
        Snspd = None, # Takes snspd_with_pads()
        L_snspd_num_squares = 1000,
    ):
    D = Device('snspd_with_inductor')

    snspd_layer = Snspd.info['snspd_layer']
    snspd_wire_width = Snspd.info['wire_width']
    inductor_width = max(snspd_wire_width*INDUCTOR_WIDTH_RATIO, INDUCTOR_MIN_WIDTH)
    L = inductor_with_pads(Pad = Pad,
                    wire_width = inductor_width,
                    num_squares = L_snspd_num_squares,
                    )

    s = D.add_ref(Snspd)
    l = D.add_ref(L)
    l.move(l.ports[2], s.ports[1])
#
    D.add_port(port = l.ports[1], name = 'wiring1')
    D.add_port(port = s.ports[1], name = 'wiring2')
    D.add_port(port = s.ports[2], name = 'wiring3')
    return D

def ntron_with_ntron(
        Ntron_with_L_1 = None, # Takes output from ntron_with_inductor()
        Ntron_with_L_2 = None, # Takes output from ntron_with_inductor()
        Resistor = None, # Takes output from resistor()
        wiring_layer = 240,
        ):

    D = Device('ntron_with_ntron')

    n1 = D << Ntron_with_L_1
    n2 = D << Ntron_with_L_2
    r = D << Resistor

    # Align the nTrons vertically
    pad_width = n1.ports['wiring2'].width
    n1.ymin = n2.ymin
    r.move(r.ports[1], n1.ports['wiring3']).movex(pad_width*2)
#
    # Place the ntron
    n2.xmin = max(r.xmax, n1.xmax) + pad_width

    # Wire together
    D << wire_basic(p1 = n1.ports['wiring3'], p2 = r.ports[1], directions = 'xyx', width = WIRING_WIDTH, layer = wiring_layer)
    D << wire_basic(p1 = n2.ports['wiring2'], p2 = r.ports[2], directions = 'xyx', width = WIRING_WIDTH, layer = wiring_layer)
    D << wire_basic(p1 = n1.ports['wiring4'], p2 = n2.ports['wiring4'], directions = 'xyx', width = WIRING_WIDTH, layer = wiring_layer)

    # Add ports
    D.add_port(port = n1.ports['wiring2'], name = 'wiring1')
    D.add_port(port = n1.ports['wiring1'], name = 'wiring2')
    D.add_port(port = n2.ports['wiring1'], name = 'wiring3')
    p4 = D.add_port(port = n2.ports['wiring3'], name = 'wiring4')
    p4.orientation = 0
    D.add_port(port = n1.ports['wiring4'], name = 'wiring5')

    return D


def snspd_with_ntron(
        Snspd_with_L = None, # Takes output from snspd_with_inductor()
        Ntron_with_L = None, # Takes output from ntron_with_inductor()
        Resistor = None, # Takes output from resistor()
        ):

    D = Device('snspd_with_ntron')

    wiring_layer = lys['m5_wiring']

    r = D << Resistor
    s = D << Snspd_with_L
    n = D << Ntron_with_L

    # Align the nTron and SNSPD vertically
    pad_width = s.ports['wiring2'].width
    s.ymin = n.ymin
    r.move(r.ports[1], s.ports['wiring2']).movex(pad_width*2)

    # Place the ntron
    n.xmin = max(r.xmax, s.xmax) + pad_width

    # Wire together
    D << wire_basic(p1 = s.ports['wiring2'], p2 = r.ports[1], directions = 'xyx', width = WIRING_WIDTH, layer = wiring_layer)
    D << wire_basic(p1 = n.ports['wiring2'], p2 = r.ports[2], directions = 'xyx', width = WIRING_WIDTH, layer = wiring_layer)
    D << wire_basic(p1 = s.ports['wiring3'], p2 = n.ports['wiring4'], directions = 'xyx', width = WIRING_WIDTH, layer = wiring_layer)

    # Add ports
    D.add_port(port = s.ports['wiring1'], name = 'wiring1')
    D.add_port(port = n.ports['wiring1'], name = 'wiring2')
    p = D.add_port(port = n.ports['wiring3'], name = 'wiring3')
    p.orientation = 0
    D.add_port(port = s.ports['wiring3'], name = 'wiring4')


    return D


def ntron_with_htron(
        Ntron_with_L = None, # Takes output from ntron_with_inductor()
        Htron = None, # Takes output from htron_with_pads()
        wiring_layer = 241,
        ):

    D = Device('ntron_with_htron')


    n = D << Ntron_with_L
    h = D << Htron

    pad_width = n.ports['wiring3'].width
    h.rotate(90)
    h.move(h.ports[1], n.ports['wiring3'])
    h.xmin = n.xmax + pad_width

    # Wire together
    D << wire_basic(p1 = n.ports['wiring3'], p2 = h.ports[1], directions = 'xyx', width = WIRING_WIDTH, layer = wiring_layer)
    D << wire_basic(p1 = n.ports['wiring4'], p2 = h.ports[3], directions = 'xy', width = WIRING_WIDTH, layer = wiring_layer)
    D << wire_basic(p1 = n.ports['wiring4'], p2 = h.ports[2], directions = 'xy', width = WIRING_WIDTH, layer = wiring_layer)

    # Add ports
    D.add_port(port = n.ports['wiring2'], name  = 'wiring1')
    D.add_port(port = n.ports['wiring1'], name  = 'wiring2')
    D.add_port(port = h.ports[4], name  = 'wiring3')
    midpoint = (h.ports[3].midpoint[0], n.ports['wiring4'].midpoint[1])
    portwidth = n.ports['wiring1'].width
    D.add_port(name  = 'wiring4', width = portwidth, midpoint = midpoint, orientation = -90)
    return D


def ktron_with_htron(
        Ktron_with_L = None, # Takes output from ntron_with_inductor()
        Htron = None, # Takes output from htron_with_pads()
        ):

    D = Device('ktron_with_htron')

    wiring_layer = lys['m5_wiring']

    k = D << Ktron_with_L
    h = D << Htron

    pad_width = k.ports['wiring3'].width
    h.rotate(90)
    h.move(h.ports[1], k.ports['wiring3'])
    h.xmin = k.xmax + pad_width

    # Wire together
    D << wire_basic(p1 = k.ports['wiring3'], p2 = h.ports[1], directions = 'xyx', width = WIRING_WIDTH, layer = wiring_layer)
    D << wire_basic(p1 = k.ports['wiring4'], p2 = h.ports[3], directions = 'xy', width = WIRING_WIDTH, layer = wiring_layer)
    D << wire_basic(p1 = k.ports['wiring4'], p2 = h.ports[2], directions = 'xy', width = WIRING_WIDTH, layer = wiring_layer)

    # Add ports
    D.add_port(port = k.ports['wiring1'], name  = 'wiring1')
    D.add_port(port = k.ports['wiring2'], name  = 'wiring2')
    D.add_port(port = k.ports['wiring3'], name  = 'wiring3')
    D.add_port(port = h.ports[4], name  = 'wiring4')
    midpoint = (h.ports[3].midpoint[0], k.ports['wiring4'].midpoint[1])
    portwidth = k.ports['wiring1'].width
    D.add_port(name  = 'wiring5', width = portwidth, midpoint = midpoint, orientation = -90)
    return D


def ntron_with_ytron(
        Ntron_with_L = None, # Takes output from ntron_with_inductor()
        Ytron = None, # Takes output from ytron_with_pads()
        Resistor = None, # Takes output from resistor()
        wiring_layer = 241,
        ):

    D = Device('ntron_with_ytron')

    r = D << Resistor
    n = D << Ntron_with_L
    y = D << Ytron

    # Align the nTron and yTron vertically
    pad_width = n.ports['wiring2'].width
    r.move(r.ports[1], n.ports['wiring3']).movex(pad_width*2)

#     Place the ytron
    y.ymin = n.ymin
    y.xmin = max(r.xmax, n.xmax) + pad_width

    # Wire together
    D << wire_basic(p1 = n.ports['wiring4'], p2 = y.ports['wiring3'], directions = 'xyx', width = WIRING_WIDTH, layer = wiring_layer)
    D << wire_basic(p1 = n.ports['wiring3'], p2 = r.ports[1], directions = 'xyx', width = WIRING_WIDTH, layer = wiring_layer)
    D << wire_basic(p1 = y.ports['wiring1'], p2 = r.ports[2], directions = 'yx', width = WIRING_WIDTH, layer = wiring_layer)

    # Add ports
    D.add_port(port = n.ports['wiring2'], name = 'wiring1')
    D.add_port(port = n.ports['wiring1'], name = 'wiring2')
    D.add_port(port = y.ports['wiring2'], name = 'wiring3')
    D.add_port(port = n.ports['wiring4'], name = 'wiring4')

    return D

def ic_wire(wire_widths = [0.25, 0.5,1,2,4], wire_widths_wide = [0.75, 1.5, 3, 4, 6], Pad = None,
            pitch = 50, layerset = lys):
    """
    Usage:

    Call ic_test_structure() with either a list of widths for the thickest part of each wire to test and a list for the
    thinnest parts of each wire. Alternatively, specify a list of widths for the thinnest part of each wire and ignore the
    wire_widths parameter. Instead you should specify the width_growth_factor which indicates by what factor the thick
    part of the wire will be larger than the thin part.
    Ex:
        ic_test_structure(wire_widths = [5,10,10,10,10], thin_width=[0.5,1,2,3,4])
        - or -
        ic_test_structure(width_growth_factor = 5, thin_width=[0.5,1,2,3,4])
    """
    D = Device('ic_wire')

    pitch =  pitch
    first_wire_center = [0,0]
    sc_layer = lys['m2_nw']
    length_test = max(wire_widths_wide)*4
    n = np.size(wire_widths)

    for i, x in enumerate(wire_widths_wide):

        Wire = pg._test_ic_wire_step(wire_widths_wide[i], wire_widths[i], wire_layer=sc_layer)
        Wire.center = [0,0]
        Wire.add_port(name = 'a', midpoint = [Wire.xmin, Wire.y], width = Wire.size[1], orientation = 180)
        Wire.add_port(name = 'b', midpoint = [Wire.xmax, Wire.y], width = Wire.size[1], orientation = 0)
        Wire.add_port(name = 1, midpoint = [-1*length_test, Wire.y], width = max(wire_widths_wide), orientation = 0)
        Wire.add_port(name = 2, midpoint = [length_test, Wire.y], width = max(wire_widths_wide), orientation = 180)
        Wire.add_ref(pr.route_basic(port1 = Wire.ports[1], port2 = Wire.ports['a'], layer = sc_layer))
        Wire.add_ref(pr.route_basic(port1 = Wire.ports[2], port2 = Wire.ports['b'], layer = sc_layer))
        Wire.ports[1].orientation = 180
        Wire.ports[2].orientation = 0
        wire = D.add_ref(Wire)

        wire.rotate(90)
        wire.x = first_wire_center[0] + i*pitch
        wire.y = 0

        f = D << port2pad(port = wire.ports[2], Pad = Pad, pad_offset = [50,0], wiring_layer = sc_layer)
        g = D << port2pad(port = wire.ports[1], Pad = Pad, pad_offset = [50,0], wiring_layer = sc_layer)
        p = D.add_port(port = f.ports[1], name = i)
        p.orientation = p.orientation+180
        p = D.add_port(port = g.ports[1], name = i+np.size(wire_widths))
        p.orientation = p.orientation + 180

        #D.info['wire_width'] =
        D.info['expected_resistance'] = EXPECTED_RSQ_WSI

    return D

def straight_wire(length = 50, width = 5,  nsquares = None, Pad = None, layerset = lys):

    D = Device('straight_wire')

    if nsquares is None:
        nsquares = length/width
    elif length is None:
        length = nsquares*width
    elif width is None:
        width = nsquares/length
    else:
        print('Two of length, width, nsquares should be defined')
        return

    wire_layer = layerset['m2_nw']
    wire = D.add_ref(pg.compass(size = (length, width), layer = wire_layer))
    taper = pg.optimal_step(start_width = wire.size[1], end_width = min(wire.size[1]*4, Pad.size[0]), num_pts = 50, width_tol = 1e-3,
                 anticrowding_factor = 1.2, layer = layerset['m2_nw'])
    taper1 = D << taper
    taper2 = D << taper
#   taper2.rotate(180)

    taper1.connect(port = 1, destination = wire.ports['E'])
    taper2.connect(port = 1, destination = wire.ports['W'])

    padinfo = {'num_squares': nsquares}
    padinfo['expected_resistance'] = padinfo['num_squares']*EXPECTED_RSQ_WSI
    padinfo['wire_width'] = width
    padinfo['length'] = padinfo['num_squares']*padinfo['wire_width']
    p1 = D << port2pad(port = taper1.ports[2], Pad = Pad, pad_offset = [Pad.size[0], 0], wiring_layer = wire_layer)
    p2 = D << port2pad(port = taper2.ports[2], Pad = Pad, pad_offset = [Pad.size[0], 0], wiring_layer = wire_layer)
    D.info.update(padinfo)
    D.add_port(p1.ports[1])
    D.add_port(port = p2.ports[1], name = 2)
    return D


####
# Alex/Adam vortex detectors
####
def shielded_nwarc(channel_width=1, terminal_width=20, terminal_length=100, radius=None, theta=30, notch_depth=None, overhang=1):
    ''' Default radius is 10x width
        Notch depth None means no notch
    '''
    D = Device()
    if radius is None:
        radius = 10 * channel_width
    nw_arc = pg.arc(radius=radius, width=channel_width, theta=theta, layer=lys['m2_nw'])
    shield_arc = pg.arc(radius=radius, width=channel_width + 2*overhang, theta=theta)
    if notch_depth is not None:
        notch_theta = theta - 10
        notch_width = overhang + notch_depth + .1
        notch_radius = radius - channel_width/2 - notch_width/2 + notch_depth
        notch = pg.arc(radius=notch_radius, width=notch_width, theta=notch_theta, start_angle=(theta-notch_theta)/2)
        shield_arc = pg.boolean(shield_arc, notch, 'A-B')
    D << nw_arc
    D << shield_arc
    D.add_port('nw_1', port=nw_arc[1])
    D.add_port('nw_2', port=nw_arc[2])
    D.add_port('shield_1', port=shield_arc[1])
    D.add_port('shield_2', port=shield_arc[2])
    return D


def shielded_nwarc_with_pad(nwarc):
    ''' Also puts on an inductor.


    '''
    # put NW_PAD_DEVs on ports 'nw_1' and 'nw_2'
    # wire the one on 'nw_1' to an inductor
    # wire that inductor to a WB_PAD_DEV
    # wire the one on 'nw_2' to WB_GND_PAD_DEV
    pass
    # this function version 2 will also:
    # put big m2_nw squares around nw_pads
    # extend coverage of shield by using 'shield_1' ports


###############################################################################
#
# SAEED KHAN DEVICES
# PACKAGING
#
###############################################################################
def Rect(length=10, width = 1, layer = 0):
    D = Device('Rec')
    xpts = [-length/2, length/2, length/2, -length/2]
    ypts = [width/2, width/2, -width/2, -width/2]

    D.add_polygon([xpts,ypts], layer = layer)
    D.center=[0,0]
    return D


def Waveguide(length=10, width = 1, layer = 0):
    D = Device('waveguide')
    D=Rect(length=length,width=width,layer=layer)
    D.center=[0,0]
    return D

def Taper(length = 10, width1 = 5, width2 = None, RelPos=0, layer = lys['wg_deep']):
    if width2 is None: width2 = width1
    xpts = [-length/2, length/2, length/2, -length/2]
    ypts = [-width1/2, -width2/2+RelPos, width2/2+RelPos, width1/2]

    D = Device('taper')
    D.add_polygon([xpts,ypts], layer = layer)
    D.center=[0,0]
    return D

def Grating(length = 25, width = 10, period = 1.24, fFactor = 0.5, layer = lys['wg_shallow']):
    D = Device('grating')
    w=period*fFactor
    total=int(np.ceil(length/period))
    #remain=length-total*period
    for i in range(total):
        xpts = [-length/2+i*period, -length/2+i*period+w, -length/2+i*period+w, -length/2+i*period]
        ypts = [-width/2, -width/2, width/2, width/2]
        D.add_polygon([xpts,ypts], layer = layer)

    #xpts = [-length/2+total*period, -length/2+total*period+remain, -length/2+total*period+remain, -length/2+total*period]
    #ypts = [-width/2, -width/2, width/2, width/2]
    D.add_polygon([xpts,ypts], layer = layer)

    D.center=[0,0]

    return D




def FatLine(Points, Width=1):

    x,y=Points
    nPoints = len(x)
    xb1=np.arange(nPoints, dtype=float)
    xb2=np.arange(nPoints, dtype=float)
    yb1=np.arange(nPoints, dtype=float)
    yb2=np.arange(nPoints, dtype=float)
    theta=np.arctan((y[1]-y[0])/(x[1]-x[0]))

    xb1[0]=x[0]+Width/2*np.sin(theta)
    yb1[0]=y[0]-Width/2*np.cos(theta)
    xb2[0]=x[0]-Width/2*np.sin(theta)
    yb2[0]=y[0]+Width/2*np.cos(theta)

    for i in range(1,nPoints):
        theta=np.arctan((y[i]-y[i-1])/(x[i]-x[i-1]))

        xb1[i]=x[i]+Width/2*np.sin(theta)
        yb1[i]=y[i]-Width/2*np.cos(theta)
        xb2[i]=x[i]-Width/2*np.sin(theta)
        yb2[i]=y[i]+Width/2*np.cos(theta)

    xb2= xb2[::-1]
    yb2=yb2[::-1]
    xb=xb1
    xb=np.append(xb,xb2)

    yb=yb1
    yb=np.append(yb,yb2)




    return [xb, yb]



def ShapeS(Points, StepRes=1000):

    p=[[Points[0][0],Points[1][0]],[Points[0][1],Points[1][1]],[Points[0][2],Points[1][2]],[Points[0][3],Points[1][3]]]
    sigma=[1,3,3,1]

    UB=np.zeros((StepRes,4))
    frac=np.zeros(StepRes)
    #u=np.ones(res)
    for i in range(StepRes):
        frac[i]=float(i)/StepRes
    P=np.zeros((StepRes,2))
    for uu in range(StepRes):
        u=frac[uu]
        dd=0
        for dd in range(4):
            d=dd+1
            UB[uu][dd]=sigma[dd]*((1-u)**(4-d))*(u**(d-1))



    for i in range(len(UB)):
        for j in range(len(p[0])):
            for k in range(len(p)):
                P[i][j] += UB[i][k]  * p[k][j]


    xvals = P[:,0]
    yvals = P[:,1]

    xvals =np.append(xvals,Points[0][3])
    yvals =np.append(yvals,Points[1][3])


    return xvals, yvals


def Bowtie(MinAx=100,MajAx=200,lChannel=30,wChannel=50,l1=300, l2=600, layer=lys['su8']):

    D = Device('Bow')

    d1 = Device('d1')

    Cir=pg.ellipse((MinAx,MajAx) ,layer=layer)
    d1.add_ref(Cir)


    r=Rect(2*lChannel+2*MinAx,wChannel,layer=layer)
    d1.add_ref(r)

    theta=np.arctan(l2/2/l1)
    #print(theta*180/np.pi)
    x=wChannel/2/np.tan(theta)
    e1=Taper(length=l1, width1=0,width2=l2,layer=layer)
    e1.movex(lChannel+MinAx+l1/2-x)
    d1.add_ref(e1)

    e2=Taper(length=l1, width1=l2,width2=0,layer=layer)
    e2.movex(-lChannel-MinAx-l1/2+x)
    d1.add_ref(e2)

    d2=Rect(2*lChannel+2*MinAx+2*l1+200,l2+200, layer=layer)

    D=pg.boolean( A = d1, B = d2, operation = 'B-A',layer=layer)


    D.rotate(90)
    D.center=[0,0]

    #d1.add_ref(d2)
    return D


#: experimentally optimized
#: note the Bowtie defaults have already been changed to reflect optimum
wWG1220 = 0.35#tech.waveguides('Strip').get_by_layer('wg_deep')[0].width  # width of 1220 sm WG
wWG1550 = 0.45#tech.waveguides('Strip1550').get_by_layer('wg_deep')[0].width  # width of 1550 sm WG
khan_nominals_1220 = dict(pGrating=0.510, wWG=wWG1220, include_mode_filter=True, label='Khan-1220-v1')
khan_nominals_1550 = dict(pGrating=0.670, wWG=wWG1550, include_mode_filter=True, label='Khan-1550-v1')


def gratingKhan(pGrating=1.24, lGrating=25, wGrating=11, fFactor=0.5, lTaper=300, wWG=0.35, include_mode_filter=False, label=None):
    # calculations
    wGratingBox=wGrating-1;
    lGratingBox=lGrating*2;

    D = Device()
    # Add grating
    r1 = D << pg.compass((lGratingBox, wGratingBox), layer=lys['wg_deep'])
    r1.center = (0, 0)
    gr1 = D << Grating(length=lGrating,width=wGrating,period=pGrating,fFactor=fFactor,layer=lys['wg_shallow'])
    gr1.center = r1.center

    # Add Taper
    tp1 = D << pg.taper(lTaper, width1=wGratingBox, width2=wWG, layer=lys['wg_deep'])
    tp1.connect(1, r1.ports['E'])

    # Mode filter (copied from grating3)
    if include_mode_filter==True:
        mf = D << mode_filter(width=wWG, sm_width=wWG-.05, mf_radius=8.5*wWG)
        # width=width, sm_width=sm_width, mf_radius=mf_radius, layer=layer
        mf.connect(2, tp1.ports[2])
        D.add_port(name=1, port=mf.ports[1])
    else:
        D.add_port(name=1, port=tp1.ports[2])
    D.add_port('wg_in_1', port=D.ports[1])

    if label is not None:
        D.label(text = label, position = (D.xmin, D.ymax + D.ysize))

    return D


def bowtie_GC(**grating_kwargs):
    ''' no thin step '''
    D = Device()
    bow = D << Bowtie()
    grating = D << gratingKhan(**grating_kwargs)
    # No alignment because they to come out just right centered at 0
    # place an exit waveguide
    covered_distance = bow.xmax - grating.xmax + 15
    if covered_distance > 0:
        extension = D << pg.compass((covered_distance, grating.ports[1].width), layer=lys['wg_deep'])
        extension.connect('E', grating.ports[1])
        D.add_port(name=1, port=extension.ports['W'])
    else:
        D.add_port(name=1, port=grating.ports[1])
    return D


def SimpleDevice(lWG=2000,wWG=.35,pGrating=1.24,lGrating=25,wGrating=11,fFactor=0.5,lTaper=300,xb=9, lChannel=30,wChannel=50,l1=300, l2=600):
    # These were used in Saeed's initial development
    D = Device('Dev')
    lShape=1000
    wGratingBox=wGrating-1;
    lGratingBox=lGrating*2;
    #Add Waveguide
    rwg=9
    lwg1=rwg*lWG/16-lShape/2
    wg1=Waveguide(length=lwg1, width = wWG,layer=lys['wg_deep'])
    wg1.movex(-lwg1/2)
    D.add_ref(wg1)
    #Add S shape
    S = Device('ShapeS')
    x1pts = [0, lShape/2, lShape/2,lShape]
    y1pts = [0,0, lShape/2,lShape/2]

    x,y=ShapeS([x1pts,y1pts])
    xf, yf=FatLine([x,y],wWG)
    S.add_polygon([xf,yf], layer = lys['wg_deep'])

    #S.movex(lShape/2)
    D.add_ref(S)

     #Add Waveguide
    lwg2=(16-rwg)*lWG/16-lShape/2
    wg2=Waveguide(length=lwg2, width = wWG,layer=lys['wg_deep'])
    wg2.movex(lShape+lwg2/2)
    wg2.movey(lShape/2)

    D.add_ref(wg2)




    # Add Input Taper
    tp1=Taper(length=lTaper,width1=wGratingBox,width2=wWG,layer=lys['wg_deep'])
    tp1.movex(-lwg1-lTaper/2)
    D.add_ref(tp1)


    # Add Output Taper
    tp2=Taper(length=lTaper,width1=wWG,width2=wGratingBox,layer=lys['wg_deep'])
    tp2.movex(lwg2+lTaper/2+lShape)
    tp2.movey(lShape/2)
    D.add_ref(tp2)


    # Add Input grating

    r1=Rect(lGratingBox,wGratingBox,layer=lys['wg_deep'])
    r1.movex(-lwg1-lTaper-lGratingBox/2)
    D.add_ref(r1)

    gr1=Grating(length=lGrating,width=wGrating,period=pGrating,fFactor=fFactor,layer=lys['wg_shallow'])
    gr1.movex(-lwg1-lTaper-lGratingBox/2)
    D.add_ref(gr1)

    # Add output grating

    r2=Rect(lGratingBox,wGratingBox,layer=lys['wg_deep'])
    r2.movex(lwg2+lTaper++lShape+lGratingBox/2)
    r2.movey(lShape/2)
    D.add_ref(r2)

    gr2=Grating(length=lGrating,width=wGrating,period=pGrating,fFactor=fFactor,layer=lys['wg_shallow'])
    gr2.movex(lwg2+lTaper++lShape+lGratingBox/2)
    gr2.movey(lShape/2)
    D.add_ref(gr2)


    # Bowtie
    MajAx=200 #major radius
    MinAx=100

    bt1=Bowtie(MinAx=MinAx,MajAx=MajAx,lChannel=lChannel,wChannel=wChannel,l1=l1, l2=l2,layer=lys['su8'])
    bt1.movex(-lwg1-lTaper-lGratingBox/2)
    #bt1.movey(-(2*lChannel+2*radius1+2*l1+200)/2)
    D.add_ref(bt1)


    bt2=Bowtie(MinAx=MinAx,MajAx=MajAx,lChannel=lChannel,wChannel=wChannel,l1=l1, l2=l2,layer=lys['su8'])
    bt2.movex(lwg2+lTaper+lShape+lGratingBox/2)
    bt2.movey(lShape/2)
    D.add_ref(bt2)


    # Bowtie for Pedestal

    #from grating simulations
    #xb = 10  #coordinate where beam intersects grating for maximum efficiency
    theta = 15*np.pi/180 #angle of fiber
    rf = 62.5 #radius of fiber

    #from fab to get position right
    h1 = 52    #height of SU-8 fiber collar
    h2 = 21.5  #height of SU-8 fiber pedestal
    hOx = 1.5 #thickness of oxide cladding

    lbc = rf/np.cos(theta)+hOx*np.tan(theta)
    xc = xb+lbc
    lb = h1*np.tan(theta)
    #lblah=2*rf*np.cos(theta)

    l3 = 2*rf*np.cos(theta)+h1*np.tan(theta);
    xCenter = xc+lb-l3/2 #Center of the ellipse
    MajAx=l3/2+5 #major radius
    MinAx=rf+5

    xp = xb+rf/np.cos(theta)+hOx*np.tan(theta)-h2/np.tan(theta)

    wPedestal = 130;
    lPedestal = 1.618*wPedestal;


    bt3=Bowtie(MinAx=MinAx,MajAx=MajAx,lChannel=lChannel,wChannel=wChannel,l1=l1, l2=l2,layer=lys['su8'])
    numPeriods=int(np.ceil(lGrating/pGrating));

    bt3.movex(-lwg1-lTaper-(lGratingBox/2-(numPeriods*pGrating/2-pGrating/2)))
    bt3.movex(-xCenter)

    #bt1.movey(-(2*lChannel+2*radius1+2*l1+200)/2)
    D.add_ref(bt3)
    r3=Rect(lPedestal,wPedestal,layer=lys['su8_thin'])
    r3.movex(-lwg1-lTaper-(lGratingBox/2-(numPeriods*pGrating/2-pGrating/2+pGrating/4))+lPedestal/2)
    r3.movex(-xp)
    D.add_ref(r3)

    bt4=Bowtie(MinAx=MinAx,MajAx=MajAx,lChannel=lChannel,wChannel=wChannel,l1=l1, l2=l2,layer=lys['su8'])

    bt4.movex(lwg2+lTaper+lShape+(lGratingBox/2-(numPeriods*pGrating/2-pGrating/2)))
    bt4.movex(xCenter)
    bt4.movey(lShape/2)
    #bt1.movey(-(2*lChannel+2*radius1+2*l1+200)/2)
    D.add_ref(bt4)
    r4=Rect(lPedestal,wPedestal,layer=lys['su8_thin'])
    r4.movex(lwg2+lTaper+lShape+(lGratingBox/2-(numPeriods*pGrating/2-pGrating/2+pGrating/4))-lPedestal/2)
    r4.movex(xp)
    r4.movey(lShape/2)
    D.add_ref(r4)


    D.center=[0,0]
    return D


###############################################################################
#
# JEFF CHILES DEVICES
# PASSIVE PHOTONICS
#
###############################################################################


def grating2(num_periods=20, period=0.75, fill_factor=0.5, ar_offset=False, width_grating=5,
             length_taper=10, width=0.4,sm_width=0.3, layer=lys['wg_deep'], layer_shallow=lys['wg_shallow'],
             partial_etch=False, include_mode_filter=False, packaging = False, mf_radius=3):


   # returns a fiber grating
   G = Device('grating')
   Temp = Device()
   if partial_etch==False:
       # make the grating teeth
       for i in range(num_periods):
           cgrating = G.add_ref(pg.rectangle(size=[period * fill_factor, width_grating], layer=layer))
           cgrating.x += i * period
   else:
       for i in range(num_periods):
           cgrating = G.add_ref(pg.rectangle(size=[period * (1-fill_factor), width_grating], layer=layer_shallow))
           cgrating.x += i * period
       Block=pg.taper(width1=width_grating,width2=width_grating,length=period*num_periods+1,layer=layer)
       block=G.add_ref(Block)
       block.xmax=cgrating.xmax
       block.y=cgrating.y

   if packaging == True:
       gluebox_dev = G << gluebox(size = (GLUEBOX_SIZE,GLUEBOX_SIZE), corner_cut = GLUEBOX_CORNER, wall_width = GLUEBOX_WALL, layer = lys['su8'])
       gluebox_dev.center = cgrating.center
       gluebox_dev.movex(-50)
   # make the taper
   tgrating = G.add_ref(pg.taper(length=length_taper, width1=width_grating, width2=width, port=None, layer=layer))
   if ar_offset==True:
       tgrating.xmin = cgrating.xmax+period
   else:
       tgrating.xmin=cgrating.xmax
   tgrating.y = cgrating.y

   if include_mode_filter==True:
       mf = G << mode_filter(width=width, sm_width=sm_width, mf_radius=mf_radius, layer=layer)
       mf.connect(2, tgrating.ports[2])
       p = G.add_port(name=1, port=mf.ports[1])
   else:
       p = G.add_port(port=tgrating.ports[2], name=1)

   return G

#same as grating2 but Shainline's AR tooth instead of Chiles'

def grating3(num_periods=20, period=0.55, fill_factor=0.81, ar_tooth = True, ar_tooth_width = 0.36, width_grating=5,
             length_taper=10, width=0.4,sm_width=0.3, layer=lys['wg_deep'],
             include_mode_filter=False, packaging = False, mf_radius=3):


   # returns a fiber grating
   G = Device('grating')
   Temp = Device()
   oldx=0

   for i in range(num_periods):
       if (ar_tooth==True) & (i == num_periods-2):
          cgrating = G.add_ref(pg.rectangle(size=[ar_tooth_width, width_grating], layer=layer))
          cgrating.xmin = oldx + period*(1-fill_factor)
          oldx = cgrating.xmax
       else:
          cgrating = G.add_ref(pg.rectangle(size=[period * fill_factor, width_grating], layer=layer))
          cgrating.xmin = oldx + period*(1-fill_factor)
          oldx = cgrating.xmax

   if packaging == True:
       gluebox_dev = G << gluebox(size = (GLUEBOX_SIZE,GLUEBOX_SIZE), corner_cut = GLUEBOX_CORNER, wall_width = GLUEBOX_WALL, layer = lys['su8'])
       gluebox_dev.center = cgrating.center
       gluebox_dev.movex(-50)
   # make the taper
   tgrating = G.add_ref(pg.taper(length=length_taper, width1=width_grating, width2=width, port=None, layer=layer))

   tgrating.xmin=cgrating.xmax
   tgrating.y = cgrating.y

   if include_mode_filter==True:
       mf = G << mode_filter(width=width, sm_width=sm_width, mf_radius=mf_radius, layer=layer)
       mf.connect(2, tgrating.ports[2])
       p = G.add_port(port=mf.ports[1], name=1)
   else:
       p = G.add_port(port=tgrating.ports[2], name=1)

   return G


def mode_filter(width=0.4, sm_width=0.3, mf_radius=3, layer=lys['wg_deep']):
    G = Device()
    Taper_down=pg.taper(length=5,width1=width,width2=sm_width,layer=layer)
    taper_down=G.add_ref(Taper_down)
    G.add_port(name=2, port=taper_down.ports[1])
    # taper_down.connect(port=1,destination=tgrating.ports[2])

    MF1=pg.arc(radius=mf_radius, width=sm_width, theta=90, angle_resolution=0.5, layer=layer)
    MF2 = pg.arc(radius=mf_radius, width=sm_width, theta=-90, angle_resolution=0.5, layer=layer)
    mf1=G.add_ref(MF1)
    mf1.connect(port=1,destination=taper_down.ports[2])

    mf2=G.add_ref(MF2)
    mf2.connect(port=1,destination=mf1.ports[2])

    taper_up=G.add_ref(Taper_down)
    taper_up.connect(port=2,destination=mf2.ports[2])
    taper_up.xmin=taper_up.xmin+mf_radius*5
    taper_up.y=taper_down.y

    R=pr.route_basic(port1=mf2.ports[2],port2=taper_up.ports[2],layer=layer)
    r=G.add_ref(R)
    p = G.add_port(port=taper_up.ports[1], name=1)
    return G

#standalone component: bendy polarizer
def bendy_polarizer(bend_radius=4,num_sections=10,wg_width=SM_WG_WIDTH,length_straight=2.0, layerset=lys):
    D=Device()
    layer=layerset['wg_deep']
    A_cw_init=pg.arc(radius=bend_radius,width=wg_width,theta=90.0,start_angle=-90,angle_resolution=0.25,layer=layer)
    A_cw=pg.arc(radius=bend_radius,width=wg_width,theta=180.0,start_angle=-90,angle_resolution=0.25,layer=layer)
    A_ccw=pg.arc(radius=bend_radius,width=wg_width,theta=-180.0,start_angle=-90,angle_resolution=0.25,layer=layer)
    A_ccw_out=pg.arc(radius=bend_radius,width=wg_width,theta=-90.0,start_angle=-90,angle_resolution=0.25,layer=layer)

    B=pg.taper(length=length_straight,width1=wg_width,width2=wg_width,layer=layer)

    a1=D.add_ref(A_cw_init)

    b1=D.add_ref(B)
    b1.connect(port=1,destination=a1.ports[2])

    prev_port=b1.ports[2]
    D.add_port(name=1,port=a1.ports[1])

    #keep connecting bends together for num_sections
    for x in range(num_sections):
        a = D.add_ref(A_ccw)
        a.connect(port=1,destination=prev_port)
        b = D.add_ref(B)
        b.connect(port=1,destination=a.ports[2])

        a = D.add_ref(A_cw)
        a.connect(port=1,destination=b.ports[2])
        b = D.add_ref(B)
        b.connect(port=1,destination=a.ports[2])
        prev_port=b.ports[2]

    a2=D.add_ref(A_ccw_out)
    a2.connect(port=1,destination=prev_port)
    D.add_port(name=2,port=a2.ports[2])


    return D


#into and out of the shallow etch mode
#basically useless, just a helper for a test structure
def shallow_deep_transition(wg_width=SM_WG_WIDTH,max_width=3,taper_length=10,layerset=lys):
    deep_layer=layerset['wg_deep']
    shallow_layer=layerset['wg_shallow']
    D=Device()
    T=Device()
    F=Device()
    Taper=pg.taper(length=taper_length,width1=wg_width,width2=max_width,layer=deep_layer)
    t1=T.add_ref(Taper)
    t2=T.add_ref(Taper)
    t2.connect(port=2,destination=t1.ports[2])

    Wg=pg.taper(length=taper_length*2,width1=wg_width,width2=wg_width,layer=deep_layer)
    wg=F.add_ref(Wg)
    T.xmin=wg.xmin
    D.add_ref(T)
    Shallow_mask=pg.boolean(A=T,B=wg,operation='A-B',layer=shallow_layer)
    shallow_mask=D.add_ref(Shallow_mask)
    D.add_port(name=2,port=t2.ports[1])
    D.add_port(name=1,port=t1.ports[1])

    return D




#essentially a different type of straight routing which blows out the waveguide to a larger width for most of the length
def longhaul_waveguide(
        port1,
        port2,
        taper_length=15,
        max_wg_width=MM_WG_WIDTH,
        layerset=lys
):
    layer=layerset['wg_deep']
    D=Device()

    #test if there is room for a routed waveguide at all between the two tapers
    if taper_length*2 >= np.sqrt((port2.x-port1.x) ** 2 + (port2.y-port1.y) ** 2):
        print('You tried to make a longhaul waveguide in too short a space.  Reduce the taper length or move ports farther apart!')
        return


    Taper1=pg.taper(length=taper_length,width1=port1.width,width2=max_wg_width,layer=layer)
    Taper2=pg.taper(length=taper_length,width2=port2.width,width1=max_wg_width,layer=layer)
    taper1=D.add_ref(Taper1)
    taper2=D.add_ref(Taper2)
    taper1.connect(port=1,destination=port1)
    taper2.connect(port=2,destination=port2)
    R=pr.route_basic(port1=taper1.ports[2],port2=taper2.ports[1],layer=layer)
    r=D.add_ref(R)
    D.add_port(name=1,port=taper1.ports[1])
    D.add_port(name=2,port=taper2.ports[2])

    return D

#beam-tap with nanospike that turns out 90 deg after the coupling portion
#this version calculates the coupling length and gap for a given coupling coefficient
def nano_tap_90(
        bend_radius=NANOTAP_BEND_RADIUS,
        wg_width=SM_WG_WIDTH,
        cc=0.01,
        rate_multiplier=1.0,
        layerset=lys
):
    layer=layerset['wg_deep']
    # port 0 is input on bus
    # port 1 is output on bus
    # port 2 is tap port, which receives a power coupling coefficient 'cc'
    D = Device()
    taper_length = 2
    min_taper_width = 0.05

    Taper = pg.taper(width1=min_taper_width, width2=wg_width, length=taper_length, layer=layer)
    taper = D.add_ref(Taper)

    #good minimum coupling length based on FDTD sims gives 3 um (220 nm si at 1220 nm)

    #the a1,a2,a3 values are fitting parameters for a sin^2 coupling coefficient
    #a1 is assumed to be 1 though in all cases.
    if cc<0.0032:
        gap=0.35
        a1=1
        a2=.0059
        a3=.011
    elif cc<0.0316:
        gap=0.3
        a1=1
        a2=.0107
        a3=.0281
    else:
        gap=0.15
        a1=1
        a2=.0666
        a3=.0991

    a2=a2*rate_multiplier

    coupling_length = (np.arcsin(np.sqrt(cc))-a3)/a2
    if coupling_length > 25:
        print('Warning: required coupling length goes beyond currently experimentally measured data.')

    Coupler = pg.taper(width1=wg_width, width2=wg_width, length=coupling_length, layer=layer)
    coupler = D.add_ref(Coupler)
    coupler.connect(port=1, destination=taper.ports[2])



    Arc = pr._gradual_bend(radius=bend_radius/0.9,width=wg_width,layer=layer,angular_coverage=3,direction='cw',num_steps=2)
    arc = D.add_ref(Arc)
    arc.connect(port=1,destination=coupler.ports[2])


    Bus = pg.taper(width1=wg_width, width2=wg_width, length=taper_length + coupling_length + bend_radius + 0.1, layer=layer)
    bus = D.add_ref(Bus)

    bus.ymin = taper.ymax + gap
    bus.xmin = taper.xmin

    D.add_port(name=1, port=bus.ports[1])
    D.add_port(name=2, port=bus.ports[2])
    D.add_port(name=3, port=arc.ports[2])

    return D

#beam-tap with nanospike that turns out 90 deg after the coupling portion
#this version accepts gap and coupling length as inputs
def nano_tap_90_simple(
        bend_radius=NANOTAP_BEND_RADIUS,
        wg_width=0.35,
        gap=0.3,
        coupling_length=10,
        layerset=lys
):
    layer=layerset['wg_deep']
    # port 0 is input on bus
    # port 1 is output on bus
    # port 2 is tap port, which receives a power coupling coefficient 'cc'
    D = Device()
    taper_length = 2
    min_taper_width = 0.05

    Taper = pg.taper(width1=min_taper_width, width2=wg_width, length=taper_length, layer=layer)
    taper = D.add_ref(Taper)

    Coupler = pg.taper(width1=wg_width, width2=wg_width, length=coupling_length, layer=layer)
    coupler = D.add_ref(Coupler)
    coupler.connect(port=1, destination=taper.ports[2])

    Arc = pg.arc(radius=bend_radius, theta=-90, width=wg_width, layer=layer)
    arc = D.add_ref(Arc)

    arc.connect(port=1, destination=coupler.ports[2])

    Bus = pg.taper(width1=wg_width, width2=wg_width, length=taper_length + coupling_length + bend_radius + 0.1,
                   layer=layer)
    bus = D.add_ref(Bus)

    bus.ymin = taper.ymax + gap
    bus.xmin = taper.xmin

    D.add_port(name=1, port=bus.ports[1])
    D.add_port(name=2, port=bus.ports[2])
    D.add_port(name=3, port=arc.ports[2])

    return D


#simple y-junction splitter consisting of two s-bend waveguides curving away from each other.
def y_junction_curved(wg_width=0.35,length_straight=0.5,y_offset=2,length_split=40,put_block=True,layerset=lys):
    layer=layerset['wg_deep']
    D=Device()
    B=pg.taper(length=length_straight,width1=wg_width,width2=wg_width,layer=layer)
    bin=D.add_ref(B)

    bup=D.add_ref(B)

    bup.xmin=bin.xmax+length_split
    bup.y=bin.y+y_offset/2

    bdown=D.add_ref(B)

    bdown.y=bup.y-y_offset
    bdown.x=bup.x

    R=pr.route_basic(port1=bin.ports[2],port2=bup.ports[1],layer=layer)
    r=D.add_ref(R)
    R=pr.route_basic(port1=bin.ports[2],port2=bdown.ports[1],layer=layer)
    r=D.add_ref(R)

    D.add_port(name=1,port=bin.ports[1])
    D.add_port(name=2,port=bup.ports[2])
    D.add_port(name=3,port=bdown.ports[2])

    if put_block:
        block_wid = .11  # determined by fabrication
        orig_port_coords = (D.xmin, D.y)
        bbox = pg.rectangle((D.xsize, D.ysize))
        bbox.xmin = orig_port_coords[0]
        bbox.y = orig_port_coords[1]
        inverse = pg.boolean(bbox, D, 'A-B')
        # remove the ones on the outside
        for poly in inverse.polygons:
            if poly.xmin - orig_port_coords[0] < .01:
                inverse.remove(poly)
        clipper = pg.rectangle((D.xsize, block_wid))
        clipper.xmin = orig_port_coords[0]
        clipper.y = orig_port_coords[1]
        central = pg.boolean(inverse, clipper, 'and')
        interstitial = pg.boolean(inverse, clipper, 'A-B')
        interstitial.remove(interstitial.polygons[0])

        yet_another = pg.rectangle((interstitial.xsize, bbox.ysize))
        yet_another.y = central.y
        yet_another.xmin = interstitial.xmin
        central = pg.boolean(central, yet_another, 'A-B')

        circular_cut = pg.circle(radius=block_wid/2)
        circular_cut.y = central.y
        circular_cut.x = central.xmax
        final_block = D << pg.boolean(central, circular_cut, 'A-B', layer=layer)
        # final_block.xmin = orig_port_coords[0]
        # final_block.y = orig_port_coords[1]

    return D


from phidl import device_lru_cache

mmi1x2_nominals = dict(wg_width=0.35, length_port=0.2, length_mmi=2.8, width_mmi=1.55, gap_mmi=0.4)

# @device_lru_cache
def mmi1x2(wg_width=0.35,length_port=0.2,length_mmi=2.8,width_mmi=1.55,gap_mmi=0.4,layer=lys['wg_deep']):

    D=Device()

    Port_wg=pg.taper(length=length_port,width1=wg_width,width2=wg_width,layer=layer)
    port_in=D.add_ref(Port_wg)
    MMI=pg.taper(length=length_mmi,width1=width_mmi,width2=width_mmi,layer=layer)
    mmi=D.add_ref(MMI)

    mmi.connect(port=1,destination=port_in.ports[2])

    port_up=D.add_ref(Port_wg)
    port_up.connect(port=1,destination=mmi.ports[2])
    port_up.movey(gap_mmi)

    port_down=D.add_ref(Port_wg)
    port_down.connect(port=1,destination=mmi.ports[2])
    port_down.movey(-gap_mmi)

    D.add_port(name=1,port=port_in.ports[1])
    D.add_port(name=2,port=port_up.ports[2])
    D.add_port(name=3,port=port_down.ports[2])

    D.flatten()
    return D

#s-bend nano beam tap that accepts coupling coefficient and determines appropriate gap and length for coupler
def nano_tap_s(
        s_offset=4,
        s_length=25,
        wg_width=SM_WG_WIDTH,
        bend_radius=NANOTAP_BEND_RADIUS,
        cc=0.01,
        rate_multiplier=1.0,
        layerset=lys
):
    layer=layerset['wg_deep']
    #port 0 is input on bus
    #port 1 is output on bus
    #port 2 is tap port, which receives a power coupling coefficient 'cc'
    D=Device()
    taper_length=2
    min_taper_width=0.05
    gap=0.3

    Taper=pg.taper(width1=min_taper_width,width2=wg_width,length=taper_length,layer=layer)
    taper=D.add_ref(Taper)

    #good minimum coupling length based on FDTD sims gives 3 um (220 nm si at 1220 nm)

    if cc<0.0032:
        gap=0.35
        a1=1
        a2=.0059
        a3=.011
    elif cc<0.0316:
        gap=0.3
        a1=1
        a2=.0107
        a3=.0281
    else:
        gap=0.15
        a1=1
        a2=.0666
        a3=.0991

    a2=a2*rate_multiplier

    coupling_length = (np.arcsin(np.sqrt(cc))-a3)/a2
    if coupling_length > 25:
        print('Warning: required coupling length goes beyond currently experimentally measured data.')

    Coupler=pg.taper(width1=wg_width,width2=wg_width,length=coupling_length,layer=layer)
    coupler=D.add_ref(Coupler)
    coupler.connect(port=1,destination=taper.ports[2])

    Out_wg=pg.taper(width1=wg_width,width2=wg_width,length=0.1,layer=layer)
    out_wg=D.add_ref(Out_wg)
    out_wg.xmin=coupler.xmax+s_length

    Bus=pg.taper(width1=wg_width,width2=wg_width,length=taper_length+coupling_length+s_length+0.1,layer=layer)
    bus=D.add_ref(Bus)

    bus.ymin=taper.ymax+gap
    bus.xmin=taper.xmin
    out_wg.y = bus.y - s_offset

    S=pr.route_basic(port1=coupler.ports[2],port2=out_wg.ports[1],num_path_pts=30,layer=layer)
    s=D.add_ref(S)

    D.add_port(name=1,port=bus.ports[1])
    D.add_port(name=2,port=bus.ports[2])
    D.add_port(name=3,port=out_wg.ports[2])

    return D

#s-bend nano beam tap that accepts coupling length and gap
def nano_tap_s_simple(
        s_offset=4,
        s_length=25,
        wg_width=SM_WG_WIDTH,
        gap=0.3,
        bend_radius=NANOTAP_BEND_RADIUS,
        coupling_length=10,
        layerset=lys
):
    layer=layerset['wg_deep']
    #port 0 is input on bus
    #port 1 is output on bus
    #port 2 is tap port, which receives a power coupling coefficient 'cc'
    D=Device()
    taper_length=2
    min_taper_width=0.05

    Taper=pg.taper(width1=min_taper_width,width2=wg_width,length=taper_length,layer=layer)
    taper=D.add_ref(Taper)

    Coupler=pg.taper(width1=wg_width,width2=wg_width,length=coupling_length,layer=layer)
    coupler=D.add_ref(Coupler)
    coupler.connect(port=1,destination=taper.ports[2])

    Out_wg=pg.taper(width1=wg_width,width2=wg_width,length=0.1,layer=layer)
    out_wg=D.add_ref(Out_wg)
    out_wg.xmin=coupler.xmax+s_length

    Bus=pg.taper(width1=wg_width,width2=wg_width,length=taper_length+coupling_length+s_length+0.1,layer=layer)
    bus=D.add_ref(Bus)

    bus.ymin=taper.ymax+gap
    bus.xmin=taper.xmin
    out_wg.y = bus.y - s_offset

    S=pr.route_basic(port1=coupler.ports[2],port2=out_wg.ports[1],num_path_pts=30,layer=layer)
    s=D.add_ref(S)

    # Arc = pg.arc(radius=bend_radius, theta=-90, width=wg_width, layer=0)
    # arc = D.add_ref(Arc)
    # arc.connect(port=1,destination=out_wg.ports[2])

    D.add_port(name=1,port=bus.ports[1])
    D.add_port(name=2,port=bus.ports[2])
    D.add_port(name=3,port=out_wg.ports[2])

    return D

#component allowing arbitrary number of consecutive in-plane crossings
#ports 1 and 2 are the main input/output respectively (left, right)
#successive ports: the top and bottom of each add'l crossing, going from left to right. (need to double check that though)
def in_plane_crossing(
        wg_width=SM_WG_WIDTH,
        mm_width=MM_WG_WIDTH,
        taper_length=1.1,
        offset=1.45,
        period=2.55,
        num_crossings=5,
        layerset=lys

):
    D=Device()
    layer=layerset['wg_deep']
    #first, taper from sm to mm waveguide nonadiabatically
    Taper=pg.taper(length=taper_length,width1=wg_width,width2=mm_width,layer=layer)
    taper_in=D.add_ref(Taper)
    taper_in.xmin=0
    taper_in.y=0

    #add the horizontal mm bus going through all the vertical crossings
    Main_bus=pg.taper(length=2*offset+mm_width+(num_crossings-1)*period,width1=mm_width,width2=mm_width,layer=layer)
    main_bus=D.add_ref(Main_bus)
    main_bus.connect(port=1,destination=taper_in.ports[2])

    #add the out taper (same as in taper but on other side and reversed)
    taper_out=D.add_ref(Taper)
    taper_out.connect(port=2,destination=main_bus.ports[2])

    D.add_port(name=1,port=taper_in.ports[1])
    D.add_port(name=2,port=taper_out.ports[1])


    x_position=taper_in.xmax+offset+mm_width/2
    #fill in the space between the tapers with the right number of crossings
    for x in range(num_crossings):
        V_bus=pg.taper(length=2*offset+mm_width,width1=mm_width,width2=mm_width,layer=layer)
        v_bus=D.add_ref(V_bus)
        v_bus.rotate(angle=90)
        v_bus.y=0
        v_bus.x=x_position
        taper_current=D.add_ref(Taper)
        taper_current.connect(port=2,destination=v_bus.ports[2])
        D.add_port(name=2*x+3,port=taper_current.ports[1])

        taper_current=D.add_ref(Taper)
        taper_current.connect(port=2,destination=v_bus.ports[1])
        D.add_port(name=2*x+4,port=taper_current.ports[1])
        x_position=x_position+period
    return D


#find a coupling coefficient for a specific output given a distribution pattern and number of outputs
#it is a little silly to recalculate the whole distribution every time, but who cares
def get_cc(
        index=0,
        num_neurons=10,
        function='uniform',
        fwhm=5
):
    if function=='uniform':
        cc=1/(num_neurons-index)
        #print(cc)
    elif function=='gaussian':
        #simulate the vals first...
        desired_power=[]
        cc=[]
        ptap=[]
        for i in range(num_neurons):
            desired_power.append(np.exp(-4*(np.log(2)*((i)-num_neurons/2)**2)/((fwhm)**2)))
        area_under_gaussian=np.sum(desired_power)
        #print(area_under_gaussian)
        p=1
        #normalization and calculation of new cc
        for i in range(num_neurons):
            #normalize this value first
            desired_power[i]=desired_power[i]/area_under_gaussian

            #calc the cc given the current power fraction remaining
            cc.append(desired_power[i]/p)

            #next, find remaining power after this tap.
            ptap.append(cc[i]*p)
            p=p-cc[i]*p
        #print(cc)
        cc=cc[index]
        #print(ptap)

    return cc

def bs_tree(
    tree_depth = 4,
    output_separation = 20,
    layer_separation = 50,
    wg_width = SM_WG_WIDTH,
    splitter=None):

    D = Device()

    if splitter is None:
        y = mmi1x2()  # use the defaults
    else:
        y = splitter
    Treelayer = Device()
    Treelayer << y
    Treelayer.add_port(port = y.ports[1], name = 0)
    Treelayer.add_port(port = y.ports[2], name = 1)
    Treelayer.add_port(port = y.ports[3], name = 2)

    D << Treelayer
    D.add_port(port = Treelayer.ports[0], name = 1)

    for i in range(1,tree_depth):
        num_port_out = 0
        num_port_in = 1
        T = Device()
        yprev = 0
        for j in range(0, 2**i):
            ycurrent = T << y
            ycurrent.xmin = 0
            ycurrent.y = yprev - 2**(tree_depth-i)*output_separation
            T.add_port(port = ycurrent.ports[1], name = -1*num_port_in)
            T.add_port(port = ycurrent.ports[2], name = num_port_out + 1)
            T.add_port(port = ycurrent.ports[3], name = num_port_out + 2)
            num_port_out = num_port_out + 2
            num_port_in = num_port_in + 1
            j = j+1
            yprev = ycurrent.y
        t = D << T
        for j in range(0, 2**i):
            t.y = Treelayer.y
            t.xmin = Treelayer.xmax + layer_separation#*(tree_depth-i-1)
            D << pr.route_manhattan(port1 = t.ports[-1*(j+1)], port2 = Treelayer.ports[j+1], bendType = 'gradual', layer = lys['wg_deep'])
        Treelayer = t

    midpoint_p = [t.xmax+layer_separation,t.ports[1].midpoint[1]+output_separation/2]

    for j in range(0, 2**tree_depth):
        numport = j+2
        D.add_port(name = numport, orientation = 180, width = wg_width, midpoint = midpoint_p)
        midpoint_p[1] = midpoint_p[1] - output_separation
        D << pr.route_manhattan(port1=D.ports[numport], port2 = t.ports[j+1], bendType = 'gradual', layer = lys['wg_deep'])
        D.ports[numport].orientation = 0

    return D

def gradual_bend(
    radius = 20,
    width = 1.0,
    theta = 90,
    angular_coverage=15,
    num_steps=10,
    angle_resolution=0.1,
    start_angle=0,
    direction='ccw',
    layer=lys['wg_deep'],
    ):

    """
    creates a theta-degree bent waveguide
    the bending radius is gradually increased until it reaches the minimum
    value of the radius at the "angular coverage" angle.
    it essentially creates a smooth transition to a bent waveguide mode.
    user can control number of steps provided.
    direction determined by start angle and cw or ccw switch
    ############
    with the default 10 "num_steps" and 15 degree coverage, effective radius is about 1.5*radius.
    """
    circular_arcangle = theta / 2 - angular_coverage
    if circular_arcangle < 0:
        raise ValueError('Cannot make the gradual section bigger than the whole bend')

    angular_coverage=np.deg2rad(angular_coverage)
    D = Device()

    #determines the increment in radius through its inverse from 0 to 1/r
    inc_rad = (radius**-1)/(num_steps)
    angle_step = angular_coverage/num_steps

    #construct a series of sub-arcs with equal angles but gradually decreasing bend radius
    arcs = []
    for x in range(num_steps):
        A = pr._arc(radius=1/((x+1)*inc_rad), width=width, theta=np.rad2deg(angle_step),
                    angle_resolution=angle_resolution, layer=layer)
        a = D.add_ref(A)
        arcs.append(a)
        if x>0:
            a.connect(port=1,destination=prevPort)
        prevPort=a.ports[2]
    D.add_port(name=1,port=arcs[0].ports[1])

    #now connect a regular bend for the normal curved portion
    B = pr._arc(radius=radius,width=width,theta=circular_arcangle,
                angle_resolution=angle_resolution,layer=layer)
    b = D.add_ref(B)
    b.connect(port=1,destination=prevPort)
    D.add_port(name=2,port=b.ports[2])

    #now create the overall structure
    Total = Device()

    #clone the half-curve into two objects and connect for a 90 deg bend.
    D1 = Total.add_ref(D)
    D2 = Total.add_ref(D)
    D2.reflect()
    D2.connect(port=2,destination=D1.ports[2])
    Total.xmin=0
    Total.ymin=0

    #orient to default settings...
    Total.reflect(p1=[0,0],p2=[1,1])
    Total.reflect(p1=[0,0],p2=[1,0])

    #orient to user-provided settings
    if direction == 'cw':
        Total.reflect(p1=[0,0],p2=[1,0])
    Total.rotate(angle=start_angle,center=Total.center)
    Total.center=[0,0]
    Total.add_port(name=1,port=D1.ports[1])
    Total.add_port(name=2,port=D2.ports[1])

    return Total

###############################################################################
#
# SONIA BUCKLEY DEVICES
# OPTOELECTRONICS
#
###############################################################################


def contact_pads(size = WIREBOND_SIZE, pad_overlap = PAD_OVERHANG,
    metal_layers = ['m5_wiring'],layerset = lys, packaging = True, device_info = {}):

    # width = pad width
    # length = pad length

    D = Device()

    Pad = Device()
    for i,metal_layer in enumerate(metal_layers):
        pad = Pad.add_ref(pg.compass(size = size, layer = layerset[metal_layer]))

    width = size[0]
    length = size[1]

    if metal_layers == ['m5_wiring']:
        po2 = pg.compass(size = (length-pad_overlap*2, width-pad_overlap*2), layer = layerset['v5'])
        D.add_ref(po2)

    pad = D << Pad

    p = D.add_port(width=length, orientation =0, midpoint = D.center, name = 1)

    if packaging == True:
        glue_dev = D << gluebox(size = (WIREBOND_SIZE[0]+25,WIREBOND_SIZE[1]+25), corner_cut = 70, wall_width = GLUEBOX_WALL, layer = lys['su8'])
        glue_dev.center = pad.center

    #p.info.update(device_info)
    #D.info.update(device_info)
    #D.flatten()
    return D

def led_pads(size = LED_PAD_SIZE, pad_overlap = PAD_OVERHANG, layerset = lys):

    # width = pad width
    # length = pad length
    D = Device()

    width = size[0]
    length = size[1]
    full_etch_layer_overlap = [-pad_overlap, pad_overlap, pad_overlap, pad_overlap, pad_overlap]#W, E, N, S, the full etch is outside the TiAu
    wiring_layer_overlap = [-pad_overlap, pad_overlap, pad_overlap, pad_overlap]#W, E, N, S, the wiring layer is outside the full etch
    partial_etch_layer_overlap = [0, 2*pad_overlap, pad_overlap, pad_overlap]#W, E, N, S, the partial etch is outside the wiring layer

    size_wg_full = (length+full_etch_layer_overlap[0]+full_etch_layer_overlap[1],
                    width+full_etch_layer_overlap[2]+full_etch_layer_overlap[3])

    size_wiring_pad = (size_wg_full[0]+wiring_layer_overlap[0]+wiring_layer_overlap[1],
                       size_wg_full[1]+wiring_layer_overlap[2]+wiring_layer_overlap[3])

    size_wg_partial = (size_wiring_pad[0]+partial_etch_layer_overlap[0]+partial_etch_layer_overlap[1],
                       size_wiring_pad[1]+partial_etch_layer_overlap[2]+partial_etch_layer_overlap[3])

    pad_m4 = D.add_ref(pg.compass(size = (length, width), layer = layerset['m4_ledpad']))
    pad_wg_full = D.add_ref(pg.compass(size = size_wg_full, layer = layerset['wg_deep']))
    pad_m5 = D.add_ref(pg.compass(size = size_wiring_pad, layer = layerset['m5_wiring']))
    pad_wg_partial = D.add_ref(pg.compass(size = size_wg_partial, layer = layerset['wg_shallow']))

    pad_m5.xmin = pad_m4.xmin-wiring_layer_overlap[0]
    pad_wg_full.xmin = pad_m4.xmin
    pad_wg_partial.xmin = pad_m4.xmin

    D.add_port(midpoint = pad_m5.center, width = pad_m5.size[1], orientation = 0, name = 1)

    return D

def nw_pads(via_size = VIA_SIZE, pad_overlap = PAD_OVERHANG, is_gnd = False, layerset = lys):

    # does not include deep etch layer for now

    D = Device()

    wiring_overlap = [pad_overlap, pad_overlap, pad_overlap, pad_overlap]
    wiring_size = [via_size[0]+wiring_overlap[0]+wiring_overlap[1], via_size[1]+wiring_overlap[2]+wiring_overlap[3]]

    tiau_overlap = [pad_overlap*2, pad_overlap*2, pad_overlap*2, pad_overlap*2]#W,E,N,S
    tiau_size = (wiring_size[0]+tiau_overlap[0]+tiau_overlap[1],
        wiring_size[1]+tiau_overlap[2]+tiau_overlap[3])

    wsi_overlap = [pad_overlap, pad_overlap, pad_overlap, pad_overlap]
    wsi_size = (tiau_size[0]+wsi_overlap[0]+wsi_overlap[1],
                tiau_size[1]+wsi_overlap[2]+wsi_overlap[3])

    via = D.add_ref(pg.compass(size = via_size, layer = layerset['v3']))
    wiring = D.add_ref(pg.compass(size = wiring_size, layer = layerset['m5_wiring']))
    tiau = D.add_ref(pg.compass(size = tiau_size, layer = layerset['m1_nwpad']))
    wsi = D.add_ref(pg.compass(size = wsi_size, layer = layerset['m2_nw']))

    D.add_port(name = 1, midpoint = wiring.center, width = wiring.size[1], orientation = 0)

    return D

def chip_marks(size = (10000,10000),
               width = 30,
               length = 200,
               center = [0,0],
               layerset = lys,
               metal_layer_str = 'm5_wiring'):

    D = Device()
    cornermark_template = Device()
    rectA = cornermark_template.add_ref(pg.compass(size = (length, width), layer = layerset[metal_layer_str]))
    rectB = cornermark_template.add_ref(pg.compass(size = (width, length-width), layer = layerset[metal_layer_str]))
    rectA.xmin = rectB.xmin
    rectA.ymin = rectB.ymax

    coords = [[-size[0]/2, size[1]/2],[size[0]/2, size[1]/2],
               [size[0]/2, -size[1]/2],[-size[0]/2, -size[1]/2]]
    rotations = [0, -90, -180, -270]

    for i in range(0,4):
        cornermark = D.add_ref(cornermark_template)
        cornermark.rotate(rotations[i])
        if coords[i][0]>0:
            cornermark.xmin = coords[i][0]+center[0]
        elif coords[i][0]<0:
            cornermark.xmax = coords[i][0]+center[0]
        else:
            print('you put something weird in for the x coordinates')
            return

        if coords[i][1]>0:
            cornermark.ymax = coords[i][1]+center[1]
        elif coords[i][1]< 0:
            cornermark.ymin = coords[i][1]+center[1]
        else:
            print('you put something weird in the y coordinates')


    return D


def contactdopant_resistance_test(test_length = 20.,
             width = 4.,
             nsquares = None,
             pad_overlap = PAD_OVERHANG,
             led_pad_device = None,
             contact_dopant = 'junk',
             layerset = lys):

    # test_length :             The length of dopant (not including under the contacts)
    # width :                    Width of the dopant region to be tested
    # pad_offset :              Distance of wirebonding pads to device
    # led_pad_overhang :             Overhang of small contact pads with shallow etch/dopant contact region
    # contact_length :          The length of the small pad contact
    # contact_dopant:     Type of dopant
    # pad_device:               Phidl pad device

    D = Device()
    Contact = Device()
    dopant_layer = layerset[contact_dopant]

    # pads to make contact with device
    smallpad = Contact.add_ref(led_pad_device)

    if nsquares is not None:
        if test_length is None or test_length == 0:
            test_length = nsquares*width
        elif width is None or width == 0:
            width = nsquares/test_length
        else:
            print('overdefined; only use two of test_length, width and nsquares')
            return
    elif nsquares is None and test_length is not None and width is not None:
        nsquares = length*width
    else:
        print('Underdefined')

    # add the contact dopant under the pad
    contact_dopant_overlap = [0, pad_overlap, pad_overlap, pad_overlap] #W,E,N,S
    contact_dopant_size = (smallpad.size[0] + contact_dopant_overlap[0]+contact_dopant_overlap[1],
                           smallpad.size[1] + contact_dopant_overlap[2]+contact_dopant_overlap[3])

    contact_area = Contact.add_ref(pg.compass(size = contact_dopant_size, layer = layerset[contact_dopant]))
    smallpad.xmin = contact_area.xmin
    Contact.add_port(name = 'dopant_end', midpoint = contact_area.ports['E'].midpoint, width = contact_area.ports['E'].width, orientation = contact_area.ports['E'].orientation)
    Contact.add_port(smallpad.ports[1])

    contact1 = D.add_ref(Contact)
    contact1.rotate(180)
    contact2 = D.add_ref(Contact)


    # make the resistance test dopant, deep and shallow etches

    partial_etch_offset = [0, 0, pad_overlap, pad_overlap]
    partial_etch_size = (test_length+partial_etch_offset[0]+partial_etch_offset[1],
                         width+partial_etch_offset[2]+partial_etch_offset[3])

    dopant_offset = [0,0,pad_overlap,pad_overlap]
    dopant_size = (partial_etch_size[0]+partial_etch_offset[0]+partial_etch_offset[1],
                         partial_etch_size[1]+partial_etch_offset[2]+partial_etch_offset[3])


    Res_test = pg.compass(size = dopant_size, layer = dopant_layer)
    res_test = D.add_ref(Res_test)
    Deep_etch= pg.compass(size = (test_length, width), layer = layerset['wg_deep'])
    deep_etch = D.add_ref(Deep_etch)
    Shallow_etch= pg.compass(size = partial_etch_size, layer = layerset['wg_shallow'])
    shallow_etch = D.add_ref(Shallow_etch)

    contact2.xmin = res_test.xmax
    contact1.xmax = res_test.xmin

    D.add_port(name = 'wiring1', port = contact1.ports[1])
    D.add_port(name = 'wiring2', port = contact2.ports[1])

    # Calculate meta-data

    contact_area = led_pad_device.size[0]*led_pad_device.size[1]

    if contact_dopant == 'dp_p+':
        expected_resistance = EXPECTED_RCONTACT_PPLUS/contact_area+EXPECTED_RSQ_PPLUS*nsquares
    elif contact_dopant == 'dp_n+':
        expected_resistance = EXPECTED_RCONTACT_NPLUS/contact_area+EXPECTED_RSQ_NPLUS*nsquares
    else:
        expected_resistance = 0

    D.info['nsquares'] = nsquares
    D.info['contact_dopant_type'] = contact_dopant
    D.info['contact_area'] = contact_area
    D.info['expected_resistance'] = expected_resistance

    return D

def wgdopant_resistance_test(test_length = 100.,
                        width = 4.,
                        nsquares = None,
                        led_pad_device = None,
                        pad_overlap = PAD_OVERHANG,
                        test_dopant = 'junk',
                        contact_dopant = 'junk',
                        layerset = lys
                        ):

    # test_length :             The length of dopant (not including under the contacts)
    # width :                   Width of the dopant region to be tested
    # pad_offset :              Distance of wirebonding pads to device
    # led_pad_overhang :             Overhang of small contact pads with shallow etch/dopant contact region
    # contact_length :          The length of the small pad contact
    # contact_dopant:           Type of dopant
    # test_dopant:              Type of dopant
    # pad_device:               Phidl pad device

    D = Device()
    if nsquares is not None:
        if test_length is None or test_length == 0:
            test_length = nsquares*width
        elif width is None or width== 0:
            width = nsquares/test_length
        else:
            print('use only two of test_length, nsquares and width')
    elif (nsquares is None) and ((test_length is not None) and (width is not None)):
        nsquares = test_length*width
    else:
        print('Underdefined')

    #make the contact area of dopant
    Plus_contact = Device()
    smallpad = Plus_contact.add_ref(led_pad_device)
    smallpad.rotate(180)

    contact_dopant_overlap = [0, pad_overlap, pad_overlap, pad_overlap] #W,E,N,S
    contact_dopant_size = (smallpad.size[0] + contact_dopant_overlap[0]+contact_dopant_overlap[1],
                           smallpad.size[1] + contact_dopant_overlap[2]+contact_dopant_overlap[3])

    plus_area = Plus_contact.add_ref(pg.compass(size = contact_dopant_size, layer = layerset[contact_dopant]))
    smallpad.xmax = plus_area.xmax
    Plus_contact.add_port(name = 'dopant_end', midpoint = plus_area.ports['E'].midpoint, width = plus_area.ports['E'].width, orientation = plus_area.ports['E'].orientation)
    Plus_contact.add_port(smallpad.ports[1])

    plus_contact1 = D.add_ref(Plus_contact)
    plus_contact1.rotate(180)
    plus_contact2 = D.add_ref(Plus_contact)

    # Make the deep etch area
    Deep_etch= pg.compass(size = (test_length, width), layer = layerset['wg_deep'])
    deep_etch = D.add_ref(Deep_etch)
    plus_contact1.connect(port = 'dopant_end', destination = deep_etch.ports['E'])
    plus_contact2.connect(port = 'dopant_end', destination = deep_etch.ports['W'])

    if test_dopant is not None:
        test_region_overlap = [0, 0, pad_overlap, pad_overlap]#W,E,N,S
        test_region_size = (deep_etch.size[0]+test_region_overlap[0]+test_region_overlap[1],
            deep_etch.size[1]+test_region_overlap[2]+test_region_overlap[3])
        test_region = D.add_ref(pg.compass(size = test_region_size, layer = layerset[test_dopant]))
        test_region.center = deep_etch.center

    D.add_port(name = 'wiring1', port = plus_contact2.ports[1])
    D.add_port(name = 'wiring2', port = plus_contact1.ports[1])

    # Calculate meta-data

    contact_area = led_pad_device.size[0]*led_pad_device.size[1]

    if test_dopant == 'dp_p+':
        expected_resistance_slope = nsquares*EXPECTED_RSQ_PPLUS
    elif test_dopant == 'dp_n+':
        expected_resistance_slope = nsquares*EXPECTED_RSQ_NPLUS
    elif test_dopant == 'dp_p':
        expected_resistance_slope = nsquares*EXPECTED_RSQ_NPLUS
    elif test_dopant == 'dp_n':
        expected_resistance_slope = nsquares*EXPECTED_RSQ_NPLUS
    elif test_dopant == 'dp_e':
        expected_resistance_slope = nsquares*EXPECTED_RSQ_E
    elif test_dopant is None:
        expected_resistance_slope = nsquares*EXPECTED_RSQ_I
    else:
        expected_resistance_slope = 0

    if contact_dopant == 'dp_p+':
        expected_resistance_contact = EXPECTED_RCONTACT_PPLUS/contact_area
    elif contact_dopant == 'dp_n+':
        expected_resistance_contact = EXPECTED_RCONTACT_NPLUS/contact_area
    else:
        expected_resistance_contact = 0

    expected_resistance = 2*expected_resistance_contact + expected_resistance_slope

    D.info['nsquares'] = nsquares
    D.info['test_dopant_type'] = test_dopant
    D.info['contact_dopant_type'] = contact_dopant
    D.info['contact_area'] = contact_area
    D.info['expected_resistance'] = expected_resistance

    return D

def diode_test(contact_dopant_length = 5.,
               wg_dopant_length = 5.,
               e_length = 10.,
               width = 4.,
               pad_offset = 100.,
               led_pad_device = None,
               pad_overlap = PAD_OVERHANG,
               exists_e_region = False,
               label=None,
               layerset = lys):

    # contact_length :      contact length
    # contact_dopant_length :   distance of pad from LED waveguide ridge
    # wg_dopant_length :       length of p+ or n+ dopant that has gold on it
    # e_length :            the length of the emissive center region. If exists_e_region = False this doesn't mean anything.
    #                       If exists_pn_region = False this overrides to e_length = waveguide_width
    # width :               width of the diode
    # led_pad_overhang:          Overhang of the pad past the dopant regions
    # pad_offset:           Distance of pad from the rest of the device
    # exists_e_region:          True means emissive centers present

    D = Device()

    #make the contact area of dopant
    #make the contact area of dopant


    Pplus_contact = Device()
    smallpad = Pplus_contact.add_ref(led_pad_device)
    smallpad.rotate(180)

    contact_dopant_overlap = [0, pad_overlap, pad_overlap, pad_overlap] #W,E,N,S
    contact_dopant_size = (smallpad.size[0] + contact_dopant_overlap[0]+contact_dopant_overlap[1],
                           smallpad.size[1] + contact_dopant_overlap[2]+contact_dopant_overlap[3])

    pplus_area = Pplus_contact.add_ref(pg.compass(size = contact_dopant_size, layer = layerset['dp_p+']))
    smallpad.xmax = pplus_area.xmax
    Pplus_contact.add_port(name = 'dopant_end', midpoint = pplus_area.ports['E'].midpoint, width = pplus_area.ports['E'].width, orientation = pplus_area.ports['E'].orientation)
    Pplus_contact.add_port(smallpad.ports[1])

    #
    Nplus_contact = Device()
    smallpad = Nplus_contact.add_ref(led_pad_device)
    nplus_area = Nplus_contact.add_ref(pg.compass(size = contact_dopant_size, layer = layerset['dp_n+']))
    smallpad.xmin = nplus_area.xmin
    Nplus_contact.add_port(name = 'dopant_end', midpoint = nplus_area.ports['W'].midpoint, width = nplus_area.ports['W'].width, orientation = nplus_area.ports['W'].orientation)
    Nplus_contact.add_port(smallpad.ports[1])

    #
    P = Device()
    N = Device()
    #
    pplus_contact = P.add_ref(Pplus_contact)
    nplus_contact = N.add_ref(Nplus_contact)
    ##

    P.add_port(port = Pplus_contact.ports[1], name = 'wiring1')
    N.add_port(port = Nplus_contact.ports[1], name = 'wiring2')


    # Make the shallow etch region of the pplus and nplus
    overlap_pe_plus_region = [0,0,pad_overlap, pad_overlap]#W,E,N,S
    size_pe_plus_region = (contact_dopant_length+overlap_pe_plus_region[0]+overlap_pe_plus_region[1],
                         width+overlap_pe_plus_region[2]+overlap_pe_plus_region[3])
    pe_pplus_region = P.add_ref(pg.compass(size = size_pe_plus_region, layer = layerset['wg_shallow']))
    pe_nplus_region = N.add_ref(pg.compass(size = size_pe_plus_region, layer = layerset['wg_shallow']))

    #Now add the pplus and nplus regions that aren't under the tiau pad
    overlap_plus_region = [0,0,pad_overlap, pad_overlap]#W,E,N,S
    size_plus_region = (size_pe_plus_region[0]+overlap_plus_region[0]+overlap_plus_region[1],
                         size_pe_plus_region[1]+overlap_plus_region[2]+overlap_plus_region[3])
    pplus_region = P.add_ref(pg.compass(size = size_plus_region, layer = layerset['dp_p+']))
    nplus_region = N.add_ref(pg.compass(size = size_plus_region, layer = layerset['dp_n+']))


    # postion them
    pe_pplus_region.y = pplus_contact.y
    pplus_region.y = pplus_contact.y
    pe_pplus_region.x = pplus_contact.xmax+contact_dopant_length/2
    pplus_region.x = pplus_contact.xmax+contact_dopant_length/2

    pe_nplus_region.y = nplus_contact.y
    nplus_region.y = nplus_contact.y
    pe_nplus_region.x = nplus_contact.xmin-contact_dopant_length/2
    nplus_region.x = nplus_contact.xmin-contact_dopant_length/2


    # make the p and n regions
    if (wg_dopant_length == 0 or wg_dopant_length is None):
        P.add_port(port = pplus_region.ports['E'], name = 'W')
        N.add_port(port = nplus_region.ports['W'], name = 'E')
        waveguide_width = e_length
    else:
        wgdopant_overlap = [0,0,pad_overlap, pad_overlap]
        size_wgdopant = [wg_dopant_length+wgdopant_overlap[0]+wgdopant_overlap[1], width+wgdopant_overlap[2]+wgdopant_overlap[3]]
        p_region = P.add_ref(pg.compass(size = size_wgdopant, layer = layerset['dp_p']))
        n_region = N.add_ref(pg.compass(size = size_wgdopant, layer = layerset['dp_n']))
        p_region.connect(port = 'W', destination = pplus_region.ports['E'])
        n_region.connect(port = 'E', destination = nplus_region.ports['W'])
    #
        P.add_port(p_region.ports['E'])
        N.add_port(n_region.ports['W'])
    #
        waveguide_width = wg_dopant_length * 2 + e_length
    #
    #    #
    #    ## Now add P and N to the main device
    #    #
    p = D.add_ref(P)
    n = D.add_ref(N)
    ##
    p.xmax = n.xmin
    p.movex(-e_length/2)
    n.movex(e_length/2)

    ####
    #
    ##
    # make the deep etch region to define the waveguide
    de_region = D.add_ref(pg.compass(size = (2*contact_dopant_length+waveguide_width, width), layer = layerset['wg_deep']))
    p.y = de_region.y
    n.y = de_region.y
    de_region.x = (p.xmax+n.xmin)/2

    ## p and n dopanats are in the ridge except for the emissive center region
    D.add_port(port = de_region.ports['N'], name = 'optical1')
    D.add_port(port = de_region.ports['S'], name = 'optical2')

    if exists_e_region:
        if e_length > 0:
            e_region_overlap = [0,0,pad_overlap, pad_overlap]
            e_region_size = (e_length+e_region_overlap[0]+e_region_overlap[1],
                             width+e_region_overlap[2]+e_region_overlap[3])
            e_region = D.add_ref(pg.compass(size = e_region_size, layer = layerset['dp_e']))
            e_region.x = (p.xmax+n.xmin)/2
    #
    D.add_port(port = p.ports['wiring1'])
    D.add_port(port = n.ports['wiring2'])

        # Calculate meta-data

    contact_area = contact_dopant_size[0]*contact_dopant_size[1]

    contact_dopant_length = 5.
    wg_dopant_length = 5
    e_length = 10
    width = 4.

    nsquares_plus= contact_dopant_length/width
    nsquares_wg_dopant = wg_dopant_length/width
    nsquares_e = e_length/width

    contact_resistance = contact_area*(EXPECTED_RCONTACT_PPLUS + EXPECTED_RCONTACT_NPLUS)
    plus_resistance = nsquares_plus*(EXPECTED_RSQ_PPLUS+EXPECTED_RSQ_NPLUS)
    wg_dopant_resistance = nsquares_wg_dopant*(EXPECTED_RSQ_P+EXPECTED_RSQ_N)

    expected_series_resistance = contact_resistance+plus_resistance+wg_dopant_resistance

    if exists_e_region:
        e_region_resistance = nsquares_e*EXPECTED_RSQ_E
        expected_series_resistance += e_region_resistance

    D.info['contact_area'] = contact_area
    D.info['expected_series_resistance'] = expected_series_resistance

    return D

# add labelys, width, include instead of is,
def led(waveguide_width = 0.35,
        pad_distance = 0.4,
        e_length = 0.05,
        LED_length = 1.,
        exists_e_region = True,
        exists_pn_region = True,
        led_pad_device = None,
        label=None,
        layerset = lys):

    # waveguide_width :     width of the LED waveguide ridge
    # pad_distance :        distance of pad from LED waveguide ridge
    # contact length :      length of p+ or n+ dopant that has gold on it
    # e_length :            the length of the emissive center region. If exists_e_region = False this doesn't mean anything.
    #                       If exists_pn_region = False this overrides to e_length = waveguide_width
    # LED_length :          The length of the LED
    # exists_e_region:          True means emissive centers present
    # exists_pn_region:         True means p,n dopants present as well as p+,n+.

    D = Device()

    contact_dopant_length = pad_distance
    if exists_pn_region:
        wg_dopant_length = (waveguide_width-e_length)/2
    else:
        wg_dopant_length = None

    width = LED_length

    led = D.add_ref(diode_test(contact_dopant_length = contact_dopant_length,
                          wg_dopant_length = wg_dopant_length,
                          e_length = e_length,
                          width = LED_length,
                          pad_offset = PAD_OVERHANG,
                          exists_e_region = exists_e_region,
                          led_pad_device = led_pad_device,
                          label=label,
                          layerset = lys
            ))

    D.info = led.info
    D.info['LED length'] = LED_length

    return D

def led_with_wb(waveguide_width = 0.35,
        pad_distance = 0.4,
        contact_length = 5.,
        e_length = 0.05,
        LED_length = 1.,
        exists_e_region = True,
        exists_pn_region = True,
        wb_pad_device = None,
        wb_gnd_pad_device = None,
        layerset = lys):

    pass

def led2(waveguide_width = 0.35,
        pad_distance = 0.4,
        contact_length = 5.,
        e_length = 0.05,
        LED_length = 1.,
        exists_e_region = True,
        exists_pn_region = True,
        layerset = lys):

    pass

def led2_with_wb(waveguide_width = 0.35,
        pad_distance = 0.4,
        contact_length = 5.,
        e_length = 0.05,
        LED_length = 1.,
        exists_e_region = True,
        exists_pn_region = True,
        wb_pad_device = None,
        wb_gnd_pad_device = None,
        layerset = lys):

    pass

def wgnw(meander_width = 0.4, num_squares = 5000.0,
         wgnw_width = 0.1, wgnw_distance = 0.2, wgnw_gap = 0.15,
         wgnw_length = 100, wg_total_length = None,
         nw_pad_device = None):
    ''' The length and width of the meander are chosen so that it is approximately square.
        If nw_pad_device is None, does not add pads. Instead, ports are placed at the meander ends
    '''

    D = Device('wgnw')
    # If given number of squares as input, calculate number of squares in meander.
    meander_pitch = meander_width/MEANDER_FILL_FACTOR

    numsquares_meander = num_squares - 2*wgnw_length/wgnw_width
    if numsquares_meander<1000:
        numsquares_meander=1000

    wgnw_pitch = wgnw_width+wgnw_gap

    meander_length = np.sqrt(numsquares_meander*meander_width*meander_pitch)

    Snspd = pg.snspd(wire_width = meander_width, wire_pitch = meander_pitch, terminals_same_side = False, size = (meander_length,None),
                              num_squares=numsquares_meander, layer = lys['m2_nw'])
    meander = D.add_ref(Snspd)
    numsquares_meander = meander.info['num_squares']
    meander.reflect(p1=(0,0), p2=(1,0))

    wgnw = D.add_ref(pg.optimal_hairpin(width = wgnw_width, pitch = wgnw_pitch, length = wgnw_length, layer = lys['m2_nw']))
    wgnw.reflect(p1 = wgnw.ports[1].midpoint, p2 = wgnw.ports[2].midpoint)
    numsquares_wgnw = 2*wgnw_length/wgnw_width
    #
    Taper = pg.optimal_step(start_width = wgnw_width, end_width = meander_width, num_pts = 50, width_tol = 1e-3,
                     anticrowding_factor = 1.2, layer = lys['m2_nw'])
    taper = D.add_ref(Taper)
    #
    taper.connect(port = 2, destination = meander.ports[1])
    #
    wgnw.xmax = meander.xmin
    wgnw.connect(port = 1, destination = taper.ports[1])

    taper2 = D.add_ref(Taper)
    taper2.reflect(taper2.ports[1].midpoint, taper2.ports[2].midpoint)
    taper2.connect(port = 1, destination = wgnw.ports[2])

    taper3 = D.add_ref(pg.optimal_step(start_width = meander_width, end_width = meander_width*4, num_pts = 50, width_tol = 1e-3,
                     anticrowding_factor = 1.2, layer = lys['m2_nw']))
    taper3.connect(port=1, destination = meander.ports[2])
    bend = D << pg.optimal_90deg(width = meander_width, num_pts = 15,
                                 length_adjust = 1, layer = lys['m2_nw'])
    bend.connect(port=2, destination=taper2.ports[2])
    ##
    # wg layer
    # wg layer wg
    wg_length = wgnw_length + wgnw_distance
    if wg_total_length is not None:
        wg_length = max(wg_length, wg_total_length)
    wg_width = wgnw_width + wgnw_pitch + wgnw_distance * 2
    wg = D.add_ref(pg.compass(size = [wg_length, wg_width], layer = lys['wg_deep']))
    wg.xmax = wgnw.xmax
    wg.y = wgnw.y
    # return D

    # pad breakouts and electrical ports
    if nw_pad_device is not None:
        route1 = D.add_ref(port2pad(port = taper3.ports[2], Pad = nw_pad_device, pad_offset = [nw_pad_device.size[0],0], wiring_layer = lys['m2_nw']))
        route2 = D.add_ref(port2pad(port = bend.ports[1], Pad = nw_pad_device, pad_offset = [nw_pad_device.size[0],nw_pad_device.size[0]/2], wiring_layer = lys['m2_nw']))
        D.add_port(midpoint = route1.ports[1].midpoint, width = route1.ports[1].width, orientation = 0, name = 'wiring1')
        p = D.add_port(port = route2.ports[1], name = 'wiring2')
        p.orientation = 0
    else:
        D.add_port(port = taper3.ports[2], name = 'wiring1')
        D.add_port(port = bend.ports[1], name = 'wiring2')

    # optical port
    D.add_port(name = 'de_edge', port = wg.ports['E'])
    D.add_port(name = 'optical', port = wg.ports['W'])
    D.ports['de_edge'].info['is_waveguide_edge'] = True
    D.move(D.ports['optical'].midpoint, destination=[0, 0])


    numsquares_taper = 3
    D.info['num_squares'] = numsquares_meander+numsquares_wgnw-meander_length/meander_width + 3*numsquares_taper
    D.info['expected_resistance'] = D.info['num_squares']*EXPECTED_RSQ_WSI
    D.info['wire_width'] = wgnw_width
    D.info['length']= wgnw_length
    return D

def wgnw_with_wb(meander_width = 0.4, num_squares = 1000,
            wgnw_width = 0.1, wgnw_gap = 0.15, wgnw_distance = 0.2, wgnw_length = 100, wg_total_length = None,
            nw_pad_device = None, wb_pad_device = None, wb_gnd_pad_device = None, pad_offset = [500,300], gnd_pad_offset = [0,-300], gnd_pad_orientation = 180,
            ):


    D = Device('wgnw_with_wb')

    wgnw_dev = wgnw(meander_width = meander_width, num_squares = num_squares,
            wgnw_width = wgnw_width, wgnw_gap = wgnw_gap, wgnw_distance = wgnw_distance,
            wgnw_length = wgnw_length, wg_total_length = wg_total_length, nw_pad_device = nw_pad_device)

    wgnw_dev.ports['wiring2'].orientation = gnd_pad_orientation
    wgnw_ref = D << wgnw_dev

    D.add_ref(port2pad(port = wgnw_ref.ports['wiring1'], Pad = wb_pad_device, pad_offset = pad_offset, wiring_layer =lys['m5_wiring']))
    D.add_ref(port2pad(port = wgnw_ref.ports['wiring2'], Pad = wb_gnd_pad_device, pad_offset = gnd_pad_offset, wiring_layer = lys['m5_wiring']))
    D.add_port(wgnw_ref.ports['optical'])

    return D

def optimal_90deg(width = 100.0, num_pts = 15, length_adjust = 1, layer = 0):

    D = Device()
    a = 2*width
    dl = 0.1
    v = 0.1

    # Get points of ideal curve
    v = np.logspace(-length_adjust,length_adjust,num_pts)
    xi = a/2.0*((1+2/np.pi*np.arcsinh(1/v)) + 1j*(1+2/np.pi*np.arcsinh(v)))
    xpts = list(np.real(xi)); ypts = list(np.imag(xi))

    # Add points for the rest of curve
    d = 2*xpts[0] # Farthest point out * 2, rounded to nearest 100
    xpts.append(width); ypts.append(d)
    xpts.append(0); ypts.append(d)
    xpts.append(0); ypts.append(0)
    xpts.append(d); ypts.append(0)
    xpts.append(d); ypts.append(width)
    xpts.append(xpts[0]); ypts.append(ypts[0])

    D.add_polygon([xpts, ypts], layer = layer)

    D.add_port(name = 1, midpoint = [a/4,d], width = a/2, orientation = 90)
    D.add_port(name = 2, midpoint = [d,a/4], width = a/2, orientation = 0)

    return D


###############################################################################
#
# SONIA - CHILES INTEGRATED DEVICES
#
###############################################################################

# device for testing optical absorpiton of hairpin WSi on waveguide

def wgnw_cutback(length_wg = 0,  # length of the single mode waveguide portion
                 length_total = None, # overrides length_wg if specified
                 length_nw = 20, # length of the wgnw
                 width_nw = 0.12,
                 gap_nw = 0.15,
                 length_taper = 20,
                 Grating = None,
                 sm_wg_width = SM_WG_WIDTH,
                 ):

    D = Device()
    width_wg = 0.2*2+2*width_nw+gap_nw # the width of the waveguide the nw sits on
    grating = D << Grating
    pitch_nw = width_nw+gap_nw
    wgnw = D.add_ref(pg.compass(size = (length_nw+0.2, width_wg), layer = lys['wg_deep']))
    nw = D.add_ref(pg.optimal_hairpin(width = width_nw, pitch = pitch_nw, length = length_nw, layer = lys['m2_nw']))
    nw.rotate(180)
    wgnw.xmax = nw.xmax

    if length_total is not None:
        length_wg = length_total - length_nw
    wg = D << pg.compass(size = (length_wg, sm_wg_width), layer = lys['wg_deep'])
    wg.xmax = wgnw.xmin - length_taper
    D << pr.route_basic(port1 = wg.ports['E'], port2 = wgnw.ports['W'], layer = lys['wg_deep'])
    grating.connect(port = 1, destination = wgnw.ports['E'])
    grating.movex(length_taper)
    D << pr.route_basic(port1 = grating.ports[1], port2 = wgnw.ports['E'], layer = lys['wg_deep'])
    D.add_port(port = wg.ports['W'], name = 1)

    return D

# device for testing optical absorpiton of hairpin WSi on waveguide, but no WSi for reference

def wgnw_cutback_reference(length_wg = 0,  # length of the single mode waveguide portion
                           length_total = None, # overrides length_wg if specified
                           length_nw = 20, # length of the wgnw
                           width_nw = 0.12,
                           gap_nw = 0.15,
                           length_taper = 20,
                           Grating = None,
                           sm_wg_width = SM_WG_WIDTH,
                           ):

    D = Device()
    width_wg = 0.2*2+gap_nw+2*width_nw # the width of the waveguide the nw sits on
    grating = D << Grating
    wgnw = D.add_ref(pg.compass(size = (length_nw+0.2, width_wg), layer = lys['wg_deep']))

    if length_total is not None:
        length_wg = length_total - length_nw
    wg = D << pg.compass(size = (length_wg, sm_wg_width), layer = lys['wg_deep'])
    wg.xmax = wgnw.xmin - length_taper
    D << pr.route_basic(port1 = wg.ports['E'], port2 = wgnw.ports['W'], layer = lys['wg_deep'])
    grating.connect(port = 1, destination = wgnw.ports['E'])
    grating.movex(length_taper)
    D << pr.route_basic(port1 = grating.ports[1], port2 = wgnw.ports['E'], layer = lys['wg_deep'])
    D.add_port(port = wg.ports['W'], name = 1)

    return D


###############################################################################
#
# ALEX TAIT DEVICES
#
# Generics from technology definition
# - relies heavily on klayout_technology
# - it will update if SOEN_PDK updates
#
# De facto two-port optical device protocol
# - input: 'wg_in_1' and output: 'wg_out_1'
#
###############################################################################

def optical_connect(*args):
    ''' Moves the first waveguide-like Device onto the output of the second one
        Needs wg_in_1 (in) and wg_out_1 (out).

        Contrast with phidl (0.8.7) connect: it works on Devices and DeviceReferences

        args are in this order (moved_third, moved_second, moved_first, does_not_move)
    '''
    if len(args) < 2:
        raise ValueError('Must specify at least two DeviceReferences to optically connect.')
    def connect_two(moving, target):
        mov_po = moving.ports['wg_in_1']
        tar_po = target.ports['wg_out_1']
        moving.rotate(angle = 180 + tar_po.orientation - mov_po.orientation, center = mov_po.midpoint)
        moving.move(mov_po.midpoint, tar_po.midpoint)
    for iMove in range(len(args)):
        if iMove == 0: continue
        connect_two(args[-iMove-1], args[-iMove])



# @lru_cache would be good
def taper_xsections(wg_src, wg_dest, taper_len, min_tip_width=0.1, keep_layernames=['wg_shallow'],
                    route_basic_options = dict(path_type='sine', width_type='sine')):
    ''' Generic between any two waveguides.
        If any layer is not present in one, it is ignored,
        unless it is in keep_layernames, in which case its width is set to min_tip_width

        TODO: document better
    '''
    # Find all of the relevant layers
    all_layer_names = set()
    for c in wg_src.components + wg_dest.components:
        all_layer_names.add(c.layer)
    D = Device()
    for lName in all_layer_names:
        src_comps = wg_src.get_by_layer(lName)
        dest_comps = wg_dest.get_by_layer(lName)

        # prune it if only one is present, unless it is explicitly in keep_layernames
        if len(src_comps) == 0 or len(dest_comps) == 0:
            if lName in keep_layernames:
                if len(src_comps) == 0:
                    new_comp = dest_comps[0].copy()
                    new_comp.width = min_tip_width
                    src_comps.append(new_comp)
                else:
                    new_comp = src_comps[0].copy()
                    new_comp.width = min_tip_width
                    dest_comps.append(new_comp)
            else:
                continue

        # If there is not the same number, pick off the remainder
        # An alternate behavior would be setting that width to zero
        while len(src_comps) != len(dest_comps):
            if len(src_comps) < len(dest_comps):
                dest_comps.pop(0 if dest_comps[0].offset > 0 else -1)
            else:
                src_comps.pop(0 if src_comps[0].offset > 0 else -1)

        # We have to make temporary ports because of route_basic.
        # TODO: Would be better to separate out the geometry part of that.
        for sco, dco in zip(src_comps, dest_comps):
            src_po = Port(midpoint=(0, sco.offset),
                          width=sco.width, orientation=0)
            dest_po = Port(midpoint=(taper_len, dco.offset),
                           width=dco.width, orientation=180)
            D << pr.route_basic(src_po, dest_po, layer=lys[lName], **route_basic_options)
    D.flatten()  # effectively removes all these ports


    # figure out the port labels
    port_extrema = np.zeros((2, 2))
    for i, wg in enumerate([wg_src, wg_dest]):
        for co in wg.components:
            port_extrema[i, 0] = min(port_extrema[i, 0], co.min)
            port_extrema[i, 1] = max(port_extrema[i, 1], co.max)
    port_ys = np.mean(port_extrema, axis=1)
    port_widths = np.diff(port_extrema, axis=1)[:,0]
    D.add_port(name = 'wg_in_1', midpoint = [0, port_ys[0]],
               width = port_widths[0], orientation = 180)
    D.add_port(name = 'wg_out_1', midpoint = [taper_len, port_ys[1]],
               width = port_widths[1], orientation = 0)
    return D


def device_from_sections(sections):
    ''' Place given sections into a parent device and connect wg_in_1's to wg_out_1's
        The order matters.

        Args:
            sections (list[Devices]): each having a wg_in_1 and wg_out_1 port
    '''
    DEV = Device()
    for iSec, sec in enumerate(sections):
        this_sec = DEV << sec
        if iSec > 0:
            optical_connect(this_sec, last_sec)
        last_sec = this_sec
        if iSec == 0:
            DEV.add_port(port = this_sec.ports['wg_in_1'])
        if iSec == len(sections) - 1:
            DEV.add_port(port = this_sec.ports['wg_out_1'])
    return DEV


def tapered_wg_straight(length, ends_wg, middle_wg, taper_len, cap_both=True, keep_layernames=['wg_shallow'], width_type='straight'):
    ''' Straight waveguide with tapers on the ends.

        Args:
            length (float): total length between port faces
            ends_wg (WGXSection): port faces have this type
            middle_wg (WGXSection): central segment will have this cross-section
            taper_len (float): for each taper
            cap_both (bool): put on the second taper on the end or not
            width_type (str): feeds through eventually to phidl route_basic
    '''
    long_len = length - taper_len * (2 if cap_both else 1)
    if long_len < 0:
        raise ValueError('Tapered waveguide not long enough to accomadate tapers')
    trans1 = taper_xsections(ends_wg, middle_wg, taper_len, route_basic_options=dict(width_type=width_type), keep_layernames=keep_layernames)
    long_sect = middle_wg.cell_straight(long_len)
    sections = [trans1, long_sect]
    if cap_both:
        trans2 = taper_xsections(middle_wg, ends_wg, taper_len, route_basic_options=dict(width_type=width_type), keep_layernames=keep_layernames)
        sections.append(trans2)
    return device_from_sections(sections)


def route_wg_straight(length, wg_name=None):
    ''' A straight waveguide with lowest loss available on the platform.
        Incorporates propagation and taper loss data to determine whether to taper out to a long haul waveguide.
        This is used just for routing because the waveguide type is not guaranteed (unless explicit)

        As of now, wg_name can also be a WGXSection

        Args:
            length (float): length of waveguide
            wg_name (str, None, WGXSection): None means it can use a mixture, depending on which is lowest loss
    '''
    if wg_name is not None:
        if isinstance(wg_name, str):
            wg_xsection = tech.waveguides(wg_name)
        else:
            wg_xsection = wg_name
        return wg_xsection.cell_straight(length)
    else:
        # Figure out which to use
        trans_obj = tech.transitions('Strip to Long haul')
        transition_loss = trans_obj.loss
        long_len = length - 2 * trans_obj.length
        strip_option_loss = length / 1e4 * tech.waveguides('Strip').loss
        tapered_option_loss = long_len / 1e4 * tech.waveguides('Long haul').loss + 2 * transition_loss

        if strip_option_loss <= tapered_option_loss:
            return tech.waveguides('Strip').cell_straight(length)
        else:
            return tapered_wg_straight(length, tech.waveguides('Strip'), tech.waveguides('Long haul'), trans_obj.length)


def route_wg_points(ptlist, wg_name=None, radius=None):
    ''' Puts a connected sequence of sections of this waveguide.

        The wg_name argument can be a string in the technology (will be deprecated)

        TODO: the argument name is misleading. Should be wg_xsection
        or a WGXSection object.
    '''
    # Compute angles and lengths (not accounting for arcs) of segments
    points = np.asarray(ptlist)
    dxdy = np.diff(points, axis=0)
    angles = (np.arctan2(dxdy[:,1], dxdy[:,0])).tolist()
    angles = np.array(angles + [angles[-1]]) * 180 / np.pi
    turns = ((angles[1:] - angles[:-1]) + 180) % 360 - 180
    if any(abs(turns) > 165):
        print('Warning: very sharp turns')

    lengths = np.sqrt(np.sum(dxdy ** 2, axis=1))
    nz = np.nonzero(lengths)
    lengths = lengths[nz]  # Check for repeated points
    turns = turns[nz]
    angles = angles[nz]

    # Straight WG type is determined above, but now we need bends
    if wg_name is not None:
        if isinstance(wg_name, str):
            bend_wg = tech.waveguides(wg_name)
        else:
            bend_wg = wg_name
    else:
        bend_wg = tech.waveguides('Strip')

    if radius is None:
        radius = bend_wg.radius

    WG = Device()
    next_point = points[0]
    for iSegment in range(len(lengths)):
        # Adjust straight section length relative to distance between anchor points
        adj_len = lengths[iSegment]
        if iSegment > 0:
            adj_len -= radius * abs(np.tan(turns[iSegment-1] / 2 * np.pi / 180))
        if iSegment < len(lengths) - 1:
            adj_len -= radius * abs(np.tan(turns[iSegment] / 2 * np.pi / 180))
        if adj_len < 0:
            raise ValueError('Length was negative. Points are too close together or turns are too sharp')

        # Start a sub-cell for the straight section plus one bend
        straight = route_wg_straight(adj_len, wg_name)
        if iSegment < len(lengths)-1:
            bent = bend_wg.cell_bend(theta=turns[iSegment], radius=radius)
            optical_connect(bent, straight)
            bent_ref = straight << bent
        straight.rotate(angles[iSegment])
        # Insert into top WG and connect to the previous bend
        straight_ref = WG << straight
        if iSegment > 0:
            optical_connect(straight, prev_bent_ref)
        # Note that this has to happen after everything is placed so that the port coordinates end up in the top WG frame
        if iSegment < len(lengths) - 1:
            prev_bent_ref = bent_ref
        # Handle ports
        if iSegment == 0:
            WG.add_port(port=straight.ports['wg_in_1'])
        elif iSegment == len(lengths) - 1:
            WG.add_port(port=straight.ports['wg_out_1'])
    WG.move(ptlist[0]).flatten()
    return WG


def route_ports(port1, port2, wg_name=None):
    ''' Puts a waveguide with one central joint at the interection of the port lines '''
    import routelib_dev.vector_math as vecmath
    ray1 = vecmath.Ray.from_port(port1)
    ray2 = vecmath.Ray.from_port(port2)
    joint = vecmath.ray_intersect(ray1, ray2).as_array()
    points = np.array([port1.midpoint, joint, port2.midpoint])
    return route_wg_points(points, wg_name)
    # if wg_name != 'Strip':
    #     raise ValueError('Currently only Strip supported for routing between ports (will change soon)')
    # return pr.route_manhattan(port1, port2, bendType='circular', layer=lys['wg_deep'], radius=tech.waveguides('Strip').radius)


def bragg_grating(length, wg_xsection=None, period=0.2, duty=0.5, dWidth=.05, perturbation='centered', cap_to='length', apodization=None, phase=np.pi):
    '''
        This expects waveguides of the NOT converted type; only one partial etch component.
        The wg_xsection's 'wg_deep' layer will be modulated **inwards**.
        You can change this behavior by enlarging the one that is passed in (or maybe this function needs a new argument)

        Args:
            length (float): total length (approximate or exact depends on cap)
            wg_xsection (float): base X section that must contain a 'wg_deep' component to be modulated
            period (float): grating period
            duty (float): ratio of wide segments
            dWidth (float): difference between min and max full widths.
            perturbation (str): 'inwards', 'outwards', or 'centered' (default)
            cap_to (str): the way to end the grating. It always starts big. Can be:
                * 'length': Lengths are guaranteed. A stub is made with the width of the regular waveguide
                * 'full': Ends on full period, meaning the end will be skinny
                * 'half': Ends on half period, meaning the end will be fat
            apodization (float, None): length of regions on the ends that apodize linearly
            phase (float): skew between left and right sides of gratings for precisely controlling coupling coefficient
    '''
    assert cap_to in ['length', 'half', 'full']
    if wg_xsection is None:
        wg_xsection = tech.waveguides('Strip')
    phase = phase % (2 * np.pi)
    phase_dependent_offset = (phase / np.pi - 1) * period / 2

    assert duty < 1
    WGBG = Device()
    nPeriods = int(np.floor(length / period))

    # Place everything except the waveguide ridge
    temp_wg = wg_xsection.copy()
    temp_wg.unconvert_from_masklike()
    portwidth = 0
    for comp in temp_wg.components:
        portwidth = max(portwidth, comp.width)
    ridge_component = temp_wg.get_by_layer('wg_deep')[0]
    assert ridge_component.offset == 0  # for now, we assume this
    temp_wg.components.remove(ridge_component)
    if len(temp_wg.components) > 0:
        temp_wg_ref = WGBG << temp_wg.cell_straight(length)

    # Some calculations on the sections
    if perturbation == 'inwards':
        min_half_width = (ridge_component.width - dWidth) / 2
    elif perturbation == 'outwards':
        min_half_width = (ridge_component.width) / 2
    elif perturbation == 'centered':
        min_half_width = (ridge_component.width - dWidth / 2) / 2
    else:
        raise ValueError('Invalid perturbation direction {}. Must be inwards, outwards, or centered.'.format(perturbation))
    max_half_width = min_half_width + dWidth / 2

    # The section generators
    @lru_cache(maxsize=2)
    def r_big(dWfrac=1.):
        half_width = (dWfrac * max_half_width +
                      (1 - dWfrac) * ridge_component.width / 2)
        return pg.rectangle((period * duty, half_width), layer=lys[ridge_component.layer])
    @lru_cache(maxsize=2)
    def r_small(dWfrac=1.):
        half_width = (dWfrac * min_half_width +
                      (1 - dWfrac) * ridge_component.width / 2)
        return pg.rectangle((period * (1 - duty), half_width), layer=lys[ridge_component.layer])

    def fractional_dW(x):
        ''' A function of dW(x) / dWidth vs. x from 0 to 1 '''
        if apodization is None or apodization == 0:
            return 1.0
        else:
            apodization_fraction = apodization / length
            if x < apodization_fraction:
                return x / apodization_fraction
            elif 1 - x < apodization_fraction:
                return (1 - x) / apodization_fraction
            else:
                return 1.0

    # Place it
    # Note this is not efficient, but it is flexible in case we later want to chirp
    ridge = Device()
    for iPeriod in range(nPeriods):
        r_big_left = ridge << r_big(fractional_dW(iPeriod / nPeriods))
        r_big_left.ymin = 0
        r_big_left.xmin = iPeriod * period
        r_big_right = ridge << r_big(fractional_dW(iPeriod / nPeriods))
        r_big_right.ymax = 0
        r_big_right.xmin = iPeriod * period + phase_dependent_offset
        r_small_left = ridge << r_small(fractional_dW(iPeriod / nPeriods))
        r_small_left.ymin = 0
        r_small_left.xmin = r_big_left.xmax
        r_small_right = ridge << r_small(fractional_dW(iPeriod / nPeriods))
        r_small_right.ymax = 0
        r_small_right.xmin = r_big_right.xmax
    if cap_to == 'half':
        r_big_left = ridge << r_big(fractional_dW(1))
        r_big_left.ymin = 0
        r_big_left.xmin = nPeriods * period
        r_big_right = ridge << r_big(fractional_dW(1))
        r_big_right.ymax = 0
        r_big_right.xmin = nPeriods * period + phase_dependent_offset

    # Fill out the right side if its phase is off
    # TODO: these are not correct if you are using apodization and phase offset together
    if phase_dependent_offset != 0:
        if cap_to == 'half':
            r_fill = ridge << pg.rectangle((r_big_left.xmax - r_big_right.xmax, max_half_width), layer=lys[ridge_component.layer])
            r_fill.ymax = 0
            r_fill.xmin = r_big_right.xmax
        else:
            r_fill = ridge << pg.rectangle((r_small_left.xmax - r_small_right.xmax, min_half_width), layer=lys[ridge_component.layer])
            r_fill.ymax = 0
            r_fill.xmin = r_small_right.xmax
    ridge.flatten()
    WGBG << ridge

    # Capping
    if cap_to == 'length':
        remainder_len = length - period * nPeriods
        if remainder_len != 0:
            r_remainder = WGBG << pg.rectangle((remainder_len, ridge_component.width), layer=lys[ridge_component.layer])
            r_remainder.xmin = ridge.xmax
            r_remainder.y = ridge.y

    # Handle ports
    WGBG.add_port(name='wg_in_1', midpoint=[0,0], width=portwidth, orientation=180)
    WGBG.add_port(name='wg_out_1', midpoint=[WGBG.xmax, 0], width=portwidth, orientation=0)

    return WGBG


def wg_terminator(wg_to_terminate=None, taper_len=20, doped=True):
    ''' Taper all the WG components to nothing and put on some dopants
    '''
    if wg_to_terminate is None:
        wg_to_terminate = tech.waveguides('Strip')
    terminus = tech.properties.WGXSection(components=[])
    taper_section = taper_xsections(wg_to_terminate, terminus, taper_len, min_tip_width=.05,
                                    keep_layernames=['wg_deep', 'wg_shallow'])
    if doped:
        dopant_section = pg.rectangle((taper_len, 3), layer=lys['dp_n+'])
        dopant_section.x = taper_section.xmax
        dopant_section.y = taper_section.y
        taper_section << dopant_section
    return taper_section


def loop_mirror_terminator(y_splitter=None):
    ''' A loop mirror (or Sagnac interferometer) consisting of a splitter with connected outputs.
        Performance should be pretty broad band, determined by the MMI.

        The only port is 'wg_in_1'
    '''
    if y_splitter is None:
        y_splitter = mmi1x2(gap_mmi=.5)
    D = Device()
    split = D << y_splitter
    # split = D << mmi1x2(length_mmi=5.75,width_mmi=2.05,gap_mmi=gap_mmi, layer=lys['wg_deep'])
    radius = tech.waveguides('Strip').radius
    width = tech.waveguides('Strip').components[0].width
    arc = D << pr._arc(radius=radius, width=width, layer=lys['wg_deep'],
                   theta=180, start_angle=-90, angle_resolution=.25)
    # Place and connect them
    arc.y = split.y
    arc.xmin = split.xmax + 10
    D << pr.route_basic(split.ports[2], arc.ports[2], layer=lys['wg_deep'])
    D << pr.route_basic(split.ports[3], arc.ports[1], layer=lys['wg_deep'])
    D.flatten()
    D.add_port(port=split.ports[1], name='wg_in_1')
    return D


def dreidel_reflector(wg_width=0.35, width_mmi=1.55, length_mmi=2.8):
    D = Device()
    mmi_section = D << pg.compass((length_mmi, width_mmi), layer=lys['wg_deep'])
    wg_approach = D << pg.compass((2, wg_width), layer=lys['wg_deep'])
    wg_approach.connect('E', mmi_section.ports['W'])
    tip_angle = 90
    triangle_length = width_mmi / 2 / (np.tan(np.pi / 180 * tip_angle / 2))
    pts = [(mmi_section.xmax, mmi_section.ymax),
           (mmi_section.xmax + triangle_length, mmi_section.y),
           (mmi_section.xmax, mmi_section.ymin)]
    D.add_polygon(pts, layer=lys['wg_deep'])
    D.add_port('wg_in_1', port=wg_approach.ports['W'])
    return D


def SOA(gain_len, taper_len=10, mid_wid=0):
    ''' A semiconductor optical amplifier that is just Strip tapered into an Active section that can be widened.
        The gain_len does not count tapers

        Args:
            gain_len (float): active section
            taper_len (float): tapers between
            mid_wid (float): additional widening applied to the gain section
    '''
    wider_xs = tech.waveguides('Active').copy()
    for comp in wider_xs.components:
        if comp.offset == 0:
            comp.width += mid_wid
        else:
            comp.offset += np.sign(comp.offset) * mid_wid / 2
    return tapered_wg_straight(gain_len + 2 * taper_len, tech.waveguides('Strip'), wider_xs, taper_len=taper_len, width_type='sine')


def _DBR_sections(wbg1_len=40, wbg2_len=10, etalon_len=40, taper_len=10, wbg_params=None):
    ''' Make a list of sections typically used for a DBR laser
    '''
    if wbg_params is None:
        wbg_params = dict()
    TAPER1 = taper_xsections(tech.waveguides('Strip'), tech.waveguides('Active'), taper_len, keep_layernames=['wg_shallow'])
    BG1 = bragg_grating(wbg1_len, tech.waveguides('Active'), **wbg_params)
    ETALON = tech.waveguides('Active').cell_straight(etalon_len)
    BG2 = bragg_grating(wbg2_len, tech.waveguides('Active'), **wbg_params)
    TAPER2 = taper_xsections(tech.waveguides('Active'), tech.waveguides('Strip'), taper_len, keep_layernames=['wg_shallow'])
    sections = [TAPER1, BG1, ETALON, BG2, TAPER2]
    return sections


def DBR_laser(wbg1_len=40, wbg2_len=10, etalon_len=40, taper_len=10, wbg_params=None):
    ''' Standard DBR that feeds through sections to the section concatenator '''
    sections = _DBR_sections(wbg1_len, wbg2_len, etalon_len, taper_len, wbg_params)
    return device_from_sections(sections)


def widened_DBR(mid_wid, midtaper_len=10, dbr_params=None):
    ''' Modified DBR where the central gain section is widened by some amount.
        The variable etalon_len does not count the length of the tapers.

        Args:
            mid_wid (float): additional width
            midtaper_len (float): length of tapers into widened section
            dbr_params (float): passed through to _DBR_sections
    '''
    # Calculate a widened cross section
    wider_xs = tech.waveguides('Active').copy()
    for comp in wider_xs.components:
        if comp.offset == 0:
            comp.width += mid_wid
        else:
            comp.offset += np.sign(comp.offset) * mid_wid / 2
    # Get the sections of the basic dbr
    if dbr_params is None:
        dbr_params = dict()
    sections = _DBR_sections(**dbr_params)
    mid_len = sections[2].xsize + 2 * midtaper_len
    # Insert the new section with tapers
    new_mid_section = tapered_wg_straight(mid_len, tech.waveguides('Active'), wider_xs, midtaper_len, width_type='sine')
    sections[2] = new_mid_section
    return device_from_sections(sections)


# def DFB_laser(*args, **kwargs):
#     ''' This is where it would REALLY help to have devices/PCells as classes '''
#     try:
#         bg_period = kwargs['wbg_params']['period']
#     except KeyError:
#         bg_period = BraggGrating.defaults.period  # right HERE!
#     kwargs['etalon_len'] = bg_period / 4
#     return DBR_laser(*args, **kwargs)


arb_via_dopant_inclusion = 0.5
arb_via_m4_inclusion = 0.8
def autometal_by_layer(D, layers_to_contact=['dp_n+', 'dp_p+']):
    ''' Find the dopants automatically and put metal on them
        Modifies device in place.
        Returns a dictionary: keys are layernames, values are lists of polygons for the metal
    '''
    pad_polys = dict()
    for heavy_dopant in layers_to_contact:
        contact = pg.extract(D, layers=[lys[heavy_dopant]])
        pad_polys[heavy_dopant] = list()
        via = pg.offset(contact, -arb_via_dopant_inclusion, layer=lys['v5'])
        D << via
        pad = pg.offset(via, arb_via_m4_inclusion, layer=lys['m4_ledpad'])
        D << pad
        pad_polys[heavy_dopant].append(pad)
    return pad_polys


def silvered_LM(gap=.2, tap_type='LM', loop_mirror_device=None):
    ''' We want to test if dobule LM gives us double efficiency.
        We also generally want to know the coupling vs. gap when used in LM lasers.
        tap_type in ['LM', 'straight', 'round']
    '''
    D = Device()
    if loop_mirror_device is None:
        loop_mirror_device = loop_mirror_terminator()


    loop1 = D << loop_mirror_device
    if tap_type == 'LM':
        loop2 = D << loop_mirror_device
        loop2.rotate(180)
    elif tap_type == 'straight':
        STAP = Device()
        straight = STAP << tech.waveguides('Strip').cell_straight(loop_mirror_device.ysize)
        term = STAP << wg_terminator(doped=False)
        turn1 = STAP << tech.waveguides('Strip').cell_bend(90)
        turn2 = STAP << tech.waveguides('Strip').cell_bend(90)
        optical_connect(turn2, straight, turn1)
        term.connect('wg_in_1', turn1.ports['wg_in_1'])
        STAP.add_port('wg_in_1', port=turn2.ports['wg_out_1'])
        loop2 = D << STAP
        loop2.rotate(180)
    elif tap_type == 'round':
        RTAP = Device()
        term = RTAP << wg_terminator(doped=False)
        turn = RTAP << tech.waveguides('Strip').cell_bend(180)
        term.connect('wg_in_1', turn.ports['wg_in_1'])
        RTAP.add_port('wg_in_1', port=turn.ports['wg_out_1'])
        loop2 = D << RTAP
        loop2.rotate(180)
    loop2.y = loop1.y
    loop2.xmin = loop1.xmax + gap
    D.add_port('wg_in_1', port=loop1.ports['wg_in_1'])
    D.add_port('wg_out_1', port=loop2.ports['wg_in_1'])
    return D


def loopMirror_laser(straight_len=40, taper_len=10, coupling_gap=.2, active_xsection=None, metalize=True):
    ''' default xsection is from the tech.waveguides, but perhaps you wanted to modify it yourself.
        straight_len is here defined as between the tapers. should
    '''
    # this is just randomly floating here. It has to be invoked after "contact_pads" is declared
    from nc_constant_devices import WB_PAD_EPOXY_DEV, WB_PAD_DEV, WB_GND_PAD_DEV, NW_PAD_DEV

    if active_xsection is None:
        active_xsection = tech.waveguides('Active')


    # Make the taper-mirror end cap
    CAPS = []
    for iCap in range(2):
        CAP = Device()
        taper = CAP << taper_xsections(active_xsection, tech.waveguides('Strip'),
                                       taper_len, keep_layernames=['wg_shallow'])#, 'dp_p+', 'dp_p', 'dp_e', 'dp_n', 'dp_n+'])
        if iCap == 0:
            mirror = CAP << loop_mirror_terminator()
            optical_connect(mirror, taper)
            CAP.add_port('wg_out_1', port=taper.ports['wg_in_1'])
            CAP.rotate(180)
        else:
            mirror = CAP << silvered_LM(gap=coupling_gap, tap_type='LM')
            optical_connect(mirror, taper)
            CAP.add_port('wg_out_1', port=mirror.ports['wg_out_1'])
        CAP.add_port('wg_in_1', port=taper.ports['wg_in_1'])
        CAPS.append(CAP)

    # Put together the cavity
    Cavity = Device()
    straight_ref = Cavity << active_xsection.cell_straight(straight_len)
    left_cap = Cavity << CAPS[0]
    right_cap = Cavity << CAPS[1]
    optical_connect(right_cap, straight_ref, left_cap)
    Cavity.add_port('wg_out_1', port=right_cap.ports['wg_out_1'])
    Cavity.flatten()
    if not metalize: return Cavity

    # and add it to the top device
    D = Device()
    cavity_ref = D << Cavity
    D.add_port('wg_out_1', port=cavity_ref.ports['wg_out_1'])


    # Metalization: device via level
    pad_polys = autometal_by_layer(D)

    # Metalization: small level (L1)
    L1_breakout_type = [None, 'tapered', 'linear'][1]
    if L1_breakout_type == None:
        pass
    elif L1_breakout_type == 'tapered':
        L1_pad_device = pg.compass(size=(20, 20), layer=lys['m4_ledpad'])
    elif L1_breakout_type == 'linear':
        L1_pad_device = pg.compass(size=(pad_polys['dp_n+'][0].xsize, 20), layer=lys['m4_ledpad'])

    alongside_ridge_n = pad_polys['dp_n+'][0]
    alongside_ridge_n.add_port(name = 'N', midpoint = [alongside_ridge_n.x, alongside_ridge_n.ymax],  width = alongside_ridge_n.xsize, orientation = 90)
    alongside_ridge_p = pad_polys['dp_p+'][0]
    alongside_ridge_p.add_port(name = 'S', midpoint = [alongside_ridge_p.x, alongside_ridge_p.ymin],  width = alongside_ridge_p.xsize, orientation = -90)

    if L1_breakout_type is not None:
        L1rect_n = D << L1_pad_device
        L1rect_n.ymin = alongside_ridge_n.ymax + 20
        L1rect_n.x = alongside_ridge_n.x
        D << pr.route_basic(L1rect_n.ports['S'], alongside_ridge_n.ports['N'], width_type='sine', layer=lys['m4_ledpad'])
        L2via_n = D << pg.offset(L1rect_n, -arb_via_m4_inclusion, layer=lys['v5'])
        L2pad_n = D << pg.offset(L1rect_n, arb_via_m4_inclusion, layer=lys['m5_wiring'])
        D.add_port(name = 'el_dn', midpoint = [L2pad_n.x, L2pad_n.ymax],  width = L2pad_n.xsize, orientation = 90)

        L1rect_p = D << L1_pad_device
        L1rect_p.ymax = alongside_ridge_p.ymin - 20
        L1rect_p.x = alongside_ridge_p.x
        D << pr.route_basic(L1rect_p.ports['N'], alongside_ridge_p.ports['S'], width_type='sine', layer=lys['m4_ledpad'])
        L2via_p = D << pg.offset(L1rect_p, -arb_via_m4_inclusion, layer=lys['v5'])
        L2pad_p = D << pg.offset(L1rect_p, arb_via_m4_inclusion, layer=lys['m5_wiring'])
        D.add_port(name = 'el_dp', midpoint = [L2pad_p.x, L2pad_p.ymin],  width = L2pad_p.xsize, orientation = -90)

    # Metalization: medium level (L2)
    L2_breakout_type = [None, 'vertical', 'leftwards'][0]
    if L2_breakout_type is None:
        pass
    elif L2_breakout_type == 'vertical':
        D << port2pad(port = D.ports['nS'], Pad = WB_GND_PAD_DEV, pad_offset = [WB_GND_PAD_DEV.size[0]+50,0], wiring_layer = lys['m5_wiring'])

        D << port2pad(port = D.ports['pN'], Pad = WB_PAD_EPOXY_DEV, pad_offset = [300,100], wiring_layer = lys['m5_wiring'])
    elif L2_breakout_type == 'leftwards':
        for L1rect in [L1rect_n, L1rect_p]:
            L2via = D << pg.offset(L1rect, -arb_via_m4_inclusion, layer=lys['v5'])
            L2pad = D << pg.offset(L1rect, arb_via_m4_inclusion, layer=lys['m5_wiring'])
            L2wire_dev = pg.rectangle(size=(L2pad.xmin - cavity_ref.xmin, L2pad.ysize), layer=lys['m5_wiring'])
            L2wire_dev.add_port(name = 'W', midpoint = [L2wire_dev.xmin, L2wire_dev.y],  width = L2wire_dev.ysize, orientation = 180)
            L2wire = D << L2wire_dev
            L2wire.y = L2pad.y
            L2wire.xmax = L2pad.xmin

    D.flatten()
    return D


def get_wg_under_nw_xsect(wgnw_width, wgnw_gap, wgnw_distance, **kwargs):
    ''' **kwargs are just used to trash additional geometry bits '''
    wg_under_nanowire_width = 2 * wgnw_width + wgnw_gap + wgnw_distance * 2
    strip_component = tech.WGXSectionComponent(layer='wg_deep',
                                               width=wg_under_nanowire_width,
                                               offset=0.0)
    wg_under_nanowire_xsect = tech.WGXSection(components=[strip_component])
    return wg_under_nanowire_xsect


def HiDRA(n_stages=5, pad_separation=50, snspd_standoff=200,
          broken=False, wg_xsection=None, tap_dev=None, mirrored=False, **nanotap_kwargs):
    # this is just randomly floating here. It has to be invoked after "contact_pads" is declared
    from nc_constant_devices import WB_PAD_EPOXY_DEV, WB_PAD_DEV, WB_GND_PAD_DEV, NW_PAD_DEV

    if wg_xsection is None:
        wg_xsection = tech.waveguides('Strip')  # or strip1550
    wg_width = wg_xsection.get_by_layer('wg_deep')[0].width
    pad_pitch = WB_PAD_EPOXY_DEV.xsize + pad_separation
    break_gap = 50 if broken else 0
    snspd_standoff -= break_gap

    if tap_dev is None:
        typical_nanotap = dict(bend_radius=wg_xsection.radius,
                             s_offset=4, s_length=25,
                             wg_width=wg_width)
        typical_nanotap.update(nanotap_kwargs)
        tap_dev = nano_tap_s(**typical_nanotap)
        # FINAL REVIEW: need to look at gap/length of broken 1550 one

    # All of the subcells that will be reused
    standoff_dev = wg_xsection.cell_straight(snspd_standoff)

    if broken:
        new_standoff_dev = Device()
        standoff_ref = new_standoff_dev << standoff_dev
        turn_away = new_standoff_dev << wg_xsection.cell_bend(-180, radius=break_gap/2)
        optical_connect(turn_away, standoff_ref)
        taper_away = new_standoff_dev << wg_terminator(wg_to_terminate=wg_xsection, doped=False)
        optical_connect(taper_away, turn_away)
        standoff_dev = new_standoff_dev
        standoff_dev.add_port(name='wg_in_1', port=standoff_ref.ports['wg_in_1'])  # this is a fake not connected to anything
        standoff_dev.add_port(name='wg_out_1', port=standoff_ref.ports['wg_out_1'])  # this is a fake not connected to anything
    pad_kwargs = dict(nw_pad_device = NW_PAD_DEV,
                      wb_pad_device = WB_PAD_EPOXY_DEV, wb_gnd_pad_device = WB_GND_PAD_DEV,
                      pad_offset = [300,0],
                      gnd_pad_offset = [100,0], gnd_pad_orientation = -90)
    nw_geom_kwargs = dict(num_squares=6000, wgnw_length=200,
                          wgnw_width=0.14, wgnw_gap=0.15, wgnw_distance=0.2)
    snspd_dev = wgnw_with_wb(**nw_geom_kwargs, **pad_kwargs)

    # Top cell devices
    wg_under_nanowire_xsect = get_wg_under_nw_xsect(**nw_geom_kwargs)
    taper_dev = taper_xsections(wg_xsection, wg_under_nanowire_xsect,
                                taper_len=50,
                                route_basic_options=dict(path_type='straight', width_type='straight'))
    # Top cell construction
    D = Device('HiDRA')
    input_wg = D << wg_xsection.cell_straight(WB_PAD_EPOXY_DEV.xsize / 2)
    prev_tap = None
    for iTap in range(n_stages):
        # The commented stuff has to do with leaving off the tap on the last one
        # if iTap == n_stages - 1:
        #     fake_tap_dev = wg_xsection.cell_bend(-90)
        #     fake_tap_dev.add_port(1, port=fake_tap_dev.ports['wg_in_1'])
        #     fake_tap_dev.add_port(2, port=fake_tap_dev.ports['wg_out_1'])
        #     this_tap = D << fake_tap_dev
        #     turn_ref = D << wg_xsection.cell_straight(1)
        # else:
        this_tap = D << tap_dev
        this_tap.reflect()
        turn_ref = D << wg_xsection.cell_bend(-90)
        if iTap == 0:
            this_tap.connect(1, input_wg.ports['wg_out_1'])
        else:
            this_tap.connect(1, prev_tap.ports[2], overlap=-pad_pitch)
            D << pr.route_basic(this_tap.ports[1], prev_tap.ports[3], layer=lys['wg_deep'])
        if iTap == n_stages - 1:
            terminus = D << wg_terminator(wg_to_terminate=wg_xsection, taper_len=40, doped=False)
            terminus.connect('wg_in_1', this_tap.ports[3])

        turn_ref.connect('wg_in_1', this_tap.ports[2])
        standoff_ref = D << standoff_dev
        standoff_ref.connect('wg_in_1', turn_ref.ports['wg_out_1'])
        taper_ref = D << taper_dev
        taper_ref.connect('wg_in_1', standoff_ref.ports['wg_out_1'], overlap=-break_gap)
        snspd_ref = D << snspd_dev
        snspd_ref.connect('optical', taper_ref.ports['wg_out_1'])
        prev_tap = this_tap

    D.add_port(port=input_wg.ports['wg_in_1'])
    D.add_port(name='optical', port=D.ports['wg_in_1'])  # redundant in order to work with bs_tree
    if mirrored:
        D.reflect()
    return D


def insert_snspd_pair(parent_device, port_up, port_down, pad_separation=200, snspd_standoff=None, wg_xsection=None):
    # currently this only works if both orientations are 0 or if port_up is 90 and port_down is -90
    if abs((port_up.orientation - port_down.orientation) % 360) == 0:
        route_fun = pr.route_manhattan180
    elif abs((port_up.orientation - port_down.orientation) % 360) == 180:
        route_fun = pr.route_manhattan90
    else:
        raise ValueError('Ports must either be parallel or antiparallel')

    # this is just randomly floating here. It has to be invoked after "contact_pads" is declared
    from nc_constant_devices import WB_PAD_EPOXY_DEV, WB_GND_PAD_DEV, NW_PAD_DEV
    if wg_xsection is None:
        wg_xsection = tech.waveguides('Strip')
    wg_width = wg_xsection.get_by_layer('wg_deep')[0].width

    pad_pitch = WB_PAD_EPOXY_DEV.xsize + pad_separation
    pad_kwargs = dict(nw_pad_device = NW_PAD_DEV,
                      wb_pad_device = WB_PAD_EPOXY_DEV, wb_gnd_pad_device = WB_GND_PAD_DEV,
                      pad_offset = [300,0],
                      gnd_pad_offset = [100,0], gnd_pad_orientation = -90)
    nw_geom_kwargs = dict(num_squares=6000, wgnw_length = 200,
                          wgnw_width = 0.14, wgnw_gap = 0.15, wgnw_distance = 0.2)
    snspd_dev = wgnw_with_wb(**nw_geom_kwargs, **pad_kwargs)

    wg_under_nanowire_xsect = get_wg_under_nw_xsect(**nw_geom_kwargs)
    taper_dev = taper_xsections(wg_xsection, wg_under_nanowire_xsect,
                                taper_len=20,
                                route_basic_options=dict(path_type='straight', width_type='straight'))

    route_radius = 20
    if snspd_standoff is None:
        snspd_standoff = route_radius * 4 + 2

    for port, side_sign in zip([port_up, port_down], [1, -1]):
        taper = parent_device << taper_dev
        taper.move(taper.ports['wg_in_1'].midpoint, port.midpoint)
        taper.move((snspd_standoff, side_sign * pad_pitch / 2))
        parent_device << route_fun(taper.ports['wg_in_1'], port,
                                   bendType='gradual', radius=route_radius,
                                   layer=lys['wg_deep'])
        snspd = parent_device << snspd_dev
        if side_sign == -1:
            snspd.reflect()
        snspd.connect('optical', taper.ports['wg_out_1'])


def HongOuMandel(pad_separation=200, snspd_standoff=None, wg_xsection=None, in_5050=None, **coupler_kwargs):
    if wg_xsection is None:
        wg_xsection = tech.waveguides('Strip')  # or strip1550
    wg_width = wg_xsection.get_by_layer('wg_deep')[0].width

    D = Device('HOM')

    if in_5050 is None:
        in_5050 = adiabatic5050_coupler(out_wg_width=wg_width, **coupler_kwargs)
    in_5050 = D << in_5050
    insert_snspd_pair(D, in_5050.ports['wg_out_1'], in_5050.ports['wg_out_2'],
                      pad_separation=pad_separation, snspd_standoff=snspd_standoff,
                      wg_xsection=wg_xsection)

    D.add_port(port=in_5050.ports['wg_in_1'])
    try:
        D.add_port(port=in_5050.ports['wg_in_2'])
    except KeyError:
        # This is a 1x2 not 2x2
        pass
    return D


def HanburyBrownTwiss(pad_separation=200, snspd_standoff=400, wg_xsection=None):
    if wg_xsection is None:
        wg_xsection = tech.waveguides('Strip')  # or strip1550
    wg_width = wg_xsection.get_by_layer('wg_deep')[0].width

    D = Device('HBT')

    in_yjunction = D << y_junction_curved(wg_width=wg_width)
    insert_snspd_pair(D, in_yjunction.ports[2], in_yjunction.ports[3],
                      pad_separation=pad_separation, snspd_standoff=snspd_standoff, wg_xsection=wg_xsection)

    D.add_port(port=in_yjunction.ports[1], name='wg_in_1')
    return D


def waveguide_hairpin(inout_spacing=None, wg_xsection=None):
    if wg_xsection is None:
        wg_xsection = tech.waveguides('Strip')
    if inout_spacing is None:
        inout_spacing = wg_xsection.radius / 2
    if inout_spacing*2 < wg_xsection.radius:
        raise ValueError('Pitch is too tight for our hairpin algorithm')
    hairpin = Device('hairpin')
    large_piece = hairpin << wg_xsection.cell_euler_bend(theta=180, radius=inout_spacing * 2, angular_coverage=45)
    target_y = abs(large_piece.ports['wg_in_1'].y - large_piece.ports['wg_out_1'].y) - inout_spacing

    small_piece = hairpin << wg_xsection.cell_s_bend_by_lateral_offset(target_y, radius=inout_spacing * 2)
    optical_connect(large_piece, small_piece)

    diff_x = abs(small_piece.ports['wg_in_1'].x - large_piece.ports['wg_out_1'].x)
    difference_wg = hairpin << wg_xsection.cell_straight(diff_x)
    optical_connect(difference_wg, large_piece)

    hairpin.add_port('wg_in_1', port=small_piece.ports['wg_in_1'])
    hairpin.add_port('wg_out_1', port=difference_wg.ports['wg_out_1'])
    hairpin.flatten()

    return hairpin


def waveguide_hairpin_tweaked(inout_spacing=None, wg_xsection=None):
    if wg_xsection is None:
        wg_xsection = tech.waveguides('Strip')
    if inout_spacing is None:
        inout_spacing = wg_xsection.radius / 2
    # if inout_spacing*2 < wg_xsection.radius:
    #     raise ValueError('Pitch is too tight for our hairpin algorithm')
    hairpin = Device('hairpin')
    tight_corner = hairpin << wg_xsection.cell_euler_bend(theta=-125)
    tight_corner.rotate(90)

    large_turn = hairpin << wg_xsection.cell_euler_bend(theta=215)
    optical_connect(large_turn, tight_corner)

    fake_target = hairpin << pg.compass()
    fake_target.ymax = tight_corner.ymin
    fake_target.x = tight_corner.ports['wg_in_1'].x - inout_spacing
    # loose_corner = hairpin << route_ports(large_turn.ports['wg_out_1'], fake_target.ports['N'], wg_name='Strip')
    loose_corner = hairpin << wg_xsection.cell_euler_bend(90)
    optical_connect(loose_corner, large_turn)
    loose_corner.movex(fake_target.ports['N'].x - loose_corner.ports['wg_out_1'].x)
    little_bridge1 = hairpin << wg_xsection.cell_straight(large_turn.ports['wg_out_1'].x - loose_corner.ports['wg_in_1'].x)
    optical_connect(little_bridge1, large_turn)
    little_bridge2 = hairpin << wg_xsection.cell_straight(loose_corner.ports['wg_out_1'].y - fake_target.ports['N'].y)
    optical_connect(little_bridge2, loose_corner)
    hairpin.remove(fake_target)

    hairpin.add_port('wg_in_1', port=tight_corner.ports['wg_in_1'])
    hairpin.add_port('wg_out_1', port=little_bridge2.ports['wg_out_1'])
    hairpin.flatten()

    return hairpin


def waveguide_meander(pitch=2.5, n_turns=4, implant_radius=None, implant_width=None, wg_xsection=None):
    ''' waveguide_meander stands for switchback
        If pitch=None, it will make it as tight as possible, given the current algorithm.
        The pitch is approximate within 10 nanometers
        implant_width: if None, puts a 1um overhang
    '''
    if wg_xsection is None:
        wg_xsection = tech.waveguides('Strip')
    hairpin_xsection = tech.waveguides('Strip')
    if wg_xsection is hairpin_xsection:
        def central_wgstraight_generator(extent, halfie=False):
            return wg_xsection.cell_straight(extent)
    else:
        def central_wgstraight_generator(extent, halfie=False):
            magic_taper_len = 10
            try:
                return tapered_wg_straight(extent, hairpin_xsection, wg_xsection, taper_len=magic_taper_len, width_type='sine', cap_both=not halfie)
            except ValueError as err:
                if 'not long enough' in err.args[0]:
                    return hairpin_xsection.cell_straight(extent)



    D = Device('waveguide_meander')

    typical_hairpin = waveguide_hairpin_tweaked(inout_spacing=pitch, wg_xsection=hairpin_xsection)
    hairpin_extent = typical_hairpin.ysize

    Half_of_it = Device()
    enterer = Half_of_it << hairpin_xsection.cell_euler_bend(theta=90)
    enterer.reflect()
    center_extent = 0
    prev_leftturn = enterer
    for iTurn in range(n_turns):
        hairpin_right = Half_of_it << typical_hairpin
        hairpin_left = Half_of_it << typical_hairpin
        hairpin_left.reflect()
        tapered_wg_straight
        central_straight1 = Half_of_it << central_wgstraight_generator(center_extent)
        optical_connect(central_straight1, prev_leftturn)
        optical_connect(hairpin_right, central_straight1)
        central_straight2 = Half_of_it << central_wgstraight_generator(center_extent + hairpin_extent)
        center_extent += 2*hairpin_extent
        optical_connect(central_straight2, hairpin_right)
        optical_connect(hairpin_left, central_straight2)
        prev_leftturn = hairpin_left
        if iTurn == n_turns - 1:
            exiter = Half_of_it << central_wgstraight_generator(center_extent / 2 - hairpin_extent/2, halfie=True)
            optical_connect(exiter, prev_leftturn)
    Half_of_it.add_port('wg_out_1', port=exiter.ports['wg_out_1'])

    if False:  # too risky just to look pretty
        enterer_Sbend = Half_of_it << hairpin_xsection.cell_s_bend_by_lateral_offset((enterer.ports['wg_in_1'].y - Half_of_it.y),
                                                                            radius=hairpin_xsection.radius*2, tol=1e-1)
    else:
        enterer_Sbend = Half_of_it << hairpin_xsection.cell_straight(1)
    enterer_Sbend.connect('wg_out_1', enterer.ports['wg_in_1'])
    enterer_stick = Half_of_it << hairpin_xsection.cell_straight(abs(enterer_Sbend.ports['wg_in_1'].x - Half_of_it.xmax))
    enterer_stick.connect('wg_out_1', enterer_Sbend.ports['wg_in_1'])
    Half_of_it.add_port('wg_in_1', port=enterer_stick.ports['wg_in_1'])

    half1 = D << Half_of_it
    half2 = D << Half_of_it
    half2.connect('wg_out_1', half1.ports['wg_out_1'])
    D.add_port('wg_in_1', port=half2.ports['wg_in_1'])
    D.add_port('wg_out_1', port=half1.ports['wg_in_1'])

    if implant_radius is not None:
        center_orig = D.center
        Implant_shape = pg.circle(radius=implant_radius, layer=lys['dp_e'])
        Implant_shape.center = center_orig
        And_implant = pg.boolean(D, Implant_shape, 'AND')
        D.info['implanted_length'] = And_implant.area() / wg_xsection.rib_width
        if implant_width is None:
            target_offset = 1.0
        else:
            target_offset = (implant_width - wg_xsection.rib_width) / 2
        D << pg.offset(And_implant, target_offset, layer=lys['dp_e'])

    D.flatten()
    return D


def waveguide_spiral(pitch=3, full_turns=10, ports_same_side=True, implant_radius=None, implant_width=None, wg_xsection=None):
    ''' Implant parameters similar to waveguide_meander
    '''
    if wg_xsection is None:
        wg_xsection = tech.waveguides('Strip')

    D = Device('Spiral')
    core1 = D << wg_xsection.cell_bend(theta=180, radius=wg_xsection.radius)
    core2 = D << wg_xsection.cell_bend(theta=180, radius=wg_xsection.radius)
    core1.connect('wg_in_1', core2.ports['wg_in_1'])
    center_orig = core1.ports['wg_in_1'].midpoint
    start_radius = abs(core1.ports['wg_out_1'].y - core2.ports['wg_out_1'].y) / 2

    # keep track of total length, starting with center twist
    wg_length = 2 * np.pi * wg_xsection.radius

    this_radius = start_radius - pitch/2
    for iTurn in range(full_turns):
        this_radius += pitch
        this_piece1 = D << wg_xsection.cell_bend(theta=180, radius=this_radius)
        optical_connect(this_piece1, core1)
        core1 = this_piece1
        this_piece2 = D << wg_xsection.cell_bend(theta=180, radius=this_radius)
        optical_connect(this_piece2, core2)
        core2 = this_piece2
        wg_length += 2 * np.pi * this_radius

    if ports_same_side:
        this_radius += pitch
        this_piece1 = D << wg_xsection.cell_bend(theta=180, radius=this_radius)
        optical_connect(this_piece1, core1)
        core1 = this_piece1
        wg_length += np.pi * this_radius
    D.add_port('wg_in_1', port=core1.ports['wg_out_1'])
    D.add_port('wg_out_1', port=core2.ports['wg_out_1'])
    D.info['wg_length'] = wg_length

    if implant_radius is not None:
        Implant_shape = pg.circle(radius=implant_radius, layer=lys['dp_e'])
        Implant_shape.center = center_orig
        And_implant = pg.boolean(D, Implant_shape, 'AND')
        D.info['implanted_length'] = And_implant.area() / wg_xsection.rib_width
        if implant_width is None:
            target_offset = 1.0
        else:
            target_offset = (implant_width - wg_xsection.rib_width) / 2
        D << pg.offset(And_implant, target_offset, layer=lys['dp_e'])

    D.rotate(-D.ports['wg_out_1'].orientation)

    D.flatten()

    return D



#####
# Microrings and MZIs and other 4-ports
#
# De facto four-port optical device protocol
# - IN input: 'wg_in_1', THRU/BAR output: 'wg_out_1'
# - ADD input: 'wg_in_2', DROP/CROSS output: 'wg_out_2'
#####

def directional_coupler(access_spacing=[4, 4], gap=0.15, coupling_length=10, wg_xsection=None):
    '''
        The straight part of the coupler is slightly less than coupling length, which counts half of the exit arcs.
        Retains a consistent 45 degree exit from coupling region, also allows any spacing for the access waveguides.
        Strip waveguide is only option.

        Args:
            access_spacing (float, list): spacing of in/out waveguides that will have ports on them
            gap (float): coupler gap
            coupling_length (float): length of coupling section
    '''
    if np.isscalar(access_spacing):
        access_spacing = [access_spacing] * 2
    if wg_xsection is None:
        wg_xsection = tech.waveguides('Strip')
    # wg_width = wg_xsection.get_by_layer('wg_deep').width
    wg_radius = wg_xsection.radius

    D=Device()

    ptlist = [np.array([0, 0])]
    def go_by(trans, y=None):
        if y is not None:
            trans = np.array([trans, y])
        ptlist.append(ptlist[-1] + np.array(trans))

    len45 = wg_radius
    go_by(wg_radius, 0)
    go_by(2 * wg_radius, -access_spacing[0] / 2 + len45)
    go_by(len45, -len45)
    go_by(coupling_length, 0)
    go_by(len45, len45)
    go_by(2 * wg_radius, access_spacing[1] / 2 - len45)
    go_by(wg_radius, 0)

    ptlist = np.array(ptlist)
    top_wg = D << route_wg_points(ptlist, wg_name=wg_xsection)
    ptlist[:, 1] = - ptlist[:, 1]
    bottom_wg = D << route_wg_points(ptlist, wg_name=wg_xsection)
    bottom_wg.ymax = top_wg.ymin - gap

    D.add_port(port = top_wg.ports['wg_in_1'])
    D.add_port(port = top_wg.ports['wg_out_1'])
    D.add_port(port = bottom_wg.ports['wg_in_1'], name = 'wg_in_2')
    D.add_port(port = bottom_wg.ports['wg_out_1'], name = 'wg_out_2')
    D.flatten()

    return D


def taperize(port, length=10, out_wg_width=1, offset=0, layer=lys['wg_deep']):
    ''' Small taper that senses the width of the port, makes an S-bend, and connects to it.
        This modifies the device so that the old port is replaced by the one at the end of the taper.
        Positive offset is to the right.
    '''
    # First make it inwards pointing so we can use route_basic
    offset_vec = np.array([length, -offset])
    theta = np.radians(port.orientation)
    c, s = np.cos(theta), np.sin(theta)
    rotated_offset_vec = np.array([[c,-s], [s, c]]).dot(offset_vec)
    new_midpoint = np.array(port.midpoint) + rotated_offset_vec
    port_temp = Port(name=port.name, midpoint=new_midpoint,
                     width=out_wg_width, orientation=port.orientation + 180)

    true_parent = port.parent
    if isinstance(true_parent, phidl.device_layout.DeviceReference):
        print('Be careful. taperize is not fully tested with DeviceReferences')
        parent_device = true_parent.parent
    else:
        parent_device = true_parent
    parent_device << pr.route_basic(port_temp, port, layer=layer)
    # Then flip it so it can override the previous non-released ports
    port_temp.orientation = port.orientation
    true_parent.ports[port.name] = port_temp
    port_temp.parent = true_parent


def adiabatic5050_coupler(is1550=False, out_wg_width=0.35, do_release=True, **hard_overrides):
    ''' This supposedly guarantees a 50:50 split ratio as long as it is longer than a certain amount
        The do_release option puts in tapered S-bends so all the outputs end up at out_wg_width.
    '''
    if is1550:
        wg_widths = np.array([[0.4, 0.5], [0.4, 0.4]])  # port 2 (in), 1 (in), 3 (out), 4 (out)
        gap1 = 1.0
        gap2 = 0.1
        taper_len = 200
    else:
        wg_widths = np.array([[0.3, 0.4], [0.3, 0.3]])  # port 2 (in), 1 (in), 3 (out), 4 (out)
        gap1 = 0.6
        gap2 = 0.1
        taper_len = 144
    wg_widths = hard_overrides.get('wg_widths', wg_widths)
    gap1 = hard_overrides.get('gap1', gap1)
    gap2 = hard_overrides.get('gap2', gap2)
    taper_len = hard_overrides.get('taper_len', taper_len)

    # These are cross sections at the beginning and end, done this way so we can use taper_xsections for the hard part
    xsection_start = tech.WGXSection(
        components=[tech.WGXSectionComponent(width=wg_widths[0][0], offset=0,
                                             layer='wg_deep'),
                    tech.WGXSectionComponent(width=wg_widths[0][1],
                                             offset=np.mean(wg_widths[0]) + gap1,
                                             layer='wg_deep'),
                    ],
                                     )
    xsection_end = tech.WGXSection(
        components=[tech.WGXSectionComponent(width=wg_widths[1][0], offset=0,
                                             layer='wg_deep'),
                    tech.WGXSectionComponent(width=wg_widths[1][1],
                                             offset=np.mean(wg_widths[1]) + gap2,
                                             layer='wg_deep'),
                    ],
                                   )

    D = Device()
    coupling_region = D << taper_xsections(xsection_start, xsection_end, taper_len=taper_len,
                                           route_basic_options=dict(path_type='straight', width_type='straight'))

    midpoints = np.array([[(0, 0), (0, np.mean(wg_widths[0]) + gap1)],
                         [(taper_len, 0), (taper_len, np.mean(wg_widths[1]) + gap2)]])

    # These are temporary if you are going to release it,
    # but they help layout the release sections in that case
    D.add_port(name='wg_in_1', midpoint=midpoints[0][1],
               orientation=180, width=wg_widths[0][1])
    D.add_port(name='wg_in_2', midpoint=midpoints[0][0],
               orientation=180, width=wg_widths[0][0])
    D.add_port(name='wg_out_1', midpoint=midpoints[1][1],
               orientation=0, width=wg_widths[1][1])
    D.add_port(name='wg_out_2', midpoint=midpoints[1][0],
               orientation=0, width=wg_widths[1][0])

    if do_release:
        for port_name, port in D.ports.items():
            offset = -15
            if port.name.endswith('_2') ^ ('in' in port.name):
                offset *= -1
            taperize(port, length=20, out_wg_width=out_wg_width, offset=offset)

    # Duplicate naming to fit other conventions. Ordering doesnt really make sense, but too bad
    D.add_port(1, port=D.ports['wg_in_1'])
    D.add_port(2, port=D.ports['wg_out_1'])
    D.add_port(3, port=D.ports['wg_out_2'])
    D.add_port(4, port=D.ports['wg_in_2'])

    return D


def clip_to(device, coord, cut_angle=0, ignore_layers=[]):
    ''' it modifies the input device and does not return '''
    if np.isscalar(coord) or len(coord) == 1:
        if cut_angle % 360 not in [0, 90, 180, 270]:
            raise ValueError('Scalar cut coordinates only work with axis-aligned angles')
        else:
            if cut_angle % 360 in [0, 180]:
                coord = [device.x, coord]
            else:
                coord = [coord, device.y]
    cut_maxlen = np.sqrt(np.sum(np.asarray(device.size) ** 2)) * 1.1
    cut_box = pg.rectangle(size=(cut_maxlen, cut_maxlen))
    cut_box.x = coord[0]
    cut_box.ymax = coord[1]
    cut_box.rotate(cut_angle, center=coord)

    device.flatten()
    for layer in device.layers:
        if layer in map(lambda lay: getattr(lay, 'gds_layer'), ignore_layers):
            continue
        new_dev = Device()
        for p in filter(lambda p: p.layers[0] == layer, device.polygons):
            device.remove(p)
            new_dev.add(p)
        clipped_dev = pg.boolean(new_dev, cut_box, operation='A-B', layer=layer)
        device.add_ref(clipped_dev)
    device.flatten()


def microracetrack(radius, gap=.15, parallel_len=0, perpendicular_len=0, bus_type='add_drop',
              ring_xsect=None, bus_xsect=None, clip_dopants=False, clip_offset=1):
    '''
        Strip waveguide is only option. However, you can add stuff down on top later (like dopants)

        Args:
            access_spacing (float, list): spacing of in/out waveguides that will have ports on them
            gap (float): coupler gap
            coupling_length (float): length of coupling section
            bus_type (str): can be 'add_drop', 'all_pass', or 'add_drop_terminated'
            clip_dopants (bool): If the ring has a doped profile, these dopants will be taken away near the coupling regions
            clip_offset (float): distance away from bus to clip the dopants. Relative to the center of the MRR waveguide
    '''
    if ring_xsect is None:
        ring_xsect = tech.waveguides('Strip')
    if bus_xsect is None:
        bus_xsect = tech.waveguides('Strip')

    ring_itself = Device()
    last_element = None
    def ring_append(dev):
        nonlocal last_element
        dev_ref = ring_itself << dev
        if last_element is not None:
            optical_connect(dev_ref, last_element)
        last_element = dev_ref

    # Turn up the resolution
    fracture_arc = 1
    arc_partial = ring_xsect.cell_bend(90 / fracture_arc, radius=radius, angle_resolution=0.05)
    def ring_append_arc():
        for i in range(fracture_arc):
            ring_append(arc_partial)
    ring_append_arc()
    ring_append(ring_xsect.cell_straight(perpendicular_len))
    ring_append_arc()
    ring_append(ring_xsect.cell_straight(parallel_len))
    ring_append_arc()
    ring_append(ring_xsect.cell_straight(perpendicular_len))
    ring_append_arc()
    ring_append(ring_xsect.cell_straight(parallel_len))
    ring_itself.flatten()
    del ring_append

    # Put in bus waveguides. One of them might be removed later down
    D = Device('microring')
    bottom_bus = D.add_ref(bus_xsect.cell_straight(ring_itself.xsize), alias='bottom_bus')
    ring_ref = D.add_ref(ring_itself, alias='ring')
    top_bus = D.add_ref(bus_xsect.cell_straight(ring_itself.xsize), alias='top_bus')

    dopant_adjusted_gap = gap - (ring_xsect.full_width - ring_xsect.rib_width) / 2
    ring_ref.ymin = bottom_bus.ymax + dopant_adjusted_gap
    ring_ref.xmin = bottom_bus.xmin
    top_bus.ymin = ring_ref.ymax + dopant_adjusted_gap

    # Handle the ports and busses
    D.add_port(port = top_bus.ports['wg_in_1'])
    D.add_port(port = top_bus.ports['wg_out_1'])
    if bus_type == 'add_drop':
        D.add_port(port = bottom_bus.ports['wg_in_1'], name = 'wg_in_2')
        D.add_port(port = bottom_bus.ports['wg_out_1'], name = 'wg_out_2')
    elif bus_type == 'add_drop_terminated':
        for turn_angle, bus_port in ((-90, 'wg_out_1'), (90, 'wg_in_1')):
            to_terminus = D << bus_xsect.cell_bend(turn_angle)
            to_terminus.connect('wg_in_1', bottom_bus.ports[bus_port])
            terminus = D << wg_terminator(wg_to_terminate=bus_xsect)
            optical_connect(terminus, to_terminus)
    elif bus_type == 'all_pass':
        D.references.remove(bottom_bus)
    else:
        raise ValueError('bus_type not accepted: {}'.format(bus_type))

    # Clip dopants in the ring away from the bus(es)
    if clip_dopants:
        clip_to(ring_itself, ring_itself.ymax - ring_xsect.full_width/2 - clip_offset, cut_angle=180, ignore_layers=[lys[l] for l in ['wg_deep', 'wg_shallow']])
        if bus_type != 'all_pass':
            clip_to(ring_itself, ring_itself.ymin + ring_xsect.full_width/2 + clip_offset, cut_angle=0, ignore_layers=[lys[l] for l in ['wg_deep', 'wg_shallow']])

        Cutter = Device()
        cut = Cutter << pg.rectangle((radius / 2, 10*radius))  # primarily for DRC purposes
        cut.x = ring_itself.x
        if bus_type != 'all_pass':
            cut.y = ring_itself.y
        else:
            cut.ymin = ring_itself.y

        new_dopants = Device()
        for comp in ring_xsect.components:
            if comp.layer not in ['wg_deep', 'wg_shallow']:
                spec = (lys[comp.layer].gds_layer, lys[comp.layer].gds_datatype)
                try:
                    implant_polygons = ring_itself.get_polygons(by_spec = True)[spec]
                except KeyError:
                    continue
                implant_ring = Device()
                implant_ring.add_polygon(implant_polygons)
                # implant_ring.center = ring_itself.center
                # implant_ring = pg.offset(D['ring'], distance=implant_overhang)
                new_dopants << pg.boolean(implant_ring, Cutter, 'A-B', layer=lys[comp.layer])
        ring_itself.remove_layers([lys['wg_deep'], lys['wg_shallow']], invert_selection=True)
        ring_itself << new_dopants

        ring_itself.flatten()
    return D


def trimmed_microracetrack(*mrr_args, implant_proportion=0.5, implant_layer=None, **mrr_kwargs):
    # note it's a proportion of the curved sections. the proportion will be wrong if you have parallel segments.
    base_mrr = microracetrack(*mrr_args, **mrr_kwargs)
    if implant_layer is None:
        implant_layer = lys['dp_e']
    ring_xsect = mrr_kwargs.get('ring_xsect', tech.waveguides('Strip'))
    implant_annular_width = ring_xsect.rib_width * 2
    implant_angle = implant_proportion * 180
    if len(mrr_args) >= 1:
        radius = mrr_args[0]
    else:
        radius = mrr_kwargs['radius']
    og_maxx = base_mrr.xmax
    og_minx = base_mrr.xmin

    if implant_proportion != 0:
        og_y = base_mrr.y
        implant_shape = pg.arc(radius=radius, width=implant_annular_width,
                               theta=implant_angle, start_angle=-implant_angle/2,
                               layer=implant_layer)
        implant1 = base_mrr << implant_shape
        implant1.y = og_y
        implant1.xmax = og_maxx - ring_xsect.rib_width / 2 + implant_annular_width / 2
        implant2 = base_mrr << implant_shape
        implant2.rotate(180)
        implant2.y = og_y
        implant2.xmin = og_minx + ring_xsect.rib_width / 2 - implant_annular_width / 2
    return base_mrr


def long_microresonator(excursion=100, coupling_length=5, gap=0.15, wg_xsection=None):
    '''
        A long microring meant to measure index in different waveguides

        wg_name is what will be passed to route_wg_points,
        so if None, it will pick the long haul

        This is not ideal behavior. You should be able to explicitly tell a waveguide and either way have the bend be in a strip
    '''
    bend_wg = tech.waveguides('Strip')
    if wg_xsection is None:
        wg_xsection = tech.waveguides('Strip')

    # Start off with the couplers and short connection
    D = Device()
    dc_shape = directional_coupler(access_spacing=10, coupling_length=coupling_length, gap=gap)
    dc1 = D << dc_shape
    dc2 = D << dc_shape
    one_eighty = D << bend_wg.cell_bend(180)
    one_eighty.connect('wg_in_1', dc1.ports['wg_in_2'])
    dc2.connect('wg_in_1', one_eighty.ports['wg_out_1'])

    # Then the excursion loop
    magic_taper_len = 7
    long_arm_shape = tapered_wg_straight(excursion, tech.waveguides('Strip'), wg_xsection, taper_len=magic_taper_len)
    long_arm_top = D << long_arm_shape
    long_arm_bottom = D << long_arm_shape
    one_eighty_far = D << bend_wg.cell_bend(-180)
    long_arm_top.connect('wg_in_1', dc1.ports['wg_out_2'])
    optical_connect(one_eighty_far, long_arm_top)
    optical_connect(long_arm_bottom, one_eighty_far)

    for turn_angle, bus_port in ((-90, 'wg_out_2'), (90, 'wg_in_2')):
        to_terminus = D << bend_wg.cell_bend(turn_angle)
        to_terminus.connect('wg_in_1', dc2.ports[bus_port])
        terminus = D << wg_terminator(wg_to_terminate=bend_wg)
        optical_connect(terminus, to_terminus)

    D.add_port(port=dc1.ports['wg_in_1'])
    D.add_port(port=dc1.ports['wg_out_1'])
    return D


def demux(n_chan=3, dL=.1, spacing=100, doped_tips=True, **microring_kwargs):
    ''' Microring demultiplexer.
        Neighboring rings have slightly different perimeters due to perpendicular length differences, but same radius

        Args:
            radius (float): of arcs used in the ring
            n_chan (int): number of MRRs and channels
            dL (float): difference in perimeter between neighbors (in microns)
            spacing (float): physical distance between placement of neighbors
            microring_kwargs (/*/*kwargs): fed through to microring. perpendicular_length here sets that of the first MRR only
                - anything that is a list will be distributed, since there are no list arguments to microring
    '''
    D = Device()
    mrr_list = []
    initial_perpendicular_length = microring_kwargs.pop('perpendicular_len', 0)
    microring_kwargs['radius'] = microring_kwargs.pop('radius', 8)
    distributed_kwargs = dict()
    for k, v in microring_kwargs.items():
        if isinstance(v, (list, tuple, np.ndarray)):
            distributed_kwargs[k] = v
    # Instantiate the references
    for iRing in range(n_chan):
        this_perp_len = iRing * dL / 2 + initial_perpendicular_length
        full_kwargs = microring_kwargs.copy()
        for k, vlist in distributed_kwargs.items():
            full_kwargs[k] = vlist[iRing]
        mrr_list.append(D << microracetrack(perpendicular_len=this_perp_len, bus_type='add_drop', **full_kwargs))
    # Arrange and connect them
    bus_xsect = microring_kwargs.get('bus_xsect', tech.waveguides('Strip'))
    Interbus = bus_xsect.cell_straight(spacing - 2 * microring_kwargs['radius'])
    for mrr_from, mrr_to in zip(mrr_list[:-1], mrr_list[1:]):
        interbus = D << Interbus
        interbus.connect('wg_in_1', mrr_from.ports['wg_out_1'])
        mrr_to.connect('wg_in_1', interbus.ports['wg_out_1'])
        # This would be for weight bank
        # D << pr.route_basic(mrr_from.ports['wg_out_2'], mrr_to.ports['wg_in_2'], layer=lys['wg_deep'])
    # Promote ports
    D.add_port(port=mrr_list[0].ports['wg_in_1'])
    D.add_port(port=mrr_list[-1].ports['wg_out_1'])
    for iRing, mrr in enumerate(mrr_list):
        D.add_port(port=mrr.ports['wg_in_2'], name='drop' + str(iRing + 1))
        to_terminus = D << bus_xsect.cell_bend(-90)
        to_terminus.connect('wg_in_1', mrr.ports['wg_out_2'])
        terminus = D << wg_terminator(wg_to_terminate=bus_xsect, doped=doped_tips)
        optical_connect(terminus, to_terminus)
    return D


def loss_rings(bus_width=0.4, ring_wg_width=SM_WG_WIDTH, ext_wg_xsect=None, gap_vector=[.2] * 4, R0=18, dR=0.008, spacing=100):
    ''' Like a demux but you can specify different kinds of parameters '''
    D = Device('Loss_ring')
    n_rings = len(gap_vector)
    wg_width_vector = n_rings * [ring_wg_width]
    # if np.isscalar(wg_width_vector):
    #     wg_width_vector = [wg_width_vector] * n_rings
    # if np.isscalar(gap_vector):
    #     gap_vector = [gap_vector] * n_rings
    bus_xsect = tech.WGXSection(components=[tech.WGXSectionComponent(width=bus_width,
                                                                     offset=0, layer='wg_deep')])

    mrr_list = []
    # Instantiate the references
    for iRing in range(n_rings):
        ring_xsect = tech.WGXSection(components=[tech.WGXSectionComponent(width=wg_width_vector[iRing],
                                                                          offset=0, layer='wg_deep')])
        mrr_list.append(D << microracetrack(radius=(R0 + dR * iRing), gap=gap_vector[iRing],
                                       ring_xsect=ring_xsect, bus_xsect=bus_xsect,
                                       bus_type='add_drop'))
    # Arrange and connect them
    for mrr_from, mrr_to in zip(mrr_list[:-1], mrr_list[1:]):
        mrr_to.xmin = mrr_from.xmin + spacing
        mrr_to.y = mrr_from.y
        D << pr.route_basic(mrr_from.ports['wg_out_1'], mrr_to.ports['wg_in_1'], layer=lys['wg_deep'])
    # Termination and promote overall ADD/THRU
    in_port = mrr_list[0].ports['wg_in_1']
    out_port = mrr_list[-1].ports['wg_out_1']
    if ext_wg_xsect is not None:
        out_taper = D << taper_xsections(bus_xsect, ext_wg_xsect, taper_len=10)
        out_taper.connect('wg_in_1', out_port)
        out_port = out_taper.ports['wg_out_1']
        in_taper = D << taper_xsections(ext_wg_xsect, bus_xsect, taper_len=10)
        in_taper.connect('wg_out_1', in_port)
        in_port = in_taper.ports['wg_in_1']
    D.add_port(port=in_port)
    D.add_port(port=out_port)

    # Promote DROPs
    drop_interface = taper_xsections(bus_xsect, ext_wg_xsect, taper_len=10)
    for iRing, mrr in enumerate(mrr_list):
        this_drop_port = mrr.ports['wg_in_2']
        # this_drop_port.orientation += 180
        if ext_wg_xsect is not None:
            interfacer = D << drop_interface
            interfacer.connect('wg_in_1', this_drop_port)
            this_drop_port = interfacer.ports['wg_out_1']
        D.add_port('drop' + str(iRing + 1), port=this_drop_port)
        to_terminus = D << bus_xsect.cell_bend(-90, radius=10)
        to_terminus.connect('wg_in_1', mrr.ports['wg_out_2'])
        terminus = D << wg_terminator(wg_to_terminate=bus_xsect)
        optical_connect(terminus, to_terminus)
    return D


def TOAD():
    ''' SOAs in balanced arms of a Mach-Zehnder
    '''
    D = Device()
    straight_sec = tech.waveguides('Strip').cell_straight(50)
    amp_sec = SOA(30, mid_wid=3)
    arm_top = D << device_from_sections([straight_sec, amp_sec])
    arm_bottom = D << device_from_sections([amp_sec, straight_sec])
    d1 = D << directional_coupler(access_spacing=[4, 10])
    d2 = D << directional_coupler(access_spacing=[10, 4])

    arm_top.move(arm_top.ports['wg_in_1'].midpoint,
                 d1.ports['wg_out_1'].midpoint)
    arm_bottom.move(arm_bottom.ports['wg_in_1'].midpoint,
                    d1.ports['wg_out_2'].midpoint)
    d2.move(d2.ports['wg_in_1'].midpoint,
            arm_top.ports['wg_out_1'].midpoint)
    return D


def focusing_gc(pitch=None, duty=None, profile=None, num_periods=50, offset=10, arc_angle=40, wg_width=None):
    ''' The profile starts with the presence of oxide
    '''
    if pitch is not None or duty is not None:
        if profile is not None:
            raise ValueError('You can only specify either pitch and duty or the profile')
    elif profile is not None:
        if pitch is not None or duty is not None:
            raise ValueError('You can only specify either pitch and duty or the profile')
    # defaults
    if profile is None:
        if pitch is None:
            pitch = .2
        if duty is None:
            duty = .5
        profile = [pitch * duty, pitch * (1 - duty)]
    pitch = sum(profile)

    if wg_width is None:
        wg_xsection = tech.waveguides('Strip')
        wg_width = wg_xsection.get_by_layer('wg_deep')[0].width

    D = Device()

    radius = offset
    for _ in range(num_periods):
        for iTooth in range(len(profile)):
            if iTooth % 2 == 1:
                tooth = pg.arc(radius=radius+profile[iTooth]/2, width=profile[iTooth], theta=arc_angle, start_angle=-arc_angle/2,
                               angle_resolution=0.5, layer=lys['wg_deep'])
                D << tooth
            radius += profile[iTooth]
    angle_rad = arc_angle/2 * np.pi / 180
    port_shift = wg_width / 2 / np.tan(angle_rad)
    central_triangle = Device()
    central_triangle.add_polygon([[port_shift, -wg_width/2],
                                  [port_shift, wg_width/2],
                                  [radius * np.cos(angle_rad), radius * np.sin(angle_rad)],
                                  [radius * np.cos(angle_rad), -radius * np.sin(angle_rad)]],
                                 layer=255)
    central_disk = pg.circle(offset, angle_resolution=0.5, layer=255)
    bulkish = pg.boolean(central_disk, central_triangle, 'and', layer=lys['wg_deep'])
    D << bulkish

    # ports
    stick_length = 0.555
    exit_stick = D << pg.compass((stick_length, wg_width), layer=lys['wg_deep'])
    exit_stick.move(origin=exit_stick.ports['E'].midpoint, destination=(port_shift, 0))
    D.add_port(name='wg_in_1', port=exit_stick.ports['W'])
    D = D.rotate(180)
    D.flatten()

    return D


connor_nominals_0 = dict(pitch = .815,
                         duty = 0.526,
                         subwavelength_tooth_width = 0.2,
                         arc_angle=35, num_periods=18, offset=19.033)

connor_nominals_23 = dict(pitch = 0.511,
                          duty = 0.553,
                          subwavelength_tooth_width = 0.1,
                          arc_angle=35, num_periods=29, offset=19.033)


def subwavelength_gc(pitch=.5, duty=.55, subwavelength_tooth_width=.1, **kwargs):
    ''' arguments are passed to focusing_gc '''
    gap_width = (pitch * (1 - duty) - subwavelength_tooth_width) / 2
    profile = [gap_width, subwavelength_tooth_width, gap_width, pitch * duty]
    GC = focusing_gc(profile=profile, **kwargs)
    left_side = (GC.xmin, GC.y)
    XCLUDE = GC << pg.rectangle((GC.xsize, GC.ysize), layer=lys['DRC_exclude'])
    XCLUDE.xmin, XCLUDE.y = left_side
    return GC



def bullseye_coupler(pitch=None, duty=None, profile=None, num_periods=20, offset=1):
    ''' The profile starts with the presence of oxide
    '''
    if pitch is not None or duty is not None:
        if profile is not None:
            raise ValueError('You can only specify either pitch and duty or the profile')
    elif profile is not None:
        if pitch is not None or duty is not None:
            raise ValueError('You can only specify either pitch and duto or the profile')
    # defaults
    if profile is None:
        if pitch is None:
            pitch = .544
        if duty is None:
            duty = .45
        profile = [pitch * duty, pitch * (1 - duty)]
    pitch = sum(profile)

    D = Device()
    radius = offset
    for _ in range(num_periods):
        for iTooth in range(len(profile)):
            if iTooth % 2 == 1:
                tooth = pg.arc(radius=radius+profile[iTooth]/2, width=profile[iTooth], theta=360,
                               angle_resolution=0.5, layer=lys['wg_deep'])
                D << tooth
            radius += profile[iTooth]
    central_disk = pg.circle(offset, angle_resolution=0.5, layer=lys['wg_deep'])
    D << central_disk
    return D


def route_manhattan0(
        port1,
        port2,
        awayness = 50,
        layer=0):
    assert port1.orientation == port2.orientation
    D = Device()
    bend = pr._arc(theta=90, width=port1.width, radius=awayness, layer=layer)
    top_downwards = D << bend
    top_downwards.connect(1, port1)
    bottom_upwards = D << bend
    bottom_upwards.connect(2, port2)
    D << pr.route_basic(top_downwards.ports[2], bottom_upwards.ports[1], layer=layer)
    return D


#####
# Electrical
#####

def SGS_pads(pad_width=120, pitch=150):
    D = Device()
    pad = pg.compass((pad_width, pad_width), layer=lys['m5_wiring'])
    pad_refs = []
    for iPad in range(5):
        pad_refs.append(D << pad)
        if iPad > 0:
            pad_refs[iPad].x = pad_refs[iPad - 1].x
            pad_refs[iPad].y = pad_refs[iPad - 1].y - pitch
    D << route_manhattan0(pad_refs[0].ports['W'], pad_refs[2].ports['W'], awayness=80, layer=lys['m5_wiring'])
    D << route_manhattan0(pad_refs[2].ports['W'], pad_refs[4].ports['W'], awayness=80, layer=lys['m5_wiring'])
    D.add_port(port=pad_refs[0].ports['E'], name='top_G')
    D.add_port(port=pad_refs[1].ports['E'], name='top_S')
    D.add_port(port=pad_refs[2].ports['E'], name='middle_G')
    D.add_port(port=pad_refs[3].ports['E'], name='bottom_S')
    D.add_port(port=pad_refs[4].ports['E'], name='bottom_G')
    return D


def fourPoint_parallels(length, width, gap, layer=lys['m5_wiring']):
    ''' Two wires that run parallel
    '''
    D = Device()
    fil_top_dev, cap_NW, cap_NE = filament(length, width, layer=layer)
    fil_top = D << fil_top_dev
    fil_bot_dev, cap_SW, cap_SE = filament(length, width, layer=layer)
    fil_bot_dev.ymax = fil_top.ymin - gap
    fil_bot = D << fil_bot_dev

    # Placement of pads
    pad_backoff = 100
    pad = pg.compass(WIREBOND_SIZE, layer=lys['m5_wiring'])

    sPad_left = D << pad
    sPad_left.xmax = cap_SW.xmin - pad_backoff
    sPad_left.ymax = cap_SW.y - 2 * pad_backoff
    nPad_left = D << pad
    nPad_left.xmax = cap_NW.xmin - pad_backoff
    nPad_left.ymin = cap_NW.y + 2 * pad_backoff

    sPad_right = D << pad
    sPad_right.xmin = cap_SE.xmax + pad_backoff
    sPad_right.ymax = cap_SE.y - 2 * pad_backoff
    nPad_right = D << pad
    nPad_right.xmin = cap_NE.xmax + pad_backoff
    nPad_right.ymin = cap_NE.y + 2 * pad_backoff

    # Breakout to pads
    for a, b in [(cap_SW.ports['S'], sPad_left.ports['N']),
                 (cap_NW.ports['N'], nPad_left.ports['S']),
                 (cap_SE.ports['S'], sPad_right.ports['N']),
                 (cap_NE.ports['N'], nPad_right.ports['S'])]:
        D << pr.route_basic(a, b, path_type='straight', width_type='sine', layer=lys['m5_wiring'])

    # Label
    pl = D << nc_text(text='W = {}\nG = {}'.format(width, gap), size=80, layer=lys['m5_wiring'])
    pl.y = nPad_right.y
    pl.x = fil_top.x

    return D


def fourPoint_wire(length, width, layer=lys['m5_wiring']):
    ''' One wire with two connections on each end. For doing 4-point resistance measurements
    '''
    D = Device()
    fil_device, fil_left, fil_right = filament(length, width, layer=layer)
    fil_ref = D << fil_device

    # Placement of pads
    pad_backoff = 100
    pad = pg.compass(WIREBOND_SIZE, layer=lys['m5_wiring'])

    iPad_left = D << pad
    iPad_left.xmax = fil_left.xmin - pad_backoff
    iPad_left.ymax = fil_left.y - 2 * pad_backoff
    vPad_left = D << pad
    vPad_left.xmax = fil_left.xmin - pad_backoff
    vPad_left.ymin = fil_left.y + 2 * pad_backoff

    iPad_right = D << pad
    iPad_right.xmin = fil_right.xmax + pad_backoff
    iPad_right.ymax = fil_right.y - 2 * pad_backoff
    vPad_right = D << pad
    vPad_right.xmin = fil_right.xmax + pad_backoff
    vPad_right.ymin = fil_right.y + 2 * pad_backoff

    # Breakout to pads
    for a, b in [(fil_left.ports['N'], vPad_left.ports['S']),
                 (fil_left.ports['S'], iPad_left.ports['N']),
                 (fil_right.ports['N'], vPad_right.ports['S']),
                 (fil_right.ports['S'], iPad_right.ports['N'])]:
        D << pr.route_basic(a, b, path_type='straight', width_type='sine', layer=lys['m5_wiring'])

    # Label
    pl = D << nc_text(text='W = {}'.format(width), size=80, layer=lys['m5_wiring'])
    pl.y = vPad_right.y
    pl.x = fil_ref.x


def vanDerPauw(radius=250, layer=lys['m5_wiring']):
    D = Device()
    circ = pg.circle(radius=radius, layer=0)
    snip = pg.rectangle(size=(radius * 4/5, radius / 5))
    snip.y = circ.y
    snip.xmax = circ.xmax
    circ = pg.boolean(circ, snip, 'A-B', layer=layer)
    snip.xmin = circ.xmin
    circ = pg.boolean(circ, snip, 'A-B', layer=layer)

    snip.rotate(90)
    snip.x = circ.x
    snip.ymax = circ.ymax
    circ = pg.boolean(circ, snip, 'A-B', layer=layer)
    snip.ymin = circ.ymin
    circ = pg.boolean(circ, snip, 'A-B', layer=layer)
    D << circ
    return D


def filament(length, width, layer=lys['m5_wiring']):
    ''' A wire with caps so that connections can be made perpecdicular at the ends '''
    D = Device()
    filament = D << pg.rectangle((length, width), layer=layer)
    cap_left = D << pg.compass((width, width), layer=layer)
    cap_left.xmin = filament.xmin
    cap_left.y = filament.y
    cap_right = D << pg.compass((width, width), layer=layer)
    cap_right.xmax = filament.xmax
    cap_right.y = filament.y
    return D, cap_left, cap_right



#####
# Generic stuff
#####

def clip(D_in, box):
    D_out = phidl.Device()
    D_db = D_in.get_polygons(by_spec=True)
    box_db = box.get_polygons(by_spec=True)
    if not len(box_db) == 1:
        raise ValueError('Only one layer at a time')
    box_polys_by_layer = list(box_db.values())
    if not len(box_polys_by_layer) == 1:

        # to be continued


        raise ValueError('Only one box at a time for now')
    box_operand = box
    for layer in D_polys.keys():
        pass


def blinders(dev_to_blind, pitch=6):
    cutter = Device()
    center = np.array([dev_to_blind.x, dev_to_blind.ymin])
    for yval in np.arange(0, dev_to_blind.ysize, pitch):
        center[1] = yval
        cut = cutter << pg.rectangle((dev_to_blind.xsize, 3), layer=lys['dp_e'])
        cut.center = center
    blinded = pg.boolean(dev_to_blind, cutter, 'and', layer=dev_to_blind.layers.pop())
    return blinded


###############################################################################
#
# SONIA + CHILES DEVICES
# Integrated passive photonics and active photonics
#
###############################################################################

def bs_tree_add_devices(D_list, Input_grating = None, Output_grating = None,
                       taper_length = 30,
                       layer_separation = 100,
                       stagger_neighbors = 0,
                       width_wg = SM_WG_WIDTH,
                       splitter=None,
                       parallel_fibers=False, rotate90_fibers=False):

    # makes a beamsplitter tree and adds all of the devices in a list
    D = Device()

    # beam splitter tree

    tree_depth = int(np.ceil(np.sqrt(np.size(D_list)+1)))
    input_grating = D << Input_grating
    output_grating = D << Output_grating
    maxsize=np.max([D_list[i].size[1] for i in range(0,len(D_list))])
    bs_tree_ref = D << bs_tree(output_separation = maxsize*1.1, tree_depth = tree_depth, layer_separation = layer_separation, wg_width=width_wg, splitter=splitter)
    if not parallel_fibers:  # old style
        input_grating.rotate(90)
        input_grating.x = bs_tree_ref.ports[1].midpoint[0]
        input_grating.y= bs_tree_ref.ports[1].midpoint[1]
        input_grating.movey(-1000)
        input_grating.movex(-500)
        D << pr.route_manhattan(port1 = input_grating.ports[1], port2 = bs_tree_ref.ports[1], bendType = 'gradual',layer = lys['wg_deep'])
        output_grating.rotate(-90)
        output_grating.xmin = bs_tree_ref.ports[2].midpoint[0]
        output_grating.y = bs_tree_ref.ports[2].midpoint[1]
        output_grating.movey(Output_grating.size[1] + 2*SM_BEND_RADIUS*10)
        output_grating.x = input_grating.x + 50
        D << pr.route_manhattan(port1 = output_grating.ports[1], port2 = bs_tree_ref.ports[2], bendType='gradual', layer = lys['wg_deep'])
    else:  # sheep style for big glueboxes
        if not rotate90_fibers:  # aligned E-W, perpendicular to the line going through the row of output ports
            input_grating.connect(1, bs_tree_ref.ports[1])
            # output_grating.xmin = input_grating.xmin
            # output_grating.y = bs_tree_ref.ports[2].midpoint[1]
            # output_grating.y += SM_BEND_RADIUS * 10  # could be as low as *2 but this only works for 1220 waveguides, and we also want 1550
            output_grating.rotate(180)
            output_grating.xmin = bs_tree_ref.xmax + SM_BEND_RADIUS * 10
            # output_grating.ymin = bs_tree_ref.ymax
            # D << pr.route_manhattan(port1 = output_grating.ports[1], port2 = bs_tree_ref.ports[2], bendType='gradual', layer = lys['wg_deep'])
            output_grating.connect(1, bs_tree_ref.ports[2])
        else:  # parallel
            candy_cane = Device()
            bent_sec = candy_cane << pg.arc(theta=-90, radius=20, width=input_grating.ports[1].width, layer=lys['wg_deep'])
            straight_sec = candy_cane << pg.compass((input_grating.ysize / 2 + 20, input_grating.ports[1].width), layer=lys['wg_deep'])
            straight_sec.connect('W', bent_sec.ports[2])
            candy_cane.add_port('wg_in_1', port=bent_sec.ports[1])
            candy_cane.add_port('wg_out_1', port=straight_sec.ports['E'])

            candy_in = D << candy_cane
            candy_in.connect('wg_out_1', bs_tree_ref.ports[1])
            input_grating.connect(1, candy_in.ports['wg_in_1'])
            candy_out = D << candy_cane
            candy_out.connect('wg_out_1', bs_tree_ref.ports[2])
            output_grating.connect(1, candy_out.ports['wg_in_1'])
            # output_grating.rotate(90)


    for i in range(0, np.size(D_list)):
        dev = D << D_list[i]
        dev.connect(port = 'optical', destination = bs_tree_ref.ports[i+3])
        dev.movex(taper_length)
        if i % 2 == 1:
            dev.movex(stagger_neighbors)
        D << pr.route_basic(port1 = bs_tree_ref.ports[i+3], port2 = dev.ports['optical'], layer = lys['wg_deep'])

    return D


###############################################################################
#
# PACKAGING
#
###############################################################################

def rectangle_cut_corners(size = (100,100), corner_cut = 70, layer = lys['su8']):
    D = Device()
    xpts = [0,0,corner_cut, size[0]-corner_cut, size[0], size[0], size[0]-corner_cut, corner_cut]
    ypts = [corner_cut, size[1]-corner_cut, size[1], size[1], size[1]-corner_cut, corner_cut, 0, 0]
    p = [xpts, ypts]
    D.add_polygon(p, layer = layer)

    return D

def gluebox(size = (400,400), corner_cut = 70, wall_width = 100, layer = lys['su8']):
    D = Device()
    G = Device()
    size_big_rect = (size[0]+2*wall_width, size[1]+2*wall_width)
    outerbox = G << rectangle_cut_corners(size = size_big_rect, corner_cut = corner_cut)
    innerbox = G << rectangle_cut_corners(size = size, corner_cut = corner_cut*size[0]/size_big_rect[0])
    innerbox.center = outerbox.center
    gluebox = D << pg.boolean(outerbox, innerbox, 'a-b', layer = layer)

    return D


###############################################################################
#
# JJ / JMS
#
###############################################################################

#from nc_library__vt_jj import jj_circle, jj_4wire, jj_series_array, jj_2stitch
#from nc_library__vt_sq import sq_washer, sq_input_coil_one_layer, sq_two_layer, sq_8wire, sq_inductance_test
#from nc_library__vt_util import vt_arg_helper, corner_fixer, vt_fill, vt_fill_circles, vt_endpoint_boxes, vt_label_maker, vt_inline_test, vt_litho_tests
#from nc_library__vt_pads_vias_wires import pad_locations, jj_pad, jj_pad_variable, jj_pad_variable_multi, via_general, jj_jj1_m4_via, via_2stitch, via_series_array, wire_tri_seg
#from nc_library__vt_res import vt_resistor_test_structures, res_stitch, res_series_array

###############################################################################
#
# end JJ / JMS
#
###############################################################################