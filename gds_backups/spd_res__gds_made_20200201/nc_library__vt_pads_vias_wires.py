import copy
import numpy as np
import gdspy
import phidl, phidl.geometry as pg, phidl.routing as pr
from phidl import Device, Layer, LayerSet, Port, quickplot2 as qp
from phidl import make_device
from nc_library import wire_basic
from vt_util import vt_layers
vt_lyrs,layer_data = vt_layers()
from nc_library__vt_util import vt_arg_helper, vt_label_maker, corner_fixer


def jj_pad(params = dict()):
    
    pad_opening_inset = vt_arg_helper(params,'pad_opening_inset',15)    
    pad_w_wire = vt_arg_helper(params,'pad_w_wire',20)
    pad_size = vt_arg_helper(params,'pad_size',[207,344])
    pad_label_size = vt_arg_helper(params,'pad_label_size',10)
    pad_gp_gap = vt_arg_helper(params,'pad_gp_gap',15)
    pad_pkg_gap = vt_arg_helper(params,'pad_pkg_gap',25)
    is_ground_pad = vt_arg_helper(params,'is_ground_pad',False)
    via_width = vt_arg_helper(params,'via_width',2)
    pad_collar_w_wall = vt_arg_helper(params,'pad_collar_w_wall',47)
    res_pad_olap = vt_arg_helper(params,'res_pad_olap',4.7)
    layers = vt_arg_helper(params,'layers',vt_lyrs)    
    if is_ground_pad == True:
        pad_size = vt_arg_helper(params,'pad_size_ground',[20,20])
    elif is_ground_pad == False:
        pad_size = vt_arg_helper(params,'pad_size',[250,350])
    
    D_pad = Device('pad')
    
    #ground plane
    if is_ground_pad == False:
        rect1 = pg.rectangle(size = [pad_size[0]+2*pad_opening_inset,pad_size[1]+2*pad_opening_inset], layer = layers['m2'])
        rect2 = pg.outline(rect1, distance = pad_gp_gap, precision = 1e-6, layer = layers['m2'])
        rect3 = D_pad.add_ref(rect2)
        rect3.center = [0,0]
        #SU8 barriers        
        rect4 = pg.rectangle(size = [pad_size[0]+2*(pad_opening_inset+pad_gp_gap+pad_pkg_gap),pad_size[1]+2*(pad_opening_inset+pad_gp_gap+pad_pkg_gap)], layer = layers['pkg'])
        rect5 = pg.outline(rect4, distance = pad_collar_w_wall, precision = 1e-6, layer = layers['pkg'])
        rect6 = D_pad.add_ref(rect5)
        rect6.center = [0,0]        
    
    #metal layers
    layer_list = ['m1p','stfp','r1p','jj1','m3p','r2p','v4']
    for ii in range(len(layer_list)):
        rect1 = pg.rectangle(size = [pad_size[0]-2*ii*pad_opening_inset,pad_size[1]-2*ii*pad_opening_inset], layer = layers[layer_list[ii]])        
        rect1.center = [0,0]
        D_pad.add_ref(rect1)
        if layer_list[ii] == 'r1p':
            tn = res_pad_olap
        elif layer_list[ii] != 'r1p':
            tn = 0
        port_name = layer_list[ii]+'_north'
        D_pad.add_port(name = port_name, midpoint = [0,pad_size[1]/2-ii*pad_opening_inset-tn], width = pad_w_wire, orientation = 90)
        port_name = layer_list[ii]+'_south'
        D_pad.add_port(name = port_name, midpoint = [0,-pad_size[1]/2+ii*pad_opening_inset+tn], width = pad_w_wire, orientation = 270)
        port_name = layer_list[ii]+'_east'
        D_pad.add_port(name = port_name, midpoint = [pad_size[0]/2-ii*pad_opening_inset-tn,0], width = pad_w_wire, orientation = 0)
        port_name = layer_list[ii]+'_west'
        D_pad.add_port(name = port_name, midpoint = [-pad_size[0]/2+ii*pad_opening_inset+tn,0], width = pad_w_wire, orientation = 180)
        
    #vias
    D_vias = Device('pad_vias')
    via_layer_list = ['v1','v2','v3']
    num_vias = len(via_layer_list)
    via_buffer = 2
    l_unit_cell_x = 2*(via_buffer+via_width)
    l_unit_cell_y = num_vias*(via_buffer+via_width)
    p_s_min = [pad_size[0]-(len(layer_list)-2)*2*pad_opening_inset,pad_size[1]-(len(layer_list)-2)*2*pad_opening_inset]
    n_x = np.floor(p_s_min[0]/l_unit_cell_x)
    n_y = np.floor(p_s_min[1]/l_unit_cell_y)
    
    Unit_cell = Device('pad_vias__unit_cell')
    for ii in range(num_vias):
        Rect1 = pg.rectangle(size = [via_width,via_width], layer = layers[via_layer_list[ii]])
        rect1 = Unit_cell.add_ref(Rect1)
        rect1.movey(ii*(via_buffer+via_width))
        if (ii % 2) != 0:
            rect1.movex(via_buffer+via_width)
        
    Row = Device('pad_vias__row')
    for ii in range(n_x.astype(int)):
        rect1 = Row.add_ref(Unit_cell)
        rect1.movex(ii*l_unit_cell_x)
            
    for ii in range(n_y.astype(int)):
        rect1 = D_vias.add_ref(Row)
        if ii % 2 != 0: 
            rect1.reflect(p1 = (0,0), p2 = (0,1))
        rect1.center = [0,0]
        rect1.movey(ii*l_unit_cell_y)
        
    D_vias.center = [0,0]
    
    if is_ground_pad == False:
        D_pad.add_port(name = 'text', midpoint = [-pad_size[0]/2-pad_opening_inset-pad_gp_gap-pad_label_size,-pad_size[1]/4], width = 0, orientation = 180)
    else:
        D_pad.add_port(name = 'text', midpoint = [-pad_size[0]/2-pad_label_size,-pad_size[1]/4], width = 0, orientation = 180)
    
    D_pad.add_ref(D_vias)
    
    D_pad.add_port(name = 'pad_anchor', midpoint = D_pad.center, width = 0, orientation = 90)
    
    return D_pad


def jj_pad_variable(params = dict()):
        
    pad_size = vt_arg_helper(params,'pad_size',[200,250])
    pad_opening_inset = vt_arg_helper(params,'pad_opening_inset',15)    
    pad_w_wire = vt_arg_helper(params,'pad_w_wire',20)
    pad_variable_bottom_layer = vt_arg_helper(params,'pad_variable_bottom_layer','jj1')
    pad_variable_opening_layer = vt_arg_helper(params,'pad_variable_opening_layer','v2')
    layers = vt_arg_helper(params,'layers',vt_lyrs) 
    
    D_pad = Device('pad')
    
    rect1 = pg.rectangle(size = pad_size, layer = layers[pad_variable_bottom_layer])
    rect2 = pg.rectangle(size = [pad_size[0]-2*pad_opening_inset,pad_size[1]-2*pad_opening_inset], layer = layers[pad_variable_opening_layer])
    rect2.center = rect1.center
    D_pad.add_ref(rect1)
    D_pad.add_ref(rect2)
    D_pad.add_port(name = 'north', midpoint = [pad_size[0]/2,pad_size[1]], width = pad_w_wire, orientation = 90)
    D_pad.add_port(name = 'south', midpoint = [pad_size[0]/2,0], width = pad_w_wire, orientation = 270)
    D_pad.add_port(name = 'east', midpoint = [pad_size[0],pad_size[1]/2], width = pad_w_wire, orientation = 0)
    D_pad.add_port(name = 'west', midpoint = [0,pad_size[1]/2], width = pad_w_wire, orientation = 180)
    D_pad.add_port(name = 'north_west', midpoint = [0,pad_size[1]], width = pad_w_wire, orientation = 90)
    D_pad.add_port(name = 'south_west', midpoint = [0,0], width = pad_w_wire, orientation = 270)
    D_pad.add_port(name = 'text', midpoint = [0,0], width = 0, orientation = 180)
    
    return D_pad


def jj_pad_variable_multi(params = dict()):
        
    pad_size = vt_arg_helper(params,'pad_size',[200,250])
    pad_opening_inset = vt_arg_helper(params,'pad_opening_inset',15)    
    pad_w_wire = vt_arg_helper(params,'pad_w_wire',20)
    pad_variable_layers = vt_arg_helper(params,'pad_variable_layers',['m3','r2','v4'])
    layers = vt_arg_helper(params,'layers',vt_lyrs) 
    
    D_pad = Device('pad')
    
    for ii in range(len(pad_variable_layers)):
        rect1 = pg.rectangle(size = [pad_size[0]-2*ii*pad_opening_inset,pad_size[1]-2*ii*pad_opening_inset], layer = layers[pad_variable_layers[ii]])
        rect1.center = [0,0]
        D_pad.add_ref(rect1)
    D_pad.add_port(name = 'north', midpoint = [0,pad_size[1]/2], width = pad_w_wire, orientation = 90)
    D_pad.add_port(name = 'north_west', midpoint = [-pad_size[0]/2,pad_size[1]/2], width = pad_w_wire, orientation = 90)
    D_pad.add_port(name = 'south', midpoint = [0,-pad_size[1]/2], width = pad_w_wire, orientation = 270)
    D_pad.add_port(name = 'south_west', midpoint = [-pad_size[0]/2,-pad_size[1]/2], width = pad_w_wire, orientation = 270)
    D_pad.add_port(name = 'east', midpoint = [pad_size[0]/2,0], width = pad_w_wire, orientation = 0)
    D_pad.add_port(name = 'west', midpoint = [-pad_size[0]/2,0], width = pad_w_wire, orientation = 180)
    D_pad.add_port(name = 'text', midpoint = [-pad_size[0]/2,-pad_size[1]/2], width = 0, orientation = 180)
    
    return D_pad


def jj_jj1_m3_via(via_width = 4, jj_bottom_contact_outset = 1, jj_top_contact_outset = 0.6, w_wire = 1, layers = vt_lyrs,**kwargs):
    
    D_via = Device('via')
    
    w_m3 = via_width+2*jj_top_contact_outset
    sq_m3 = pg.rectangle(size = [w_m3,w_m3], layer = layers['m3'])
    sq_via = pg.rectangle(size = [via_width,via_width], layer = layers['v3'])
    w_jj1 = via_width+2*jj_bottom_contact_outset
    sq_jj1 = pg.rectangle(size = [w_jj1,w_jj1], layer = layers['jj1'])
    sq_jj1.center = [0,0]
    sq_m3.center = sq_jj1.center
    sq_via.center = sq_jj1.center
    
    D_via.add_ref(sq_m3)
    D_via.add_ref(sq_via)
    D_via.add_ref(sq_jj1)
    
    D_via.add_port(name = 'jj1_north', midpoint = [0,w_jj1/2], width = w_wire, orientation = 90)
    D_via.add_port(name = 'jj1_east', midpoint = [w_jj1/2,0], width = w_wire, orientation = 0)
    D_via.add_port(name = 'jj1_south', midpoint = [0,-w_jj1/2], width = w_wire, orientation = 270)
    D_via.add_port(name = 'jj1_west', midpoint = [-w_jj1/2,0], width = w_wire, orientation = 180)
    
    D_via.add_port(name = 'm3_north', midpoint = [0,w_m3/2], width = w_wire, orientation = 90)
    D_via.add_port(name = 'm3_east', midpoint = [w_m3/2,0], width = w_wire, orientation = 0)
    D_via.add_port(name = 'm3_south', midpoint = [0,-w_m3/2], width = w_wire, orientation = 270)
    D_via.add_port(name = 'm3_west', midpoint = [-w_m3/2,0], width = w_wire, orientation = 180)
    
    return D_via


def via_general(params = dict()):
    
    via_width = vt_arg_helper(params,'via_width',4)
    via_bottom_contact_outset = vt_arg_helper(params,'via_bottom_contact_outset',1)
    via_top_contact_outset = vt_arg_helper(params,'via_top_contact_outset',0.6)
    via_layer = vt_arg_helper(params,'via_layer','v2')
    layer_below = vt_arg_helper(params,'layer_below','m2o')
    layer_above = vt_arg_helper(params,'layer_above','jj1')
    layers = vt_arg_helper(params,'layers',vt_lyrs)
    
    D_via = Device('via')
    
    w_above = via_width+2*via_top_contact_outset
    sq_above = pg.rectangle(size = [w_above,w_above], layer = layers[layer_above])
    sq_via = pg.rectangle(size = [via_width,via_width], layer = layers[via_layer])
    w_below = via_width+2*via_bottom_contact_outset
    sq_below = pg.rectangle(size = [w_below,w_below], layer = layers[layer_below])
    sq_below.center = [0,0]
    sq_above.center = sq_below.center
    sq_via.center = sq_below.center
    D_via.add_ref(sq_above)
    D_via.add_ref(sq_via)
    D_via.add_ref(sq_below)
    if layer_below == 'jj2':
        w_below2 = via_width+4*via_bottom_contact_outset
        sq_below2 = pg.rectangle(size = [w_below2,w_below2], layer = layers['jj1'])
        sq_below2.center = sq_via.center
        D_via.add_ref(sq_below2)        
        
    w_wire = via_width+2*via_bottom_contact_outset    
    D_via.add_port(name = 'below_north', midpoint = [0,w_below/2], width = w_wire, orientation = 90)
    D_via.add_port(name = 'below_east', midpoint = [w_below/2,0], width = w_wire, orientation = 0)
    D_via.add_port(name = 'below_south', midpoint = [0,-w_below/2], width = w_wire, orientation = 270)
    D_via.add_port(name = 'below_west', midpoint = [-w_below/2,0], width = w_wire, orientation = 180)
    
    w_wire = via_width+2*via_top_contact_outset    
    D_via.add_port(name = 'above_north', midpoint = [0,w_above/2], width = w_wire, orientation = 90)
    D_via.add_port(name = 'above_east', midpoint = [w_above/2,0], width = w_wire, orientation = 0)
    D_via.add_port(name = 'above_south', midpoint = [0,-w_above/2], width = w_wire, orientation = 270)
    D_via.add_port(name = 'above_west', midpoint = [-w_above/2,0], width = w_wire, orientation = 180)
    
    return D_via


def via_2stitch(params = dict()):
    
    inter_via_gap = vt_arg_helper(params,'inter_via_gap',1)
    inter_via_gap = vt_arg_helper(params,'inter_via_gap',1)
    w_wire = vt_arg_helper(params,'w_wire',1)
    layers = vt_arg_helper(params,'layers',vt_lyrs)
    layer_below = vt_arg_helper(params,'layer_below','jj1')
    layer_above = vt_arg_helper(params,'layer_above','m3') 
    layers = vt_arg_helper(params,'layers',vt_lyrs) 
    
    D_vias = Device('via_2stitch')
                      
    via = via_general(params)
    pitch = via.ysize+inter_via_gap
    via1 = D_vias.add_ref(via)
    via2 = D_vias.add_ref(via)
    via2.movey(pitch)
    wire = wire_basic(via1.ports['below_north'],via2.ports['below_south'],'y',w_wire,layers[layer_below])
    D_vias.add_ref(wire)
    if layer_below == 'jj2':
        wireb = wire_basic(via1.ports['below_north'],via2.ports['below_south'],'y',w_wire,layers['jj1'])
        wireb = D_vias.add_ref(wireb)
    D_vias.add_port(name = 'below_north', midpoint = via1.ports['below_north'].midpoint, width = w_wire, orientation = 90)
    D_vias.add_port(name = 'below_south', midpoint = via1.ports['below_south'].midpoint, width = w_wire, orientation = 270)
    wire_length = (via2.ports['above_south'].midpoint[1]-via1.ports['above_north'].midpoint[1])/2
    via1_wire = pg.rectangle(size = [w_wire,wire_length], layer = layers[layer_above])
    via1_wire.add_port(name = 'wire_north', midpoint = [w_wire/2,wire_length], width = w_wire, orientation = 90)
    via1_wire.add_port(name = 'wire_south', midpoint = [w_wire/2,0], width = w_wire, orientation = 270)
    wire1 = D_vias.add_ref(via1_wire)
    wire1.connect(port = 'wire_north', destination = via1.ports['above_south'])
    wire2 = D_vias.add_ref(via1_wire)
    wire2.connect(port = 'wire_south', destination = via2.ports['above_north'])
    
    D_vias.add_port(name = 'north', midpoint = wire2.ports['wire_north'].midpoint, width = w_wire, orientation = 90)
    D_vias.add_port(name = 'south', midpoint = wire1.ports['wire_south'].midpoint, width = w_wire, orientation = 270)
        
    return D_vias


def via_series_array(params = dict()):
    
    via_width = vt_arg_helper(params,'via_width',4)
    inter_via_gap = vt_arg_helper(params,'inter_via_gap',1)
    num_vias_per_col = vt_arg_helper(params,'num_vias_per_col',20)
    num_cols_vias = vt_arg_helper(params,'num_cols_vias',10)
    w_wire = vt_arg_helper(params,'w_wire',1)
    layers = vt_arg_helper(params,'layers',vt_lyrs)
    layer_below = vt_arg_helper(params,'layer_below','jj1')
    via_layer = vt_arg_helper(params,'via_layer','v2')
    layer_above = vt_arg_helper(params,'layer_above','m3')    
    pad_size = vt_arg_helper(params,'pad_size',[200,250])
    pad_size_ground = vt_arg_helper(params,'pad_size_ground',[44,33])
    pad_w_wire = vt_arg_helper(params,'pad_w_wire',20)
    pad_label_size = vt_arg_helper(params,'pad_label_size',10)
    pad_pitch = vt_arg_helper(params,'pad_pitch',[200,375])
    pad_y_backset = vt_arg_helper(params,'pad_y_backset',100)
    pad_gnd_y_backset = vt_arg_helper(params,'pad_gnd_y_backset',50)
    device_label_size = vt_arg_helper(params,'device_label_size',100)
    label_layer = vt_arg_helper(params,'label_layer','m3l')
    layers = vt_arg_helper(params,'layers',vt_lyrs)        
    
    D_vias = Device('via_series')
                      
    vias = via_2stitch(params)
    Cf = corner_fixer(w_wire,layers[layer_above])
    
    Column = Device('via_sub_series')
    Column.add_port(name = 'south', midpoint = vias.ports['south'].midpoint, width = w_wire, orientation = 270)    
    for ii in range(round(num_vias_per_col/2)):
        vias_inst1 = Column.add_ref(vias)
        vias_inst1.movey(ii*vias.ysize)    
    Column.add_port(name = 'north', midpoint = vias_inst1.ports['north'].midpoint, width = w_wire, orientation = 90)

    cf1 = Column.add_ref(Cf)
    cf1.connect(port = 'south', destination = Column.ports['north'])
    Column.add_port(name = 'connector_north_north', midpoint = cf1.ports['north'].midpoint, width = w_wire, orientation = 90)
    Column.add_port(name = 'connector_north_east', midpoint = cf1.ports['east'].midpoint, width = w_wire, orientation = 0)
    Column.add_port(name = 'connector_north_west', midpoint = cf1.ports['west'].midpoint, width = w_wire, orientation = 180)
    cf2 = Column.add_ref(Cf)
    cf2.connect(port = 'north', destination = Column.ports['south'])
    Column.add_port(name = 'connector_south_south', midpoint = cf2.ports['south'].midpoint, width = w_wire, orientation = 270)
    Column.add_port(name = 'connector_south_east', midpoint = cf2.ports['east'].midpoint, width = w_wire, orientation = 0)
    Column.add_port(name = 'connector_south_west', midpoint = cf2.ports['west'].midpoint, width = w_wire, orientation = 180)
    
    D_vias.add_port(name = 'south_west_south', midpoint = cf2.ports['south'].midpoint, width = w_wire, orientation = 270)
    D_vias.add_port(name = 'south_west_west', midpoint = cf2.ports['west'].midpoint, width = w_wire, orientation = 180)
    pitch = vias.xsize+inter_via_gap    
    col_inst1 = D_vias.add_ref(Column)
    for ii in range(num_cols_vias-1):
        col_inst2 = D_vias.add_ref(Column)
        col_inst2.movex((ii+1)*pitch)
        if ((ii+1)%2) == 0:            
            wire_inst = wire_basic(col_inst1.ports['connector_south_east'],col_inst2.ports['connector_south_west'],'x',w_wire,layers[layer_above])
        else:
            wire_inst = wire_basic(col_inst1.ports['connector_north_east'],col_inst2.ports['connector_north_west'],'x',w_wire,layers[layer_above])
        D_vias.add_ref(wire_inst)
        col_inst1 = col_inst2
    D_vias.add_port(name = 'south_east_south', midpoint = col_inst2.ports['connector_south_south'].midpoint, width = w_wire, orientation = 270)
    D_vias.add_port(name = 'south_east_east', midpoint = col_inst2.ports['connector_south_east'].midpoint, width = w_wire, orientation = 0)
    
    text_coords = [D_vias.xmin-4*w_wire,D_vias.y]
    vias_width = D_vias.xsize
   
    #pads       
    params['is_ground_pad'] = False
    D_pad = jj_pad(params)
    params['is_ground_pad'] = True
    D_pad_gnd = jj_pad(params)
   
    port_label = layer_above+'_north' 
    port_label_gnd = layer_above+'_west'
    wire_layer = layer_above
    if layer_above == 'm2' or layer_above == 'm2o':
        port_label = 'm1_north'
        port_label_gnd = 'm1_west'
        wire_layer = 'm1'
        Via = via_general(params)
        via1 = D_vias.add_ref(Via)
        via1.connect(port = 'above_north', destination = D_vias.ports['south_west_south'])
        via2 = D_vias.add_ref(Via)
        via2.connect(port = 'above_north', destination = D_vias.ports['south_east_south'])
    pad1 = D_vias.add_ref(D_pad)
    pad1.center = [D_vias.ports['south_west_west'].midpoint[0]-vias_width/2,D_vias.ports['south_west_west'].midpoint[1]-pad_size[1]/2-pad_y_backset]
    if layer_above != 'm2o':                                                                         
        wire1 = wire_tri_seg(p1 = pad1.ports[port_label],p2 = D_vias.ports['south_west_west'],initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.5,0.3], directions = 'yx', layer = layers[wire_layer])
    elif layer_above == 'm2o':
        wire1 = wire_tri_seg(p1 = pad1.ports[port_label],p2 = via1.ports['below_west'],initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.5,0.3], directions = 'yx', layer = layers[wire_layer])
         
    pad2 = D_vias.add_ref(D_pad_gnd)
    if layer_above != 'm2o':
        pad2.center = [D_vias.ports['south_east_east'].midpoint[0]+pad_size_ground[0]/2+pad_gnd_y_backset,D_vias.ports['south_east_east'].midpoint[1]]
        wire2 = wire_tri_seg(p1 = pad2.ports[port_label_gnd],p2 = D_vias.ports['south_east_east'],initial_width = pad_w_wire, final_width = w_wire, length_factors = [4,4,0.5], directions = 'xx', layer = layers[wire_layer])
    elif layer_above == 'm2o':
        pad2.center = [via2.ports['below_east'].midpoint[0]+pad_size_ground[0]/2+pad_gnd_y_backset,via2.ports['below_east'].midpoint[1]]
        wire2 = wire_tri_seg(p1 = pad2.ports[port_label_gnd],p2 = via2.ports['below_east'],initial_width = pad_w_wire, final_width = w_wire, length_factors = [4,4,0.5], directions = 'xx', layer = layers[wire_layer])
    
    D_vias.add_ref(wire1)
    D_vias.add_ref(wire2)
    
    D_vias.add_port(name = 'pad_anchor', midpoint = [pad1.ports['m3_north'].midpoint[0],pad1.ports['m3_north'].midpoint[1]-pad_size[1]/2], width = pad_size[0], orientation = 90)
            
    #label
    text_string = 'via series\n'+layer_below+' / '+via_layer+' / '+layer_above+'\nw = '+str(via_width)+'\nnum_vias = '+str(num_cols_vias*num_vias_per_col)  
    Text_label = vt_label_maker(text_string = text_string, text_size = device_label_size, justify = 'right', layer = layers[label_layer])
    text_label = D_vias.add_ref(Text_label)
    text_label.center = [text_coords[0]-text_label.xsize/2,text_coords[1]]
#    
#    #pad label
#    Text_label = vt_label_maker(text_string = 'gnd', text_size = pad_label_size, justify = 'left', layer = layers[label_layer])
#    text_label = D_vias.add_ref(Text_label)
#    text_label.rotate(90)
#    text_label.center = pad2.ports['text'].midpoint
    
    return D_vias


def vt_m1_v3_via_sequence(params = dict()):
    
    w_wire = vt_arg_helper(params,'w_wire',2.02)
    layers = vt_arg_helper(params,'layers',vt_lyrs)    
    
    D_vias = Device('m1_m3_via_sequence')
        
    params_mod = copy.deepcopy(params)
    params_mod['via_layer'] = 'v1'
    params_mod['layer_below'] = 'm1'
    params_mod['layer_above'] = 'm2o'
    via1 = D_vias.add_ref(via_general(params_mod))  
    params_mod['via_layer'] = 'v2'
    params_mod['layer_below'] = 'm2o'
    params_mod['layer_above'] = 'jj1'
    via2 = D_vias.add_ref(via_general(params_mod)) 
    params_mod['via_layer'] = 'v3'
    params_mod['layer_below'] = 'jj1'
    params_mod['layer_above'] = 'm3'
    via3 = D_vias.add_ref(via_general(params_mod))
    via2.connect(port = 'below_west', destination = via1.ports['above_east'])
    via2.movex(w_wire)
    D_vias.add_ref(wire_basic(via1.ports['above_east'],via2.ports['below_west'],'x',w_wire,layers['m2o']))
    via3.connect(port = 'below_west', destination = via2.ports['above_east'])
    via3.movex(w_wire)
    D_vias.add_ref(wire_basic(via2.ports['above_east'],via3.ports['below_west'],'x',w_wire,layers['jj1']))

    D_vias.add_port(name = 'm1_west', midpoint = via1.ports['below_west'].midpoint, width = w_wire, orientation = 180)
    D_vias.add_port(name = 'm1_south', midpoint = via1.ports['below_south'].midpoint, width = w_wire, orientation = 270)
    D_vias.add_port(name = 'm1_north', midpoint = via1.ports['below_north'].midpoint, width = w_wire, orientation = 90)
    D_vias.add_port(name = 'm3_east', midpoint = via3.ports['above_east'].midpoint, width = w_wire, orientation = 0)
    D_vias.add_port(name = 'm3_south', midpoint = via3.ports['above_south'].midpoint, width = w_wire, orientation = 270)
    D_vias.add_port(name = 'm3_north', midpoint = via3.ports['above_north'].midpoint, width = w_wire, orientation = 90)
    
    return D_vias


def vt_m1_jj1_via_sequence(params = dict()):
    
    w_wire = vt_arg_helper(params,'w_wire',2.02)
    layers = vt_arg_helper(params,'layers',vt_lyrs)    
    
    D_vias = Device('m1_jj1_via_sequence')
        
    params_mod = copy.deepcopy(params)
    params_mod['via_layer'] = 'v1'
    params_mod['layer_below'] = 'm1'
    params_mod['layer_above'] = 'm2o'
    via1 = D_vias.add_ref(via_general(params_mod))  
    params_mod['via_layer'] = 'v2'
    params_mod['layer_below'] = 'm2o'
    params_mod['layer_above'] = 'jj1'
    via2 = D_vias.add_ref(via_general(params_mod)) 
    via2.connect(port = 'below_west', destination = via1.ports['above_east'])
    via2.movex(w_wire)
    D_vias.add_ref(wire_basic(via1.ports['above_east'],via2.ports['below_west'],'x',w_wire,layers['m2o']))

    D_vias.add_port(name = 'm1_west', midpoint = via1.ports['below_west'].midpoint, width = w_wire, orientation = 180)
    D_vias.add_port(name = 'm1_south', midpoint = via1.ports['below_south'].midpoint, width = w_wire, orientation = 270)
    D_vias.add_port(name = 'm1_north', midpoint = via1.ports['below_north'].midpoint, width = w_wire, orientation = 90)
    D_vias.add_port(name = 'jj1_east', midpoint = via2.ports['above_east'].midpoint, width = w_wire, orientation = 0)
    D_vias.add_port(name = 'jj1_south', midpoint = via2.ports['above_south'].midpoint, width = w_wire, orientation = 270)
    D_vias.add_port(name = 'jj1_north', midpoint = via2.ports['above_north'].midpoint, width = w_wire, orientation = 90)
    
    return D_vias


def wire_tri_seg(p1, p2, initial_width = 20, final_width = 1, length_factors = [0.8,0.1], directions = 'yx', layer = 240):

    D_wire = Device(name = 'wire_tri_seg')

    if isinstance(p1, phidl.device_layout.Port): p1 = p1.midpoint
    if isinstance(p2, phidl.device_layout.Port): p2 = p2.midpoint
    
    if directions == 'yx' or directions == 'yy':
        
        if p2[0] > p1[0] and p2[1] >= p1[1]:
            case_label = 'regular_x'
            p1_x = p1[0]
            p2_x = p2[0]
            p1_y = p1[1]
            p2_y = p2[1]
        if p2[0] > p1[0] and p2[1] <= p1[1]:
            case_label = 'mirror_y'
            p1_x = p1[0]
            p2_x = p2[0]
            p1_y = p2[1]
            p2_y = p1[1]
        elif p2[0] < p1[0] and p2[1] >= p1[1]:
            case_label = 'mirror_x'
            p1_x = p1[0]
            p2_x = p1[0]+(p1[0]-p2[0])
            p1_y = p1[1]
            p2_y = p2[1]
        elif p2[0] == p1[0]:
            case_label = 'vertical'
            p1_x = p1[0]
            p2_x = p2[0]
            p1_y = p1[1]
            p2_y = p2[1]
        elif p2[0] < p1[0] and p2[1] <= p1[1]:
            case_label = 'rotation_180'
            p1_x = p2[0]
            p2_x = p1[0]
            p1_y = p2[1]
            p2_y = p1[1]            
    
    if directions == 'xx' or directions == 'xy':
        
        if p2[0] > p1[0] and p2[1] > p1[1]:
            case_label = 'regular_x'
            p1_x = p1[0]
            p2_x = p2[0]
            p1_y = p1[1]
            p2_y = p2[1]
        if p2[0] > p1[0] and p2[1] <= p1[1]:
            case_label = 'mirror_y'
            p1_x = p1[0]
            p2_x = p2[0]
            p1_y = p2[1]
            p2_y = p1[1]
        elif p2[0] < p1[0] and p2[1] >= p1[1]:
            case_label = 'mirror_x'
            p1_x = p1[0]
            p2_x = p1[0]+(p1[0]-p2[0])
            p1_y = p1[1]
            p2_y = p2[1]
        elif p2[1] == p1[1]:
            case_label = 'horizontal'
            p1_x = p1[0]
            p2_x = p2[0]
            p1_y = p1[1]
            p2_y = p2[1]
        elif p2[0] < p1[0] and p2[1] <= p1[1]:
            case_label = 'rotation_180'
            p1_x = p2[0]
            p2_x = p1[0]
            p1_y = p2[1]
            p2_y = p1[1]            
        
    total_length_x = p2_x-p1_x
    total_length_y = p2_y-p1_y
    
    w1 = initial_width
    w2 = final_width
    if directions == 'yx':
        
        len1_x = (1-length_factors[1])*total_length_x
        len2_x = length_factors[1]*total_length_x
        len1_y = length_factors[0]*total_length_y
        len2_y = (1-length_factors[0])*total_length_y
                
        x0 = p1_x
        x1 = x0+w1/2
        x2 = x1
        x3 = x2-w1/2+len1_x
        x4 = x3+len2_x
        x5 = x4
        x6 = x3-w2
        x7 = x0-w1/2
        x8 = x7
        
        y0 = p1_y
        y1 = y0
        y2 = y1+len1_y
        y3 = y2+len2_y-w2/2
        y4 = y3
        y5 = y4+w2
        y6 = y5
        y7 = y2+w1
        y8 = y0
            
        D_wire.add_polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4),(x5,y5),(x6,y6),(x7,y7),(x8,y8)],layer)
        
    elif directions == 'yy':
        
        total_length_factor = np.sum(length_factors)
        length_factors = length_factors/total_length_factor
        len1_y = length_factors[0]*total_length_y
        len2_y = length_factors[1]*total_length_y
        len3_y = length_factors[2]*total_length_y
                
        if p1_x-w1/2 <= p2_x and p1_x+w1/2 >= p2_x:
            
            x0 = p1_x
            x1 = x0+w1/2
            x2 = x1
            x3 = p2_x+w2/2
            x4 = x3
            x5 = p2_x-w2/2
            x6 = x5
            x7 = x0-w1/2
            x8 = x7
            
            y0 = p1_y
            y1 = y0
            y2 = y1+len1_y
            y3 = y2+len2_y#-w2
            y4 = y3+len3_y#+w2
            y5 = y4
            y6 = y5-len3_y
            y7 = y2#+w1
            y8 = y1
            
            D_wire.add_polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4),(x5,y5),(x6,y6),(x7,y7),(x8,y8)],layer)
            
        else:
            
            x0 = p1_x
            x1 = x0+w1/2
            x2 = x1
            x3 = p2_x+w2/2
            x4 = x3
            x5 = p2_x-w2/2
            x6 = x5
            x7 = x0-w1/2
            x8 = x7
            
            y0 = p1_y
            y1 = y0
            y2 = y1+len1_y
            y3 = y2+len2_y-w2
            y4 = y3+w2+len3_y
            y5 = y4
            y6 = y5-len3_y
            y7 = y2+w1
            y8 = y1
            
            D_wire.add_polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4),(x5,y5),(x6,y6),(x7,y7),(x8,y8)],layer)
            
    elif directions == 'xy':
        
        len1_x = length_factors[0]*total_length_x
        len2_x = (1-length_factors[0])*total_length_x
        len1_y = length_factors[0]*total_length_y
        len2_y = (1-length_factors[0])*total_length_y
                
        x0 = p1_x
        x1 = x0
        x2 = x1+len1_x+w1
        x3 = x0+len1_x+len2_x+w2/2
        x4 = x3
        x5 = x4-w2
        x6 = x5
        x7 = x0+len1_x
        x8 = x0
        
        y0 = p1_y
        y1 = y0-w1/2
        y2 = y1
        y4 = y0+len1_y+len2_y
        y3 = y4-len2_y-w2
        y5 = y4
        y6 = y5-len2_y
        y7 = y0+w1/2
        y8 = y7
            
        D_wire.add_polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4),(x5,y5),(x6,y6),(x7,y7),(x8,y8)],layer)
        
    elif directions == 'xx':
        
        total_length_factor = np.sum(length_factors)
        length_factors = length_factors/total_length_factor
        len1_x = length_factors[0]*total_length_x
        len2_x = length_factors[1]*total_length_x
        len3_x = length_factors[2]*total_length_x
                
        if p1_y-w1/2 <= p2_y and p1_y+w1/2 >= p2_y:
            
            x0 = p1_x
            x1 = x0
            x2 = x1+len1_x
            x3 = x2+len2_x
            x4 = x0+total_length_x
            x5 = x4
            x6 = x5-len3_x
            x7 = x0+len1_x
            x8 = x0
            
            y0 = p1_y
            y1 = y0-w1/2
            y2 = y1
            y3 = p2_y-w2/2
            y4 = y3
            y5 = y4+w2
            y6 = y5
            y7 = p1_y+w1/2
            y8 = y7
            
            D_wire.add_polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4),(x5,y5),(x6,y6),(x7,y7),(x8,y8)],layer) 
            
        else:
            
            x0 = p1_x
            x1 = x0
            x2 = x1+len1_x+w1
            x3 = p2_x-len3_x
            x4 = p2_x
            x5 = x4
            x6 = x5-len3_x-w2
            x7 = x0+len1_x
            x8 = x0
            
            y0 = p1_y
            y1 = y0-w1/2
            y2 = y1
            y3 = p2_y-w2/2
            y4 = y3
            y5 = y4+w2
            y6 = y5
            y7 = p1_y+w1/2
            y8 = y7
            
            D_wire.add_polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4),(x5,y5),(x6,y6),(x7,y7),(x8,y8)],layer)
        
    if case_label == 'mirror_x':
        D_wire.reflect(p1 = [x0,y0], p2 = [x0,y0+1])
    if case_label == 'mirror_y':
        D_wire.reflect(p1 = [x0,y0], p2 = [x0+1,y0])
        D_wire.movey(total_length_y)
        
    if case_label == 'rotation_180':
        D_wire.rotate(180)
        if directions == 'yx':
            D_wire.move([p2[0]-D_wire.xmin,p1[1]-D_wire.ymax])
        elif directions == 'yy':
            D_wire.move([p2[0]-D_wire.xmin-final_width/2,p1[1]-D_wire.ymax])
            
    return D_wire


def pad_locations(params = dict()):
    
    chip_size = vt_arg_helper(params,'chip_size',[5000,5000])
    pad_x_inset = vt_arg_helper(params,'pad_x_inset',600)
    pad_y_inset = vt_arg_helper(params,'pad_y_inset',25)
    pad_pitch = vt_arg_helper(params,'pad_pitch',[200,375])
    pad_size = vt_arg_helper(params,'pad_size',[180,220])
    pad_size_ground = vt_arg_helper(params,'pad_size_ground',[53,67])
    pad_size_ground_master = vt_arg_helper(params,'pad_size_ground_master',[411,412])
    ground_master_inset = vt_arg_helper(params,'ground_master_inset',[411,412])
    pad_row_spacing_vec = vt_arg_helper(params,'pad_row_spacing_vec',[pad_pitch[1],pad_pitch[1],pad_pitch[1],pad_pitch[1]])
    die_type = vt_arg_helper(params,'die_type','vt02')
    layers = vt_arg_helper(params,'layers',vt_lyrs)
    
    D_pads = Device('pad_locations')  
#    chip_edge = pg.rectangle(size = chip_size, layer = layers['pl']) 

#    ce = D_pads.add_ref(chip_edge)
#    ce.center = [0,0]
    cwx = chip_size[0]/2
    cwy = chip_size[1]/2
    
    marker_box = pg.rectangle(size = [pad_size[0],pad_size[1]], layer = layers['pl'])
    marker_box_ground = pg.rectangle(size = [pad_size_ground[0],pad_size_ground[1]], layer = layers['pl'])
    marker_box_ground_master = pg.rectangle(size = [pad_size_ground_master[0],pad_size_ground_master[1]], layer = layers['pl'])
    
    if die_type == 'c3':
        for ii in range(10):
            
            #pad south a
            temp_str_s_a = 'p_s_a_'+str(ii+1)
            coords_s_a = [-cwx+pad_x_inset+2*ii*pad_pitch[0]+pad_size[0]/2,-cwy+pad_y_inset+pad_size[1]/2]
            D_pads.add_port(name = temp_str_s_a, midpoint = coords_s_a, width = pad_size[0], orientation = 270)
            mb = D_pads.add_ref(marker_box)
            mb.center = coords_s_a
            
            #pad east a
            temp_str_e_a = 'p_e_a_'+str(ii+1)
            coords_e_a = [cwx-pad_y_inset-pad_size[1]/2,-cwy+pad_x_inset+2*ii*pad_pitch[0]+pad_size[0]/2]
            D_pads.add_port(name = temp_str_e_a, midpoint = coords_e_a, width = pad_size[0], orientation = 0)
            mb = D_pads.add_ref(marker_box)
            mb.rotate(90)
            mb.center = coords_e_a
            
            #pad north a
            temp_str_n_a = 'p_n_a_'+str(ii+1)
            coords_n_a = [cwx-pad_x_inset-pad_size[0]/2-2*ii*pad_pitch[0],cwy-pad_y_inset-pad_size[1]/2]
            D_pads.add_port(name = temp_str_n_a, midpoint = coords_n_a, width = pad_size[0], orientation = 90)
            mb = D_pads.add_ref(marker_box)
            mb.center = coords_n_a
            
            #pad west a
            temp_str_w_a = 'p_w_a_'+str(ii+1)
            coords_w_a = [-cwx+pad_y_inset+pad_size[1]/2,cwy-pad_x_inset-2*ii*pad_pitch[0]-pad_size[0]/2]
            D_pads.add_port(name = temp_str_w_a, midpoint = coords_w_a, width = pad_size[0], orientation = 180)
            mb = D_pads.add_ref(marker_box)
            mb.rotate(90)
            mb.center = coords_w_a           
            
        for ii in range(9):
            
            #pad south b
            temp_str_s_b = 'p_s_b_'+str(ii+1)
            coords_s_b = [-cwx+pad_x_inset+(2*ii+1)*pad_pitch[0]+pad_size[0]/2,-cwy+pad_y_inset+pad_size[1]/2+pad_pitch[1]]
            D_pads.add_port(name = temp_str_s_b, midpoint = coords_s_b, width = pad_size[0], orientation = 270)
            mb = D_pads.add_ref(marker_box)
            mb.center = coords_s_b
            
            #pad east b
            temp_str_e_b = 'p_e_b_'+str(ii+1)
            coords_e_a = [cwx-pad_y_inset-pad_size[1]/2-pad_pitch[1],-cwy+pad_x_inset+(2*ii+1)*pad_pitch[0]+pad_size[0]/2]
            D_pads.add_port(name = temp_str_e_b, midpoint = coords_e_a, width = pad_size[0], orientation = 0)
            mb = D_pads.add_ref(marker_box)
            mb.rotate(90)
            mb.center = coords_e_a
            
            #pad north b
            temp_str_n_b = 'p_n_b_'+str(ii+1)
            coords_n_b = [cwx-pad_x_inset-pad_size[0]/2-(2*ii+1)*pad_pitch[0],cwy-pad_y_inset-pad_size[1]/2-pad_pitch[1]]
            D_pads.add_port(name = temp_str_n_b, midpoint = coords_n_b, width = pad_size[0], orientation = 90)
            mb = D_pads.add_ref(marker_box)
            mb.center = coords_n_b
            
            #pad west b
            temp_str_w_b = 'p_w_b_'+str(ii+1)
            coords_w_b = [-cwx+pad_y_inset+pad_size[1]/2+pad_pitch[1],cwy-pad_x_inset-(2*ii+1)*pad_pitch[0]-pad_size[0]/2]
            D_pads.add_port(name = temp_str_w_b, midpoint = coords_w_b, width = pad_size[0], orientation = 180)
            mb = D_pads.add_ref(marker_box)
            mb.rotate(90)
            mb.center = coords_w_b
    
    if die_type == 'vt02_perimeter':
        for ii in range(10):
            
            #pad south a
            temp_str_s_a = 'p_s_a_'+str(ii+1)
            coords_s_a = [-cwx+pad_x_inset+ii*pad_pitch[0]+pad_size[0]/2,-cwy+pad_y_inset+pad_size[1]/2]
            D_pads.add_port(name = temp_str_s_a, midpoint = coords_s_a, width = pad_size[0], orientation = 270)
            mb = D_pads.add_ref(marker_box)
            mb.center = coords_s_a
            
            #pad east a
            temp_str_e_a = 'p_e_a_'+str(ii+1)
            coords_e_a = [cwx-pad_y_inset-pad_size[1]/2,-cwy+pad_x_inset+ii*pad_pitch[0]+pad_size[0]/2]
            D_pads.add_port(name = temp_str_e_a, midpoint = coords_e_a, width = pad_size[0], orientation = 0)
            mb = D_pads.add_ref(marker_box)
            mb.rotate(90)
            mb.center = coords_e_a
            
            #pad north a
            temp_str_n_a = 'p_n_a_'+str(ii+1)
            coords_n_a = [cwx-pad_x_inset-pad_size[0]/2-ii*pad_pitch[0],cwy-pad_y_inset-pad_size[1]/2]
            D_pads.add_port(name = temp_str_n_a, midpoint = coords_n_a, width = pad_size[0], orientation = 90)
            mb = D_pads.add_ref(marker_box)
            mb.center = coords_n_a
            
            #pad west a
            temp_str_w_a = 'p_w_a_'+str(ii+1)
            coords_w_a = [-cwx+pad_y_inset+pad_size[1]/2,cwy-pad_x_inset-ii*pad_pitch[0]-pad_size[0]/2]
            D_pads.add_port(name = temp_str_w_a, midpoint = coords_w_a, width = pad_size[0], orientation = 180)
            mb = D_pads.add_ref(marker_box)
            mb.rotate(90)
            mb.center = coords_w_a           
    
            
    if die_type == 'vt02_two-sides':
        str_list = ['a','b','c']
        aa = np.cumsum(pad_row_spacing_vec[:])
        bb = np.cumsum(pad_row_spacing_vec[2:])
        for ii in range(9):
            for jj in range(3):
                
                #south
                temp_str_s_a = 'p_s_'+str_list[jj]+'_'+str(ii+1)
#                coords_s_a = [-cwx+x_inset+ii*pad_pitch[0]+pad_size[0]/2,-cwy+y_inset+pad_size[1]/2+jj*pad_pitch[1]]
                if jj > 0:
                    tn_s = aa[jj-1]
                else:
                    tn_s = 0
                coords_s_a = [-cwx+pad_x_inset+ii*pad_pitch[0]+pad_size[0]/2,-cwy+pad_y_inset+pad_size[1]/2+tn_s]
                D_pads.add_port(name = temp_str_s_a, midpoint = coords_s_a, width = pad_size[0], orientation = 270)
                mb = D_pads.add_ref(marker_box)
                mb.center = coords_s_a
                            
                #north
                if jj > 0:
                    tn_n = bb[jj-1]
                else:
                    tn_n = 0
                temp_str_n_a = 'p_n_'+str_list[jj]+'_'+str(ii+1)
                coords_n_a = [cwx-pad_x_inset-ii*pad_pitch[0]-pad_size[0]/2,cwy-pad_y_inset-pad_size[1]/2-tn_n]
                D_pads.add_port(name = temp_str_n_a, midpoint = coords_n_a, width = pad_size[0], orientation = 90)
                mb = D_pads.add_ref(marker_box)
                mb.center = coords_n_a
                        
    return D_pads


def vt_inductor(params = dict()):

    induct_w_wire = vt_arg_helper(params,'induct_w_wire',2.4)
    induct_wire_pitch = vt_arg_helper(params,'induct_wire_pitch',4.4)
    induct_num_squares = vt_arg_helper(params,'induct_num_squares',245)
    induct_l_pre = vt_arg_helper(params,'induct_l_pre',9.7)
    induct_turn_ratio = vt_arg_helper(params,'induct_turn_ratio',2.7)
    induct_layer = vt_arg_helper(params,'induct_layer','m3')
    induct_include_inductex_ports = vt_arg_helper(params,'induct_include_inductex_ports',False)
    layers = vt_arg_helper(params,'layers',vt_lyrs)
    
    D_ind = Device('inductor')
    
    meander = D_ind.add_ref(pg.snspd(wire_width = induct_w_wire, wire_pitch = induct_wire_pitch, size = None, num_squares = induct_num_squares, turn_ratio = induct_turn_ratio, terminals_same_side = False, layer = layers[induct_layer]))
    D_ind.add_ref(wire_basic([meander.ports[1].midpoint[0]-induct_l_pre,meander.ports[1].midpoint[1]],[meander.ports[1].midpoint[0],meander.ports[1].midpoint[1]],'x',induct_w_wire,layers[induct_layer]))
    D_ind.add_ref(wire_basic([meander.ports[2].midpoint[0],meander.ports[2].midpoint[1]],[meander.ports[2].midpoint[0]+induct_l_pre,meander.ports[2].midpoint[1]],'x',induct_w_wire,layers[induct_layer]))
    
    D_ind.add_port(name = 'east', midpoint = [meander.ports[2].midpoint[0]+induct_l_pre,meander.ports[2].midpoint[1]], width = induct_w_wire, orientation = 0)
    D_ind.add_port(name = 'west', midpoint = [meander.ports[1].midpoint[0]-induct_l_pre,meander.ports[1].midpoint[1]], width = induct_w_wire, orientation = 180)
    
    #inductor ports
    if induct_include_inductex_ports == True:
        if induct_layer == 'm1':
            lyr = 'ipm1'
        if induct_layer == 'jj1':
            lyr = 'ipj1'
        if induct_layer == 'm3':
            lyr = 'ipm3'
        Inductor_port = pg.rectangle(size = [induct_w_wire,2], layer = layers[lyr])
        in_po1 = D_ind.add_ref(Inductor_port)
        in_po2 = D_ind.add_ref(Inductor_port)
        in_po1.center = D_ind.ports['west'].midpoint
        in_po2.center = D_ind.ports['east'].midpoint
        D_ind.label(text = 'P1+ '+induct_layer, position = D_ind.ports['west'].midpoint, layer = layers['ipl'])
        D_ind.label(text = 'P1- '+induct_layer, position = D_ind.ports['east'].midpoint, layer = layers['ipl'])
    
    return D_ind