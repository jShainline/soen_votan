import copy
import numpy as np
import gdspy
import phidl, phidl.geometry as pg
from phidl import Device, Layer, LayerSet, Port
from phidl import make_device

from nc_library import wire_basic
from vt_util import vt_layers

from nc_library__vt_util import vt_arg_helper, corner_fixer, vt_label_maker
from nc_library__vt_pads_vias_wires import jj_pad, jj_pad_variable_multi, wire_tri_seg
#from nc_library__vt_res import res_stitch

vt_lyrs,layer_data = vt_layers()

def vt_resistor_meander(params = dict()):
    
    res_layer = vt_arg_helper(params,'res_layer','r1')
    res_w_wire = vt_arg_helper(params,'res_w_wire',2)
    res_num_squares_meander = vt_arg_helper(params,'res_num_squares_meander',5000)
    pad_size = vt_arg_helper(params,'pad_size',[200,250])
    pad_y_backset = vt_arg_helper(params,'pad_y_backset',100)
    pad_w_wire = vt_arg_helper(params,'pad_w_wire',20)
    device_label_size = vt_arg_helper(params,'device_label_size',5)
    label_layer = vt_arg_helper(params,'label_layer','m4l')
    layers = vt_arg_helper(params,'layers',vt_lyrs)
    
    D_res = Device('resistor_meander') 
    
    #res meander
    wire_pitch = 2*res_w_wire
    meander = D_res.add_ref(pg.snspd(wire_width = res_w_wire, wire_pitch = wire_pitch, size = None, num_squares = res_num_squares_meander, turn_ratio = 4, terminals_same_side = False, layer = layers[res_layer])).rotate(90)

    #pads and wires    
    params_mod = copy.deepcopy(params)   
    Pad = jj_pad(params_mod)
    params_mod['is_ground_pad'] = True    
    Pad_gnd = jj_pad(params_mod)
    
    pad1 = D_res.add_ref(Pad)
    pad1.connect(port = 'm1p_north', destination = meander.ports[1])
    pad1.move([0,-pad_y_backset])   
    pad_port = res_layer+'_north'
    if res_layer == 'r1' or res_layer == 'r2':
        pad_port = 'pad_anchor'
    if res_layer == 'm2i':
        pad_port = 'm1p_north'
    wire1 = wire_tri_seg(p1 = pad1.ports[pad_port], p2 = meander.ports[1], initial_width = pad_w_wire, final_width = res_w_wire, length_factors = [0.75,0.1,0.1], directions = 'yy', layer = layers[res_layer])
    D_res.add_ref(wire1)
    
    pad2 = D_res.add_ref(Pad_gnd)
    pad2.connect(port = 'm1p_south', destination = meander.ports[2])
    pad2.move([0,pad_y_backset/2])
    pad_port = res_layer+'_south'
    if res_layer == 'r1' or res_layer == 'r2':
        pad_port = 'pad_anchor'
    if res_layer == 'm2i':
        pad_port = 'm1p_south'
    wire2 = wire_tri_seg(p1 = pad2.ports[pad_port], p2 = meander.ports[2], initial_width = pad_w_wire, final_width = res_w_wire, length_factors = [0.75,0.1,0.1], directions = 'yy', layer = layers[res_layer])
    D_res.add_ref(wire2)
    
    #text label
    string = 'resistor meander\nw_wire = '+str(res_w_wire)+' um\npitch = '+str(wire_pitch)+' um\nnum_sq = '+str(int(round(res_num_squares_meander)))+'\nlayer = '+res_layer
    Text_label = vt_label_maker(text_string = string, text_size = device_label_size, justify = 'right', layer = layers[label_layer])
    text_label = D_res.add_ref(Text_label)
    text_label.center = ([meander.xmin-text_label.xsize/2-10,meander.y])
    
    #pad anchor
    D_res.add_port(name = 'pad_anchor', midpoint = [pad1.ports['m1p_north'].midpoint[0],pad1.ports['m1p_west'].midpoint[1]], width = pad_size[0], orientation = 90)
        
    return D_res

def vt_resistor_4_wire(params = dict()):
    
    res_layer = vt_arg_helper(params,'res_layer','r1')
    res_w_wire = vt_arg_helper(params,'res_w_wire',2)
    res_resistance = vt_arg_helper(params,'res_resistance',50)
    res_per_sq_r1 = vt_arg_helper(params,'res_per_sq_r1',0.1)
    res_per_sq_r2 = vt_arg_helper(params,'res_per_sq_r2',2)
    pad_size = vt_arg_helper(params,'pad_size',[200,250])
    pad_size_ground = vt_arg_helper(params,'pad_size_ground',[200,250])
    pad_pitch = vt_arg_helper(params,'pad_pitch',[200,250])
    pad_y_backset = vt_arg_helper(params,'pad_y_backset',100)
    pad_w_wire = vt_arg_helper(params,'pad_w_wire',20)
    device_label_size = vt_arg_helper(params,'device_label_size',5)
    label_layer = vt_arg_helper(params,'label_layer','m4l')
    layers = vt_arg_helper(params,'layers',vt_lyrs)
    
    D_res = Device('resistor_4_wire') 
        
    #connection points
    Cf = corner_fixer(res_w_wire,layers[res_layer])
    cf_left = D_res.add_ref(Cf)
    cf_right = D_res.add_ref(Cf)
    if res_layer == 'r1':
        res_per_sq = res_per_sq_r1
    elif res_layer == 'r2':
        res_per_sq = res_per_sq_r2        
    res_num_squares = res_resistance/res_per_sq
    cf_right.movex(res_num_squares*res_w_wire+res_w_wire)
    
    #resistor
    wire1 = wire_basic(cf_left.ports['east'],cf_right.ports['west'],'x',res_w_wire,layers[res_layer])            
    D_res.add_ref(wire1)
    
    #pads 
    pad_y_backset_local = 100
    params_mod = copy.deepcopy(params)   
    Pad = jj_pad(params_mod)
    params_mod['is_ground_pad'] = True    
    Pad_gnd = jj_pad(params_mod)   
    
    pad_I_in = D_res.add_ref(Pad)
    pad_I_gnd = D_res.add_ref(Pad_gnd)
    pad_V_in = D_res.add_ref(Pad)
    pad_V_gnd = D_res.add_ref(Pad_gnd)
    
    pad_V_in.connect(port = 'm1p_north', destination = cf_left.ports['south'])
    pad_V_in.move([pad_size[0]/2,-pad_y_backset_local])
    
    pad_I_in.connect(port = 'm1p_east', destination = pad_V_in.ports['m1p_west'])
    pad_I_in.movex(-pad_pitch[0]+pad_size[0])
    
    pad_V_gnd.connect(port = 'm1p_south', destination = cf_right.ports['north'])
    pad_V_gnd.move([-pad_size_ground[0]/2,2*pad_size_ground[1]])
    
    pad_I_gnd.connect(port = 'm1p_west', destination = pad_V_gnd.ports['m1p_east'])
    pad_I_gnd.movex(pad_size_ground[0])
    
    #wires
    wire_I_in = wire_tri_seg(p1 = pad_I_in.ports['r1p_north'], p2 = cf_left.ports['west'], initial_width = pad_w_wire, final_width = res_w_wire, length_factors = [0.75,0.05], directions = 'yx', layer = layers[res_layer])
    D_res.add_ref(wire_I_in)
    wire_V_in = wire_tri_seg(p1 = pad_V_in.ports['r1p_north'], p2 = cf_left.ports['south'], initial_width = pad_w_wire, final_width = res_w_wire, length_factors = [0.4,0.3,0.05], directions = 'yy', layer = layers[res_layer])
    D_res.add_ref(wire_V_in)
    wire_V_gnd = wire_tri_seg(p1 = pad_V_gnd.ports['r1p_south'], p2 = cf_right.ports['north'], initial_width = pad_w_wire, final_width = res_w_wire, length_factors = [0.3,0.3,0.05], directions = 'yy', layer = layers[res_layer])
    D_res.add_ref(wire_V_gnd)    
    wire_I_gnd = wire_tri_seg(p1 = pad_I_gnd.ports['r1p_south'], p2 = cf_right.ports['east'], initial_width = pad_w_wire, final_width = res_w_wire, length_factors = [0.3,0.1], directions = 'yx', layer = layers[res_layer])
    D_res.add_ref(wire_I_gnd)
    

    #text label
    string = 'resistor 4 wire\nw_wire = '+str(res_w_wire)+'um\nnum_squares = '+str(np.around(res_num_squares,decimals = 2))+'\nlayer = '+res_layer
    Text_label = vt_label_maker(text_string = string, text_size = device_label_size, layer = layers[label_layer])
    text_label = D_res.add_ref(Text_label)
    text_label.connect(port = 'south', destination = pad_V_gnd.ports['m1p_north'])
    
    #pad anchor
    D_res.add_port(name = 'pad_anchor', midpoint = [pad_I_in.x,pad_I_in.y], width = pad_size[0], orientation = 90)
    
    return D_res


def vt_res_stitch_4_wire(params = dict()):
    
    res_layer = vt_arg_helper(params,'res_layer','r1')
    res_w_wire = vt_arg_helper(params,'res_w_wire',2)
    res_resistance = vt_arg_helper(params,'res_resistance',50)
    res_per_sq_r1 = vt_arg_helper(params,'res_per_sq_r1',0.1)
    res_per_sq_r2 = vt_arg_helper(params,'res_per_sq_r2',2)
    pad_size = vt_arg_helper(params,'pad_size',[200,250])
    pad_size_ground = vt_arg_helper(params,'pad_size_ground',[200,250])
    pad_pitch = vt_arg_helper(params,'pad_pitch',[200,250])
    pad_w_wire = vt_arg_helper(params,'pad_w_wire',20)
    device_label_size = vt_arg_helper(params,'device_label_size',5)
    label_layer = vt_arg_helper(params,'label_layer','m4l')
    layers = vt_arg_helper(params,'layers',vt_lyrs)
    
    D_res = Device('resistor_4_wire') 
        
    #connection points
    if res_layer == 'r1':
        res_per_sq = res_per_sq_r1
        metal_layer = 'm1'
        metal_layer_pad = 'm1p'
    elif res_layer == 'r2':
        res_per_sq = res_per_sq_r2
        metal_layer = 'm3'
        metal_layer_pad = 'm3p'
    Cf = corner_fixer(res_w_wire,layers[metal_layer])
    cf_left = D_res.add_ref(Cf)
    cf_right = D_res.add_ref(Cf)
    res_num_squares = res_resistance/res_per_sq
#    cf_right.movex(res_num_squares*res_w_wire+res_w_wire)
    
    #resistor
    params_mod = copy.deepcopy(params)
    params_mod['res_include_label'] = False
    Res = res_stitch_simp(params_mod)
    res = D_res.add_ref(Res)
    res.connect(port = 'west', destination = cf_left.ports['east'])
    cf_right.connect(port = 'west', destination = res.ports['east'])
    
    #pads 
    pad_y_backset_local = 100
    params_mod = copy.deepcopy(params)   
    Pad = jj_pad(params_mod)
    params_mod['is_ground_pad'] = True    
    Pad_gnd = jj_pad(params_mod)   
    
    pad_I_in = D_res.add_ref(Pad)
    pad_I_gnd = D_res.add_ref(Pad_gnd)
    pad_V_in = D_res.add_ref(Pad)
    pad_V_gnd = D_res.add_ref(Pad_gnd)
    
    pad_V_in.connect(port = 'm1p_north', destination = cf_left.ports['south'])
    pad_V_in.move([pad_size[0]/2,-pad_y_backset_local])
    
    pad_I_in.connect(port = 'm1p_east', destination = pad_V_in.ports['m1p_west'])
    pad_I_in.movex(-pad_pitch[0]+pad_size[0])
    
    pad_V_gnd.connect(port = 'm1p_south', destination = cf_right.ports['north'])
    pad_V_gnd.move([-pad_size_ground[0]/2,2*pad_size_ground[1]])
    
    pad_I_gnd.connect(port = 'm1p_west', destination = pad_V_gnd.ports['m1p_east'])
    pad_I_gnd.movex(pad_size_ground[0])
    
    #wires
    wire_I_in = wire_tri_seg(p1 = pad_I_in.ports[metal_layer_pad+'_north'], p2 = cf_left.ports['west'], initial_width = pad_w_wire, final_width = res_w_wire, length_factors = [0.75,0.05], directions = 'yx', layer = layers[metal_layer])
    D_res.add_ref(wire_I_in)
    wire_V_in = wire_tri_seg(p1 = pad_V_in.ports[metal_layer_pad+'_north'], p2 = cf_left.ports['south'], initial_width = pad_w_wire, final_width = res_w_wire, length_factors = [0.4,0.3,0.05], directions = 'yy', layer = layers[metal_layer])
    D_res.add_ref(wire_V_in)
    wire_V_gnd = wire_tri_seg(p1 = pad_V_gnd.ports[metal_layer_pad+'_south'], p2 = cf_right.ports['north'], initial_width = pad_w_wire, final_width = res_w_wire, length_factors = [0.3,0.3,0.05], directions = 'yy', layer = layers[metal_layer])
    D_res.add_ref(wire_V_gnd)    
    wire_I_gnd = wire_tri_seg(p1 = pad_I_gnd.ports[metal_layer_pad+'_south'], p2 = cf_right.ports['east'], initial_width = pad_w_wire, final_width = res_w_wire, length_factors = [0.3,0.1], directions = 'yx', layer = layers[metal_layer])
    D_res.add_ref(wire_I_gnd)
    

    #text label
    string = 'resistor 4 wire\nw_wire = '+str(res_w_wire)+'um\nnum_squares = '+str(np.around(res_num_squares,decimals = 2))+'\nlayer = '+res_layer
    Text_label = vt_label_maker(text_string = string, text_size = device_label_size, layer = layers[label_layer])
    text_label = D_res.add_ref(Text_label)
    text_label.connect(port = 'south', destination = pad_V_gnd.ports['m1p_north'])
    
    #pad anchor
    D_res.add_port(name = 'pad_anchor', midpoint = [pad_I_in.x,pad_I_in.y], width = pad_size[0], orientation = 90)
    
    return D_res


def vt_resistor_series_array(params = dict()):
        
    res_num_squares = vt_arg_helper(params,'res_num_squares',5)
    res_w_wire = vt_arg_helper(params,'res_w_wire',2)
    res_num_in_column = vt_arg_helper(params,'res_num_in_column',10)
    res_num_columns = vt_arg_helper(params,'res_num_columns',10)
    res_space = vt_arg_helper(params,'res_space',2)
    w_wire = vt_arg_helper(params,'w_wire',1)
    pad_size = vt_arg_helper(params,'pad_size',[150,200])
    pad_size_ground = vt_arg_helper(params,'pad_size_ground',[15,20])
    pad_pitch = vt_arg_helper(params,'pad_pitch',[343,445])
    pad_w_wire = vt_arg_helper(params,'pad_w_wire',20)
    pad_y_backset = vt_arg_helper(params,'pad_y_backset',20)
    pad_gnd_y_backset = vt_arg_helper(params,'pad_gnd_y_backset',23)
    device_label_size = vt_arg_helper(params,'device_label_size',10)
    label_layer = vt_arg_helper(params,'label_layer','m4l')
    layers = vt_arg_helper(params,'layers',vt_lyrs)
    
    D_res = Device('resistor_series_array')
    
    Cf = corner_fixer(w_wire,layers['m3'])
    
    Res_col = Device('resistor_column')
    
    params_mod = copy.deepcopy(params)
    params_mod['res_lower_layer'] = 'm3'
    params_mod['res_stitch_drop_via_hack'] = True
    Res = res_stitch(params_mod)
    y_pitch = Res.ysize+res_space
    num = np.round(res_num_in_column/2)
    for ii in range(num.astype(int)):
        res1 = Res_col.add_ref(Res)
        res2 = Res_col.add_ref(Res)
        res1.movey(2*ii*y_pitch)
        res2.movey((2*ii+1)*y_pitch)
        wire1 = wire_basic(res1.ports['m3_east_north'],res2.ports['m3_east_south'],'y',w_wire,layers['m3'])
        wire2 = wire_basic(res2.ports['m3_west_north'],[res2.ports['m3_west_north'].midpoint[0],res2.ports['m3_west_north'].midpoint[1]+wire1.ysize],'y',w_wire,layers['m3'])
        Res_col.add_ref(wire1)
        Res_col.add_ref(wire2)
        if ii == 0:
            Res_col.add_port(name = 'south', midpoint = res1.ports['m3_west'].midpoint, width = w_wire, orientation = 180)
        if ii == num-1:
            cf = Res_col.add_ref(Cf)
            cf.connect(port = 'south',destination = res2.ports['m3_west_north'])
            cf.movey(wire1.ysize)
            Res_col.add_port(name = 'north', midpoint = cf.ports['east'].midpoint, width = w_wire, orientation = 0)
    
    Res_cols = Device('resistor_columns')
    x_pitch = Res_col.xsize+res_space
    num = np.round(res_num_columns/2)
    for ii in range(num.astype(int)):
        col1 = Res_cols.add_ref(Res_col)
        col2 = Res_cols.add_ref(Res_col)
        col1.movex(2*ii*x_pitch)
        col2.reflect([0,0],[1,0])
        col2.move([(2*ii+1)*x_pitch,col1.ports['north'].midpoint[1]-col2.ports['south'].midpoint[1]])
        wire1 = wire_basic(col1.ports['north'],col2.ports['south'],'x',w_wire,layers['m3'])
        wire2 = wire_basic(col2.ports['north'],[col2.ports['north'].midpoint[0]+wire1.xsize,col2.ports['north'].midpoint[1]],'x',w_wire,layers['m3'])
        Res_cols.add_ref(wire1)
        Res_cols.add_ref(wire2)
        if ii == 0:
            Res_cols.add_port(name = 'west', midpoint = col1.ports['south'].midpoint, width = w_wire, orientation = 180)
        if ii == num-1:
            Res_cols.add_port(name = 'east', midpoint = [col2.ports['north'].midpoint[0]+wire1.xsize,col2.ports['north'].midpoint[1]], width = w_wire, orientation = 0)
      
    D_res.add_ref(Res_cols)
    
    #pads
    x_mid = Res_cols.ports['west'].midpoint[0]+(Res_cols.ports['east'].midpoint[0]-Res_cols.ports['west'].midpoint[0])/2
    x_min = D_res.xmin
    y_mid = D_res.y
    y_min = D_res.ymin
    y_max = D_res.ymax
    x_width = Res_cols.xsize
    
    params_mod = copy.deepcopy(params)  
    D_pad = jj_pad(params_mod)    
   
    cf2 = D_res.add_ref(Cf)
    cf2.connect(port = 'east', destination = Res_cols.ports['west'])
    cf3 = D_res.add_ref(Cf)
    cf3.connect(port = 'east', destination = cf2.ports['west'])
    pad1 = D_res.add_ref(D_pad)                                                                         
    pad1.center = [x_mid,y_min-pad_size[1]/2-pad_y_backset]
    pad1.movex(-pad_pitch[0]-x_width/2-pad_w_wire/2-2*w_wire)
    D_res.add_ref(wire_tri_seg(p1 = pad1.ports['m3p_north'],p2 = cf3.ports['west'],initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.4,0.1], directions = 'yx', layer = layers['m3']))
    pad2 = D_res.add_ref(D_pad)                                                                         
    pad2.center = [pad1.x+pad_pitch[0],pad1.y]
    D_res.add_ref(wire_tri_seg(p1 = pad2.ports['m3p_north'],p2 = cf3.ports['south'],initial_width = pad_w_wire, final_width = w_wire, length_factors = [4,4,4], directions = 'yy', layer = layers['m3']))
         
    params_mod = copy.deepcopy(params)
    params_mod['is_ground_pad'] = True
    D_pad_gnd = jj_pad(params_mod)
    pad_gnd = D_res.add_ref(D_pad_gnd)
    pad_gnd.center = [Res_cols.ports['east'].midpoint[0]+pad_size_ground[0]/2+pad_gnd_y_backset,Res_cols.ports['east'].midpoint[1]]
    D_res.add_ref( wire_tri_seg(p1 = pad_gnd.ports['m3p_west'],p2 = Res_cols.ports['east'],initial_width = pad_w_wire, final_width = w_wire, length_factors = [4,4,0.5], directions = 'xx', layer = layers['m3']))
    
    #label
    D_res.add_port(name = 'label_port', midpoint = [x_mid,y_max+device_label_size], width = 0, orientation = 90)
    text_string = 'low res\nw_res = '+str(res_w_wire)+' um\nres = '+str(res_num_squares)+' sq\nnum_res_tot = '+str(res_num_columns*res_num_in_column)
    Text_label = vt_label_maker(text_string = text_string, text_size = device_label_size, justify = 'right', layer = layers[label_layer])
    text_label = D_res.add_ref(Text_label)
    text_label.connect(port = 'south', destination = D_res.ports['label_port'])
    
    #pad anchor
    D_res.add_port(name = 'pad_anchor', midpoint = [pad1.x,pad1.y], width = pad_size[0], orientation = 90)
        
    return D_res


def res_stitch(params = dict()):
    
    jj_bottom_contact_outset = vt_arg_helper(params,'jj_bottom_contact_outset',2)
    jj_top_contact_outset = vt_arg_helper(params,'jj_top_contact_outset',1)
    res_num_squares = vt_arg_helper(params,'res_num_squares',5)
    res_w_wire = vt_arg_helper(params,'res_w_wire',2)
    res_w_cntct = vt_arg_helper(params,'res_w_cntct',1)
    res_outset = vt_arg_helper(params,'res_outset',1)
    res_l_lead = vt_arg_helper(params,'res_l_lead',2)
    res_w_lead = vt_arg_helper(params,'res_w_lead',3)
    via_width = vt_arg_helper(params,'via_width',1)
    res_layer = vt_arg_helper(params,'res_layer','r2') 
    res_lower_layer = vt_arg_helper(params,'res_lower_layer','v3') 
    w_wire = vt_arg_helper(params,'w_wire',1)
    layers = vt_arg_helper(params,'layers',vt_lyrs)
    res_stitch_drop_via_hack = vt_arg_helper(params,'res_stitch_drop_via_hack',False)
    
    D_res = Device('resistor_stitch')    
    
    D_via = Device('resistor_stitch_via')
    if res_lower_layer == 'jj1':
        r0 = pg.rectangle(size = [via_width+2*jj_bottom_contact_outset,via_width+2*jj_bottom_contact_outset], layer = layers['jj1'])
        r0.center = [0,0]
        D_via.add_ref(r0)
        half_width = (via_width+2*jj_bottom_contact_outset)/2
        D_via.add_port(name = 'jj1_west', midpoint = [-half_width,0], width = w_wire, orientation = 180)
        D_via.add_port(name = 'jj1_east', midpoint = [half_width,0], width = w_wire, orientation = 0)
        D_via.add_port(name = 'jj1_north', midpoint = [0,half_width], width = w_wire, orientation = 90)
        D_via.add_port(name = 'jj1_south', midpoint = [0,-half_width], width = w_wire, orientation = 270)        
    if res_lower_layer == 'v3' or res_lower_layer == 'jj1':
        if res_stitch_drop_via_hack != True:
            r1 = pg.rectangle(size = [via_width,via_width], layer = layers['v2'])
            r1.center = [0,0]
            D_via.add_ref(r1)
        r2 = pg.rectangle(size = [via_width+2*jj_top_contact_outset,via_width+2*jj_top_contact_outset], layer = layers['m4'])        
        r2.center = [0,0]
        D_via.add_ref(r2)
        half_width_x = (via_width+2*jj_top_contact_outset)/2
        half_width_y = half_width_x
        D_via.add_port(name = 'm3_west', midpoint = [-half_width_x,0], width = w_wire, orientation = 180)
        D_via.add_port(name = 'm3_east', midpoint = [half_width_x,0], width = w_wire, orientation = 0)
        D_via.add_port(name = 'm3_north', midpoint = [0,half_width_y], width = w_wire, orientation = 90)
        D_via.add_port(name = 'm3_south', midpoint = [0,-half_width_y], width = w_wire, orientation = 270)
        res_length = res_num_squares*res_w_wire+2*(res_w_cntct+res_outset)
    if res_lower_layer == 'm3': 
        r1 = pg.rectangle(size = [res_w_cntct+w_wire+res_outset,res_w_wire-2*res_outset], layer = layers['m3'])
        r1.center = [0,0]
        D_via.add_ref(r1)
        half_width_x = (res_w_cntct+w_wire+res_outset)/2        
        half_width_y = (res_w_wire-2*res_outset)/2
        D_via.add_port(name = 'm3_west', midpoint = [-half_width_x,0], width = w_wire, orientation = 180)        
        D_via.add_port(name = 'm3_east', midpoint = [half_width_x,0], width = w_wire, orientation = 0)
        D_via.add_port(name = 'm3_north_east', midpoint = [half_width_x-res_w_cntct/2,half_width_y], width = w_wire, orientation = 90)
        D_via.add_port(name = 'm3_north_east_wire', midpoint = [half_width_x-w_wire/2,half_width_y], width = w_wire, orientation = 90)
        D_via.add_port(name = 'm3_south_east', midpoint = [half_width_x-res_w_cntct/2,-half_width_y], width = w_wire, orientation = 270)
        D_via.add_port(name = 'm3_north_west', midpoint = [-half_width_x+res_w_cntct/2,half_width_y], width = w_wire, orientation = 90)
        D_via.add_port(name = 'm3_north_west_wire', midpoint = [-half_width_x+w_wire/2,half_width_y], width = w_wire, orientation = 90)
        D_via.add_port(name = 'm3_south_west', midpoint = [-half_width_x+res_w_cntct/2,-half_width_y], width = w_wire, orientation = 270)   
        res_length = res_num_squares*res_w_wire+2*(res_w_cntct)
    
    via_1 = D_res.add_ref(D_via)
    via_2 = D_res.add_ref(D_via)
    via_2.connect(port = 'm3_west', destination = via_1.ports['m3_east'])
    if res_lower_layer == 'v3' or res_lower_layer == 'jj1':
        tn = res_length+2*res_l_lead-2*res_outset
    if res_lower_layer == 'm3':
        tn = res_length+2*res_l_lead            
    via_2.movex(tn)
    
    Res = pg.rectangle(size = [res_length,res_w_wire], layer = layers[res_layer])
    if res_lower_layer == 'v3' or res_lower_layer == 'jj1':
        tn = res_outset+res_w_cntct/2
    if res_lower_layer == 'm3':
        tn = res_w_cntct/2
    Res.add_port(name = 'west', midpoint = [tn,Res.ymax-res_outset], width = res_w_cntct, orientation = 270)
    Res.add_port(name = 'east', midpoint = [Res.xmax-tn,Res.ymax-res_outset], width = res_w_cntct, orientation = 270)
    res = D_res.add_ref(Res)
    res.center = [via_1.center[0]+(via_2.center[0]-via_1.center[1])/2,via_1.center[1]]
    
    if res_lower_layer == 'v3' or res_lower_layer == 'jj1':
        res_contact = pg.rectangle(size = [res_w_cntct,res_w_wire-2*res_outset], layer = layers['m3'])
        res_contact.add_port(name = 'north', midpoint = [res_w_cntct/2,res_contact.ymax], width = res_w_cntct, orientation = 90)
        res_contact.add_port(name = 'west', midpoint = [0,res_contact.ymin+res_contact.size[1]/2], width = w_wire, orientation = 180)
        res_contact.add_port(name = 'east', midpoint = [res_contact.xmax,res_contact.ymin+res_contact.size[1]/2], width = w_wire, orientation = 0)
        west_contact = D_res.add_ref(res_contact)
        west_contact.connect(port = 'north', destination = res.ports['west'])
        east_contact = D_res.add_ref(res_contact)
        east_contact.connect(port = 'north', destination = res.ports['east'])
        wire1 = wire_basic(via_1.ports['m3_east'],west_contact.ports['west'],'x',res_w_lead,layers['m3'])
        wire2 = wire_basic(east_contact.ports['east'],via_2.ports['m3_west'],'x',res_w_lead,layers['m3'])   
        D_res.add_ref(wire1)
        D_res.add_ref(wire2)
    elif res_lower_layer == 'm3':
        via_1.connect(port = 'm3_north_east', destination = res.ports['west'])
        via_2.connect(port = 'm3_north_west', destination = res.ports['east'])
           
    D_res.add_port(name = 'stitch_center', midpoint = [D_res.xmin+D_res.xsize/2,D_res.ymin+D_res.ysize/2], width = 0, orientation = 0)
    D_res.add_port(name = 'm3_west', midpoint = via_1.ports['m3_west'].midpoint, width = w_wire, orientation = 180)
    D_res.add_port(name = 'm3_east', midpoint = via_2.ports['m3_east'].midpoint, width = w_wire, orientation = 0)
    if res_lower_layer == 'v2' or res_lower_layer == 'jj1':
        D_res.add_port(name = 'm3_west_north', midpoint = via_1.ports['m3_north'].midpoint, width = w_wire, orientation = 90)
        D_res.add_port(name = 'm3_west_south', midpoint = via_1.ports['m3_south'].midpoint, width = w_wire, orientation = 270)        
        D_res.add_port(name = 'm3_east_north', midpoint = via_2.ports['m3_north'].midpoint, width = w_wire, orientation = 90)
        D_res.add_port(name = 'm3_east_south', midpoint = via_2.ports['m3_south'].midpoint, width = w_wire, orientation = 270)        
    if res_lower_layer == 'm3':
        D_res.add_port(name = 'm3_west_north', midpoint = via_1.ports['m3_north_west_wire'].midpoint, width = w_wire, orientation = 90)
        D_res.add_port(name = 'm3_west_south', midpoint = via_1.ports['m3_south_west'].midpoint, width = w_wire, orientation = 270)
        D_res.add_port(name = 'm3_east_north', midpoint = via_2.ports['m3_north_east_wire'].midpoint, width = w_wire, orientation = 90)
        D_res.add_port(name = 'm3_east_south', midpoint = via_2.ports['m3_south_east'].midpoint, width = w_wire, orientation = 270)
    if res_lower_layer == 'jj1':
        D_res.add_port(name = 'jj1_west', midpoint = via_1.ports['jj1_west'].midpoint, width = w_wire, orientation = 180)
        D_res.add_port(name = 'jj1_east', midpoint = via_2.ports['jj1_east'].midpoint, width = w_wire, orientation = 0)
            
    return D_res


def res_50_ohm(params = dict()):
    
    res_50_ohm_w_wire = vt_arg_helper(params,'res_50_ohm_w_wire',23)
    res_50_ohm_outset = vt_arg_helper(params,'res_50_ohm_outset',23)
    
    params_mod = copy.deepcopy(params)
    params_mod['res_layer'] = 'r2'
    params_mod['res_resistance'] = 50
    params_mod['res_w_wire'] = res_50_ohm_w_wire
    params_mod['res_outset'] = res_50_ohm_outset
    
    D_res = res_stitch_simp(params_mod)
            
    return D_res


def res_stitch_simp(params = dict()):
    
    res_resistance = vt_arg_helper(params,'res_resistance',50)
    res_w_wire = vt_arg_helper(params,'res_w_wire',2)
    res_outset = vt_arg_helper(params,'res_outset',1)
    res_wire_olap = vt_arg_helper(params,'res_wire_olap',1)
    res_l_lead = vt_arg_helper(params,'res_l_lead',2)
    res_include_label = vt_arg_helper(params,'res_include_label',True)
    w_wire = vt_arg_helper(params,'w_wire',3)
    res_layer = vt_arg_helper(params,'res_layer','r2')
    res_per_sq_r1 = vt_arg_helper(params,'res_per_sq_r1',0.1)
    res_per_sq_r2 = vt_arg_helper(params,'res_per_sq_r2',2)
    res_label_size = vt_arg_helper(params,'res_label_size',5)
    label_layer = vt_arg_helper(params,'label_layer','m3l')
    layers = vt_arg_helper(params,'layers',vt_lyrs)
    
    D_res = Device('resistor_stitch')
    
    if res_layer == 'r1':
        lead_layer = 'm1'
        res_num_squares = np.around(res_resistance/res_per_sq_r1, decimals = 2)
    if res_layer == 'r2':
        lead_layer = 'm3'
        res_num_squares = np.around(res_resistance/res_per_sq_r2, decimals = 2)
            
    res_length = res_num_squares*res_w_wire
    Res = pg.rectangle(size = [res_length,res_w_wire], layer = layers[res_layer])
    Res.add_port(name = 'east', midpoint = [Res.xmax,Res.y], width = res_w_wire, orientation = 0)
    Res.add_port(name = 'west', midpoint = [Res.xmin,Res.y], width = res_w_wire, orientation = 180)
    res = D_res.add_ref(Res)

    Tab = pg.rectangle(size = [res_wire_olap+res_outset,max(w_wire+2*res_outset,res_w_wire+2*res_outset)], layer = layers[res_layer])
    Tab.add_port(name = 'east', midpoint = [Tab.xmax,Tab.y], width = w_wire, orientation = 0)
    Tab.add_port(name = 'west', midpoint = [Tab.xmin,Tab.y], width = w_wire, orientation = 180)
    tab1 = D_res.add_ref(Tab)
    tab2 = D_res.add_ref(Tab)
    tab1.connect(port = 'east', destination = res.ports['west'])
    tab2.connect(port = 'west', destination = res.ports['east'])

    Lead = pg.rectangle(size = [res_l_lead+res_wire_olap,max(w_wire,res_w_wire)], layer = layers[lead_layer])
    Lead.add_port(name = 'east', midpoint = [Lead.xmax,Lead.y], width = w_wire, orientation = 0)
    Lead.add_port(name = 'west', midpoint = [Lead.xmin,Lead.y], width = w_wire, orientation = 180)
    lead1 = D_res.add_ref(Lead)
    lead2 = D_res.add_ref(Lead)
    lead1.connect(port = 'east', destination = tab1.ports['west'])
    lead1.movex(res_wire_olap)
    lead2.connect(port = 'west', destination = tab2.ports['east'])
    lead2.movex(-res_wire_olap)
    
    D_res.add_port(name = 'east', midpoint = lead2.ports['east'].midpoint, width = w_wire, orientation = 0)
    D_res.add_port(name = 'west', midpoint = lead1.ports['west'].midpoint, width = w_wire, orientation = 180)
    
    #text label
    if res_include_label == True:
        string = 'r = '+str(np.around(res_resistance,decimals = 3))+' ohm'
        Text_label = vt_label_maker(text_string = string, text_size = res_label_size, justify = 'center', layer = layers[label_layer])
        text_label = D_res.add_ref(Text_label)
        text_label.center = [res.x,res.ymax+res_label_size]
    
    return D_res