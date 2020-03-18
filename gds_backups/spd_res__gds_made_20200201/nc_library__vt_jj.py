import copy
import numpy as np
import gdspy
import phidl, phidl.geometry as pg
from phidl import Device, Layer, LayerSet, Port
from phidl import make_device

from nc_library import wire_basic
from vt_util import vt_layers

from nc_library__vt_util import vt_arg_helper, corner_fixer, vt_label_maker
from nc_library__vt_pads_vias_wires import jj_pad, jj_jj1_m3_via, wire_tri_seg
from nc_library__vt_res import res_50_ohm

vt_lyrs,layer_data = vt_layers()

def jj_circle(params = dict()):
    
    jj_junc_diam = vt_arg_helper(params,'jj_junc_diam',2)
    jj_num_pts = vt_arg_helper(params,'jj_num_pts',40)
    jj_via_inset = vt_arg_helper(params,'jj_via_inset',0.4)
    jj_top_contact_outset = vt_arg_helper(params,'jj_top_contact_outset',1)
    jj_bottom_contact_outset = vt_arg_helper(params,'jj_bottom_contact_outset',1.1)
    via_width = vt_arg_helper(params,'via_width',4)
    jj_include_shunt = vt_arg_helper(params,'jj_include_shunt',True)
    jj_shunt_w_wire = vt_arg_helper(params,'jj_shunt_w_wire',2)
    jj_shunt_resistance = vt_arg_helper(params,'jj_shunt_resistance',4)
    res_per_sq_r2 = vt_arg_helper(params,'res_per_sq_r2',2)
    shunt_w_cntct = vt_arg_helper(params,'shunt_w_cntct',2)
    shunt_outset = vt_arg_helper(params,'shunt_outset',1)
    shunt_l_lead = vt_arg_helper(params,'shunt_l_lead',1)
    include_fine_moats = vt_arg_helper(params,'include_fine_moats',True)
    ground_plane_moat_width_fine = vt_arg_helper(params,'ground_plane_moat_width_fine',1)
    ground_plane_buffer_fine = vt_arg_helper(params,'ground_plane_buffer_fine',2)
    jj_include_flux_moats = vt_arg_helper(params,'jj_include_flux_moats',True)
    w_wire = vt_arg_helper(params,'w_wire',1)
    layers = vt_arg_helper(params,'layers',vt_lyrs)
        
    D_jj = Device('jj')
    
    #draw shapes and center    
    angle_resolution = 360/jj_num_pts
    via_diam = jj_junc_diam-2*jj_via_inset
    sq1_w = via_diam+2*jj_top_contact_outset
    sq2_w = jj_junc_diam+2*jj_bottom_contact_outset
    lower_contact = pg.rectangle(size = [sq2_w,sq2_w], layer = layers['jj1'])
    lower_contact.add_port(name = 'east', midpoint = [sq2_w,sq2_w/2], width = w_wire, orientation = 0)
    upper_contact = pg.rectangle(size = [sq1_w,sq1_w], layer = layers['m3'])
    upper_contact.add_port(name = 'east', midpoint = [sq1_w,sq1_w/2], width = w_wire, orientation = 0)
    jj = pg.circle(radius = jj_junc_diam/2, angle_resolution = angle_resolution, layer = layers['jj2'])
    via = pg.circle(radius = via_diam/2, angle_resolution = angle_resolution, layer = layers['v3'])
    l_c = D_jj.add_ref(lower_contact)
    u_c = D_jj.add_ref(upper_contact)
    lower_contact.center = jj.center
    upper_contact.center = jj.center
    D_jj.add_ref(jj) 
    D_jj.add_ref(via)
    
    #shunt resistor
    if jj_include_shunt == True:    
        #calculate num squares
        shunt_num_squares = jj_shunt_resistance/res_per_sq_r2
        
        Shunt_partial = Device('jj_shunt_resistor')
        Shunt_full = Device('jj_shunt_assembly')
        shunt_resistor_length = shunt_num_squares*jj_shunt_w_wire+2*(shunt_w_cntct+shunt_outset)
        shunt_resistor = pg.rectangle(size = [shunt_resistor_length,jj_shunt_w_wire], layer = layers['r2'])
        shunt_resistor.add_port(name = 'west', midpoint = [shunt_outset+shunt_w_cntct/2,shunt_resistor.ymax-shunt_outset], width = shunt_w_cntct, orientation = 270)
        shunt_resistor.add_port(name = 'east', midpoint = [shunt_resistor.xmax-shunt_outset-shunt_w_cntct/2,shunt_resistor.ymax-shunt_outset], width = shunt_w_cntct, orientation = 270)
        Shunt_partial.add_ref(shunt_resistor)
        
        Shunt_contact = pg.rectangle(size = [shunt_w_cntct,jj_shunt_w_wire-2*shunt_outset], layer = layers['m3'])
        Shunt_contact.add_port(name = 'north', midpoint = [shunt_w_cntct/2,Shunt_contact.ymax], width = shunt_w_cntct, orientation = 90)
        Shunt_contact.add_port(name = 'west', midpoint = [0,Shunt_contact.ymin+Shunt_contact.size[1]/2], width = w_wire, orientation = 180)
        Shunt_contact.add_port(name = 'east', midpoint = [Shunt_contact.xmax,Shunt_contact.ymin+Shunt_contact.size[1]/2], width = w_wire, orientation = 0)
        west_contact = Shunt_partial.add_ref(Shunt_contact)
        west_contact.connect(port = 'north', destination = shunt_resistor.ports['west'])
        east_contact = Shunt_partial.add_ref(Shunt_contact)
        east_contact.connect(port = 'north', destination = shunt_resistor.ports['east'])
        
        Shunt_partial.add_port(name = 'west', midpoint = west_contact.ports['west'].midpoint, width = w_wire, orientation = 180)
        Shunt_partial.add_port(name = 'east', midpoint = east_contact.ports['east'].midpoint, width = w_wire, orientation = 0)
        shunt_partial = D_jj.add_ref(Shunt_partial)
        shunt_partial.connect(port = 'west', destination = u_c.ports['east'])
        shunt_partial.movex(shunt_l_lead+shunt_outset+(sq2_w-sq1_w)/2)
        wire1 = wire_basic(u_c.ports['east'],shunt_partial.ports['west'], 'x', w_wire, layers['m3'] )
        Shunt_full.add_ref(wire1)     
        
        Via = jj_jj1_m3_via(via_width,jj_bottom_contact_outset,jj_top_contact_outset,w_wire,layers)
        via = Shunt_full.add_ref(Via)
        via.connect(port = 'm3_west', destination = shunt_partial.ports['east'])
        via.movex(shunt_outset+shunt_l_lead+np.max([jj_bottom_contact_outset-jj_top_contact_outset,0]))
        wire2 = wire_basic(shunt_partial.ports['east'],via.ports['m3_west'], 'x', w_wire, layers['m3'] )
        Shunt_full.add_ref(wire2)
        w_wire3 = min([via.size[1],l_c.size[1]])
        wire3 = wire_basic(l_c.ports['east'],via.ports['jj1_west'], 'xyx', w_wire3, layers['jj1'] )
        Shunt_full.add_ref(wire3)  
    
        D_jj.add_ref(Shunt_full)      
    
    #add ports
    D_jj.add_port(name = 'm3_north', midpoint = [upper_contact.center[0],upper_contact.ymax], width = w_wire, orientation = 90)
    D_jj.add_port(name = 'm3_south', midpoint = [upper_contact.center[0],upper_contact.ymin], width = w_wire, orientation = 270)
    D_jj.add_port(name = 'm3_east', midpoint = [upper_contact.xmax,upper_contact.center[1]], width = w_wire, orientation = 0)
    D_jj.add_port(name = 'm3_west', midpoint = [upper_contact.xmin,upper_contact.center[1]], width = w_wire, orientation = 180)
    D_jj.add_port(name = 'jj1_north', midpoint = [lower_contact.center[0],lower_contact.ymax], width = w_wire, orientation = 90)
    D_jj.add_port(name = 'jj1_south', midpoint = [lower_contact.center[0],lower_contact.ymin], width = w_wire, orientation = 270)
    D_jj.add_port(name = 'jj1_east', midpoint = [lower_contact.xmax,lower_contact.center[1]], width = w_wire, orientation = 0)
    D_jj.add_port(name = 'jj1_west', midpoint = [lower_contact.xmin,lower_contact.center[1]], width = w_wire, orientation = 180)
    
    #flux moats
    x_min = D_jj.xmin
    x_max = D_jj.xmax
    y_min = D_jj.ymin
    y_max = D_jj.ymax
    extra_length = 1
    extra_pitch = 2
    num_extra = 3
    if include_fine_moats == True:
        if jj_include_flux_moats == True:
            for kk in range(num_extra):
                coord_pairs = [[x_min-ground_plane_buffer_fine-ground_plane_moat_width_fine/2-kk*ground_plane_moat_width_fine-kk*extra_pitch,y_min-kk*extra_length],
                               [x_min-ground_plane_buffer_fine-ground_plane_moat_width_fine/2-kk*ground_plane_moat_width_fine-kk*extra_pitch,y_max+kk*extra_length]]    
                moat_points = []
                for ii in range(len(coord_pairs)):
                    moat_points.append(np.array(coord_pairs[ii])) 
                moat_path = gdspy.PolyPath(moat_points, width = ground_plane_moat_width_fine)    
                D_jj.add_polygon(moat_path.polygons, layer = layers['m2m'])
            
        if jj_include_flux_moats == True or jj_include_flux_moats == 'only_via_side':
            for kk in range(num_extra):
                coord_pairs = [[x_max+ground_plane_buffer_fine+ground_plane_moat_width_fine/2+kk*ground_plane_moat_width_fine+kk*extra_pitch,y_min-kk*extra_length],
                               [x_max+ground_plane_buffer_fine+ground_plane_moat_width_fine/2+kk*ground_plane_moat_width_fine+kk*extra_pitch,y_max+kk*extra_length]]    
                moat_points = []
                for ii in range(len(coord_pairs)):
                    moat_points.append(np.array(coord_pairs[ii])) 
                moat_path = gdspy.PolyPath(moat_points, width = ground_plane_moat_width_fine)    
                D_jj.add_polygon(moat_path.polygons, layer = layers['m2m'])    
         
    return D_jj

def jj_4wire(params = dict()):    
    
    jj_junc_diam = vt_arg_helper(params,'jj_junc_diam',2)
    jj_shunt_resistance = vt_arg_helper(params,'jj_shunt_resistance',4)
    w_wire = vt_arg_helper(params,'w_wire',1)
    l_lead = vt_arg_helper(params,'l_lead',10)
    include_50_ohm = vt_arg_helper(params,'include_50_ohm',True)
    pad_size = vt_arg_helper(params,'pad_size',[200,250])
    pad_size_ground = vt_arg_helper(params,'pad_size_ground',[100,150])
    pad_pitch = vt_arg_helper(params,'pad_pitch',[200,375])
    pad_y_backset = vt_arg_helper(params,'pad_y_backset',100)
    pad_gnd_y_backset = vt_arg_helper(params,'pad_gnd_y_backset',100)
    pad_w_wire = vt_arg_helper(params,'pad_w_wire',20)
    pad_label_size = vt_arg_helper(params,'pad_label_size',10)
    device_label_size = vt_arg_helper(params,'device_label_size',5)
    label_layer = vt_arg_helper(params,'label_layer','m4l')
    include_big_moats = vt_arg_helper(params,'include_big_moats',False)
    ground_plane_moat_width = vt_arg_helper(params,'ground_plane_moat_width',5)
    ground_plane_buffer = vt_arg_helper(params,'ground_plane_buffer',15)
    layers = vt_arg_helper(params,'layers',vt_lyrs)
    
    D_jj = Device('jj')
            
    D_jj1 = jj_circle(params)
    D_jj.add_ref(D_jj1)
        
    #initial leads
    Upper_lead_wire = pg.rectangle(size = [w_wire,l_lead], layer = layers['m3'])
    Upper_lead_wire.add_port(name = 'top', midpoint = [0,l_lead-w_wire/2], width = w_wire, orientation = 180)
    Upper_lead_wire.add_port(name = 'bottom', midpoint = [w_wire/2,0], width = w_wire, orientation = 270)
    Upper_lead_wire.add_port(name = 'side', midpoint = [w_wire,l_lead-w_wire/2], width = w_wire, orientation = 0)
    u_wire = D_jj.add_ref(Upper_lead_wire)
    u_wire.connect(port = 'bottom', destination = D_jj1.ports['m3_north'])
    
    Lower_lead_wire = pg.rectangle(size = [w_wire,l_lead], layer = layers['jj1'])
    Lower_lead_wire.add_port(name = 'top', midpoint = [w_wire/2,l_lead], width = w_wire, orientation = 90)
    Lower_lead_wire.add_port(name = 'bottom', midpoint = [w_wire/2,0], width = w_wire, orientation = 270)
    l_wire = D_jj.add_ref(Lower_lead_wire)
    l_wire.connect(port = 'top', destination = D_jj1.ports['jj1_south'])
    
    #corner fixer instead of via
    Cf = corner_fixer(w_wire,layers['jj1'])
    cf = D_jj.add_ref(Cf)
    cf.connect(port = 'north', destination = l_wire.ports['bottom'])
    
    #50 ohm in series
    if include_50_ohm == True:
        Res = res_50_ohm(params)
        res1 = D_jj.add_ref(Res)
        res1.connect(port = 'east',destination = u_wire.ports['top'])
        move_amount = 10
        res1.movex(-move_amount)
        wire_res1 = wire_basic(res1.ports['east'],u_wire.ports['top'],'x',w_wire,layers['m3'])
        D_jj.add_ref(wire_res1)
        res2 = D_jj.add_ref(Res)
        res2.connect(port = 'west',destination = u_wire.ports['side'])
        res2.movex(move_amount)
        wire_res2 = wire_basic(res2.ports['west'],u_wire.ports['side'],'x',w_wire,layers['m3'])
        D_jj.add_ref(wire_res2)
        
    #pads and pad wiring
    params_mod = copy.deepcopy(params)
    params_mod['is_ground_pad'] = False
    D_pad = jj_pad(params_mod)
    params_mod['is_ground_pad'] = True
    D_pad_gnd = jj_pad(params_mod)
    pad_space_x = pad_pitch[0]-pad_size[0]
    x0 = D_jj1.x
    x1 = x0-pad_size[0]/2-pad_space_x/2
    y0 = D_jj1.y
    pad_y_backset_extra = 0
    y1 = y0-pad_size[1]/2-pad_y_backset-pad_y_backset_extra

    pad1 = D_jj.add_ref(D_pad)
    pad1.center = [x1,y1]
    if include_50_ohm == True:
        wire1 = wire_tri_seg(p1 = pad1.ports['m3_north'], p2 = res1.ports['west'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.75,0.1], directions = 'yx', layer = layers['m3'])
    elif include_50_ohm == False:
        wire1 = wire_tri_seg(p1 = pad1.ports['m3_north'], p2 = u_wire.ports['top'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.75,0.1], directions = 'yx', layer = layers['m3'])
    D_jj.add_ref(wire1)

    pad2 = D_jj.add_ref(D_pad)
    pad2.center = pad1.center+np.array([pad_pitch[0],0])
    if include_50_ohm == True:
        D_jj.add_ref(wire_tri_seg(p1 = pad2.ports['m3_north'], p2 = res2.ports['east'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.8,0.1], directions = 'yx', layer = layers['m3']))
    elif include_50_ohm == False:
        D_jj.add_ref(wire_tri_seg(p1 = pad2.ports['m3_north'], p2 = u_wire.ports['side'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.8,0.1], directions = 'yx', layer = layers['m3']))
    
    pad1_gnd = D_jj.add_ref(D_pad_gnd)
    pad1_gnd.center = cf.ports['south']
    pad1_gnd.movey(-pad_size_ground[1]/2-pad_gnd_y_backset)
    D_jj.add_ref(wire_tri_seg(p1 = pad1_gnd.ports['jj1_north'], p2 = cf.ports['south'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.1,0.1,0.05], directions = 'yy', layer = layers['jj1']))    

#    pad2_gnd = D_jj.add_ref(D_pad_gnd)
#    pad2_gnd.center = [x_list[3],y_list[3]]
#    D_jj.add_ref(wire_tri_seg(p1 = pad2_gnd.ports['jj1_north'], p2 = cf.ports['east'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.4,0.2], directions = 'yx', layer = layers['jj1'], layers = layers))
    
    D_jj.add_port(name = 'pad_anchor', midpoint = [pad1.ports['m3_north'].midpoint[0],pad1.ports['m3_west'].midpoint[1]], width = pad_size[0], orientation = 90)
    
    #labels
    text_coords = [u_wire.ports['top'].midpoint[0],D_jj.ymax+3*device_label_size]
    text_string = 'jj_4wire\nd = '+str(jj_junc_diam)+' um\nr = '+str(jj_shunt_resistance)+' ohm'
    Text_label = vt_label_maker(text_string = text_string, text_size = device_label_size, justify = 'center', layer = layers[label_layer])
    text_label = D_jj.add_ref(Text_label)
    text_label.center = text_coords
    
    #pad labels
#    text_strings = ['I+','I- / gnd','V+','V- / gnd']
#    text_coords = [pad1.ports['text'].midpoint,pad2.ports['text'].midpoint,pad3.ports['text'].midpoint,pad4.ports['text'].midpoint]
#    for ii in range(len(text_strings)):
#        Text_label = vt_label_maker(text_string = text_strings[ii], text_size = pad_label_size, justify = 'left', layer = layers[label_layer])
#        text_label = D_jj.add_ref(Text_label)
#        text_label.rotate(90)
#        text_label.center = text_coords[ii]
        
    #flux moats  
#    if include_big_moats == True:
#        x1 = pad1.x-pad_w_wire/2-ground_plane_buffer-ground_plane_moat_width/2
#        x2 = pad4.x+pad_w_wire/2+ground_plane_buffer+ground_plane_moat_width/2
#        y1 = pad1.ymax+ground_plane_buffer+ground_plane_moat_width/2
#        y2 = D_jj.ymax+ground_plane_moat_width/2+ground_plane_buffer
#        #outer
#    #    [pad1.xmin-ground_plane_buffer-ground_plane_moat_width/2,pad1.ymin],
#        coord_pairs = [[pad1.xmin,y1],
#                       [x1,y1],[x1,y2-18],[pad2.x,y2],[pad3.x,y2],[x2,y2-18],[x2,pad4.ymax+ground_plane_buffer+ground_plane_moat_width/2],
#                       [pad4.xmax+ground_plane_buffer+ground_plane_moat_width/2,pad4.ymax+ground_plane_buffer+ground_plane_moat_width/2],
#                       [pad4.xmax+ground_plane_buffer+ground_plane_moat_width/2,pad4.ymin-ground_plane_buffer-ground_plane_moat_width/2],        
#                       [pad3.xmax+ground_plane_buffer+ground_plane_moat_width/2,pad4.ymin-ground_plane_buffer-ground_plane_moat_width/2],
#                       [pad3.xmax+ground_plane_buffer+ground_plane_moat_width/2,pad3.ymin]]
#        moat_points = []
#        for ii in range(len(coord_pairs)):
#            moat_points.append(np.array(coord_pairs[ii])) 
#        moat_path = gdspy.PolyPath(moat_points, width = ground_plane_moat_width)    
#        D_jj.add_polygon(moat_path.polygons, layer = layers['m2m'])
#        #inner
#        coord_pairs = [[pad1.x+pad_w_wire/2+ground_plane_buffer,y1+ground_plane_moat_width/2],[pad3.x-pad_w_wire/2-ground_plane_buffer,y1+ground_plane_moat_width/2]]
#        moat_points = []
#        for ii in range(len(coord_pairs)):
#            moat_points.append(np.array(coord_pairs[ii])) 
#        moat_path = gdspy.PolyPath(moat_points, width = ground_plane_moat_width)    
#        D_jj.add_polygon(moat_path.polygons, layer = layers['m2m'])
            
    return D_jj


def jj_series_array(params = dict()):
    
    jj_junc_diam = vt_arg_helper(params,'jj_junc_diam',2.26)
    jj_include_shunt = vt_arg_helper(params,'jj_include_shunt',True)
    jj_shunt_resistance = vt_arg_helper(params,'jj_shunt_resistance',4.1)
    num_jjs_per_col = vt_arg_helper(params,'num_jjs_per_col',10)
    num_cols_jjs = vt_arg_helper(params,'num_cols_jjs',4)
    inter_jj_gap = vt_arg_helper(params,'inter_jj_gap',1)
    w_wire = vt_arg_helper(params,'w_wire',1)    
    pad_size = vt_arg_helper(params,'pad_size',[200,250])
    pad_size_ground = vt_arg_helper(params,'pad_size_ground',[100,150])
    pad_w_wire = vt_arg_helper(params,'pad_w_wire',20)
    pad_label_size = vt_arg_helper(params,'pad_label_size',10)
    pad_y_backset = vt_arg_helper(params,'pad_y_backset',100)
    device_label_size = vt_arg_helper(params,'device_label_size',100) 
    include_big_moats = vt_arg_helper(params,'include_big_moats',False)
    include_fine_moats = vt_arg_helper(params,'include_fine_moats',True)
    ground_plane_moat_width = vt_arg_helper(params,'ground_plane_moat_width',5)
    ground_plane_buffer = vt_arg_helper(params,'ground_plane_buffer',15)
    ground_plane_moat_width_fine = vt_arg_helper(params,'ground_plane_moat_width_fine',1)
    ground_plane_buffer_fine = vt_arg_helper(params,'ground_plane_buffer_fine',10)
    label_layer = vt_arg_helper(params,'label_layer','m4l')
    layers = vt_arg_helper(params,'layers',vt_lyrs)   
    
    D_jjs = Device('jj_series')
                    
    jjs = jj_2stitch(params)    
    Cf = corner_fixer(w_wire,layers['m3'])
    
    Column = Device('jj_sub_series')
    Column.add_port(name = 'south', midpoint = jjs.ports['south'].midpoint, width = w_wire, orientation = 270)
    for ii in range(round(num_jjs_per_col/2)):
        jj_inst1 = Column.add_ref(jjs)
        jj_inst1.movey(ii*jjs.ysize)
    Column.add_port(name = 'north', midpoint = jj_inst1.ports['north'].midpoint, width = w_wire, orientation = 90)        
    
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
    
    D_jjs.add_port(name = 'south_west_south', midpoint = cf2.ports['south'].midpoint, width = w_wire, orientation = 270)
    D_jjs.add_port(name = 'south_west_west', midpoint = cf2.ports['west'].midpoint, width = w_wire, orientation = 180)
    pitch = jjs.xsize+inter_jj_gap
    col_inst1 = D_jjs.add_ref(Column)
    for ii in range(num_cols_jjs-1):
        col_inst2 = D_jjs.add_ref(Column)
        col_inst2.movex((ii+1)*pitch)
        if ((ii+1)%2) == 0:            
            wire_inst = wire_basic(col_inst1.ports['connector_south_east'],col_inst2.ports['connector_south_west'],'x',w_wire,layers['m3'])
        else:
            wire_inst = wire_basic(col_inst1.ports['connector_north_east'],col_inst2.ports['connector_north_west'],'x',w_wire,layers['m3'])
        D_jjs.add_ref(wire_inst)
        col_inst1 = col_inst2
    D_jjs.add_port(name = 'south_east_south', midpoint = col_inst2.ports['connector_south_south'].midpoint, width = w_wire, orientation = 270)
    D_jjs.add_port(name = 'south_east_east', midpoint = col_inst2.ports['connector_south_east'].midpoint, width = w_wire, orientation = 0)
    jjs_size = D_jjs.size
    
    #fine flux moats
    x_min = D_jjs.xmin
    x_max = D_jjs.xmax
    y_min = D_jjs.ymin
    y_max = D_jjs.ymax
    if include_fine_moats == True:
        #left
        coord_pairs = [[x_min-ground_plane_buffer_fine-ground_plane_moat_width_fine/2,y_min+jjs.size[1]/2],
                       [x_min-ground_plane_buffer_fine-ground_plane_moat_width_fine/2,y_max]]
        moat_points = []
        for ii in range(len(coord_pairs)):
            moat_points.append(np.array(coord_pairs[ii])) 
        moat_path = gdspy.PolyPath(moat_points, width = ground_plane_moat_width_fine)    
        D_jjs.add_polygon(moat_path.polygons, layer = layers['m2m'])
        #right
        coord_pairs = [[x_max+ground_plane_buffer_fine+ground_plane_moat_width_fine/2,y_min+jjs.size[1]/2],
                       [x_max+ground_plane_buffer_fine+ground_plane_moat_width_fine/2,y_max]]
        moat_points = []
        for ii in range(len(coord_pairs)):
            moat_points.append(np.array(coord_pairs[ii])) 
        moat_path = gdspy.PolyPath(moat_points, width = ground_plane_moat_width_fine)    
        D_jjs.add_polygon(moat_path.polygons, layer = layers['m2m'])
        #bottom
        coord_pairs = [[x_min+ground_plane_buffer_fine+ground_plane_moat_width_fine/2,y_min-ground_plane_buffer_fine-ground_plane_moat_width_fine/2],
                       [D_jjs.ports['south_east_east'].midpoint[0]-ground_plane_buffer_fine-ground_plane_moat_width_fine/2,y_min-ground_plane_buffer_fine-ground_plane_moat_width_fine/2]]
        moat_points = []
        for ii in range(len(coord_pairs)):
            moat_points.append(np.array(coord_pairs[ii])) 
        moat_path = gdspy.PolyPath(moat_points, width = ground_plane_moat_width_fine)    
        D_jjs.add_polygon(moat_path.polygons, layer = layers['m2m'])
        #top
        coord_pairs = [[x_min+ground_plane_buffer_fine+ground_plane_moat_width_fine/2,y_max+ground_plane_buffer_fine+ground_plane_moat_width_fine/2],
                       [x_max-ground_plane_buffer_fine-ground_plane_moat_width_fine/2,y_max+ground_plane_buffer_fine+ground_plane_moat_width_fine/2]]
        moat_points = []
        for ii in range(len(coord_pairs)):
            moat_points.append(np.array(coord_pairs[ii])) 
        moat_path = gdspy.PolyPath(moat_points, width = ground_plane_moat_width_fine)    
        D_jjs.add_polygon(moat_path.polygons, layer = layers['m2m'])

    #label
    if jj_include_shunt == True:
        text_string = 'jj_series\nn_jjs = '+str(num_jjs_per_col*num_cols_jjs)+'\nd = '+str(jj_junc_diam)+' um\nr = '+str(jj_shunt_resistance)+' ohm'
    elif jj_include_shunt == False:
        text_string = 'jj_series\nn_jjs = '+str(num_jjs_per_col*num_cols_jjs)+'\nd = '+str(jj_junc_diam)+' um'    
    
    Text_label = vt_label_maker(text_string = text_string, text_size = device_label_size, justify = 'center', layer = layers[label_layer])
    text_coords = [D_jjs.center[0],D_jjs.ymax+Text_label.ysize/2+device_label_size]#[pad1.ports['m3_north'].midpoint[0]+(pad2.ports['m3_north'].midpoint[0]-pad1.ports['m3_north'].midpoint[0])/2,D_jjs.ymax+tn*device_label_size]
    text_label = D_jjs.add_ref(Text_label)
    text_label.center = text_coords
    text_label_size_x = text_label.size[0]
    
    #pads
    params['is_ground_pad'] = False    
    D_pad = jj_pad(params)
    params['is_ground_pad'] = True
    D_pad_gnd = jj_pad(params)

    x1 = D_jjs.ports['south_west_south'].midpoint[0]
    y1 = D_jjs.ports['south_west_south'].midpoint[1]
    pad1 = D_jjs.add_ref(D_pad)
    pad1.center = [x1-pad_size[0]/2,y1-pad_size[1]/2-pad_y_backset]
    tn = 0.65
    if jj_include_shunt == False:
        tn = 0.7
    wire1 = wire_tri_seg(p1 = pad1.ports['m3_north'], p2 = D_jjs.ports['south_west_west'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [tn,0.1], directions = 'yx', layer = layers['m3'])
    
    x1 = D_jjs.ports['south_east_south'].midpoint[0]
    pad2 = D_jjs.add_ref(D_pad_gnd)
    pad2.center = [x1+pad_size_ground[0],y1-pad_y_backset+pad_size_ground[1]/2]
    wire2 = wire_tri_seg(p1 = pad2.ports['m3_north'], p2 = D_jjs.ports['south_east_east'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.35,0.1], directions = 'yx', layer = layers['m3'])
    
    D_jjs.add_ref(wire1)
    D_jjs.add_ref(wire2)
    
    D_jjs.add_port(name = 'pad_anchor', midpoint = [pad1.ports['m3_north'].midpoint[0],pad1.ports['m3_west'].midpoint[1]], width = pad_size[0], orientation = 90)
        
    #pad label
#    Text_label = vt_label_maker(text_string = 'gnd', text_size = pad_label_size, layer = layers[label_layer])
#    text_label = D_jjs.add_ref(Text_label)
#    text_label.rotate(90)
#    text_label.center = pad2.ports['text'].midpoint
                
    #flux moats 
    if include_big_moats == True:
        x1 = pad1.x-pad_w_wire/2-ground_plane_buffer-ground_plane_moat_width/2
        x2 = min(D_jjs.ports['south_west_west'].midpoint[0]-jjs.size[0]-2*ground_plane_buffer-3*ground_plane_moat_width/2,D_jjs.x-max(jjs_size[0]/2,text_label_size_x/2))
        x3 = max(D_jjs.ports['south_east_east'].midpoint[0]+jjs.size[0]+2*ground_plane_buffer+3*ground_plane_moat_width/2,D_jjs.x+max(jjs_size[0]/2,text_label_size_x/2))
        x4 = pad2.x+pad_w_wire/2+ground_plane_buffer+ground_plane_moat_width/2
        y1 = pad1.ymax+ground_plane_buffer
        y2 = D_jjs.ports['south_west_west'].midpoint[1]
        y3 = D_jjs.ymax+ground_plane_moat_width/2+ground_plane_buffer
        y4 = pad2.ymax+ground_plane_buffer
        #outer
        coord_pairs = [[pad1.xmin-ground_plane_buffer-ground_plane_moat_width/2,pad1.ymin],
                       [pad1.xmin-ground_plane_buffer-ground_plane_moat_width/2,y1],
                       [x1,y1],[x1,y2],[x2,y3],[x3,y3],[x4,y2],[x4,y4]]
        moat_points = []
        for ii in range(len(coord_pairs)):
            moat_points.append(np.array(coord_pairs[ii])) 
        moat_path = gdspy.PolyPath(moat_points, width = ground_plane_moat_width)    
        D_jjs.add_polygon(moat_path.polygons, layer = layers['m2m'])
        #inner
        if pad1.ymax+ground_plane_buffer+ground_plane_moat_width/2 > pad2.ymin-ground_plane_buffer-ground_plane_moat_width/2:
            coord_pairs = [[pad1.xmax+ground_plane_buffer+ground_plane_moat_width/2,pad1.ymin],
                           [pad1.xmax+ground_plane_buffer+ground_plane_moat_width/2,pad2.ymin-ground_plane_buffer-ground_plane_moat_width/2],
                           [pad2.xmax,pad2.ymin-ground_plane_buffer-ground_plane_moat_width/2]]
        else:
            coord_pairs = [[pad1.xmax+ground_plane_buffer+ground_plane_moat_width/2,pad1.ymin],
                           [pad1.xmax+ground_plane_buffer+ground_plane_moat_width/2,pad1.ymax+ground_plane_buffer+ground_plane_moat_width/2],
                           [pad1.x+pad_w_wire/2+ground_plane_buffer+ground_plane_moat_width/2,pad1.ymax+ground_plane_buffer+ground_plane_moat_width/2],
                           [pad1.x+pad_w_wire/2+ground_plane_buffer+ground_plane_moat_width/2,pad2.ymin-ground_plane_buffer-ground_plane_moat_width/2],
                           [pad2.xmax,pad2.ymin-ground_plane_buffer-ground_plane_moat_width/2]]
        moat_points = []
        for ii in range(len(coord_pairs)):
            moat_points.append(np.array(coord_pairs[ii])) 
        moat_path = gdspy.PolyPath(moat_points, width = ground_plane_moat_width)    
        D_jjs.add_polygon(moat_path.polygons, layer = layers['m2m'])
        
        coord_pairs = [[pad1.x+pad_w_wire/2+ground_plane_buffer+ground_plane_moat_width/2,pad2.ymax+ground_plane_buffer+ground_plane_moat_width/2],
                       [pad2.x-pad_w_wire/2-ground_plane_buffer-ground_plane_moat_width/2,pad2.ymax+ground_plane_buffer+ground_plane_moat_width/2]]
        moat_points = []
        for ii in range(len(coord_pairs)):
            moat_points.append(np.array(coord_pairs[ii])) 
        moat_path = gdspy.PolyPath(moat_points, width = ground_plane_moat_width)    
        D_jjs.add_polygon(moat_path.polygons, layer = layers['m2m'])
#    #right
#    coord_pairs = [[pad3.xmax+ground_plane_buffer+ground_plane_moat_width/2,pad3.ymin],[pad3.xmax+ground_plane_buffer+ground_plane_moat_width/2,pad4.ymin-ground_plane_buffer]]
#    moat_points = []
#    for ii in range(len(coord_pairs)):
#        moat_points.append(np.array(coord_pairs[ii])) 
#    moat_path = gdspy.PolyPath(moat_points, width = ground_plane_moat_width)    
#    D_jjs.add_polygon(moat_path.polygons, layer = layers['m2m'])
        
    return D_jjs

def jj_2stitch(params = dict()):
        
    inter_jj_gap = vt_arg_helper(params,'inter_jj_gap',1)
    w_wire = vt_arg_helper(params,'w_wire',1)
    layers = vt_arg_helper(params,'layers',vt_lyrs)
    
    D_jjs = Device('jj_2stitch')
                  
    jj = jj_circle(params)         

    pitch = jj.ysize+inter_jj_gap
    jj1 = D_jjs.add_ref(jj)
    jj2 = D_jjs.add_ref(jj)
    jj2.movey(pitch)
    wire = wire_basic(jj1.ports['jj1_north'],jj2.ports['jj1_south'],'y',w_wire,vt_lyrs['jj1'])
    D_jjs.add_ref(wire)
    D_jjs.add_port(name = 'jj1_north', midpoint = jj1.ports['jj1_north'].midpoint, width = w_wire, orientation = 90)
    D_jjs.add_port(name = 'jj1_south', midpoint = jj1.ports['jj1_south'].midpoint, width = w_wire, orientation = 270)
    wire_length = (jj2.ports['m3_south'].midpoint[1]-jj1.ports['m3_north'].midpoint[1])/2
    jj1_wire = pg.rectangle(size = [w_wire,wire_length], layer = layers['m3'])
    jj1_wire.add_port(name = 'wire_north', midpoint = [w_wire/2,wire_length], width = w_wire, orientation = 90)
    jj1_wire.add_port(name = 'wire_south', midpoint = [w_wire/2,0], width = w_wire, orientation = 270)
    wire1 = D_jjs.add_ref(jj1_wire)
    wire1.connect(port = 'wire_north', destination = jj1.ports['m3_south'])
    wire2 = D_jjs.add_ref(jj1_wire)
    wire2.connect(port = 'wire_south', destination = jj2.ports['m3_north'])
    
    D_jjs.add_port(name = 'north', midpoint = wire2.ports['wire_north'].midpoint, width = w_wire, orientation = 90)
    D_jjs.add_port(name = 'south', midpoint = wire1.ports['wire_south'].midpoint, width = w_wire, orientation = 270)
        
    return D_jjs


def jj_with_leads(params = dict()):
    
    jj_xcntct_length = vt_arg_helper(params,'jj_xcntct_length',8.6)
    w_wire = vt_arg_helper(params,'w_wire',2.6)
    layers = vt_arg_helper(params,'layers',vt_lyrs)
    
    D_jj = Device('jj_pads_and_leads')
    
    #jj
    jj = D_jj.add_ref(jj_circle(params))
    
    #xcntct
    Xcntct = Device('jj_xcntct')
    Xcntct.add_ref(pg.union(pg.cross(length = 2*jj_xcntct_length, width = w_wire, layer = layers['m3']), by_layer = False, layer = layers['m3']))
    Xcntct.add_port(name = 'north', midpoint = [0,jj_xcntct_length], width = w_wire, orientation = 90)
    Xcntct.add_port(name = 'south', midpoint = [0,-jj_xcntct_length], width = w_wire, orientation = 270)
    Xcntct.add_port(name = 'east', midpoint = [jj_xcntct_length,0], width = w_wire, orientation = 0)
    Xcntct.add_port(name = 'west', midpoint = [-jj_xcntct_length,0], width = w_wire, orientation = 180)
    xcntct = D_jj.add_ref(Xcntct)
    xcntct.connect(port = 'south', destination = jj.ports['m3_north'])
    
    #50 ohm
    res = D_jj.add_ref(res_50_ohm(params))
    res.connect(port = 'west', destination = xcntct.ports['north'])
    
    #lower cntct
    D_jj.add_ref(wire_basic(jj.ports['jj1_south'],[jj.ports['jj1_south'].midpoint[0],jj.ports['jj1_south'].midpoint[1]-jj_xcntct_length],'y',w_wire,layers['jj1']))
    
    #add ports
    D_jj.add_port(name = 'xcntct_north', midpoint = res.ports['east'].midpoint, width = w_wire, orientation = 90)
    D_jj.add_port(name = 'xcntct_east', midpoint = xcntct.ports['east'].midpoint, width = w_wire, orientation = 0)
    D_jj.add_port(name = 'xcntct_west', midpoint = xcntct.ports['west'].midpoint, width = w_wire, orientation = 180)
    D_jj.add_port(name = 'jj1_south', midpoint = [jj.ports['jj1_south'].midpoint[0],jj.ports['jj1_south'].midpoint[1]-jj_xcntct_length], width = w_wire, orientation = 270)
    
    return D_jj