import copy
import numpy as np
import gdspy
import phidl, phidl.geometry as pg
from phidl import Device, Layer, LayerSet, Port
from phidl import make_device

from nc_library import wire_basic
from vt_util import vt_layers

from nc_library__vt_util import vt_arg_helper, corner_fixer, vt_label_maker
from nc_library__vt_pads_vias_wires import jj_jj1_m3_via, jj_pad, wire_tri_seg, via_general, vt_m1_v3_via_sequence
from nc_library__vt_res import res_stitch, res_stitch_simp, res_50_ohm
from nc_library__vt_jj import jj_circle

vt_lyrs,layer_data = vt_layers()

def sq_washer(params = dict()):
        
    jj_junc_diam = vt_arg_helper(params,'jj_junc_diam',1.15)
    jj_via_inset = vt_arg_helper(params,'jj_via_inset',0.4)
    jj_top_contact_outset = vt_arg_helper(params,'jj_top_contact_outset',1)
    jj_bottom_contact_outset = vt_arg_helper(params,'jj_bottom_contact_outset',1)
    sq_wash_slit_gap = vt_arg_helper(params,'sq_wash_slit_gap',1)
    sq_wash_w_wire_wide = vt_arg_helper(params,'sq_wash_w_wire_wide',30)
    sq_wash_w_in = vt_arg_helper(params,'sq_wash_w_in',13)
    sq_washer_include_inductor_ports = vt_arg_helper(params,'sq_washer_include_inductor_ports',True)    
    sq_ground_plane_hole = vt_arg_helper(params,'sq_ground_plane_hole',False)
    sq_gph_outset = vt_arg_helper(params,'sq_gph_outset',10)
    via_width = vt_arg_helper(params,'via_width',2)
    jj_or_via_termination = vt_arg_helper(params,'jj_or_via_termination','jj')
    w_wire = vt_arg_helper(params,'w_wire',1)
    res_w_wire = vt_arg_helper(params,'res_w_wire',2)
    layers = vt_arg_helper(params,'layers',vt_lyrs)
    
    w1 = sq_wash_slit_gap
    w2 = sq_wash_w_wire_wide
    w3 = sq_wash_w_in
    D_sq = Device('squid washer')
    D_half = Device('half squid washer')
    
    D_half.add_polygon( [(-w1/2, 0), (-w3/2-w2,0), (-w3/2-w2,2*w2+w3), (0,2*w2+w3), (0,w2+w3), (-w3/2,w2+w3), (-w3/2,w2), (-w1/2,w2)],layers['jj1'])
    D_sq.add_ref(D_half)
    right_half = D_sq.add_ref(D_half)
    right_half.reflect(p1 = [0,0], p2 = [0,1])
    
    if jj_or_via_termination == 'jj':
        via_diam = jj_junc_diam-2*jj_via_inset
        sq1_w = via_diam+2*jj_top_contact_outset
        sq2_w = jj_junc_diam+2*jj_bottom_contact_outset        
    elif jj_or_via_termination == 'via':
        via_diam = via_width
        sq1_w = via_diam+2*jj_top_contact_outset
        sq2_w = via_diam+2*jj_bottom_contact_outset
    
    D_sq.add_port(name = 'jj_left', midpoint = [-w1/2-sq2_w/2,0], width = sq2_w, orientation = 270)
    D_sq.add_port(name = 'jj_right', midpoint = [w1/2+sq2_w/2,0], width = sq2_w, orientation = 270)    
    D_sq.add_port(name = 'm3_left', midpoint = [-w1/2-sq1_w/2,0], width = sq2_w, orientation = 270)
    D_sq.add_port(name = 'm3_right', midpoint = [w1/2+sq1_w/2,0], width = sq2_w, orientation = 270)    
    
    y_coord = 0
    D_sq.add_port(name = 'bias_north', midpoint = [D_sq.xsize/2,y_coord+w_wire/2], width = w_wire, orientation = 0)#this is the mobile bias point. it is no longer north. sorry. moved it due to conversation with Ben and Malcolm
    D_sq.add_port(name = 'washer_center', midpoint = [0,D_sq.ysize/2], width = 0, orientation = 180)
    D_sq.add_port(name = 'shunt', midpoint = [0,D_sq.ymin+max(via_width/2+jj_top_contact_outset,res_w_wire/2)], width = 0, orientation = 180)
    x_size = D_sq.xsize
            
    #ground plane hole
    if sq_ground_plane_hole == True:
        Ground_plane_hole = pg.rectangle(size = D_sq.size+2*sq_gph_outset, layer = layers['m2'])
        gph = D_sq.add_ref(Ground_plane_hole)
        gph.center = D_sq.ports['washer_center'].midpoint 
        
    #inductor ports
    if sq_washer_include_inductor_ports == True:
        w4 = x_size/2-w1/2-2
        Inductor_port = pg.rectangle(size = [w4,2], layer = layers['ip1'])
        in_po1 = D_sq.add_ref(Inductor_port)
        in_po2 = D_sq.add_ref(Inductor_port)
        in_po1.center = [-w1/2-1-w4/2,0]
        in_po2.center = [w1/2+1+w4/2,0]
        D_sq.label(text = 'P1+ JJ1', position = [-w1/2-1-w4/2,0.1], layer = layers['ipt'])
        D_sq.label(text = 'P1- JJ1', position = [w1/2+1+w4/2,0.1], layer = layers['ipt'])
           
    return D_sq

def sq_input_coil_one_layer(params = dict()):
    
    sq_incoil_w_wire = vt_arg_helper(params,'sq_incoil_w_wire',1)
    sq_incoil_w_coil = vt_arg_helper(params,'sq_incoil_w_coil',10)
    sq_incoil_l_leads = vt_arg_helper(params,'sq_incoil_l_leads',5)
    sq_incoil_lead_gap = vt_arg_helper(params,'sq_incoil_lead_gap',2)
    sq_incoil_include_inductor_ports = vt_arg_helper(params,'sq_incoil_include_inductor_ports',True)
    layers = vt_arg_helper(params,'layers',vt_lyrs) 
    
    D_sq = Device('squid input coil _ one layer')
    
    w1 = sq_incoil_w_wire
    w2 = sq_incoil_w_coil
    l1 = sq_incoil_l_leads
    l2 = l1+w2
    g = sq_incoil_lead_gap
    
    D_half = Device('half squid input coil')
    D_half.add_polygon( [(0,g/2),(0,g/2+w1),(l1,g/2+w1),(l1,w2/2),(l2,w2/2),(l2,0),(l2-w1,0),(l2-w1,w2/2-w1),(l1+w1,w2/2-w1),(l1+w1,g/2)],layers['m3'])
    D_sq.add_ref(D_half)
    lower_half = D_sq.add_ref(D_half)
    lower_half.reflect(p1 = [0,0], p2 = [1,0])
    
    D_sq.add_port(name = 'upper', midpoint = [0,g/2+w1/2], width = w1, orientation = 180)
    D_sq.add_port(name = 'lower', midpoint = [0,-g/2-w1/2], width = w1, orientation = 180)
    D_sq.add_port(name = 'coil_center', midpoint = [l1+w2/2,0], width = w1, orientation = 0)
    
    #inductor ports
    if sq_incoil_include_inductor_ports == True:
        Inductor_port = pg.rectangle(size = [2,w1], layer = layers['ipm3'])
        in_po1 = D_sq.add_ref(Inductor_port)
        in_po2 = D_sq.add_ref(Inductor_port)
        in_po1.center = D_sq.ports['upper'].midpoint
        in_po2.center = D_sq.ports['lower'].midpoint
        D_sq.label(text = 'P2+ M3', position = [D_sq.ports['upper'].midpoint[0]+0.1,D_sq.ports['upper'].midpoint[1]], layer = layers['ipl'])
        D_sq.label(text = 'P2- M3', position = [D_sq.ports['lower'].midpoint[0]+0.1,D_sq.ports['lower'].midpoint[1]], layer = layers['ipl'])
    
    return D_sq

def sq_two_layer(params = dict()):
    
    jj_junc_diam = vt_arg_helper(params,'jj_junc_diam',2)
    jj_top_contact_outset = vt_arg_helper(params,'jj_top_contact_outset',1)
    jj_bottom_contact_outset = vt_arg_helper(params,'jj_bottom_contact_outset',1.1)
    sq_wash_w_in = vt_arg_helper(params,'sq_wash_w_in',16.5)
    sq_incoil_outset = vt_arg_helper(params,'sq_incoil_outset',5)
    sq_incoil_w_wire = vt_arg_helper(params,'sq_incoil_w_wire',1)
    sq_l_leads = vt_arg_helper(params,'sq_l_leads',5)
    sq_include_via = vt_arg_helper(params,'sq_include_via',False)  
    sq_include_shunt = vt_arg_helper(params,'sq_include_shunt',True)
    sq_squid_include_inductor_ports = vt_arg_helper(params,'sq_squid_include_inductor_ports',True)
    sq_res_num_squares = vt_arg_helper(params,'sq_res_num_squares',5)
    sq_res_w_wire = vt_arg_helper(params,'sq_res_w_wire',2)
    sq_res_w_cntct = vt_arg_helper(params,'sq_res_w_cntct',1)
    sq_res_outset = vt_arg_helper(params,'sq_res_outset',1)
    sq_res_l_lead = vt_arg_helper(params,'sq_res_l_lead',2)
    sq_res_w_lead = vt_arg_helper(params,'sq_res_w_lead',3)
    sq_washer_jj_offset = vt_arg_helper(params,'sq_washer_jj_offset',5)
    via_width = vt_arg_helper(params,'via_width',4)
    w_wire = vt_arg_helper(params,'w_wire',1)
    layers = vt_arg_helper(params,'layers',vt_lyrs)
                 
    params['res_num_squares'] = sq_res_num_squares                   
    params['res_w_wire'] = sq_res_w_wire                  
    params['res_w_cntct'] = sq_res_w_cntct
    params['res_outset'] = sq_res_outset
    params['res_l_lead'] = sq_res_l_lead
    params['res_w_lead'] = sq_res_w_lead
    
    D_sq = Device('squid')
   
    #washer
    D_washer = sq_washer(params)    
    washer = D_sq.add_ref(D_washer)
    
    #jjs
    w_temp = jj_junc_diam+2*jj_bottom_contact_outset
    if sq_squid_include_inductor_ports == False:
        
        D_jj = jj_circle(params)         
    
        jj1 = D_sq.add_ref(D_jj)
        jj2 = D_sq.add_ref(D_jj)
        jj1.reflect([0,0],[0,1])
        jj1.connect(port = 'jj1_north', destination = washer.ports['jj_left'])    
        jj2.connect(port = 'jj1_north', destination = washer.ports['jj_right'])
        #move back one square for litho
        jj1.movey(-sq_washer_jj_offset)
        jj_wire1 = wire_basic(jj1.ports['jj1_north'],washer.ports['jj_left'],'y',w_temp,layers['jj1'])
        jj2.movey(-sq_washer_jj_offset)
        jj_wire2 = wire_basic(jj2.ports['jj1_north'],washer.ports['jj_right'],'y',w_temp,layers['jj1'])
        D_sq.add_ref(jj_wire1)
        D_sq.add_ref(jj_wire2)
    
    #inductor ports
    if sq_squid_include_inductor_ports == True:
        Inductor_port = pg.rectangle(size = [w_temp-1,2], layer = layers['ipj1'])
        in_po1 = D_sq.add_ref(Inductor_port)
        in_po2 = D_sq.add_ref(Inductor_port)
        in_po1.center = washer.ports['jj_left'].midpoint
        in_po2.center = washer.ports['jj_right'].midpoint
        D_sq.label(text = 'P1+ JJ1', position = washer.ports['jj_left'].midpoint, layer = layers['ipl'])
        D_sq.label(text = 'P1- JJ1', position = washer.ports['jj_right'].midpoint, layer = layers['ipl'])    
    
    if sq_squid_include_inductor_ports == False:
        
        jj_wire = wire_basic(jj1.ports['m3_west'],jj2.ports['m3_west'],'x',w_wire,vt_lyrs['m3'])
        mp = [jj1.ports['m3_west'].midpoint[0]+(jj2.ports['m3_west'].midpoint[0]-jj1.ports['m3_west'].midpoint[0])/2,jj1.ports['m3_south'].midpoint[1]+(jj1.ports['m3_north'].midpoint[1]-jj1.ports['m3_south'].midpoint[1])/2-w_wire/2]    
        D_sq.add_ref(jj_wire)
        D_sq.add_port(name = 'bias_point', midpoint = mp, width = w_wire, orientation = 270)
        
        #via to pad leads  
        if  sq_include_via == True:
            Via_sq = jj_jj1_m3_via(via_width,jj_bottom_contact_outset,jj_top_contact_outset,w_wire,layers)
            via_sq = D_sq.add_ref(Via_sq)
            via_sq.connect(port = 'jj1_south', destination = washer.ports['bias_north'])
        
        #wire for two connections
        l_wire = sq_l_leads
        if sq_include_via == False:
            contact_wire = pg.rectangle(size = [w_wire,l_wire], layer = layers['jj1'])
            contact_wire.add_port(name = 'north', midpoint = [w_wire/2,l_wire], width = w_wire, orientation = 90)
            contact_wire.add_port(name = 'south', midpoint = [w_wire/2,0], width = w_wire, orientation = 270)
            contact_wire.add_port(name = 'north_east', midpoint = [w_wire,l_wire-w_wire/2], width = w_wire, orientation = 0)
            contact_wire.add_port(name = 'north_west', midpoint = [0,l_wire-w_wire/2], width = w_wire, orientation = 180)
            contact_wire.add_port(name = 'south_east', midpoint = [w_wire,l_wire/2], width = w_wire, orientation = 0)
            contact_wire.add_port(name = 'south_west', midpoint = [0,l_wire/2], width = w_wire, orientation = 180)
        elif sq_include_via == True:
            contact_wire = pg.rectangle(size = [w_wire,l_wire], layer = layers['m3'])
            contact_wire.add_port(name = 'north', midpoint = [w_wire/2,l_wire], width = w_wire, orientation = 90)
            contact_wire.add_port(name = 'south', midpoint = [w_wire/2,0], width = w_wire, orientation = 270)
            contact_wire.add_port(name = 'north_east', midpoint = [w_wire,l_wire-w_wire/2], width = w_wire, orientation = 0)
            contact_wire.add_port(name = 'north_west', midpoint = [0,l_wire-w_wire/2], width = w_wire, orientation = 180)
            contact_wire.add_port(name = 'south_east', midpoint = [w_wire,l_wire/2], width = w_wire, orientation = 0)
            contact_wire.add_port(name = 'south_west', midpoint = [0,l_wire/2], width = w_wire, orientation = 180)
            
        cw1 = D_sq.add_ref(contact_wire)
        if sq_include_via == False:
            cw3 = D_sq.add_ref(contact_wire)
            cw1.connect(port = 'south', destination = washer.ports['bias_north'])
            cw3.connect(port = 'south', destination = cw1.ports['north_east'])
        elif sq_include_via == True:
            cw3 = D_sq.add_ref(contact_wire)
            cw1.connect(port = 'south', destination = via_sq.ports['m3_north']) 
            cw3.connect(port = 'south', destination = cw1.ports['north_east'])
        
        contact_wire = pg.rectangle(size = [w_wire,l_wire], layer = layers['m3'])
        contact_wire.add_port(name = 'north', midpoint = [w_wire/2,l_wire], width = w_wire, orientation = 90)
        contact_wire.add_port(name = 'south', midpoint = [w_wire/2,0], width = w_wire, orientation = 270)
        contact_wire.add_port(name = 'north_east', midpoint = [w_wire,l_wire-w_wire/2], width = w_wire, orientation = 0)
        contact_wire.add_port(name = 'north_west', midpoint = [0,l_wire-w_wire/2], width = w_wire, orientation = 180)
        contact_wire.add_port(name = 'south_east', midpoint = [w_wire,l_wire/2], width = w_wire, orientation = 0)
        contact_wire.add_port(name = 'south_west', midpoint = [0,l_wire/2], width = w_wire, orientation = 180)
        cw2 = D_sq.add_ref(contact_wire)
        cw2.connect(port = 'north', destination = D_sq.ports['bias_point'])
        
        D_sq.add_port(name = 'north', midpoint = washer.ports['bias_north'].midpoint, width = w_wire, orientation = 0)
        D_sq.add_port(name = 'south', midpoint = cw2.ports['south'].midpoint, width = w_wire, orientation = 270)
        if sq_include_via == False:
            D_sq.add_port(name = 'north_east', midpoint = cw3.ports['north_east'].midpoint, width = w_wire, orientation = 0)
            D_sq.add_port(name = 'north_west', midpoint = cw3.ports['north_west'].midpoint, width = w_wire, orientation = 180)
            D_sq.add_port(name = 'north_south', midpoint = cw3.ports['north'].midpoint, width = w_wire, orientation = 270)
        if sq_include_via == True:
            D_sq.add_port(name = 'north_east', midpoint = cw3.ports['north_east'].midpoint, width = w_wire, orientation = 0)
            D_sq.add_port(name = 'north_west', midpoint = cw3.ports['north_west'].midpoint, width = w_wire, orientation = 180)
        D_sq.add_port(name = 'south_east', midpoint = cw2.ports['south_east'].midpoint, width = w_wire, orientation = 0)
        D_sq.add_port(name = 'south_west', midpoint = cw2.ports['south_west'].midpoint, width = w_wire, orientation = 180)
                
    #shunt resistor
    if sq_include_shunt == True:
        params_mod = copy.deepcopy(params)
        params_mod['res_include_via'] = True
        params_mod['res_lower_layer'] = 'm3'
        Res = res_stitch(params)
        res = D_sq.add_ref(Res)
        res.connect(port = 'stitch_center', destination = washer.ports['shunt'])
    
    #input coil
    params['sq_incoil_w_coil'] = sq_wash_w_in+2*(sq_incoil_outset+sq_incoil_w_wire)
    D_incoil = sq_input_coil_one_layer(params)
    d_incoil = D_sq.add_ref(D_incoil)
    d_incoil.connect(port = 'coil_center', destination = washer.ports['washer_center']) 
    
    D_sq.add_port(name = 'center', midpoint = washer.ports['washer_center'].midpoint, width = w_wire, orientation = 180)
    D_sq.add_port(name = 'incoil_north', midpoint = d_incoil.ports['upper'].midpoint, width = sq_incoil_w_wire, orientation = 180)
    D_sq.add_port(name = 'incoil_south', midpoint = d_incoil.ports['lower'].midpoint, width = sq_incoil_w_wire, orientation = 180)
            
    return D_sq

def sq_8wire(params = dict()):

    jj_junc_diam = vt_arg_helper(params,'jj_junc_diam',2)
    sq_wash_w_in = vt_arg_helper(params,'sq_wash_w_in',16.5)
    sq_l_leads = vt_arg_helper(params,'sq_l_leads',20)
    sq_incoil_w_wire = vt_arg_helper(params,'sq_incoil_w_wire',1)
    sq_incoil_outset = vt_arg_helper(params,'sq_incoil_outset',5)
    sq_addflux_include = vt_arg_helper(params,'sq_addflux_include',True)
    sq_addflux_w_wire = vt_arg_helper(params,'sq_addflux_w_wire',2)
    sq_addflux_outset = vt_arg_helper(params,'sq_addflux_outset',2)
    sq_addflux_extent = vt_arg_helper(params,'sq_addflux_extent',20)    
    sq_include_via = vt_arg_helper(params,'sq_include_via',False)
    sq_include_shunt = vt_arg_helper(params,'sq_include_shunt',True)
    shunt_num_squares = vt_arg_helper(params,'shunt_num_squares',1.15)#for jj shunts
    sq_res_num_squares = vt_arg_helper(params,'sq_res_num_squares',0.13)#for resonance-dampling shunt
    w_wire = vt_arg_helper(params,'w_wire',1)
    include_50_ohm = vt_arg_helper(params,'include_50_ohm',True)
    res_50_ohm_w_wire = vt_arg_helper(params,'res_50_ohm_w_wire',1)
    res_50_ohm_num_squares = vt_arg_helper(params,'res_50_ohm_num_squares',25)
    res_50_ohm_outset = vt_arg_helper(params,'res_50_ohm_outset',0.2)
    sq_res_50_ohm_lateral_gap = vt_arg_helper(params,'sq_res_50_ohm_lateral_gap',10)
    pad_size = vt_arg_helper(params,'pad_size',[200,250])
    pad_size_ground = vt_arg_helper(params,'pad_size_ground',[100,150])
    pad_pitch = vt_arg_helper(params,'pad_pitch',[200,375])
    pad_y_backset = vt_arg_helper(params,'pad_y_backset',100)
    pad_w_wire = vt_arg_helper(params,'pad_w_wire',20)
    pad_label_size = vt_arg_helper(params,'pad_label_size',10)
    device_label_size = vt_arg_helper(params,'device_label_size',5)    
    include_big_moats = vt_arg_helper(params,'include_big_moats',False)
    ground_plane_moat_width = vt_arg_helper(params,'ground_plane_moat_width',5)
    ground_plane_buffer = vt_arg_helper(params,'ground_plane_buffer',15)
    ground_plane_moat_width_fine = vt_arg_helper(params,'ground_plane_moat_width_fine',1)
    ground_plane_buffer_fine = vt_arg_helper(params,'ground_plane_buffer_fine',10)
    label_layer = vt_arg_helper(params,'label_layer','m4l')
    layers = vt_arg_helper(params,'layers',vt_lyrs)
    
    D_sq = Device('squid')
                      
    Sq = sq_two_layer(params)
    sq = D_sq.add_ref(Sq)
    
    #corner fixer 
    Cf = corner_fixer(w_wire,layers['m3'])
    cf2 = D_sq.add_ref(Cf)
    cf2.connect(port = 'north', destination = sq.ports['south'])

    sq_size = D_sq.size
    
    #50 ohm in series
    if include_50_ohm == True:
        l_wire = (sq_res_50_ohm_lateral_gap-cf2.xsize)/2
        w_wire_res_cntct = res_50_ohm_w_wire-2*res_50_ohm_outset
        Wire_lateral = pg.rectangle(size = [l_wire,w_wire], layer = layers['m3'])
        Wire_lateral.add_port(name = 'west', midpoint = [0,w_wire/2], width = w_wire, orientation = 180)
        Wire_lateral.add_port(name = 'west_south', midpoint = [w_wire_res_cntct/2,0], width = w_wire_res_cntct, orientation = 270)
        Wire_lateral.add_port(name = 'east', midpoint = [l_wire,w_wire/2], width = w_wire, orientation = 0)
        Wire_lateral.add_port(name = 'east_south', midpoint = [l_wire-w_wire_res_cntct/2,0], width = w_wire_res_cntct, orientation = 270)
        wire_lateral_1 = D_sq.add_ref(Wire_lateral)
        wire_lateral_1.connect(port = 'east',destination = cf2.ports['west'])
        wire_lateral_2 = D_sq.add_ref(Wire_lateral)
        wire_lateral_2.connect(port = 'west',destination = cf2.ports['east'])
        params_mod = copy.deepcopy(params)
        params_mod['res_w_wire'] = res_50_ohm_w_wire
        params_mod['res_num_squares'] = res_50_ohm_num_squares
        params_mod['res_outset'] = res_50_ohm_outset
        params_mod['res_lower_layer'] = 'm3'
        Res = res_stitch(params_mod)
        res1 = D_sq.add_ref(Res)
        res1.connect(port = 'm3_east',destination = wire_lateral_1.ports['west_south'])
        Res = res_stitch(params_mod)
        res2 = D_sq.add_ref(Res)
        res2.connect(port = 'm3_west',destination = wire_lateral_2.ports['east_south'])
        #moats
        Wire_moat = pg.rectangle(size = [ground_plane_moat_width_fine,res_50_ohm_num_squares*res_50_ohm_w_wire], layer = layers['m2m'])
        for ii in range(3):
            wire_moat_1 = D_sq.add_ref(Wire_moat)
            wire_moat_1.center = res1.center
            wire_moat_1.movex(-res1.xsize/2-(ii+1)*(ground_plane_buffer_fine+ground_plane_moat_width_fine)+ground_plane_moat_width_fine/2)
            wire_moat_2 = D_sq.add_ref(Wire_moat)
            wire_moat_2.center = res2.center
            wire_moat_2.movex(res1.xsize/2+(ii+1)*(ground_plane_buffer_fine+ground_plane_moat_width_fine)-ground_plane_moat_width_fine/2)

    #pads for squid
    x_mid = sq.ports['south'].midpoint[0]
    pad_space_x = pad_pitch[0]-pad_size[0]    
    params['is_ground_pad'] = False
    D_pad = jj_pad(params)
    params['is_ground_pad'] = True
    D_pad_gnd = jj_pad(params)
   
    pad3 = D_sq.add_ref(D_pad)
    sq_extra_offset = Res.xsize/2+sq_l_leads+50
    pad3_x = x_mid-pad_size[0]/2-pad_space_x/2
    pad5_x = pad3_x+pad_pitch[0]
    pad7_x = pad3_x+2*pad_pitch[0]
    pad1_x = pad3_x-pad_pitch[0]
    pad_y = sq.ymin-pad_size[1]/2-pad_y_backset-sq_extra_offset
    pad3.center = [pad3_x,pad_y]
    pad_y_gap = 60
    if include_50_ohm == True:
        wire3 = wire_tri_seg(p1 = pad3.ports['m3_north'], p2 = res1.ports['m3_west_north'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.7,0.05], directions = 'yx', layer = layers['m3'], layers = layers)
    elif include_50_ohm == False:
        wire3 = wire_tri_seg(p1 = pad3.ports['m3_north'], p2 = cf2.ports['west'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.8,0.1], directions = 'yx', layer = layers['m3'], layers = layers)
         
    pad4 = D_sq.add_ref(D_pad_gnd)
    pad4.center = [sq.ports['north_west'].midpoint[0],pad_y+pad_size[1]/2+pad_size_ground[1]/2+pad_y_gap]#[pad3_x+pad_size[0]/2+pad_size_ground[0]/2,pad_y+pad_size[1]/2+pad_size_ground[1]/2+pad_y_gap]
    wire4 = wire_tri_seg(p1 = pad4.ports['jj1_north'], p2 = sq.ports['north_south'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.4,0.2,0.1], directions = 'yy', layer = layers['jj1'], layers = layers)

    pad5 = D_sq.add_ref(D_pad)
    pad5.center = [pad5_x,pad_y]
        
    pad6 = D_sq.add_ref(D_pad_gnd)
    pad6.center = [pad5_x+pad_size[0]/2+pad_size_ground[0]/2,pad_y+pad_size[1]/2+pad_size_ground[1]/2+pad_y_gap]
    
    if sq_include_via == False:
        if include_50_ohm == True:
#            wire5 = wire_tri_seg(p1 = pad5.ports['jj1_north'], p2 = res2.ports['jj1_east'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.85,0.15], directions = 'yx', layer = layers['jj1'], layers = layers)
            wire5 = wire_tri_seg(p1 = pad5.ports['m3_north'], p2 = res2.ports['m3_east_north'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.72,0.10], directions = 'yx', layer = layers['m3'], layers = layers)
        elif include_50_ohm == False:
#            wire5 = wire_tri_seg(p1 = pad5.ports['jj1_north'], p2 = sq.ports['north_west'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.85,0.15], directions = 'yx', layer = layers['jj1'], layers = layers)
            wire5 = wire_tri_seg(p1 = pad5.ports['m3_north'], p2 = sq.ports['south_east'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.85,0.15], directions = 'yx', layer = layers['jj1'], layers = layers)
        wire6 = wire_tri_seg(p1 = pad6.ports['jj1_north'], p2 = sq.ports['north'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.3,0.25], directions = 'yx', layer = layers['jj1'], layers = layers)  
    elif sq_include_via == True:
        wire5 = wire_tri_seg(p1 = pad5.ports['m3_north'], p2 = sq.ports['north_west'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.85,0.15], directions = 'yx', layer = layers['m3'], layers = layers)
        wire6 = wire_tri_seg(p1 = pad6.ports['m3_north'], p2 = sq.ports['north'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.4,0.15], directions = 'yx', layer = layers['m3'], layers = layers)  
    
    D_sq.add_ref(wire3)
    D_sq.add_ref(wire4)
    D_sq.add_ref(wire5)
    D_sq.add_ref(wire6)
    
    #pads for input coil
    pad1 = D_sq.add_ref(D_pad)
    pad1.center = [pad1_x,pad_y]
    wire1 = wire_tri_seg(p1 = pad1.ports['m3_north'], p2 = sq.ports['incoil_north'], initial_width = pad_w_wire, final_width = sq_incoil_w_wire, length_factors = [0.9,0.1], directions = 'yx', layer = layers['m3'], layers = layers)
    
    pad2 = D_sq.add_ref(D_pad_gnd)
    pad2.center = [pad1_x+pad_size[0]/2+pad_size_ground[0]/2,pad_y+pad_size[1]/2+pad_size_ground[1]/2+pad_y_gap]
    wire2 = wire_tri_seg(p1 = pad2.ports['m3_north'], p2 = sq.ports['incoil_south'], initial_width = pad_w_wire, final_width = sq_incoil_w_wire, length_factors = [0.4,0.05], directions = 'yx', layer = layers['m3'], layers = layers)
    
    D_sq.add_ref(wire1)
    D_sq.add_ref(wire2) 
    
    D_sq.add_port(name = 'pad_anchor', midpoint = [pad1.ports['m3_north'].midpoint[0],pad1.ports['m3_west'].midpoint[1]], width = pad_size[0], orientation = 90)
    
    #linearizing input line for flux-locked loop operation
    if sq_addflux_include == True:
        D_addflux = Device('squid_flux_adder')
        x_coord = sq.ports['center'].midpoint[0]+sq_wash_w_in/2+sq_incoil_outset+sq_incoil_w_wire+sq_addflux_outset+sq_addflux_w_wire/2
        y_coord = sq.ports['center'].midpoint[1]
        wire_vert = wire_basic([x_coord,y_coord-sq_addflux_extent/2],[x_coord,y_coord+sq_addflux_extent/2],'y',sq_addflux_w_wire,layers['m3'])
        D_addflux.add_ref(wire_vert)
        D_addflux.add_port(name = 'addflux_north', midpoint = [x_coord+sq_addflux_w_wire/2,y_coord+sq_addflux_extent/2-sq_addflux_w_wire/2], width = sq_addflux_w_wire, orientation = 0)
        D_addflux.add_port(name = 'addflux_south', midpoint = [x_coord+sq_addflux_w_wire/2,y_coord-sq_addflux_extent/2+sq_addflux_w_wire/2], width = sq_addflux_w_wire, orientation = 0)
        pad7 = D_sq.add_ref(D_pad)
        pad7.center= [pad7_x,pad_y]
        pad8 = D_sq.add_ref(D_pad_gnd)
        pad8.center = [pad7_x+pad_size[0]/2+pad_size_ground[0]/2,pad_y+pad_size[1]/2+pad_size_ground[1]/2+pad_y_gap]
        wire7 = wire_tri_seg(p1 = pad7.ports['m3_north'], p2 = D_addflux.ports['addflux_south'], initial_width = pad_w_wire, final_width = sq_addflux_w_wire, length_factors = [0.77,0.05], directions = 'yx', layer = layers['m3'], layers = layers)
        wire8 = wire_tri_seg(p1 = pad8.ports['m3_north'], p2 = D_addflux.ports['addflux_north'], initial_width = pad_w_wire, final_width = sq_addflux_w_wire, length_factors = [0.75,0.1], directions = 'yx', layer = layers['m3'], layers = layers)
        D_addflux.add_ref(wire7)
        D_addflux.add_ref(wire8)
    
    D_sq.add_ref(D_addflux)   
    
    #labels        
    if sq_include_shunt == True:    
        text_string = 'sq_8wire\nd_jjs = '+str(jj_junc_diam)+' um\nr_jj = '+str(shunt_num_squares)+' sq\nr_sq = '+str(sq_res_num_squares)+' sq\nw_in = '+str(sq_wash_w_in)+' um'
    elif sq_include_shunt == False:
        text_string = 'sq_8wire\nd_jjs = '+str(jj_junc_diam)+' um\nr_jj = '+str(shunt_num_squares)+' sq\nw_in = '+str(sq_wash_w_in)+' um'
    text_coords = [sq.ports['south'].midpoint[0],D_sq.ymax+4.5*device_label_size]    
    Text_label = vt_label_maker(text_string = text_string, text_size = device_label_size, justify = 'center', layer = layers[label_layer])
    text_label = D_sq.add_ref(Text_label)
    text_label.center = text_coords
    
    #pad labels
    text_strings = ['I+','I- / gnd','V+','V- / gnd','incoil+','incoil- / gnd','add_flux I+','add_flux I-']
    text_coords = [pad3.ports['text'].midpoint,pad4.ports['text'].midpoint,pad5.ports['text'].midpoint,pad6.ports['text'].midpoint,pad1.ports['text'].midpoint,pad2.ports['text'].midpoint,pad7.ports['text'].midpoint,pad8.ports['text'].midpoint]
    for ii in range(len(text_strings)):
        Text_label = vt_label_maker(text_string = text_strings[ii], text_size = pad_label_size, layer = layers[label_layer])
        text_label = D_sq.add_ref(Text_label)
        text_label.rotate(90)
        text_label.center = text_coords[ii]
        
    #flux moats
    if include_big_moats == True:
        #outer
        coord_pairs = [[pad1.x-pad_w_wire/2-ground_plane_buffer-ground_plane_moat_width/2,pad1.ymax+3*ground_plane_buffer+ground_plane_moat_width/2],
                       [pad1.x-pad_w_wire/2-ground_plane_buffer-ground_plane_moat_width/2,sq.ports['incoil_north'].midpoint[1]],
                       [sq.ports['south'].midpoint[0]-sq_size[1]/2,D_sq.ymax+ground_plane_buffer+ground_plane_moat_width/2],
                       [sq.ports['south'].midpoint[0]+sq_size[1]/2,D_sq.ymax+ground_plane_buffer+ground_plane_moat_width/2],
                       [pad8.x+pad_w_wire/2+ground_plane_buffer+ground_plane_moat_width/2,sq.ports['incoil_north'].midpoint[1]+ground_plane_buffer+ground_plane_moat_width/2],
                       [pad8.x+pad_w_wire/2+ground_plane_buffer+ground_plane_moat_width/2,pad8.ymax+ground_plane_buffer+ground_plane_moat_width/2],
                       ]
        moat_points = []
        for ii in range(len(coord_pairs)):
            moat_points.append(np.array(coord_pairs[ii])) 
        moat_path = gdspy.PolyPath(moat_points, width = ground_plane_moat_width)    
        D_sq.add_polygon(moat_path.polygons, layer = layers['m2m'])
        #inner
        for kk in range(4):
            coord_pairs = [[pad1.x+pad_w_wire/2+ground_plane_buffer+ground_plane_moat_width/2+kk*2*pad_pitch[0],pad1.ymax+ground_plane_buffer+ground_plane_moat_width/2],
                           [pad3.x-pad_w_wire/2-ground_plane_buffer-ground_plane_moat_width/2+kk*2*pad_pitch[0],pad1.ymax+ground_plane_buffer+ground_plane_moat_width/2]]
            moat_points = []
            for ii in range(len(coord_pairs)):
                moat_points.append(np.array(coord_pairs[ii])) 
            moat_path = gdspy.PolyPath(moat_points, width = ground_plane_moat_width)    
            D_sq.add_polygon(moat_path.polygons, layer = layers['m2m'])
    
    return D_sq

def sq_washer_quad(params = dict()):
    
    sq_jj_junc_diam = vt_arg_helper(params,'sq_jj_junc_diam',1.15)
    jj_bottom_contact_outset = vt_arg_helper(params,'jj_bottom_contact_outset',1)
    sq_wash_slit_gap = vt_arg_helper(params,'sq_wash_slit_gap',1)
    sq_wash_w_wide = vt_arg_helper(params,'sq_wash_w_wide',30)
    sq_wash_jj_offset = vt_arg_helper(params,'sq_wash_jj_offset',6)
    sq_wash_w_bridge = vt_arg_helper(params,'sq_wash_w_bridge',5)
    sq_wash_w_in1 = vt_arg_helper(params,'sq_wash_w_in1',13)
    sq_wash_w_in2 = vt_arg_helper(params,'sq_wash_w_in2',17)
    sq_wash_l_bridge = vt_arg_helper(params,'sq_wash_l_bridge',17)
    sq_incoil_w_wire = vt_arg_helper(params,'sq_incoil_w_wire',2)
    sq_incoil_via_backset = vt_arg_helper(params,'sq_incoil_via_backset',1)
    sq_bias_y_offset = vt_arg_helper(params,'sq_bias_y_offset',7.2)
    sq_shunt_w_lead = vt_arg_helper(params,'sq_shunt_w_lead',5.2)
    w_wire = vt_arg_helper(params,'w_wire',2.2)
    layers = vt_arg_helper(params,'layers',vt_lyrs)
    
    D_sq = Device('squid washer')
    D_half = Device('half squid washer')
    
    jj_junc_diam = sq_jj_junc_diam
    w1 = sq_wash_slit_gap+sq_wash_w_in1+sq_wash_w_in2
    jj_size = jj_junc_diam+2*jj_bottom_contact_outset
    
    x1 = 0
    x2 = x1
    x3 = x2+sq_wash_w_in2
    x4 = x3
    x5 = x4-w1
    x6 = x5
    x7 = x6+sq_wash_w_in1
    x8 = x7
    x9 = x8
    x10 = x9-jj_size
    x11 = x10
    x12 = x8-sq_wash_w_in1-sq_wash_w_wide
    x13 = x12
    x14 = x13+2*sq_wash_w_wide+w1
    x15 = x14
    x16 = x1+sq_wash_w_bridge
    x17 = x16
    
    y1 = 0
    y2 = y1+sq_wash_l_bridge+jj_size+sq_wash_w_wide
    y3 = y2
    y4 = y3+w1
    y5 = y4
    y6 = y3
    y7 = y6
    y8 = y1+sq_wash_l_bridge+jj_size
    y9 = y8-sq_wash_jj_offset
    y10 = y9
    y11 = y8
    y12 = y8
    y13 = y12+2*sq_wash_w_wide+w1
    y14 = y13
    y15 = y8
    y16 = y15
    y17 = y1
    
    D_half.add_polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4),(x5,y5),(x6,y6),(x7,y7),(x8,y8),(x9,y9),(x10,y10),(x11,y11),(x12,y12),(x13,y13),(x14,y14),(x15,y15),(x16,y16),(x17,y17)],layers['jj1'])
    D_half.add_port(name = 'washer_half_port', midpoint = [x1+(x17-x1)/2,y1], width = 1, orientation = 270)
    D_half.add_port(name = 'washer_jj_port', midpoint = [x8-jj_size/2,y8-sq_wash_jj_offset], width = 1, orientation = 270)
    D_half.add_port(name = 'incoil_bridge', midpoint = [D_half.ports['washer_half_port'].midpoint[0],D_half.ports['washer_half_port'].midpoint[1]+sq_wash_l_bridge+jj_size+sq_wash_w_wide+sq_incoil_via_backset], width = sq_incoil_w_wire, orientation = 270)
    one_half = D_sq.add_ref(D_half)
    two_half = D_sq.add_ref(D_half)
    two_half.connect(port = 'washer_half_port', destination = one_half.ports['washer_half_port'])    

    D_sq.add_port(name = 'jj_port_1', midpoint = [one_half.ports['washer_jj_port'].midpoint[0],one_half.ports['washer_jj_port'].midpoint[1]], width = jj_size, orientation = 270)
    D_sq.add_port(name = 'jj_port_2', midpoint = [two_half.ports['washer_jj_port'].midpoint[0],two_half.ports['washer_jj_port'].midpoint[1]], width = jj_size, orientation = 90)
    D_sq.add_port(name = 'incoil_bridge_1', midpoint = [one_half.ports['incoil_bridge'].midpoint[0],one_half.ports['incoil_bridge'].midpoint[1]], width = jj_size, orientation = 90)    
    D_sq.add_port(name = 'incoil_bridge_2', midpoint = [two_half.ports['incoil_bridge'].midpoint[0],two_half.ports['incoil_bridge'].midpoint[1]], width = jj_size, orientation = 270) 
    D_sq.add_port(name = 'washer_hole_ne', midpoint = [x4,y4], width = jj_size, orientation = 90)
    D_sq.add_port(name = 'washer_hole_nw', midpoint = [x5,y5], width = jj_size, orientation = 180)
    D_sq.add_port(name = 'washer_hole_sw', midpoint = [x6,y6], width = jj_size, orientation = 270)
    D_sq.add_port(name = 'washer_hole_se', midpoint = [x3,y3], width = jj_size, orientation = 0)
    D_sq.add_port(name = 'bias_point', midpoint = [two_half.xmax,D_sq.ports['jj_port_2'].midpoint[1]-sq_wash_jj_offset-w_wire/2-sq_bias_y_offset], width = w_wire, orientation = 0)
    shunt_coords = [[one_half.xmin,y6-sq_wash_w_wide+sq_shunt_w_lead/2],[two_half.xmin,two_half.ymin+2*sq_wash_w_wide+w1-sq_shunt_w_lead/2]]
    D_sq.add_port(name = 'shunt_1', midpoint = shunt_coords[0], width = sq_shunt_w_lead, orientation = 180)
    D_sq.add_port(name = 'shunt_2', midpoint = shunt_coords[1], width = sq_shunt_w_lead, orientation = 180)
           
    return D_sq


def sq_quad(params = dict()):
    
    sq_jj_junc_diam = vt_arg_helper(params,'sq_jj_junc_diam',1.15)
    sq_jj_shunt_resistance = vt_arg_helper(params,'sq_jj_shunt_resistance',0.93)
    sq_jj_shunt_w_wire = vt_arg_helper(params,'sq_jj_shunt_w_wire',5.5)
    sq_l_leads = vt_arg_helper(params,'sq_l_leads',11.4)
    sq_lead_gap = vt_arg_helper(params,'sq_lead_gap',6.3)
    sq_bias_y_coord = vt_arg_helper(params,'sq_bias_y_coord',7.2)
    sq_incoil_w_wire = vt_arg_helper(params,'sq_incoil_w_wire',2)
    sq_incoil_w_pitch = vt_arg_helper(params,'sq_incoil_w_pitch',2)
    sq_incoil_outset = vt_arg_helper(params,'sq_incoil_outset',0.6)
    sq_incoil_numturns = vt_arg_helper(params,'sq_incoil_numturns',5)    
    sq_incoil_distance_outside_washer = vt_arg_helper(params,'sq_incoil_distance_outside_washer',19.9)
    sq_washer_include_inductor_ports = vt_arg_helper(params,'sq_washer_include_inductor_ports',True)    
    sq_ground_plane_hole = vt_arg_helper(params,'sq_ground_plane_hole',False)
    sq_gph_outset = vt_arg_helper(params,'sq_gph_outset',10)
    sq_addflux_include = vt_arg_helper(params,'sq_addflux_include',True)
    sq_addflux_outset = vt_arg_helper(params,'sq_addflux_outset',2.7)
    sq_addflux_w_wire = vt_arg_helper(params,'sq_addflux_w_wire',2.01)
    sq_addflux_extent = vt_arg_helper(params,'sq_addflux_extent',21.1)
    sq_include_shunt = vt_arg_helper(params,'sq_include_shunt',True)
    sq_shunt_resistance = vt_arg_helper(params,'sq_shunt_resistance',0.1)
    sq_shunt_w_wire = vt_arg_helper(params,'sq_shunt_w_wire',6.4)
    sq_shunt_l_lead = vt_arg_helper(params,'sq_shunt_l_lead',1.4)
    sq_shunt_w_lead = vt_arg_helper(params,'sq_shunt_w_lead',9.4)
    sq_shunt_outset = vt_arg_helper(params,'sq_shunt_outset',1.3)
    via_width = vt_arg_helper(params,'via_width',2)
    via_top_contact_outset = vt_arg_helper(params,'via_top_contact_outset',2)
    w_wire = vt_arg_helper(params,'w_wire',2.03)
    device_label_size = vt_arg_helper(params,'device_label_size',2.03)
    induct_include_inductex_ports = vt_arg_helper(params,'induct_include_inductex_ports',False)
    label_layer = vt_arg_helper(params,'label_layer','m3')
    layers = vt_arg_helper(params,'layers',vt_lyrs)

    D_sq = Device('squid heart')
    
    D_wash = sq_washer_quad(params)
    washer = D_sq.add_ref(D_wash)
    
    #jjs
    if induct_include_inductex_ports == False:
        params_mod = copy.deepcopy(params)
        params_mod['jj_include_flux_moats'] = False
        params_mod['jj_shunt_resistance'] = sq_jj_shunt_resistance
        params_mod['jj_junc_diam'] = sq_jj_junc_diam
        params_mod['jj_shunt_w_wire'] = sq_jj_shunt_w_wire
        D_jj = jj_circle(params_mod)
        jj1 = D_sq.add_ref(D_jj)
        jj2 = D_sq.add_ref(D_jj)
        jj1.connect(port = 'jj1_south', destination = washer.ports['jj_port_1'])
        jj2.connect(port = 'jj1_south', destination = washer.ports['jj_port_2'])
        D_sq.add_ref(wire_basic(jj1.ports['m3_west'],jj2.ports['m3_west'],'xyx',w_wire,layers['m3']))
    
    #incoil bridge
    D_incoil = Device('sq_incoil')
    vias1 = D_incoil.add_ref(vt_m1_v3_via_sequence(params))
    vias1.connect(port = 'm1_west', destination = D_wash.ports['incoil_bridge_1'])
    vias2 = D_incoil.add_ref(vt_m1_v3_via_sequence(params))
    vias2.connect(port = 'm1_west', destination = D_wash.ports['incoil_bridge_2'])
    D_incoil.add_ref(wire_basic(vias1.ports['m1_west'],vias2.ports['m1_west'],'y',sq_incoil_w_wire,layers['m1']))
    
    #incoil
    p1 = vias1.ports['m3_south'].midpoint
    x_east = D_wash.ports['washer_hole_ne'].midpoint[0]+sq_incoil_outset+sq_incoil_w_wire/2
    y_north = D_wash.ports['washer_hole_ne'].midpoint[1]+sq_incoil_outset+sq_incoil_w_wire/2
    x_west = D_wash.ports['washer_hole_nw'].midpoint[0]-sq_incoil_outset-sq_incoil_w_wire/2
    y_south = D_wash.ports['washer_hole_sw'].midpoint[1]-sq_incoil_outset-sq_incoil_w_wire/2    
    p2 = [x_east,p1[1]]
    p3 = [x_east,y_north]
    p4 = [x_west,y_north]    
    points = [p1,p2,p3,p4]
    for ii in range(sq_incoil_numturns):
        y_north += sq_incoil_w_pitch
        points.append([x_west,y_south])#southwest
        x_west += -sq_incoil_w_pitch
        x_east += sq_incoil_w_pitch
        points.append([x_east,y_south])#southeast
        y_south += -sq_incoil_w_pitch
        points.append([x_east,y_north])#northeast
        points.append([x_west,y_north])#northwest    
    route_path = gdspy.PolyPath(points, width = sq_incoil_w_wire)
    D_incoil_half = Device('sq_incoil_half')
    D_incoil_half.add_polygon(route_path.polygons, layer = layers['m3'])
    D_incoil_half.add_port(name = 'incoil_half_inner', midpoint = p1, width = sq_incoil_w_wire, orientation = 180)
    D_incoil_half.add_port(name = 'incoil_half_outer', midpoint = points[-1], width = sq_incoil_w_wire, orientation = 0)
    
    incoil_north = D_incoil.add_ref(D_incoil_half)
    incoil_north.connect(port = 'incoil_half_inner', destination = vias1.ports['m3_south'])
    temp_length = washer.xmin-sq_l_leads-incoil_north.ports['incoil_half_outer'].midpoint[0]
    p1 = np.array([incoil_north.ports['incoil_half_outer'].midpoint[0],incoil_north.ports['incoil_half_outer'].midpoint[1]])
    D_incoil.add_ref(wire_basic(p1,p1+np.array([temp_length,0]),'xyx',sq_incoil_w_wire,layers['m3']))
    D_sq.add_port(name = 'incoil_1', midpoint = p1+np.array([temp_length,0]),width = sq_incoil_w_wire, orientation = 180)
    
    incoil_south = D_incoil.add_ref(D_incoil_half)
    incoil_south.connect(port = 'incoil_half_inner', destination = vias2.ports['m3_south'])
    p1 = incoil_south.ports['incoil_half_outer'].midpoint
    p2 = D_sq.ports['incoil_1'].midpoint
    temp_length = washer.xmax+sq_incoil_distance_outside_washer+sq_incoil_w_wire/2-p1[0]
    points_wrap_back = np.array([p1,[p1[0]+temp_length,p1[1]],[p1[0]+temp_length,p2[1]+sq_lead_gap+sq_incoil_w_wire],[p2[0],p2[1]+sq_lead_gap+sq_incoil_w_wire]])
    route_path = gdspy.PolyPath(points_wrap_back, width = sq_incoil_w_wire)
    D_incoil.add_polygon(route_path.polygons, layer = layers['m3'])    
    D_sq.add_port(name = 'incoil_2', midpoint = [p2[0],p2[1]+sq_lead_gap+sq_incoil_w_wire],width = sq_incoil_w_wire, orientation = 180)
    
    D_sq.add_ref(D_incoil)
    
    #addflux
    if sq_addflux_include == True:
        D_addflux = Device('squid_flux_adder')
        x_coord = washer.x
        y_coord = D_incoil.ymin-sq_addflux_outset-sq_addflux_w_wire/2
        wire_vert = wire_basic([x_coord-sq_addflux_extent/2,y_coord],[x_coord+sq_addflux_extent/2,y_coord],'x',sq_addflux_w_wire,layers['m3'])
        D_sq.add_ref(wire_vert)
        D_sq.add_port(name = 'addflux_east', midpoint = [x_coord+sq_addflux_extent/2,y_coord], width = sq_addflux_w_wire, orientation = 0)
        D_sq.add_port(name = 'addflux_southeast', midpoint = [x_coord+sq_addflux_extent/2-sq_addflux_w_wire/2,y_coord-sq_addflux_w_wire/2], width = sq_addflux_w_wire, orientation = 270)
        D_sq.add_port(name = 'addflux_west', midpoint = [x_coord-sq_addflux_extent/2,y_coord], width = sq_addflux_w_wire, orientation = 180)    
        D_sq.add_port(name = 'addflux_southwest', midpoint = [x_coord-sq_addflux_extent/2+sq_addflux_w_wire/2,y_coord--sq_addflux_w_wire/2], width = sq_addflux_w_wire, orientation = 270)
        D_sq.add_ref(D_addflux)
    
    #bias leads
    if induct_include_inductex_ports == False:
        Cf = corner_fixer(w_wire,layers['m3'])
        cf = D_sq.add_ref(Cf)
        cf.connect(port = 'south', destination = jj2.ports['m3_north'])
        p1 = cf.ports['north'].midpoint
        p2 = [D_wash.ports['bias_point'].midpoint[0]+sq_l_leads,jj2.ymax+2*w_wire+w_wire/2]
        D_sq.add_ref(wire_basic(p1,p2,'yx',w_wire,layers['m3']))
        D_sq.add_port(name = 'bias_1', midpoint = p2, width = w_wire, orientation = 0)
        D_sq.add_ref(wire_basic(D_wash.ports['bias_point'],[D_wash.ports['bias_point'].midpoint[0]+sq_l_leads,D_wash.ports['bias_point'].midpoint[1]],'x',w_wire,layers['jj1']))
        D_sq.add_port(name = 'bias_2', midpoint = [D_wash.ports['bias_point'].midpoint[0]+sq_l_leads,D_wash.ports['bias_point'].midpoint[1]],width = w_wire, orientation = 0)
    
    #shunt resistor
    if sq_include_shunt == True:
        D_res = Device('sq_shunt_resistor')
        params_mod = copy.deepcopy(params)
        params_mod['res_resistance'] = sq_shunt_resistance
        params_mod['res_w_wire'] = sq_shunt_w_wire
        params_mod['res_l_lead'] = sq_shunt_l_lead
        params_mod['w_wire'] = sq_shunt_w_lead
        params_mod['res_layer'] = 'r2'
        params_mod['res_include_label'] = False
        res = D_res.add_ref(res_stitch_simp(params_mod))
        res.rotate(90)    
        x_coord = min(washer.ports['shunt_1'].midpoint[0],washer.ports['shunt_2'].midpoint[0])-res.xsize/2-sq_shunt_outset
        y_coord = washer.ports['shunt_2'].midpoint[1]+(washer.ports['shunt_1'].midpoint[1]-washer.ports['shunt_2'].midpoint[1])/2
        res.center = [x_coord,y_coord]
        params_mod = copy.deepcopy(params)
        params_mod['via_layer'] = 'v3'
        params_mod['layer_below'] = 'jj1'
        params_mod['layer_above'] = 'm3'
        Via = via_general(params_mod)
        via1 = D_sq.add_ref(Via) 
        via1.connect(port = 'above_south', destination = res.ports['east']) 
        via2 = D_sq.add_ref(Via) 
        via2.connect(port = 'above_north', destination = res.ports['west'])
        D_res.add_ref(wire_basic(via1.ports['below_east'],washer.ports['shunt_1'],'xyx',sq_shunt_w_lead,layers['jj1']))       
        D_res.add_ref(wire_basic(via2.ports['below_east'],washer.ports['shunt_2'],'xyx',sq_shunt_w_lead,layers['jj1']))       
        D_sq.add_ref(D_res)
    
    #label
    if sq_include_shunt == True:
        text_string = 'sq_quad'
    elif sq_include_shunt == False:
        text_string = 'sq_quad'
    
    Text_label = vt_label_maker(text_string = text_string, text_size = device_label_size, justify = 'center', layer = layers[label_layer])
    text_coords = [D_sq.xmin-Text_label.ysize,washer.ports['shunt_2'].midpoint[1]+(washer.ports['shunt_1'].midpoint[1]-washer.ports['shunt_2'].midpoint[1])/2]
    text_label = D_sq.add_ref(Text_label)
    text_label.rotate(270)
    text_label.center = text_coords
    
    #inductor ports
    if induct_include_inductex_ports == True:
        
        Inductor_port_bias_jj1 = pg.rectangle(size = [w_wire+1,2], layer = layers['ipj1'])
        Inductor_port_bias_m3 = pg.rectangle(size = [w_wire+1,2], layer = layers['ipm3'])
        in_po1 = D_sq.add_ref(Inductor_port_bias_jj1)
        in_po2 = D_sq.add_ref(Inductor_port_bias_m3)
        in_po1.center = washer.ports['jj_port_1'].midpoint
        in_po2.center = washer.ports['jj_port_2'].midpoint
        D_sq.label(text = 'P1+ jj1', position = washer.ports['jj_port_1'].midpoint, layer = layers['ipl'])
        D_sq.label(text = 'P1- jj1', position = washer.ports['jj_port_2'].midpoint, layer = layers['ipl'])
        
        Inductor_port_incoil = pg.rectangle(size = [sq_incoil_w_wire+1,2], layer = layers['ipm3'])
        in_po1 = D_sq.add_ref(Inductor_port_incoil)
        in_po2 = D_sq.add_ref(Inductor_port_incoil)
        in_po1.center = D_sq.ports['incoil_1'].midpoint
        in_po2.center = D_sq.ports['incoil_2'].midpoint
        D_sq.label(text = 'P2+ m3', position = D_sq.ports['incoil_1'].midpoint, layer = layers['ipl'])
        D_sq.label(text = 'P2- m3', position = D_sq.ports['incoil_2'].midpoint, layer = layers['ipl'])  
                
        Inductor_port_addflux = pg.rectangle(size = [sq_addflux_w_wire+1,2], layer = layers['ipm3'])
        in_po1 = D_sq.add_ref(Inductor_port_addflux)
        in_po2 = D_sq.add_ref(Inductor_port_addflux)
        in_po1.center = D_sq.ports['addflux_east'].midpoint
        in_po2.center = D_sq.ports['addflux_west'].midpoint
        D_sq.label(text = 'P3+ m3', position = D_sq.ports['addflux_east'].midpoint, layer = layers['ipl'])
        D_sq.label(text = 'P3- m3', position = D_sq.ports['addflux_west'].midpoint, layer = layers['ipl'])        
    
    return D_sq

def sq_quad_8_wire(params = dict()):
    
    sq_incoil_w_wire = vt_arg_helper(params,'sq_incoil_w_wire',2.02)
    sq_addflux_w_wire = vt_arg_helper(params,'sq_addflux_w_wire',2.3)
    sq_res_50_ohm_lateral_gap = vt_arg_helper(params,'sq_res_50_ohm_lateral_gap',12.3)
    pad_size = vt_arg_helper(params,'pad_size',[75,123])
    pad_size_ground = vt_arg_helper(params,'pad_size_ground',[44,33])
    pad_pitch = vt_arg_helper(params,'pad_pitch',[206,398])
    pad_y_backset = vt_arg_helper(params,'pad_y_backset',77)
    pad_gnd_y_backset = vt_arg_helper(params,'pad_gnd_y_backset',33)
    pad_w_wire = vt_arg_helper(params,'pad_w_wire',2.1)
    w_wire = vt_arg_helper(params,'w_wire',1.9)
    layers = vt_arg_helper(params,'layers',vt_lyrs)

    D_sq = Device('squid heart')
    
    D_squid_heart = sq_quad(params)
    sq = D_sq.add_ref(D_squid_heart)
    sq.rotate(90)
    
    #lead contacts
    Cf1 = corner_fixer(sq_incoil_w_wire,layers['m3'])
    Cf2 = corner_fixer(w_wire,layers['m3'])
    Cf3 = corner_fixer(w_wire,layers['jj1'])
    cf_incoil1 = D_sq.add_ref(Cf1)
    cf_incoil1.connect(port = 'north', destination = sq.ports['incoil_2'])
    cf_bias1 = D_sq.add_ref(Cf3)
    cf_bias1.connect(port = 'north', destination = sq.ports['bias_2'])
    
    #50-ohm resistors
    Res = res_50_ohm(params)
    res_in = D_sq.add_ref(Res)
    res_in.connect(port = 'east', destination = cf_incoil1.ports['west'])
    resistor_length = res_in.xsize
    res_addflux = D_sq.add_ref(Res)
    res_addflux.connect(port = 'west', destination = sq.ports['addflux_southeast'])
    
    #wires for I-V
    D_res = Device('sq_res_contact_structure')
    params_mod = copy.deepcopy(params)
    params_mod['via_layer'] = 'v3'
    params_mod['via_layer_below'] = 'jj1'
    params_mod['via_layer_above'] = 'm3'
    via1 = D_sq.add_ref(via_general(params_mod)) 
    via1.connect(port = 'above_south', destination = sq.ports['bias_1'])
    wire_r1 = D_res.add_ref(wire_basic(via1.ports['below_west'],[sq.xmin-2*w_wire,via1.ports['below_west'].midpoint[1]],'x',w_wire,layers['jj1']))
    via2 = D_sq.add_ref(via_general(params_mod))
    via2.connect(port = 'below_east', destination = via1.ports['below_west'])
    via2.movex(-wire_r1.xsize)
    res_sq1 = D_res.add_ref(Res)
    res_sq1.connect(port = 'east', destination = via2.ports['above_west'])
    res_sq2 = D_res.add_ref(Res)
    res_sq2.connect(port = 'east', destination = via2.ports['above_west'])
    res_sq2.movey(sq_res_50_ohm_lateral_gap)
    D_res.add_ref(wire_basic(via2.ports['above_north'],res_sq2.ports['east'],'yx',w_wire,layers['m3']))
    D_sq.add_ref(D_res)
    
    #pads
    D_pad = jj_pad(params)
    params_mod = copy.deepcopy(params)
    params_mod['is_ground_pad'] = True
    D_pad_gnd = jj_pad(params_mod)    
    center_coord = [sq.x,sq.y]
    
    pad_x_space = pad_pitch[0]-pad_size[0]
    
    #squid bias 1
    pad1 = D_sq.add_ref(D_pad)
    extra_backset_local = 20#20+resistor_length
    pad1.center = [center_coord[0]-pad_x_space/2-5*pad_size[0]/2-pad_pitch[0],sq.ymin-pad_size[1]/2-pad_y_backset-extra_backset_local]
    D_sq.add_ref(wire_tri_seg(p1 = pad1.ports['m3_north'], p2 = res_sq2.ports['west'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.4,0.1], directions = 'yx', layer = layers['m3']))
      
    pad2 = D_sq.add_ref(D_pad)
    pad2.center = [pad1.center[0]+pad_pitch[0],pad1.center[1]]
    D_sq.add_ref(wire_tri_seg(p1 = pad2.ports['m3_north'], p2 = res_sq1.ports['west'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.2,0.1], directions = 'yx', layer = layers['m3']))
        
    pad3 = D_sq.add_ref(D_pad)    
    pad3.center = [pad2.center[0]+pad_pitch[0],pad2.center[1]]
    cf_sq_in = D_sq.add_ref(Cf2)
    cf_sq_in.connect(port = 'east', destination = res_in.ports['west'])
    D_sq.add_ref(wire_tri_seg(p1 = pad3.ports['m3_north'], p2 = cf_sq_in.ports['south'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.2,0.1], directions = 'yx', layer = layers['m3']))

    pad4 = D_sq.add_ref(D_pad)
    pad4.center = [pad3.center[0]+pad_pitch[0],pad3.center[1]]
    D_sq.add_ref(wire_tri_seg(p1 = pad4.ports['m3_north'], p2 = res_addflux.ports['east'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.7,0.1], directions = 'yx', layer = layers['m3']))

    #ground pads
    pad1_gnd = D_sq.add_ref(D_pad_gnd)
    pad1_gnd.center = sq.ports['incoil_1'].midpoint
    pad1_gnd.move([0,-pad_size_ground[1]/2-pad_gnd_y_backset])
    D_sq.add_ref(wire_tri_seg(p1 = pad1_gnd.ports['m3_north'], p2 = sq.ports['incoil_1'], initial_width = pad_w_wire, final_width = sq_incoil_w_wire, length_factors = [0.1,0.2,0.05], directions = 'yy', layer = layers['m3']))
    pad2_gnd = D_sq.add_ref(D_pad_gnd)
    pad2_gnd.center = cf_bias1.ports['south'].midpoint
    pad2_gnd.move([0,pad_size_ground[1]/2+pad_gnd_y_backset])
    D_sq.add_ref(wire_tri_seg(p1 = pad2_gnd.ports['jj1_north'], p2 = cf_bias1.ports['north'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.1,0.2,0.05], directions = 'yy', layer = layers['jj1']))
    pad3_gnd = D_sq.add_ref(D_pad_gnd)
    pad3_gnd.center = [sq.ports['addflux_west'].midpoint[0],pad1_gnd.center[1]]
    D_sq.add_ref(wire_tri_seg(p1 = pad3_gnd.ports['m3_north'], p2 = sq.ports['addflux_west'], initial_width = pad_w_wire, final_width = sq_addflux_w_wire, length_factors = [0.1,0.2,0.05], directions = 'yy', layer = layers['m3']))
    
    D_sq.add_port(name = 'pad_anchor', midpoint = pad1.center, width = pad_w_wire, orientation = 90)
    
    return D_sq

def sq_inductance_test(params = dict()):
    
    D_sq = Device('squid inductance test')
    
    params['sq_incoil_w_coil'] = params['sq_wash_w_in']+2*(params['sq_incoil_outset']+params['sq_incoil_w_wire'])    
    washer = D_sq.add_ref(sq_washer(params))
    incoil = D_sq.add_ref(sq_input_coil_one_layer(params))
    washer.connect(port = 'washer_center', destination = incoil.ports['coil_center'])
    
    return D_sq