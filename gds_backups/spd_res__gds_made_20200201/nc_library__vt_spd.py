import copy
import numpy as np
import gdspy
import phidl, phidl.geometry as pg
from phidl import Device, Layer, LayerSet, Port
from phidl import make_device

from nc_library import wire_basic
from vt_util import vt_layers

from nc_library__vt_util import vt_arg_helper, corner_fixer, vt_label_maker, vt_fiber_collar
from nc_library__vt_pads_vias_wires import jj_pad, jj_jj1_m3_via, wire_tri_seg, vt_m1_v3_via_sequence, vt_m1_jj1_via_sequence, via_general
from nc_library__vt_res import res_stitch_simp, res_50_ohm
from nc_library__vt_sq import sq_quad

vt_lyrs,layer_data = vt_layers()

def vt_spd_meander(params = dict()):
        
    spd_w_wire = vt_arg_helper(params,'spd_w_wire',1.0)
    spd_wire_pitch = vt_arg_helper(params,'spd_wire_pitch',2.0)
    spd_inductance = vt_arg_helper(params,'spd_inductance',100e-9)
    spd_turn_ratio = vt_arg_helper(params,'spd_turn_ratio',4)
    spd_ground_plane_buffer = vt_arg_helper(params,'spd_ground_plane_buffer',9.6)
    spd_pad_w_wire = vt_arg_helper(params,'spd_pad_w_wire',21.6)
    inductance_per_sq_stf = vt_arg_helper(params,'inductance_per_sq_stf',200e-12)
    pad_size = vt_arg_helper(params,'pad_size',[200,250])
    pad_y_backset = vt_arg_helper(params,'pad_y_backset',100)
    pad_w_wire = vt_arg_helper(params,'pad_w_wire',20)
    device_label_size = vt_arg_helper(params,'device_label_size',5)
    label_layer = vt_arg_helper(params,'label_layer','m4l')
    layers = vt_arg_helper(params,'layers',vt_lyrs)
    
    D_spd = Device('spd')
    
    #spd meander
    spd_num_squares = spd_inductance/inductance_per_sq_stf
    meander = D_spd.add_ref(pg.snspd(wire_width = spd_w_wire, wire_pitch = spd_wire_pitch, size = None, num_squares = spd_num_squares, turn_ratio = spd_turn_ratio, terminals_same_side = False, layer = layers['stf'])).rotate(90)

    #pads and wires    
    params_mod = copy.deepcopy(params)   
    Pad = jj_pad(params_mod)
    params_mod['is_ground_pad'] = True    
    Pad_gnd = jj_pad(params_mod)
    
    pad1 = D_spd.add_ref(Pad)
    pad1.connect(port = 'm1p_north', destination = meander.ports[1])
    pad1.move([0,-pad_y_backset])    
    wire1 = wire_tri_seg(p1 = pad1.ports['stfp_north'], p2 = meander.ports[1], initial_width = spd_pad_w_wire, final_width = spd_w_wire, length_factors = [0.75,0.1,0.1], directions = 'yy', layer = layers['stf'])
    D_spd.add_ref(wire1)
    
    pad2 = D_spd.add_ref(Pad_gnd)
    pad2.connect(port = 'm1p_south', destination = meander.ports[2])
    pad2.move([0,pad_y_backset/2])
    wire2 = wire_tri_seg(p1 = pad2.ports['stfp_south'], p2 = meander.ports[2], initial_width = spd_pad_w_wire, final_width = spd_w_wire, length_factors = [0.75,0.1,0.1], directions = 'yy', layer = layers['stf'])
    D_spd.add_ref(wire2)
    
    #ground plane clearout
    Gp_clearout = pg.rectangle(size = [meander.xsize+2*spd_ground_plane_buffer,meander.ysize+2*spd_ground_plane_buffer], layer = layers['m2i'])
    gpc = D_spd.add_ref(Gp_clearout)
    gpc.center = meander.center
    
    #fiber core
    fiber_core = D_spd.add_ref(pg.circle(radius = 10, angle_resolution = 360/100, layer = layers['pkfc']))
    fiber_core.center = meander.center
                       
    #text label
    string = 'spd meander\nw_wire = '+str(spd_w_wire)+'um\npitch = '+str(spd_wire_pitch)+'um\nnum_sq = '+str(int(round(spd_num_squares)))
    Text_label = vt_label_maker(text_string = string, text_size = device_label_size, justify = 'right', layer = layers[label_layer])
    text_label = D_spd.add_ref(Text_label)
    text_label.center = ([meander.xmin-text_label.xsize/2-10,gpc.y])
    
    #pad anchor
    D_spd.add_port(name = 'pad_anchor', midpoint = [pad1.ports['m1p_north'].midpoint[0],pad1.ports['m1p_west'].midpoint[1]], width = pad_size[0], orientation = 90)

    return D_spd


def vt_sy_spd(params = dict()):
    
    spd_w_wire = vt_arg_helper(params,'spd_w_wire',1.0)
    spd_wire_pitch = vt_arg_helper(params,'spd_wire_pitch',2.0)
    spd_inductance = vt_arg_helper(params,'spd_inductance',100e-9)
    spd_turn_ratio = vt_arg_helper(params,'spd_turn_ratio',4)
    spd_ground_plane_buffer = vt_arg_helper(params,'spd_ground_plane_buffer',9.6)
    spd_sy_l_lead = vt_arg_helper(params,'spd_sy_l_lead',10.1)
    spd_sy_w_lead = vt_arg_helper(params,'spd_sy_w_lead',5.2)
    spd_sy_m1_olap = vt_arg_helper(params,'spd_sy_m1_olap',4.9)
    spd_sy_m1_outset = vt_arg_helper(params,'spd_sy_m1_outset',1.1)
    spd_sy_m1_res_l_extra = vt_arg_helper(params,'spd_sy_m1_res_l_extra',11.1)
    spd_sy_include_res = vt_arg_helper(params,'spd_sy_include_res',True) 
    spd_res_y_offset = vt_arg_helper(params,'spd_res_y_offset',19.7)
    spd_contact_taper_length = vt_arg_helper(params,'spd_contact_taper_length',11.7)     
    sfg_has_sfq = vt_arg_helper(params,'sfg_has_sfq',False)
    jj_res = vt_arg_helper(params,'jj_res',4.125)
    w_wire = vt_arg_helper(params,'w_wire',2.7)
    sy_res_w_wire = vt_arg_helper(params,'sy_res_w_wire',1.9)
    spd_tau = vt_arg_helper(params,'spd_tau',51e-9)
    inductance_per_sq_stf = vt_arg_helper(params,'inductance_per_sq_stf',200e-12)
    device_label_size = vt_arg_helper(params,'device_label_size',5)
    label_layer = vt_arg_helper(params,'label_layer','m4l')
    layers = vt_arg_helper(params,'layers',vt_lyrs)
    
    D_spd = Device('synapse_spd')
        
    #spd meander
    spd_num_squares = spd_inductance/inductance_per_sq_stf
    meander = D_spd.add_ref(pg.snspd(wire_width = spd_w_wire, wire_pitch = spd_wire_pitch, size = None, num_squares = spd_num_squares, turn_ratio = spd_turn_ratio, terminals_same_side = False, layer = layers['stf'])).rotate(90)
    D_spd.add_port(name = 'fiber_port', midpoint = meander.center, width = meander.xsize, orientation = 0)
    fiber_core = D_spd.add_ref(pg.circle(radius = 10, angle_resolution = 360/100, layer = layers['pkfc']))
    fiber_core.center = meander.center
    
    #spd contacts
    Cntct = Device('spd_sy_cntct')
    Cntct.add_ref(wire_tri_seg(p1 = [0,0], p2 = [0,spd_sy_l_lead+spd_sy_m1_olap], initial_width = spd_w_wire, final_width = spd_sy_w_lead, length_factors = [spd_sy_l_lead/2,spd_sy_l_lead/2,spd_sy_m1_olap], directions = 'yy', layer = layers['stf']))
    rect = Cntct.add_ref(pg.rectangle(size = [spd_sy_w_lead+2*spd_sy_m1_outset,2*spd_res_y_offset], layer = layers['m1']))
    rect.move([-spd_sy_w_lead/2-spd_sy_m1_outset,spd_sy_l_lead])
    Cntct.add_port(name = 'stf', midpoint = [0,0], width = spd_w_wire, orientation = 270)
    Cntct.add_port(name = 'm1', midpoint = [0,spd_sy_l_lead+2*spd_res_y_offset], width = spd_sy_w_lead+2*spd_sy_m1_outset, orientation = 90)
    Cntct.add_port(name = 'res_port', midpoint = [spd_sy_w_lead/2+spd_sy_m1_outset,spd_sy_l_lead+spd_res_y_offset], width = w_wire, orientation = 0)
    cntct1 = D_spd.add_ref(Cntct)
    cntct2 = D_spd.add_ref(Cntct)
    cntct1.connect(port = 'stf', destination = meander.ports[2])
    cntct2.connect(port = 'stf', destination = meander.ports[1])
    
    if spd_sy_include_res == True:
        
        #taper
        tpr = D_spd.add_ref(pg.taper(length = spd_contact_taper_length, width1 = spd_sy_w_lead+2*spd_sy_m1_outset, width2 = w_wire, port = None, layer = layers['m1'])).rotate(90)
        tpr.connect(port = 1, destination = cntct1.ports['m1'])
        
        #vias to room temp
        vias1 = D_spd.add_ref(vt_m1_v3_via_sequence(params))
        vias1.connect(port = 'm1p_west', destination = tpr.ports[2])
        
        #resistor 50 ohms from room temp
        res1 = D_spd.add_ref(res_50_ohm(params))
        res1.connect(port = 'west', destination = vias1.ports['m3_east'])
     
        #resistor to jj/sq
        res_resistance = spd_inductance/spd_tau-jj_res
        params_mod = copy.deepcopy(params)
        params_mod['res_w_wire'] = sy_res_w_wire
        params_mod['res_resistance'] = res_resistance
        params_mod['res_layer'] = 'r1'
        res2 = D_spd.add_ref(res_stitch_simp(params_mod))
        res2.connect(port = 'west', destination = cntct1.ports['res_port'])
        res2.movex(spd_sy_m1_res_l_extra)
        D_spd.add_ref(wire_basic(cntct1.ports['res_port'],res2.ports['west'],'x',w_wire,layers['m1']))
        
        #vias to jj
        if sfg_has_sfq == False:
            vias2 = D_spd.add_ref(vt_m1_v3_via_sequence(params))
        elif sfg_has_sfq == True:
            vias2 = D_spd.add_ref(vt_m1_jj1_via_sequence(params))
        vias2.connect(port = 'm1p_west', destination = res2.ports['east'])
        

    #ground plane clearout
    Gp_clearout = pg.rectangle(size = [meander.xsize+2*spd_ground_plane_buffer,meander.ysize+2*spd_ground_plane_buffer], layer = layers['m2i'])
    gpc = D_spd.add_ref(Gp_clearout)
    gpc.center = meander.center   
    
    #fiber collar
    min_coords = [D_spd.xmin,D_spd.ymin]
    max_coords = [D_spd.xmax,D_spd.ymax]
    fc = D_spd.add_ref(vt_fiber_collar(params))
    fc.connect(port = 'fiber_center', destination = D_spd.ports['fiber_port'])
    
    if spd_sy_include_res == True:
        D_spd.add_port(name = 'pad', midpoint = res1.ports['east'].midpoint, width = w_wire, orientation = 90)
    elif spd_sy_include_res == False:
        D_spd.add_port(name = 'pad', midpoint = cntct1.ports['m1'].midpoint, width = spd_sy_w_lead+2*spd_sy_m1_outset, orientation = 90)
    D_spd.add_port(name = 'pad_gnd', midpoint = cntct2.ports['m1'].midpoint, width = spd_sy_w_lead+2*spd_sy_m1_outset, orientation = 270)
    if spd_sy_include_res == True:
        if sfg_has_sfq == False:
            D_spd.add_port(name = 'wire_port', midpoint = vias2.ports['m3_east'].midpoint, width = w_wire, orientation = 0)
        elif sfg_has_sfq == True:
            D_spd.add_port(name = 'wire_port', midpoint = vias2.ports['jj1_east'].midpoint, width = w_wire, orientation = 0)
    D_spd.add_port(name = 'min_coords_without_collar', midpoint = min_coords, width = w_wire, orientation = 0)
    D_spd.add_port(name = 'max_coords_without_collar', midpoint = max_coords, width = w_wire, orientation = 0)
    D_spd.add_port(name = 'fiber_center', midpoint = fc.ports['fiber_center'].midpoint, width = w_wire, orientation = 270)
                       
    #text label
    string = 'spd meander\nw_wire = '+str(spd_w_wire)+' um\npitch = '+str(spd_wire_pitch)+' um\nnum_sq = '+str(int(round(spd_num_squares)))+'\ninductance = '+str(int(round(spd_inductance*1e9)))+' nH'
    Text_label = vt_label_maker(text_string = string, text_size = device_label_size, justify = 'right', layer = layers[label_layer])
    text_label = D_spd.add_ref(Text_label)
    text_label.center = [gpc.xmin-text_label.xsize/2-device_label_size,gpc.y]    
    
    return D_spd


def vt_spd_sq(params = dict()):
    
    w_wire = vt_arg_helper(params,'w_wire',1.3)
    spd_sq_extra_x = vt_arg_helper(params,'spd_sq_extra_x',37.1)
    spd_tau = vt_arg_helper(params,'spd_tau',100e-9)
    spd_cntrl_inductance = vt_arg_helper(params,'spd_cntrl_inductance',1e-6) 
    sq_incoil_w_wire = vt_arg_helper(params,'sq_incoil_w_wire',1.9)
    sq_addflux_w_wire = vt_arg_helper(params,'sq_addflux_w_wire',2.5)
    pad_w_wire = vt_arg_helper(params,'pad_w_wire',21.3)
    pad_y_backset = vt_arg_helper(params,'pad_y_backset',27.83)
    pad_pitch = vt_arg_helper(params,'pad_pitch',[343,445])
    pad_size = vt_arg_helper(params,'pad_size',[44,33])
    pad_size_ground = vt_arg_helper(params,'pad_size_ground',[44,33])
    pad_gnd_y_backset = vt_arg_helper(params,'pad_gnd_y_backset',29.83)
    device_label_size = vt_arg_helper(params,'device_label_size',2.03)
    label_layer = vt_arg_helper(params,'label_layer','m3')
    layers = vt_arg_helper(params,'layers',vt_lyrs)    
    
    D_spd = Device('sfg')
    
    #spd
    params_mod = copy.deepcopy(params)
    params_mod['fiber_collar_orientation'] = 0
    spd = D_spd.add_ref(vt_sy_spd(params_mod))
    
    #squid
    sq = D_spd.add_ref(sq_quad(params))
    sq.connect(port = 'incoil_2', destination = spd.ports['wire_port'])
    sq.movex(spd_sq_extra_x)
    D_spd.add_ref(wire_basic(sq.ports['incoil_2'],spd.ports['wire_port'],'x',sq_incoil_w_wire,layers['m3']))
    params_mod = copy.deepcopy(params)
    params_mod['via_layer'] = 'v3'
    params_mod['layer_below'] = 'jj1'
    params_mod['layer_above'] = 'm3'
    Via_sq = via_general(params_mod)
    sq_via1 = D_spd.add_ref(Via_sq)
    sq_via2 = D_spd.add_ref(Via_sq)
    sq_via1.connect(port = 'above_west', destination = sq.ports['bias_1'])
    sq_via2.center = [sq_via1.center[0],sq.ymax+4*sq_incoil_w_wire]
    D_spd.add_ref(wire_basic(sq_via1.ports['below_north'],sq_via2.ports['below_south'],'y',w_wire,layers['jj1']))
    res_50_1 = D_spd.add_ref(res_50_ohm(params))
    res_50_1.connect(port = 'west', destination = sq_via2.ports['above_north'])
    res_50_2 = D_spd.add_ref(res_50_ohm(params))
    res_50_2.connect(port = 'west', destination = sq_via2.ports['above_east'])
    
    #####
    #pads
    #####
    Pad = jj_pad(params)
    params_mod = copy.deepcopy(params)
    params_mod['is_ground_pad'] = True
    Pad_gnd = jj_pad(params_mod)
    
    #spd pads
    pad_spd = D_spd.add_ref(Pad)
    initial_x_offset = -250
    initial_y_offset = 25
    pad_spd.connect(port = 'm1p_south', destination = spd.ports['pad'])
    pad_spd.move([initial_x_offset,initial_y_offset+pad_y_backset])
    D_spd.add_ref(wire_tri_seg(p1 = pad_spd.ports['m3_south'], p2 = spd.ports['pad'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [5,4,1], directions = 'yy', layer = layers['m3']))
    pad_spd_gnd = D_spd.add_ref(Pad_gnd)
    pad_spd_gnd.connect(port = 'm1p_north', destination = spd.ports['pad_gnd'])            

    #mi loop ground
    pad_mi_gnd = D_spd.add_ref(Pad_gnd)
    pad_mi_gnd.center = (sq.ports['incoil_1'].midpoint)
    pad_mi_gnd.move([-pad_size_ground[0]/2,-pad_size_ground[1]/2-4*sq_incoil_w_wire])
    D_spd.add_ref(wire_basic(pad_mi_gnd.ports['m3_north'],sq.ports['incoil_1'],'yx',sq_incoil_w_wire,layers['m3']))    
    
    #sq
    pad_sq_I = D_spd.add_ref(Pad)
    pad_sq_I.center = [pad_spd.center[0]+pad_pitch[0],pad_spd.center[1]]
    D_spd.add_ref(wire_tri_seg(p1 = pad_sq_I.ports['m3_south'], p2 = res_50_1.ports['east'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.2,0.5,0.2], directions = 'yy', layer = layers['m3']))
    pad_sq_V = D_spd.add_ref(Pad)
    pad_sq_V.center = [pad_sq_I.center[0]+pad_pitch[0],pad_spd.center[1]]
    D_spd.add_ref(wire_tri_seg(p1 = pad_sq_V.ports['m3_south'], p2 = res_50_2.ports['east'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.5,0.1], directions = 'yx', layer = layers['m3']))    
    pad_sq_gnd = D_spd.add_ref(Pad_gnd)
    pad_sq_gnd.center = sq.ports['bias_2'].midpoint
    pad_sq_gnd.movex(pad_size_ground[0]/2+pad_gnd_y_backset)
    D_spd.add_ref(wire_tri_seg(p1 = pad_sq_gnd.ports['jj1_west'], p2 = sq.ports['bias_2'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.4,0.6,0.1], directions = 'xx', layer = layers['jj1']))    
    
    #addflux
    pad_flux = D_spd.add_ref(Pad)
    pad_flux.center = [pad_sq_V.center[0]+pad_pitch[0],pad_spd.center[1]]
    D_spd.add_ref(wire_tri_seg(p1 = pad_flux.ports['m3_south'], p2 = sq.ports['addflux_east'], initial_width = pad_w_wire, final_width = sq_addflux_w_wire, length_factors = [0.6,0.05], directions = 'yx', layer = layers['m3']))    
    pad_flux_gnd = D_spd.add_ref(Pad_gnd)
    pad_flux_gnd.center = sq.ports['addflux_west'].midpoint
    pad_flux_gnd.move([-pad_size_ground[0]/2,-pad_size_ground[1]/2-4*sq_incoil_w_wire])
    D_spd.add_ref(wire_basic(pad_flux_gnd.ports['m3_north'],sq.ports['addflux_west'],'yx',sq_addflux_w_wire,layers['m3']))
    
    #control spd
    params_mod = copy.deepcopy(params)
    params_mod['fiber_collar_orientation'] = 180
    params_mod['spd_inductance'] = spd_cntrl_inductance
    params_mod['spd_sy_include_res'] = False
    spd_cntrl = D_spd.add_ref(vt_sy_spd(params_mod))
    pad_spd_cntrl = D_spd.add_ref(Pad)
    pad_spd_cntrl.center = [pad_flux.center[0]+pad_pitch[0],pad_spd.center[1]]
    spd_cntrl.center = pad_spd_cntrl.center
    spd_cntrl.move([spd_cntrl.center[0]-spd_cntrl.ports['fiber_center'].midpoint[0],-pad_size[1]/2-spd_cntrl.ysize/2-pad_y_backset/2])
    D_spd.add_ref(wire_tri_seg(p1 = pad_spd_cntrl.ports['m1p_south'], p2 = spd_cntrl.ports['pad'], initial_width = pad_w_wire, final_width = spd_cntrl.ports['pad'].width, length_factors = [2,5,1], directions = 'yy', layer = layers['m1']))    
    pad_spd_cntrl_gnd = D_spd.add_ref(Pad_gnd)
    pad_spd_cntrl_gnd.connect(port = 'm1p_north', destination = spd_cntrl.ports['pad_gnd'])
    
    #pad anchor
    D_spd.add_port(name = 'pad_anchor', midpoint = [pad_spd.ports['m3_north'].midpoint[0],pad_spd.ports['m3_west'].midpoint[1]], width = pad_w_wire, orientation = 270)
    
    #label
    Text_label = Device('spd_sq_text_label')
    text_label1 = Text_label.add_ref(vt_label_maker(text_string = 'SPD_SQ', text_size = 1.4*device_label_size, justify = 'center', layer = layers[label_layer]))
    text_string = '\ntau_spd = '+str(int(round(spd_tau*1e9)))+' ns'
    text_label2 = Text_label.add_ref(vt_label_maker(text_string = text_string, text_size = device_label_size, justify = 'center', layer = layers[label_layer]))
    text_coords = [pad_spd_gnd.xmax+(pad_flux_gnd.xmin-pad_spd_gnd.xmax)/2,pad_flux_gnd.y+(pad_spd_gnd.y-pad_flux_gnd.y)/2]
    text_label1.connect(port = 'south', destination = text_label2.ports['north'])
    text_label = D_spd.add_ref(Text_label)
    text_label.center = text_coords
        
    return D_spd