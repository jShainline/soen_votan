import copy
import numpy as np
import gdspy
import phidl, phidl.geometry as pg
from phidl import Device, Layer, LayerSet, Port
from phidl import make_device

from nc_library import wire_basic
from vt_util import vt_layers

from nc_library__vt_util import vt_arg_helper, corner_fixer, vt_label_maker
from nc_library__vt_pads_vias_wires import jj_pad, wire_tri_seg, vt_inductor, vt_m1_v3_via_sequence, via_general
from nc_library__vt_spd import vt_sy_spd
from nc_library__vt_jj import jj_with_leads, jj_circle
from nc_library__vt_sq import sq_quad
from nc_library__vt_res import res_stitch_simp, res_50_ohm

vt_lyrs,layer_data = vt_layers()


def vt_sfg(params = dict()):
    
    sfg_spd_inductance = vt_arg_helper(params,'sfg_spd_inductance',133e-9)
    spd_cntrl_inductance = vt_arg_helper(params,'spd_cntrl_inductance',1e-6)
    spd_sy_w_lead = vt_arg_helper(params,'spd_sy_w_lead',5.3)
    spd_sy_m1_outset = vt_arg_helper(params,'spd_sy_m1_outset',1.3)
    spd_jj_extra_x = vt_arg_helper(params,'spd_jj_extra_x',21.3)
    jj_inductor_extra_x = vt_arg_helper(params,'jj_inductor_extra_x',31.3)
    w_wire = vt_arg_helper(params,'w_wire',1.3)
    induct_w_wire = vt_arg_helper(params,'induct_w_wire',1.5)
    induct_num_squares_sfq1 = vt_arg_helper(params,'induct_num_squares_sfq1',234)
    induct_num_squares_sfq2 = vt_arg_helper(params,'induct_num_squares_sfq2',334)
    spd_sy_l_lead = vt_arg_helper(params,'spd_sy_l_lead',10.1)
    spd_sy_w_lead = vt_arg_helper(params,'spd_sy_w_lead',5.2)
    spd_sy_m1_olap = vt_arg_helper(params,'spd_sy_m1_olap',4.9)
    spd_sy_m1_outset = vt_arg_helper(params,'spd_sy_m1_outset',1.1)
    spd_contact_taper_length = vt_arg_helper(params,'spd_contact_taper_length',11.1)
    sfg_si_w_wire = vt_arg_helper(params,'sfg_si_w_wire',1.5)
    sfg_si_inductance = vt_arg_helper(params,'sfg_si_inductance',7.4e-9)
    sfg_si_tau = vt_arg_helper(params,'sfg_si_tau',1.34e-6)
    sfg_w_wire_sq_res_connect = vt_arg_helper(params,'sfg_w_wire_sq_res_connect',26)
    sfg_si_turn_ratio = vt_arg_helper(params,'sfg_si_turn_ratio',3)
    sfg_w_wire_tau = vt_arg_helper(params,'sfg_w_wire_tau',25)
    sfg_w_wire_flux_purge = vt_arg_helper(params,'sfg_w_wire_flux_purge',1.19)
    sfg_w_res_flux_purge = vt_arg_helper(params,'sfg_w_res_flux_purge',1.09)
    sfg_resistance_flux_purge = vt_arg_helper(params,'sfg_resistance_flux_purge',20.1)
    sfg_flux_purge_jj_offset = vt_arg_helper(params,'sfg_flux_purge_jj_offset',20.1)
    sfg_flux_purge_layer = vt_arg_helper(params,'sfg_flux_purge_layer','r2')
    sfg_include_tau_leak = vt_arg_helper(params,'sfg_include_tau_leak',True)
    sfg_has_spd = vt_arg_helper(params,'sfg_has_spd',True)
    sfg_has_sfq = vt_arg_helper(params,'sfg_has_sfq',False)
    sfg_unique_label = vt_arg_helper(params,'sfg_unique_label','dummy')
    induct_flux_purge = vt_arg_helper(params,'induct_flux_purge',100e-9)    
    inductance_per_sq_stf = vt_arg_helper(params,'inductance_per_sq_stf',200e-12)
    inductance_per_sq_m1 = vt_arg_helper(params,'inductance_per_sq_m1',200e-15)
    inductance_per_sq_m3 = vt_arg_helper(params,'inductance_per_sq_m3',200e-15)
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
    
    D_sfg = Device('sfg')
    
    #spd
    if sfg_has_spd == True:
        params_mod = copy.deepcopy(params)
        params_mod['spd_inductance'] = sfg_spd_inductance
        params_mod['fiber_collar_orientation'] = 0
        spd = D_sfg.add_ref(vt_sy_spd(params_mod))
    elif sfg_has_spd == False:
        cf_init = D_sfg.add_ref(corner_fixer(w_wire = w_wire, layer = layers['jj1']))
    
    #jj_sfq
    if sfg_has_sfq == True:
        jj_sfq = D_sfg.add_ref(jj_circle(params))
        if sfg_has_spd == True:
            jj_sfq.connect(port = 'jj1_north', destination = spd.ports['wire_port'])
            jj_sfq.movex(w_wire)
            D_sfg.add_ref(wire_basic(jj_sfq.ports['jj1_north'],spd.ports['wire_port'],'x',w_wire,layers['jj1']))
        elif sfg_has_spd == False:
            jj_sfq.connect(port = 'jj1_north', destination = cf_init.ports['east'])
            jj_sfq.movex(w_wire)
            D_sfg.add_ref(wire_basic(jj_sfq.ports['jj1_north'],cf_init.ports['east'],'x',w_wire,layers['jj1']))
        params_mod = copy.deepcopy(params)
        params_mod['induct_num_squares'] = induct_num_squares_sfq1
        params_mod['induct_layer'] = 'jj1'
        ind_sfq1 = D_sfg.add_ref(vt_inductor(params_mod))
        ind_sfq1.connect(port = 'west', destination = jj_sfq.ports['jj1_west'])        
    
    #jj1
    jj1 = D_sfg.add_ref(jj_with_leads(params))
    if sfg_has_sfq == False:
        if sfg_has_spd == True:
            jj1.connect(port = 'xcntct_west', destination = spd.ports['wire_port'])
            jj1.movex(spd_jj_extra_x)
            D_sfg.add_ref(wire_basic(spd.ports['wire_port'],jj1.ports['xcntct_west'],'x',w_wire,layers['m3']))
        elif sfg_has_spd == False:
            jj1.connect(port = 'xcntct_west', destination = cf_init.ports['east'])
            jj1.movex(spd_jj_extra_x)
            D_sfg.add_ref(wire_basic(cf_init.ports['east'],jj1.ports['xcntct_west'],'x',w_wire,layers['m3']))
    if sfg_has_sfq == True:
        jj1.connect(port = 'xcntct_west', destination = jj_sfq.ports['m3_south'])
        jj1.movex(spd_jj_extra_x)
        D_sfg.add_ref(wire_basic(jj_sfq.ports['m3_south'],jj1.ports['xcntct_west'],'x',w_wire,layers['m3']))
        tpr_sfq = D_sfg.add_ref(pg.taper(length = 3*w_wire, width1 = w_wire, width2 = induct_w_wire, port = None, layer = layers['jj1']))
        tpr_sfq.connect(port = 1, destination = jj1.ports['jj1_south'])
        params_mod['induct_num_squares'] = induct_num_squares_sfq2
        params_mod['induct_layer'] = 'jj1'
        ind_sfq2 = D_sfg.add_ref(vt_inductor(params_mod))
        ind_sfq2.connect(port = 'west', destination = tpr_sfq.ports[2])
    Tpr = pg.taper(length = jj_inductor_extra_x, width1 = w_wire, width2 = induct_w_wire, port = None, layer = layers['m3'])
    tpr1 = D_sfg.add_ref(Tpr)
    tpr1.connect(port = 1, destination = jj1.ports['xcntct_east'])
    
    #inductor 1
    params_mod = copy.deepcopy(params)
    params_mod['induct_num_squares'] = induct_num_squares_sfq1
    params_mod['induct_layer'] = 'm3'
    Ind = vt_inductor(params_mod)
    ind1 = D_sfg.add_ref(Ind)
    ind1.connect(port = 'west', destination = tpr1.ports[2])
    
    #jj2
    tpr2 = D_sfg.add_ref(Tpr)
    tpr2.connect(port = 2, destination = ind1.ports['east'])
    jj2 = D_sfg.add_ref(jj_with_leads(params))
    jj2.connect(port = 'xcntct_west', destination = tpr2.ports[1])
    tpr3 = D_sfg.add_ref(Tpr)
    tpr3.connect(port = 1, destination = jj2.ports['xcntct_east'])
    
    #inductor 2
    ind2 = D_sfg.add_ref(Ind)
    ind2.connect(port = 'west', destination = tpr3.ports[2])
    
    #jj3 
    tpr4 = D_sfg.add_ref(Tpr)
    tpr4.connect(port = 2, destination = ind2.ports['east'])
    jj3 = D_sfg.add_ref(jj_with_leads(params))
    jj3.connect(port = 'xcntct_west', destination = tpr4.ports[1])
    tpr5 = D_sfg.add_ref(Tpr)
    tpr5.connect(port = 1, destination = jj3.ports['xcntct_east'])
    
    #vias in si loop
    Cntct = Device('L_si_cntct_ind')
    vias1 = Cntct.add_ref(vt_m1_v3_via_sequence(params))
    vias1.connect(port = 'm3_east', destination = tpr5.ports[2])
    Cntct.add_port(name = 'via', midpoint = vias1.ports['m3_east'].midpoint, width = w_wire, orientation = 180)
    
    #L_si determine layer and num_squares
    sfg_induct_num_squares_stf = np.around(sfg_si_inductance/inductance_per_sq_stf,decimals = 2)
    if sfg_induct_num_squares_stf < 40:        
        inductor_layer = 'm1'
        sfg_induct_num_squares = np.around(sfg_si_inductance/inductance_per_sq_m1, decimals = 2)
    elif sfg_induct_num_squares_stf >= 40:
        sfg_induct_num_squares = sfg_induct_num_squares_stf
        inductor_layer = 'stf'
        
    #L_si contact
    if inductor_layer == 'stf':
        tpr6 = Cntct.add_ref(pg.taper(length = spd_contact_taper_length, width1 = w_wire, width2 = spd_sy_w_lead+2*spd_sy_m1_outset, port = None, layer = layers['m1']))
        tpr6.connect(port = 1, destination = vias1.ports['m1_west'])
        rect_m1 = Cntct.add_ref(pg.rectangle(size = [spd_sy_m1_olap,spd_sy_w_lead+2*spd_sy_m1_outset], layer = layers['m1']))
        rect_m1.center = [tpr6.ports[2].midpoint[0]+spd_sy_m1_olap/2,tpr6.ports[2].midpoint[1]]
        rect_stf = Cntct.add_ref(pg.rectangle(size = [spd_sy_m1_olap,spd_sy_w_lead], layer = layers['stf']))
        rect_stf.center = [rect_m1.center[0],tpr6.ports[2].midpoint[1]]
        tpr_stf = Cntct.add_ref(pg.taper(length = spd_sy_l_lead, width1 = spd_sy_w_lead, width2 = sfg_si_w_wire, port = None, layer = layers['stf']))
        tpr_stf.connect(port = 1, destination = vias1.ports['m1_west'])
        tpr_stf.movex(spd_contact_taper_length+spd_sy_m1_olap)
        Cntct.add_port(name = 'stf', midpoint = tpr_stf.ports[2].midpoint, width = sfg_si_w_wire, orientation = 0)
        cntct1 = D_sfg.add_ref(Cntct)
        cntct1.connect(port = 'via', destination = tpr5.ports[2])
    elif inductor_layer == 'm1':
        tpr6 = Cntct.add_ref(pg.taper(length = spd_contact_taper_length, width1 = w_wire, width2 = sfg_si_w_wire, port = None, layer = layers['m1']))
        tpr6.connect(port = 1, destination = vias1.ports['m1_west'])
        Cntct.add_port(name = 'stf', midpoint = tpr6.ports[2].midpoint, width = sfg_si_w_wire, orientation = 0)
        cntct1 = D_sfg.add_ref(Cntct)
        cntct1.connect(port = 'via', destination = tpr5.ports[2])
        
    #L_si
    if inductor_layer == 'stf':
        induct_si = D_sfg.add_ref(pg.snspd(wire_width = sfg_si_w_wire, wire_pitch = 2*sfg_si_w_wire, size = None, num_squares = sfg_induct_num_squares, turn_ratio = sfg_si_turn_ratio, terminals_same_side = False, layer = layers['stf']))
        induct_si.connect(port = 1, destination = cntct1.ports['stf'])
        cntct2 = D_sfg.add_ref(Cntct)
        cntct2.connect(port = 'stf', destination = induct_si.ports[2])
    elif inductor_layer == 'm1':
        induct_si = D_sfg.add_ref(pg.snspd(wire_width = sfg_si_w_wire, wire_pitch = 2*sfg_si_w_wire, size = None, num_squares = sfg_induct_num_squares, turn_ratio = sfg_si_turn_ratio, terminals_same_side = False, layer = layers['m1']))
        induct_si.connect(port = 1, destination = cntct1.ports['stf'])
        cntct2 = D_sfg.add_ref(Cntct)
        cntct2.connect(port = 'stf', destination = induct_si.ports[2])
    
    #squid
    if inductor_layer == 'stf':
        sq = D_sfg.add_ref(sq_quad(params))
        sq.connect(port = 'incoil_2', destination = cntct2.ports['via'])
    if inductor_layer == 'm1':
        sq = D_sfg.add_ref(sq_quad(params))
        sq.connect(port = 'incoil_2', destination = cntct2.ports['via'])
    params_mod = copy.deepcopy(params)
    params_mod['via_layer'] = 'v2'
    params_mod['layer_below'] = 'jj1'
    params_mod['layer_above'] = 'm3'
    sq_via1 = D_sfg.add_ref(via_general(params_mod))
    sq_via2 = D_sfg.add_ref(via_general(params_mod))
    sq_via1.connect(port = 'above_west', destination = sq.ports['bias_1'])
    sq_via2.connect(port = 'below_south', destination = sq_via1.ports['below_north'])
    sq_via2.center = [sq_via2.center[0],induct_si.ymax]
    D_sfg.add_ref(wire_basic(sq_via1.ports['below_north'],sq_via2.ports['below_south'],'y',w_wire,layers['jj1']))
    res_50_1 = D_sfg.add_ref(res_50_ohm(params))
    res_50_1.connect(port = 'west', destination = sq_via2.ports['above_north'])
    res_50_2 = D_sfg.add_ref(res_50_ohm(params))
    res_50_2.connect(port = 'west', destination = sq_via2.ports['above_east'])
    
    #r_si
    if sfg_include_tau_leak == True:
        vias2 = D_sfg.add_ref(vt_m1_v3_via_sequence(params))
        vias2.connect(port = 'm3_east', destination  = sq.ports['incoil_1'])
        params_mod = copy.deepcopy(params)
        params_mod['res_w_wire'] = sfg_w_wire_tau
        params_mod['res_resistance'] = sfg_si_inductance/sfg_si_tau
        params_mod['res_layer'] = 'r1'
        params_mod['w_wire'] = sfg_w_wire_sq_res_connect
        params_mod['res_include_label'] = True
        res = D_sfg.add_ref(res_stitch_simp(params_mod))
        res.center = [jj3.xmax+(induct_si.xmin-jj3.xmax)/2,sq.ports['bias_1'].midpoint[1]]
        D_sfg.add_ref(wire_tri_seg(p1 = res.ports['east'], p2 = vias2.ports['m1_west'], initial_width = sfg_w_wire_sq_res_connect, final_width = sq_incoil_w_wire, length_factors = [0.1,2,0.1], directions = 'xx', layer = layers['m1']))
    
    #flux_purge
    D_fp = Device('sfg_flux_purge')
    params_mod = copy.deepcopy(params)
    params_mod['res_w_wire'] = sfg_w_res_flux_purge
    params_mod['res_resistance'] = sfg_resistance_flux_purge
    params_mod['res_layer'] = sfg_flux_purge_layer
    params_mod['w_wire'] = sfg_w_wire_flux_purge
    Res_fp = res_stitch_simp(params_mod)
#    x_coords = np.array([jj1.x,jj2.x,jj3.x,induct_si.x,sq.x])
    x_coords = np.array([jj1.x,jj2.x,jj3.x,induct_si.x])
    temp_y = sfg_flux_purge_jj_offset+sfg_w_wire_flux_purge/2
#    y_coords = np.array([jj1.ymin-temp_y,jj2.ymin-temp_y,jj3.ymin-temp_y,induct_si.y,sq.ymax+4*sfg_w_wire_flux_purge])
    y_coords = np.array([jj1.ymin-temp_y,jj2.ymin-temp_y,jj3.ymin-temp_y,induct_si.y])
    res_fp1 = D_fp.add_ref(Res_fp)
    res_fp1.center = [x_coords[0],y_coords[0]]
    res_fp2 = D_fp.add_ref(Res_fp)
    res_fp2.center = [x_coords[1],y_coords[1]]
    res_fp3 = D_fp.add_ref(Res_fp)
    res_fp3.center = [x_coords[2],y_coords[2]]
    res_fp4 = D_fp.add_ref(Res_fp)
    res_fp4.rotate(90)
    res_fp4.center = [x_coords[3],y_coords[3]] 
#    res_fp5 = D_fp.add_ref(Res_fp)
#    res_fp5.center = [x_coords[4],y_coords[4]] 
    Cf_fp = corner_fixer(sfg_w_wire_flux_purge,layers['m3'])
    cf_fp1 = D_fp.add_ref(Cf_fp)
    cf_fp1.connect(port = 'north', destination = res_fp4.ports['west'])
    cf_fp2 = D_fp.add_ref(Cf_fp)
    cf_fp2.connect(port = 'south', destination = res_fp4.ports['east'])    
    #wires
    D_fp.add_ref(wire_basic(res_fp1.ports['east'],res_fp2.ports['west'],'xyxxxxxx',sfg_w_wire_flux_purge,layers['m3']))
    D_fp.add_ref(wire_basic(res_fp2.ports['east'],res_fp3.ports['west'],'xyxxxxxx',sfg_w_wire_flux_purge,layers['m3']))
    D_fp.add_ref(wire_basic(res_fp3.ports['east'],cf_fp1.ports['west'],'xxxxxyx',sfg_w_wire_flux_purge,layers['m3']))
#    D_fp.add_ref(wire_basic(cf_fp2.ports['east'],res_fp5.ports['west'],'xyx',sfg_w_wire_flux_purge,layers['m3']))
 
    
    D_sfg.add_ref(D_fp)
    
    #####
    #pads
    #####
    Pad = jj_pad(params)
    params_mod = copy.deepcopy(params)
    params_mod['is_ground_pad'] = True
    Pad_gnd = jj_pad(params_mod)
    
    #spd pads
    Cf = corner_fixer(w_wire,layers['m3'])
    pad_spd = D_sfg.add_ref(Pad)
    if sfg_has_spd == True:
        initial_x_offset = -500
        initial_y_offset = pad_spd.ysize/2+spd.ysize/2#+20
        cf1 = D_sfg.add_ref(Cf)
        cf1.connect(port = 'south', destination = spd.ports['pad'])
        pad_spd.connect(port = 'pad_anchor', destination = spd.ports['fiber_center'])
        pad_spd.move([initial_x_offset,initial_y_offset])#+pad_y_backset
        D_sfg.add_ref(wire_tri_seg(p1 = pad_spd.ports['m1_south'], p2 = cf1.ports['west'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.2,0.1], directions = 'yx', layer = layers['m3']))
        pad_spd_gnd = D_sfg.add_ref(Pad_gnd)
        pad_spd_gnd.connect(port = 'm1_north', destination = spd.ports['pad_gnd'])        
    elif sfg_has_spd == False:
        initial_x_offset = -500
        initial_y_offset = pad_spd.ysize/2
        pad_spd.connect(port = 'jj1_south', destination = cf_init.ports['north'])
        pad_spd.move([initial_x_offset,initial_y_offset])#+pad_y_backset
        D_sfg.add_ref(wire_tri_seg(p1 = pad_spd.ports['jj1_south'], p2 = cf_init.ports['west'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.2,0.1], directions = 'yx', layer = layers['jj1']))

    #jj_sfq ground pad
    if sfg_has_sfq == True:
        pad_jjsfq_gnd = D_sfg.add_ref(Pad_gnd) 
        pad_jjsfq_gnd.connect(port = 'jj1_north', destination = ind_sfq1.ports['east'])
        pad_jjsfq_gnd.movey(-pad_gnd_y_backset)
        D_sfg.add_ref(wire_tri_seg(p1 = pad_jjsfq_gnd.ports['m3_north'], p2 = ind_sfq1.ports['east'], initial_width = pad_w_wire, final_width = induct_w_wire, length_factors = [4,4,0.5], directions = 'yy', layer = layers['jj1']))
        
    #jj1 pads
    pad_jj1 = D_sfg.add_ref(Pad)
    pad_jj1_gnd = D_sfg.add_ref(Pad_gnd) 
    pad_jj1.center = [pad_spd.center[0]+pad_pitch[0],pad_spd.center[1]]
    D_sfg.add_ref(wire_tri_seg(p1 = pad_jj1.ports['m3_south'], p2 = jj1.ports['xcntct_north'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [1,4,3], directions = 'yy', layer = layers['m3']))
    if sfg_has_sfq == False:
        pad_jj1_gnd.connect(port = 'jj1_north', destination = jj1.ports['jj1_south'])
        pad_jj1_gnd.movey(-pad_gnd_y_backset)
        D_sfg.add_ref(wire_tri_seg(p1 = pad_jj1_gnd.ports['jj1_north'], p2 = jj1.ports['jj1_south'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [4,4,2], directions = 'yy', layer = layers['jj1']))
    if sfg_has_sfq == True:
        pad_jj1_gnd.connect(port = 'jj1_north', destination = ind_sfq2.ports['east'])
        pad_jj1_gnd.movey(-pad_gnd_y_backset)
        D_sfg.add_ref(wire_tri_seg(p1 = pad_jj1_gnd.ports['jj1_north'], p2 = ind_sfq2.ports['east'], initial_width = pad_w_wire, final_width = induct_w_wire, length_factors = [4,4,0.5], directions = 'yy', layer = layers['jj1']))

    #jj2 pads
    pad_jj2 = D_sfg.add_ref(Pad)
    pad_jj2_gnd = D_sfg.add_ref(Pad_gnd) 
    pad_jj2.center = [pad_jj1.center[0]+pad_pitch[0],pad_spd.center[1]]
    D_sfg.add_ref(wire_tri_seg(p1 = pad_jj2.ports['m3_south'], p2 = jj2.ports['xcntct_north'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [1.4,4,1], directions = 'yy', layer = layers['m3']))
    pad_jj2_gnd.connect(port = 'm1_north', destination = jj2.ports['jj1_south'])
    pad_jj2_gnd.movey(-pad_gnd_y_backset)
    D_sfg.add_ref(wire_tri_seg(p1 = pad_jj2_gnd.ports['jj1_north'], p2 = jj2.ports['jj1_south'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [4,4,1], directions = 'yy', layer = layers['jj1']))    

    #jj3 pads
    cf2 = D_sfg.add_ref(Cf)
    cf2.connect(port = 'south', destination = jj3.ports['xcntct_north'])
    pad_jj3 = D_sfg.add_ref(Pad)
    pad_jj3_gnd = D_sfg.add_ref(Pad_gnd) 
    pad_jj3.center = [pad_jj2.center[0]+pad_pitch[0],pad_spd.center[1]]
    D_sfg.add_ref(wire_tri_seg(p1 = pad_jj3.ports['m3_south'], p2 = cf2.ports['east'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.4,0.2], directions = 'yx', layer = layers['m3']))
    pad_jj3_gnd.connect(port = 'm1_north', destination = jj3.ports['jj1_south'])
    pad_jj3_gnd.movey(-pad_gnd_y_backset)
    D_sfg.add_ref(wire_tri_seg(p1 = pad_jj3_gnd.ports['jj1_north'], p2 = jj3.ports['jj1_south'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [4,4,1], directions = 'yy', layer = layers['jj1']))    
    
    #si loop ground
    pad_si_gnd = D_sfg.add_ref(Pad_gnd)
    if sfg_include_tau_leak == True:
        pad_si_gnd.connect(port = 'm1_east', destination = res.ports['west'])
        pad_si_gnd.movex(-pad_gnd_y_backset)
        D_sfg.add_ref(wire_tri_seg(p1 = pad_si_gnd.ports['m1_east'], p2 = res.ports['west'], initial_width = pad_w_wire, final_width = sfg_w_wire_sq_res_connect, length_factors = [1,2,1], directions = 'xx', layer = layers['m1']))    
    elif sfg_include_tau_leak == False:
        pad_si_gnd.connect(port = 'm1_east', destination = sq.ports['incoil_1'])
        pad_si_gnd.move([-pad_size_ground[0]/2-pad_gnd_y_backset,-pad_size_ground[1]/2-pad_gnd_y_backset])
        D_sfg.add_ref(wire_tri_seg(p1 = pad_si_gnd.ports['m3_east'], p2 = sq.ports['incoil_1'], initial_width = pad_w_wire, final_width = sq_incoil_w_wire, length_factors = [1,2,1], directions = 'xx', layer = layers['m3']))    
    
    #sq
    pad_sq_I = D_sfg.add_ref(Pad)
    pad_sq_I.center = [pad_jj3.center[0]+pad_pitch[0],pad_spd.center[1]]
    D_sfg.add_ref(wire_tri_seg(p1 = pad_sq_I.ports['m3_south'], p2 = res_50_1.ports['east'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.2,0.5,0.2], directions = 'yy', layer = layers['m3']))
    pad_sq_V = D_sfg.add_ref(Pad)
    pad_sq_V.center = [pad_sq_I.center[0]+pad_pitch[0],pad_spd.center[1]]
    D_sfg.add_ref(wire_tri_seg(p1 = pad_sq_V.ports['m3_south'], p2 = res_50_2.ports['east'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.2,0.5], directions = 'yx', layer = layers['m3']))    
    pad_sq_gnd = D_sfg.add_ref(Pad_gnd)
    pad_sq_gnd.center = sq.ports['bias_2'].midpoint
    pad_sq_gnd.movex(pad_size_ground[0]/2+pad_gnd_y_backset)
    D_sfg.add_ref(wire_tri_seg(p1 = pad_sq_gnd.ports['jj1_west'], p2 = sq.ports['bias_2'], initial_width = pad_w_wire, final_width = w_wire, length_factors = [0.4,0.6,0.1], directions = 'xx', layer = layers['jj1']))    
    
    #addflux
    pad_flux = D_sfg.add_ref(Pad)
    pad_flux.center = [pad_sq_V.center[0]+pad_pitch[0],pad_spd.center[1]]
    D_sfg.add_ref(wire_tri_seg(p1 = pad_flux.ports['m3_south'], p2 = sq.ports['addflux_east'], initial_width = pad_w_wire, final_width = sq_addflux_w_wire, length_factors = [0.2,0.5], directions = 'yx', layer = layers['m3']))    
    pad_flux_gnd = D_sfg.add_ref(Pad_gnd)
    pad_flux_gnd.center = [induct_si.x,sq.ports['addflux_west'].midpoint[1]+20]
    D_sfg.add_ref(wire_tri_seg(p1 = pad_flux_gnd.ports['m3_east'], p2 = sq.ports['addflux_west'], initial_width = pad_w_wire, final_width = sq_addflux_w_wire, length_factors = [0.4,0.6,0.1], directions = 'xx', layer = layers['m3']))    
    
    #control spd
    if sfg_has_spd == True:
        params_mod = copy.deepcopy(params)
        params_mod['spd_sy_include_res'] = False
        params_mod['spd_inductance'] = spd_cntrl_inductance
        params_mod['fiber_collar_orientation'] = 180
        spd_cntrl = D_sfg.add_ref(vt_sy_spd(params_mod))
        pad_spd_cntrl = D_sfg.add_ref(Pad)
        pad_spd_cntrl.center = [pad_flux.center[0]+pad_pitch[0],pad_spd.center[1]]
        spd_cntrl.center = pad_spd_cntrl.center
        spd_cntrl.move([spd_cntrl.center[0]-spd_cntrl.ports['fiber_center'].midpoint[0],-pad_spd_cntrl.ysize/2-spd_cntrl.ysize/2])
        D_sfg.add_ref(wire_tri_seg(p1 = pad_spd_cntrl.ports['m3_south'], p2 = spd_cntrl.ports['pad'], initial_width = pad_w_wire, final_width = spd_cntrl.ports['pad'].width, length_factors = [2,5,1], directions = 'yy', layer = layers['m1']))    
        pad_spd_cntrl_gnd = D_sfg.add_ref(Pad_gnd)
        pad_spd_cntrl_gnd.connect(port = 'm1_north', destination = spd_cntrl.ports['pad_gnd'])
    
    #flux_purge
    pad_fp = D_sfg.add_ref(Pad)
    pad_fp.center = [pad_spd.center[0]-pad_pitch[0],pad_spd.center[1]]
    pad_fp_gnd = D_sfg.add_ref(Pad_gnd)
    spd_fp_outset = 4*sfg_w_wire_flux_purge
    pad_fp_gnd.center = [sq.ports['bias_1'].midpoint[0]-pad_fp_gnd.xsize,induct_si.ymax+pad_fp_gnd.ysize]
    
    #flux_purge_wires
    p1 = res_fp1.ports['west'].midpoint
    if sfg_has_sfq == False:
        x_east = min([spd.ports['max_coords_without_collar'].midpoint[0]+4*sq_addflux_w_wire,pad_jj1_gnd.xmin-4*sq_addflux_w_wire])
    elif sfg_has_sfq == True: 
        x_east = pad_jjsfq_gnd.xmax+(pad_jj1_gnd.xmin-pad_jjsfq_gnd.xmax)/2
    y_north = res_fp1.ports['west'].midpoint[1] 
    #flux_purge extra inductor       
    induct_num_squares = induct_flux_purge/inductance_per_sq_m3
    params_mod = copy.deepcopy(params)
    params_mod['induct_num_squares'] = induct_num_squares
    params_mod['induct_layer'] = 'm3'
    params_mod['induct_w_wire'] = sfg_w_wire_flux_purge
    params_mod['induct_wire_pitch'] = 2*sfg_w_wire_flux_purge
    ind_fp = D_sfg.add_ref(vt_inductor(params_mod))    
    if sfg_has_spd == True:
        x_west = spd.ports['min_coords_without_collar'].midpoint[0]-spd_fp_outset
        y_south = pad_spd_gnd.ymin-spd_fp_outset
        p2 = [x_east,y_north]
        p3 = [x_east,y_south]
        p4 = [x_west,y_south]    
        points = [p1,p2,p3,p4]        
        route_path = gdspy.PolyPath(points, width = sfg_w_wire_flux_purge)
        Wire = Device('sffg_flux_purge_input_wire')
        Wire.add_polygon(route_path.polygons, layer = layers['m3'])
        Wire.add_port(name = 'flux_purge_input', midpoint = points[-1], width = sfg_w_wire_flux_purge, orientation = 180)
        wire = D_sfg.add_ref(Wire)        
        ind_fp.center = [pad_spd.x,spd.y]
        D_sfg.add_ref(wire_basic(ind_fp.ports['east'],wire.ports['flux_purge_input'],'xyx',sfg_w_wire_flux_purge,layers['m3']))
        D_sfg.add_ref(wire_tri_seg(p1 = pad_fp.ports['m3_south'], p2 = ind_fp.ports['west'], initial_width = pad_w_wire, final_width = sfg_w_wire_flux_purge, length_factors = [0.3,0.01], directions = 'yx', layer = layers['m3']))    
    elif sfg_has_spd == False:
        ind_fp.center = [pad_spd.x,ind_sfq1.y]
        D_sfg.add_ref(wire_basic(ind_fp.ports['east'],res_fp1.ports['west'],'xyx',sfg_w_wire_flux_purge,layers['m3']))
        D_sfg.add_ref(wire_tri_seg(p1 = pad_fp.ports['m3_south'], p2 = ind_fp.ports['west'], initial_width = pad_w_wire, final_width = sfg_w_wire_flux_purge, length_factors = [0.3,0.01], directions = 'yx', layer = layers['m3']))    
#    D_sfg.add_ref(wire_tri_seg(p1 = pad_fp_gnd.ports['m3_south'], p2 = res_fp5.ports['east'], initial_width = pad_w_wire, final_width = sfg_w_wire_flux_purge, length_factors = [0.3,0.2], directions = 'yx', layer = layers['m3']))    
    D_sfg.add_ref(wire_tri_seg(p1 = pad_fp_gnd.ports['m3_south'], p2 = cf_fp2.ports['east'], initial_width = pad_w_wire, final_width = sfg_w_wire_flux_purge, length_factors = [0.3,0.2], directions = 'yx', layer = layers['m3']))    
        
    #pad anchor
    D_sfg.add_port(name = 'pad_anchor', midpoint = [pad_spd.ports['m3_north'].midpoint[0],pad_spd.ports['m3_west'].midpoint[1]], width = pad_w_wire, orientation = 270)
    
    #label
    Text_label = Device('sfg_text_label')
    text_label1 = Text_label.add_ref(vt_label_maker(text_string = 'SFG: '+sfg_unique_label, text_size = 1.2*device_label_size, justify = 'center', layer = layers[label_layer]))
    if sfg_include_tau_leak == True:
        text_string = 'si_inductance = '+str(np.around(sfg_si_inductance*1e9,decimals = 2))+' nH\nnum_sq_induct = '+str(sfg_induct_num_squares)+'\ntau_si = '+str(np.around(sfg_si_tau*1e6, decimals = 2))+' us\nr_si = '+str(np.around(sfg_si_inductance/sfg_si_tau, decimals = 3))+' ohm'
    if sfg_include_tau_leak == False:    
        text_string = 'si_inductance = '+str(np.around(sfg_si_inductance*1e9,decimals = 2))+' nH\nnum_sq_induct = '+str(sfg_induct_num_squares)
    text_label2 = Text_label.add_ref(vt_label_maker(text_string = text_string, text_size = device_label_size, justify = 'center', layer = layers[label_layer]))
    text_coords = [jj1.x+(jj2.x-jj1.x)/2,pad_si_gnd.y]
    text_label1.connect(port = 'south', destination = text_label2.ports['north'])
    text_label = D_sfg.add_ref(Text_label)
    text_label.center = text_coords
        
    return D_sfg