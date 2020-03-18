from __future__ import division, print_function, absolute_import
import numpy as np

import phidl
import phidl.geometry as pg
from phidl import make_device, Device, Layer, LayerSet, quickplot2 as qp
from nc_devices import *

from vt_util import *
from vt_params import vt_parameters
from vt_post_processing import *

import copy
import time

#%%
print('\n\nbuilding gds...')
t_tot = time.time()

D_jj_a = Device('votan1__jj_a')#jj chip with 40uA Ic
D_jj_b = Device('votan1__jj_b')#jj chip with 100uA Ic
D_sq_a = Device('votan1__sq_a')#squid chip with 40uA Ic
D_sq_b = Device('votan1__sq_b')#squid chip with 100uA Ic

vt_lyrs,layer_data = vt_layers()
vt_lyrs_post,layer_data_post = vt_layers_post()
write_lyp(vt_lyrs,layer_data,'vt')
write_lyp(vt_lyrs_post,layer_data_post,'vt_post')
params  = vt_parameters(vt_lyrs)
params_mod = copy.deepcopy(params)

#%% geometry setup

x_inset = 600#how far lower left pad is inset from chip edge
y_inset = 25

P_loc = pad_locations(**params)
p_loc_jj_a = D_jj_a.add_ref(P_loc)
p_loc_jj_b = D_jj_b.add_ref(P_loc)
p_loc_sq_a = D_sq_a.add_ref(P_loc)
p_loc_sq_b = D_sq_b.add_ref(P_loc)
p_loc_jj_a.center = [0,0] 
p_loc_jj_b.center = [0,0]
p_loc_sq_a.center = [0,0] 
p_loc_sq_b.center = [0,0]

num_jjs = len(params['jj_junc_diam_vec_040'])

#%% which to run

which_to_run = dict(run_all = True,
                    jj_single = False,
                    jj_series = False,
                    jj_series_dense = False,
                    sq_single = False,
                    via_series = False,
                    inline_test = False,
                    resistance_test = False,
                    isolation_test = False,
                    litho_test = False,
                    endpoint_boxes = False,
                    extra_pads = False,
                    labels = False,
                    post_processing_jj_a = True,
                    post_processing_jj_b = True,
                    post_processing_sq_a = True,
                    post_processing_sq_b = True)

chip_to_run = dict(jj_a = True,
                   jj_b = True,
                   sq_a = True,
                   sq_b = True)

#%% single josephson junctions

if which_to_run['run_all'] == True or which_to_run['jj_single'] == True:

    port_list = ['p_s_a_1','p_s_a_3','p_s_a_5','p_s_a_7','p_e_a_1','p_e_a_3','p_e_a_5','p_e_a_7']
    
    if chip_to_run['jj_a'] == True:
        D_jj = Device('jj_four_wires')
        params_mod = copy.deepcopy(params)
        for pp in range(len(params['shunt_num_squares_vec_040'])):
            params_mod['shunt_num_squares'] = params['shunt_num_squares_vec_040'][pp]
            for ii in range(len(params['jj_junc_diam_vec_040'])):
                params_mod['jj_junc_diam'] = params['jj_junc_diam_vec_040'][ii]
                jj = jj_4wire(params_mod)
                jj_4wire_instance = D_jj.add_ref(jj)
                jj_4wire_instance.connect(port = 'pad_anchor', destination = p_loc_jj_a.ports[port_list[pp*num_jjs+ii]])
            
        D_jj_a.add_ref(D_jj)    

    if chip_to_run['jj_b'] == True:
        D_jj = Device('jj_four_wires')
        params_mod = copy.deepcopy(params)
        for pp in range(len(params['shunt_num_squares_vec_100'])):
            params_mod['shunt_num_squares'] = params['shunt_num_squares_vec_100'][pp]
            for ii in range(len(params['jj_junc_diam_vec_100'])):
                params_mod['jj_junc_diam'] = params['jj_junc_diam_vec_100'][ii]
                jj = jj_4wire(params_mod)
                jj_4wire_instance = D_jj.add_ref(jj)
                jj_4wire_instance.connect(port = 'pad_anchor', destination = p_loc_jj_b.ports[port_list[pp*num_jjs+ii]])
            
        D_jj_b.add_ref(D_jj) 
        
#%% jj series arrays
    
if which_to_run['run_all'] == True or which_to_run['jj_series'] == True:    

    port_list = ['p_n_a_1','p_n_a_2','p_n_a_3','p_n_a_4','p_n_a_5','p_n_a_6','p_n_a_7','p_n_a_8','p_w_a_1','p_w_a_2','p_w_a_3','p_w_a_4','p_w_a_5','p_w_a_6','p_w_a_7','p_w_a_8']
    if chip_to_run['jj_a'] == True:
                
        num_groups = len(params['num_cols_jjs_vec'])
        D_jjs = Device('jj_series_arrays')
        params_mod = copy.deepcopy(params)
        params_mod['jj_include_flux_moats'] = False
        for pp in range(len(params['shunt_num_squares_vec_040'])):
            params_mod['shunt_num_squares'] = params['shunt_num_squares_vec_040'][pp]
            for ii in range(num_jjs):
                params_mod['jj_junc_diam'] = params['jj_junc_diam_vec_040'][ii]
                for kk in range(num_groups):
                    params_mod['num_cols_jjs'] = params['num_cols_jjs_vec'][kk]
                    jjs = jj_series_array(params_mod)
                    jj_series_array_instance = D_jjs.add_ref(jjs)
                    jj_series_array_instance.connect(port = 'pad_anchor', destination = p_loc_jj_a.ports[port_list[pp*num_jjs*num_groups+ii*num_groups+kk]])
                
        D_jj_a.add_ref(D_jjs)
        
    if chip_to_run['jj_b'] == True:
                
        num_groups = len(params['num_cols_jjs_vec'])
        D_jjs = Device('jj_series_arrays')
        params_mod = copy.deepcopy(params)
        params_mod['jj_include_flux_moats'] = False
        for pp in range(len(params['shunt_num_squares_vec_100'])):
            params_mod['shunt_num_squares'] = params['shunt_num_squares_vec_100'][pp]
            for ii in range(num_jjs):
                params_mod['jj_junc_diam'] = params['jj_junc_diam_vec_100'][ii]
                for kk in range(num_groups):
                    params_mod['num_cols_jjs'] = params['num_cols_jjs_vec'][kk]
                    jjs = jj_series_array(params_mod)
                    jj_series_array_instance = D_jjs.add_ref(jjs)
                    jj_series_array_instance.connect(port = 'pad_anchor', destination = p_loc_jj_b.ports[port_list[pp*num_jjs*num_groups+ii*num_groups+kk]])
                
        D_jj_b.add_ref(D_jjs)        
  
#%% jj series arrays dense

if which_to_run['run_all'] == True or which_to_run['jj_series_dense'] == True:

    if chip_to_run['jj_a'] == True:
        
        num_jjs_dense = len(params['jj_junc_diam_vec_dense_040'])
        D_jjs_dense = Device('jj_series_arrays_dense')
        params_mod = copy.deepcopy(params)
        params_mod['jj_include_shunt'] = False
        params_mod['pad_size'] = params['inline_test_pad_size']
        params_mod['pad_pitch'] = [2*params_mod['pad_size'][0],params_mod['pad_size'][1]]
        params_mod['jj_include_flux_moats'] = False
        num_cols_jjs_vec = params['num_cols_jjs_vec_dense']
        for kk in range(len(num_cols_jjs_vec)):
            D_jjs_dense_row = Device('jj_series_arrays_dense_row')
            params_mod['num_cols_jjs'] = num_cols_jjs_vec[kk]
            for ii in range(num_jjs_dense):
                params_mod['jj_junc_diam'] = params['jj_junc_diam_vec_dense_040'][ii]
                jjs = jj_series_array(params_mod)
                jj_series_array_instance = D_jjs_dense_row.add_ref(jjs) 
                jj_series_array_instance.movex(4*ii*params_mod['pad_size'][0]+60)
            row_instance = D_jjs_dense.add_ref(D_jjs_dense_row)
            row_instance.movey(kk*460)
         
        D_jjs_dense_text = Device('jj_series_dense_text')
        Text = pg.text(text = 'jj series arrays', size = params['block_label_size'], layer = vt_lyrs[params['label_layer']])
        text = D_jjs_dense_text.add_ref(Text)
        text.rotate(90)
        text.center = D_jjs_dense.center
        text.move([-D_jjs_dense.xsize/2-params['block_label_size'],0])
        D_jjs_dense.add_ref(D_jjs_dense_text)
               
        jjs_dense = D_jj_a.add_ref(D_jjs_dense)            
        jjs_dense.center = [40,-970]   

    if chip_to_run['jj_b'] == True:
        
        num_jjs_dense = len(params['jj_junc_diam_vec_dense_100'])
        D_jjs_dense = Device('jj_series_arrays_dense')
        params_mod = copy.deepcopy(params)
        params_mod['jj_include_shunt'] = False
        params_mod['pad_size'] = params['inline_test_pad_size']
        params_mod['pad_pitch'] = [2*params_mod['pad_size'][0],params_mod['pad_size'][1]]
        params_mod['jj_include_flux_moats'] = False
        num_cols_jjs_vec = params['num_cols_jjs_vec_dense']
        for kk in range(len(num_cols_jjs_vec)):
            D_jjs_dense_row = Device('jj_series_arrays_dense_row')
            params_mod['num_cols_jjs'] = num_cols_jjs_vec[kk]
            for ii in range(num_jjs_dense):
                params_mod['jj_junc_diam'] = params['jj_junc_diam_vec_dense_100'][ii]
                jjs = jj_series_array(params_mod)
                jj_series_array_instance = D_jjs_dense_row.add_ref(jjs) 
                jj_series_array_instance.movex(4*ii*params_mod['pad_size'][0]+60)
            row_instance = D_jjs_dense.add_ref(D_jjs_dense_row)
            row_instance.movey(kk*485)
                     
        D_jjs_dense_text = Device('jj_series_dense_text')
        Text = pg.text(text = 'jj series arrays', size = params['block_label_size'], layer = vt_lyrs[params['label_layer']])
        text = D_jjs_dense_text.add_ref(Text)
        text.rotate(90)
        text.center = D_jjs_dense.center
        text.move([-D_jjs_dense.xsize/2-params['block_label_size'],0])
        D_jjs_dense.add_ref(D_jjs_dense_text)            
                
        jjs_dense = D_jj_b.add_ref(D_jjs_dense)
        jjs_dense.center = [40,-970] 
        
#%% single squids

if which_to_run['run_all'] == True or which_to_run['sq_single'] == True:
        
    port_list = ['p_s_a_1','p_s_a_5','p_e_a_1','p_e_a_5','p_n_a_1','p_n_a_5','p_w_a_1','p_w_a_5']
    
    if chip_to_run['sq_a'] == True:
        
        D_single_squids = Device('single_squids')
        params_mod = copy.deepcopy(params)
        params_mod['sq_ground_plane_hole'] = True
        params_mod['sq_wash_w_wire_wide'] = params['sq_wash_w_wire_wide_040']
        params_mod['sq_wash_w_in'] = params['sq_wash_w_in_040']
        params_mod['shunt_num_squares'] = params['sq_jj_shunt_num_squares_040']
        params_mod['jj_include_flux_moats'] = 'only_via_side'
        params_mod['sq_incoil_outset'] = params['sq_incoil_outset_040']
        params_mod['sq_washer_jj_offset'] = params['sq_washer_jj_offset_040']
        for pp in range(2):
            if pp == 0:
                params_mod['sq_include_shunt'] = True
                params_mod['sq_res_num_squares'] = params['sq_res_num_squares_040']
                params_mod['sq_res_w_wire'] = params['sq_res_w_wire_040']
            elif pp == 1:
                params_mod['sq_include_shunt'] = False
            for ii in range(len(params['sq_jj_junc_diam_vec_040'])):
                params_mod['jj_junc_diam'] = params['sq_jj_junc_diam_vec_040'][ii]
                sq = sq_8wire(params_mod)
                sq_inst = D_single_squids.add_ref(sq)
                sq_inst.connect(port = 'pad_anchor', destination = p_loc_sq_a.ports[port_list[pp*len(params['sq_jj_junc_diam_vec_040'])+ii]])
            
        D_sq_a.add_ref(D_single_squids)
    
    if chip_to_run['sq_b'] == True:
        
        D_single_squids = Device('single_squids')
        params_mod = copy.deepcopy(params)
        params_mod['sq_ground_plane_hole'] = False
        params_mod['sq_wash_w_wire_wide'] = params['sq_wash_w_wire_wide_100']
        params_mod['sq_wash_w_in'] = params['sq_wash_w_in_100']
        params_mod['shunt_num_squares'] = params['sq_jj_shunt_num_squares_100']
        params_mod['jj_include_flux_moats'] = 'only_via_side'
        params_mod['sq_incoil_outset'] = params['sq_incoil_outset_100']
        params_mod['sq_washer_jj_offset'] = params['sq_washer_jj_offset_100']
        for pp in range(2):
            if pp == 0:
                params_mod['sq_include_shunt'] = True
                params_mod['sq_res_num_squares'] = params['sq_res_num_squares_100']
                params_mod['sq_res_w_wire'] = params['sq_res_w_wire_100']
            elif pp == 1:
                params_mod['sq_include_shunt'] = False
            for ii in range(len(params['sq_jj_junc_diam_vec_100'])):
                params_mod['jj_junc_diam'] = params['sq_jj_junc_diam_vec_100'][ii]
                sq = sq_8wire(params_mod)
                sq_inst = D_single_squids.add_ref(sq)
                sq_inst.connect(port = 'pad_anchor', destination = p_loc_sq_a.ports[port_list[pp*len(params['sq_jj_junc_diam_vec_100'])+ii]])
            
        D_sq_b.add_ref(D_single_squids)         

#%% via series arrays 
    
if which_to_run['run_all'] == True or which_to_run['via_series'] == True:    
              
    if chip_to_run['sq_a'] == True:
        
        num_groups = len(params['num_cols_vias_vec'])
        D_via_series = Device('via_series_arrays')
        params_mod = copy.deepcopy(params)
        params_mod['pad_pitch'] = [175,200]
        params_mod['pad_size'] = [150,150]
        for kk in range(len(params['layer_below_vec'])):
            params_mod['layer_below'] = params['layer_below_vec'][kk]
            params_mod['via_layer'] = params['via_layer_vec'][kk]
            params_mod['layer_above'] = params['layer_above_vec'][kk]
            if params['layer_below_vec'][kk] == 'jj2':
                params_mod['inter_via_gap'] = 0
            for ii in range(num_groups):
                params_mod['num_cols_vias'] = params['num_cols_vias_vec'][ii]
                vias = via_series_array(params_mod)
                via_series_array_instance = D_via_series.add_ref(vias)
                via_series_array_instance.movex(2*(kk*num_groups+ii)*0.85*params['pad_pitch'][0])    
        
        D_via_text = Device('via_series_text')
        Text = pg.text(text = 'via series arrays', size = params['block_label_size'], layer = vt_lyrs[params['label_layer']])
        text = D_via_text.add_ref(Text)
        text.rotate(90)
        text.center = D_via_series.center
        text.move([-D_via_series.xsize/2-params['block_label_size'],0])
        D_via_series.add_ref(D_via_text)
            
        via_block = D_sq_a.add_ref(D_via_series)
        via_block.center = [0,0]
        via_block.move([0,-via_block.size[1]/2-1035])
              
    if chip_to_run['sq_b'] == True:
        
        num_groups = len(params['num_cols_vias_vec'])
        D_via_series = Device('via_series_arrays')
        params_mod = copy.deepcopy(params)
        params_mod['pad_pitch'] = [175,200]
        params_mod['pad_size'] = [150,150]
        for kk in range(len(params['layer_below_vec'])):
            params_mod['layer_below'] = params['layer_below_vec'][kk]
            params_mod['via_layer'] = params['via_layer_vec'][kk]
            params_mod['layer_above'] = params['layer_above_vec'][kk]
            if params['layer_below_vec'][kk] == 'jj2':
                params_mod['inter_via_gap'] = 0
            for ii in range(num_groups):
                params_mod['num_cols_vias'] = params['num_cols_vias_vec'][ii]
                vias = via_series_array(params_mod)
                via_series_array_instance = D_via_series.add_ref(vias)
                via_series_array_instance.movex(2*(kk*num_groups+ii)*0.85*params['pad_pitch'][0])    
        
        D_via_text = Device('via_series_text')
        Text = pg.text(text = 'via series arrays', size = params['block_label_size'], layer = vt_lyrs[params['label_layer']])
        text = D_via_text.add_ref(Text)
        text.rotate(90)
        text.center = D_via_series.center
        text.move([-D_via_series.xsize/2-params['block_label_size'],0])
        D_via_series.add_ref(D_via_text)
            
        via_block = D_sq_b.add_ref(D_via_series)
        via_block.center = [0,0]
        via_block.move([0,-via_block.size[1]/2-1035])        
    
#%% inline test
    
if which_to_run['run_all'] == True or which_to_run['inline_test'] == True:              
        
        inline_test = vt_inline_test(copy.deepcopy(params))
        
        ilt_jj_a = D_jj_a.add_ref(inline_test)    
        ilt_sq_a = D_sq_a.add_ref(inline_test)
        ilt_jj_b = D_jj_b.add_ref(inline_test)
        ilt_sq_b = D_sq_b.add_ref(inline_test)
        
        ilt_jj_a.center = [0,0]  
        ilt_jj_a.movex(230)
        ilt_jj_a.movey(inline_test.size[1]+170)
        ilt_jj_b.center = [0,0]  
        ilt_jj_b.movex(350)
        ilt_jj_b.movey(inline_test.size[1]+170)
        
        ilt_sq_a.rotate(-90)
        ilt_sq_a.center = [0,0]  
        ilt_sq_a.move([120,inline_test.size[1]+85])         
        ilt_sq_b.rotate(-90)
        ilt_sq_b.center = [0,0]  
        ilt_sq_b.move([-10,inline_test.size[1]+135])

#%% resistance tests

if which_to_run['run_all'] == True or which_to_run['resistance_test'] == True:
    
    if chip_to_run['sq_a'] == True:
        
        #just straight wires
        rt_coords = [0,-600]
        params_mod = copy.deepcopy(params)
        params_mod['pad_size'] = [120,120]
        res_test = vt_resistor_test_structures(params_mod)
        rtb = D_sq_a.add_ref(res_test)
        rtb.center = [0,0]
        rtb.move(rt_coords)
        
        #large numbers of resistors with small numbers of squares
        params_mod = copy.deepcopy(params)
        params_mod['pad_size'] = [120,120]
    #    params_mod['pad_pitch'] = [120,100]
        x_start = 780
        x_pitch = 320
        y_start = 400
        coords = [[x_start,y_start],[x_start,y_start+500],[x_start,y_start+1000],
                  [x_start+x_pitch,y_start-280],[x_start+x_pitch,y_start+260],[x_start+x_pitch,y_start+910],
                  [x_start+2*x_pitch+30,y_start-90],[x_start+2*x_pitch+30,y_start+790]]
        for ii in range(len(params['res_num_squares_vec_low_r'])):
            params_mod['res_num_squares'] = params['res_num_squares_vec_low_r'][ii]
            params_mod['res_w_wire'] = params['res_w_wire_vec'][ii]    
            params_mod['res_num_in_column'] = params['res_num_in_column_vec'][ii]                           
            D_res = res_series_array(params_mod)
            res = D_sq_a.add_ref(D_res)
            res.center = coords[ii]
            
        Text = pg.text(text = 'low\nresistance\ntests', size = params['block_label_size'], layer = vt_lyrs[params['label_layer']])    
        text = D_sq_a.add_ref(Text)         
        text.center = [x_start,y_start-280]
    
    if chip_to_run['sq_b'] == True:
        
        #just straight wires
        rt_coords = [0,-600]
        params_mod = copy.deepcopy(params)
        params_mod['pad_size'] = [120,120]
        res_test = vt_resistor_test_structures(params_mod)
        rtb = D_sq_b.add_ref(res_test)
        rtb.center = [0,0]
        rtb.move(rt_coords)
        
        #large numbers of resistors with small numbers of squares
        params_mod = copy.deepcopy(params)
        params_mod['pad_size'] = [120,120]
    #    params_mod['pad_pitch'] = [120,100]
        x_start = 780
        x_pitch = 320
        y_start = 400
        coords = [[x_start,y_start],[x_start,y_start+500],[x_start,y_start+1000],
                  [x_start+x_pitch,y_start-280],[x_start+x_pitch,y_start+260],[x_start+x_pitch,y_start+910],
                  [x_start+2*x_pitch+30,y_start-90],[x_start+2*x_pitch+30,y_start+790]]
        for ii in range(len(params['res_num_squares_vec_low_r'])):
            params_mod['res_num_squares'] = params['res_num_squares_vec_low_r'][ii]
            params_mod['res_w_wire'] = params['res_w_wire_vec'][ii]    
            params_mod['res_num_in_column'] = params['res_num_in_column_vec'][ii]                           
            D_res = res_series_array(params_mod)
            res = D_sq_b.add_ref(D_res)
            res.center = coords[ii]
            
        Text = pg.text(text = 'low\nresistance\ntests', size = params['block_label_size'], layer = vt_lyrs[params['label_layer']])    
        text = D_sq_b.add_ref(Text)         
        text.center = [x_start,y_start-280]        
    
#%% electrical isolation test
    
if which_to_run['run_all'] == True or which_to_run['isolation_test'] == True: 
    
    x_start = 1180
    y_start = 600

    if chip_to_run['jj_a'] == True:
    
        Combs = Device('electrical_isolation_tests')
        comb = Combs.add_ref(pg.test_comb(pad_size = [100,100], wire_width = 1, wire_gap = 3,
                      comb_layer = params['layers']['m0i'], overlap_zigzag_layer = params['layers']['jj1'],
                      comb_pad_layer = params['layers']['v0'], comb_gnd_layer = params['layers']['v0'], overlap_pad_layer = params['layers']['v2']))
        comb.center = [x_start,y_start]
        
        comb = Combs.add_ref(pg.test_comb(pad_size = [100,100], wire_width = 1, wire_gap = 3,
                      comb_layer = params['layers']['jj1'], overlap_zigzag_layer = params['layers']['m4'],
                      comb_pad_layer = params['layers']['v4'], comb_gnd_layer = params['layers']['v4'], overlap_pad_layer = params['layers']['v4']))
        comb.center = [x_start,y_start+400]
        
        D_jj_a.add_ref(Combs)

    if chip_to_run['jj_b'] == True:
    
        Combs = Device('electrical_isolation_tests')
        comb = Combs.add_ref(pg.test_comb(pad_size = [100,100], wire_width = 1, wire_gap = 3,
                      comb_layer = params['layers']['m0i'], overlap_zigzag_layer = params['layers']['jj1'],
                      comb_pad_layer = params['layers']['v0'], comb_gnd_layer = params['layers']['v0'], overlap_pad_layer = params['layers']['v2']))
        comb.center = [x_start,y_start]
        
        comb = Combs.add_ref(pg.test_comb(pad_size = [100,100], wire_width = 1, wire_gap = 3,
                      comb_layer = params['layers']['jj1'], overlap_zigzag_layer = params['layers']['m4'],
                      comb_pad_layer = params['layers']['v4'], comb_gnd_layer = params['layers']['v4'], overlap_pad_layer = params['layers']['v4']))
        comb.center = [x_start,y_start+400]
        
        D_jj_b.add_ref(Combs)
        
#%% litho tests
    
if which_to_run['run_all'] == True or which_to_run['litho_test'] == True:
        
    if chip_to_run['jj_a'] == True:
        
        params_mod = copy.deepcopy(params)
        params_mod['jj_junc_diam_vec'] = params['jj_junc_diam_vec_040']
        D_lt,D_lt_inv = vt_litho_tests(params_mod)
        litho_test_coords = [-1050,1170]
        d_lt = D_jj_a.add_ref(D_lt)
        d_lt.rotate(-90)
        d_lt.center = [litho_test_coords[0],litho_test_coords[1]]
        d_lt_inv = D_jj_a.add_ref(D_lt_inv)
        d_lt_inv.rotate(-90)
        d_lt_inv.center = [litho_test_coords[0],litho_test_coords[1]-615]

    if chip_to_run['jj_b'] == True:

        params_mod = copy.deepcopy(params)
        params_mod['jj_junc_diam_vec'] = params['jj_junc_diam_vec_100']
        D_lt,D_lt_inv = vt_litho_tests(params_mod)        
        litho_test_coords = [-900,820]
        d_lt = D_jj_b.add_ref(D_lt)
        d_lt.rotate(-90)
        d_lt.center = [litho_test_coords[0],litho_test_coords[1]]
#        d_lt_inv = D_jj_b.add_ref(D_lt_inv)
#        d_lt_inv.rotate(-90)
#        d_lt_inv.center = [litho_test_coords[0],litho_test_coords[1]-615]
        
    if chip_to_run['sq_a'] == True:        

        params_mod = copy.deepcopy(params)
        params_mod['jj_junc_diam_vec'] = params['jj_junc_diam_vec_040']
        D_lt,D_lt_inv = vt_litho_tests(params_mod)        
        litho_test_coords = [-1300,800]
        d_lt = D_sq_a.add_ref(D_lt)
        d_lt.center = [litho_test_coords[0],litho_test_coords[1]]
        d_lt_inv = D_sq_a.add_ref(D_lt_inv)
        d_lt_inv.center = [litho_test_coords[0]+615,litho_test_coords[1]]
        
    if chip_to_run['sq_b'] == True:        

        params_mod = copy.deepcopy(params)
        params_mod['jj_junc_diam_vec'] = params['jj_junc_diam_vec_100']
        D_lt,D_lt_inv = vt_litho_tests(params_mod)        
        litho_test_coords = [-1000,960]
        d_lt = D_sq_b.add_ref(D_lt)
        d_lt.center = [litho_test_coords[0],litho_test_coords[1]]
#        d_lt_inv = D_sq_b.add_ref(D_lt_inv)
#        d_lt_inv.center = [litho_test_coords[0]+615,litho_test_coords[1]]        
        
    
#%% endpoint boxes
        
if which_to_run['run_all'] == True or which_to_run['endpoint_boxes'] == True:
    
    params_mod = copy.deepcopy(params)
    params_mod['label_layer'] = 'm0l'
    endpoint_boxes = vt_endpoint_boxes(params_mod)
    epb_jj_a = D_jj_a.add_ref(endpoint_boxes)
    epb_jj_b = D_jj_b.add_ref(endpoint_boxes)
    epb_jj_a.center = [0,50]
    epb_jj_b.center = [0,50]
    epb_sq_a = D_sq_a.add_ref(endpoint_boxes)
    epb_sq_b = D_sq_b.add_ref(endpoint_boxes)
    epb_sq_a.center = [-450,50]
    epb_sq_b.center = [-450,50]
    
#%% extra pads

if which_to_run['run_all'] == True or which_to_run['extra_pads'] == True:
    
    D_extra_pads_a = Device('extra_pads')
    D_extra_pads_b = Device('extra_pads')
    params_mod = copy.deepcopy(params)
    params_mod['is_ground_pad'] = False
    D_pad = jj_pad(params_mod)
    params_mod['is_ground_pad'] = True
    D_pad_gnd = jj_pad(params_mod)
    
    port_list_a = ['p_s_a_9','p_s_a_10','p_e_a_9','p_e_a_10','p_n_a_9','p_n_a_10','p_w_a_9','p_w_a_10']
    for ii in range(len(port_list_a)):
        pad = D_extra_pads_a.add_ref(D_pad)
        pad.connect(port = 'pad_anchor', destination = p_loc_jj_a.ports[port_list_a[ii]])
    port_list_a = ['p_s_b_9','p_e_b_9','p_n_b_9','p_w_b_9']
    for ii in range(len(port_list_a)):
        pad = D_extra_pads_a.add_ref(D_pad_gnd)
        pad.connect(port = 'pad_anchor', destination = p_loc_jj_a.ports[port_list_a[ii]])        
    D_jj_a.add_ref(D_extra_pads_a)
    D_jj_b.add_ref(D_extra_pads_a)
    
    port_list_b = ['p_s_a_9','p_s_a_10','p_e_a_9','p_e_a_10','p_n_a_9','p_n_a_10','p_w_a_9','p_w_a_10']#['p_s_a_10','p_e_a_10','p_n_a_10','p_w_a_10']
    for ii in range(len(port_list_b)):
        pad = D_extra_pads_b.add_ref(D_pad_gnd)
        pad.connect(port = 'pad_anchor', destination = p_loc_sq_a.ports[port_list_b[ii]])
    port_list_b = ['p_s_b_9','p_e_b_9','p_n_b_9','p_w_b_9']#['p_s_a_10','p_e_a_10','p_n_a_10','p_w_a_10']
    for ii in range(len(port_list_b)):
        pad = D_extra_pads_b.add_ref(D_pad_gnd)
        pad.connect(port = 'pad_anchor', destination = p_loc_sq_a.ports[port_list_b[ii]])        
    D_sq_a.add_ref(D_extra_pads_b)
    D_sq_b.add_ref(D_extra_pads_b)

#%% nist logo / phi labels
    
if which_to_run['run_all'] == True or which_to_run['labels'] == True:
    
    nist_coords1 = [-2200,2325]
    nist_coords2 = [2175,2325]
    D_nist = pg.import_gds('nist_small.gds', cellname = None, flatten = True)
    D_nist.remap_layers(layermap = {0:vt_lyrs[params['label_layer']]})
    D_nist.add_port(name = 'south', midpoint = [D_nist.center[0], D_nist.ymin], width = 0, orientation = 270)
    D_nist.add_port(name = 'south_west', midpoint = [D_nist.xmin, D_nist.ymin], width = 0, orientation = 270)
    D_nist.add_port(name = 'south_east', midpoint = [D_nist.xmax, D_nist.ymin], width = 0, orientation = 270)
    
    nist_a1 = D_jj_a.add_ref(D_nist)
    nist_a1.center = nist_coords1
    nist_a2 = D_jj_a.add_ref(D_nist)
    nist_a2.center = nist_coords2
    
    nist_a1 = D_jj_b.add_ref(D_nist)
    nist_a1.center = nist_coords1
    nist_a2 = D_jj_b.add_ref(D_nist)
    nist_a2.center = nist_coords2    
    
    nist_b1 = D_sq_a.add_ref(D_nist)
    nist_b1.center = nist_coords1
    nist_b2 = D_sq_a.add_ref(D_nist)
    nist_b2.center = nist_coords2
    
    nist_b1 = D_sq_b.add_ref(D_nist)
    nist_b1.center = nist_coords1
    nist_b2 = D_sq_b.add_ref(D_nist)
    nist_b2.center = nist_coords2
    
    Phi_1 = Device('phi_label_1')
    Phi_2 = Device('phi_label_2')
    phi_text_size = params['block_label_size']
    y_space = phi_text_size/2
    text_str1 = 'physics and \nhardware for \ninformation'
    text_str2 = 'physics and \nhardware for \nintelligence'
    
    Tex1 = pg.text(text = text_str1, size = phi_text_size, justify = 'left', layer = vt_lyrs[params['label_layer']])
    Tex2 = pg.text(text = text_str2, size = phi_text_size, justify = 'right', layer = vt_lyrs[params['label_layer']])
    Phi_1.add_ref(Tex1)
    Phi_1.add_port(name = 'north', midpoint = [Phi_1.center[0], Phi_1.ymax], width = 0, orientation = 90) 
    Phi_1.add_port(name = 'north_west', midpoint = [Phi_1.xmin, Phi_1.ymax], width = 0, orientation = 90)
    Phi_2.add_ref(Tex2)
    Phi_2.add_port(name = 'north', midpoint = [Phi_2.center[0], Phi_2.ymax], width = 0, orientation = 90) 
    Phi_2.add_port(name = 'north_east', midpoint = [Phi_2.xmax, Phi_2.ymax], width = 0, orientation = 90)
    
    tex_a_1 = D_jj_a.add_ref(Phi_1)         
    tex_a_1.connect(port = 'north_west', destination = nist_a1.ports['south_west'])
    tex_a_1.movey(-y_space)
    
    tex_a_2 = D_jj_a.add_ref(Phi_2)         
    tex_a_2.connect(port = 'north_east', destination = nist_a2.ports['south_east'])
    tex_a_2.movey(-y_space)    
     
    tex_a_1 = D_jj_b.add_ref(Phi_1)         
    tex_a_1.connect(port = 'north_west', destination = nist_a1.ports['south_west'])
    tex_a_1.movey(-y_space)
    
    tex_a_2 = D_jj_b.add_ref(Phi_2)         
    tex_a_2.connect(port = 'north_east', destination = nist_a2.ports['south_east'])
    tex_a_2.movey(-y_space) 
    
    tex_b_1 = D_sq_a.add_ref(Phi_1)         
    tex_b_1.connect(port = 'north_west', destination = nist_b1.ports['south_west'])
    tex_b_1.movey(-y_space)
    
    tex_b_2 = D_sq_a.add_ref(Phi_2)         
    tex_b_2.connect(port = 'north_east', destination = nist_b2.ports['south_east'])
    tex_b_2.movey(-y_space)    
    
    tex_b_1 = D_sq_b.add_ref(Phi_1)         
    tex_b_1.connect(port = 'north_west', destination = nist_b1.ports['south_west'])
    tex_b_1.movey(-y_space)
    
    tex_b_2 = D_sq_b.add_ref(Phi_2)         
    tex_b_2.connect(port = 'north_east', destination = nist_b2.ports['south_east'])
    tex_b_2.movey(-y_space) 
    
#%% votan labels
    
if which_to_run['run_all'] == True or which_to_run['labels'] == True:
    
    votan_coords = [0,1450]
    votan_text_size = 100
    
    text_str = 'phi : vt01 jj a'
    Tex = pg.text(text = text_str, size = votan_text_size, layer = vt_lyrs[params['label_layer']])
    tex = D_jj_a.add_ref(Tex)            
    tex.center = votan_coords
    tex.movex(175)
    
    text_str = 'phi : vt01 jj b'
    Tex = pg.text(text = text_str, size = votan_text_size, layer = vt_lyrs[params['label_layer']])    
    tex = D_jj_b.add_ref(Tex)            
    tex.center = votan_coords
    tex.movex(345)    
    
    text_str = 'phi : vt01 sq a'
    Tex = pg.text(text = text_str, size = votan_text_size, layer = vt_lyrs[params['label_layer']])
    tex = D_sq_a.add_ref(Tex)            
    tex.center = votan_coords
    tex.movex(-330)
        
    text_str = 'phi : vt01 sq b'
    Tex = pg.text(text = text_str, size = votan_text_size, layer = vt_lyrs[params['label_layer']])
    tex = D_sq_b.add_ref(Tex)            
    tex.center = votan_coords  
    tex.move([-10,70])

#%% shainline label
    
if which_to_run['run_all'] == True or which_to_run['labels'] == True:
    
    shainline_text_size = params['block_label_size']
    text_str = params['chip_signature']
    Tex = pg.text(text = text_str, size = shainline_text_size, justify = 'right', layer = vt_lyrs[params['label_layer']])
    tex = D_jj_a.add_ref(Tex)            
    tex.center = [params['chip_size'][0]/2-tex.xsize/2-45,-params['chip_size'][1]/2+tex.ysize/2+45]
    tex = D_jj_b.add_ref(Tex)            
    tex.center = [params['chip_size'][0]/2-tex.xsize/2-45,-params['chip_size'][1]/2+tex.ysize/2+45]    
    tex = D_sq_a.add_ref(Tex)            
    tex.center = [params['chip_size'][0]/2-tex.xsize/2-45,-params['chip_size'][1]/2+tex.ysize/2+45]   
    tex = D_sq_b.add_ref(Tex)            
    tex.center = [params['chip_size'][0]/2-tex.xsize/2-45,-params['chip_size'][1]/2+tex.ysize/2+45]

#%% make die
if chip_to_run['jj_a'] == True:
    make_die_vt(D_jj_a,params['cleave_street_width'],params['cleave_street_length'],vt_lyrs[params['label_layer']],vt_lyrs['ce'],chip_label = params['chip_label_jj_a'],label_size = params['block_label_size'],die_size = params['chip_size'])
    D_jj_a.write_gds('vt01jja.gds')

if chip_to_run['jj_b'] == True:    
    make_die_vt(D_jj_b,params['cleave_street_width'],params['cleave_street_length'],vt_lyrs[params['label_layer']],vt_lyrs['ce'],chip_label = params['chip_label_jj_b'],label_size = params['block_label_size'],die_size = params['chip_size'])
    D_jj_b.write_gds('vt01jjb.gds')

if chip_to_run['sq_a'] == True:
    make_die_vt(D_sq_a,params['cleave_street_width'],params['cleave_street_length'],vt_lyrs[params['label_layer']],vt_lyrs['ce'],chip_label = params['chip_label_sq_a'],label_size = params['block_label_size'],die_size = params['chip_size'])
    D_sq_a.write_gds('vt01sqa.gds')

if chip_to_run['sq_b'] == True:
    make_die_vt(D_sq_b,params['cleave_street_width'],params['cleave_street_length'],vt_lyrs[params['label_layer']],vt_lyrs['ce'],chip_label = params['chip_label_sq_b'],label_size = params['block_label_size'],die_size = params['chip_size'])
    D_sq_b.write_gds('vt01sqb.gds')

elapsed = time.time() - t_tot
print('initial gds build duration = '+str(elapsed)+' s')

#%% post-processing

if which_to_run['post_processing_jj_a'] == True:
    D_jj_a_post, layer_list = vt_data_prep(D_jj_a, new_device_name = 'votan_1__jj_a', params = params)
    D_jj_a_post.write_gds('vt01jja_post.gds')
    
if which_to_run['post_processing_jj_b'] == True:    
    D_jj_b_post, layer_list = vt_data_prep(D_jj_b, new_device_name = 'votan_1__jj_b', params = params)
    D_jj_b_post.write_gds('vt01jjb_post.gds')

if which_to_run['post_processing_sq_a'] == True:
    D_sq_a_post, layer_list = vt_data_prep(D_sq_a, new_device_name = 'votan_1__sq_a', params = params)
    D_sq_a_post.write_gds('vt01sqa_post.gds')
    
if which_to_run['post_processing_sq_b'] == True:    
    D_sq_b_post, layer_list = vt_data_prep(D_sq_b, new_device_name = 'votan_1__sq_b', params = params)
    D_sq_b_post.write_gds('vt01sqb_post.gds')    