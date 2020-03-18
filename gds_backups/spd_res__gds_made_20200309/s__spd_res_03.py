from __future__ import division, print_function, absolute_import
import numpy as np

#import phidl
#import phidl.geometry as pg
#C:\Users\jms4\AppData\Local\Continuum\anaconda3\Lib\site-packages\phidl
from phidl import Device
from f__physical_constants import physical_constants
#from phidl import quickplot2 as qp

from vt_util import vt_layers, vt_layers_post, write_lyp
from vt_params__spd_res import vt_parameters
from vt_post_processing__spd_res import vt_data_prep

import copy
import time

from nc_library__vt_util import vt_make_die, vt_inline_test, vt_litho_tests, vt_endpoint_boxes, vt_fill_array
from nc_library__vt_spd import vt_spd_meander, vt_spd_sq
from nc_library__vt_res import vt_resistor_meander, vt_resistor_4_wire, vt_res_stitch_4_wire, vt_resistor_series_array, vt_res_interdigitated_3_wire

    
#%% misc setup
print('\n\nbuilding gds...')
t_tot = time.time()

p = physical_constants()

vt_lyrs,layer_data = vt_layers()
vt_lyrs_post,layer_data_post = vt_layers_post()
write_lyp(vt_lyrs,layer_data,'vt')
write_lyp(vt_lyrs_post,layer_data_post,'vt_post')
params  = vt_parameters()

#%% spd_res
                
t_spd_start = time.time()
print('building spd_res ...')

D_spd_res_full = Device('votan2__spd_res')

#initialize die
params_mod = copy.deepcopy(params)
params_mod['die_type'] = 'vt02_two-sides'
params_mod['chip_label'] = params['chip_label']+' spd / res'
params_mod['chip_label_small'] = params['chip_label_small']
params_mod['ground_pad_coords'] = [[2250,1530],[2250,-1500],[-2250,1530],[-2250,-1500]]
params_mod['ground_pad_size'] = [[400,1100],[400,1300],[400,1100],[400,1300]]
params_mod['pad_row_spacing_vec'] = [580,520,780,560]
D_spd_res_full, p_loc_spd_res = vt_make_die(D_spd_res_full,params_mod)

#spd meanders

#vary w_wire
port_list = ['p_s_a_1','p_s_a_2','p_s_a_3','p_s_a_4','p_s_a_5','p_s_a_6','p_s_a_7','p_s_a_8','p_s_a_9']
params_mod = copy.deepcopy(params)
for ii in range(len(params_mod['spd_w_wire_vec'])):            
    params_mod['spd_w_wire'] = np.around(params_mod['spd_w_wire_vec'][ii], decimals = 2)
    params_mod['spd_wire_pitch'] = np.around(2*params_mod['spd_w_wire_vec'][ii], decimals = 2)
    D_spd = vt_spd_meander(params_mod)
    spd = D_spd_res_full.add_ref(D_spd)
    spd.connect(port = 'pad_anchor', destination = p_loc_spd_res.ports[port_list[ii]]) 
#vary inductance
port_list = ['p_s_b_1','p_s_b_2','p_s_b_3','p_s_b_4','p_s_b_5','p_s_b_6','p_s_b_7','p_s_b_8','p_s_b_9']
params_mod = copy.deepcopy(params)
for ii in range(len(params_mod['spd_inductance_vec'])):
    params_mod['spd_inductance'] = params_mod['spd_inductance_vec'][ii]
    D_spd = vt_spd_meander(params_mod)
    spd = D_spd_res_full.add_ref(D_spd)
    spd.connect(port = 'pad_anchor', destination = p_loc_spd_res.ports[port_list[ii]])        
#consistency tests
#port_list = ['p_s_c_1','p_s_c_2','p_s_c_3','p_s_c_4','p_s_c_5','p_s_c_6','p_s_c_7','p_s_c_8','p_s_c_9']
#params_mod = copy.deepcopy(params)
#for ii in range(len(port_list)):
#    D_spd = vt_spd_meander(params_mod)
#    spd = D_spd_res_full.add_ref(D_spd)
#    spd.connect(port = 'pad_anchor', destination = p_loc_spd_res.ports[port_list[ii]])  
                    
#resistance tests      
    
#resistor meanders
x_start = -2065
y_start = -1400
x_delta = 410
x_coord = x_start
y_coord = y_start
x_extra = 140
params_mod = copy.deepcopy(params)
for jj in range(len(params['res_layer_list'])):
    params_mod['res_layer'] = params['res_layer_list'][jj]
    x_coord = x_coord+x_extra
    for ii in range(3):
        params_mod['res_num_squares_meander'] = params['res_num_squares_meander_vec'][jj*3+ii]
        D_res = vt_resistor_meander(params_mod)
        res = D_spd_res_full.add_ref(D_res)
        res.xmin = x_coord
        res.ymin = y_coord
        x_coord = x_coord+x_delta
        
#res-stitch series arrays
params_mod = copy.deepcopy(params)
params_mod['pad_pitch'] = [320,0]
res_layer_list = ['r1','r2','stf']
x_start = -2300
y_start = -770
x_delta = 690
y_delta = 410
x_coord = x_start
y_coord = y_start
for jj in range(len(res_layer_list)):       
    params_mod['res_layer'] = res_layer_list[jj]
    for ii in range(len(params['res_num_squares_vec_low_r'])): 
        if res_layer_list[jj] == 'r1' or 'r2' or 'm3':
            params_mod['res_num_squares'] = params['res_num_squares_vec_low_r'][ii]
        if res_layer_list[jj] == 'stf':
            params_mod['res_num_squares'] = params['spd_num_sqrs_cntct_array']
            params_mod['res_num_columns'] = params['spd_num_columns_vec'][ii]
            params_mod['res_num_in_column'] = params['spd_num_in_column_vec'][ii]
        D_res = vt_resistor_series_array(params_mod)
        res = D_spd_res_full.add_ref(D_res)
        res.xmin = x_coord
        res.ymin = y_coord
        x_coord = x_coord+x_delta
    x_coord = x_start
    y_coord = y_coord+y_delta    
    
#resistors interdigitated    
params_mod = copy.deepcopy(params)
res_layer_list = ['r1','r2'] 
params_mod['pad_pitch'] = [320,0]
x_start = -2030
y_start = 450
x_delta = 640
y_delta = 465
x_coord = x_start
y_coord = y_start
for jj in range(len(res_layer_list)):       
    params_mod['res_intdig_layer'] = res_layer_list[jj]        
    for ii in range(len(params['res_intdig_num_digits_vec'])):
        params_mod['res_intdig_num_digits'] = params['res_intdig_num_digits_vec'][ii] 
        params_mod['res_intdig_l_digit'] = params['res_intdig_l_digit_vec'][ii]          
        D_res = vt_res_interdigitated_3_wire(params_mod)
        res = D_spd_res_full.add_ref(D_res)
        res.xmin = x_coord
        res.ymin = y_coord
        y_coord = y_coord+y_delta
    x_coord = x_start+x_delta
    y_coord = y_start

##resistor 4 wire
params_mod = copy.deepcopy(params)
res_layer_list = ['r1','r2','m3','stf'] 
params_mod['pad_pitch'] = [310,0]  
x_start = -700
y_start = 680
x_delta = 640
y_delta = 450
x_coord = x_start
y_coord = y_start
for jj in range(len(res_layer_list)):       
    params_mod['res_layer'] = res_layer_list[jj]  
    temp_vec = params['res_resistance_vec__4wire__'+params_mod['res_layer']]
    for ii in range(len(temp_vec)):#range(len(params_mod['res_layer_list'])): 
        params_mod['res_resistance'] = temp_vec[ii]           
        D_res = vt_resistor_4_wire(params_mod)
        res = D_spd_res_full.add_ref(D_res)
        res.xmin = x_coord
        res.ymin = y_coord
        x_coord = x_coord+x_delta
    x_coord = x_start
    y_coord = y_coord+y_delta

#resistor stitch 4 wire
params_mod = copy.deepcopy(params)
res_layer_list = ['r1','r2','m3','stf'] 
params_mod['pad_pitch'] = [310,0]  
x_start = 690
y_start = 680
x_delta = 640
y_delta = 450
x_coord = x_start
y_coord = y_start
for jj in range(len(res_layer_list)):       
    params_mod['res_layer'] = res_layer_list[jj]
    temp_vec = params['res_resistance_vec__4wire__'+params_mod['res_layer']]      
    for ii in range(len(temp_vec)):
        params_mod['res_resistance'] = temp_vec[ii]           
        D_res = vt_res_stitch_4_wire(params_mod)
        res = D_spd_res_full.add_ref(D_res)
        res.xmin = x_coord
        res.ymin = y_coord
        x_coord = x_coord+x_delta
    x_coord = x_start
    y_coord = y_coord+y_delta
       
#endpoint boxes
vt_lyrs,layer_data = vt_layers()
import phidl.geometry as pg
D_ept = pg.rectangle(size = [300,300], layer = vt_lyrs['m1e'])
ept = D_spd_res_full.add_ref(D_ept)
ept.center = [2200,-300]
D_ept = pg.rectangle(size = [300,300], layer = vt_lyrs['stfe'])
ept = D_spd_res_full.add_ref(D_ept)
ept.center = [2200,430]        
                   
# write gds
D_spd_res_full.write_gds('spd_res_03.gds') 

t_spd_stop = time.time()
print('time to build spd_res = '+str(t_spd_stop-t_spd_start)+' s')       
    

#%% total time    
elapsed = time.time() - t_tot
print('initial gds build duration = '+str(elapsed)+' s')

#%% post-processing
D_spd_res_full_post = vt_data_prep(D_spd_res_full, new_device_name = 'spd_res', params = params)
D_spd_res_full_post.write_gds('spd_res_03_post.gds')
