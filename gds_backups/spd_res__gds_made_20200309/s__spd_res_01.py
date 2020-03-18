from __future__ import division, print_function, absolute_import
import numpy as np

#import phidl
#import phidl.geometry as pg
from phidl import Device
from f__physical_constants import physical_constants
#from phidl import quickplot2 as qp

from vt_util import vt_layers, vt_layers_post, write_lyp
from vt_params__spd_res import vt_parameters
from vt_post_processing__spd_res import vt_data_prep

import copy
import time

from nc_library__vt_util import vt_make_die, vt_inline_test, vt_litho_tests, vt_endpoint_boxes
from nc_library__vt_spd import vt_spd_meander, vt_spd_sq
from nc_library__vt_res import vt_resistor_meander, vt_resistor_4_wire, vt_res_stitch_4_wire, vt_resistor_series_array

    
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
params_mod['chip_label_small'] = params['chip_label_small']+'\nSPD / Res'
params_mod['ground_pad_coords'] = [[2250,1430],[2250,-1450],[-2250,1430],[-2250,-1450]]
params_mod['ground_pad_size'] = [[400,1300],[400,1300],[400,1300],[400,1300]]
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
port_list = ['p_s_c_1','p_s_c_2','p_s_c_3','p_s_c_4','p_s_c_5','p_s_c_6','p_s_c_7','p_s_c_8','p_s_c_9']
params_mod = copy.deepcopy(params)
for ii in range(len(port_list)):
    D_spd = vt_spd_meander(params_mod)
    spd = D_spd_res_full.add_ref(D_spd)
    spd.connect(port = 'pad_anchor', destination = p_loc_spd_res.ports[port_list[ii]])  
                    
#resistance tests    
    
#resistor meanders
port_list = ['p_n_a_1','p_n_a_2','p_n_a_3','p_n_a_4','p_n_a_6','p_n_a_7','p_n_a_8','p_n_a_9']
params_mod = copy.deepcopy(params)
res_layer_list = ['r1','r1','r1','r1','r2','r2','r2','r2']
for ii in range(len(res_layer_list)):
    params_mod['res_layer'] = res_layer_list[ii]
    params_mod['res_num_squares_meander'] = params['res_num_squares_meander_vec'][ii]
    D_res = vt_resistor_meander(params_mod)
    res = D_spd_res_full.add_ref(D_res)
    res.connect(port = 'pad_anchor', destination = p_loc_spd_res.ports[port_list[ii]])

#resistor 4 wire
params_mod = copy.deepcopy(params)
port_list = ['p_n_b_1','p_n_b_3','p_n_c_1','p_n_c_3']
res_layer_list = ['r1','r1','r2','r2']        
for ii in range(len(res_layer_list)):#range(len(params_mod['res_layer_list'])):
    params_mod['res_layer'] = res_layer_list[ii]
    params_mod['res_resistance'] = params['res_resistance_vec__4wire'][ii]           
    D_res = vt_resistor_4_wire(params_mod)
    res = D_spd_res_full.add_ref(D_res)
    res.connect(port = 'pad_anchor', destination = p_loc_spd_res.ports[port_list[ii]])

#resistor stitch 4 wire
params_mod = copy.deepcopy(params)
port_list = ['p_n_b_6','p_n_b_8','p_n_c_6','p_n_c_8']
res_layer_list = ['r1','r1','r2','r2']        
for ii in range(len(res_layer_list)):#range(len(params_mod['res_layer_list'])):
    params_mod['res_layer'] = res_layer_list[ii]
    params_mod['res_resistance'] = params['res_resistance_vec__4wire'][ii]           
    D_res = vt_res_stitch_4_wire(params_mod)
    res = D_spd_res_full.add_ref(D_res)
    res.connect(port = 'pad_anchor', destination = p_loc_spd_res.ports[port_list[ii]])
    
#low-resistance series arrays (large numbers of resistors with small numbers of squares)
params_mod = copy.deepcopy(params)
params_mod['pad_pitch'] = [350,0]
res_layer_list = ['r1','r1','r1','r1','r1','r1','r2','r2','r2','r2','r2','r2','stf','stf','stf','stf','stf','stf']
init_coords = [1950,250]
tn = 0
x_coord = 2000
y_coord = 360
for ii in range(len(params['res_num_squares_vec_low_r'])):
    if ii == 6:
        x_coord = 2000
        y_coord = -80
    if ii == 12:
        x_coord = 2000
        y_coord = -520
    params_mod['res_layer'] = res_layer_list[ii]
    params_mod['res_num_squares'] = params['res_num_squares_vec_low_r'][ii]                         
    D_res = vt_resistor_series_array(params_mod)
    res = D_spd_res_full.add_ref(D_res)
    res.rotate(180)
    res.center = [x_coord,y_coord]
    x_coord = x_coord-750
                   
# write gds
D_spd_res_full.write_gds('spd_res.gds') 

t_spd_stop = time.time()
print('time to build spd_res = '+str(t_spd_stop-t_spd_start)+' s')       
    

#%% total time    
elapsed = time.time() - t_tot
print('initial gds build duration = '+str(elapsed)+' s')

#%% post-processing
D_spd_res_full_post = vt_data_prep(D_spd_res_full, new_device_name = 'spd_res', params = params)
D_spd_res_full_post.write_gds('spd_res_post.gds')
