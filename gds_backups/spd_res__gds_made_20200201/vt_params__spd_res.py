# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 16:18:51 2019

@author: jms4
"""

from phidl import LayerSet
import numpy as np


#%%
def vt_parameters():
    
    spd_params = dict(spd_w_wire = 0.8,
                      spd_w_wire_vec = np.linspace(0.5,1.3,9),
                      spd_wire_pitch = 1.6,
                      spd_inductance = 1.4e-6,
                      spd_inductance_vec = np.linspace(0.7e-6,2e-6,9),
                      spd_sq_inductance_vec = [500e-9,825e-9],
                      spd_turn_ratio = 3,                      
                      spd_ground_plane_buffer = 10,
                      spd_tau = 50e-9,
                      spd_res_y_offset = 21,#how far down current bias wire to make contact to resistor to Jsf
                      spd_pad_w_wire = 20,
                      spd_sy_l_lead = 10,
                      spd_sy_w_lead = 5,
                      spd_sy_m1_olap = 5,
                      spd_sy_m1_outset = 1,
                      spd_sy_m1_res_l_extra = 20,
                      spd_sy_include_res = True,
                      spd_contact_taper_length = 10,                      
                      spd_sq_extra_x = 35,#for spd-to-squid
                      spd_sy_res_w_wire_vec = [2,1.2],
                      spd_sq_tau_vec = [30e-9,40e-9,50e-9,60e-9,100e-9,200e-9],                      
                      spd_cntrl_inductance = 1e-6,                     
                      )
    
    resistor_params = dict(shunt_w_cntct = 2,#jj shunt resistor
                           shunt_outset = 0.2,#jj shunt resistor
                           shunt_l_lead = 1,#jj shunt resistorres_layer = 'r1',#for resistance tests
                           res_layer_list = ['m1','stf','r1','jj1','m3','r2'],#for resistance tests
                           res_w_wire = 2,#for resistance tests                       
                           res_num_squares = 13,#for resistance tests
                           res_num_squares_meander = 4000,#for resistance tests
                           res_num_squares_meander_vec = [2500,5000,10000,20000,250,500,1000,2000],#for resistance tests
                           res_num_squares_vec_low_r = [0.5,1,2,4,8,16,0.5,1,2,4,8,16],#for low-resistance series arrays 
                           res_num_in_column = 4,#for low-resistance series arrays, divisible by 2
                           res_num_columns = 4,#for low-resistance series arrays, divisible by 2
                           res_space = 2,#for low-resistance series arrays   
                           res_w_cntct = 1,#for low-resistance series arrays
                           res_outset = 0.2,
                           res_wire_olap = 3,
                           res_l_lead = 1,
                           res_50_ohm_w_wire = 1.4,
                           res_50_ohm_outset = 0.2,
                           res_per_sq_r1 = 0.2,
                           res_per_sq_r2 = 2,
                           sy_res_w_wire = 1.4,
                           res_label_size = 4,
                           res_resistance_vec__4wire = [0.2,2,2,20],
                           res_pad_olap = 5,
                           )


    pad_params = dict(pad_size = [200,220],
                      pad_size_ground = [32,32],
                      pad_opening_inset = 0.5,
                      pad_pitch = [450,830],
                      pad_y_backset = 80,
                      pad_gnd_y_backset = 30,
                      pad_w_wire = 10,
                      pad_gp_gap = 15,
                      pad_pkg_gap = 25,
                      inline_test_pad_size = [204,204],
                      pad_x_inset = 600,#how far left and right pads are inset from chip edge
                      pad_y_inset = 70.5,#how far lower and upper pads are inset from chip edge
                      )
    
    packaging_params = dict(fiber_collar_num_pts = 40,
                            fiber_collar_diam = 135,
                            fiber_collar_w_wall = 50,
                            fiber_collar_w_channel = 50,
                            fiber_collar_l_channel = 50,
                            fiber_collar_w_box = 350,
                            fiber_collar_l_box = 700,
                            fiber_collar_outer_corner_length = 50,
                            fiber_collar_inner_corner_length = 20,
                            fiber_collar_layer = 'pkg',
                            fiber_collar_orientation = 180,
                            pad_collar_w_wall = 30,
                            )
    
    misc_params = dict(block_label_size = 40,
                       pad_label_size = 20,
                       device_label_size = 12,
                       label_layer = 'm1',
                       ground_plane_hole_size = 2,
                       ground_plane_hole_pitch = 20,
                       ground_plane_moat_width = 5,
                       ground_plane_buffer = 5,
                       ground_plane_buffer_fine = 2,
                       ground_plane_moat_width_fine = 1,
                       include_big_moats = False,
                       include_fine_moats = True,
                       )
    
    chip_params = dict(chip_size = [5000,5000],
                       chip_edge_layer = 'ce',
                       die_type = 'vt02_perimeter',
                       cleave_street_width = 25,
                       cleave_street_length = 125,
                       chip_label = 'phi :',
                       chip_label_small = 'NIST PHI\nVotan2',
                       group_label_1 = 'physics and \nhardware for \ninformation',
                       group_label_2 = 'physics and \nhardware for \nintelligence',
                       chip_signature = 'Shainline\n20200201',
                       endpoint_box_width = 250,
                       endpoint_box_outline = 10,
                       fill_w_widest = 20, 
                       fill_gap = 20,
                       fill_w_delta = 1, 
                       fill_layer_list = ['stf','r1','r2'],
                       fill_exclude_offset = 15,                       
                       post_processing_min_dim = 0.4,
                       post_processing_precision = 1e-6,
                       post_processing_elim_small_features = True,
                       )
    
    params = {**spd_params,**pad_params,**chip_params,**resistor_params,**misc_params,**packaging_params}
    
    return params