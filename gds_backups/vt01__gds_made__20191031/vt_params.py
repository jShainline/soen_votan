# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 16:18:51 2019

@author: jms4
"""

from phidl import LayerSet
import numpy as np


#%%
def vt_parameters(layers):
    
    jc_vec = np.array([2e7,1.5e7,1e7,0.5e7])#np.flipud(np.linspace(0.5e7,2e7,8))
    d_vec_040 = np.around(1e6*np.sqrt(4*40e-6/(np.pi*jc_vec)), decimals = 2)#all entries ending in _040 are for 40uA Ic junctions
    d_vec_100 = np.around(1e6*np.sqrt(4*100e-6/(np.pi*jc_vec)), decimals = 2)#all entries ending in _100 are for 100uA Ic junctions
    jj_params = dict(jj_junc_diam = 2.26, #targeting 40uA Ic with 1kA/cm^2 jc
                     jj_num_pts = 40,
                     jj_junc_diam_vec_040 = d_vec_040,#these values are chosen to try to hit 40uA Ic if jc ranges from 500A/cm^2 to 2kA/cm^2
                     jj_junc_diam_vec_100 = d_vec_100,#these values are chosen to try to hit 100uA Ic if jc ranges from 500A/cm^2 to 2kA/cm^2
                     jj_junc_diam_vec_dense_040 = np.around(np.linspace(d_vec_040[0],d_vec_040[-1],8), decimals = 2),
                     jj_junc_diam_vec_dense_100 = np.around(np.linspace(d_vec_100[0],d_vec_100[-1],8), decimals = 2),
                     jj_via_inset = 0.3,#chosen to make shorting to jj1 unlikely
                     jj_top_contact_outset = 1, 
                     jj_bottom_contact_outset = 1.2,
                     jj_include_shunt = True,
                     jj_include_flux_moats = True,
                     include_50_ohm = True,
                     shunt_w_wire = 2,#jj shunt resistor
                     shunt_num_squares = 2.05,#jj shunt resistor
                     shunt_num_squares_vec_040 = [1.15,2.05],#jj shunt resistor, first number is for beta_c = 0.3, second number for beta_c = 0.95
                     shunt_num_squares_vec_100 = [0.45,0.8],#jj shunt resistor, first number is for beta_c = 0.3, second number for beta_c = 0.95
                     shunt_w_cntct = 2,#jj shunt resistor
                     shunt_outset = 0.2,#jj shunt resistor
                     shunt_l_lead = 1,#jj shunt resistor
                     inter_jj_gap = 1.3,#for series arrays
                     num_jjs_per_col = 20,#for series arrays, divisible by 2
                     num_cols_jjs = 4,#for series arrays, divisible by 2
                     num_cols_jjs_vec = [4,6],#for series arrays, divisible by 2
                     num_cols_jjs_vec_dense = [2,4,8],#for series arrays, divisible by 2
                     w_wire = 2,#misc wiring
                     via_width = 2,                                      
                     l_lead = 10,                     
                     layers = layers)
    
    squid_params = dict(sq_jj_junc_diam_vec_040 = d_vec_040[:],
                        sq_jj_junc_diam_vec_100 = d_vec_100[:],
                        sq_jj_shunt_num_squares_040 = 1.15,
                        sq_jj_shunt_num_squares_100 = 0.45,
                        sq_wash_w_wire_wide_040 = 20,                        
                        sq_wash_w_wire_wide_100 = 9,
                        sq_wash_w_in_040 = 10,#16.5,#we need L = 26 pH; L = 1.25 mu0 d where d is sq_wash_w_in (inner diameter of square washer) when width of wire is greater than d. This gives d = Phi0*beta_L/(2.5*mu0*Ic), which gives 16.5um. I'm using 13 to try to be safe with beta_L definitely less than 1 even if the slit adds inductance
                        sq_wash_w_in_100 = 40,
                        sq_wash_slit_gap = 2,
                        sq_l_leads = 20,
                        sq_incoil_w_wire = 2,#input coil
                        sq_incoil_outset = 0.5,#input coil
                        sq_incoil_outset_040 = 5,#input coil
                        sq_incoil_outset_100 = 0.5,#input coil
                        sq_washer_jj_offset = 5,#how far back to place washer from jjs
                        sq_washer_jj_offset_040 = 5,#how far back to place washer from jjs
                        sq_washer_jj_offset_100 = 5,#how far back to place washer from jjs
                        sq_incoil_l_leads = 10,#input coil
                        sq_incoil_lead_gap = 4,#input coil
                        sq_addflux_include = True,#additional input for operating in flux-locked loop
                        sq_addflux_w_wire = 2,#additional input for operating in flux-locked loop
                        sq_addflux_outset = 2,#outset beyond incoil
                        sq_addflux_extent = 20,#length of line
                        sq_include_via = False,#True or False. need yes for series array. for single four wire can be either
                        sq_res_num_squares = 0.1,#for resonance-damping resistor (this one is kind of complicated, but I'm just shooting for an L/r long-pass filter with cutoff around 2 GHz (f_c = R/2piL with L the self-inductance of the washer (the jjs aren't in the filter)). I don't know exactly where the LC resonances will be)
                        sq_res_num_squares_040 = 0.1,#for resonance-damping resistor (this one is kind of complicated, but I'm just shooting for an L/r long-pass filter with cutoff around 2 GHz (f_c = R/2piL with L the self-inductance of the washer (the jjs aren't in the filter)). I don't know exactly where the LC resonances will be)
                        sq_res_num_squares_100 = 0.2,#for resonance-damping resistor (this one is kind of complicated, but I'm just shooting for an L/r long-pass filter with cutoff around 2 GHz (f_c = R/2piL with L the self-inductance of the washer (the jjs aren't in the filter)). I don't know exactly where the LC resonances will be)
                        sq_res_w_wire = 15,#for resonance-damping resistor
                        sq_res_w_wire_040 = 10,#for resonance-damping resistor
                        sq_res_w_wire_100 = 5,#for resonance-damping resistor
                        sq_res_w_cntct = 1,#for resonance-damping resistor 
                        sq_res_outset = 0.5,#for resonance-damping resistor 
                        sq_res_l_lead = 2,#for resonance-damping resistor
                        sq_res_w_lead = 3,#for resonance-damping resistor
                        sq_res_50_ohm_lateral_gap = 10,#gap between two vertically running 50 ohm resistors
                        sq_include_shunt = True,#for resonance-damping resistor
                        sq_ground_plane_hole = False,#make hole in ground plane under squid to decrease inductance
                        sq_gph_outset = 15,#extent of ground plane hole beyond squid washer
                        sq_washer_include_inductor_ports = False,#for inductex/fasthenry
                        sq_incoil_include_inductor_ports = False,#for inductex/fasthenry
                        sq_squid_include_inductor_ports = False,#for inductex/fasthenry
                        )
    
    resistor_params = dict(w_wire_vec = [1,2],#for resistance tests
                           res_num_squares_vec = [1,2,4,8,16,32,64],#for resistance tests
                           res_num_squares_vec_meander = [1000,4000],#for resistance tests
                           res_num_squares_vec_meander_Nb = [1000,4000],#for resistance tests
                           res_num_squares = 0.1,#for low-resistance series arrays
                           res_num_squares_vec_low_r = [0.1,0.2,0.3,0.1,0.1,0.1,0.01,0.02],#for low-resistance series arrays
                           res_w_wire = 10,#for low-resistance series arrays
                           res_w_wire_vec = [10,10,10,10,10,10,60,50],#for low-resistance series arrays
                           res_w_cntct = 1,#for low-resistance series arrays
                           res_outset = 0.2,#for low-resistance series arrays
                           res_l_lead = 1,#for low-resistance series arrays
                           res_space = 2,#for low-resistance series arrays
                           res_num_in_column = 10,#divisible by 2, for low-resistance series arrays
                           res_num_in_column_vec = [10,10,10,20,30,40,10,10],#divisible by 2, for low-resistance series arrays
                           res_num_columns = 10,#divisible by 2, for low-resistance series arrays
                           res_50_ohm_w_wire = 2,
                           res_50_ohm_num_squares = 25,
                           res_50_ohm_outset = 0.2,
                           )
    
    misc_params = dict(layer_below = 'jj1',#for via series arrays
                       layer_below_vec = ['m0i','jj1','jj2'],#for via series arrays
                       via_layer = 'v2',#for via series arrays
                       via_layer_vec = ['v0','v2','v2'],#for via series arrays
                       layer_above = 'm4',#for via series arrays
                       layer_above_vec = ['jj1','m4','m4'],#for via series arrays
                       num_vias_per_col = 20,#for viaseries arrays
                       num_cols_vias = 10,#for series arrays, divisible by 2
                       num_cols_vias_vec = [2,4,8],#for via series arrays
                       inter_via_gap = 1,#for series arrays   
                       block_label_size = 40,
                       pad_label_size = 20,
                       device_label_size = 6,
                       label_layer = 'm4l',
                       ground_plane_hole_size = 5,
                       ground_plane_hole_pitch = 20,
                       ground_plane_moat_width = 5,
                       ground_plane_buffer = 10,
                       ground_plane_moat_width_fine = 1,
                       ground_plane_buffer_fine = 10,
                       include_big_moats = False,
                       include_fine_moats = True
                       )
    
    pad_params = dict(pad_size = [200,220],
                     pad_opening_inset = 1,
                     pad_pitch = [200,375],
                     pad_y_backset = 75,
                     pad_w_wire = 30,
                     pad_w_ground_strap = 35,
                     inline_test_pad_size = [90,90],
                     pad_gp_gap = 15)
    
    chip_params = dict(chip_size = [5000,5000],
                       cleave_street_width = 25,
                       cleave_street_length = 125,
                       chip_label_jj_a = 'NIST PHI Votan1\n40uA JJs',
                       chip_label_jj_b = 'NIST PHI Votan1\n100uA JJs',
                       chip_label_sq_a = 'NIST PHI Votan1\n40uA SQUIDs',
                       chip_label_sq_b = 'NIST PHI Votan1\n100uA SQUIDs',
                       chip_signature = 'Shainline\n20191031',
                       endpoint_box_width = 250,
                       endpoint_box_outline = 10,
                       fill_w_widest = 20, 
                       fill_gap = 20,
                       fill_w_delta = 1, 
                       fill_layer_list = ['jj1f','jj2f','m4f','rf'],
                       fill_exclude_offset = 15,
                       x_inset = 600,#how far lower left pad is inset from chip edge
                       y_inset = 25,
                       post_processing_min_dim = 0.01,
                       post_processing_precision = 1e-6,
                       post_processing_elim_small_features = False
                       )
    
    params = {**jj_params,**squid_params,**pad_params,**chip_params,**resistor_params,**misc_params}
    
    return params