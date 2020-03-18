# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 16:18:51 2019

@author: jms4
"""

from phidl import LayerSet
import numpy as np


#%%
def vt_parameters(layers):
    
#    jc_vec = np.array([2e7,1.5e7,1e7,0.75e7,0.5e7])#np.flipud(np.linspace(0.5e7,2e7,8))
#    d_vec = np.around(1e6*np.sqrt(4*40e-6/(np.pi*jc_vec)), decimals = 2)#all entries ending in _040 are for 40uA Ic junctions  
    jc_target = 1e7#amps per m^2
    ic_target_fq = 40e-6
    ic_target_fq_vec = np.array([20e-6,35e-6,40e-6,45e-6,60e-6])
    ic_target_sq = 100e-6
    ic_target_sq_vec = np.array([80e-6,95e-6,100e-6,105e-6,120e-6])
    
    d_vec_fq = np.around(1e6*np.sqrt(4*ic_target_fq_vec/(np.pi*jc_target)), decimals = 2)
    d_vec_sq = np.around(1e6*np.sqrt(4*ic_target_sq_vec/(np.pi*jc_target)), decimals = 2)
    
    jj_params = dict(jj_ic_target = ic_target_fq,
                     jj_junc_diam = np.around(1e6*np.sqrt(4*ic_target_fq/(np.pi*jc_target)), decimals = 2), #targeting 40uA Ic with 1kA/cm^2 jc
                     jj_cap_density = 11.5e-15,#farads per microamp
                     jj_jc_target = jc_target,
                     jj_shunt_resistance = 4.11,
                     jj_shunt_w_wire = 2,#jj shunt resistor
                     jj_num_pts = 40,
                     jj_junc_diam_vec = d_vec_fq,
                     jj_via_inset = 0.3,#chosen to make shorting to jj1 unlikely
                     jj_top_contact_outset = 1, 
                     jj_bottom_contact_outset = 1.2,
                     jj_include_shunt = True,
                     jj_include_flux_moats = False,
                     jj_xcntct_length = 7,
                     jj_res = 4.125,#jj resistance in normal state
                     jj_beta_c = 0.95,
                     include_50_ohm = True,                     
                     inter_jj_gap = 1.3,#for series arrays
                     num_jjs_per_col = 20,#for series arrays, divisible by 2
                     num_cols_jjs = 4,#for series arrays, divisible by 2
                     w_wire = 2.04,#misc wiring                                    
                     l_lead = 30,                     
                     )
    
    squid_params = dict(sq_jj_ic_target = ic_target_fq,
                        sq_jj_junc_diam = np.around(1e6*np.sqrt(4*ic_target_sq/(np.pi*jc_target)), decimals = 2),
                        sq_jj_junc_diam_vec = d_vec_sq,
                        sq_jj_shunt_resistance = 0.93,
                        sq_jj_shunt_w_wire = 4,#jj shunt resistor
                        sq_jj_beta_c = 0.3,
                        sq_wash_w_wide = 20, 
                        sq_wash_w_bridge = 6, 
                        sq_wash_w_in1 = 6,#16.5,#we need L = 26 pH; L = 1.25 mu0 d where d is sq_wash_w_in (inner diameter of square washer) when width of wire is greater than d. This gives d = Phi0*beta_L/(2.5*mu0*Ic), which gives 16.5um. I'm using 13 to try to be safe with beta_L definitely less than 1 even if the slit adds inductance                        
                        sq_wash_w_in2 = 12,
                        sq_wash_l_bridge = 3,
                        sq_wash_slit_gap = 2,
                        sq_wash_jj_offset = 5,#how far back to place washer from jjs
                        sq_l_leads = 10,#must be less than sq_incoil_distance_outside_washer
                        sq_lead_gap = 11.1,
                        sq_bias_y_offset = 3.7,#how far up washer to make bias contact
                        sq_incoil_w_wire = 1.4,#input coil
                        sq_incoil_w_pitch = 3,#input coil
                        sq_incoil_numturns = 6,#input coil
                        sq_incoil_via_backset = 1,#how far into the washer opening is the via?
                        sq_incoil_outset = 1,#input coil extent outside washer opening                        
                        sq_incoil_distance_outside_washer = 15,#input coil extent wire coming back is spaced away from washer edge
                        sq_addflux_include = True,#additional input for operating in flux-locked loop
                        sq_addflux_w_wire = 2.3,#additional input for operating in flux-locked loop
                        sq_addflux_outset = 4.1,#outset beyond incoil
                        sq_addflux_extent = 65,#length of line
                        sq_include_via = False,#True or False. need yes for series array. for single four wire can be either
                        sq_shunt_resistance = 0.06,#for resonance-damping resistor (this one is kind of complicated, but I'm just shooting for an L/r long-pass filter with cutoff around 2 GHz (f_c = R/2piL with L the self-inductance of the washer (the jjs aren't in the filter)). I don't know exactly where the LC resonances will be)                        
                        sq_shunt_w_wire = 15,#for resonance-damping resistor 
                        sq_shunt_l_lead = 2,#for resonance-damping resistor
                        sq_shunt_w_lead = 3,#for resonance-damping resistor
                        sq_shunt_outset = 2,#for resonance-damping resistor
                        sq_include_shunt = True,#for resonance-damping resistor
                        sq_res_50_ohm_lateral_gap = 10,#gap between two vertically running 50 ohm resistors
                        sq_ground_plane_hole = False,#make hole in ground plane under squid to decrease inductance
                        sq_gph_outset = 15,#extent of ground plane hole beyond squid washer
                        sq_washer_include_inductor_ports = False,#for inductex/fasthenry
                        sq_incoil_include_inductor_ports = False,#for inductex/fasthenry
                        sq_squid_include_inductor_ports = False,#for inductex/fasthenry
                        )
    
    spd_params = dict(spd_w_wire = 1.0,
                      spd_w_wire_vec = np.linspace(0.6,1.5,10),
                      spd_wire_pitch = 2.0,
                      spd_inductance = 1e-6,
                      spd_inductance_vec = [1.75e-6,2e-6,2.25e-6,2.5e-6,3e-6],
                      spd_sq_inductance_vec = [500e-9,825e-9],
                      spd_turn_ratio = 3,                      
                      spd_ground_plane_buffer = 10,
                      spd_tau = 50e-9,
                      spd_res_y_offset = 21,#how far down current bias wire to make contact to resistor to Jsf
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
                           res_num_squares_vec_low_r = [0.1,0.1,0.1,0.01,0.01,0.01,0.005,0.005,0.005],#for low-resistance series arrays 
                           res_w_wire_vec_low_r = [10,10,10,50,50,50,200,200,200],#for low-resistance series arrays
                           res_num_in_column = 6,#for low-resistance series arrays, divisible by 2
                           res_num_in_column_vec = [20,20,20,6,6,6,6,6,6],#for low-resistance series arrays, divisible by 2
                           res_num_columns = 4,#for low-resistance series arrays, divisible by 2
                           res_num_columns_vec = [8,16,32,8,16,32,8,16,32],#for low-resistance series arrays, divisible by 2
                           res_space = 2,#for low-resistance series arrays   
                           res_w_cntct = 1,#for low-resistance series arrays
                           res_outset = 0.2,
                           res_wire_olap = 3,
                           res_l_lead = 1,
                           res_50_ohm_w_wire = 1.4,
                           res_50_ohm_outset = 0.2,
                           res_per_sq_r1 = 0.1,
                           res_per_sq_r2 = 2,
                           sy_res_w_wire = 1.4,
                           res_label_size = 4,
                           )
    
    inductor_params = dict(induct_w_wire = 1,
                           induct_wire_pitch = 2,
                           induct_num_squares = 390,
                           induct_num_squares_sfq1 = 244,#left inductor in DCSFQ
                           induct_num_squares_sfq2 = 144,#right inductor in DCSFQ
                           induct_num_squares_jtl = 244,#jtls in sffg
                           induct_flux_purge = 10e-9,
                           induct_l_pre = 10,
                           induct_turn_ratio = 2,
                           induct_layer = 'm3',
                           inductance_per_sq_stf = 200e-12,
                           inductance_per_sq_m1 = 200e-15,
                           inductance_per_sq_m3 = 200e-15,
                           induct_include_inductex_ports = False)
    
    sfg_params = dict(sfg_spd_inductance = 825e-9,
                      spd_jj_extra_x = 35,
                      sfg_spd_tau = 50e-9,
                      jj_inductor_extra_x = 25,
                      sfg_si_w_wire = 2,                      
                      sfg_si_inductance = 8e-9,
                      sfg_si_turn_ratio = 3,
                      
                      sfg_tau_unique_label_vec = ['md_L   lg_tau','md_L   md_tau','md_L   sm_tau','sm_L   md_tau','md_L   md_tau','lg_L   md_tau'],
                      sfg_tau_si_inductance_vec = [247e-9,247e-9,247e-9,8.2e-9,82.3e-9,823e-9],
                      sfg_tau_si_w_wire_vec = [2,2,2,2,2,2],
                      sfg_si_tau_vec = [20e-6,2e-6,200e-9,2e-6,2e-6,2e-6],
                      sfg_w_wire_tau_vec = [20,4,2,4,4,4],
                      sfg_w_wire_sq_res_connect_vec = [18,3,1.8,3,3,3],
                      sfg_tau_si_turn_ratio_vec = [3,3,3,2,3,3],
                      sfg_w_wire_tau = 20,
                      sfg_w_wire_sq_res_connect = 2.1,
                      sfg_si_tau = 8e-6,
                      
                      sfg_inf_unique_label_vec = ['sm_L   tau_inf','md_L   tau_inf','lg_L  tau_inf','sfq sm_L tau_inf','sfq md_L tau_inf','sfq lg_L tau_inf'],
                      sfg_inf_si_w_wire_vec = [2,2,2,2,2,2],
                      sfg_inf_si_inductance_vec = [8.2e-9,82.3e-9,823e-9,800e-12,8.2e-9,82.3e-9],
                      sfg_inf_si_turn_ratio_vec = [2,3,3,2,2,3],
                      
                      sfg_w_wire_flux_purge = 1,#1.01,
                      sfg_w_res_flux_purge = 1.92,
                      sfg_resistance_flux_purge = 22,
                      sfg_flux_purge_jj_offset = 12,
                      sfg_flux_purge_layer = 'r2',
                      sfg_include_tau_leak = True,
                      sfg_has_spd = True,#if false, spd is replaced by an electrical port for driving dcsfq
                      sfg_has_sfq = False,#toggles between standard sfg and version where spd goes into dcsfq
                      )
    
    sfq_params = dict(sfq_induct_w_wire_vec = [2,1.5,5,2,1.5,5],
                      sfq_induct_inductance_vec = [80e-12,800e-12,8e-9,80e-12,800e-12,8e-9],                      
                      sfq_induct_turn_ratio_vec = [1,1,3,1,1,3],                     
                      sfq_include_tau_leak = False,
                      sfq_has_spd_vec = [True,True,True,False,False,False],#if false, spd is replaced by an electrical port for driving dcsfq
                      sfq_has_sfq = True,#toggles between standard sfg and version where spd goes into dcsfq
                      sfq_unique_label_vec = ['spd-sfq _ sm_L','spd-sfq _ md_L','spd-sfq _ lg_L','dcsfq _ sm_L','dcsfq _ md_L','dcsfq _ lg_L',]
                      )
    
    via_params = dict(via_width = 2,
                      via_bottom_contact_outset = 1,
                      via_top_contact_outset = 0.6,
                      via_layer = 'v1',
                      via_layer_below = 'm1',
                      via_layer_above = 'm2o',
                      via_layer_vec = ['v1','v2','v3','v3'],#for via series arrays
                      via_layer_below_vec = ['m1','m2o','jj1','jj2'],#for via series arrays
                      via_layer_above_vec = ['m2o','jj1','m3','m3'],#for via series arrays
                      num_cols_vias = 10,#for via series arrays, divisible by 2
                      num_vias_per_col = 20,#for via series arrays
                      inter_via_gap = 1,#for series arrays                        
                      )
    
    pad_params = dict(pad_size = [200,220],
                      pad_size_ground = [32,32],
                      pad_opening_inset = 0.5,
                      pad_pitch = [400,830],
                      pad_y_backset = 100,
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
                       device_label_size = 6,
                       label_layer = 'm3l',
                       ground_plane_hole_size = 5,
                       ground_plane_hole_pitch = 20,
                       ground_plane_moat_width = 5,
                       ground_plane_buffer = 10,
                       ground_plane_moat_width_fine = 1,
                       ground_plane_buffer_fine = 10,
                       include_big_moats = False,
                       include_fine_moats = True,
                       )
    
    chip_params = dict(chip_size = [5000,5000],
                       chip_edge_layer = 'ce',
                       die_type = 'vt02_perimeter',
                       cleave_street_width = 25,
                       cleave_street_length = 125,
                       chip_label = 'phi : vt02',
                       chip_label_small = 'NIST PHI\nVotan2',
                       group_label_1 = 'physics and \nhardware for \ninformation',
                       group_label_2 = 'physics and \nhardware for \nintelligence',
                       chip_signature = 'Shainline\n20200124',
                       endpoint_box_width = 250,
                       endpoint_box_outline = 10,
                       fill_w_widest = 20, 
                       fill_gap = 20,
                       fill_w_delta = 1, 
                       fill_layer_list = ['jj1f','jj2f','m3f','rf'],
                       fill_exclude_offset = 15,                       
                       post_processing_min_dim = 0.01,
                       post_processing_precision = 1e-6,
                       post_processing_elim_small_features = False,
                       )
    
    params = {**jj_params,**squid_params,**spd_params,**via_params,**pad_params,**chip_params,**resistor_params,**misc_params,**sfg_params,**sfq_params,**inductor_params,**packaging_params}
    
    return params