#%% init
import numpy as np
import time
import copy

import phidl.geometry as pg
from phidl import make_device, Device, Layer, LayerSet, quickplot2 as qp
#from nc_library import *

from vt_util import vt_layers
from vt_params__spd_res import vt_parameters
from nc_library__vt_util import vt_fill, vt_fill_array

vt_lyrs,layer_data = vt_layers()
params  = vt_parameters()

#%%
def vt_data_prep(D_in, new_device_name = 'spd_res', params = params):
    
    layers = vt_lyrs
    
    print('\n\npost processing...\n')
    
    precision = params['post_processing_precision']
    fill_exclude_offset = params['fill_exclude_offset']
    
    t_tot = time.time()
    
    D_out = Device(new_device_name)
    
    #chip_edge
    chip_edge = pg.extract(D_in,[layers['ce']])
    D_out.add_ref(chip_edge)

    #make fill
    print('generating fill...')
    t_fill = time.time()
    D_fill = vt_fill_array(params)
    D_fill.center = D_in.center
    print('fill generation time = '+str(time.time()-t_fill)+' s\n')
    
    #fill exclude  
    print('calculating fill exclusions...')
    t_fill = time.time()
    fill_exclude_layers_list = ['m1','m1e','m2','m3','m3cs','stf','stfe','stfp','r1','r1p','r2','r2p']
    layer_list = []
    for ii in range(len(fill_exclude_layers_list)):
        layer_list.append(layers[fill_exclude_layers_list[ii]])
    fill_exclude_region = pg.extract(D_in,layer_list)
    fill_exclude_region.flatten(single_layer = layers['dp'])
    fill_exclude_region = pg.boolean(fill_exclude_region,fill_exclude_region,'or',precision = precision,layer = layers['dp'])
    fill_exclude_region = pg.offset(fill_exclude_region, distance = fill_exclude_offset, join_first = True, precision = precision, layer = layers['dp'])

#    fill_region_m1 = pg.extract(D_fill,[layers['m1']])
#    fill_region_m1 = pg.boolean(fill_region_m1,pg.offset(fill_exclude_region,distance = 0*params['fill_w_delta'], join_first = True, precision = precision, layer = layers['dp']),'A-B',precision = precision,layer = layers['m0f'])            
    fill_region_stf = pg.extract(D_fill,[layers['stf']])
    fill_region_stf = pg.boolean(fill_region_stf,pg.offset(fill_exclude_region,distance = 0*params['fill_w_delta'], join_first = True, precision = precision, layer = layers['dp']),'A-B',precision = precision,layer = layers['stff'])            
    fill_region_r1 = pg.extract(D_fill,[layers['r1']])
    fill_region_r1 = pg.boolean(fill_region_r1,pg.offset(fill_exclude_region,distance = 1*params['fill_w_delta'], join_first = True, precision = precision, layer = layers['dp']),'A-B',precision = precision,layer = layers['r1f'])        
    fill_region_r2 = pg.extract(D_fill,[layers['r2']])
    fill_region_r2 = pg.boolean(fill_region_r2,pg.offset(fill_exclude_region,distance = 2*params['fill_w_delta'], join_first = True, precision = precision, layer = layers['dp']),'A-B',precision = precision,layer = layers['r2f'])
    print('fill exclusion calculation time = '+str(time.time()-t_fill)+' s\n')

    #m1 / gound plane
    print('m1/gp...')
    t_m1 = time.time()
    m1 = pg.extract(D_in,[layers['m1'],layers['m1e'],layers['m2'],layers['m3cs']])#,layers['m1p']
    m1invert = pg.extract(D_in,[layers['m1i']])
    m3 = pg.extract(D_in,[layers['m3']])
    device_layers = pg.extract(D_in,[layers['stf'],layers['r1'],layers['r2'],layers['m3']])
    device_layers_offset = pg.offset(device_layers, distance = params['ground_plane_buffer'], precision = precision, layer = layers['dp'])
    device_layers_pads = pg.extract(D_in,[layers['stfp'],layers['r1p'],layers['r2p']])
    device_layers_offset = pg.boolean(device_layers_offset,device_layers_pads,'A-B',precision = precision,layer = layers['m1'])
    m1 = pg.boolean(m1,device_layers_offset,'A+B',precision = precision,layer = layers['m1'])
    m1 = pg.boolean(m1,m3,'A-B',precision = precision,layer = layers['m1'])
    m1 = pg.boolean(m1,m1invert,'A+B',precision = precision,layer = layers['m1'])
    m1i = pg.copy_layer(m1, layer = layers['m1'], new_layer = layers['m1i'])

    #gp holes
    params_mod = copy.deepcopy(params)
    params_mod['fill_w_widest'] = params['ground_plane_hole_size']
    params_mod['fill_gap'] = params['ground_plane_hole_pitch']-params['ground_plane_hole_size']
    params_mod['fill_layer_list'] = ['dp']
    D_holes_gp = vt_fill(params_mod)
    D_holes_gp.center = D_in.center
    D_holes_gp = pg.boolean(D_holes_gp,fill_exclude_region,'A-B',precision = precision,layer = layers['dp'])
    m1 = pg.boolean(m1,D_holes_gp,'A+B',precision = precision,layer = layers['m1'])
    
    if params['post_processing_elim_small_features'] == True:
        D_out.add_ref(vt_eliminate_small_features(m1,params['post_processing_min_dim'],precision,layers['m1']))
        m1i = pg.copy_layer(D_out, layer = layers['m1'], new_layer = layers['m1i'])
        D_out.add_ref(m1i)
    else:
        D_out.add_ref(m1)
        m1i = pg.copy_layer(D_out, layer = layers['m1'], new_layer = layers['m1i'])
        D_out.add_ref(m1i)
    print('gp calculation time = '+str(time.time()-t_m1)+' s\n')
    
    #stf
    print('stf...')
    invert_stf = False
    t_stf = time.time()
    stf = pg.extract(D_in,[layers['stf'],layers['stfp']])
    if invert_stf == True:
        stf_inverted = pg.boolean(chip_edge,stf,'A-B',precision = precision,layer = layers['stf'])
        stf_inverted = pg.boolean(stf_inverted,fill_region_stf,'A-B',precision = precision,layer = layers['stf'])
        if params['post_processing_elim_small_features'] == True:
            D_out.add_ref(vt_eliminate_small_features(stf_inverted,params['post_processing_min_dim'],precision,layers['stf']))
        else:
            D_out.add_ref(stf_inverted)
    elif invert_stf == False:
        stf = pg.boolean(stf,fill_region_stf,'A+B',precision = precision,layer = layers['stf'])
        if params['post_processing_elim_small_features'] == True:
            D_out.add_ref(vt_eliminate_small_features(stf,params['post_processing_min_dim'],precision,layers['stf']))
        else:
            D_out.add_ref(stf)    
    print('stf calculation time = '+str(time.time()-t_stf)+' s\n')
    
    #r1
    print('r1...')
    t_res = time.time()
    r1 = pg.extract(D_in,[layers['r1'],layers['r1p']])
    r1 = pg.boolean(r1,r1,'A+B',precision = precision,layer = layers['r1'])
    r1 = pg.boolean(r1,fill_region_r1,'A+B',precision = precision,layer = layers['r1'])
    if params['post_processing_elim_small_features'] == True:
        D_out.add_ref(vt_eliminate_small_features(r1,params['post_processing_min_dim'],precision,layers['r1']))
    else:
        D_out.add_ref(r1)
    print('resistor calculation time = '+str(time.time()-t_res)+' s\n')
    
    #r2
    print('r2...')
    t_res = time.time()
    r2 = pg.extract(D_in,[layers['r2'],layers['r2p']])
    r2 = pg.boolean(r2,r2,'A+B',precision = precision,layer = layers['r2'])
    r2 = pg.boolean(r2,fill_region_r2,'A+B',precision = precision,layer = layers['r2'])
    if params['post_processing_elim_small_features'] == True:
        D_out.add_ref(vt_eliminate_small_features(r2,params['post_processing_min_dim'],precision,layers['r2']))
    else:
        D_out.add_ref(r2)
    print('resistor calculation time = '+str(time.time()-t_res)+' s\n')    
    
    elapsed = time.time() - t_tot
    print('total post processing duration = '+str(elapsed)+' s')
    
    return D_out

def vt_eliminate_small_features(element, feature_size = 0.1, precision = 1e-6, layer = 1):
    
    #remove small bits
    element = pg.offset(element, distance = -feature_size/2, join_first = True, precision = precision, layer = layer)
    element = pg.boolean(element,element,'or',precision = precision,layer = layer)
    element = pg.offset(element, distance = feature_size/2, join_first = True, precision = precision, layer = layer)
    
    #remove small gaps
    element = pg.offset(element, distance = feature_size/2, join_first = True, precision = precision, layer = layer)
    element = pg.boolean(element,element,'or',precision = precision,layer = layer)
    element = pg.offset(element, distance = -feature_size/2, join_first = True, precision = precision, layer = layer)
    
    return element

def vt_add_via_layers(main_device,via_layer_list,layers):
    
    Vias = Device('vias') 
    for ii in range(len(via_layer_list)):
        via = pg.extract(main_device,[layers[via_layer_list[ii]],layers[via_layer_list[ii]+'e']])
        via.flatten(single_layer = layers[via_layer_list[ii]])
        Vias.add_ref(via)
    
    return Vias