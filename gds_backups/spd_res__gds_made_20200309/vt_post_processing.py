#%% init
import numpy as np
import time
import copy

import phidl.geometry as pg
from phidl import make_device, Device, Layer, LayerSet, quickplot2 as qp
#from nc_library import *

from vt_util import *
from vt_params import vt_parameters

vt_lyrs,layer_data = vt_layers()
params  = vt_parameters(vt_lyrs)

#%%
def vt_data_prep(D_in, new_device_name = 'vt01', params = params):
    
    layers = params['layers']
    
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
    D_fill = vt_fill(**params)
    D_fill.center = D_in.center
    print('fill generation time = '+str(time.time()-t_fill)+' s\n')
    
    #fill exclude  
    print('calculating fill exclusions...')
    t_fill = time.time()
    fill_exclude_layers_list = ['m0','m0e','m0l','m0i','v0','v0e','jj1','jj1e','jj2','jj2e','v2','v2e','res','m4','m4e','m4l','v4','v4e']
    layer_list = []
    for ii in range(len(fill_exclude_layers_list)):
        layer_list.append(layers[fill_exclude_layers_list[ii]])
    fill_exclude_region = pg.extract(D_in,layer_list)
    fill_exclude_region.flatten(single_layer = layers['dp'])
    fill_exclude_region = pg.boolean(fill_exclude_region,fill_exclude_region,'or',precision = precision,layer = layers['dp'])
    fill_exclude_region = pg.offset(fill_exclude_region, distance = fill_exclude_offset, join_first = True, precision = precision, layer = layers['dp'])

#    fill_region_m0 = pg.extract(D_fill,[layers['m0']])
#    fill_region_m0 = pg.boolean(fill_region_m0,pg.offset(fill_exclude_region,distance = 0*params['fill_w_delta'], join_first = True, precision = precision, layer = layers['dp']),'A-B',precision = precision,layer = layers['m0f'])            
    fill_region_jj1 = pg.extract(D_fill,[layers['jj1f']])
    fill_region_jj1 = pg.boolean(fill_region_jj1,pg.offset(fill_exclude_region,distance = 0*params['fill_w_delta'], join_first = True, precision = precision, layer = layers['dp']),'A-B',precision = precision,layer = layers['jj1f'])            
    fill_region_jj2 = pg.extract(D_fill,[layers['jj2f']])
    fill_region_jj2 = pg.boolean(fill_region_jj2,pg.offset(fill_exclude_region,distance = 1*params['fill_w_delta'], join_first = True, precision = precision, layer = layers['dp']),'A-B',precision = precision,layer = layers['jj2f'])        
    fill_region_m4 = pg.extract(D_fill,[layers['m4f']])
    fill_region_m4 = pg.boolean(fill_region_m4,pg.offset(fill_exclude_region,distance = 2*params['fill_w_delta'], join_first = True, precision = precision, layer = layers['dp']),'A-B',precision = precision,layer = layers['m4f'])
    fill_region_res = pg.extract(D_fill,[layers['rf']])
    fill_region_res = pg.boolean(fill_region_res,pg.offset(fill_exclude_region,distance = 3*params['fill_w_delta'], join_first = True, precision = precision, layer = layers['dp']),'A-B',precision = precision,layer = layers['rf'])    
    print('fill exclusion calculation time = '+str(time.time()-t_fill)+' s\n')

    #m0
    print('ground plane...')
    t_gp = time.time()
    m0 = pg.extract(D_in,[layers['m0'],layers['m0e']])
    jj1_e = pg.extract(D_in,[layers['jj1e']])
    jj1_e = pg.outline(jj1_e,params['endpoint_box_outline'],precision = precision, layer = layers['m0'])
    jj2_e = pg.extract(D_in,[layers['jj2e']])
    jj2_e = pg.outline(jj2_e,params['endpoint_box_outline'],precision = precision, layer = layers['m0'])
    m4_e = pg.extract(D_in,[layers['m4e']])
    m4_e = pg.outline(m4_e,params['endpoint_box_outline'],precision = precision, layer = layers['m0'])
    m0 = pg.boolean(m0,jj1_e,'A+B',precision = precision,layer = layers['m0'])
    m0 = pg.boolean(m0,jj2_e,'A+B',precision = precision,layer = layers['m0'])
    m0 = pg.boolean(m0,m4_e,'A+B',precision = precision,layer = layers['m0'])
    params_mod = copy.deepcopy(params)
    params_mod['fill_w_widest'] = params['ground_plane_hole_size']
    params_mod['fill_gap'] = params['ground_plane_hole_pitch']-params['ground_plane_hole_size']
    params_mod['fill_layer_list'] = ['dp']
    D_holes_gp = vt_fill(**params_mod)
    D_holes_gp.center = D_in.center 
#    D_holes_gp = pg.boolean(chip_edge,D_holes_gp,'A-B',precision,layers['dp'])
    D_holes_gp = pg.boolean(D_holes_gp,fill_exclude_region,'A-B',precision = precision,layer = layers['dp'])
    D_holes_gp = pg.boolean(D_holes_gp,pg.extract(D_in,[layers['m0l']]),'A+B',precision = precision,layer = layers['dp'])
    m0 = pg.boolean(m0,D_holes_gp,'A+B',precision = precision,layer = layers['m0'])
    
    #ground plane inverted regions
    m0i = pg.extract(D_in,[layers['m0i']])
    m0i_offset = pg.offset(m0i,distance = params['ground_plane_buffer'], join_first = True, precision = precision, layer = layers['m0i'])
    m0i_outline = pg.boolean(m0i_offset,m0i,'A-B',precision = precision,layer = layers['m0'])
    m0 = pg.boolean(m0,m0i_outline,'A+B',precision = precision,layer = layers['m0'])
    
#    m0i = pg.extract(D_in,[layers['m0i']])
##    m0_invert = pg.deepcopy(m0)
##    m0_invert = pg.boolean(m0_invert,m0i,'A-B',precision = precision,layer = layers['m0'])
#    m0_inverted = pg.invert(m0i, border = params['ground_plane_buffer'], precision = precision, layer = layers['m0'])
#    m0 = pg.boolean(m0,m0_inverted,'A+B',precision = precision,layer = layers['m0'])
    
    m0m = pg.extract(D_in,[layers['m0m']])#m0 moats
    m0 = pg.boolean(m0,m0m,'A+B',precision = precision,layer = layers['m0'])
    if params['post_processing_elim_small_features'] == True:
        D_out.add_ref(vt_eliminate_small_features(m0,params['post_processing_min_dim'],precision,layers['m0']))
    else:
        D_out.add_ref(m0)
    print('gp calculation time = '+str(time.time()-t_gp)+' s\n')
    
    #jj1
    print('jj1...')
    t_jj1 = time.time()
    jj1 = pg.extract(D_in,[layers['jj1']])
    jj1_inverted = pg.boolean(chip_edge,jj1,'A-B',precision = precision,layer = layers['jj1'])
    jj1_inverted_fill_clearout = pg.boolean(jj1_inverted,fill_region_jj1,'A-B',precision = precision,layer = layers['jj1'])
    if params['post_processing_elim_small_features'] == True:
        D_out.add_ref(vt_eliminate_small_features(jj1_inverted_fill_clearout,params['post_processing_min_dim'],precision,layers['jj1']))
    else:
        D_out.add_ref(jj1_inverted_fill_clearout)
    print('jj1 calculation time = '+str(time.time()-t_jj1)+' s\n')
    
    #jj2
    print('jj2...')
    t_jj2 = time.time()
    jj2 = pg.extract(D_in,[layers['jj2']])
    jj2_inverted = pg.boolean(chip_edge,jj2,'A-B',precision = precision,layer = layers['jj2'])
    jj2_inverted_fill_clearout = pg.boolean(jj2_inverted,fill_region_jj2,'A-B',precision = precision,layer = layers['jj2'])
    if params['post_processing_elim_small_features'] == True:
        D_out.add_ref(vt_eliminate_small_features(jj2_inverted_fill_clearout,params['post_processing_min_dim'],precision,layers['jj2']))
    else:
        D_out.add_ref(jj2_inverted_fill_clearout)
    print('jj2 calculation time = '+str(time.time()-t_jj2)+' s\n')
    
    #resistor
    print('resistor...')
    t_res = time.time()
    res = pg.extract(D_in,[layers['res']])
    res = pg.boolean(res,fill_region_res,'A+B',precision = precision,layer = layers['res'])
    if params['post_processing_elim_small_features'] == True:
        D_out.add_ref(vt_eliminate_small_features(res,params['post_processing_min_dim'],precision,layers['res']))
    else:
        D_out.add_ref(res)
    print('resistor calculation time = '+str(time.time()-t_res)+' s\n')
    
    #m4
    print('m4...')
    t_m4 = time.time()
#    m4 = pg.extract(D_in,[layers['m4'],layers['m4l']])
    m4 = pg.extract(D_in,[layers['m4']])
    m4_inverted = pg.boolean(chip_edge,m4,'A-B',precision = precision,layer = layers['m4'])
    m4_inverted_fill_clearout = pg.boolean(m4_inverted,fill_region_m4,'A-B',precision = 1e-6,layer = layers['m4'])
    if params['post_processing_elim_small_features'] == True:
        m4_inverted_fill_clearout = vt_eliminate_small_features(m4_inverted_fill_clearout,params['post_processing_min_dim'],precision,layers['m4'])
    m4_inverted_fill_clearout = pg.boolean(m4_inverted_fill_clearout,pg.extract(D_in,[layers['m4l']]),'A-B',precision = precision, layer = layers['m4'])
    D_out.add_ref(m4_inverted_fill_clearout)
    print('m4 calculation time = '+str(time.time()-t_m4)+' s\n')
    
    #vias
    print('vias...')
    t_vias = time.time()
    via_layer_list = ['v0','v2','v4']
    D_out.add_ref(vt_add_via_layers(D_in,via_layer_list,layers))
    print('via calculation time = '+str(time.time()-t_vias)+' s\n')
    
    elapsed = time.time() - t_tot
    print('total post processing duration = '+str(elapsed)+' s')
    
    return D_out, layer_list

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