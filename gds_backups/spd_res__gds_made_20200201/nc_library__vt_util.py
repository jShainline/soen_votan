import copy
import numpy as np
import gdspy
import phidl, phidl.geometry as pg, phidl.routing as pr
from phidl import Device, Layer, LayerSet, Port, quickplot2 as qp
from phidl import make_device
#from nc_library import *

from vt_util import vt_layers
vt_lyrs,layer_data = vt_layers()

def vt_arg_helper(params,parameter_name,default_value):
    
    if parameter_name in params.keys(): 
        value = params[parameter_name] 
    else: 
        value = default_value
    
    return value


def corner_fixer(w_wire = 1, layer = vt_lyrs['m3']):
    
    D_cf = Device('corner fixer')
    
    Sq = pg.rectangle(size = [w_wire,w_wire], layer = layer)
    D_cf.add_ref(Sq)    
    D_cf.add_port(name = 'north', midpoint = [w_wire/2,w_wire], width = w_wire, orientation = 90)
    D_cf.add_port(name = 'south', midpoint = [w_wire/2,0], width = w_wire, orientation = 270)
    D_cf.add_port(name = 'east', midpoint = [w_wire,w_wire/2], width = w_wire, orientation = 0)
    D_cf.add_port(name = 'west', midpoint = [0,w_wire/2], width = w_wire, orientation = 180)
    
    return D_cf


def vt_fill(params = dict()):
    
    fill_w_widest = vt_arg_helper(params,'fill_w_widest',40)
    fill_w_delta = vt_arg_helper(params,'fill_w_delta',4)
    fill_gap = vt_arg_helper(params,'fill_gap',10)
    fill_layer_list = vt_arg_helper(params,'fill_layer_list',['jj1f','jj2f','m3f'])
    chip_size = vt_arg_helper(params,'chip_size',[5000,5000])
    layers = vt_arg_helper(params,'layers',vt_lyrs)    

    D_fill = Device('fill')
    
    Boxes = Device('boxes')
    for ii in range(len(fill_layer_list)):
        Box = pg.rectangle(size = [fill_w_widest-ii*2*fill_w_delta,fill_w_widest-ii*2*fill_w_delta], layer = layers[fill_layer_list[ii]])
        box = Boxes.add_ref(Box)
        box.center = [0,0]
    
    period = fill_w_widest+fill_gap
    num_x = np.floor((chip_size[0]+fill_gap)/period)
    num_y = np.floor((chip_size[1]+fill_gap)/period)
    
    Row = Device('fill_row')
    for ii in range(num_x.astype(int)):
        row_element = Row.add_ref(Boxes)
        row_element.movex(ii*period)
    
    for ii in range(num_y.astype(int)):
        col_element = D_fill.add_ref(Row)
        col_element.movey(ii*period)
        
    D_fill.center = [0,0]
    
    return D_fill


def vt_fill_circles(fill_w_widest = 40, fill_w_delta = 4, fill_gap = 10, fill_layer_list = ['jj1f','jj2f','m3f'], chip_size = [5000,5000], layers = vt_lyrs, **kwargs):

    D_fill = Device('fill')
    
    Boxes = Device('circles')
    for ii in range(len(fill_layer_list)):
        Box = pg.circle(radius = (fill_w_widest-ii*2*fill_w_delta)/2, angle_resolution = 3.6, layer = layers[fill_layer_list[ii]])
        box = Boxes.add_ref(Box)
        box.center = [0,0]
    
    period = fill_w_widest+fill_gap
    num_x = np.floor((chip_size[0]+fill_gap)/period)
    num_y = np.floor((chip_size[1]+fill_gap)/period)
    
    Row = Device('fill_row')
    for ii in range(num_x.astype(int)):
        row_element = Row.add_ref(Boxes)
        row_element.movex(ii*period)
    
    for ii in range(num_y.astype(int)):
        col_element = D_fill.add_ref(Row)
        col_element.movey(ii*period)
        
    D_fill.center = [0,0]
    
    return D_fill


def vt_endpoint_boxes(params = dict()):
    
    endpoint_box_width = vt_arg_helper(params,'endpoint_box_width',300)
    pad_label_size = vt_arg_helper(params,'pad_label_size',10)
    block_label_size = vt_arg_helper(params,'block_label_size',20)
    label_layer = vt_arg_helper(params,'label_layer','m3l')
    layers = vt_arg_helper(params,'layers',vt_lyrs)
        
    D_epb = Device('endpoint boxes')
    
    outer_layer_list = ['m1e','stfe','m1' ,'m2e','m2' ,'jj2e','jj1e','jj1'     ,'m3e','r2']
    inner_layer_list = [''   ,''    ,'v1e',''   ,'v2e',''    ,''    ,'v3e'     ,''   ,'v4e']
    label_list =       ['m1','stf'  ,'v1' ,'m2' ,'v2' ,'jj2' ,'jj1' ,'v3 / jj1','m3' ,'v4 / r2']
    inner_box_inset = 2
    w = endpoint_box_width
    x_shift_factor = 1.2
    for ii in range(len(outer_layer_list)):
        Epb = pg.rectangle(size = [w,w], layer = layers[outer_layer_list[ii]])
        epb = D_epb.add_ref(Epb)
        epb.movex(ii*x_shift_factor*w)
        Text_label = vt_label_maker(text_string = label_list[ii], text_size = pad_label_size, layer = layers[label_layer])
        text_label = D_epb.add_ref(Text_label)
        text_label.center = [ii*x_shift_factor*w+w/2,w+2*pad_label_size]       
        if inner_layer_list[ii] != '':
            Epb = pg.rectangle(size = [w-2*inner_box_inset,w-2*inner_box_inset], layer = layers[inner_layer_list[ii]])
            epb = D_epb.add_ref(Epb)
            epb.move([ii*x_shift_factor*w+inner_box_inset,inner_box_inset])                
    D_epb.center = [0,0]
    
    #label
#    D_epb.add_port(name = 'label_port', midpoint = [0,D_epb.ymax+0.5*block_label_size], width = 0, orientation = 90)
#    Text_label = vt_label_maker(text_string = 'endpoint', text_size = block_label_size, layer = layers[label_layer])
#    text_label = D_epb.add_ref(Text_label)
#    text_label.connect(port = 'south', destination = D_epb.ports['label_port'])

    return D_epb


def vt_inline_test(params = dict()):
    
    inline_test_pad_size = vt_arg_helper(params,'inline_test_pad_size',[200,300])
    pad_label_size = vt_arg_helper(params,'pad_label_size',10)
    block_label_size = vt_arg_helper(params,'block_label_size',40)
    pad_label_size = vt_arg_helper(params,'pad_label_size',20)
    layers = vt_arg_helper(params,'layers',vt_lyrs)
    
    from nc_library import wire_basic
    from nc_library__vt_pads_vias_wires import jj_pad_variable_multi

    D_inline = Device('inline_tests')   
    
    wire_length = 100
    wire_width_vec = [2,20]
    pad_opening_inset = 2
    pad_label_size = 20
    
    pad_layer_sets = [['m1','v1'],['m1','v1'],['m1','v1'],
                      ['m2i','v2'],
                      ['jj1','v3'],['jj1','jj2','v3'],
                      ['m3','v4'],['m3','r2','v4']]
    wire_layers = ['m1','stf','r1',
                   'm2i',
                   'jj1','jj1',
                   'm3','r2']
        
    metals_label_list = ['m1','stf','r1','m2','jj1','jj1','m3','m3']
    label_layers = ['m1l','m1l','m1l','m2l','m2l','m2l','m3l','m3l']
    for ii in range(len(pad_layer_sets)):        
        for kk in range(len(wire_width_vec)):
            Test_device = Device('inline_test')
            params_mod = copy.deepcopy(params)            
            params_mod['pad_size'] = inline_test_pad_size
            params_mod['pad_opening_inset'] = pad_opening_inset
            params_mod['pad_w_wire'] = wire_width_vec[kk]
            params_mod['pad_variable_layers'] = pad_layer_sets[ii]
            Pad = jj_pad_variable_multi(params_mod)
            pad1 = Test_device.add_ref(Pad)
            pad2 = Test_device.add_ref(Pad)
            pad2.connect(port = 'north', destination = pad1.ports['south'])
            pad2.movey(-wire_length)
            if wire_layers[ii] == 'r1' or wire_layers[ii] == 'stf' or wire_layers[ii] == 'r2':
                wire1 = wire_basic([pad2.ports['north'].midpoint[0],pad2.ports['north'].midpoint[1]-10],[pad1.ports['south'].midpoint[0],pad1.ports['south'].midpoint[1]+10],'y',wire_width_vec[kk],layers[wire_layers[ii]])            
            else:
                wire1 = wire_basic(pad2.ports['north'],pad1.ports['south'],'y',wire_width_vec[kk],layers[wire_layers[ii]])            
            Test_device.add_ref(wire1)
            string = 'metal: '+metals_label_list[ii]+'; via: '+pad_layer_sets[ii][-1]+'; wire: '+wire_layers[ii]
            Text_label = vt_label_maker(text_string = string, text_size = pad_label_size, layer = layers[label_layers[ii]])
            text_label = Test_device.add_ref(Text_label)
            text_label.connect(port = 'south_west', destination = pad2.ports['text'])
            test1 = D_inline.add_ref(Test_device)
            test1.move([ii*2*inline_test_pad_size[0],kk*3530]) 
#            test1.move([ii*1.5*inline_test_pad_size[0],kk*(2*inline_test_pad_size[1]+1.5*wire_length)]) 
    
    #label
    D_inline.center = [0,0]
#    D_inline.add_port(name = 'label_port', midpoint = [0,D_inline.ymax+block_label_size], width = 0, orientation = 90)
#    Text_label = vt_label_maker(text_string = 'in-line tests', text_size = block_label_size, layer = layers['m1l'])
#    text_label = D_inline.add_ref(Text_label)
#    text_label.connect(port = 'south', destination = D_inline.ports['label_port'])

    return D_inline


def vt_litho_tests(params = dict()):
    
    jj_junc_diam_vec = vt_arg_helper(params,'jj_junc_diam_vec',10)
    layers = vt_arg_helper(params,'layers',vt_lyrs)
    
    D_litho_tests = Device('litho_tests')
    D_litho_tests_inv = Device('litho_tests_tone_inverted')

    size_list = jj_junc_diam_vec
    text_size = 4
    
    x_padding = 12
    y_padding = 10
    
    #isolated squares    
    layer_list = ['stf','r1','jj1','jj2','r2'] 
    x = 0
    y = 0    
    D_litho_tests_squares = Device('litho_tests_squares')
    for ii in range(len(layer_list)):
        Rec_col = Device('litho_rectangle_column')
        for jj in range(len(size_list)):             
            Rec_inst = Device('litho_rectangle_instance')
            
            Rec = pg.rectangle(size = [size_list[jj],size_list[jj]], layer = layers[layer_list[ii]])
            Rec.add_port(name = 'west', midpoint = [Rec.xmin,Rec.y], width = 0, orientation = 180)
            rec = Rec_inst.add_ref(Rec)
            rec.center = [x,y]
            
            text_str = layer_list[ii]+' '+str(size_list[jj])
            Tex = pg.text(text = text_str, size = text_size, justify = 'right', layer = layers[layer_list[ii]])
            Tex.add_port(name = 'east', midpoint = [Tex.xmax,Tex.y], width = 0, orientation = 0)
            tex = Rec_inst.add_ref(Tex)            
            tex.connect(port = 'east', destination = rec.ports['west'])            
            tex.movex(-text_size)
            
            Rec_col.add_ref(Rec_inst)
            
            y = y+Rec_inst.ysize+y_padding
        
        x = x+Rec_col.xsize+x_padding            
        y = 0
        D_litho_tests_squares.add_ref(Rec_col)
    isolated_rectangles = D_litho_tests.add_ref(D_litho_tests_squares)
    
    #isolated circles
    layer_list = ['stf','r1','jj1','jj2','r2']
    x = 0
    y = 0    
    D_litho_tests_circles = Device('litho_tests_circles') 
    for ii in range(len(layer_list)):
        Cir_col = Device('litho_circle_column')
        for jj in range(len(size_list)):              
            Cir_inst = Device('litho_circle_instance')
            
            Cir = pg.circle(radius = size_list[jj]/2, angle_resolution = 3.6, layer = layers[layer_list[ii]])
            Cir.add_port(name = 'west', midpoint = [Cir.xmin,Cir.y], width = 0, orientation = 180)
            cir = Cir_inst.add_ref(Cir)
            cir.center = [x,y]
                        
            text_str = layer_list[ii]+' '+str(size_list[jj])
            Tex = pg.text(text = text_str, size = text_size, justify = 'right', layer = layers[layer_list[ii]])
            Tex.add_port(name = 'east', midpoint = [Tex.xmax,Tex.y], width = 0, orientation = 0)
            tex = Cir_inst.add_ref(Tex)            
            tex.connect(port = 'east', destination = cir.ports['west'])
            tex.movex(-text_size)
            
            Cir_col.add_ref(Cir_inst)
            
            y = y+Cir_inst.ysize+y_padding
            
        x = x+Cir_col.xsize+x_padding
        y = 0
        D_litho_tests_circles.add_ref(Cir_col)
    isolated_circles = D_litho_tests.add_ref(D_litho_tests_circles)    
        
    #arrays of squares
    array_size = 8
    D_litho_tests_square_arrays = Device('litho_tests_square_arrays') 
    x = 0
    y = 0
    for ii in range(len(layer_list)):
        Rec_col = Device('litho_rect_array_column')
        for jj in range(len(size_list)):
            Rec_inst = Device('litho_rect_array_inst')
            
            Rec = vt_fill(fill_w_widest = size_list[jj], fill_w_delta = 0, fill_gap = size_list[jj], fill_layer_list = [layer_list[ii]], chip_size = [2*array_size*size_list[jj],2*array_size*size_list[jj]], layers = layers)
            Rec.add_port(name = 'west', midpoint = [Rec.xmin,Rec.y], width = 0, orientation = 180)
            rec = Rec_inst.add_ref(Rec)            
            rec.center = [x,y]
            
            text_str = layer_list[ii]+' '+str(size_list[jj])
            Tex = pg.text(text = text_str, size = text_size, justify = 'right', layer = layers[layer_list[ii]])
            Tex.add_port(name = 'east', midpoint = [0,Tex.size[1]/2], width = 0, orientation = 0)
            tex = Rec_inst.add_ref(Tex)            
            tex.connect(port = 'east', destination = rec.ports['west'])
            tex.movex(-text_size)
            
            Rec_col.add_ref(Rec_inst)
            
            y = y+Rec_inst.ysize+y_padding
            
        x = x+Rec_col.xsize+x_padding
        y = 0
        D_litho_tests_square_arrays.add_ref(Rec_col)
    square_arrays = D_litho_tests.add_ref(D_litho_tests_square_arrays)
    
    #arrays of circles
    D_litho_tests_circle_arrays = Device('litho_tests_circle_arrays') 
    x = 0
    y = 0
    for ii in range(len(layer_list)):
        Cir_col = Device('litho_circle_array_column')
        for jj in range(len(size_list)):
            Cir_inst = Device('litho_circle_array_inst')
            
            Cir = vt_fill_circles(fill_w_widest = size_list[jj], fill_w_delta = 0, fill_gap = size_list[jj], fill_layer_list = [layer_list[ii]], chip_size = [2*array_size*size_list[jj],2*array_size*size_list[jj]], layers = layers)
            Cir.add_port(name = 'west', midpoint = [Cir.xmin,Cir.y], width = 0, orientation = 180)
            cir = Cir_inst.add_ref(Cir)            
            cir.center = [x,y]
            
            text_str = layer_list[ii]+' '+str(size_list[jj])
            Tex = pg.text(text = text_str, size = text_size, justify = 'right', layer = layers[layer_list[ii]])
            Tex.add_port(name = 'east', midpoint = [0,Tex.size[1]/2], width = 0, orientation = 0)
            tex = Cir_inst.add_ref(Tex)            
            tex.connect(port = 'east', destination = cir.ports['west'])
            tex.movex(-text_size)
            
            Cir_col.add_ref(Cir_inst)
            
            y = y+Cir_inst.ysize+y_padding
            
        x = x+Cir_col.xsize+x_padding
        y = 0
        D_litho_tests_circle_arrays.add_ref(Cir_col)
    circle_arrays = D_litho_tests.add_ref(D_litho_tests_circle_arrays)

    #isolated lines
    D_litho_tests_isolated_lines = Device('litho_tests_isolated_lines') 
    x = 0
    y = 0
    layer_list = ['m1','stf','r1','m3','r2']
    y_extent = 10*size_list[-1]
    x_extent = 100
    for ii in range(len(layer_list)):
        Line_col = Device('litho_isolated_line_column')
        for jj in range(len(size_list)):
            Line_inst = Device('litho_isolated_line_inst')
            
            Line = pg.rectangle(size = [x_extent,size_list[jj]], layer = layers[layer_list[ii]])
            Line.add_port(name = 'west', midpoint = [Line.xmin,Line.y], width = 0, orientation = 180)
            line = Line_inst.add_ref(Line)            
            line.center = [x,y]
            
            text_str = layer_list[ii]+' '+str(size_list[jj])
            Tex = pg.text(text = text_str, size = text_size, justify = 'right', layer = layers[layer_list[ii]])
            Tex.add_port(name = 'east', midpoint = [0,Tex.size[1]/2], width = 0, orientation = 0)
            tex = Line_inst.add_ref(Tex)            
            tex.connect(port = 'east', destination = line.ports['west'])
            tex.movex(-text_size)
            
            Line_col.add_ref(Line_inst)
            
            y = y+Line_inst.ysize+y_padding
            
        x = x+Line_col.xsize+x_padding
        y = 0
        D_litho_tests_isolated_lines.add_ref(Line_col)
    isolated_lines = D_litho_tests.add_ref(D_litho_tests_isolated_lines)

    #line arrays
    D_litho_tests_line_arrays = Device('litho_tests_line_arrays') 
    x = 0
    y = 0
    layer_list = ['m1','stf','r1','m3','r2']
    for ii in range(len(layer_list)):
        Line_col = Device('litho_isolated_line_column')
        for jj in range(len(size_list)):
            Line_inst = Device('litho_line_array_inst')
            
            Line = pg.rectangle(size = [x_extent,size_list[jj]], layer = layers[layer_list[ii]])
            Line_bundle = Device('litho_test_line_bundle')            
            for kk in range(array_size):
                line = Line_bundle.add_ref(Line)
                line.center = [x,y+2*kk*size_list[jj]]
            Line_bundle.add_port(name = 'west', midpoint = [Line_bundle.xmin,Line_bundle.y], width = 0, orientation = 180)            
            line_bundle = Line_inst.add_ref(Line_bundle)    
#                rec = Rec_bundle.add_ref(Rec)
#                rec.movex(2*kk*size_list[jj])
            
            text_str = layer_list[ii]+' '+str(size_list[jj])
            Tex = pg.text(text = text_str, size = text_size, justify = 'right', layer = layers[layer_list[ii]])
            Tex.add_port(name = 'east', midpoint = [0,Tex.size[1]/2], width = 0, orientation = 0)
            tex = Line_inst.add_ref(Tex)            
            tex.connect(port = 'east', destination = line_bundle.ports['west'])
            tex.movex(-text_size)
            
            Line_col.add_ref(Line_inst)
            
            y = y+Line_inst.ysize+y_padding
            
        x = x+Line_col.xsize+x_padding
        y = 0
        D_litho_tests_line_arrays.add_ref(Line_col)
    line_arrays = D_litho_tests.add_ref(D_litho_tests_line_arrays)    


    #placement of sub-structures
    square_arrays.center = isolated_rectangles.center
    square_arrays.movey(-isolated_rectangles.ysize/2-square_arrays.ysize/2-4*y_padding)
    circle_arrays.center = square_arrays.center
    circle_arrays.movex(square_arrays.xsize/2+circle_arrays.xsize/2+4*x_padding)
    isolated_circles.center = circle_arrays.center
    isolated_circles.movey(circle_arrays.ysize/2+isolated_circles.ysize/2+4*y_padding)
    isolated_lines.center = [square_arrays.xmax+(circle_arrays.xmin-square_arrays.xmax)/2,isolated_rectangles.ymax]
    isolated_lines.movey(isolated_lines.ysize/2+4*y_padding) 
    line_arrays.center = [square_arrays.xmax+(circle_arrays.xmin-square_arrays.xmax)/2,square_arrays.ymin]
    line_arrays.movey(-line_arrays.ysize/2-4*y_padding)

        
    #invert
    D_litho_tests.center = [0,0]
    bounding_box_size = [D_litho_tests.size[0]+10,D_litho_tests.size[1]+10]
    for ii in range(len(layer_list)):
        layer = layers[layer_list[ii]]
        D_single_layer = pg.extract(D_litho_tests, [layer])
        B_box = pg.rectangle(size = bounding_box_size, layer = layer)
        B_box.center = [0,0]
        D_inv = pg.boolean(B_box, D_single_layer, 'A-B', precision = 1e-6, layer = layer)
        D_litho_tests_inv.add_ref(D_inv)
    D_litho_tests_inv.center = [0,0]
        
    return D_litho_tests, D_litho_tests_inv


def vt_label_maker(text_string = 'votan', text_size = 10, justify = 'left', layer = 94):
    
    D_text = Device('text_label')
    
    Text = pg.text(text = text_string, size = text_size, justify = justify, layer = layer)
    D_text.add_ref(Text)
    D_text.add_port(name = 'north', midpoint = [D_text.center[0], D_text.ymax+text_size/2+1], width = 0, orientation = 90) 
    D_text.add_port(name = 'east', midpoint = [D_text.xmax+text_size/2, D_text.center[1]], width = 0, orientation = 0)
    D_text.add_port(name = 'south', midpoint = [D_text.center[0], D_text.ymin-text_size/2-1], width = 0, orientation = 270)
    D_text.add_port(name = 'west', midpoint = [D_text.xmin-text_size/2, D_text.center[1]], width = 0, orientation = 180)
    D_text.add_port(name = 'north_west', midpoint = [D_text.xmin-text_size/2, D_text.ymax+text_size/2+1], width = 0, orientation = 90)
    D_text.add_port(name = 'north_east', midpoint = [D_text.xmax+text_size/2, D_text.ymax+text_size/2+1], width = 0, orientation = 90)
    D_text.add_port(name = 'south_west', midpoint = [D_text.xmin-text_size/2, D_text.ymin-text_size/2-1], width = 0, orientation = 270)
    D_text.add_port(name = 'south_east', midpoint = [D_text.xmax+text_size/2, D_text.ymin-text_size/2-1], width = 0, orientation = 270)
    
    return D_text


def vt_label_chips(params = dict()):
    
    D_text = Device('chip_labels')

    # nist logo / phi labels
    nist_coords1 = [-2270,2395]
    nist_coords2 = [2255,2395]
    D_nist = pg.import_gds('nist_small.gds', cellname = None, flatten = True)
    D_nist.remap_layers(layermap = {0:vt_lyrs[params['label_layer']]})
    D_nist.add_port(name = 'south', midpoint = [D_nist.center[0], D_nist.ymin], width = 0, orientation = 270)
    D_nist.add_port(name = 'south_west', midpoint = [D_nist.xmin, D_nist.ymin], width = 0, orientation = 270)
    D_nist.add_port(name = 'south_east', midpoint = [D_nist.xmax, D_nist.ymin], width = 0, orientation = 270)
    
    nist_a1 = D_text.add_ref(D_nist)
    nist_a1.center = nist_coords1
    nist_a2 = D_text.add_ref(D_nist)
    nist_a2.center = nist_coords2 
    
    Phi_1 = Device('phi_label_1')
    Phi_2 = Device('phi_label_2')
    phi_text_size = params['block_label_size']
    y_space = phi_text_size/2
    text_str1 = params['group_label_1']
    text_str2 = params['group_label_2']
    
    Tex1 = pg.text(text = text_str1, size = phi_text_size, justify = 'left', layer = vt_lyrs[params['label_layer']])
    Tex2 = pg.text(text = text_str2, size = phi_text_size, justify = 'right', layer = vt_lyrs[params['label_layer']])
    Phi_1.add_ref(Tex1)
    Phi_1.add_port(name = 'north', midpoint = [Phi_1.center[0], Phi_1.ymax], width = 0, orientation = 90) 
    Phi_1.add_port(name = 'north_west', midpoint = [Phi_1.xmin, Phi_1.ymax], width = 0, orientation = 90)
    Phi_2.add_ref(Tex2)
    Phi_2.add_port(name = 'north', midpoint = [Phi_2.center[0], Phi_2.ymax], width = 0, orientation = 90) 
    Phi_2.add_port(name = 'north_east', midpoint = [Phi_2.xmax, Phi_2.ymax], width = 0, orientation = 90)

    tex_a_1 = D_text.add_ref(Phi_1)         
    tex_a_1.connect(port = 'north_west', destination = nist_a1.ports['south_west'])
    tex_a_1.movey(-y_space)
    
    tex_a_2 = D_text.add_ref(Phi_2)         
    tex_a_2.connect(port = 'north_east', destination = nist_a2.ports['south_east'])
    tex_a_2.movey(-y_space)
    
    #votan labels    
    votan_coords = [-2300,0]
    votan_text_size = 100
    
    text_str = params['chip_label']
    Tex = pg.text(text = text_str, size = votan_text_size, layer = vt_lyrs[params['label_layer']])
    tex = D_text.add_ref(Tex) 
    tex.rotate(90)           
    tex.center = votan_coords    
    
    #small chip labels
    text_size = params['block_label_size']
    text_str = params['chip_label_small']
    Tex = pg.text(text = text_str, size = text_size, justify = 'left', layer = vt_lyrs[params['label_layer']])
    tex = D_text.add_ref(Tex)            
    tex.center = [-params['chip_size'][0]/2+tex.xsize/2+45,-params['chip_size'][1]/2+tex.ysize/2+45]
    
    # signature
    signature_text_size = params['block_label_size']
    text_str = params['chip_signature']
    Tex = pg.text(text = text_str, size = signature_text_size, justify = 'right', layer = vt_lyrs[params['label_layer']])
    tex = D_text.add_ref(Tex)            
    tex.center = [params['chip_size'][0]/2-tex.xsize/2-45,-params['chip_size'][1]/2+tex.ysize/2+45]
    
    return D_text


def vt_make_die(D,params = dict()):

    chip_size = vt_arg_helper(params,'chip_size',[5000,5000])
    cleave_street_width = vt_arg_helper(params,'cleave_street_width',50)
    cleave_street_length = vt_arg_helper(params,'cleave_street_length',500)
    cleave_street_layer = vt_arg_helper(params,'cleave_street_layer','m3cs')
    label_layer = vt_arg_helper(params,'label_layer','m2l')
    chip_edge_layer = vt_arg_helper(params,'chip_edge_layer','ce')
#    chip_label = vt_arg_helper(params,'chip_label','vt02')
    layers = vt_arg_helper(params,'layers',vt_lyrs)
    
#    chip_label = params['chip_label_small_jj'],label_size = params['block_label_size'],die_size = params['chip_size']
     
    #chip edge
    D.center = (0,0)
    Chip_edge = pg.rectangle(size = chip_size, layer = layers[chip_edge_layer])
    chip_edge = D.add_ref(Chip_edge)
    chip_edge.center = [0,0]
    
    #cleave streets
    D_cs = Device('cleave_streets')
    x1 = 0
    x2 = cleave_street_length
    x3 = x2
    x4 = cleave_street_width
    x5 = x4
    x6 = x1
    y1 = 0
    y2 = y1
    y3 = -cleave_street_width
    y4 = y3
    y5 = -cleave_street_length
    y6 = y5
    D_cs.add_polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4),(x5,y5),(x6,y6)],layers[cleave_street_layer])
    x_shifts = [-chip_size[0]/2,chip_size[0]/2,chip_size[0]/2,-chip_size[0]/2]
    y_shifts = [chip_size[1]/2,chip_size[1]/2,-chip_size[0]/2,-chip_size[0]/2]
    rotations = [0,270,180,90]
    for ii in range(4):
        cs = D.add_ref(D_cs)
        cs.rotate(rotations[ii])
        cs.move([x_shifts[ii],y_shifts[ii]])
        
    #pad locations    
    from nc_library__vt_pads_vias_wires import pad_locations, jj_pad
    P_loc = pad_locations(params)
    p_loc = D.add_ref(P_loc)
    p_loc.center = [0,0] 
        
    params_mod = copy.deepcopy(params)
    params_mod['is_ground_pad'] = True
    for ii in range(len(params['ground_pad_coords'])):    
        params_mod['pad_size_ground'] = params_mod['ground_pad_size'][ii]
        Pad_ground_master = jj_pad(params_mod)
        pad_ground_master = D.add_ref(Pad_ground_master)
        pad_ground_master.center = params_mod['ground_pad_coords'][ii]
       
    #labels
    D_text = vt_label_chips(params)
    D.add_ref(D_text)

    return D, p_loc


def vt_fiber_collar(params = dict()):
        
    fiber_collar_num_pts = vt_arg_helper(params,'fiber_collar_num_pts',40) 
    fiber_collar_diam = vt_arg_helper(params,'fiber_collar_diam',135) 
    fiber_collar_w_wall = vt_arg_helper(params,'fiber_collar_w_wall',135) 
    fiber_collar_w_channel = vt_arg_helper(params,'fiber_collar_w_channel',40) 
    fiber_collar_l_channel = vt_arg_helper(params,'fiber_collar_l_channel',200)     
    fiber_collar_w_box = vt_arg_helper(params,'fiber_collar_w_box',300)
    fiber_collar_l_box = vt_arg_helper(params,'fiber_collar_l_box',500)
    fiber_collar_outer_corner_length = vt_arg_helper(params,'fiber_collar_corner_length',45)
    fiber_collar_inner_corner_length = vt_arg_helper(params,'fiber_collar_inner_corner_length',25)
    fiber_collar_layer = vt_arg_helper(params,'fiber_collar_layer','pkg')
    fiber_collar_orientation = vt_arg_helper(params,'fiber_collar_orientation',180)
    layers = vt_arg_helper(params,'layers',vt_lyrs)    
    
    D_fc = Device('fiber_collar')
    
    w_wall = fiber_collar_w_wall
    l_tot = 2*w_wall+fiber_collar_diam+fiber_collar_l_channel+fiber_collar_l_box
    w_tot = 2*w_wall+max([fiber_collar_diam,fiber_collar_w_box])
    cl = fiber_collar_outer_corner_length
    l_box = fiber_collar_l_box
    w_box = fiber_collar_w_box
    cl_i = fiber_collar_inner_corner_length
        
    #glue box outer
    Gbo = Device('temp_struct')
    points = [[-cl,0],[0,cl],[0,w_tot-cl],[-cl,w_tot],[-l_tot+cl,w_tot],[-l_tot,w_tot-cl],[-l_tot,cl],[-l_tot+cl,0]]    
    glue_box = Gbo.add_polygon(points,layers[fiber_collar_layer])
    glue_box.center = [0,0]
    
    #glue box inner
    Gbi = Device('temp_struct')
    points = [[-cl_i,0],[0,cl_i],[0,w_box-cl_i],[-cl_i,w_box],[-l_box+cl_i,w_box],[-l_box,w_box-cl_i],[-l_box,cl_i],[-l_box+cl_i,0]]    
    glue_box_i = Gbi.add_polygon(points,layers[fiber_collar_layer])
    glue_box_i.center = [glue_box.xmin,glue_box.y]
    glue_box_i.movex(l_box/2+w_wall)    
    
    #fiber collar
    fc = pg.circle(radius = fiber_collar_diam/2, angle_resolution = 360/fiber_collar_num_pts, layer = layers[fiber_collar_layer])
    fc.center = [glue_box.xmax,glue_box.x]
    fc.movex(-fiber_collar_diam/2-w_wall)
    
    #glue channel
    channel = pg.rectangle(size = [fiber_collar_l_channel+fiber_collar_diam/2,fiber_collar_w_channel], layer = layers[fiber_collar_layer])
    channel.center = fc.center
    channel.movex(-fiber_collar_diam/2-fiber_collar_l_channel/2)
    
    #booleans
    inner_region_1 = pg.boolean(A = Gbi, B = channel, operation = 'A+B', precision = 1e-6, num_divisions = [1,1], layer = layers[fiber_collar_layer])
    inner_region_2 = pg.boolean(A = inner_region_1, B = fc, operation = 'A+B', precision = 1e-6, num_divisions = [1,1], layer = layers[fiber_collar_layer])
    inner_region_3 = pg.boolean(A = Gbo, B = inner_region_2, operation = 'A-B', precision = 1e-6, num_divisions = [1,1], layer = layers[fiber_collar_layer])
    
    D_fc.add_ref(inner_region_3)
#    D_fc.add_ref(Gbo)
#    D_fc.add_ref(Gbi)
#    D_fc.add_ref(fc)
#    D_fc.add_ref(channel)
    
    D_fc.add_port(name = 'fiber_center', midpoint = [D_fc.xmax-w_wall-fiber_collar_diam/2,D_fc.y], width = fiber_collar_diam, orientation = -fiber_collar_orientation)
    
    return D_fc