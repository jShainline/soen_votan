# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 16:18:51 2019

@author: jms4
"""

from phidl import LayerSet
import numpy as np

#%%
def vt_layers():
    vt_lyrs = LayerSet()
    
    color_mat = gds_colors()
    
    #gds layer name, layer number, data type, description, color from color list, dither
    layer_data = [['m0',40,0,'ground plane',32,'I1'],
                  ['m0l',40,1,'ground plane label',32,'I0'],
                  ['m0i',40,2,'ground plane invert',32,'I1'],
                  ['m0e',40,3,'ground plane endpoint',32,'I1'],
                  ['m0m',40,4,'ground plane moats',32,'I2'],
                  ['v0',50,0,'via from jj1 to ground plane',11,'I5'],
                  ['v0e',50,1,'via from jj1 to ground plane',11,'I5'],
                  ['jj2',22,0,'jj top contact',13,'I5'],
                  ['jj2e',22,1,'jj top contact endpoint',13,'I2'],
                  ['jj2f',22,2,'jj top contact fill',13,'I1'],
                  ['jj1',21,0,'jj bottom contact / m3',14,'I4'],
                  ['jj1e',21,1,'jj bottom contact endpoint',14,'I2'],
                  ['jj1f',21,2,'jj bottom contact fill',14,'I1'],                  
                  ['v2',52,0,'via to jj top and bottom',12,'I1'],
                  ['v2e',52,1,'via to jj endpoint',11,'I5'],
                  ['m4',42,0,'wiring / m4',34,'I9'],
                  ['m4e',42,1,'wiring endpoint',34,'I2'],
                  ['m4f',42,2,'wiring fill',34,'I1'],
                  ['m4l',42,3,'m4 / label layer',34,'I0'],
                  ['res',32,0,'resistor',19,'I9'],
                  ['rf',32,1,'resistor fill',16,'I1'],                  
                  ['v4',54,0,'via to wiring / pad opening',16,'I1'],
                  ['v4e',54,1,'v4 endpoint',11,'I2'],
                  ['ip1',19,0,'inductor ports _ jj1',11,'I1'],
                  ['ip2',19,0,'inductor ports _ m4',11,'I1'],
                  ['ipt',18,0,'inductor port text labels',11,'I1'],
                  ['pl',95,0,'pad locations',44,'I1'],
                  ['ce',96,0,'chip edge',41,'I1'],
                  ['dp',99,10,'data prep dummy',8,'I1'],
                  ]
 

#                  ['m1',41,0,'small pads / m1',1,'I1'],
#                  ['spd',11,0,'nanowire / inductor / m2',2,'I2'],
#                  ['v1',51,0,'via to nw',3,'I3'],
#                  ['v1e',51,1,'via to nw',3,'I2'],                 
    
    num_layers = len(layer_data)    
    for ii in range(num_layers):        
        color_number = layer_data[ii][4]
        color_hex = '#{0:02x}{1:02x}{2:02x}'.format(clamp(color_mat[:,color_number-1][0]*256), clamp(color_mat[:,color_number-1][1]*256), clamp(color_mat[:,color_number-1][2]*256))
        vt_lyrs.add_layer(name = layer_data[ii][0], gds_layer = layer_data[ii][1], gds_datatype = layer_data[ii][2],description = layer_data[ii][3], color = color_hex, inverted = False,alpha = 0.6, dither = layer_data[ii][5])
    
    return vt_lyrs, layer_data

#%%
def vt_layers_post():
    vt_lyrs = LayerSet()
    
    color_mat = gds_colors()
    
    #gds layer name, layer number, data type, description, color from color list, dither
    layer_data = [['m0',40,0,'ground plane',32,'I1'],
                  ['v0',50,0,'via from jj1 to ground plane',11,'I5'],
                  ['jj2',22,0,'jj top contact',13,'I5'],
                  ['jj1',21,0,'jj bottom contact / m3',14,'I4'],                  
                  ['v2',52,0,'via to jj top and bottom',12,'I1'],
                  ['m4',42,0,'wiring / m4',34,'I9'],
                  ['res',32,0,'resistor',19,'I9'],                 
                  ['v4',54,0,'via to wiring / pad opening',16,'I1'],
                  ['ce',96,0,'chip edge',41,'I1']
                  ]              
    
    num_layers = len(layer_data)    
    for ii in range(num_layers):        
        color_number = layer_data[ii][4]
        color_hex = '#{0:02x}{1:02x}{2:02x}'.format(clamp(color_mat[:,color_number-1][0]*256), clamp(color_mat[:,color_number-1][1]*256), clamp(color_mat[:,color_number-1][2]*256))
        vt_lyrs.add_layer(name = layer_data[ii][0], gds_layer = layer_data[ii][1], gds_datatype = layer_data[ii][2],description = layer_data[ii][3], color = color_hex, inverted = False,alpha = 0.6, dither = layer_data[ii][5])
    
    return vt_lyrs, layer_data

#%%
def write_lyp(lyrs,layer_data,lyp_file_name):

    num_layers = len(lyrs._layers)
    color_mat = gds_colors()
    
#    gds_layers = np.zeros([num_layers,1])
#    for ii in range(num_layers):
#        gds_layers[ii] = layer_data[ii][1]
#        
#    index_array,gds_layers_sorted = np.argsort(gds_layers)
    
    A = '<?xml version="1.0" encoding="utf-8"?>\n<layer-properties>'
    
    for kk in range(num_layers):
        
        ii = kk#index_array[kk]
        color_number = layer_data[ii][4]
        color_hex = '#{0:02x}{1:02x}{2:02x}'.format(clamp(color_mat[:,color_number-1][0]*256), clamp(color_mat[:,color_number-1][1]*256), clamp(color_mat[:,color_number-1][2]*256))
        A = A + '<properties>\n<frame-color>'+color_hex+'</frame-color>\n'
        A = A + '<fill-color>'+color_hex+'</fill-color>\n'
        A = A + '<frame-brightness>0</frame-brightness>\n'
        A = A + '<fill-brightness>0</fill-brightness>\n'
        A = A + '<dither-pattern>'+layer_data[ii][5]+'</dither-pattern>\n'
        A = A + '<visible>true</visible>\n'
        A = A + '<transparent>false</transparent>\n'
        A = A + '<width>1</width>\n'
        A = A + '<marked>false</marked>\n'
        A = A + '<animation>0</animation>\n'
        A = A + '<name>'+str(layer_data[ii][1])+'/'+str(layer_data[ii][2])+': '+str(layer_data[ii][0])+'; '+str(layer_data[ii][3])+'</name>\n'
        A = A + '<source>'+str(layer_data[ii][1])+'/'+str(layer_data[ii][2])+'@1'+'</source>\n'
        A = A + '</properties>\n'
    
    A = A + '</layer-properties>'    
    
    print(A,file=open(lyp_file_name+'.lyp','w'))
#    with open('vt.lyp','w') as text_file:
#        text_file.write(A)
  
    return

#%%
def gds_colors():
    
    ## define colors
    #blues  lightest to darkest
    blueVec1 = np.array([145,184,219]); blue1 = blueVec1/256;
    blueVec2 = np.array([96,161,219]); blue2 = blueVec2/256;
    blueVec3 = np.array([24,90,149]); blue3 = blueVec3/256;
    blueVec4 = np.array([44,73,100]); blue4 = blueVec4/256;
    blueVec5 = np.array([4,44,80]); blue5 = blueVec5/256;
    #reds  lightest to darkest
    redVec1 = np.array([246,177,156]); red1=redVec1/256;
    redVec2 = np.array([246,131,98]); red2 = redVec2/256;
    redVec3 = np.array([230,69,23]); red3 = redVec3/256;
    redVec4 = np.array([154,82,61]); red4 = redVec4/256;
    redVec5 = np.array([123,31,4]); red5 = redVec5/256;
    #greens  lightest to darkest
    greenVec1 = np.array([142,223,180]); green1 = greenVec1/256;
    greenVec2 = np.array([89,223,151]); green2 = greenVec2/256;
    greenVec3 = np.array([16,162,84]); green3 = greenVec3/256;
    greenVec4 = np.array([43,109,74]); green4 = greenVec4/256;
    greenVec5 = np.array([3,87,42]); green5 = greenVec5/256;
    #yellows  lightest to darkest
    yellowVec1 = np.array([246,204,156]); yellow1 = yellowVec1/256;
    yellowVec2 = np.array([246,185,98]); yellow2 = yellowVec2/256;
    yellowVec3 = np.array([230,144,23]); yellow3 = yellowVec3/256;
    yellowVec4 = np.array([154,115,61]); yellow4 = yellowVec4/256;
    yellowVec5 = np.array([123,74,4]); yellow5 = yellowVec5/256;
    
    #blue grays
    gBlueVec1 = np.array([197,199,202]); gBlue1 = gBlueVec1/256;
    gBlueVec2 = np.array([195,198,202]); gBlue2 = gBlueVec2/256;
    gBlueVec3 = np.array([142,145,149]); gBlue3 = gBlueVec3/256;
    gBlueVec4 = np.array([108,110,111]); gBlue4 = gBlueVec4/256;
    gBlueVec5 = np.array([46,73,97]); gBlue5 = gBlueVec5/256;
    #red grays
    gRedVec1 = np.array([242,237,236]); gRed1 = gRedVec1/256;
    gRedVec2 = np.array([242,235,233]); gRed2 = gRedVec2/256;
    gRedVec3 = np.array([230,231,218]); gRed3 = gRedVec3/256;
    gRedVec4 = np.array([172,167,166]); gRed4 = gRedVec4/256;
    gRedVec5 = np.array([149,88,71]); gRed5 = gRedVec5/256;
    #green grays
    gGreenVec1 = np.array([203,209,206]); gGreen1 = gGreenVec1/256;
    gGreenVec2 = np.array([201,209,204]); gGreen2 = gGreenVec2/256;
    gGreenVec3 = np.array([154,162,158]); gGreen3 = gGreenVec3/256;
    gGreenVec4 = np.array([117,122,119]); gGreen4 = gGreenVec4/256;
    gGreenVec5 = np.array([50,105,76]); gGreen5 = gGreenVec5/256;
    #yellow grays
    gYellowVec1 = np.array([242,240,236]); gYellow1 = gYellowVec1/256;
    gYellowVec2 = np.array([242,239,233]); gYellow2 = gYellowVec2/256;
    gYellowVec3 = np.array([230,225,218]); gYellow3 = gYellowVec3/256;
    gYellowVec4 = np.array([172,169,166]); gYellow4 = gYellowVec4/256;
    gYellowVec5 =np.array( [149,117,71]); gYellow5 = gYellowVec5/256;
    
    #pure grays (white to black)
    gVec1 = np.array([256,256,256]); g1 = gVec1/256;
    gVec2 = np.array([242,242,242]); g2 = gVec2/256;
    gVec3 = np.array([230,230,230]); g3 = gVec3/256;
    gVec4 = np.array([204,204,204]); g4 = gVec4/256;
    gVec5 = np.array([179,179,179]); g5 = gVec5/256;
    gVec6 = np.array([153,153,153]); g6 = gVec6/256;
    gVec7 = np.array([128,128,128]); g7 = gVec7/256;
    gVec8 = np.array([102,102,102]); g8 = gVec8/256;
    gVec9 = np.array([77,77,77]); g9 = gVec9/256;
    gVec10 = np.array([51,51,51]); g10 = gVec10/256;
    gVec11 = np.array([26,26,26]); g11 = gVec11/256;
    gVec12 = np.array([0,0,0]); g12 = gVec12/256;
    
    color_mat = np.column_stack((blue1,blue2,blue3,blue4,blue5,red1,red2,red3,red4,red5,green1,green2,green3,green4,green5,yellow1,yellow2,yellow3,yellow4,yellow5,
    gBlue1,gBlue2,gBlue3,gBlue4,gBlue5,gRed1,gRed2,gRed3,gRed4,gRed5,gGreen1,gGreen2,gGreen3,gGreen4,gGreen5,gYellow1,gYellow2,gYellow3,gYellow4,gYellow5,
    g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12))
                    
    return color_mat

#%%
def clamp(x): return int(max(0, min(x, 255)))