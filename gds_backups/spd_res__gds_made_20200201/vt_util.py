from phidl import LayerSet
import numpy as np

#%%
def vt_layers():
    vt_lyrs = LayerSet()
    
    color_mat = gds_colors()
    
    #gds layer name, layer number, data type, description, color from color list, dither
    layer_data = [['m1',40,0,'wiring for stf and r1',32,'I3'],
                  ['m1e',40,1,'wiring for stf and r1 endpoint',32,'I2'],
                  ['m1f',40,2,'wiring for stf and r1 fill',32,'I1'],
                  ['m1l',40,3,'wiring for stf and r1 label',32,'I1'],
                  ['m1p',40,4,'wiring for stf and r1 pad',32,'I4'],
                  ['stf',10,0,'superconducting thin film for spds and inductors',3,'I5'],
                  ['stfe',10,1,'superconducting thin film endpoint',3,'I1'],
                  ['stff',10,2,'superconducting thin film fill',3,'I1'],
                  ['stfp',10,3,'superconducting thin film pad',3,'I3'],
                  ['r1',30,0,'spd and loop resistor',18,'I5'],
                  ['r1f',30,2,'spd and loop resistor fill',18,'I1'],
                  ['r1p',30,3,'spd and loop resistor pad',18,'I3'],
                  ['v1',50,0,'via from m1 to m2',11,'I5'],
                  ['v1e',50,1,'v1 endpoint',11,'I1'],                  
                  ['m2',41,0,'ground plane',32,'I1'],
                  ['m2e',41,1,'m2 endpoint',32,'I1'],
                  ['m2l',41,3,'m2 label',32,'I0'],
                  ['m2i',41,4,'m2 invert',32,'I1'],
                  ['m2o',41,5,'m2 offset',32,'I2'],#postprocessing will trace around this layer to allow vias
                  ['m2m',41,6,'m2 moats',32,'I2'],                  
                  ['v2',51,0,'via from m2 to jj1, jj2, or m1',11,'I6'],
                  ['v2e',51,1,'v2 endpoint',11,'I1'],
                  ['jj1',21,0,'jj bottom contact / m3',14,'I4'],
                  ['jj1e',21,1,'jj bottom contact endpoint',14,'I2'],
                  ['jj1f',21,2,'jj bottom contact fill',14,'I1'],                  
                  ['jj2',22,0,'jj top contact',13,'I5'],
                  ['jj2e',22,1,'jj top contact endpoint',13,'I2'],
                  ['jj2f',22,2,'jj top contact fill',13,'I1'],
                  ['v3',52,0,'via to jj top and bottom contacts',11,'I8'],
                  ['v3e',52,1,'v3 endpoint',11,'I5'],
                  ['m3',42,0,'jj contact metal',34,'I9'],
                  ['m3e',42,1,'jj contact metal endpoint',34,'I1'],
                  ['m3f',42,2,'jj contact metal fill',34,'I1'],
                  ['m3l',42,3,'jj contact metal label',34,'I1'],
                  ['m3p',42,4,'jj contact metal pad',34,'I4'],
                  ['m3cs',43,6,'m3 label',34,'I0'],
                  ['r2',31,0,'jj shunt resistor',8,'I9'],
                  ['r2f',31,2,'jj shunt resistor fill',8,'I1'],
                  ['r2p',31,3,'jj shunt resistor pad',8,'I3'],
                  ['v4',54,0,'via to r2 / pad opening',16,'I1'],
                  ['v4e',54,1,'v4 endpoint',16,'I2'],                  
                  ['pkg',60,0,'SU8 packaging layer',1,'I1'],                  
                  ['pkfc',60,1,'fiber core dummy layer',21,'I1'],
                  ['ipm1',19,0,'inductor port m1',11,'I1'],
                  ['ipj1',19,0,'inductor port jj1',11,'I1'],
                  ['ipm3',19,0,'inductor port m3',11,'I1'],
                  ['ipl',19,0,'inductor port labels',11,'I1'],
                  ['pl',95,0,'pad locations',47,'I1'],
                  ['ce',96,0,'chip edge',48,'I1'],
                  ['dp',99,10,'data prep dummy',8,'I1'],
                  ]
    
#                  ['v4',53,0,'via from m4 to m3',11,'I9'],
#                  ['v4e',53,1,'via from m4 to m3 endpoint',11,'I1'],                  
#                  ['m4',43,0,'upper metal wiring',34,'I9'],
#                  ['m4e',43,1,'m4 endpoint',34,'I2'],
#                  ['m4f',43,2,'m4 fill',34,'I1'],
#                  ['m4l',43,3,'m4 label',34,'I0'],
#                  ['m4cs',43,6,'m4 label',34,'I0'],
#                  ['r3',32,0,'resistor / pad cap',19,'I9'],
#                  ['r3f',32,1,'resistor / pad cap fill',16,'I1'],
#                  ['v5',54,0,'via to r3 / pad opening',16,'I1'],
#                  ['v5e',54,1,'v4 endpoint',11,'I2'],
                
    
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
    layer_data = [['m1',40,0,'wiring for stf and r1',32,'I9'],
                  ['stf',10,0,'superconducting thin film for spds and inductors',3,'I5'],
                  ['r1',30,0,'spd and loop resistor',18,'I5'],
                  ['v1',50,0,'via from m1 to m2',11,'I5'],                 
                  ['m2',41,0,'ground plane',32,'I1'],                 
                  ['v2',51,0,'via from m2 to jj1, jj2, or m1',11,'I6'],
                  ['jj1',21,0,'jj bottom contact / m3',14,'I4'],                 
                  ['jj2',22,0,'jj top contact',13,'I5'],
                  ['v3',52,0,'via to jj top and bottom contacts',11,'I8'],
                  ['m3',42,0,'jj contact metal',34,'I9'],
                  ['r2',31,0,'jj shunt resistor',8,'I9'],
                  ['v4',54,0,'via to r2 / pad opening',16,'I1'],             
                  ['pkg',60,0,'SU8 packaging layer',1,'I1'], 
                  ['ce',96,0,'chip edge',48,'I1'],
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