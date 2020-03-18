#%% initialization
import numpy as np
from pylab import *
import os
import shutil
import glob


#%% find all .dat files 
src = os.getcwd()
iv_data_file_list = []
sq_data_file_list = []
for file_name in glob.glob('*_IV_*.dat'):
   iv_data_file_list.append(file_name)

for file_name in glob.glob('*_incoil_*.dat'):
   sq_data_file_list.append(file_name)
   
#%% plot I-Vs
  
for file_name in iv_data_file_list:
    I = []
    V = []
    with open(file_name,'r') as file:
        counter = 0
        for line in file:
            text = line.split(',',1)
            if counter > 0:
                I.append(float(text[0]))
                V.append(float(text[1]))
            counter += 1
        fig, axes = plt.subplots(1,1)
        axes.plot(np.asarray(V)*1e3, np.asarray(I)*1e6, 'o-', linewidth = 1, markersize = 3)#, label = 'I_sq = {0} uA'.format(I_sq_bias[ii]*1e6)
        axes.set_xlabel(r'Voltage across SQUID (mV)', fontsize=20)
        axes.set_ylabel(r'SQUID bias current current (uA)', fontsize=20)
        axes.legend(loc='best')
        grid(True,which='both')
        title(file_name)
        plt.show()
        plt.savefig(file_name[0:-4]+'.png', bbox_inches = 'tight')
        plt.savefig(file_name[0:-4]+'.pdf', bbox_inches = 'tight')
        

#%% plot SQUIDs
id_string = 'sq bias current'
for file_name in sq_data_file_list:
    sq_bias_lines = []
    with open(file_name,'r') as file:
        for num, line in enumerate(file,1):
            if id_string in line:
                sq_bias_lines.append(num)
    num_biases = len(sq_bias_lines)
    num_per_bias = sq_bias_lines[1]-sq_bias_lines[0]-3
    I_sq = np.zeros((num_biases,1))
    I_incoil = np.zeros((num_per_bias,num_biases))
    V_sq = np.zeros((num_per_bias,num_biases))
    with open(file_name,'r') as file:
        lines = file.readlines()
        for ii in range(num_biases):
            I_sq[ii] = lines[sq_bias_lines[ii]]
            for jj in range(num_per_bias):
                text = lines[sq_bias_lines[ii]+2+jj].split(',',1)
                I_incoil[jj,ii] = float(text[0])
                V_sq[jj,ii] = float(text[1])
    fig, axes = plt.subplots(1,1)
    for ii in range(num_biases):
        I = I_incoil[:,ii]
        V = V_sq[:,ii]
        axes.plot(I*1e6, V*1e3, 'o-', linewidth = 1, markersize = 3, label = 'I_sq = {0} uA'.format(float(I_sq[ii])*1e6))
    axes.set_ylabel(r'Voltage across SQUID (mV)', fontsize=20)
    axes.set_xlabel(r'Incoil current current (uA)', fontsize=20)
    axes.legend(loc='best')
    grid(True,which='both')
    title(file_name)
    plt.show()
    plt.savefig(file_name[0:-4]+'.png', bbox_inches = 'tight')
    plt.savefig(file_name[0:-4]+'.pdf', bbox_inches = 'tight')
       
